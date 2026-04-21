#include "report_writer.h"
#include <iostream>
#include <stdexcept>
#include <cstdio>
#include <sys/stat.h>
#include <errno.h>
#include <cstring>

bool ReportWriter::ensure_dir(const std::string& path) {
    struct stat st{};
    if (stat(path.c_str(), &st) == 0) return S_ISDIR(st.st_mode);
    if (mkdir(path.c_str(), 0755) == 0) return true;
    if (errno == EEXIST) return true;
    return false;
}

ReportWriter::SplitSink* ReportWriter::get_or_open_split(const std::string& sample) {
    std::lock_guard<std::mutex> lock(split_map_mtx_);
    auto it = split_files_.find(sample);
    if (it != split_files_.end()) return it->second.get();

    std::string path = config_.split_dir + "/" + sample + ".fastq.gz";
    // Compression level 1 favors throughput when split output is enabled.
    gzFile gz = gzopen(path.c_str(), "wb1");
    if (!gz) throw std::runtime_error("Cannot open split file: " + path);
    gzbuffer(gz, 1 << 20);

    auto sink = std::make_unique<SplitSink>();
    sink->gz = gz;
    SplitSink* ptr = sink.get();
    split_files_[sample] = std::move(sink);
    return ptr;
}

void ReportWriter::write_fastq_gz(gzFile gz, const FastqRecord& rec) {
    gzputc(gz, '@');
    gzwrite(gz, rec.name.data(),     static_cast<unsigned>(rec.name.size()));
    gzputc(gz, '\n');
    gzwrite(gz, rec.sequence.data(), static_cast<unsigned>(rec.sequence.size()));
    gzwrite(gz, "\n+\n", 3);
    gzwrite(gz, rec.quality.data(),  static_cast<unsigned>(rec.quality.size()));
    gzputc(gz, '\n');
}

ReportWriter::ReportWriter(const ReporterConfig& cfg) : config_(cfg) {
    if (!cfg.split_dir.empty())
        if (!ensure_dir(cfg.split_dir))
            throw std::runtime_error("Cannot create split directory: " + cfg.split_dir);
    if (!open_tsv_output())
        throw std::runtime_error("Cannot open TSV output: " + cfg.output_file);
    write_tsv_header();
}

ReportWriter::~ReportWriter() {
    std::lock_guard<std::mutex> map_lock(split_map_mtx_);
    for (auto& [name, sink] : split_files_)
        if (sink && sink->gz) gzclose(sink->gz);
}

bool ReportWriter::open_tsv_output() {
    if (config_.output_file.empty() || config_.output_file == "-") {
        tsv_out_ = &std::cout;
    } else {
        tsv_fout_.open(config_.output_file);
        if (!tsv_fout_.is_open()) return false;
        tsv_fout_.rdbuf()->pubsetbuf(tsv_buf_, sizeof(tsv_buf_));
        tsv_out_ = &tsv_fout_;
    }
    return true;
}

void ReportWriter::write_tsv_header() {
    *tsv_out_ <<
        "read_id\t"
        "status\t"
        "strand\t"
        "anchor1_pos\t"
        "anchor2_pos\t"
        "extracted_cb\t"
        "extracted_umi\t"
        "matched_cb\t"
        "edit_distance\t"
        "n_candidates\t"
        "sample\n";
}

void ReportWriter::append_tsv_line(std::string& out, const ReportEntry& e) {
    static constexpr const char* UNASSIGNED = "unassigned";

    out.append(e.read_id);
    out.push_back('\t');
    out.append(status_str(e.status));
    out.push_back('\t');
    out.append(e.strand);
    out.push_back('\t');
    out.append(std::to_string(e.anchor1_pos));
    out.push_back('\t');
    out.append(std::to_string(e.anchor2_pos));
    out.push_back('\t');
    out.append(e.extracted_cb);
    out.push_back('\t');
    out.append(e.extracted_umi);
    out.push_back('\t');
    out.append(e.matched_cb.empty() ? UNASSIGNED : e.matched_cb);
    out.push_back('\t');
    out.append(std::to_string(e.edit_distance));
    out.push_back('\t');
    out.append(std::to_string(e.n_candidates));
    out.push_back('\t');
    out.append(e.sample.empty() ? UNASSIGNED : e.sample);
    out.push_back('\n');
}

void ReportWriter::add_entries(std::vector<ReportChunkItem>&& items) {
    if (items.empty()) return;

    int anchors = 0;
    int matched = 0;
    int ambig   = 0;

    std::string tsv_chunk;
    tsv_chunk.reserve(items.size() * 128);
    for (const auto& item : items) {
        const auto& e = item.entry;
        if (e.anchors_found) ++anchors;
        if (e.cb_matched)    ++matched;
        if (e.ambiguous)     ++ambig;
        append_tsv_line(tsv_chunk, e);
    }

    reads_with_anchors_.fetch_add(anchors, std::memory_order_relaxed);
    reads_matched_.fetch_add(matched, std::memory_order_relaxed);
    ambiguous_calls_.fetch_add(ambig, std::memory_order_relaxed);

    {
        std::lock_guard<std::mutex> lock(tsv_mtx_);
        tsv_out_->write(tsv_chunk.data(), static_cast<std::streamsize>(tsv_chunk.size()));
    }

    // Only stream unambiguous, matched reads with a known sample to split files.
    if (!config_.split_dir.empty()) {
        for (const auto& item : items) {
            const auto& e = item.entry;
            if (e.status != CallStatus::unique || e.sample.empty()) continue;
            SplitSink* sink = get_or_open_split(e.sample);
            std::lock_guard<std::mutex> sample_lock(sink->mtx);
            write_fastq_gz(sink->gz, item.record);
        }
    }
}

bool ReportWriter::finalize() {
    {
        std::lock_guard<std::mutex> lock(tsv_mtx_);
        tsv_out_->flush();
    }

    std::lock_guard<std::mutex> map_lock(split_map_mtx_);
    for (auto& [name, sink] : split_files_) {
        if (!sink || !sink->gz) continue;
        gzflush(sink->gz, Z_FINISH);
        gzclose(sink->gz);
        sink->gz = nullptr;
    }
    split_files_.clear();
    if (config_.verbose) {
        std::cerr << "Reads with anchors: " << reads_with_anchors_.load() << "\n";
        std::cerr << "CBs matched:        " << reads_matched_.load()      << "\n";
        std::cerr << "Ambiguous calls:    " << ambiguous_calls_.load()    << "\n";
        if (!config_.split_dir.empty())
            std::cerr << "Split FASTQ.gz written to: " << config_.split_dir << "/\n";
    }
    return tsv_out_->good();
}
