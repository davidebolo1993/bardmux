#include "report_writer.h"
#include <iostream>
#include <stdexcept>
#include <sys/stat.h>   // mkdir
#include <errno.h>
#include <cstring>

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

bool ReportWriter::ensure_dir(const std::string& path) {
    struct stat st{};
    if (stat(path.c_str(), &st) == 0) {
        return S_ISDIR(st.st_mode);   // already exists
    }
    // Try to create (non-recursive; one level)
    if (mkdir(path.c_str(), 0755) == 0) return true;
    // May have been created by another thread between stat and mkdir
    if (errno == EEXIST) return true;
    return false;
}

gzFile ReportWriter::get_or_open_split(const std::string& sample) {
    // Caller holds mtx_
    auto it = split_files_.find(sample);
    if (it != split_files_.end()) return it->second;

    std::string path = config_.split_dir + "/" + sample + ".fastq.gz";
    gzFile gz = gzopen(path.c_str(), "wb");
    if (!gz)
        throw std::runtime_error("Cannot open split file: " + path);
    gzbuffer(gz, 1 << 20);   // 1 MiB write buffer
    split_files_[sample] = gz;
    return gz;
}

void ReportWriter::write_fastq_gz(gzFile gz, const FastqRecord& rec) {
    // Write standard 4-line FASTQ
    // name may already contain the comment after a space — prepend '@'
    std::string header = "@";
    header += rec.name;
    header += '\n';
    gzwrite(gz, header.data(),         static_cast<unsigned>(header.size()));
    gzwrite(gz, rec.sequence.data(),   static_cast<unsigned>(rec.sequence.size()));
    gzwrite(gz, "\n+\n",              3);
    gzwrite(gz, rec.quality.data(),    static_cast<unsigned>(rec.quality.size()));
    gzwrite(gz, "\n",                 1);
}

// ---------------------------------------------------------------------------
// Construction / destruction
// ---------------------------------------------------------------------------

ReportWriter::ReportWriter(const ReporterConfig& cfg) : config_(cfg) {
    // Create split directory up front if requested
    if (!cfg.split_dir.empty()) {
        if (!ensure_dir(cfg.split_dir))
            throw std::runtime_error("Cannot create split directory: " + cfg.split_dir);
    }
    if (!open_tsv_output())
        throw std::runtime_error("Cannot open TSV output: " + cfg.output_file);
    write_tsv_header();
}

ReportWriter::~ReportWriter() {
    // Flush gzip split files (finalize() should have been called, but be safe)
    for (auto& [name, gz] : split_files_)
        if (gz) gzclose(gz);
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
    *tsv_out_ << "read_id\tanchors_found\tstrand\tanchor1_pos\tanchor2_pos\t"
                 "extracted_cb\textracted_umi\tcb_matched\tmatched_cb\t"
                 "edit_distance\tsample\tambiguous\n";
}

// ---------------------------------------------------------------------------
// Core: add_entry  (thread-safe)
// ---------------------------------------------------------------------------

void ReportWriter::add_entry(const ReportEntry& e, const FastqRecord& rec) {
    // Atomic counter updates — no lock needed
    if (e.anchors_found)  reads_with_anchors_++;
    if (e.cb_matched)     reads_matched_++;
    if (e.ambiguous)      ambiguous_calls_++;

    // Build TSV line locally (outside lock)
    std::string line;
    line.reserve(160);
    line += e.read_id;                                line += '\t';
    line += e.anchors_found ? "true" : "false";       line += '\t';
    line += e.strand;                                 line += '\t';
    line += std::to_string(e.anchor1_pos);            line += '\t';
    line += std::to_string(e.anchor2_pos);            line += '\t';
    line += e.extracted_cb;                           line += '\t';
    line += e.extracted_umi;                          line += '\t';
    line += e.cb_matched  ? "true" : "false";         line += '\t';
    line += e.matched_cb;                             line += '\t';
    line += std::to_string(e.edit_distance);          line += '\t';
    line += e.sample;                                 line += '\t';
    line += e.ambiguous   ? "true" : "false";         line += '\n';

    // Serialised section: TSV write + optional split FASTQ write
    std::lock_guard<std::mutex> lock(mtx_);

    tsv_out_->write(line.data(), static_cast<std::streamsize>(line.size()));

    // Stream to per-sample FASTQ.gz only when we have a confirmed, unambiguous match
    if (!config_.split_dir.empty() && e.cb_matched && !e.sample.empty() && !e.ambiguous) {
        gzFile gz = get_or_open_split(e.sample);
        write_fastq_gz(gz, rec);
    }
}

// ---------------------------------------------------------------------------
// Finalize
// ---------------------------------------------------------------------------

bool ReportWriter::finalize() {
    std::lock_guard<std::mutex> lock(mtx_);

    tsv_out_->flush();

    for (auto& [name, gz] : split_files_) {
        gzflush(gz, Z_FINISH);
        gzclose(gz);
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
