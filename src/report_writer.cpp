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

gzFile ReportWriter::get_or_open_split(const std::string& sample) {
    auto it = split_files_.find(sample);
    if (it != split_files_.end()) return it->second;
    std::string path = config_.split_dir + "/" + sample + ".fastq.gz";
    gzFile gz = gzopen(path.c_str(), "wb");
    if (!gz) throw std::runtime_error("Cannot open split file: " + path);
    gzbuffer(gz, 1 << 20);
    split_files_[sample] = gz;
    return gz;
}

void ReportWriter::write_fastq_gz(gzFile gz, const FastqRecord& rec) {
    // '@' + name + '\n'
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

void ReportWriter::add_entry(const ReportEntry& e, const FastqRecord& rec) {
    if (e.anchors_found)  reads_with_anchors_++;
    if (e.cb_matched)     reads_matched_++;
    if (e.ambiguous)      ambiguous_calls_++;

    // Build TSV into a stack buffer — zero heap allocation
    char buf[1024];
    int n = std::snprintf(buf, sizeof(buf),
        "%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\t%d\t%s\t%s\n",
        e.read_id.c_str(),
        e.anchors_found ? "true" : "false",
        e.strand.c_str(),
        e.anchor1_pos,
        e.anchor2_pos,
        e.extracted_cb.c_str(),
        e.extracted_umi.c_str(),
        e.cb_matched ? "true" : "false",
        e.matched_cb.c_str(),
        e.edit_distance,
        e.sample.c_str(),
        e.ambiguous ? "true" : "false");

    if (n < 0) return;  // snprintf error
    if (n >= static_cast<int>(sizeof(buf))) n = sizeof(buf) - 1;

    std::lock_guard<std::mutex> lock(mtx_);
    tsv_out_->write(buf, n);

    if (!config_.split_dir.empty() && e.cb_matched && !e.sample.empty() && !e.ambiguous) {
        gzFile gz = get_or_open_split(e.sample);
        write_fastq_gz(gz, rec);
    }
}

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
