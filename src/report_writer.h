#pragma once

#include <string>
#include <fstream>
#include <mutex>
#include <atomic>
#include <unordered_map>
#include <memory>
#include <zlib.h>
#include "fastq_reader.h"

struct ReportEntry {
    std::string read_id;
    bool        anchors_found  = false;
    std::string strand         = "FWD";
    int         anchor1_pos    = -1;
    int         anchor2_pos    = -1;
    std::string extracted_cb;
    std::string extracted_umi;
    bool        cb_matched     = false;
    std::string matched_cb;
    int         edit_distance  = -1;
    std::string sample;
    bool        ambiguous      = false;
};

struct ReporterConfig {
    std::string output_file;   // TSV output (empty = stdout)
    std::string split_dir;     // directory for per-sample FASTQ.gz (empty = disabled)
    bool        compress = true;  // gzip-compress split files (.fastq.gz)
    bool        verbose  = false;
};

class ReportWriter {
public:
    explicit ReportWriter(const ReporterConfig& config = ReporterConfig());
    ~ReportWriter();

    // Thread-safe: write TSV line + optionally stream read to per-sample FASTQ.gz
    void add_entry(const ReportEntry& entry, const FastqRecord& record);

    // Flush and close all output streams
    bool finalize();

    int reads_with_anchors() const { return reads_with_anchors_.load(); }
    int reads_matched()      const { return reads_matched_.load(); }
    int ambiguous_calls()    const { return ambiguous_calls_.load(); }

private:
    ReporterConfig config_;

    // TSV output
    std::ostream*  tsv_out_  = nullptr;
    std::ofstream  tsv_fout_;
    char           tsv_buf_[1 << 20]{};

    // Per-sample gzip split files: sample_name -> gzFile
    std::unordered_map<std::string, gzFile> split_files_;

    std::mutex mtx_;   // protects tsv_out_ writes AND split_files_ map

    std::atomic<int> reads_with_anchors_{0};
    std::atomic<int> reads_matched_     {0};
    std::atomic<int> ambiguous_calls_   {0};

    bool   open_tsv_output();
    void   write_tsv_header();
    gzFile get_or_open_split(const std::string& sample);
    void   write_fastq_gz(gzFile gz, const FastqRecord& rec);
    bool   ensure_dir(const std::string& path);
};
