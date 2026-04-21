#pragma once

#include <string>
#include <fstream>
#include <mutex>
#include <atomic>
#include <unordered_map>
#include <vector>
#include <zlib.h>
#include "fastq_reader.h"

// Possible per-read call status
enum class CallStatus {
    no_anchor,   // no valid anchor pair found
    no_cb,       // anchors found but CB region too short
    no_match,    // CB extracted, no whitelist match within max_edit_distance
    ambiguous,   // ≥2 whitelist entries at same best edit distance
    unique       // single unambiguous match
};

inline const char* status_str(CallStatus s) {
    switch (s) {
        case CallStatus::no_anchor:  return "no_anchor";
        case CallStatus::no_cb:      return "no_cb";
        case CallStatus::no_match:   return "no_match";
        case CallStatus::ambiguous:  return "ambiguous";
        case CallStatus::unique:     return "unique";
    }
    return "unknown";
}

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
    int         n_candidates   = 0;    // number of whitelist entries at best distance
    std::string sample;
    bool        ambiguous      = false;
    CallStatus  status         = CallStatus::no_anchor;
};

struct ReportChunkItem {
    ReportEntry entry;
    FastqRecord record;
};

struct ReporterConfig {
    std::string output_file;
    std::string split_dir;
    bool        verbose = false;
};

class ReportWriter {
public:
    explicit ReportWriter(const ReporterConfig& config = ReporterConfig());
    ~ReportWriter();

    void add_entries(std::vector<ReportChunkItem>&& items);
    bool finalize();

    int reads_with_anchors() const { return reads_with_anchors_.load(); }
    int reads_matched()      const { return reads_matched_.load(); }
    int ambiguous_calls()    const { return ambiguous_calls_.load(); }

private:
    ReporterConfig config_;
    std::ostream*  tsv_out_  = nullptr;
    std::ofstream  tsv_fout_;
    char           tsv_buf_[1 << 20]{};

    std::unordered_map<std::string, gzFile> split_files_;
    std::mutex mtx_;

    std::atomic<int> reads_with_anchors_{0};
    std::atomic<int> reads_matched_     {0};
    std::atomic<int> ambiguous_calls_   {0};

    bool   open_tsv_output();
    void   write_tsv_header();
    static void append_tsv_line(std::string& out, const ReportEntry& e);
    gzFile get_or_open_split(const std::string& sample);
    void   write_fastq_gz(gzFile gz, const FastqRecord& rec);
    bool   ensure_dir(const std::string& path);
};
