#pragma once

#include <string>
#include "fastq_reader.h"

struct AnchorConfig {
    std::string anchor1        = "CTACACGACGCTCTTCCGATCT";   // TruSeq R1
    std::string anchor2        = "TTTCTTATATGGG";             // TSO
    int         max_errors     = 3;
    int         cb_length      = 16;
    int         umi_length     = 12;
};

struct AnchorMatch {
    bool        found                  = false;
    bool        is_reverse_complement  = false;
    int         anchor1_end            = -1;   // right boundary of left anchor in raw read
    int         anchor2_start          = -1;   // left boundary of right anchor in raw read
    std::string extracted_cb;
    std::string extracted_umi;
};

class AnchorFinder {
public:
    explicit AnchorFinder(const AnchorConfig& cfg = AnchorConfig());

    AnchorMatch find_anchors(const FastqRecord& record) const;

    // Public utility used for testing
    static std::string reverse_complement(const std::string& seq);

private:
    AnchorConfig config_;
    std::string  anchor1_rc_;
    std::string  anchor2_rc_;

    // Find first occurrence of pattern in text[search_start..search_end] with <= max_errors mismatches
    int find_first(const std::string& text, const std::string& pattern,
                   int max_errors, int search_start = 0, int search_end = -1) const;

    // Forward: A1_fwd precedes A2_fwd  →  extract CB+UMI directly
    AnchorMatch try_fwd(const std::string& seq) const;

    // RC: A2_rc precedes A1_rc in raw read  →  RC the region to recover CB+UMI
    AnchorMatch try_rc(const std::string& seq) const;

    AnchorMatch extract_fwd(const std::string& seq, int a1_end,   int a2_start) const;
    AnchorMatch extract_rc (const std::string& seq, int a2rc_end, int a1rc_start) const;
};
