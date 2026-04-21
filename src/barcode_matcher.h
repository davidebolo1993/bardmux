#pragma once

#include <string>
#include <vector>
#include <unordered_map>
#include <cstdint>

struct BarcodeEntry {
    std::string barcode;
    std::string sample;
    bool has_sample = false;

    void set_barcode(const std::string& s) {
        barcode = s;
    }
    void set_sample(const std::string& s) {
        sample = s;
        has_sample = true;
    }
};

struct BarcodeMatch {
    bool        found          = false;
    bool        ambiguous      = false;
    std::string matched_barcode;
    int         edit_distance  = -1;
    int         num_candidates = 0;
    std::string sample;
};

struct MatcherConfig {
    int max_edit_distance = 2;
    int min_margin        = 1;
    int kmer_length       = 8;
};

class BarcodeMatcher {
public:
    explicit BarcodeMatcher(const MatcherConfig& config = MatcherConfig());

    bool         load_whitelist(const std::string& filename);
    BarcodeMatch match_barcode(const std::string& barcode) const;
    size_t       whitelist_size() const { return whitelist_.size(); }

private:
    MatcherConfig config_;
    std::vector<BarcodeEntry> whitelist_;   // contiguous array — cache-friendly scan

    // Exact match: barcode string → whitelist index (O(1))
    std::unordered_map<std::string, uint32_t> exact_index_;

    // Flat CSR k-mer index: kmer_hash → [begin, end) into flat_candidates_
    std::unordered_map<uint64_t, std::pair<uint32_t,uint32_t>> kmer_offsets_;
    std::vector<uint32_t> flat_candidates_;   // sorted, deduplicated per bucket

    void build_kmer_index();

    int edit_distance_myers(const char* a, int na, const char* b, int nb) const;

    std::vector<uint32_t> kmer_candidates(const std::string& query) const;
};
