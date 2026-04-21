#pragma once

#include <string>
#include <vector>
#include <unordered_map>
#include <cstdint>

struct BarcodeEntry {
    std::string barcode;
    std::string sample;
    bool        has_sample = false;
};

struct BarcodeMatch {
    bool        found          = false;
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
    std::vector<BarcodeEntry> whitelist_;

    // Exact match: barcode string -> whitelist index (O(1) lookup, checked first)
    std::unordered_map<std::string, uint32_t> exact_index_;

    // Approximate match: k-mer hash -> sorted list of whitelist indices
    std::unordered_map<uint64_t, std::vector<uint32_t>> kmer_index_;

    void build_kmer_index();

    // Bitwise Myers bit-vector edit distance (O(n*m/64))
    int edit_distance_myers(const std::string& s1, const std::string& s2) const;

    // Return candidate whitelist indices sharing >= 1 k-mer with query
    std::vector<uint32_t> kmer_candidates(const std::string& query) const;
};
