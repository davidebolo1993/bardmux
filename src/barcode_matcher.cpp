#include "barcode_matcher.h"
#include <fstream>
#include <algorithm>
#include <limits>
#include <cctype>
#include <cstdlib>
#include <vector>
#include <unordered_map>

namespace {
inline char nt_upper(unsigned char c) {
    if (c >= 'a' && c <= 'z') c = static_cast<unsigned char>(c - ('a' - 'A'));
    return static_cast<char>(c);
}

inline bool encode_acgt(unsigned char c, uint8_t& out) {
    switch (nt_upper(c)) {
        case 'A': out = 0; return true;
        case 'C': out = 1; return true;
        case 'G': out = 2; return true;
        case 'T': out = 3; return true;
        default: return false;
    }
}

inline void normalize_barcode(std::string& s) {
    for (char& c : s) {
        c = nt_upper(static_cast<unsigned char>(c));
        if (c != 'A' && c != 'C' && c != 'G' && c != 'T' && c != 'N') c = 'N';
    }
}

inline bool kmer_hash_checked(const char* s, int start, int k, uint64_t& h) {
    h = 0;
    for (int i = 0; i < k; ++i) {
        uint8_t b = 0;
        if (!encode_acgt(static_cast<unsigned char>(s[start + i]), b)) return false;
        h = (h << 2) | b;
    }
    return true;
}
} // namespace

// Myers bit-vector edit distance — operates on raw char* to avoid string copies.
// For sequences ≤ 64 bp: O(n·m/64), zero heap allocation.
int BarcodeMatcher::edit_distance_myers(const char* a, int n,
                                         const char* b, int m) const
{
    if (n == 0) return m;
    if (m == 0) return n;
    if (std::abs(n - m) > config_.max_edit_distance)
        return config_.max_edit_distance + 1;

    if (n > 64 || m > 64) {
        // Banded Levenshtein fallback for long strings.
        const int bw  = config_.max_edit_distance;
        const int INF = config_.max_edit_distance + 1;
        std::vector<int> prev(m + 1), curr(m + 1, INF);
        for (int j = 0; j <= m; ++j) prev[j] = j;
        for (int i = 1; i <= n; ++i) {
            std::fill(curr.begin(), curr.end(), INF);
            int lo = std::max(1, i - bw);
            int hi = std::min(m, i + bw);
            if (lo == 1) curr[0] = i;

            int row_best = INF;
            for (int j = lo; j <= hi; ++j) {
                int cost = (a[i-1] == b[j-1]) ? 0 : 1;
                int del  = prev[j] + 1;
                int ins  = curr[j-1] + 1;
                int sub  = prev[j-1] + cost;
                curr[j]  = std::min({del, ins, sub});
                row_best = std::min(row_best, curr[j]);
            }
            std::swap(prev, curr);
            if (row_best > config_.max_edit_distance)
                return config_.max_edit_distance + 1;
        }
        return prev[m];
    }

    // Full ASCII Peq table preserves exact character semantics (e.g. N != A).
    uint64_t Peq[256] = {};
    for (int j = 0; j < m; ++j)
        Peq[static_cast<unsigned char>(b[j])] |= (uint64_t(1) << j);

    const uint64_t mask = (m == 64) ? ~uint64_t(0) : (uint64_t(1) << m) - 1;
    uint64_t VP = mask, VN = 0;
    int score = m;

    for (int i = 0; i < n; ++i) {
        uint64_t eq = Peq[static_cast<unsigned char>(a[i])];
        uint64_t X  = eq | VN;
        uint64_t D0 = ((VP + (X & VP)) ^ VP) | X;
        uint64_t HN = VP & D0;
        uint64_t HP = VN | ~(VP | D0);
        X  = (HP << 1) | 1;
        VN = X & D0;
        VP = (HN << 1) | ~(X | D0) & mask;
        if (HP & (uint64_t(1) << (m-1))) ++score;
        else if (HN & (uint64_t(1) << (m-1))) --score;
    }
    return score;
}

// Build flat CSR k-mer index
void BarcodeMatcher::build_kmer_index() {
    kmer_offsets_.clear();
    flat_candidates_.clear();
    exact_index_.clear();

    const int k = config_.kmer_length;

    // Pass 1: collect all (hash → index) pairs
    std::unordered_map<uint64_t, std::vector<uint32_t>> tmp;
    tmp.reserve(whitelist_.size() * 8);

    for (uint32_t i = 0; i < static_cast<uint32_t>(whitelist_.size()); ++i) {
        const std::string& bc = whitelist_[i].barcode;
        int len = static_cast<int>(bc.size());
        exact_index_.try_emplace(bc, i); // keep first occurrence on duplicates

        if (len < k) continue;
        for (int j = 0; j <= len - k; ++j) {
            uint64_t h = 0;
            if (kmer_hash_checked(bc.c_str(), j, k, h))
                tmp[h].push_back(i);
        }
    }

    // Pass 2: sort+dedup each bucket, lay out in flat array (CSR)
    flat_candidates_.reserve(tmp.size() * 4);
    kmer_offsets_.reserve(tmp.size());

    for (auto& [h, vec] : tmp) {
        std::sort(vec.begin(), vec.end());
        vec.erase(std::unique(vec.begin(), vec.end()), vec.end());

        uint32_t begin = static_cast<uint32_t>(flat_candidates_.size());
        flat_candidates_.insert(flat_candidates_.end(), vec.begin(), vec.end());
        uint32_t end   = static_cast<uint32_t>(flat_candidates_.size());
        kmer_offsets_[h] = {begin, end};
    }
}

BarcodeMatcher::BarcodeMatcher(const MatcherConfig& config) : config_(config) {}

bool BarcodeMatcher::load_whitelist(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) return false;

    whitelist_.clear();
    std::string line;
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        if (!line.empty() && line.back() == '\r') line.pop_back();

        auto tab = line.find('\t');
        std::string barcode = (tab != std::string::npos) ? line.substr(0, tab) : line;

        // Strip 10X partition suffix -[digit]
        if (barcode.size() >= 2) {
            size_t last = barcode.size() - 1;
            if (barcode[last-1] == '-' && std::isdigit((unsigned char)barcode[last]))
                barcode.resize(last - 1);
        }
        normalize_barcode(barcode);
        if (barcode.empty()) continue;

        BarcodeEntry e;
        e.set_barcode(barcode);
        if (tab != std::string::npos) {
            auto tab2 = line.find('\t', tab + 1);
            std::string samp = line.substr(tab + 1,
                tab2 == std::string::npos ? std::string::npos : tab2 - tab - 1);
            if (!samp.empty()) e.set_sample(samp);
        }
        whitelist_.push_back(e);
    }
    build_kmer_index();
    return true;
}

std::vector<uint32_t> BarcodeMatcher::kmer_candidates(const std::string& query) const {
    const int k    = config_.kmer_length;
    const int qlen = static_cast<int>(query.size());

    if (qlen < k || kmer_offsets_.empty()) {
        std::vector<uint32_t> all(whitelist_.size());
        for (uint32_t i = 0; i < static_cast<uint32_t>(whitelist_.size()); ++i) all[i] = i;
        return all;
    }

    // Deduplicate candidates with thread-local epoch marks (no O(N) scan).
    thread_local std::vector<uint32_t> seen_epoch;
    thread_local uint32_t epoch = 1;

    if (seen_epoch.size() < whitelist_.size()) {
        seen_epoch.assign(whitelist_.size(), 0);
        epoch = 1;
    }
    ++epoch;
    if (epoch == 0) {
        std::fill(seen_epoch.begin(), seen_epoch.end(), 0);
        epoch = 1;
    }

    std::vector<uint32_t> cands;
    cands.reserve(256);
    for (int j = 0; j <= qlen - k; ++j) {
        uint64_t h = 0;
        if (!kmer_hash_checked(query.c_str(), j, k, h)) continue;

        auto it = kmer_offsets_.find(h);
        if (it != kmer_offsets_.end()) {
            auto [begin, end] = it->second;
            for (uint32_t x = begin; x < end; ++x) {
                uint32_t idx = flat_candidates_[x];
                if (seen_epoch[idx] != epoch) {
                    seen_epoch[idx] = epoch;
                    cands.push_back(idx);
                }
            }
        }
    }
    return cands;
}

BarcodeMatch BarcodeMatcher::match_barcode(const std::string& barcode) const {
    BarcodeMatch result;
    if (barcode.empty() || whitelist_.empty()) return result;

    std::string query = barcode;
    normalize_barcode(query);

    // Tier 1: exact match O(1)
    {
        auto it = exact_index_.find(query);
        if (it != exact_index_.end()) {
            uint32_t idx         = it->second;
            result.found          = true;
            result.matched_barcode = whitelist_[idx].barcode;
            result.edit_distance  = 0;
            result.num_candidates = 1;
            if (whitelist_[idx].has_sample)
                result.sample = whitelist_[idx].sample;

            // For min_margin > 1, verify there isn't another near-tie candidate.
            if (config_.min_margin > 1) {
                int best_other = std::numeric_limits<int>::max();
                auto candidates = kmer_candidates(query);
                for (uint32_t cidx : candidates) {
                    if (cidx == idx) continue;
                    const std::string& bc = whitelist_[cidx].barcode;
                    int d = edit_distance_myers(query.c_str(),
                                                static_cast<int>(query.size()),
                                                bc.c_str(),
                                                static_cast<int>(bc.size()));
                    if (d < best_other) best_other = d;
                    if (best_other < config_.min_margin) break;
                }
                if (best_other < config_.min_margin)
                    result.ambiguous = true;
            }
            return result;
        }
    }

    // Tier 2: k-mer filter
    auto candidates = kmer_candidates(query);
    if (candidates.empty()) return result;

    // Tier 3: Myers edit distance on raw char* — no string copies
    int best_dist   = std::numeric_limits<int>::max();
    int second_dist = std::numeric_limits<int>::max();
    std::vector<uint32_t> best_idx;
    const int na = static_cast<int>(query.size());

    for (uint32_t idx : candidates) {
        const std::string& bc = whitelist_[idx].barcode;
        int nb = static_cast<int>(bc.size());
        int d  = edit_distance_myers(query.c_str(), na, bc.c_str(), nb);

        if (d < best_dist) {
            second_dist = best_dist;
            best_dist   = d;
            best_idx.clear();
            best_idx.push_back(idx);
        } else if (d == best_dist) {
            best_idx.push_back(idx);
        } else if (d < second_dist) {
            second_dist = d;
        }
    }

    if (best_dist > config_.max_edit_distance) return result;

    result.found           = true;
    result.matched_barcode = whitelist_[best_idx[0]].barcode;
    result.edit_distance   = best_dist;
    result.num_candidates  = static_cast<int>(best_idx.size());
    if (result.num_candidates > 1)
        result.ambiguous = true;
    else if (second_dist != std::numeric_limits<int>::max() &&
             (second_dist - best_dist) < config_.min_margin)
        result.ambiguous = true;

    if (whitelist_[best_idx[0]].has_sample)
        result.sample = whitelist_[best_idx[0]].sample;
    return result;
}
