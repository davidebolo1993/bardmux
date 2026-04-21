#include "barcode_matcher.h"
#include <fstream>
#include <algorithm>
#include <limits>
#include <cctype>

// 2-bit DNA encoding: A=0, C=1, G=2, T=3
static inline uint64_t base2bit(char c) {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default:            return 0;
    }
}

static uint64_t kmer_hash(const std::string& s, int start, int k) {
    uint64_t h = 0;
    for (int i = 0; i < k; ++i) h = (h << 2) | base2bit(s[start + i]);
    return h;
}

// Myers (1999) bit-vector edit distance.
// For barcodes <= 64 bp: O(n*m/64) time, O(m) space, zero heap allocation.
// Falls back to banded Levenshtein DP for longer sequences.
int BarcodeMatcher::edit_distance_myers(const std::string& a,
                                         const std::string& b) const
{
    const int n = static_cast<int>(a.size());
    const int m = static_cast<int>(b.size());
    if (n == 0) return m;
    if (m == 0) return n;

    if (n > 64 || m > 64) {
        // Banded DP fallback
        std::vector<int> prev(m + 1), curr(m + 1);
        for (int j = 0; j <= m; ++j) prev[j] = j;
        for (int i = 1; i <= n; ++i) {
            curr[0] = i;
            for (int j = 1; j <= m; ++j) {
                int cost = (a[i-1] == b[j-1]) ? 0 : 1;
                curr[j] = std::min({prev[j]+1, curr[j-1]+1, prev[j-1]+cost});
            }
            std::swap(prev, curr);
        }
        return prev[m];
    }

    // Build Peq: bitmask of pattern positions equal to each character
    uint64_t Peq[256] = {};
    for (int j = 0; j < m; ++j)
        Peq[static_cast<unsigned char>(b[j])] |= (uint64_t(1) << j);

    uint64_t X, D0, HN, HP, VN, VP;
    VP = (m == 64) ? ~uint64_t(0) : (uint64_t(1) << m) - 1;
    VN = 0;
    int score = m;

    for (int i = 0; i < n; ++i) {
        X  = Peq[static_cast<unsigned char>(a[i])] | VN;
        D0 = ((VP + (X & VP)) ^ VP) | X;
        HN = VP & D0;
        HP = VN | ~(VP | D0);
        X  = (HP << 1) | 1;
        VN = X & D0;
        VP = (HN << 1) | ~(X | D0);
        if (HP & (uint64_t(1) << (m-1))) score++;
        else if (HN & (uint64_t(1) << (m-1))) score--;
    }
    return score;
}

void BarcodeMatcher::build_kmer_index() {
    kmer_index_.clear();
    exact_index_.clear();
    const int k = config_.kmer_length;

    for (uint32_t i = 0; i < static_cast<uint32_t>(whitelist_.size()); ++i) {
        const std::string& bc = whitelist_[i].barcode;

        // Exact-match index: O(1) lookup before any edit-distance work
        exact_index_[bc] = i;

        const int len = static_cast<int>(bc.size());
        if (len < k) continue;
        for (int j = 0; j <= len - k; ++j)
            kmer_index_[kmer_hash(bc, j, k)].push_back(i);
    }
    for (auto& [key, vec] : kmer_index_) {
        std::sort(vec.begin(), vec.end());
        vec.erase(std::unique(vec.begin(), vec.end()), vec.end());
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
            auto last = barcode.size() - 1;
            if (barcode[last-1] == '-' && std::isdigit((unsigned char)barcode[last]))
                barcode.resize(last - 1);
        }
        if (barcode.empty()) continue;

        BarcodeEntry e;
        e.barcode    = std::move(barcode);
        e.has_sample = (tab != std::string::npos);
        if (e.has_sample) {
            auto tab2 = line.find('\t', tab + 1);
            e.sample = line.substr(tab + 1,
                tab2 == std::string::npos ? std::string::npos : tab2 - tab - 1);
        }
        whitelist_.push_back(std::move(e));
    }
    build_kmer_index();
    return true;
}

std::vector<uint32_t> BarcodeMatcher::kmer_candidates(const std::string& query) const {
    const int k    = config_.kmer_length;
    const int qlen = static_cast<int>(query.size());

    if (qlen < k || kmer_index_.empty()) {
        std::vector<uint32_t> all(whitelist_.size());
        for (uint32_t i = 0; i < static_cast<uint32_t>(whitelist_.size()); ++i) all[i] = i;
        return all;
    }

    std::vector<uint8_t> hits(whitelist_.size(), 0);
    for (int j = 0; j <= qlen - k; ++j) {
        auto it = kmer_index_.find(kmer_hash(query, j, k));
        if (it != kmer_index_.end())
            for (uint32_t idx : it->second)
                if (hits[idx] < 255) hits[idx]++;
    }

    std::vector<uint32_t> cands;
    cands.reserve(256);
    for (uint32_t i = 0; i < static_cast<uint32_t>(whitelist_.size()); ++i)
        if (hits[i] > 0) cands.push_back(i);
    return cands;
}

BarcodeMatch BarcodeMatcher::match_barcode(const std::string& barcode) const {
    BarcodeMatch result;
    if (barcode.empty() || whitelist_.empty()) return result;

    // --- Fast path: exact match (edit distance = 0) ---
    // No need to run Myers or even touch the k-mer index.
    {
        auto it = exact_index_.find(barcode);
        if (it != exact_index_.end()) {
            uint32_t idx       = it->second;
            result.found           = true;
            result.matched_barcode = whitelist_[idx].barcode;
            result.edit_distance   = 0;
            result.num_candidates  = 1;
            if (whitelist_[idx].has_sample)
                result.sample = whitelist_[idx].sample;
            return result;
        }
    }

    // --- Slow path: approximate matching via k-mer filter + Myers ---
    auto candidates = kmer_candidates(barcode);
    if (candidates.empty()) return result;

    int best_dist   = std::numeric_limits<int>::max();
    int second_dist = std::numeric_limits<int>::max();
    std::vector<uint32_t> best_idx;

    for (uint32_t idx : candidates) {
        int d = edit_distance_myers(barcode, whitelist_[idx].barcode);
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
    if (whitelist_[best_idx[0]].has_sample)
        result.sample = whitelist_[best_idx[0]].sample;
    return result;
}
