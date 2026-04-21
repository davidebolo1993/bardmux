#include "barcode_matcher.h"
#include <fstream>
#include <algorithm>
#include <limits>
#include <cctype>
#include <vector>
#include <unordered_map>

// 2-bit DNA encoding LUT
static const uint8_t BASE2BIT[256] = {
    // A=0 C=1 G=2 T=3 (upper and lower); everything else=0
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,1,2,0,0,0,3,0,0,0,0,0,0,0,0, 0,0,0,0,0,3,0,0,0,0,0,0,0,0,0,0,  // @A..T
    0,0,1,2,0,0,0,3,0,0,0,0,0,0,0,0, 0,0,0,0,0,3,0,0,0,0,0,0,0,0,0,0,  // `a..t
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
};

static inline uint64_t kmer_hash(const char* s, int start, int k) {
    uint64_t h = 0;
    for (int i = 0; i < k; ++i)
        h = (h << 2) | BASE2BIT[static_cast<unsigned char>(s[start + i])];
    return h;
}

// Myers bit-vector edit distance — operates on raw char* to avoid string copies.
// For sequences ≤ 64 bp: O(n·m/64), zero heap allocation.
int BarcodeMatcher::edit_distance_myers(const char* a, int n,
                                         const char* b, int m) const
{
    if (n == 0) return m;
    if (m == 0) return n;

    if (n > 64 || m > 64) {
        // Banded Levenshtein fallback
        const int bw = config_.max_edit_distance + 1;
        std::vector<int> prev(m + 1), curr(m + 1);
        for (int j = 0; j <= m; ++j) prev[j] = j;
        for (int i = 1; i <= n; ++i) {
            curr[0] = i;
            int lo = std::max(1, i - bw), hi = std::min(m, i + bw);
            for (int j = lo; j <= hi; ++j) {
                int cost = (a[i-1] == b[j-1]) ? 0 : 1;
                curr[j] = std::min({prev[j]+1, curr[j-1]+1, prev[j-1]+cost});
            }
            std::swap(prev, curr);
        }
        return prev[m];
    }

    // Reduced Peq: only 4 DNA bases needed → 4-entry array indexed by 2-bit code
    uint64_t Peq[4] = {};
    for (int j = 0; j < m; ++j)
        Peq[BASE2BIT[static_cast<unsigned char>(b[j])]] |= (uint64_t(1) << j);

    const uint64_t mask = (m == 64) ? ~uint64_t(0) : (uint64_t(1) << m) - 1;
    uint64_t VP = mask, VN = 0;
    int score = m;

    for (int i = 0; i < n; ++i) {
        uint64_t eq = Peq[BASE2BIT[static_cast<unsigned char>(a[i])]];
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
        const char* bc = whitelist_[i].barcode;
        int len = static_cast<int>(std::strlen(bc));

        exact_index_[bc] = i;

        if (len < k) continue;
        for (int j = 0; j <= len - k; ++j)
            tmp[kmer_hash(bc, j, k)].push_back(i);
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

    // Count k-mer hits per candidate using the flat CSR index
    std::vector<uint8_t> hits(whitelist_.size(), 0);
    for (int j = 0; j <= qlen - k; ++j) {
        auto it = kmer_offsets_.find(kmer_hash(query.c_str(), j, k));
        if (it != kmer_offsets_.end()) {
            auto [begin, end] = it->second;
            for (uint32_t x = begin; x < end; ++x) {
                uint32_t idx = flat_candidates_[x];
                if (hits[idx] < 255) hits[idx]++;
            }
        }
    }

    std::vector<uint32_t> cands;
    cands.reserve(256);
    for (uint32_t i = 0; i < static_cast<uint32_t>(whitelist_.size()); ++i)
        if (hits[i]) cands.push_back(i);
    return cands;
}

BarcodeMatch BarcodeMatcher::match_barcode(const std::string& barcode) const {
    BarcodeMatch result;
    if (barcode.empty() || whitelist_.empty()) return result;

    // Tier 1: exact match O(1)
    {
        auto it = exact_index_.find(barcode);
        if (it != exact_index_.end()) {
            uint32_t idx         = it->second;
            result.found          = true;
            result.matched_barcode = whitelist_[idx].barcode;
            result.edit_distance  = 0;
            result.num_candidates = 1;
            if (whitelist_[idx].has_sample)
                result.sample = whitelist_[idx].sample;
            return result;
        }
    }

    // Tier 2: k-mer filter
    auto candidates = kmer_candidates(barcode);
    if (candidates.empty()) return result;

    // Tier 3: Myers edit distance on raw char* — no string copies
    int best_dist   = std::numeric_limits<int>::max();
    int second_dist = std::numeric_limits<int>::max();
    std::vector<uint32_t> best_idx;
    const int na = static_cast<int>(barcode.size());

    for (uint32_t idx : candidates) {
        const char* bc = whitelist_[idx].barcode;
        int nb = static_cast<int>(std::strlen(bc));
        int d  = edit_distance_myers(barcode.c_str(), na, bc, nb);

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
