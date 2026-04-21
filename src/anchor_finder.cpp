#include "anchor_finder.h"
#include <algorithm>
#include <cstdint>
#include <limits>
#include <vector>

namespace {
inline char nt_upper(unsigned char c) {
    if (c >= 'a' && c <= 'z') c = static_cast<unsigned char>(c - ('a' - 'A'));
    return static_cast<char>(c);
}

inline char complement_base(unsigned char c) {
    switch (nt_upper(c)) {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        default:  return 'N';
    }
}

inline void uppercase_inplace(std::string& s) {
    for (char& c : s) c = nt_upper(static_cast<unsigned char>(c));
}
} // namespace

std::string AnchorFinder::reverse_complement(const std::string& seq) {
    std::string rc(seq.size(), 'N');
    const auto n = seq.size();
    for (std::size_t i = 0; i < n; ++i)
        rc[n - 1 - i] = complement_base(static_cast<unsigned char>(seq[i]));
    return rc;
}

AnchorFinder::AnchorFinder(const AnchorConfig& cfg)
    : config_(cfg)
{
    uppercase_inplace(config_.anchor1);
    uppercase_inplace(config_.anchor2);
    anchor1_rc_ = reverse_complement(config_.anchor1);
    anchor2_rc_ = reverse_complement(config_.anchor2);
}

AnchorFinder::AnchorHit AnchorFinder::find_first(const std::string& text,
                                                 const std::string& pattern,
                                                 int search_start,
                                                 int search_end) const
{
    if (config_.max_edits >= 0)
        return find_first_edit(text, pattern, config_.max_edits, search_start, search_end);
    return find_first_hamming(text, pattern, config_.max_errors, search_start, search_end);
}

AnchorFinder::AnchorHit AnchorFinder::find_first_hamming(const std::string& text,
                                                         const std::string& pattern,
                                                         int max_errors,
                                                         int search_start,
                                                         int search_end) const
{
    AnchorHit miss;
    const int tlen = static_cast<int>(text.size());
    const int plen = static_cast<int>(pattern.size());
    if (__builtin_expect(plen == 0 || plen > tlen, 0)) return miss;

    if (search_end < 0 || search_end > tlen - plen)
        search_end = tlen - plen;
    if (__builtin_expect(search_start > search_end, 0)) return miss;

    const char* t = text.data() + search_start;
    const char* p = pattern.data();
    const int span = search_end - search_start;
    const int mid = plen >> 1;

    for (int i = 0; i <= span; ++i) {
        int mm = 0;
        // Fast sentinels: reject obvious non-matches before full scan.
        if (nt_upper(static_cast<unsigned char>(t[i])) != p[0] &&
            __builtin_expect(++mm > max_errors, 1)) {
            continue;
        }
        if (mid > 0 && mid < plen - 1 &&
            nt_upper(static_cast<unsigned char>(t[i + mid])) != p[mid] &&
            __builtin_expect(++mm > max_errors, 1)) {
            continue;
        }
        if (plen > 1 &&
            nt_upper(static_cast<unsigned char>(t[i + plen - 1])) != p[plen - 1] &&
            __builtin_expect(++mm > max_errors, 1)) {
            continue;
        }

        bool reject = false;
        for (int j = 0; j < plen; ++j) {
            if (j == 0 || j == mid || j == plen - 1) continue;
            if (nt_upper(static_cast<unsigned char>(t[i + j])) != p[j]) {
                if (__builtin_expect(++mm > max_errors, 1)) {
                    reject = true;
                    break;
                }
            }
        }
        if (reject) continue;
        AnchorHit hit;
        hit.start = search_start + i;
        hit.end   = hit.start + plen;
        hit.score = mm;
        return hit;
    }
    return miss;
}

AnchorFinder::AnchorHit AnchorFinder::find_first_edit(const std::string& text,
                                                      const std::string& pattern,
                                                      int max_edits,
                                                      int search_start,
                                                      int search_end) const
{
    AnchorHit miss;
    const int tlen = static_cast<int>(text.size());
    const int plen = static_cast<int>(pattern.size());
    if (plen == 0 || tlen == 0) return miss;

    const int min_span = std::max(1, plen - max_edits);
    const int max_span = plen + max_edits;
    if (search_end < 0 || search_end > tlen - min_span)
        search_end = tlen - min_span;
    if (search_start < 0) search_start = 0;
    if (search_start > search_end) return miss;

    const int window_end = std::min(tlen, search_end + max_span);
    const int n = window_end - search_start;
    if (n <= 0) return miss;

    const char* txt = text.data() + search_start;
    const char* pat = pattern.data();

    const int INF = max_edits + 1;
    std::vector<int> prev_d(plen + 1), curr_d(plen + 1);
    std::vector<int> prev_s(plen + 1), curr_s(plen + 1);
    for (int j = 0; j <= plen; ++j) {
        prev_d[j] = j;
        prev_s[j] = 0;
    }

    int best_start = std::numeric_limits<int>::max();
    int best_end   = -1;
    int best_score = INF;

    for (int i = 1; i <= n; ++i) {
        curr_d[0] = 0;
        curr_s[0] = i;
        int row_best = INF;

        const char tc = nt_upper(static_cast<unsigned char>(txt[i - 1]));
        for (int j = 1; j <= plen; ++j) {
            int d_del = prev_d[j] + 1;
            int s_del = prev_s[j];

            int d_ins = curr_d[j - 1] + 1;
            int s_ins = curr_s[j - 1];

            int cost  = (tc == pat[j - 1]) ? 0 : 1;
            int d_sub = prev_d[j - 1] + cost;
            int s_sub = prev_s[j - 1];

            int best_d = d_del;
            int best_s = s_del;
            if (d_ins < best_d || (d_ins == best_d && s_ins < best_s)) {
                best_d = d_ins;
                best_s = s_ins;
            }
            if (d_sub < best_d || (d_sub == best_d && s_sub < best_s)) {
                best_d = d_sub;
                best_s = s_sub;
            }

            curr_d[j] = best_d;
            curr_s[j] = best_s;
            if (best_d < row_best) row_best = best_d;
        }

        if (curr_d[plen] <= max_edits) {
            int global_start = search_start + curr_s[plen];
            int global_end   = search_start + i;
            int score        = curr_d[plen];

            if (global_start >= search_start && global_start <= search_end &&
                global_end > global_start) {
                if (global_start < best_start ||
                    (global_start == best_start && score < best_score) ||
                    (global_start == best_start && score == best_score && global_end < best_end)) {
                    best_start = global_start;
                    best_end   = global_end;
                    best_score = score;
                }
            }
        }

        // No need to continue when even the best cell in this row already
        // exceeds the threshold and we passed the latest feasible end for the
        // current best start.
        if (best_start != std::numeric_limits<int>::max() &&
            row_best > max_edits &&
            i > (best_start - search_start) + max_span) {
            break;
        }

        std::swap(prev_d, curr_d);
        std::swap(prev_s, curr_s);
    }

    if (best_start == std::numeric_limits<int>::max()) return miss;
    AnchorHit hit;
    hit.start = best_start;
    hit.end   = best_end;
    hit.score = best_score;
    return hit;
}

AnchorMatch AnchorFinder::extract_fwd(const std::string& seq, int a1_pos, int a1_end, int a2_pos) const
{
    AnchorMatch m;
    m.found                = true;
    m.is_reverse_complement = false;
    m.anchor1_pos          = a1_pos;
    m.anchor2_pos          = a2_pos;

    const int cb_len  = config_.cb_length;
    const int umi_len = config_.umi_length;

    if (a2_pos - a1_end >= cb_len)
        m.extracted_cb = seq.substr(a1_end, cb_len);

    int umi_start = a1_end + cb_len;
    int umi_end   = std::min(umi_start + umi_len, a2_pos);
    if (umi_end > umi_start)
        m.extracted_umi = seq.substr(umi_start, umi_end - umi_start);

    uppercase_inplace(m.extracted_cb);
    uppercase_inplace(m.extracted_umi);
    return m;
}

AnchorMatch AnchorFinder::extract_rc(const std::string& seq, int a2rc_pos, int a2rc_end, int a1rc_pos) const
{
    AnchorMatch m;
    m.found                = true;
    m.is_reverse_complement = true;
    m.anchor1_pos          = a1rc_pos;
    m.anchor2_pos          = a2rc_pos;

    const int cb_len  = config_.cb_length;
    const int umi_len = config_.umi_length;
    const int a1rc_start = a1rc_pos;
    const int region_len = a1rc_start - a2rc_end;

    if (__builtin_expect(region_len < cb_len, 0)) return m;

    // Use a stack buffer for the RC region (≤ 128 bp for CB+UMI)
    // Avoids heap allocation for the most common short-insert case.
    const int MAX_STACK = 128;
    if (region_len <= MAX_STACK) {
        char buf[MAX_STACK + 1];
        const char* src = seq.data() + a2rc_end;
        for (int i = 0; i < region_len; ++i)
            buf[i] = complement_base(static_cast<unsigned char>(src[region_len - 1 - i]));
        buf[region_len] = '\0';
        m.extracted_cb.assign(buf, cb_len);
        int umi_avail = region_len - cb_len;
        if (umi_avail > 0)
            m.extracted_umi.assign(buf + cb_len,
                                   std::min(umi_len, umi_avail));
    } else {
        // Fallback for unusually large regions
        std::string rc = seq.substr(a2rc_end, region_len);
        for (int i = 0, j = region_len - 1; i < j; ++i, --j)
            std::swap(rc[i], rc[j]);
        for (char& c : rc) c = complement_base(static_cast<unsigned char>(c));
        m.extracted_cb = rc.substr(0, cb_len);
        if (region_len > cb_len)
            m.extracted_umi = rc.substr(cb_len,
                std::min(umi_len, region_len - cb_len));
    }
    return m;
}

AnchorMatch AnchorFinder::try_fwd(const std::string& seq) const
{
    const int slen  = static_cast<int>(seq.size());
    const int a1len = static_cast<int>(config_.anchor1.size());
    const int a2len = static_cast<int>(config_.anchor2.size());
    const int ed    = std::max(0, config_.max_edits);
    const int min_a1_span = (config_.max_edits >= 0) ? std::max(1, a1len - ed) : a1len;
    const int min_a2_span = (config_.max_edits >= 0) ? std::max(1, a2len - ed) : a2len;
    const int expected = config_.cb_length + config_.umi_length;
    const int max_region = expected + (config_.use_gap_slack ? std::max(0, config_.gap_slack) : 0);

    int pos = 0;
    while (true) {
        AnchorHit a1 = find_first(seq, config_.anchor1, pos, slen - min_a1_span);
        if (__builtin_expect(!a1.found(), 1)) break;
        int a1_end = a1.end;
        int a2_start = a1_end + config_.cb_length;
        int a2_end = std::min(slen - min_a2_span, a1_end + max_region);
        if (a2_start <= a2_end) {
            AnchorHit a2 = find_first(seq, config_.anchor2, a2_start, a2_end);
            if (__builtin_expect(a2.found(), 0))
                return extract_fwd(seq, a1.start, a1_end, a2.start);
        }
        pos = a1.start + 1;
    }
    return AnchorMatch{};
}

AnchorMatch AnchorFinder::try_rc(const std::string& seq) const
{
    const int slen    = static_cast<int>(seq.size());
    const int a2rclen = static_cast<int>(anchor2_rc_.size());
    const int a1rclen = static_cast<int>(anchor1_rc_.size());
    const int ed      = std::max(0, config_.max_edits);
    const int min_a2_span = (config_.max_edits >= 0) ? std::max(1, a2rclen - ed) : a2rclen;
    const int min_a1_span = (config_.max_edits >= 0) ? std::max(1, a1rclen - ed) : a1rclen;
    const int expected = config_.cb_length + config_.umi_length;
    const int max_region = expected + (config_.use_gap_slack ? std::max(0, config_.gap_slack) : 0);

    int pos = 0;
    while (true) {
        AnchorHit a2 = find_first(seq, anchor2_rc_, pos, slen - min_a2_span);
        if (__builtin_expect(!a2.found(), 1)) break;
        int a2rc_end = a2.end;
        int a1_start = a2rc_end + config_.cb_length;
        int a1_end = std::min(slen - min_a1_span, a2rc_end + max_region);
        if (a1_start <= a1_end) {
            AnchorHit a1 = find_first(seq, anchor1_rc_, a1_start, a1_end);
            if (__builtin_expect(a1.found(), 0))
                return extract_rc(seq, a2.start, a2rc_end, a1.start);
        }
        pos = a2.start + 1;
    }
    return AnchorMatch{};
}

AnchorMatch AnchorFinder::find_anchors(const FastqRecord& record) const {
    {
        auto m = try_fwd(record.sequence);
        if (__builtin_expect(m.found, 0)) return m;
    }
    {
        auto m = try_rc(record.sequence);
        if (__builtin_expect(m.found, 0)) return m;
    }
    return AnchorMatch{};
}
