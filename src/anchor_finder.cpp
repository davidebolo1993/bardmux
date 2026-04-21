#include "anchor_finder.h"
#include <algorithm>
#include <cstdint>

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

int AnchorFinder::find_first(const std::string& text,
                              const std::string& pattern,
                              int max_errors,
                              int search_start,
                              int search_end) const
{
    const int tlen = static_cast<int>(text.size());
    const int plen = static_cast<int>(pattern.size());
    if (__builtin_expect(plen == 0 || plen > tlen, 0)) return -1;
    if (search_end < 0 || search_end > tlen - plen)
        search_end = tlen - plen;
    if (__builtin_expect(search_start > search_end, 0)) return -1;

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

        for (int j = 0; j < plen; ++j) {
            if (j == 0 || j == mid || j == plen - 1) continue;
            if (nt_upper(static_cast<unsigned char>(t[i + j])) != p[j]) {
                if (__builtin_expect(++mm > max_errors, 1)) goto next_pos;
            }
        }
        return search_start + i;
        next_pos:;
    }
    return -1;
}

AnchorMatch AnchorFinder::extract_fwd(const std::string& seq, int a1_pos, int a2_pos) const
{
    AnchorMatch m;
    m.found                = true;
    m.is_reverse_complement = false;
    m.anchor1_pos          = a1_pos;
    m.anchor2_pos          = a2_pos;

    const int cb_len  = config_.cb_length;
    const int umi_len = config_.umi_length;
    const int a1_end = a1_pos + static_cast<int>(config_.anchor1.size());

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

AnchorMatch AnchorFinder::extract_rc(const std::string& seq, int a2rc_pos, int a1rc_pos) const
{
    AnchorMatch m;
    m.found                = true;
    m.is_reverse_complement = true;
    m.anchor1_pos          = a1rc_pos;
    m.anchor2_pos          = a2rc_pos;

    const int cb_len  = config_.cb_length;
    const int umi_len = config_.umi_length;
    const int a2rc_end   = a2rc_pos + static_cast<int>(anchor2_rc_.size());
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
    const int expected = config_.cb_length + config_.umi_length;
    const int max_region = expected + (config_.use_gap_slack ? std::max(0, config_.gap_slack) : 0);

    int pos = 0;
    while (true) {
        int a1_pos = find_first(seq, config_.anchor1, config_.max_errors, pos, slen - a1len);
        if (__builtin_expect(a1_pos < 0, 1)) break;
        int a1_end = a1_pos + a1len;
        int a2_start = a1_end + config_.cb_length;
        int a2_end = std::min(slen - a2len, a1_end + max_region);
        if (a2_start <= a2_end) {
            int a2_pos = find_first(seq, config_.anchor2, config_.max_errors, a2_start, a2_end);
            if (__builtin_expect(a2_pos >= 0, 0))
                return extract_fwd(seq, a1_pos, a2_pos);
        }
        pos = a1_pos + 1;
    }
    return AnchorMatch{};
}

AnchorMatch AnchorFinder::try_rc(const std::string& seq) const
{
    const int slen    = static_cast<int>(seq.size());
    const int a2rclen = static_cast<int>(anchor2_rc_.size());
    const int a1rclen = static_cast<int>(anchor1_rc_.size());
    const int expected = config_.cb_length + config_.umi_length;
    const int max_region = expected + (config_.use_gap_slack ? std::max(0, config_.gap_slack) : 0);

    int pos = 0;
    while (true) {
        int a2rc_pos = find_first(seq, anchor2_rc_, config_.max_errors, pos, slen - a2rclen);
        if (__builtin_expect(a2rc_pos < 0, 1)) break;
        int a2rc_end = a2rc_pos + a2rclen;
        int a1_start = a2rc_end + config_.cb_length;
        int a1_end = std::min(slen - a1rclen, a2rc_end + max_region);
        if (a1_start <= a1_end) {
            int a1rc_pos = find_first(seq, anchor1_rc_, config_.max_errors, a1_start, a1_end);
            if (__builtin_expect(a1rc_pos >= 0, 0))
                return extract_rc(seq, a2rc_pos, a1rc_pos);
        }
        pos = a2rc_pos + 1;
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
