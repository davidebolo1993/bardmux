#include "anchor_finder.h"
#include <algorithm>

// Complement lookup table
static const char comp_lut[256] = {
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','T','N','G','N','N','N','C','N','N','N','N','N','N','N','N',
    'N','N','N','N','A','N','N','N','N','N','N','N','N','N','N','N',
    'N','t','N','g','N','N','N','c','N','N','N','N','N','N','N','N',
    'N','N','N','N','a','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
};

std::string AnchorFinder::reverse_complement(const std::string& seq) {
    std::string rc(seq.size(), 'N');
    const auto n = seq.size();
    for (std::size_t i = 0; i < n; ++i)
        rc[n - 1 - i] = comp_lut[static_cast<unsigned char>(seq[i])];
    return rc;
}

AnchorFinder::AnchorFinder(const AnchorConfig& cfg)
    : config_(cfg),
      anchor1_rc_(reverse_complement(cfg.anchor1)),
      anchor2_rc_(reverse_complement(cfg.anchor2))
{}

int AnchorFinder::find_first(const std::string& text,
                              const std::string& pattern,
                              int max_errors,
                              int search_start,
                              int search_end) const
{
    const int tlen = static_cast<int>(text.size());
    const int plen = static_cast<int>(pattern.size());
    if (plen == 0 || plen > tlen) return -1;
    if (search_end < 0 || search_end > tlen - plen)
        search_end = tlen - plen;
    if (search_start > search_end) return -1;

    const char* t = text.data() + search_start;
    const char* p = pattern.data();
    const int span = search_end - search_start;

    for (int i = 0; i <= span; ++i) {
        int mm = 0;
        for (int j = 0; j < plen; ++j) {
            if (t[i + j] != p[j]) {
                if (++mm > max_errors) goto next_pos;
            }
        }
        return search_start + i;
        next_pos:;
    }
    return -1;
}

// Extract CB+UMI from a FORWARD-orientation anchor pair.
// Layout in raw read: ... A1 ... [CB (cb_length bp)][UMI (umi_length bp)] ... A2 ...
// a1_end:    position just after anchor1 ends
// a2_start:  position where anchor2 starts
AnchorMatch AnchorFinder::extract_fwd(const std::string& seq,
                                       int a1_end, int a2_start) const
{
    AnchorMatch m;
    m.found                = true;
    m.is_reverse_complement = false;
    m.anchor1_end          = a1_end;
    m.anchor2_start        = a2_start;

    const int cb_len  = config_.cb_length;
    const int umi_len = config_.umi_length;

    int avail = a2_start - a1_end;
    if (avail >= cb_len)
        m.extracted_cb = seq.substr(a1_end, cb_len);

    int umi_start = a1_end + cb_len;
    int umi_end   = std::min(umi_start + umi_len, a2_start);
    if (umi_end > umi_start)
        m.extracted_umi = seq.substr(umi_start, umi_end - umi_start);

    return m;
}

// Extract CB+UMI from a REVERSE-COMPLEMENT anchor pair.
//
// On the original sense strand:   5' A1 - CB - UMI - A2 3'
// On the raw (antisense) read:    5' A2_rc - UMI_rc - CB_rc - A1_rc 3'
//
// So in the raw read A2_rc appears BEFORE A1_rc.
// We search for A2_rc first (left anchor), then A1_rc (right anchor).
// The region between them, when reverse-complemented, yields: CB (first cb_length bp) + UMI.
//
// a2rc_end:   position just after A2_rc ends  (left boundary of insert)
// a1rc_start: position where A1_rc starts     (right boundary of insert)
AnchorMatch AnchorFinder::extract_rc(const std::string& seq,
                                      int a2rc_end, int a1rc_start) const
{
    AnchorMatch m;
    m.found                = true;
    m.is_reverse_complement = true;
    m.anchor1_end          = a2rc_end;    // reported as the "left" boundary
    m.anchor2_start        = a1rc_start;  // reported as the "right" boundary

    const int cb_len  = config_.cb_length;
    const int umi_len = config_.umi_length;

    int region_len = a1rc_start - a2rc_end;
    if (region_len < cb_len) return m;  // truncated — cannot extract full CB

    // RC the region to restore sense-strand orientation
    std::string region_rc = reverse_complement(seq.substr(a2rc_end, region_len));

    m.extracted_cb  = region_rc.substr(0, cb_len);
    int umi_avail   = static_cast<int>(region_rc.size()) - cb_len;
    if (umi_avail > 0)
        m.extracted_umi = region_rc.substr(cb_len, std::min(umi_len, umi_avail));

    return m;
}

// Forward orientation search:
//   Scan for A1_fwd, then A2_fwd to the right of it.
AnchorMatch AnchorFinder::try_fwd(const std::string& seq) const
{
    const int slen  = static_cast<int>(seq.size());
    const int a1len = static_cast<int>(config_.anchor1.size());
    const int a2len = static_cast<int>(config_.anchor2.size());

    int pos = 0;
    while (true) {
        int a1_pos = find_first(seq, config_.anchor1, config_.max_errors, pos, slen - a1len);
        if (a1_pos < 0) break;
        int a1_end = a1_pos + a1len;
        int a2_pos = find_first(seq, config_.anchor2, config_.max_errors,
                                 a1_end + config_.cb_length, slen - a2len);
        if (a2_pos >= 0)
            return extract_fwd(seq, a1_end, a2_pos);
        pos = a1_pos + 1;
    }
    return AnchorMatch{};
}

// Reverse-complement orientation search:
//   On the raw read, A2_rc appears FIRST, A1_rc appears SECOND.
//   Scan for A2_rc, then A1_rc to the right of it.
AnchorMatch AnchorFinder::try_rc(const std::string& seq) const
{
    const int slen   = static_cast<int>(seq.size());
    const int a2rclen = static_cast<int>(anchor2_rc_.size());
    const int a1rclen = static_cast<int>(anchor1_rc_.size());

    int pos = 0;
    while (true) {
        int a2rc_pos = find_first(seq, anchor2_rc_, config_.max_errors, pos, slen - a2rclen);
        if (a2rc_pos < 0) break;
        int a2rc_end = a2rc_pos + a2rclen;
        // Minimum gap between the two anchors = cb_length (exact; UMI can be truncated)
        int a1rc_pos = find_first(seq, anchor1_rc_, config_.max_errors,
                                   a2rc_end + config_.cb_length, slen - a1rclen);
        if (a1rc_pos >= 0)
            return extract_rc(seq, a2rc_end, a1rc_pos);
        pos = a2rc_pos + 1;
    }
    return AnchorMatch{};
}

AnchorMatch AnchorFinder::find_anchors(const FastqRecord& record) const {
    // Try forward orientation first (most common for the sense strand)
    {
        auto m = try_fwd(record.sequence);
        if (m.found) return m;
    }
    // Try reverse-complement orientation
    {
        auto m = try_rc(record.sequence);
        if (m.found) return m;
    }
    return AnchorMatch{};
}
