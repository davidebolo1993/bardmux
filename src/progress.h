#pragma once
#include <atomic>
#include <chrono>
#include <cstdio>
#include <unistd.h>

// Each counter is padded to a full cache line (64 bytes) to eliminate
// false sharing between worker threads.
struct alignas(64) PaddedAtomic {
    std::atomic<long> value{0};
    char pad[64 - sizeof(std::atomic<long>)]{};
};

struct ProgressDelta {
    long total       = 0;
    long with_anchor = 0;
    long cb_matched  = 0;
    long exact       = 0;
    long ambiguous   = 0;
    long no_donor    = 0;
};

class ProgressTracker {
public:
    explicit ProgressTracker(long interval = 100000)
        : interval_(interval),
          t_start_(std::chrono::steady_clock::now()),
          tty_(isatty(fileno(stderr)))
    {}

    void record_read(bool anchor_found, bool cb_matched, bool exact_match,
                     bool ambiguous, bool no_donor)
    {
        ProgressDelta d;
        d.total       = 1;
        d.with_anchor = anchor_found ? 1 : 0;
        d.cb_matched  = cb_matched ? 1 : 0;
        d.exact       = exact_match ? 1 : 0;
        d.ambiguous   = ambiguous ? 1 : 0;
        d.no_donor    = no_donor ? 1 : 0;
        record_batch(d);
    }

    void record_batch(const ProgressDelta& d) {
        if (d.total <= 0) return;
        long cur = total_.value.fetch_add(d.total, std::memory_order_relaxed) + d.total;
        if (d.with_anchor)
            with_anchor_.value.fetch_add(d.with_anchor, std::memory_order_relaxed);
        if (d.cb_matched)
            cb_matched_.value.fetch_add(d.cb_matched, std::memory_order_relaxed);
        if (d.exact)
            exact_.value.fetch_add(d.exact, std::memory_order_relaxed);
        if (d.ambiguous)
            ambiguous_.value.fetch_add(d.ambiguous, std::memory_order_relaxed);
        if (d.no_donor)
            no_donor_.value.fetch_add(d.no_donor, std::memory_order_relaxed);

        if (interval_ <= 0) return;
        long milestone = cur / interval_;
        long prev = last_milestone_.value.load(std::memory_order_relaxed);
        if (milestone > prev)
            if (last_milestone_.value.compare_exchange_strong(prev, milestone,
                    std::memory_order_relaxed, std::memory_order_relaxed))
                print_progress(cur, false);
    }

    void finalize() { print_progress(total_.value.load(), true); }

private:
    long  interval_;
    bool  tty_;
    std::chrono::steady_clock::time_point t_start_;

    PaddedAtomic total_;
    PaddedAtomic with_anchor_;
    PaddedAtomic cb_matched_;
    PaddedAtomic exact_;
    PaddedAtomic ambiguous_;
    PaddedAtomic no_donor_;
    PaddedAtomic last_milestone_;

    void print_progress(long n, bool final) const {
        if (n == 0) return;
        auto now    = std::chrono::steady_clock::now();
        double secs = std::chrono::duration<double>(now - t_start_).count();
        double rate = n / (secs > 0.001 ? secs : 0.001);

        long anchor  = with_anchor_.value.load();
        long matched = cb_matched_.value.load();
        long ex      = exact_.value.load();
        long approx  = matched - ex;
        long ambig   = ambiguous_.value.load();
        long nodon   = no_donor_.value.load();
        long nomatch = anchor - matched - ambig;

        auto pct = [&](long x) -> double {
            return n > 0 ? 100.0 * x / n : 0.0;
        };

        if (final) {
            std::fprintf(stderr,
                "%s"
                "\n"
                "┌──────────────────────── bardmux summary ───────────────────────────┐\n"
                "│  Reads processed   : %10ld   (%.0f reads/sec,  %.1f s total)\n"
                "│  Anchors found     : %10ld   (%5.1f%%)\n"
                "│  CB matched        : %10ld   (%5.1f%%)\n"
                "│    ├─ exact  (d=0) : %10ld   (%5.1f%%)\n"
                "│    └─ approx (d>0) : %10ld   (%5.1f%%)\n"
                "│  Ambiguous calls   : %10ld   (%5.1f%%)\n"
                "│  No donor label    : %10ld   (%5.1f%%)\n"
                "│  No CB match       : %10ld   (%5.1f%%)\n"
                "└─────────────────────────────────────────────────────────────────────┘\n",
                tty_ ? "\r\033[2K" : "",
                n, rate, secs,
                anchor,  pct(anchor),
                matched, pct(matched),
                ex,      pct(ex),
                approx,  pct(approx),
                ambig,   pct(ambig),
                nodon,   pct(nodon),
                nomatch, pct(nomatch));
        } else {
            char buf[256];
            std::snprintf(buf, sizeof(buf),
                "[bardmux] %9ld reads | anchors %5.1f%% | matched %5.1f%% "
                "(d=0: %5.1f%% d>0: %5.1f%%) | ambig %4.1f%% | %7.0f rd/s | %.0fs",
                n, pct(anchor), pct(matched), pct(ex), pct(approx),
                pct(ambig), rate, secs);
            if (tty_)
                std::fprintf(stderr, "\r%-140s", buf);
            else
                std::fprintf(stderr, "%s\n", buf);
        }
        std::fflush(stderr);
    }
};
