#pragma once
#include <atomic>
#include <chrono>
#include <cstdio>
#include <unistd.h>   // isatty

class ProgressTracker {
public:
    explicit ProgressTracker(long interval = 100000)
        : interval_(interval),
          t_start_(std::chrono::steady_clock::now()),
          tty_(isatty(fileno(stderr)))
    {}

    // Called by workers (thread-safe, lock-free).
    void record_read(bool anchor_found,
                     bool cb_matched,
                     bool exact_match,
                     bool ambiguous,
                     bool no_donor)
    {
        long cur = ++total_;
        if (anchor_found) ++with_anchor_;
        if (cb_matched)   ++cb_matched_;
        if (exact_match)  ++exact_;
        if (ambiguous)    ++ambiguous_;
        if (no_donor)     ++no_donor_;

        // Fire whenever we cross the next multiple of interval_.
        // Use a threshold: if cur / interval_ > last_milestone, try to claim it.
        long milestone = cur / interval_;
        long prev = last_milestone_.load(std::memory_order_relaxed);
        if (milestone > prev) {
            // One thread wins the CAS and prints; others skip
            if (last_milestone_.compare_exchange_strong(prev, milestone,
                    std::memory_order_relaxed, std::memory_order_relaxed))
                print_progress(cur, false);
        }
    }

    void finalize() { print_progress(total_.load(), true); }

private:
    long  interval_;
    bool  tty_;
    std::chrono::steady_clock::time_point t_start_;

    std::atomic<long> total_       {0};
    std::atomic<long> with_anchor_ {0};
    std::atomic<long> cb_matched_  {0};
    std::atomic<long> exact_       {0};
    std::atomic<long> ambiguous_   {0};
    std::atomic<long> no_donor_    {0};
    std::atomic<long> last_milestone_{0};

    void print_progress(long n, bool final) const {
        if (n == 0) return;

        auto now    = std::chrono::steady_clock::now();
        double secs = std::chrono::duration<double>(now - t_start_).count();
        double rate = n / (secs > 0.001 ? secs : 0.001);

        long anchor  = with_anchor_.load();
        long matched = cb_matched_.load();
        long ex      = exact_.load();
        long approx  = matched - ex;
        long ambig   = ambiguous_.load();
        long nodon   = no_donor_.load();
        long nomatch = anchor - matched - ambig;

        auto pct = [&](long x) -> double {
            return n > 0 ? 100.0 * x / n : 0.0;
        };

        if (final) {
            std::fprintf(stderr,
                "%s"   // clear the in-place line if on TTY
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
                tty_ ? "\r\033[2K" : "",   // erase current line on TTY
                n,    rate, secs,
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
                std::fprintf(stderr, "\r%-140s", buf);  // overwrite in place
            else
                std::fprintf(stderr, "%s\n", buf);
        }
        std::fflush(stderr);
    }
};
