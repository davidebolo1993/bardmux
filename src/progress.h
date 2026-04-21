#pragma once
#include <atomic>
#include <chrono>
#include <cstdio>
#include <string>
#include <unistd.h>   // isatty

class ProgressTracker {
public:
    explicit ProgressTracker(int interval = 100000)
        : interval_(interval),
          t_start_(std::chrono::steady_clock::now()),
          tty_(isatty(fileno(stderr)))
    {}

    // Called by workers (thread-safe, lock-free)
    void record_read(bool anchor_found,
                     bool cb_matched,
                     bool exact_match,   // edit_distance == 0
                     bool ambiguous,
                     bool no_donor)      // anchor+CB found but no sample label
    {
        ++total_;
        if (anchor_found)  ++with_anchor_;
        if (cb_matched)    ++cb_matched_;
        if (exact_match)   ++exact_;
        if (ambiguous)     ++ambiguous_;
        if (no_donor)      ++no_donor_;

        // Print every `interval_` reads (only one thread wins the CAS)
        long cur = total_.load(std::memory_order_relaxed);
        if (cur % interval_ == 0) {
            long expected = cur;
            if (last_printed_.compare_exchange_strong(expected, cur,
                    std::memory_order_relaxed, std::memory_order_relaxed))
                print(cur, false);
        }
    }

    // Final summary — always printed
    void finalize() { print(total_.load(), true); }

private:
    int  interval_;
    bool tty_;
    std::chrono::steady_clock::time_point t_start_;

    std::atomic<long> total_      {0};
    std::atomic<long> with_anchor_{0};
    std::atomic<long> cb_matched_ {0};
    std::atomic<long> exact_      {0};
    std::atomic<long> ambiguous_  {0};
    std::atomic<long> no_donor_   {0};
    std::atomic<long> last_printed_{-1};

    void print(long n, bool final) const {
        if (n == 0) return;

        auto now     = std::chrono::steady_clock::now();
        double secs  = std::chrono::duration<double>(now - t_start_).count();
        double rate  = n / (secs > 0 ? secs : 1.0);

        long  anchor   = with_anchor_.load();
        long  matched  = cb_matched_.load();
        long  exact    = exact_.load();
        long  approx   = matched - exact;
        long  ambig    = ambiguous_.load();
        long  nodonor  = no_donor_.load();
        long  nomatch  = anchor - matched - ambig;

        // Percentages vs total reads
        auto pct = [&](long x) -> double {
            return n > 0 ? 100.0 * x / n : 0.0;
        };

        char buf[512];
        if (final) {
            std::snprintf(buf, sizeof(buf),
                "\n"
                "┌─────────────────────── bardmux summary ───────────────────────────┐\n"
                "│ Reads processed   : %10ld   (%.0f reads/sec, %.1f s total)  \n"
                "│ Anchors found     : %10ld   (%5.1f%%)                        \n"
                "│ CB matched        : %10ld   (%5.1f%%)                        \n"
                "│   ├─ exact (d=0)  : %10ld   (%5.1f%%)                        \n"
                "│   └─ approx (d>0) : %10ld   (%5.1f%%)                        \n"
                "│ Ambiguous calls   : %10ld   (%5.1f%%)                        \n"
                "│ No donor label    : %10ld   (%5.1f%%)                        \n"
                "│ No CB match       : %10ld   (%5.1f%%)                        \n"
                "└────────────────────────────────────────────────────────────────────┘\n",
                n, rate, secs,
                anchor,  pct(anchor),
                matched, pct(matched),
                exact,   pct(exact),
                approx,  pct(approx),
                ambig,   pct(ambig),
                nodonor, pct(nodonor),
                nomatch, pct(nomatch));
        } else {
            std::snprintf(buf, sizeof(buf),
                "[bardmux] %8ld reads | anchors %5.1f%% | CB matched %5.1f%% "
                "(exact %5.1f%% approx %5.1f%%) | ambig %5.1f%% | %.0f rd/s",
                n,
                pct(anchor),
                pct(matched),
                pct(exact),
                pct(approx),
                pct(ambig),
                rate);
        }

        if (tty_ && !final)
            std::fprintf(stderr, "\r%-120s", buf);   // overwrite in place on TTY
        else
            std::fprintf(stderr, "%s\n", buf);
    }
};
