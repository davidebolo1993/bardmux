#include <algorithm>
#include <cctype>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
#include <unistd.h>

namespace fs = std::filesystem;

struct Config {
    std::string bardmux_file;
    std::string bardmux2_file;
    std::string wf_file;
    std::string tmp_dir;
    std::string out_file;
    std::string disagreements_file;
    std::string wf_only_status_file;
    std::string bardmux_unassigned = "unassigned";
    std::string sort_mem;
    int         threads = 1;
    int         max_disagreements = 20;
    uint64_t    progress_rows = 5000000ULL;
    bool        strip_wf_suffix = true;
    bool        keep_temp = false;
    bool        progress = true;
};

struct BardmuxParseStats {
    uint64_t input_rows          = 0;
    uint64_t status_rows         = 0;
    uint64_t unique_assigned_rows = 0;
    uint64_t malformed_rows      = 0;
};

struct WfParseStats {
    uint64_t input_rows          = 0;
    uint64_t assigned_rows       = 0;
    uint64_t malformed_rows      = 0;
};

struct BardmuxAssignedRec {
    std::string id;
    std::string cb;
    int         edit = -1;
};

struct WfRec {
    std::string id;
    std::string cb;
};

struct BardmuxStatusRec {
    std::string id;
    std::string status;
};

struct JoinStats {
    uint64_t bardmux_distinct_ids = 0;
    uint64_t wf_distinct_ids      = 0;
    uint64_t both_ids             = 0;
    uint64_t bardmux_only_ids     = 0;
    uint64_t wf_only_ids          = 0;

    uint64_t strict_comparable_ids = 0;
    uint64_t strict_concordant_ids = 0;

    uint64_t bardmux_dup_ids       = 0;
    uint64_t wf_dup_ids            = 0;
    uint64_t bardmux_dup_conflict  = 0;
    uint64_t wf_dup_conflict       = 0;

    // edit_distance -> (strict_concordant, strict_discordant)
    std::map<int, std::pair<uint64_t, uint64_t>> strict_edit_hist;
    std::map<std::string, uint64_t> wf_only_status_counts;

    struct Disagreement {
        std::string id;
        std::string bardmux_cb;
        std::string wf_cb;
        int         bardmux_edit = -1;
    };
    std::vector<Disagreement> disagreement_examples;
};

struct B2JoinStats {
    uint64_t a_assigned_distinct_ids = 0;
    uint64_t b_assigned_distinct_ids = 0;
    uint64_t assigned_intersection_ids = 0;
    uint64_t a_only_assigned_ids = 0;
    uint64_t b_only_assigned_ids = 0;

    uint64_t assigned_strict_comparable_ids = 0;
    uint64_t assigned_cb_match_ids = 0;
    uint64_t assigned_cb_mismatch_ids = 0;

    uint64_t a_dup_assigned_groups = 0;
    uint64_t b_dup_assigned_groups = 0;
    uint64_t a_dup_assigned_conflicting_groups = 0;
    uint64_t b_dup_assigned_conflicting_groups = 0;

    uint64_t a_status_distinct_ids = 0;
    uint64_t b_status_distinct_ids = 0;
    uint64_t shared_status_ids = 0;
    uint64_t a_only_status_ids = 0;
    uint64_t b_only_status_ids = 0;

    uint64_t anchor_found_both_shared_ids = 0;
    uint64_t anchor_found_a_only_shared_ids = 0;
    uint64_t anchor_found_b_only_shared_ids = 0;
    uint64_t anchor_found_neither_shared_ids = 0;

    std::map<int, std::pair<uint64_t, uint64_t>> strict_edit_hist_a; // editA -> (match, mismatch)
    std::map<std::string, uint64_t> status_transition_counts; // "statusA\tstatusB" -> count

    struct Disagreement {
        std::string id;
        std::string a_cb;
        std::string b_cb;
        int         a_edit = -1;
        int         b_edit = -1;
    };
    std::vector<Disagreement> disagreement_examples;
};

static std::string fmt_u64(uint64_t x) {
    std::string s = std::to_string(x);
    int insert_pos = static_cast<int>(s.size()) - 3;
    while (insert_pos > 0) {
        s.insert(static_cast<std::size_t>(insert_pos), ",");
        insert_pos -= 3;
    }
    return s;
}

class PhaseProgress {
public:
    PhaseProgress(const std::string& phase,
                  bool enabled,
                  uint64_t every_rows,
                  uint64_t total_bytes = 0)
        : phase_(phase),
          enabled_(enabled),
          every_rows_(std::max<uint64_t>(1, every_rows)),
          total_bytes_(total_bytes),
          tty_(isatty(fileno(stderr))),
          t0_(std::chrono::steady_clock::now())
    {
        if (enabled_) {
            std::fprintf(stderr, "[compare] %s started\n", phase_.c_str());
            std::fflush(stderr);
        }
    }

    void tick_row(uint64_t bytes = 0) {
        if (!enabled_) return;
        ++rows_;
        bytes_ += bytes;
        if (rows_ - last_print_rows_ >= every_rows_) {
            print(false);
            last_print_rows_ = rows_;
        }
    }

    void print_done(const std::string& suffix = "") {
        if (!enabled_) return;
        print(true, suffix);
    }

private:
    std::string phase_;
    bool enabled_;
    uint64_t every_rows_;
    uint64_t total_bytes_;
    bool tty_;
    uint64_t rows_ = 0;
    uint64_t bytes_ = 0;
    uint64_t last_print_rows_ = 0;
    std::chrono::steady_clock::time_point t0_;

    void print(bool final, const std::string& suffix = "") const {
        auto now = std::chrono::steady_clock::now();
        double sec = std::chrono::duration<double>(now - t0_).count();
        if (sec < 1e-6) sec = 1e-6;
        double rows_s = static_cast<double>(rows_) / sec;
        double mb_s = (static_cast<double>(bytes_) / (1024.0 * 1024.0)) / sec;
        double pct = (total_bytes_ > 0) ? (100.0 * static_cast<double>(bytes_) / static_cast<double>(total_bytes_)) : -1.0;

        std::ostringstream oss;
        oss << "[compare] " << phase_ << " | rows " << fmt_u64(rows_)
            << " | " << std::fixed << std::setprecision(0) << rows_s << " rows/s";
        if (bytes_ > 0) {
            oss << " | " << std::setprecision(1) << mb_s << " MiB/s";
        }
        if (pct >= 0.0) {
            if (pct > 100.0) pct = 100.0;
            oss << " | " << std::setprecision(1) << pct << "%";
        }
        if (!suffix.empty()) oss << " | " << suffix;
        std::string line = oss.str();

        if (tty_ && !final) {
            std::fprintf(stderr, "\r%-140s", line.c_str());
        } else {
            if (tty_ && final) std::fprintf(stderr, "\r\033[2K");
            std::fprintf(stderr, "%s\n", line.c_str());
        }
        std::fflush(stderr);
    }
};

static void usage(const char* prog) {
    std::cerr
        << "Usage:\n"
        << "  " << prog << " --bardmux ASSIGN.tsv --wf WF.tsv --tmp TMP_DIR [options]\n"
        << "  " << prog << " --bardmux ASSIGN_A.tsv --bardmux2 ASSIGN_B.tsv --tmp TMP_DIR [options]\n\n"
        << "Required:\n"
        << "  --bardmux FILE      bardmux assignment TSV\n"
        << "  --wf FILE           wf-single-cell read->CB table (mode: bardmux_vs_wf)\n"
        << "  --bardmux2 FILE     second bardmux assignment TSV (mode: bardmux_vs_bardmux)\n"
        << "  --tmp DIR           Temporary directory for intermediates\n\n"
        << "Options:\n"
        << "  --out FILE          Write summary report to FILE (default: stdout)\n"
        << "  --threads INT       sort parallel threads (default: 1)\n"
        << "  --sort-mem STR      sort memory hint (GNU sort -S), e.g. 50% or 4G\n"
        << "  --bardmux-unassigned X  Token treated as unassigned in matched_cb (default: unassigned)\n"
        << "  --no-wf-strip-suffix  Keep wf read ids as-is (default strips trailing _<digits>)\n"
        << "  --disagreements FILE   Write strict discordant examples to FILE\n"
        << "  --wf-only-status FILE  In wf mode: write wf-only bardmux-status table.\n"
        << "                         In bardmux2 mode: write status-transition table.\n"
        << "  --max-disagreements INT  Max discordant examples kept (default: 20)\n"
        << "  --progress-rows INT   Progress update interval in rows (default: 5000000)\n"
        << "  --no-progress         Disable progress messages\n"
        << "  --keep-temp         Keep intermediate files\n"
        << "  -h, --help          Show this help\n";
}

static std::string shell_quote(const std::string& s) {
    std::string q = "'";
    for (char c : s) {
        if (c == '\'') q += "'\\''";
        else q += c;
    }
    q += "'";
    return q;
}

static int run_cmd(const std::string& cmd) {
    return std::system(cmd.c_str());
}

static std::string trim_copy(const std::string& s) {
    std::size_t b = 0;
    while (b < s.size() && std::isspace(static_cast<unsigned char>(s[b]))) ++b;
    std::size_t e = s.size();
    while (e > b && std::isspace(static_cast<unsigned char>(s[e - 1]))) --e;
    return s.substr(b, e - b);
}

static std::string normalize_bardmux_read_id(std::string id) {
    auto sp = id.find(' ');
    if (sp != std::string::npos) id.resize(sp);
    return id;
}

static bool is_all_digits(const std::string& s, std::size_t from) {
    if (from >= s.size()) return false;
    for (std::size_t i = from; i < s.size(); ++i) {
        if (!std::isdigit(static_cast<unsigned char>(s[i]))) return false;
    }
    return true;
}

static std::string normalize_wf_read_id(std::string id, bool strip_suffix) {
    if (!strip_suffix) return id;
    auto us = id.rfind('_');
    if (us != std::string::npos && is_all_digits(id, us + 1)) id.resize(us);
    return id;
}

static bool parse_int_safe(const std::string& s, int& out) {
    try {
        std::size_t pos = 0;
        int v = std::stoi(s, &pos);
        if (pos != s.size()) return false;
        out = v;
        return true;
    } catch (...) {
        return false;
    }
}

static bool parse_u64_safe(const std::string& s, uint64_t& out) {
    try {
        std::size_t pos = 0;
        unsigned long long v = std::stoull(s, &pos);
        if (pos != s.size()) return false;
        out = static_cast<uint64_t>(v);
        return true;
    } catch (...) {
        return false;
    }
}

static bool parse_args(int argc, char** argv, Config& cfg) {
    if (argc < 2) return false;
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if ((a == "-h") || (a == "--help")) return false;
        else if (a == "--bardmux" && i + 1 < argc) cfg.bardmux_file = argv[++i];
        else if (a == "--our" && i + 1 < argc) cfg.bardmux_file = argv[++i];
        else if (a == "--wf" && i + 1 < argc) cfg.wf_file = argv[++i];
        else if (a == "--bardmux2" && i + 1 < argc) cfg.bardmux2_file = argv[++i];
        else if (a == "--bardmux-b" && i + 1 < argc) cfg.bardmux2_file = argv[++i];
        else if (a == "--tmp" && i + 1 < argc) cfg.tmp_dir = argv[++i];
        else if (a == "--out" && i + 1 < argc) cfg.out_file = argv[++i];
        else if (a == "--threads" && i + 1 < argc) {
            int t = 1;
            if (!parse_int_safe(argv[++i], t)) return false;
            cfg.threads = std::max(1, t);
        }
        else if (a == "--sort-mem" && i + 1 < argc) cfg.sort_mem = argv[++i];
        else if (a == "--bardmux-unassigned" && i + 1 < argc) cfg.bardmux_unassigned = argv[++i];
        else if (a == "--our-unassigned" && i + 1 < argc) cfg.bardmux_unassigned = argv[++i];
        else if (a == "--no-wf-strip-suffix") cfg.strip_wf_suffix = false;
        else if (a == "--disagreements" && i + 1 < argc) cfg.disagreements_file = argv[++i];
        else if (a == "--wf-only-status" && i + 1 < argc) cfg.wf_only_status_file = argv[++i];
        else if (a == "--max-disagreements" && i + 1 < argc) {
            int d = 0;
            if (!parse_int_safe(argv[++i], d)) return false;
            cfg.max_disagreements = std::max(0, d);
        }
        else if (a == "--progress-rows" && i + 1 < argc) {
            uint64_t p = 0;
            if (!parse_u64_safe(argv[++i], p)) return false;
            cfg.progress_rows = std::max<uint64_t>(1, p);
        }
        else if (a == "--no-progress") cfg.progress = false;
        else if (a == "--keep-temp") cfg.keep_temp = true;
        else {
            std::cerr << "Unknown or incomplete argument: " << a << "\n";
            return false;
        }
    }
    if (cfg.bardmux_file.empty() || cfg.tmp_dir.empty()) return false;
    const bool has_wf = !cfg.wf_file.empty();
    const bool has_b2 = !cfg.bardmux2_file.empty();
    if (has_wf == has_b2) {
        std::cerr << "Provide exactly one of --wf or --bardmux2\n";
        return false;
    }
    return true;
}

static uint64_t try_get_file_size(const fs::path& p) {
    std::error_code ec;
    const auto sz = fs::file_size(p, ec);
    if (ec) return 0;
    return static_cast<uint64_t>(sz);
}

static bool normalize_bardmux_file(const Config& cfg,
                                   const std::string& input_file,
                                   const std::string& phase_label,
                                   const fs::path& out_assigned_norm,
                                   const fs::path& out_status_norm,
                                   BardmuxParseStats& stats) {
    std::ifstream in(input_file);
    if (!in.is_open()) {
        std::cerr << "Cannot open bardmux file: " << input_file << "\n";
        return false;
    }
    std::ofstream out_assigned(out_assigned_norm);
    if (!out_assigned.is_open()) {
        std::cerr << "Cannot write normalized assigned file: " << out_assigned_norm << "\n";
        return false;
    }
    std::ofstream out_status(out_status_norm);
    if (!out_status.is_open()) {
        std::cerr << "Cannot write normalized status file: " << out_status_norm << "\n";
        return false;
    }

    PhaseProgress prog(phase_label, cfg.progress, cfg.progress_rows, try_get_file_size(input_file));
    std::string line;
    bool first = true;
    while (std::getline(in, line)) {
        prog.tick_row(static_cast<uint64_t>(line.size() + 1));
        if (line.empty()) continue;
        ++stats.input_rows;

        if (first) {
            first = false;
            if (line.rfind("read_id\t", 0) == 0) continue;
        }

        std::string cols[9];
        int col = 0;
        std::size_t start = 0;
        for (std::size_t i = 0; i <= line.size(); ++i) {
            if (i == line.size() || line[i] == '\t') {
                if (col <= 8) cols[col] = line.substr(start, i - start);
                ++col;
                start = i + 1;
            }
        }
        if (col < 9) {
            ++stats.malformed_rows;
            continue;
        }

        const std::string status = cols[1];
        const std::string matched_cb = cols[7];
        const std::string edit_str = cols[8];

        std::string id = normalize_bardmux_read_id(cols[0]);
        if (id.empty()) {
            ++stats.malformed_rows;
            continue;
        }

        out_status << id << '\t' << status << '\n';
        ++stats.status_rows;

        if (status != "unique") continue;
        if (matched_cb.empty() || matched_cb == cfg.bardmux_unassigned) continue;

        int ed = -1;
        if (!parse_int_safe(edit_str, ed)) ed = -1;

        out_assigned << id << '\t' << matched_cb << '\t' << ed << '\n';
        ++stats.unique_assigned_rows;
    }
    std::ostringstream sfx;
    sfx << "unique_assigned=" << fmt_u64(stats.unique_assigned_rows)
        << ", malformed=" << fmt_u64(stats.malformed_rows);
    prog.print_done(sfx.str());
    return true;
}

static bool normalize_wf_file(const Config& cfg, const fs::path& out_norm, WfParseStats& stats) {
    std::ifstream in(cfg.wf_file);
    if (!in.is_open()) {
        std::cerr << "Cannot open wf file: " << cfg.wf_file << "\n";
        return false;
    }
    std::ofstream out(out_norm);
    if (!out.is_open()) {
        std::cerr << "Cannot write normalized wf file: " << out_norm << "\n";
        return false;
    }

    PhaseProgress prog("normalize wf", cfg.progress, cfg.progress_rows, try_get_file_size(cfg.wf_file));
    std::string line;
    bool first = true;
    while (std::getline(in, line)) {
        prog.tick_row(static_cast<uint64_t>(line.size() + 1));
        if (line.empty()) continue;
        ++stats.input_rows;

        std::size_t sep = line.find('\t');
        if (sep == std::string::npos) sep = line.find(' ');
        if (sep == std::string::npos) {
            ++stats.malformed_rows;
            continue;
        }

        std::string id = trim_copy(line.substr(0, sep));
        std::string cb = trim_copy(line.substr(sep + 1));

        if (first) {
            first = false;
            std::string id_low = id;
            std::string cb_low = cb;
            std::transform(id_low.begin(), id_low.end(), id_low.begin(),
                           [](unsigned char c){ return static_cast<char>(std::tolower(c)); });
            std::transform(cb_low.begin(), cb_low.end(), cb_low.begin(),
                           [](unsigned char c){ return static_cast<char>(std::tolower(c)); });
            if ((id_low == "read_id" || id_low == "readid") &&
                (cb_low.find("barcode") != std::string::npos || cb_low == "cb")) {
                continue;
            }
        }

        id = normalize_wf_read_id(id, cfg.strip_wf_suffix);
        if (id.empty() || cb.empty()) {
            ++stats.malformed_rows;
            continue;
        }

        out << id << '\t' << cb << '\n';
        ++stats.assigned_rows;
    }
    std::ostringstream sfx;
    sfx << "assigned=" << fmt_u64(stats.assigned_rows)
        << ", malformed=" << fmt_u64(stats.malformed_rows);
    prog.print_done(sfx.str());
    return true;
}

static bool sort_file(const Config& cfg, const fs::path& in, const fs::path& out) {
    auto t0 = std::chrono::steady_clock::now();
    if (cfg.progress) {
        std::cerr << "[compare] sort started: " << in << "\n";
    }

    auto run_sort = [&](bool with_tuning) -> int {
        std::ostringstream cmd;
        cmd << "LC_ALL=C sort -k1,1 -T " << shell_quote(cfg.tmp_dir) << " ";
        if (with_tuning) {
            if (!cfg.sort_mem.empty()) cmd << "-S " << shell_quote(cfg.sort_mem) << " ";
            if (cfg.threads > 1) cmd << "--parallel=" << cfg.threads << " ";
        }
        cmd << shell_quote(in.string()) << " -o " << shell_quote(out.string());
        return run_cmd(cmd.str());
    };

    int rc = run_sort(true);
    if (rc != 0 && (!cfg.sort_mem.empty() || cfg.threads > 1)) {
        std::cerr << "Sort with tuning flags failed, retrying without tuning flags...\n";
        rc = run_sort(false);
    }
    if (rc != 0) {
        std::cerr << "Sort failed for file: " << in << "\n";
        return false;
    }
    if (cfg.progress) {
        auto t1 = std::chrono::steady_clock::now();
        double sec = std::chrono::duration<double>(t1 - t0).count();
        std::cerr << "[compare] sort done: " << out << " (" << std::fixed << std::setprecision(1) << sec << " s)\n";
    }
    return true;
}

static bool read_next_bardmux_assigned(std::ifstream& in, BardmuxAssignedRec& rec) {
    std::string line;
    if (!std::getline(in, line)) return false;
    std::size_t t1 = line.find('\t');
    if (t1 == std::string::npos) return false;
    std::size_t t2 = line.find('\t', t1 + 1);
    if (t2 == std::string::npos) return false;
    rec.id = line.substr(0, t1);
    rec.cb = line.substr(t1 + 1, t2 - (t1 + 1));
    int ed = -1;
    if (!parse_int_safe(line.substr(t2 + 1), ed)) ed = -1;
    rec.edit = ed;
    return true;
}

static bool read_next_wf(std::ifstream& in, WfRec& rec) {
    std::string line;
    if (!std::getline(in, line)) return false;
    std::size_t t1 = line.find('\t');
    if (t1 == std::string::npos) return false;
    rec.id = line.substr(0, t1);
    rec.cb = line.substr(t1 + 1);
    return true;
}

static bool read_next_bardmux_status(std::ifstream& in, BardmuxStatusRec& rec) {
    std::string line;
    if (!std::getline(in, line)) return false;
    std::size_t t1 = line.find('\t');
    if (t1 == std::string::npos) return false;
    rec.id = line.substr(0, t1);
    rec.status = line.substr(t1 + 1);
    return true;
}

template <typename TRec, typename NextFn>
static std::vector<TRec> consume_group(TRec first, std::ifstream& in, NextFn next_fn, bool& has_next, TRec& next_out) {
    std::vector<TRec> g;
    g.push_back(std::move(first));
    while (true) {
        TRec tmp;
        if (!next_fn(in, tmp)) {
            has_next = false;
            return g;
        }
        if (tmp.id != g.front().id) {
            has_next = true;
            next_out = std::move(tmp);
            return g;
        }
        g.push_back(std::move(tmp));
    }
}

static bool all_same_cb_bardmux(const std::vector<BardmuxAssignedRec>& g) {
    if (g.empty()) return true;
    for (std::size_t i = 1; i < g.size(); ++i)
        if (g[i].cb != g[0].cb) return false;
    return true;
}

static bool all_same_cb_wf(const std::vector<WfRec>& g) {
    if (g.empty()) return true;
    for (std::size_t i = 1; i < g.size(); ++i)
        if (g[i].cb != g[0].cb) return false;
    return true;
}

template <typename TRec>
static int count_unique_cb(const std::vector<TRec>& g) {
    if (g.empty()) return 0;
    std::vector<std::string> cbs;
    cbs.reserve(g.size());
    for (const auto& r : g) cbs.push_back(r.cb);
    std::sort(cbs.begin(), cbs.end());
    cbs.erase(std::unique(cbs.begin(), cbs.end()), cbs.end());
    return static_cast<int>(cbs.size());
}

static int min_edit_in_group(const std::vector<BardmuxAssignedRec>& g) {
    int best = std::numeric_limits<int>::max();
    for (const auto& r : g) best = std::min(best, r.edit);
    return (best == std::numeric_limits<int>::max()) ? -1 : best;
}

static bool merge_join_assigned(const Config& cfg,
                                const fs::path& bardmux_sorted,
                                const fs::path& wf_sorted,
                                const fs::path& wf_only_ids_file,
                                JoinStats& js) {
    std::ifstream a(bardmux_sorted);
    std::ifstream b(wf_sorted);
    std::ofstream wf_only_out(wf_only_ids_file);
    if (!a.is_open() || !b.is_open() || !wf_only_out.is_open()) {
        std::cerr << "Cannot open assigned sorted intermediate files\n";
        return false;
    }
    PhaseProgress prog("merge assigned", cfg.progress, cfg.progress_rows, 0);

    BardmuxAssignedRec ar;
    WfRec  br;
    bool has_a = read_next_bardmux_assigned(a, ar);
    bool has_b = read_next_wf(b, br);

    while (has_a && has_b) {
        if (ar.id < br.id) {
            BardmuxAssignedRec next_a;
            auto ga = consume_group(std::move(ar), a, read_next_bardmux_assigned, has_a, next_a);
            ++js.bardmux_distinct_ids;
            ++js.bardmux_only_ids;
            if (ga.size() > 1) {
                ++js.bardmux_dup_ids;
                if (!all_same_cb_bardmux(ga)) ++js.bardmux_dup_conflict;
            }
            prog.tick_row();
            if (has_a) ar = std::move(next_a);
            continue;
        }
        if (br.id < ar.id) {
            WfRec next_b;
            auto gb = consume_group(std::move(br), b, read_next_wf, has_b, next_b);
            ++js.wf_distinct_ids;
            ++js.wf_only_ids;
            wf_only_out << gb.front().id << '\n';
            if (gb.size() > 1) {
                ++js.wf_dup_ids;
                if (!all_same_cb_wf(gb)) ++js.wf_dup_conflict;
            }
            prog.tick_row();
            if (has_b) br = std::move(next_b);
            continue;
        }

        BardmuxAssignedRec next_a;
        WfRec  next_b;
        auto ga = consume_group(std::move(ar), a, read_next_bardmux_assigned, has_a, next_a);
        auto gb = consume_group(std::move(br), b, read_next_wf, has_b, next_b);

        ++js.bardmux_distinct_ids;
        ++js.wf_distinct_ids;
        ++js.both_ids;

        if (ga.size() > 1) {
            ++js.bardmux_dup_ids;
            if (!all_same_cb_bardmux(ga)) ++js.bardmux_dup_conflict;
        }
        if (gb.size() > 1) {
            ++js.wf_dup_ids;
            if (!all_same_cb_wf(gb)) ++js.wf_dup_conflict;
        }

        const int uniq_bardmux = count_unique_cb(ga);
        const int uniq_wf = count_unique_cb(gb);

        if (uniq_bardmux == 1 && uniq_wf == 1) {
            ++js.strict_comparable_ids;
            const bool concordant = (ga.front().cb == gb.front().cb);
            const int ed = min_edit_in_group(ga);
            if (concordant) {
                ++js.strict_concordant_ids;
                js.strict_edit_hist[ed].first++;
            } else {
                js.strict_edit_hist[ed].second++;
                if (static_cast<int>(js.disagreement_examples.size()) < cfg.max_disagreements) {
                    JoinStats::Disagreement d;
                    d.id = ga.front().id;
                    d.bardmux_cb = ga.front().cb;
                    d.wf_cb = gb.front().cb;
                    d.bardmux_edit = ga.front().edit;
                    js.disagreement_examples.push_back(std::move(d));
                }
            }
        }

        prog.tick_row();
        if (has_a) ar = std::move(next_a);
        if (has_b) br = std::move(next_b);
    }

    while (has_a) {
        BardmuxAssignedRec next_a;
        auto ga = consume_group(std::move(ar), a, read_next_bardmux_assigned, has_a, next_a);
        ++js.bardmux_distinct_ids;
        ++js.bardmux_only_ids;
        if (ga.size() > 1) {
            ++js.bardmux_dup_ids;
            if (!all_same_cb_bardmux(ga)) ++js.bardmux_dup_conflict;
        }
        prog.tick_row();
        if (has_a) ar = std::move(next_a);
    }
    while (has_b) {
        WfRec next_b;
        auto gb = consume_group(std::move(br), b, read_next_wf, has_b, next_b);
        ++js.wf_distinct_ids;
        ++js.wf_only_ids;
        wf_only_out << gb.front().id << '\n';
        if (gb.size() > 1) {
            ++js.wf_dup_ids;
            if (!all_same_cb_wf(gb)) ++js.wf_dup_conflict;
        }
        prog.tick_row();
        if (has_b) br = std::move(next_b);
    }

    std::ostringstream sfx;
    sfx << "intersection=" << fmt_u64(js.both_ids)
        << ", wf_only=" << fmt_u64(js.wf_only_ids)
        << ", bardmux_only=" << fmt_u64(js.bardmux_only_ids);
    prog.print_done(sfx.str());
    return true;
}

static std::string collapse_status_group(const std::vector<BardmuxStatusRec>& group) {
    if (group.empty()) return "not_in_bardmux_table";
    std::vector<std::string> statuses;
    statuses.reserve(group.size());
    for (const auto& r : group) statuses.push_back(r.status);
    std::sort(statuses.begin(), statuses.end());
    statuses.erase(std::unique(statuses.begin(), statuses.end()), statuses.end());
    if (statuses.empty()) return "not_in_bardmux_table";
    if (statuses.size() == 1) return statuses.front();
    return "multiple_statuses";
}

static bool count_wf_only_statuses(const Config& cfg,
                                   const fs::path& wf_only_ids_file,
                                   const fs::path& bardmux_status_sorted,
                                   JoinStats& js) {
    std::ifstream wf_ids(wf_only_ids_file);
    std::ifstream st(bardmux_status_sorted);
    if (!wf_ids.is_open() || !st.is_open()) {
        std::cerr << "Cannot open wf-only/status files for status breakdown\n";
        return false;
    }

    std::string wf_id;
    bool has_wf = static_cast<bool>(std::getline(wf_ids, wf_id));
    PhaseProgress prog("wf-only status breakdown", cfg.progress, cfg.progress_rows, try_get_file_size(wf_only_ids_file));

    BardmuxStatusRec sr;
    bool has_st = read_next_bardmux_status(st, sr);

    while (has_wf) {
        if (!has_st || wf_id < sr.id) {
            js.wf_only_status_counts["not_in_bardmux_table"]++;
            prog.tick_row(static_cast<uint64_t>(wf_id.size() + 1));
            has_wf = static_cast<bool>(std::getline(wf_ids, wf_id));
            continue;
        }

        if (sr.id < wf_id) {
            BardmuxStatusRec next_sr;
            auto g = consume_group(std::move(sr), st, read_next_bardmux_status, has_st, next_sr);
            (void)g;
            if (has_st) sr = std::move(next_sr);
            continue;
        }

        BardmuxStatusRec next_sr;
        auto g = consume_group(std::move(sr), st, read_next_bardmux_status, has_st, next_sr);
        js.wf_only_status_counts[collapse_status_group(g)]++;
        prog.tick_row(static_cast<uint64_t>(wf_id.size() + 1));
        if (has_st) sr = std::move(next_sr);

        has_wf = static_cast<bool>(std::getline(wf_ids, wf_id));
    }

    std::ostringstream sfx;
    sfx << "wf_only_ids=" << fmt_u64(js.wf_only_ids)
        << ", status_bins=" << fmt_u64(static_cast<uint64_t>(js.wf_only_status_counts.size()));
    prog.print_done(sfx.str());
    return true;
}

static bool group_has_anchor(const std::vector<BardmuxStatusRec>& group) {
    for (const auto& r : group) {
        if (r.status != "no_anchor") return true;
    }
    return false;
}

static bool merge_join_assigned_b2(const Config& cfg,
                                   const fs::path& a_sorted,
                                   const fs::path& b_sorted,
                                   B2JoinStats& js) {
    std::ifstream a(a_sorted);
    std::ifstream b(b_sorted);
    if (!a.is_open() || !b.is_open()) {
        std::cerr << "Cannot open bardmux-vs-bardmux assigned sorted files\n";
        return false;
    }

    PhaseProgress prog("merge assigned (bardmux vs bardmux)", cfg.progress, cfg.progress_rows, 0);

    BardmuxAssignedRec ar;
    BardmuxAssignedRec br;
    bool has_a = read_next_bardmux_assigned(a, ar);
    bool has_b = read_next_bardmux_assigned(b, br);

    while (has_a && has_b) {
        if (ar.id < br.id) {
            BardmuxAssignedRec next_a;
            auto ga = consume_group(std::move(ar), a, read_next_bardmux_assigned, has_a, next_a);
            ++js.a_assigned_distinct_ids;
            ++js.a_only_assigned_ids;
            if (ga.size() > 1) {
                ++js.a_dup_assigned_groups;
                if (!all_same_cb_bardmux(ga)) ++js.a_dup_assigned_conflicting_groups;
            }
            prog.tick_row();
            if (has_a) ar = std::move(next_a);
            continue;
        }
        if (br.id < ar.id) {
            BardmuxAssignedRec next_b;
            auto gb = consume_group(std::move(br), b, read_next_bardmux_assigned, has_b, next_b);
            ++js.b_assigned_distinct_ids;
            ++js.b_only_assigned_ids;
            if (gb.size() > 1) {
                ++js.b_dup_assigned_groups;
                if (!all_same_cb_bardmux(gb)) ++js.b_dup_assigned_conflicting_groups;
            }
            prog.tick_row();
            if (has_b) br = std::move(next_b);
            continue;
        }

        BardmuxAssignedRec next_a;
        BardmuxAssignedRec next_b;
        auto ga = consume_group(std::move(ar), a, read_next_bardmux_assigned, has_a, next_a);
        auto gb = consume_group(std::move(br), b, read_next_bardmux_assigned, has_b, next_b);

        ++js.a_assigned_distinct_ids;
        ++js.b_assigned_distinct_ids;
        ++js.assigned_intersection_ids;

        if (ga.size() > 1) {
            ++js.a_dup_assigned_groups;
            if (!all_same_cb_bardmux(ga)) ++js.a_dup_assigned_conflicting_groups;
        }
        if (gb.size() > 1) {
            ++js.b_dup_assigned_groups;
            if (!all_same_cb_bardmux(gb)) ++js.b_dup_assigned_conflicting_groups;
        }

        const int uniq_a = count_unique_cb(ga);
        const int uniq_b = count_unique_cb(gb);
        if (uniq_a == 1 && uniq_b == 1) {
            ++js.assigned_strict_comparable_ids;
            const int ed_a = min_edit_in_group(ga);
            if (ga.front().cb == gb.front().cb) {
                ++js.assigned_cb_match_ids;
                js.strict_edit_hist_a[ed_a].first++;
            } else {
                ++js.assigned_cb_mismatch_ids;
                js.strict_edit_hist_a[ed_a].second++;
                if (static_cast<int>(js.disagreement_examples.size()) < cfg.max_disagreements) {
                    B2JoinStats::Disagreement d;
                    d.id = ga.front().id;
                    d.a_cb = ga.front().cb;
                    d.b_cb = gb.front().cb;
                    d.a_edit = ga.front().edit;
                    d.b_edit = gb.front().edit;
                    js.disagreement_examples.push_back(std::move(d));
                }
            }
        }

        prog.tick_row();
        if (has_a) ar = std::move(next_a);
        if (has_b) br = std::move(next_b);
    }

    while (has_a) {
        BardmuxAssignedRec next_a;
        auto ga = consume_group(std::move(ar), a, read_next_bardmux_assigned, has_a, next_a);
        ++js.a_assigned_distinct_ids;
        ++js.a_only_assigned_ids;
        if (ga.size() > 1) {
            ++js.a_dup_assigned_groups;
            if (!all_same_cb_bardmux(ga)) ++js.a_dup_assigned_conflicting_groups;
        }
        prog.tick_row();
        if (has_a) ar = std::move(next_a);
    }
    while (has_b) {
        BardmuxAssignedRec next_b;
        auto gb = consume_group(std::move(br), b, read_next_bardmux_assigned, has_b, next_b);
        ++js.b_assigned_distinct_ids;
        ++js.b_only_assigned_ids;
        if (gb.size() > 1) {
            ++js.b_dup_assigned_groups;
            if (!all_same_cb_bardmux(gb)) ++js.b_dup_assigned_conflicting_groups;
        }
        prog.tick_row();
        if (has_b) br = std::move(next_b);
    }

    std::ostringstream sfx;
    sfx << "assigned_intersection=" << fmt_u64(js.assigned_intersection_ids)
        << ", a_only=" << fmt_u64(js.a_only_assigned_ids)
        << ", b_only=" << fmt_u64(js.b_only_assigned_ids);
    prog.print_done(sfx.str());
    return true;
}

static bool compare_statuses_b2(const Config& cfg,
                                const fs::path& a_status_sorted,
                                const fs::path& b_status_sorted,
                                B2JoinStats& js) {
    std::ifstream a(a_status_sorted);
    std::ifstream b(b_status_sorted);
    if (!a.is_open() || !b.is_open()) {
        std::cerr << "Cannot open bardmux-vs-bardmux status sorted files\n";
        return false;
    }

    PhaseProgress prog("compare statuses (bardmux vs bardmux)", cfg.progress, cfg.progress_rows, 0);

    BardmuxStatusRec ar;
    BardmuxStatusRec br;
    bool has_a = read_next_bardmux_status(a, ar);
    bool has_b = read_next_bardmux_status(b, br);

    while (has_a && has_b) {
        if (ar.id < br.id) {
            BardmuxStatusRec next_a;
            auto ga = consume_group(std::move(ar), a, read_next_bardmux_status, has_a, next_a);
            ++js.a_status_distinct_ids;
            ++js.a_only_status_ids;
            (void)ga;
            prog.tick_row();
            if (has_a) ar = std::move(next_a);
            continue;
        }
        if (br.id < ar.id) {
            BardmuxStatusRec next_b;
            auto gb = consume_group(std::move(br), b, read_next_bardmux_status, has_b, next_b);
            ++js.b_status_distinct_ids;
            ++js.b_only_status_ids;
            (void)gb;
            prog.tick_row();
            if (has_b) br = std::move(next_b);
            continue;
        }

        BardmuxStatusRec next_a;
        BardmuxStatusRec next_b;
        auto ga = consume_group(std::move(ar), a, read_next_bardmux_status, has_a, next_a);
        auto gb = consume_group(std::move(br), b, read_next_bardmux_status, has_b, next_b);

        ++js.a_status_distinct_ids;
        ++js.b_status_distinct_ids;
        ++js.shared_status_ids;

        const std::string sa = collapse_status_group(ga);
        const std::string sb = collapse_status_group(gb);
        js.status_transition_counts[sa + "\t" + sb]++;

        const bool anchor_a = group_has_anchor(ga);
        const bool anchor_b = group_has_anchor(gb);
        if (anchor_a && anchor_b) ++js.anchor_found_both_shared_ids;
        else if (anchor_a) ++js.anchor_found_a_only_shared_ids;
        else if (anchor_b) ++js.anchor_found_b_only_shared_ids;
        else ++js.anchor_found_neither_shared_ids;

        prog.tick_row();
        if (has_a) ar = std::move(next_a);
        if (has_b) br = std::move(next_b);
    }

    while (has_a) {
        BardmuxStatusRec next_a;
        auto ga = consume_group(std::move(ar), a, read_next_bardmux_status, has_a, next_a);
        ++js.a_status_distinct_ids;
        ++js.a_only_status_ids;
        (void)ga;
        prog.tick_row();
        if (has_a) ar = std::move(next_a);
    }
    while (has_b) {
        BardmuxStatusRec next_b;
        auto gb = consume_group(std::move(br), b, read_next_bardmux_status, has_b, next_b);
        ++js.b_status_distinct_ids;
        ++js.b_only_status_ids;
        (void)gb;
        prog.tick_row();
        if (has_b) br = std::move(next_b);
    }

    std::ostringstream sfx;
    sfx << "shared_status_ids=" << fmt_u64(js.shared_status_ids)
        << ", both_anchor=" << fmt_u64(js.anchor_found_both_shared_ids);
    prog.print_done(sfx.str());
    return true;
}

static void write_wf_only_status_table(std::ostream& out, const JoinStats& js) {
    out << "wf_only_bardmux_status\tcount\tpct_of_wf_only\n";
    const double denom = static_cast<double>(js.wf_only_ids);

    std::vector<std::pair<std::string, uint64_t>> rows(js.wf_only_status_counts.begin(), js.wf_only_status_counts.end());
    std::sort(rows.begin(), rows.end(), [](const auto& a, const auto& b) {
        if (a.second != b.second) return a.second > b.second;
        return a.first < b.first;
    });

    for (const auto& kv : rows) {
        double pct = (denom > 0.0) ? (100.0 * static_cast<double>(kv.second) / denom) : 0.0;
        out << kv.first << '\t' << kv.second << '\t' << std::fixed << std::setprecision(4) << pct << '\n';
    }
}

static void write_summary(std::ostream& out,
                          const BardmuxParseStats& bardmux_ps,
                          const WfParseStats& wf_ps,
                          const JoinStats& js)
{
    const double eps = 1e-12;
    const uint64_t strict_discordant = js.strict_comparable_ids - js.strict_concordant_ids;

    const double strict_concordance = (js.strict_comparable_ids > 0)
        ? (100.0 * js.strict_concordant_ids / (js.strict_comparable_ids + eps)) : 0.0;
    const double jaccard = (js.bardmux_distinct_ids + js.wf_distinct_ids - js.both_ids > 0)
        ? (100.0 * js.both_ids / (js.bardmux_distinct_ids + js.wf_distinct_ids - js.both_ids + eps)) : 0.0;
    const double bardmux_vs_wf_cov = (js.wf_distinct_ids > 0)
        ? (100.0 * js.both_ids / (js.wf_distinct_ids + eps)) : 0.0;
    const double wf_vs_bardmux_cov = (js.bardmux_distinct_ids > 0)
        ? (100.0 * js.both_ids / (js.bardmux_distinct_ids + eps)) : 0.0;

    out << "metric\tvalue\n";
    out << "comparison_mode\tbardmux_vs_wf\n";
    out << "bardmux_input_rows\t" << bardmux_ps.input_rows << "\n";
    out << "bardmux_status_rows\t" << bardmux_ps.status_rows << "\n";
    out << "bardmux_unique_assigned_rows\t" << bardmux_ps.unique_assigned_rows << "\n";
    out << "bardmux_malformed_rows\t" << bardmux_ps.malformed_rows << "\n";
    out << "wf_input_rows\t" << wf_ps.input_rows << "\n";
    out << "wf_assigned_rows\t" << wf_ps.assigned_rows << "\n";
    out << "wf_malformed_rows\t" << wf_ps.malformed_rows << "\n";
    out << "bardmux_assigned_distinct_ids\t" << js.bardmux_distinct_ids << "\n";
    out << "wf_assigned_distinct_ids\t" << js.wf_distinct_ids << "\n";
    out << "intersection_assigned_ids\t" << js.both_ids << "\n";
    out << "bardmux_only_assigned_ids\t" << js.bardmux_only_ids << "\n";
    out << "wf_only_assigned_ids\t" << js.wf_only_ids << "\n";
    out << "strict_comparable_ids\t" << js.strict_comparable_ids << "\n";
    out << "strict_concordant_ids\t" << js.strict_concordant_ids << "\n";
    out << "strict_discordant_ids\t" << strict_discordant << "\n";
    out << "strict_concordance_pct\t" << std::fixed << std::setprecision(4) << strict_concordance << "\n";
    out << "assigned_id_jaccard_pct\t" << std::fixed << std::setprecision(4) << jaccard << "\n";
    out << "bardmux_vs_wf_assigned_coverage_pct\t" << std::fixed << std::setprecision(4) << bardmux_vs_wf_cov << "\n";
    out << "wf_vs_bardmux_assigned_coverage_pct\t" << std::fixed << std::setprecision(4) << wf_vs_bardmux_cov << "\n";
    out << "bardmux_duplicate_id_groups\t" << js.bardmux_dup_ids << "\n";
    out << "bardmux_duplicate_id_conflicting_cb_groups\t" << js.bardmux_dup_conflict << "\n";
    out << "wf_duplicate_id_groups\t" << js.wf_dup_ids << "\n";
    out << "wf_duplicate_id_conflicting_cb_groups\t" << js.wf_dup_conflict << "\n";

    out << "\nstrict_edit_distance\tconcordant_ids\tdiscordant_ids\n";
    for (const auto& kv : js.strict_edit_hist) {
        out << kv.first << '\t' << kv.second.first << '\t' << kv.second.second << "\n";
    }

    out << "\n";
    write_wf_only_status_table(out, js);

    out << "\n# disagreement_examples\n";
    out << "read_id\tbardmux_cb\twf_cb\tbardmux_edit\n";
    for (const auto& ex : js.disagreement_examples) {
        out << ex.id << '\t' << ex.bardmux_cb << '\t' << ex.wf_cb << '\t' << ex.bardmux_edit << "\n";
    }
}

static void write_disagreements_file(const fs::path& out_path, const JoinStats& js) {
    std::ofstream out(out_path);
    if (!out.is_open()) {
        std::cerr << "Cannot open disagreements file: " << out_path << "\n";
        return;
    }
    out << "read_id\tbardmux_cb\twf_cb\tbardmux_edit\n";
    for (const auto& ex : js.disagreement_examples) {
        out << ex.id << '\t' << ex.bardmux_cb << '\t' << ex.wf_cb << '\t' << ex.bardmux_edit << "\n";
    }
}

static void write_wf_only_status_file(const fs::path& out_path, const JoinStats& js) {
    std::ofstream out(out_path);
    if (!out.is_open()) {
        std::cerr << "Cannot open wf-only status file: " << out_path << "\n";
        return;
    }
    write_wf_only_status_table(out, js);
}

static void write_status_transition_table(std::ostream& out, const B2JoinStats& js) {
    out << "status_a\tstatus_b\tcount\tpct_of_shared_status_ids\n";
    const double denom = static_cast<double>(js.shared_status_ids);

    std::vector<std::pair<std::string, uint64_t>> rows(js.status_transition_counts.begin(),
                                                       js.status_transition_counts.end());
    std::sort(rows.begin(), rows.end(), [](const auto& a, const auto& b) {
        if (a.second != b.second) return a.second > b.second;
        return a.first < b.first;
    });

    for (const auto& kv : rows) {
        const std::string& key = kv.first;
        auto tab = key.find('\t');
        std::string sa = (tab == std::string::npos) ? key : key.substr(0, tab);
        std::string sb = (tab == std::string::npos) ? "" : key.substr(tab + 1);
        const double pct = (denom > 0.0) ? (100.0 * static_cast<double>(kv.second) / denom) : 0.0;
        out << sa << '\t' << sb << '\t' << kv.second << '\t'
            << std::fixed << std::setprecision(4) << pct << '\n';
    }
}

static void write_summary_b2(std::ostream& out,
                             const BardmuxParseStats& a_ps,
                             const BardmuxParseStats& b_ps,
                             const B2JoinStats& js) {
    const double eps = 1e-12;
    const double cb_match_pct = (js.assigned_strict_comparable_ids > 0)
        ? (100.0 * js.assigned_cb_match_ids / (js.assigned_strict_comparable_ids + eps)) : 0.0;
    const double jaccard = (js.a_assigned_distinct_ids + js.b_assigned_distinct_ids - js.assigned_intersection_ids > 0)
        ? (100.0 * js.assigned_intersection_ids /
           (js.a_assigned_distinct_ids + js.b_assigned_distinct_ids - js.assigned_intersection_ids + eps)) : 0.0;
    const double a_vs_b_cov = (js.b_assigned_distinct_ids > 0)
        ? (100.0 * js.assigned_intersection_ids / (js.b_assigned_distinct_ids + eps)) : 0.0;
    const double b_vs_a_cov = (js.a_assigned_distinct_ids > 0)
        ? (100.0 * js.assigned_intersection_ids / (js.a_assigned_distinct_ids + eps)) : 0.0;

    uint64_t exact_status_match = 0;
    for (const auto& kv : js.status_transition_counts) {
        auto tab = kv.first.find('\t');
        if (tab == std::string::npos) continue;
        if (kv.first.substr(0, tab) == kv.first.substr(tab + 1)) exact_status_match += kv.second;
    }
    const double exact_status_match_pct = (js.shared_status_ids > 0)
        ? (100.0 * exact_status_match / (js.shared_status_ids + eps)) : 0.0;
    const double both_anchor_pct = (js.shared_status_ids > 0)
        ? (100.0 * js.anchor_found_both_shared_ids / (js.shared_status_ids + eps)) : 0.0;
    const double a_only_anchor_pct = (js.shared_status_ids > 0)
        ? (100.0 * js.anchor_found_a_only_shared_ids / (js.shared_status_ids + eps)) : 0.0;
    const double b_only_anchor_pct = (js.shared_status_ids > 0)
        ? (100.0 * js.anchor_found_b_only_shared_ids / (js.shared_status_ids + eps)) : 0.0;
    const double neither_anchor_pct = (js.shared_status_ids > 0)
        ? (100.0 * js.anchor_found_neither_shared_ids / (js.shared_status_ids + eps)) : 0.0;

    out << "metric\tvalue\n";
    out << "comparison_mode\tbardmux_vs_bardmux\n";
    out << "bardmux_a_input_rows\t" << a_ps.input_rows << "\n";
    out << "bardmux_a_status_rows\t" << a_ps.status_rows << "\n";
    out << "bardmux_a_unique_assigned_rows\t" << a_ps.unique_assigned_rows << "\n";
    out << "bardmux_a_malformed_rows\t" << a_ps.malformed_rows << "\n";
    out << "bardmux_b_input_rows\t" << b_ps.input_rows << "\n";
    out << "bardmux_b_status_rows\t" << b_ps.status_rows << "\n";
    out << "bardmux_b_unique_assigned_rows\t" << b_ps.unique_assigned_rows << "\n";
    out << "bardmux_b_malformed_rows\t" << b_ps.malformed_rows << "\n";
    out << "bardmux_a_assigned_distinct_ids\t" << js.a_assigned_distinct_ids << "\n";
    out << "bardmux_b_assigned_distinct_ids\t" << js.b_assigned_distinct_ids << "\n";
    out << "assigned_intersection_ids\t" << js.assigned_intersection_ids << "\n";
    out << "bardmux_a_only_assigned_ids\t" << js.a_only_assigned_ids << "\n";
    out << "bardmux_b_only_assigned_ids\t" << js.b_only_assigned_ids << "\n";
    out << "assigned_strict_comparable_ids\t" << js.assigned_strict_comparable_ids << "\n";
    out << "assigned_cb_match_ids\t" << js.assigned_cb_match_ids << "\n";
    out << "assigned_cb_mismatch_ids\t" << js.assigned_cb_mismatch_ids << "\n";
    out << "assigned_cb_match_pct\t" << std::fixed << std::setprecision(4) << cb_match_pct << "\n";
    out << "assigned_id_jaccard_pct\t" << std::fixed << std::setprecision(4) << jaccard << "\n";
    out << "bardmux_a_vs_b_assigned_coverage_pct\t" << std::fixed << std::setprecision(4) << a_vs_b_cov << "\n";
    out << "bardmux_b_vs_a_assigned_coverage_pct\t" << std::fixed << std::setprecision(4) << b_vs_a_cov << "\n";
    out << "bardmux_a_duplicate_assigned_groups\t" << js.a_dup_assigned_groups << "\n";
    out << "bardmux_a_duplicate_assigned_conflicting_groups\t" << js.a_dup_assigned_conflicting_groups << "\n";
    out << "bardmux_b_duplicate_assigned_groups\t" << js.b_dup_assigned_groups << "\n";
    out << "bardmux_b_duplicate_assigned_conflicting_groups\t" << js.b_dup_assigned_conflicting_groups << "\n";
    out << "shared_status_ids\t" << js.shared_status_ids << "\n";
    out << "status_exact_match_ids\t" << exact_status_match << "\n";
    out << "status_exact_match_pct\t" << std::fixed << std::setprecision(4) << exact_status_match_pct << "\n";
    out << "anchor_found_both_shared_ids\t" << js.anchor_found_both_shared_ids << "\n";
    out << "anchor_found_a_only_shared_ids\t" << js.anchor_found_a_only_shared_ids << "\n";
    out << "anchor_found_b_only_shared_ids\t" << js.anchor_found_b_only_shared_ids << "\n";
    out << "anchor_found_neither_shared_ids\t" << js.anchor_found_neither_shared_ids << "\n";
    out << "anchor_found_both_shared_pct\t" << std::fixed << std::setprecision(4) << both_anchor_pct << "\n";
    out << "anchor_found_a_only_shared_pct\t" << std::fixed << std::setprecision(4) << a_only_anchor_pct << "\n";
    out << "anchor_found_b_only_shared_pct\t" << std::fixed << std::setprecision(4) << b_only_anchor_pct << "\n";
    out << "anchor_found_neither_shared_pct\t" << std::fixed << std::setprecision(4) << neither_anchor_pct << "\n";

    out << "\nstrict_edit_distance_a\tmatch_ids\tmismatch_ids\n";
    for (const auto& kv : js.strict_edit_hist_a) {
        out << kv.first << '\t' << kv.second.first << '\t' << kv.second.second << "\n";
    }

    out << "\n";
    write_status_transition_table(out, js);

    out << "\n# disagreement_examples\n";
    out << "read_id\tbardmux_a_cb\tbardmux_b_cb\tbardmux_a_edit\tbardmux_b_edit\n";
    for (const auto& ex : js.disagreement_examples) {
        out << ex.id << '\t' << ex.a_cb << '\t' << ex.b_cb << '\t'
            << ex.a_edit << '\t' << ex.b_edit << "\n";
    }
}

static void write_disagreements_file_b2(const fs::path& out_path, const B2JoinStats& js) {
    std::ofstream out(out_path);
    if (!out.is_open()) {
        std::cerr << "Cannot open disagreements file: " << out_path << "\n";
        return;
    }
    out << "read_id\tbardmux_a_cb\tbardmux_b_cb\tbardmux_a_edit\tbardmux_b_edit\n";
    for (const auto& ex : js.disagreement_examples) {
        out << ex.id << '\t' << ex.a_cb << '\t' << ex.b_cb << '\t'
            << ex.a_edit << '\t' << ex.b_edit << "\n";
    }
}

int main(int argc, char** argv) {
    Config cfg;
    if (!parse_args(argc, argv, cfg)) {
        usage(argv[0]);
        return 1;
    }

    try {
        fs::create_directories(cfg.tmp_dir);
    } catch (const std::exception& e) {
        std::cerr << "Cannot create tmp dir: " << cfg.tmp_dir << " (" << e.what() << ")\n";
        return 1;
    }

    const bool mode_wf = !cfg.wf_file.empty();
    std::vector<fs::path> temp_files;

    if (mode_wf) {
        const fs::path bardmux_assigned_norm   = fs::path(cfg.tmp_dir) / "bardmux.assigned.norm.tsv";
        const fs::path bardmux_status_norm     = fs::path(cfg.tmp_dir) / "bardmux.status.norm.tsv";
        const fs::path wf_norm                 = fs::path(cfg.tmp_dir) / "wf.assigned.norm.tsv";

        const fs::path bardmux_assigned_sorted = fs::path(cfg.tmp_dir) / "bardmux.assigned.norm.sorted.tsv";
        const fs::path bardmux_status_sorted   = fs::path(cfg.tmp_dir) / "bardmux.status.norm.sorted.tsv";
        const fs::path wf_sorted               = fs::path(cfg.tmp_dir) / "wf.assigned.norm.sorted.tsv";

        const fs::path wf_only_ids_file        = fs::path(cfg.tmp_dir) / "wf.only.assigned.ids.tsv";
        temp_files = {bardmux_assigned_norm, bardmux_status_norm, wf_norm,
                      bardmux_assigned_sorted, bardmux_status_sorted, wf_sorted, wf_only_ids_file};

        BardmuxParseStats bardmux_ps;
        WfParseStats wf_ps;

        if (!normalize_bardmux_file(cfg, cfg.bardmux_file, "normalize bardmux", bardmux_assigned_norm, bardmux_status_norm, bardmux_ps)) return 1;
        if (!normalize_wf_file(cfg, wf_norm, wf_ps)) return 1;

        if (!sort_file(cfg, bardmux_assigned_norm, bardmux_assigned_sorted)) return 1;
        if (!sort_file(cfg, bardmux_status_norm, bardmux_status_sorted)) return 1;
        if (!sort_file(cfg, wf_norm, wf_sorted)) return 1;

        JoinStats js;
        if (!merge_join_assigned(cfg, bardmux_assigned_sorted, wf_sorted, wf_only_ids_file, js)) return 1;
        if (!count_wf_only_statuses(cfg, wf_only_ids_file, bardmux_status_sorted, js)) return 1;

        if (!cfg.out_file.empty()) {
            std::ofstream ofs(cfg.out_file);
            if (!ofs.is_open()) {
                std::cerr << "Cannot open output report: " << cfg.out_file << "\n";
                return 1;
            }
            write_summary(ofs, bardmux_ps, wf_ps, js);
            std::cerr << "Summary written to: " << cfg.out_file << "\n";
        } else {
            write_summary(std::cout, bardmux_ps, wf_ps, js);
        }

        if (!cfg.disagreements_file.empty()) write_disagreements_file(cfg.disagreements_file, js);
        if (!cfg.wf_only_status_file.empty()) write_wf_only_status_file(cfg.wf_only_status_file, js);
    } else {
        const fs::path a_assigned_norm   = fs::path(cfg.tmp_dir) / "bardmux_a.assigned.norm.tsv";
        const fs::path a_status_norm     = fs::path(cfg.tmp_dir) / "bardmux_a.status.norm.tsv";
        const fs::path b_assigned_norm   = fs::path(cfg.tmp_dir) / "bardmux_b.assigned.norm.tsv";
        const fs::path b_status_norm     = fs::path(cfg.tmp_dir) / "bardmux_b.status.norm.tsv";

        const fs::path a_assigned_sorted = fs::path(cfg.tmp_dir) / "bardmux_a.assigned.norm.sorted.tsv";
        const fs::path a_status_sorted   = fs::path(cfg.tmp_dir) / "bardmux_a.status.norm.sorted.tsv";
        const fs::path b_assigned_sorted = fs::path(cfg.tmp_dir) / "bardmux_b.assigned.norm.sorted.tsv";
        const fs::path b_status_sorted   = fs::path(cfg.tmp_dir) / "bardmux_b.status.norm.sorted.tsv";

        temp_files = {a_assigned_norm, a_status_norm, b_assigned_norm, b_status_norm,
                      a_assigned_sorted, a_status_sorted, b_assigned_sorted, b_status_sorted};

        BardmuxParseStats a_ps, b_ps;
        if (!normalize_bardmux_file(cfg, cfg.bardmux_file, "normalize bardmux A", a_assigned_norm, a_status_norm, a_ps)) return 1;
        if (!normalize_bardmux_file(cfg, cfg.bardmux2_file, "normalize bardmux B", b_assigned_norm, b_status_norm, b_ps)) return 1;

        if (!sort_file(cfg, a_assigned_norm, a_assigned_sorted)) return 1;
        if (!sort_file(cfg, a_status_norm, a_status_sorted)) return 1;
        if (!sort_file(cfg, b_assigned_norm, b_assigned_sorted)) return 1;
        if (!sort_file(cfg, b_status_norm, b_status_sorted)) return 1;

        B2JoinStats js;
        if (!merge_join_assigned_b2(cfg, a_assigned_sorted, b_assigned_sorted, js)) return 1;
        if (!compare_statuses_b2(cfg, a_status_sorted, b_status_sorted, js)) return 1;

        if (!cfg.out_file.empty()) {
            std::ofstream ofs(cfg.out_file);
            if (!ofs.is_open()) {
                std::cerr << "Cannot open output report: " << cfg.out_file << "\n";
                return 1;
            }
            write_summary_b2(ofs, a_ps, b_ps, js);
            std::cerr << "Summary written to: " << cfg.out_file << "\n";
        } else {
            write_summary_b2(std::cout, a_ps, b_ps, js);
        }

        if (!cfg.disagreements_file.empty()) write_disagreements_file_b2(cfg.disagreements_file, js);
        if (!cfg.wf_only_status_file.empty()) {
            std::ofstream out(cfg.wf_only_status_file);
            if (!out.is_open()) {
                std::cerr << "Cannot open status-transition file: " << cfg.wf_only_status_file << "\n";
                return 1;
            }
            write_status_transition_table(out, js);
        }
    }

    if (!cfg.keep_temp) {
        std::error_code ec;
        for (const auto& p : temp_files) fs::remove(p, ec);
    }

    return 0;
}
