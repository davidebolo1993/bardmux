#include <algorithm>
#include <cctype>
#include <cstdint>
#include <cstdlib>
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

namespace fs = std::filesystem;

struct Config {
    std::string bardmux_file;
    std::string wf_file;
    std::string tmp_dir;
    std::string out_file;
    std::string disagreements_file;
    std::string wf_only_status_file;
    std::string bardmux_unassigned = "unassigned";
    std::string sort_mem;
    int         threads = 1;
    int         max_disagreements = 20;
    bool        strip_wf_suffix = true;
    bool        keep_temp = false;
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

static void usage(const char* prog) {
    std::cerr
        << "Usage: " << prog << " --bardmux ASSIGN.tsv --wf WF.tsv --tmp TMP_DIR [options]\n\n"
        << "Required:\n"
        << "  --bardmux FILE      bardmux assignment TSV\n"
        << "  --wf FILE           wf-single-cell read->CB table\n"
        << "  --tmp DIR           Temporary directory for intermediates\n\n"
        << "Options:\n"
        << "  --out FILE          Write summary report to FILE (default: stdout)\n"
        << "  --threads INT       sort parallel threads (default: 1)\n"
        << "  --sort-mem STR      sort memory hint (GNU sort -S), e.g. 50% or 4G\n"
        << "  --bardmux-unassigned X  Token treated as unassigned in matched_cb (default: unassigned)\n"
        << "  --no-wf-strip-suffix  Keep wf read ids as-is (default strips trailing _<digits>)\n"
        << "  --disagreements FILE   Write strict discordant examples to FILE\n"
        << "  --wf-only-status FILE  Write wf-only bardmux-status table to FILE\n"
        << "  --max-disagreements INT  Max discordant examples kept (default: 20)\n"
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

static bool parse_args(int argc, char** argv, Config& cfg) {
    if (argc < 2) return false;
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if ((a == "-h") || (a == "--help")) return false;
        else if (a == "--bardmux" && i + 1 < argc) cfg.bardmux_file = argv[++i];
        else if (a == "--our" && i + 1 < argc) cfg.bardmux_file = argv[++i];
        else if (a == "--wf" && i + 1 < argc) cfg.wf_file = argv[++i];
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
        else if (a == "--keep-temp") cfg.keep_temp = true;
        else {
            std::cerr << "Unknown or incomplete argument: " << a << "\n";
            return false;
        }
    }
    if (cfg.bardmux_file.empty() || cfg.wf_file.empty() || cfg.tmp_dir.empty()) return false;
    return true;
}

static bool normalize_bardmux_file(const Config& cfg,
                                   const fs::path& out_assigned_norm,
                                   const fs::path& out_status_norm,
                                   BardmuxParseStats& stats) {
    std::ifstream in(cfg.bardmux_file);
    if (!in.is_open()) {
        std::cerr << "Cannot open bardmux file: " << cfg.bardmux_file << "\n";
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

    std::string line;
    bool first = true;
    while (std::getline(in, line)) {
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

    std::string line;
    bool first = true;
    while (std::getline(in, line)) {
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
    return true;
}

static bool sort_file(const Config& cfg, const fs::path& in, const fs::path& out) {
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
        if (has_b) br = std::move(next_b);
    }

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

static bool count_wf_only_statuses(const fs::path& wf_only_ids_file,
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

    BardmuxStatusRec sr;
    bool has_st = read_next_bardmux_status(st, sr);

    while (has_wf) {
        if (!has_st || wf_id < sr.id) {
            js.wf_only_status_counts["not_in_bardmux_table"]++;
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
        if (has_st) sr = std::move(next_sr);

        has_wf = static_cast<bool>(std::getline(wf_ids, wf_id));
    }

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

    const fs::path bardmux_assigned_norm   = fs::path(cfg.tmp_dir) / "bardmux.assigned.norm.tsv";
    const fs::path bardmux_status_norm     = fs::path(cfg.tmp_dir) / "bardmux.status.norm.tsv";
    const fs::path wf_norm                 = fs::path(cfg.tmp_dir) / "wf.assigned.norm.tsv";

    const fs::path bardmux_assigned_sorted = fs::path(cfg.tmp_dir) / "bardmux.assigned.norm.sorted.tsv";
    const fs::path bardmux_status_sorted   = fs::path(cfg.tmp_dir) / "bardmux.status.norm.sorted.tsv";
    const fs::path wf_sorted               = fs::path(cfg.tmp_dir) / "wf.assigned.norm.sorted.tsv";

    const fs::path wf_only_ids_file        = fs::path(cfg.tmp_dir) / "wf.only.assigned.ids.tsv";

    BardmuxParseStats bardmux_ps;
    WfParseStats wf_ps;

    if (!normalize_bardmux_file(cfg, bardmux_assigned_norm, bardmux_status_norm, bardmux_ps)) return 1;
    if (!normalize_wf_file(cfg, wf_norm, wf_ps)) return 1;

    if (!sort_file(cfg, bardmux_assigned_norm, bardmux_assigned_sorted)) return 1;
    if (!sort_file(cfg, bardmux_status_norm, bardmux_status_sorted)) return 1;
    if (!sort_file(cfg, wf_norm, wf_sorted)) return 1;

    JoinStats js;
    if (!merge_join_assigned(cfg, bardmux_assigned_sorted, wf_sorted, wf_only_ids_file, js)) return 1;
    if (!count_wf_only_statuses(wf_only_ids_file, bardmux_status_sorted, js)) return 1;

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

    if (!cfg.keep_temp) {
        std::error_code ec;
        fs::remove(bardmux_assigned_norm, ec);
        fs::remove(bardmux_status_norm, ec);
        fs::remove(wf_norm, ec);
        fs::remove(bardmux_assigned_sorted, ec);
        fs::remove(bardmux_status_sorted, ec);
        fs::remove(wf_sorted, ec);
        fs::remove(wf_only_ids_file, ec);
    }

    return 0;
}
