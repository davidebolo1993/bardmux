#include <iostream>
#include <fstream>
#include <string>
#include <chrono>
#include <thread>
#include <stdexcept>
#include "fastq_reader.h"
#include "anchor_finder.h"
#include "barcode_matcher.h"
#include "report_writer.h"
#include "thread_pool.h"
#include "progress.h"

struct Config {
    std::string input_file;
    std::string barcode_file;
    std::string output_file;
    std::string split_dir;

    std::string anchor1           = "CTACACGACGCTCTTCCGATCT";
    std::string anchor2           = "TTTCTTATATGGG";
    int         max_anchor_errors = 3;
    int         max_anchor_edits  = -1;
    int         cb_length         = 16;
    int         umi_length        = 12;
    int         anchor_gap_slack  = 20;
    bool        use_anchor_gap_slack = true;

    int  max_edit_distance = 2;
    int  min_margin        = 1;
    int  kmer_length       = 8;

    int  num_threads    = std::max(1u, std::thread::hardware_concurrency());
    int  batch_size     = 4096;
    int  progress_every = 100000;

    bool verbose = false;
};

void print_usage(const char* prog) {
    std::cerr
        << "Usage: " << prog << " [options] -i INPUT -b BARCODES\n\n"
        << "Required:\n"
        << "  -i FILE     Input FASTQ / FASTQ.gz (or - for stdin)\n"
        << "  -b FILE     Barcode whitelist (1-col or 2-col TSV)\n\n"
        << "Anchor:\n"
        << "  -A SEQ         Anchor 1  [CTACACGACGCTCTTCCGATCT]\n"
        << "  --anchor2 SEQ  Anchor 2  [TTTCTTATATGGG]\n"
        << "  -E INT         Max anchor mismatches [3]\n"
        << "  --anchor-edits INT  Max anchor edit distance (indel-aware, disabled by default)\n"
        << "  -C INT         CB length bp [16]\n"
        << "  -U INT         UMI length bp [12]\n"
        << "  --gap-slack INT  Extra bp allowed between anchors [20]\n"
        << "  --no-fallback   Disable extra gap slack (strict CB+UMI window)\n\n"
        << "Barcode matching:\n"
        << "  -e INT      Max edit distance [2]\n"
        << "  -m INT      Min margin for unambiguous call [1]\n"
        << "  -k INT      K-mer length for whitelist filter [8]\n\n"
        << "Output:\n"
        << "  -o FILE     TSV assignments [stdout]\n"
        << "  -d DIR      Stream unique matched reads to DIR/<sample>.fastq.gz\n\n"
        << "Performance:\n"
        << "  -t INT      Threads [hardware_concurrency]\n"
        << "  -B INT      Reads per batch [4096]\n"
        << "  -p INT      Progress interval in reads [100000] (0 = off)\n"
        << "  -v          Verbose startup info\n";
}

bool parse_arguments(int argc, char** argv, Config& cfg) {
    if (argc < 2) { print_usage(argv[0]); return false; }
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if      (a=="-i"  && i+1<argc) cfg.input_file        = argv[++i];
        else if (a=="-b"  && i+1<argc) cfg.barcode_file      = argv[++i];
        else if (a=="-o"  && i+1<argc) cfg.output_file       = argv[++i];
        else if (a=="-d"  && i+1<argc) cfg.split_dir         = argv[++i];
        else if (a=="-A"  && i+1<argc) cfg.anchor1           = argv[++i];
        else if (a=="--anchor2"&&i+1<argc) cfg.anchor2       = argv[++i];
        else if (a=="-E"  && i+1<argc) cfg.max_anchor_errors = std::stoi(argv[++i]);
        else if (a=="--anchor-edits" && i+1<argc) cfg.max_anchor_edits = std::stoi(argv[++i]);
        else if (a=="-C"  && i+1<argc) cfg.cb_length         = std::stoi(argv[++i]);
        else if (a=="-U"  && i+1<argc) cfg.umi_length        = std::stoi(argv[++i]);
        else if (a=="--gap-slack" && i+1<argc) cfg.anchor_gap_slack = std::stoi(argv[++i]);
        else if (a=="--no-fallback") cfg.use_anchor_gap_slack = false;
        else if (a=="-e"  && i+1<argc) cfg.max_edit_distance = std::stoi(argv[++i]);
        else if (a=="-m"  && i+1<argc) cfg.min_margin        = std::stoi(argv[++i]);
        else if (a=="-k"  && i+1<argc) cfg.kmer_length       = std::stoi(argv[++i]);
        else if (a=="-t"  && i+1<argc) cfg.num_threads       = std::stoi(argv[++i]);
        else if (a=="-B"  && i+1<argc) cfg.batch_size        = std::stoi(argv[++i]);
        else if (a=="-p"  && i+1<argc) cfg.progress_every    = std::stoi(argv[++i]);
        else if (a=="-v")              cfg.verbose            = true;
        else if (a=="-h"||a=="--help") { print_usage(argv[0]); return false; }
        else { std::cerr << "Unknown argument: " << a << "\n"; return false; }
    }
    if (cfg.input_file.empty() || cfg.barcode_file.empty()) {
        std::cerr << "Error: -i and -b are required\n";
        print_usage(argv[0]); return false;
    }
    if (cfg.max_anchor_errors < 0 || cfg.max_edit_distance < 0 ||
        cfg.max_anchor_edits < -1 ||
        cfg.cb_length <= 0 || cfg.umi_length < 0 || cfg.kmer_length <= 0 ||
        cfg.num_threads <= 0 || cfg.batch_size <= 0 ||
        cfg.progress_every < 0 || cfg.min_margin < 0 || cfg.anchor_gap_slack < 0) {
        std::cerr << "Error: invalid numeric argument value\n";
        return false;
    }
    return true;
}

struct BatchItem { FastqRecord record; };

static void process_batch(std::vector<BatchItem>  batch,
                           const AnchorFinder&     anchor_finder,
                           const BarcodeMatcher&   barcode_matcher,
                           ReportWriter&           report_writer,
                           ProgressTracker*        progress)
{
    std::vector<ReportChunkItem> chunk;
    chunk.reserve(batch.size());
    ProgressDelta delta;

    for (auto& item : batch) {
        AnchorMatch am = anchor_finder.find_anchors(item.record);

        ReportEntry entry;
        entry.read_id       = item.record.name;
        entry.anchors_found = am.found;
        entry.strand        = am.is_reverse_complement ? "RC" : "FWD";
        entry.anchor1_pos   = am.anchor1_pos;
        entry.anchor2_pos   = am.anchor2_pos;
        entry.extracted_cb  = am.extracted_cb;
        entry.extracted_umi = am.extracted_umi;

        // Determine call status
        if (!am.found) {
            entry.status = CallStatus::no_anchor;
        } else if (am.extracted_cb.empty()) {
            entry.status = CallStatus::no_cb;
        } else {
            BarcodeMatch bm = barcode_matcher.match_barcode(am.extracted_cb);
            entry.edit_distance = bm.edit_distance;
            entry.n_candidates  = bm.num_candidates;

            if (!bm.found) {
                entry.status = CallStatus::no_match;
            } else if (bm.ambiguous) {
                entry.status    = CallStatus::ambiguous;
                entry.ambiguous = true;
                // Intentionally leave matched_cb and sample blank:
                // the call is not trustworthy for a unique assignment.
            } else {
                entry.status       = CallStatus::unique;
                entry.cb_matched   = true;
                entry.matched_cb   = bm.matched_barcode;
                entry.sample       = bm.sample;
            }
        }

        bool exact    = entry.cb_matched && entry.edit_distance == 0;
        bool no_donor = entry.cb_matched && entry.sample.empty();

        ++delta.total;
        if (am.found)          ++delta.with_anchor;
        if (entry.cb_matched)  ++delta.cb_matched;
        if (exact)             ++delta.exact;
        if (entry.ambiguous)   ++delta.ambiguous;
        if (no_donor)          ++delta.no_donor;

        chunk.push_back(ReportChunkItem{std::move(entry), std::move(item.record)});
    }

    report_writer.add_entries(std::move(chunk));
    if (progress) progress->record_batch(delta);
}

int main(int argc, char** argv) {
    Config cfg;
    if (!parse_arguments(argc, argv, cfg)) return 1;

    try {
        if (cfg.verbose) {
            std::cerr << "=== bardmux v2.0.0 ===\n"
                      << "Input:     " << cfg.input_file   << "\n"
                      << "Barcodes:  " << cfg.barcode_file << "\n"
                      << "Output:    " << (cfg.output_file.empty() ? "stdout" : cfg.output_file) << "\n"
                      << "Split dir: " << (cfg.split_dir.empty()   ? "(none)"  : cfg.split_dir)  << "\n"
                      << "Threads:   " << cfg.num_threads  << "\n"
                      << "Anchor mode: " << (cfg.max_anchor_edits >= 0 ? "edit distance" : "mismatch")
                      << (cfg.max_anchor_edits >= 0 ? (" (k=" + std::to_string(cfg.max_anchor_edits) + ")")
                                                    : (" (k=" + std::to_string(cfg.max_anchor_errors) + ")"))
                      << "\n";
        }

        AnchorConfig acfg;
        acfg.anchor1    = cfg.anchor1;
        acfg.anchor2    = cfg.anchor2;
        acfg.max_errors = cfg.max_anchor_errors;
        acfg.max_edits  = cfg.max_anchor_edits;
        acfg.cb_length  = cfg.cb_length;
        acfg.umi_length = cfg.umi_length;
        acfg.gap_slack  = cfg.anchor_gap_slack;
        acfg.use_gap_slack = cfg.use_anchor_gap_slack;
        AnchorFinder anchor_finder(acfg);

        MatcherConfig mcfg;
        mcfg.max_edit_distance = cfg.max_edit_distance;
        mcfg.min_margin        = cfg.min_margin;
        mcfg.kmer_length       = cfg.kmer_length;
        BarcodeMatcher barcode_matcher(mcfg);

        if (!barcode_matcher.load_whitelist(cfg.barcode_file)) {
            std::cerr << "Error: cannot load whitelist: " << cfg.barcode_file << "\n";
            return 1;
        }
        if (cfg.verbose)
            std::cerr << "Loaded " << barcode_matcher.whitelist_size() << " barcodes\n";

        ReporterConfig rcfg;
        rcfg.output_file = cfg.output_file;
        rcfg.split_dir   = cfg.split_dir;
        rcfg.verbose     = cfg.verbose;
        ReportWriter report_writer(rcfg);

        std::unique_ptr<ProgressTracker> progress_owned;
        ProgressTracker* progress = nullptr;
        if (cfg.progress_every > 0) {
            progress_owned = std::make_unique<ProgressTracker>(cfg.progress_every);
            progress = progress_owned.get();
        }

        ThreadPool pool(cfg.num_threads);
        FastqReader reader(cfg.input_file);
        FastqRecord record;

        std::vector<BatchItem> batch;
        batch.reserve(cfg.batch_size);

        auto dispatch = [&]() {
            auto b = std::make_shared<std::vector<BatchItem>>(std::move(batch));
            batch.clear();
            batch.reserve(cfg.batch_size);
            pool.submit([b, &anchor_finder, &barcode_matcher, &report_writer, progress]() {
                process_batch(std::move(*b), anchor_finder, barcode_matcher,
                              report_writer, progress);
            });
        };

        while (reader.next(record)) {
            batch.push_back({std::move(record)});
            if (static_cast<int>(batch.size()) >= cfg.batch_size) dispatch();
        }
        if (!batch.empty()) dispatch();
        pool.wait_all();

        if (!report_writer.finalize()) {
            std::cerr << "Error: failed to write report\n"; return 1;
        }

        if (progress) progress->finalize();

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n"; return 1;
    }
    return 0;
}
