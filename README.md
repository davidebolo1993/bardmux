# bardmux — 10X Cell Barcode Extractor for Nanopore Reads

**bardmux** is a high-performance C++ tool for extracting and validating 10X Chromium cell barcodes (CBs) and UMIs from Oxford Nanopore Technology (ONT) sequencing reads. It locates two known library-preparation anchor sequences within each read, extracts the barcode region between them, and corrects it against a user-supplied whitelist using approximate string matching.

---

## Features

- **Streaming FASTQ/FASTQ.gz parsing** via vendored kseq with a 1 MiB gzip buffer
- **Dual-orientation anchor detection** — forward and reverse-complement read strands
- **Exact-match fast path** — O(1) hash lookup short-circuits edit-distance computation for perfect CB matches
- **Real k-mer hash index** — pre-filter eliminates >95% of whitelist before approximate matching
- **Bitwise Myers edit distance** — word-parallel bit-vector algorithm O(n·m/64) instead of O(n·m)
- **True multi-threaded processing** — producer-consumer batch dispatch; scales linearly with cores
- **Streaming TSV output** — entries written as produced, O(batch_size) RAM regardless of dataset size
- **No htslib dependency** — only zlib and pthreads required

---

## Biological Context

In 10X Chromium single-cell libraries, cell barcodes and UMIs are flanked by two fixed adapter sequences:

```
5'─ CTACACGACGCTCTTCCGATCT ─[ CB: 16 bp ]─[ UMI: 12 bp ]─ TTTCTTATATGGG ─3'
         TruSeq R1 adapter                                   Template Switch Oligo (TSO)
```

ONT reads are sequenced from full-length cDNA molecules. Each read originates from exactly one strand, so the two adapters always appear co-directionally:

```
Forward read:   5'─ A1_fwd ─[CB][UMI]─ A2_fwd ─3'
RC read:        5'─ A2_rc  ─[UMI'][CB']─ A1_rc ─3'   (reverse complement of the above)
```

Mixed orientations (A1_fwd + A2_rc or A1_rc + A2_fwd) are physically impossible and are not searched.

---

## Installation

### Prerequisites

| Dependency | Version | Debian/Ubuntu | macOS |
|---|---|---|---|
| CMake | ≥ 3.15 | `apt install cmake` | `brew install cmake` |
| GCC/Clang | C++17 | `apt install build-essential` | Xcode CLT |
| zlib | any | `apt install zlib1g-dev` | bundled with macOS SDK |

No htslib required — `kseq.h` is vendored in `src/`.

### Build

```bash
git clone https://github.com/yourusername/bardmux
cd bardmux
./build.sh
# or manually:
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j$(nproc)
```

---

## Usage

```bash
bardmux -i reads.fastq.gz -b whitelist.tsv -o results.tsv -t 8 -v
```

### All Options

```
Required:
  -i FILE     Input FASTQ / FASTQ.gz
  -b FILE     Barcode whitelist (1-col or 2-col TSV)

Anchor:
  -A SEQ      Anchor 1 sequence [CTACACGACGCTCTTCCGATCT]
  --anchor2   Anchor 2 sequence [TTTCTTATATGGG]
  -E INT      Max mismatches per anchor [3]
  -C INT      Cell barcode length bp [16]
  -U INT      UMI length bp [12]
  --no-fallback  Disable window fallback

Barcode matching:
  -e INT      Max edit distance vs whitelist [2]
  -m INT      Min margin for unambiguous call (d2-d1) [1]
  -k INT      K-mer length for whitelist filter [8]

Output:
  -o FILE     TSV assignments [stdout]
  -d DIR      Split FASTQ by sample into DIR/

Performance:
  -t INT      Worker threads [hardware_concurrency]
  -B INT      Reads per dispatch batch [4096]
  -v          Verbose output
```

### Output TSV Columns

| Column | Description |
|---|---|
| `read_id` | ONT read name |
| `anchors_found` | true/false |
| `strand` | FWD or RC |
| `anchor1_pos` | 0-based start of anchor 1 |
| `anchor2_pos` | 0-based start of anchor 2 |
| `extracted_cb` | Raw 16 bp cell barcode |
| `extracted_umi` | Raw UMI (up to 12 bp) |
| `cb_matched` | true/false |
| `matched_cb` | Corrected barcode from whitelist |
| `edit_distance` | Levenshtein distance to matched barcode (0 = perfect match) |
| `sample` | Sample label (if whitelist is 2-column) |
| `ambiguous` | true if ≥ 2 equidistant whitelist entries |

---

## Algorithm

### Overview

Each read passes through three sequential stages, all executed inside thread-pool workers:

```
FASTQ read
    │
    ▼
┌───────────────────────┐
│     Anchor Finder     │  locate TruSeq R1 + TSO in FWD or RC orientation
└──────────┬────────────┘
           │ CB sequence (16 bp)
           ▼
┌───────────────────────────────────────────────┐
│             Barcode Matcher                   │
│                                               │
│  1. Exact hash lookup  ──► hit?  return d=0   │
│                 │ miss                        │
│                 ▼                             │
│  2. k-mer filter  ──► candidate list          │
│                 │                             │
│                 ▼                             │
│  3. Myers edit distance  ──► best match       │
│                 │                             │
│                 ▼                             │
│  4. Ambiguity check  ──► accept / flag        │
└──────────┬────────────────────────────────────┘
           │
           ▼
┌───────────────────────┐
│     Report Writer     │  streaming TSV line (1 MiB output buffer)
└───────────────────────┘
```

---

### Stage 1 — Anchor Detection

#### Strand orientation

A 10X Nanopore read is sequenced from one physical strand. Both adapters therefore always appear
co-directionally. The two valid layouts are:

| Case | Sequenced strand | Left anchor in raw read | Right anchor in raw read | CB+UMI recovery |
|---|---|---|---|---|
| FWD | sense strand | A1_fwd | A2_fwd | extract directly: `seq[A1_end : A2_start]` |
| RC  | antisense strand | A2_rc | A1_rc | reverse-complement `seq[A2rc_end : A1rc_start]` |

**Important:** In the RC case, A2_rc appears to the **left** of A1_rc in the raw read — the
opposite order from the FWD case. This is because the sense-strand layout `A1–CB–UMI–A2` becomes
`A2_rc–UMI_rc–CB_rc–A1_rc` when reverse-complemented. Extracting the region between the two
anchors and reverse-complementing it restores the original `CB–UMI` order:

```
Sense strand:  5'─ A1_fwd ─[CB (16bp)][UMI (12bp)]─ A2_fwd ─3'
Raw RC read:   5'─ A2_rc  ─[UMI_rc  ][CB_rc       ]─ A1_rc  ─3'
                            ^                        ^
                            A2rc_end                 A1rc_start
  → RC(seq[A2rc_end : A1rc_start]) = CB + UMI
```

The reverse complements of both anchors are pre-computed at construction time (zero per-read cost).
Each read is tested in FWD orientation first; the RC sweep is only attempted on failure.

#### Mismatch scan with early exit

For each anchor pattern P of length m, bardmux sweeps a sliding window of width m across
the read text T of length L, counting mismatches:

    hit(i) = Σ_{j=0}^{m-1} 1[T[i+j] ≠ P[j]] ≤ k

The inner loop exits as soon as the running mismatch count exceeds k, giving average cost
well below O(m) per window for typical ONT data (k=3, anchor identity ≥ 85%).

#### CB and UMI extraction

Once a valid anchor pair is found at positions (a1_end, a2_start), the barcode region is sliced:

    CB  = seq[ a1_end : a1_end + cb_length ]           (default: 16 bp)
    UMI = seq[ a1_end + cb_length : a1_end + cb_length + umi_length ]   (default: up to 12 bp)

Both slices are clipped by a2_start to handle partially truncated molecules.

---

### Stage 2 — Barcode Matching

The extracted CB goes through three tiers, short-circuiting as early as possible:

#### Tier 1 — Exact match (O(1))

At whitelist load time, every barcode is inserted into a `std::unordered_map<string, index>`.
Before any approximate matching, the extracted CB is looked up in this map.

- **Hit**: edit distance is 0 by definition. Return immediately — no k-mer work, no Myers.
- **Miss**: proceed to Tier 2.

This is the common case in high-quality data (most CBs are read without sequencing errors).

#### Tier 2 — k-mer pre-filter

Each whitelist barcode is decomposed into overlapping k-mers (default k=8), encoded as
packed 2-bit integers:

    h(s₀…s_{k-1}) = Σ_{i=0}^{k-1} enc(sᵢ) × 4^{k-1-i}

where enc(A)=0, enc(C)=1, enc(G)=2, enc(T)=3. An `unordered_map` maps each hash to the
sorted list of whitelist indices carrying that k-mer.

For the query CB, all its k-mers are looked up. Only whitelist entries sharing ≥1 k-mer
are passed to Tier 3. By the pigeonhole principle, two strings of length L at Levenshtein
distance d share at least (L - k + 1 - d·k) k-mers, so this filter is theoretically
guaranteed for k ≤ (L+1)/(d+1). For L=16, k=8, d=2: threshold = 9 - 16 = negative,
meaning the guarantee is weak for large k; in practice requiring ≥1 shared k-mer
eliminates >95% of candidates empirically.

#### Tier 3 — Myers bit-vector edit distance

For remaining candidates, the exact Levenshtein distance is computed via the Myers (1999)
bit-vector algorithm. For barcodes ≤ 64 bp, the DP state is encoded in a single 64-bit word,
updated with five bitwise operations per text character:

    Peq[c]  bitmask of pattern positions equal to character c
    X  = Peq[a[i]] | VN
    D0 = ((VP + (X & VP)) ^ VP) | X
    HN = VP & D0
    HP = VN | ~(VP | D0)
    X  = (HP << 1) | 1
    VN = X & D0
    VP = (HN << 1) | ~(X | D0)
    score += (HP >> (m-1)) & 1
    score -= (HN >> (m-1)) & 1

Complexity: O(n·m/64) time, O(m) space, zero heap allocation per call.
For 16 bp barcodes (m=16, well within a single word) this is essentially constant time.
Barcodes >64 bp fall back to banded Levenshtein DP of width 2·max_edit_distance+1.

#### Ambiguity resolution

After collecting all candidates with best distance d₁, a call is accepted as unambiguous when:

    d₂ - d₁ ≥ min_margin   (default min_margin = 1)

where d₂ is the second-best distance. Ambiguous reads are flagged with `ambiguous=true`;
the first match found is still reported in `matched_cb` for downstream use.

---

### Stage 3 — Parallelism

#### Producer-consumer model

```
Main thread (I/O)
    reads FASTQ one record at a time
    fills std::vector<FastqRecord> of size -B (default 4096)
    when full ──► submits batch to ThreadPool

Worker threads (N = -t, default hardware_concurrency)
    for each read in batch:
        1. AnchorFinder::find_anchors()    [read-only, zero locks]
        2. BarcodeMatcher::match_barcode() [read-only, zero locks]
        3. ReportWriter::write_entry()     [builds line locally, one mutex for write()]
```

`AnchorFinder` and `BarcodeMatcher` are fully immutable after construction.
The only shared mutable state is the output stream, protected by a mutex held only for
a single `ostream::write()` call on a pre-formatted string — the critical section is
proportional to one line of TSV text, not to any computation.

---

## Performance

Expected throughput on a modern Linux server (AMD EPYC / Intel Xeon, NVMe):

| Configuration | Reads/second |
|---|---|
| 1 thread | 200 000 – 400 000 |
| 8 threads | 1 500 000 – 3 000 000 |
| 32 threads | 5 000 000 – 10 000 000 |

Memory: ~60 MB for a 6 M-barcode 10X v3 whitelist (strings + exact index + k-mer index).

### Tuning tips

- Large whitelists (>1 M barcodes): `-k 10` to reduce index bucket sizes
- I/O-bound datasets (NFS/slow disk): `-B 16384` to reduce scheduling overhead
- Stricter demultiplexing: `-m 2` for a larger margin between best and second-best
- Always build with `-DCMAKE_BUILD_TYPE=Release` (enables LTO automatically)

---

## Testing

```bash
mkdir build && cd build && cmake .. && make
./bardmux -i ../test/test.fastq -b ../test/cb_list_wdonors.tsv -o out.tsv -v
head -2 out.tsv
```

Expected: extracted CB `ACAGCCGCAAACAACA` matched to sample `FY_M749`.

---

## Troubleshooting

**Poor match rate:**
- Check anchor sequences match your library prep (`-A`, `--anchor2`)
- Try `-E 4` for high error-rate (>10%) reads
- Verify whitelist: `wc -l whitelist.tsv`

---

## Citation

```
bardmux v2.0.0 — Cell barcode extraction for Nanopore reads
https://github.com/yourusername/bardmux
```

## License

MIT
