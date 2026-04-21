# bardmux — 10X Barcode/UMI Extraction for ONT Reads

`bardmux` extracts 10X cell barcodes (CB) and UMIs from Nanopore reads by finding two adapter anchors, then matches CBs against a whitelist with error tolerance.

It supports `FASTQ`, `FASTQ.gz`, and streamed input from `stdin` (`-i -`).

## Build

```bash
mkdir -p build
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

## Quick Start

```bash
./build/bardmux \
  -i reads.fastq.gz \
  -b whitelist.tsv \
  -o assignments.tsv \
  -t 8 -v
```

## Best Throughput Run (Parallel Decompression + Streaming + Split Output)

```bash
pigz -dc -p 24 reads.fastq.gz | \
./build/bardmux \
  -i - \
  -b whitelist.tsv \
  -o assignments.tsv \
  -d split_fastq \
  -t 24 \
  -B 16384 \
  -p 0
```

Notes:
- `pigz -p 24` parallelizes decompression outside the tool.
- `-i -` makes `bardmux` read from stdin.
- `-d split_fastq` writes unique matched reads as `split_fastq/<sample>.fastq.gz`.
- `-p 0` disables progress printing overhead.

## CLI Options

Required:
- `-i FILE` input FASTQ/FASTQ.gz, or `-` for stdin
- `-b FILE` barcode whitelist TSV (1-col barcode, or 2-col barcode+sample)

Anchor detection:
- `-A SEQ` anchor1 sequence (default `CTACACGACGCTCTTCCGATCT`)
- `--anchor2 SEQ` anchor2 sequence (default `TTTCTTATATGGG`)
- `-E INT` max mismatches per anchor (default `3`)
- `--anchor-edits INT` max anchor edit distance (indel-aware; default disabled)
- `-C INT` CB length (default `16`)
- `-U INT` UMI length (default `12`)
- `--gap-slack INT` extra allowed bp between anchors beyond `CB+UMI` (default `20`)
- `--no-fallback` strict anchor spacing: only `CB+UMI` window (no extra slack)

If `--anchor-edits` is set, anchor matching uses edit distance (indels allowed) and `-E` is not used.

Barcode matching:
- `-e INT` max edit distance to whitelist (default `2`)
- `-m INT` minimum margin for unique call (`d2-d1`, default `1`)
- `-k INT` k-mer length for prefilter (default `8`)

Execution/output:
- `-o FILE` output TSV (default stdout)
- `-d DIR` split unique reads by sample into `DIR/<sample>.fastq.gz`
- `-t INT` worker threads (default `hardware_concurrency`)
- `-B INT` reads per batch (default `4096`)
- `-p INT` progress interval in reads (default `100000`, `0` disables)
- `-v` verbose startup/summary

## Output TSV

Columns:
- `read_id`
- `status`
- `strand`
- `anchor1_pos`
- `anchor2_pos`
- `extracted_cb`
- `extracted_umi`
- `matched_cb`
- `edit_distance`
- `n_candidates`
- `sample`

### Status Meaning

- `unique`: anchors found, CB extracted, one accepted whitelist match
- `ambiguous`: best match is not unique under ambiguity rules
- `no_match`: anchors and CB found, but no whitelist entry within `-e`
- `no_cb`: anchors found but region too short to extract a full CB
- `no_anchor`: no valid anchor pair found

For unassigned calls, `matched_cb` and `sample` are emitted as `unassigned`.

## Algorithm Summary

1. Anchor detection (FWD then RC):
- Search anchor pairs with mismatch threshold `-E`, or with bounded edit distance when `--anchor-edits` is enabled.
- Enforce biological spacing constraint:
  - minimum region: `CB length`
  - maximum region: `CB + UMI + gap_slack` (or strict `CB+UMI` with `--no-fallback`)
- Extract region between anchors, reverse-complement in RC orientation, then split into CB and UMI.

2. Barcode matching:
- Tier 1: exact hash lookup.
- Tier 2: k-mer prefilter to generate candidates.
- Tier 3: edit distance:
  - Myers bit-vector for lengths `<= 64`
  - banded Levenshtein fallback for longer strings
- Ambiguity rule: unique if margin satisfies `d2 - d1 >= m`.

3. Parallel execution:
- One producer thread parses input and dispatches batches.
- Worker threads perform anchor + barcode matching.
- Output is buffered and written in batch chunks.

## References

- G. Myers, 1999. A fast bit-vector algorithm for approximate string matching based on dynamic programming.
- V. I. Levenshtein, 1966. Binary codes capable of correcting deletions, insertions, and reversals.
- E. Ukkonen, 1985. Algorithms for approximate string matching.

## Minimal Test

```bash
./build/bardmux \
  -i test/test.fastq \
  -b test/cb_list_wdonors.tsv \
  -o out.tsv -p 0
```

Expected test CB/UMI:
- CB: `ACAGCCGCAAACAACA`
- UMI: `GCAGTGTGCG`

## License

MIT
