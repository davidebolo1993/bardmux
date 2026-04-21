#!/bin/bash
set -e
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD="$DIR/build"
mkdir -p "$BUILD"
cd "$BUILD"
cmake -DCMAKE_BUILD_TYPE=Release "$DIR"
make -j"$(nproc 2>/dev/null || echo 4)"
echo ""
echo "=== Build complete: $BUILD/bardmux ==="
echo "Run: $BUILD/bardmux -i reads.fastq.gz -b barcodes.tsv -o out.tsv -t 8 -v"
