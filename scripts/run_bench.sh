#!/usr/bin/env bash
set -euo pipefail
BUILD_DIR="${1:-build}"
BENCH="$BUILD_DIR/numa_bench"
if [ ! -x "$BENCH" ]; then
    echo "ERROR: $BENCH no existe. Compila con NUMA_BUILD_BENCH=ON"
    exit 1
fi
"$BENCH" --benchmark_format=console --benchmark_min_time=0.2s
