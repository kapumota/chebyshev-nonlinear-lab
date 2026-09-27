#!/usr/bin/env bash
set -euo pipefail
BUILD_DIR="${1:-build-asan}"
cmake -S . -B "$BUILD_DIR" \
    -DCMAKE_BUILD_TYPE=Debug \
    -DNUMA_ENABLE_SANITIZERS=ON \
    -DNUMA_BUILD_BENCH=OFF
cmake --build "$BUILD_DIR" -j
cd "$BUILD_DIR"
ASAN_OPTIONS=detect_leaks=1 UBSAN_OPTIONS=print_stacktrace=1 \
    ctest --output-on-failure
