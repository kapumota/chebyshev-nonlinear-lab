#!/usr/bin/env bash
set -euo pipefail
BUILD_DIR="${1:-build-omp}"
cmake -S . -B "$BUILD_DIR" \
    -DCMAKE_BUILD_TYPE=Release \
    -DNUMA_BUILD_APPS=ON
cmake --build "$BUILD_DIR" -j
echo "Build OpenMP OK -> $BUILD_DIR"
