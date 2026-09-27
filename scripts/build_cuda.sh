#!/usr/bin/env bash
set -euo pipefail
BUILD_DIR="${1:-build-cuda}"
cmake -S . -B "$BUILD_DIR" \
    -DCMAKE_BUILD_TYPE=Release \
    -DNUMA_BUILD_CUDA=ON \
    -DNUMA_BUILD_TESTS=OFF \
    -DNUMA_BUILD_BENCH=OFF
cmake --build "$BUILD_DIR" -j
echo "Build CUDA OK -> $BUILD_DIR"
