#!/usr/bin/env bash
set -euo pipefail

BUILD_DIR="${1:-build-mpi}"

cmake -S . -B "$BUILD_DIR" \
    -DCMAKE_BUILD_TYPE=Release \
    -DNUMA_BUILD_MPI=ON \
    -DNUMA_BUILD_BENCH=OFF

cmake --build "$BUILD_DIR" -j"$(nproc)"
ctest --test-dir "$BUILD_DIR" -L mpi --output-on-failure

echo "Compilación y verificación MPI correctas en $BUILD_DIR"
