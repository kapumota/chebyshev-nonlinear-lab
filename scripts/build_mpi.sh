#!/usr/bin/env bash
set -euo pipefail
BUILD_DIR="${1:-build-mpi}"
CXX="${CXX:-mpicxx}"
cmake -S . -B "$BUILD_DIR" \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_CXX_COMPILER="$CXX" \
    -DNUMA_HAS_MPI=1
cmake --build "$BUILD_DIR" -j
echo "Build MPI OK -> $BUILD_DIR"
