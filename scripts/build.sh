#!/usr/bin/env bash
set -euo pipefail
BUILD_DIR="${1:-build}"
cmake -S . -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE=Release
cmake --build "$BUILD_DIR" -j
echo "Build OK -> $BUILD_DIR"
