#!/usr/bin/env bash
set -euo pipefail
BUILD_DIR="${1:-build}"
cd "$BUILD_DIR"
ctest --output-on-failure
