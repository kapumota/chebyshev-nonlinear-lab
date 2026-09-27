#!/usr/bin/env bash
set -euo pipefail
echo ">>> 1/4 Build secuencial"
bash scripts/build.sh build

echo ">>> 2/4 Tests"
bash scripts/run_tests.sh build

echo ">>> 3/4 Demos"
./build/demo_newton
./build/demo_broyden
./build/demo_kantorovich
./build/demo_chebyshev
./build/demo_chebyshev_omp

echo ">>> 4/4 Benchmarks"
bash scripts/run_bench.sh build

echo ">>> TODO OK"
