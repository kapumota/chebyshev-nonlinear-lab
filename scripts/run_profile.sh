#!/usr/bin/env bash
set -euo pipefail
# Perfilado con perf. Requiere linux-tools-generic y permisos.
BUILD_DIR="${1:-build}"
BIN="${2:-$BUILD_DIR/demo_chebyshev}"

if [ ! -x "$BIN" ]; then
    echo "ERROR: $BIN no existe."
    exit 1
fi

if ! command -v perf >/dev/null 2>&1; then
    echo "ERROR: perf no instalado. Ejecuta: sudo apt install linux-tools-generic"
    exit 1
fi

echo ">>> perf record sobre $BIN"
perf record -g -o perf.data "$BIN"
echo ">>> perf report"
perf report -i perf.data --stdio | head -n 60
echo ">>> Resultado guardado en perf.data"
