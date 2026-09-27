#!/usr/bin/env bash
set -euo pipefail
OUT="${1:-chebyshev-solver.zip}"
rm -f "$OUT"
zip -r "$OUT" \
    CMakeLists.txt README.md .gitignore .clang-format \
    cmake include apps tests scripts .github \
    -x "*/build/*" "*/build-*/*" "*/.git/*" "*.zip"
echo "Paquete creado: $OUT"
ls -lh "$OUT"
