#!/usr/bin/env bash
# Captura el entorno de ejecución para acompañar resultados de campañas.
set -euo pipefail

OUT="${1:-env_$(date +%Y%m%d_%H%M%S).txt}"

{
    echo "fecha: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "commit: $(git rev-parse HEAD 2>/dev/null || echo desconocido)"
    if git rev-parse --is-inside-work-tree >/dev/null 2>&1; then
       if [[ -z "$(git status --porcelain)" ]]; then
           echo "arbol_limpio: si"
       else
           echo "arbol_limpio: no"
       fi
    else
       echo "arbol_limpio: desconocido"
    fi
    echo "kernel: $(uname -sr)"
    echo "cpu: $(lscpu 2>/dev/null | awk -F: '/Model name/ {gsub(/^ +/, "", $2); print $2; exit}')"
    echo "nucleos: $(nproc)"
    echo "memoria_kb: $(awk '/MemTotal/ {print $2}' /proc/meminfo 2>/dev/null || echo desconocido)"
    echo "g++: $(g++ --version 2>/dev/null | head -1 || echo no_disponible)"
    echo "cmake: $(cmake --version 2>/dev/null | head -1 || echo no_disponible)"
    echo "mpi: $(mpirun --version 2>/dev/null | head -1 || echo no_disponible)"
    if command -v nvidia-smi >/dev/null 2>&1; then
        echo "gpu: $(nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv,noheader | paste -sd ';' -)"
        echo "nvcc: $(nvcc --version 2>/dev/null | tail -1 || echo no_disponible)"
    else
        echo "gpu: no_disponible"
    fi
} > "$OUT"

echo "Entorno guardado en $OUT"
