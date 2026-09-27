#!/usr/bin/env bash
set -euo pipefail

fail=0

printf '%s\n' 'Verificando estado del repositorio...'

if [[ -n "$(git status --porcelain)" ]]; then
    printf '%s\n' '[FALLA] El working tree no esta limpio.'
    fail=1
else
    printf '%s\n' '[ok] Working tree limpio.'
fi

if git diff --check --quiet; then
    printf '%s\n' '[ok] git diff --check no detecta errores.'
else
    printf '%s\n' '[FALLA] git diff --check detecta errores.'
    fail=1
fi

if git ls-files | grep -Eq '(^|/)(build|build-mpi|build-verify)/'; then
    printf '%s\n' '[FALLA] Existen directorios de build versionados.'
    fail=1
else
    printf '%s\n' '[ok] No hay directorios de build versionados.'
fi

if git ls-files | grep -Eq '\.(zip|tar|tar\.gz|tgz|o|obj|a|so|dylib|dll|exe)$'; then
    printf '%s\n' '[FALLA] Existen binarios o paquetes accidentales versionados.'
    fail=1
else
    printf '%s\n' '[ok] No hay binarios o paquetes accidentales versionados.'
fi

required_files=(
    README.md
    LICENSE
    CITATION.cff
    CONTRIBUTING.md
    NOTICE.md
    docs/research/BASELINE.md
    docs/research/PROVENANCE.md
    docs/research/REPRODUCIBILITY.md
    docs/research/PUBLIC_RELEASE_CHECKLIST.md
)

for file in "${required_files[@]}"; do
    if [[ ! -f "$file" ]]; then
        printf '[FALLA] Falta %s\n' "$file"
        fail=1
    fi
done

if grep -q 'family-names: "Contributors"' CITATION.cff; then
    printf '%s\n' '[FALLA] CITATION.cff conserva metadata generica. Reemplazarla por autores reales antes de publicar.'
    fail=1
else
    printf '%s\n' '[ok] CITATION.cff no conserva el autor generico inicial.'
fi

if [[ "$fail" -ne 0 ]]; then
    printf '%s\n' 'Preflight publico: FAIL'
    exit 1
fi

printf '%s\n' 'Preflight publico estructural: PASS'
printf '%s\n' 'Faltan las verificaciones cientificas de R0-G10 antes de cambiar la visibilidad.'
