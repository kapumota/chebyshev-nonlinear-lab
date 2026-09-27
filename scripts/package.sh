#!/usr/bin/env bash
set -euo pipefail

OUT="${1:-chebyshev-nonlinear-lab.zip}"
ROOT_NAME="${2:-chebyshev-nonlinear-lab}"
OUT_ABS="$(realpath -m "$OUT")"

if ! command -v git >/dev/null 2>&1; then
    printf '%s\n' 'ERROR: git es obligatorio para crear un paquete reproducible.' >&2
    exit 1
fi

if ! command -v zip >/dev/null 2>&1; then
    printf '%s\n' 'ERROR: zip es obligatorio para crear el paquete.' >&2
    exit 1
fi

if [[ -n "$(git status --porcelain)" ]]; then
    printf '%s\n' 'ERROR: el working tree debe estar limpio antes de empaquetar.' >&2
    exit 1
fi

TMP_DIR="$(mktemp -d)"
trap 'rm -rf "$TMP_DIR"' EXIT

rm -f "$OUT_ABS"

git archive \
    --format=tar \
    --prefix="${ROOT_NAME}/" \
    HEAD \
    | tar -xf - -C "$TMP_DIR"

(
    cd "$TMP_DIR"
    zip -qr "$OUT_ABS" "$ROOT_NAME"
)

printf 'Paquete creado: %s\n' "$OUT_ABS"
sha256sum "$OUT_ABS"
ls -lh "$OUT_ABS"
