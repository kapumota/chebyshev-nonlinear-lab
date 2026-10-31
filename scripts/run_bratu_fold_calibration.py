#!/usr/bin/env python3
"""Ejecuta la calibración reproducible del fold de Bratu3D para C1."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess

SIZES = (31, 63, 127)
RUNS = (
    ("coarse", 0.25),
    ("fine", 0.125),
)

FIXED = {
    "lambda_start": 0.0,
    "lambda_stop": 10.2,
    "min_step": 1e-4,
    "spectral_stop": 5e-2,
    "eigen_rtol": 1e-7,
    "linear_rtol": 1e-10,
    "threads": 1,
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def git_output(*args: str) -> str:
    return subprocess.run(
        ["git", *args], check=True, text=True, capture_output=True
    ).stdout.strip()


def require_clean_tracked_tree() -> str:
    head = git_output("rev-parse", "HEAD")
    worktree = subprocess.run(["git", "diff", "--quiet"], check=False)
    index = subprocess.run(["git", "diff", "--cached", "--quiet"], check=False)
    if worktree.returncode != 0 or index.returncode != 0:
        raise SystemExit(
            "La calibración del fold exige archivos rastreados sin cambios respecto de HEAD."
        )
    return head


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--executable",
        default="build-petsc/calibrate_b6_arclength",
    )
    parser.add_argument(
        "--output-dir",
        default="experiments/c1/fold_calibration",
    )
    parser.add_argument(
        "--manifest",
        default="experiments/c1/fold_calibration/bratu_fold_run_manifest.json",
    )
    args = parser.parse_args()

    root = Path.cwd().resolve()
    executable = (root / args.executable).resolve()
    output_dir = (root / args.output_dir).resolve()
    manifest_path = (root / args.manifest).resolve()

    if not executable.is_file():
        raise SystemExit(f"No existe calibrate_b6: {executable}")

    source_commit = require_clean_tracked_tree()
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)

    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "1"
    env["OPENBLAS_NUM_THREADS"] = "1"
    env["MKL_NUM_THREADS"] = "1"

    runs = []
    for n in SIZES:
        for label, step in RUNS:
            stem = f"bratu_n{n}_{label}"
            csv_path = output_dir / f"{stem}.csv"
            stdout_path = output_dir / f"{stem}.stdout.txt"
            stderr_path = output_dir / f"{stem}.stderr.txt"

            command = [
                str(executable),
                "--n", str(n),
                "--lambda-start", str(FIXED["lambda_start"]),
                "--lambda-stop", str(FIXED["lambda_stop"]),
                "--step", str(step),
                "--min-step", str(FIXED["min_step"]),
                "--spectral-stop", str(FIXED["spectral_stop"]),
                "--eigen-rtol", str(FIXED["eigen_rtol"]),
                "--linear-rtol", str(FIXED["linear_rtol"]),
                "--out", str(csv_path),
            ]

            completed = subprocess.run(
                command,
                cwd=root,
                env=env,
                text=True,
                capture_output=True,
                check=False,
            )
            stdout_path.write_text(completed.stdout, encoding="utf-8")
            stderr_path.write_text(completed.stderr, encoding="utf-8")

            if completed.returncode != 0:
                raise SystemExit(
                    f"calibrate_b6 falló para n={n}, run={label}, "
                    f"returncode={completed.returncode}. "
                    f"Revise {stderr_path}."
                )
            if not csv_path.is_file():
                raise SystemExit(f"No se produjo CSV para n={n}, run={label}.")

            runs.append({
                "n": n,
                "label": label,
                "step": step,
                "csv": {
                    "path": str(csv_path.relative_to(root)),
                    "sha256": sha256(csv_path),
                },
                "stdout": {
                    "path": str(stdout_path.relative_to(root)),
                    "sha256": sha256(stdout_path),
                },
                "stderr": {
                    "path": str(stderr_path.relative_to(root)),
                    "sha256": sha256(stderr_path),
                },
                "returncode": completed.returncode,
            })

    manifest = {
        "schema_version": 1,
        "artifact": "C1_BRATU_FOLD_RUNS",
        "status": "CALIBRATION_ONLY",
        "campaign": "C1",
        "source_commit": source_commit,
        "independent_of_compared_solver": True,
        "continuation_method": "natural_then_pseudo_arclength_fallback",
        "executable": {
            "path": str(executable.relative_to(root)),
            "sha256": sha256(executable),
        },
        "fixed": FIXED,
        "runs": runs,
    }
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(
        "C1_BRATU_FOLD_RUNS=PASS "
        f"runs={len(runs)} source_commit={source_commit}"
    )
    print(f"manifest={manifest_path.relative_to(root)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
