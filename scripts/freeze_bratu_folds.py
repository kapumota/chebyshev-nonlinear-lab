#!/usr/bin/env python3
"""Congela lambda_fold_h desde dos refinamientos independientes de continuation."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import tempfile

SIZES = (31, 63, 127)
LABELS = ("coarse", "fine")
MAX_RELATIVE_DIFFERENCE = 0.0025
MIN_R2 = 0.98


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def load_points(path: Path, eigen_rtol: float, spectral_stop: float) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {
            "lambda",
            "spectral_gap",
            "eigen_residual",
            "continuation_iters",
            "n_F",
            "n_JVP",
            "n_prec",
        }
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise SystemExit(
                f"{path}: faltan columnas: {','.join(sorted(missing))}"
            )
        points = list(reader)

    if len(points) < 4:
        raise SystemExit(f"{path}: se requieren al menos cuatro puntos aceptados.")

    parsed = []
    previous_lambda = -math.inf
    for row in points:
        lam = float(row["lambda"])
        gap = float(row["spectral_gap"])
        residual = float(row["eigen_residual"])
        if not all(math.isfinite(x) for x in (lam, gap, residual)):
            raise SystemExit(f"{path}: valor no finito.")
        if lam <= previous_lambda:
            raise SystemExit(f"{path}: lambda no es estrictamente creciente.")
        if gap <= 0.0:
            raise SystemExit(f"{path}: spectral_gap debe ser positivo.")
        if residual > eigen_rtol:
            raise SystemExit(
                f"{path}: eigen_residual={residual:.6g} supera {eigen_rtol:.6g}."
            )
        parsed.append({"lambda": lam, "gap": gap, "residual": residual})
        previous_lambda = lam

    if parsed[-1]["gap"] > spectral_stop:
        raise SystemExit(
            f"{path}: el último spectral_gap={parsed[-1]['gap']:.6g} "
            f"no alcanza spectral_stop={spectral_stop:.6g}."
        )
    return parsed


def fit_proxy(points: list[dict]) -> tuple[float, float, float, float]:
    tail = points[-4:]
    xs = [p["lambda"] for p in tail]
    ys = [p["gap"] ** 2 for p in tail]
    count = float(len(xs))
    sx = sum(xs)
    sy = sum(ys)
    sxx = sum(x * x for x in xs)
    sxy = sum(x * y for x, y in zip(xs, ys))
    denominator = count * sxx - sx * sx
    if abs(denominator) < 1e-14:
        raise SystemExit("Ajuste degenerado de gap^2.")

    slope = (count * sxy - sx * sy) / denominator
    intercept = (sy - slope * sx) / count
    if not slope < 0.0:
        raise SystemExit("El ajuste local de gap^2 no tiene pendiente negativa.")

    root = -intercept / slope
    if not math.isfinite(root) or root <= xs[-1]:
        raise SystemExit("El proxy espectral no extrapola un fold posterior al último punto.")

    mean_y = sy / count
    ss_tot = sum((y - mean_y) ** 2 for y in ys)
    ss_res = sum(
        (y - (intercept + slope * x)) ** 2
        for x, y in zip(xs, ys)
    )
    r2 = 1.0 if ss_tot == 0.0 else 1.0 - ss_res / ss_tot
    if r2 < MIN_R2:
        raise SystemExit(
            f"R2={r2:.6g} inferior al mínimo preregistrado {MIN_R2:.6g}."
        )
    return root, r2, slope, intercept


def validate_manifest(data: dict, root: Path) -> dict[tuple[int, str], dict]:
    if data.get("schema_version") != 1:
        raise SystemExit("Manifest de folds con schema_version inválido.")
    if data.get("artifact") != "C1_BRATU_FOLD_RUNS":
        raise SystemExit("Artifact de runs Bratu inválido.")
    if data.get("status") != "CALIBRATION_ONLY":
        raise SystemExit("Los runs Bratu no están marcados CALIBRATION_ONLY.")
    if data.get("campaign") != "C1":
        raise SystemExit("Los runs Bratu no pertenecen a C1.")
    if data.get("independent_of_compared_solver") is not True:
        raise SystemExit("Los runs deben ser independientes del solver comparado.")

    fixed = data.get("fixed", {})
    if float(fixed.get("spectral_stop", -1.0)) != 5e-2:
        raise SystemExit("spectral_stop cambió respecto del protocolo.")
    if float(fixed.get("eigen_rtol", -1.0)) != 1e-7:
        raise SystemExit("eigen_rtol cambió respecto del protocolo.")
    if int(fixed.get("threads", 0)) != 1:
        raise SystemExit("El protocolo de fold exige un hilo.")

    expected_steps = {"coarse": 0.25, "fine": 0.125}
    observed = {}
    for item in data.get("runs", []):
        n = int(item["n"])
        label = str(item["label"])
        if n not in SIZES or label not in LABELS:
            raise SystemExit("Run Bratu inesperado.")
        key = (n, label)
        if key in observed:
            raise SystemExit("Run Bratu duplicado.")
        if float(item["step"]) != expected_steps[label]:
            raise SystemExit(f"Paso inesperado en n={n}, run={label}.")
        csv_path = root / item["csv"]["path"]
        if not csv_path.is_file():
            raise SystemExit(f"No existe {csv_path}.")
        if sha256(csv_path) != item["csv"]["sha256"]:
            raise SystemExit(f"Hash de CSV no coincide para n={n}, run={label}.")
        observed[key] = item

    expected = {(n, label) for n in SIZES for label in LABELS}
    if set(observed) != expected:
        raise SystemExit("El manifest no contiene exactamente seis runs Bratu.")
    return observed


def freeze(manifest_path: Path, output_path: Path) -> dict:
    root = Path.cwd().resolve()
    data = load_json(manifest_path)
    runs = validate_manifest(data, root)
    fixed = data["fixed"]
    eigen_rtol = float(fixed["eigen_rtol"])
    spectral_stop = float(fixed["spectral_stop"])

    folds = {}
    evidence = []
    for n in SIZES:
        proxies = {}
        for label in LABELS:
            item = runs[(n, label)]
            csv_path = root / item["csv"]["path"]
            points = load_points(csv_path, eigen_rtol, spectral_stop)
            proxy, r2, slope, intercept = fit_proxy(points)
            proxies[label] = {
                "lambda_fold_proxy": proxy,
                "r2_gap2_linear": r2,
                "slope": slope,
                "intercept": intercept,
                "last_lambda": points[-1]["lambda"],
                "last_spectral_gap": points[-1]["gap"],
                "max_eigen_residual": max(p["residual"] for p in points),
                "points": len(points),
                "csv_path": item["csv"]["path"],
                "csv_sha256": item["csv"]["sha256"],
            }

        coarse = proxies["coarse"]["lambda_fold_proxy"]
        fine = proxies["fine"]["lambda_fold_proxy"]
        relative = abs(fine - coarse) / abs(fine)
        if relative > MAX_RELATIVE_DIFFERENCE:
            raise SystemExit(
                f"n={n}: discrepancia relativa coarse/fine={relative:.6g} "
                f"supera {MAX_RELATIVE_DIFFERENCE:.6g}. "
                "Se requiere pseudo-arclength antes de C1."
            )

        folds[str(n)] = {
            "lambda_fold_h": fine,
            "selection": "fine_proxy_after_stability_gate",
            "coarse_proxy": coarse,
            "fine_proxy": fine,
            "relative_difference": relative,
            "fine_r2": proxies["fine"]["r2_gap2_linear"],
        }
        evidence.append({
            "n": n,
            "coarse": proxies["coarse"],
            "fine": proxies["fine"],
        })

    frozen = {
        "schema_version": 1,
        "artifact": "C1_BRATU_FOLDS",
        "status": "FROZEN_CALIBRATION_ONLY",
        "campaign": "C1",
        "source_commit": data["source_commit"],
        "independent_of_compared_solver": True,
        "freeze_rule": {
            "coarse_step": 0.25,
            "fine_step": 0.125,
            "eigen_rtol": 1e-7,
            "spectral_stop": 5e-2,
            "min_r2_gap2_linear": MIN_R2,
            "max_relative_difference": MAX_RELATIVE_DIFFERENCE,
            "selected_estimate": "fine_proxy",
            "pseudo_arclength_if_gate_fails": True,
        },
        "run_manifest": {
            "path": str(manifest_path),
            "sha256": sha256(manifest_path),
        },
        "folds": folds,
        "evidence": evidence,
    }

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        json.dumps(frozen, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return frozen


def self_test() -> int:
    with tempfile.TemporaryDirectory(prefix="c1-fold-selftest-") as tmp:
        temp = Path(tmp)
        old_cwd = Path.cwd()
        try:
            import os
            os.chdir(temp)
            runs = []
            for n in SIZES:
                true_fold = 9.9 + n * 1e-5
                for label, step, bias in (
                    ("coarse", 0.25, -0.0010),
                    ("fine", 0.125, 0.0),
                ):
                    path = temp / f"n{n}_{label}.csv"
                    lambdas = [
                        true_fold - 0.60,
                        true_fold - 0.40,
                        true_fold - 0.20,
                        true_fold - 0.05,
                    ]
                    with path.open("w", encoding="utf-8", newline="") as handle:
                        handle.write(
                            "lambda,spectral_gap,eigen_residual,"
                            "continuation_iters,n_F,n_JVP,n_prec\n"
                        )
                        shifted_fold = true_fold + bias
                        for lam in lambdas:
                            gap2 = max(1e-12, 0.04 * (shifted_fold - lam))
                            gap = math.sqrt(gap2)
                            handle.write(
                                f"{lam:.17g},{gap:.17g},1e-9,3,4,10,10\n"
                            )
                    runs.append({
                        "n": n,
                        "label": label,
                        "step": step,
                        "csv": {
                            "path": path.name,
                            "sha256": sha256(path),
                        },
                    })

            manifest = temp / "manifest.json"
            manifest.write_text(
                json.dumps({
                    "schema_version": 1,
                    "artifact": "C1_BRATU_FOLD_RUNS",
                    "status": "CALIBRATION_ONLY",
                    "campaign": "C1",
                    "source_commit": "a" * 40,
                    "independent_of_compared_solver": True,
                    "fixed": {
                        "spectral_stop": 5e-2,
                        "eigen_rtol": 1e-7,
                        "threads": 1,
                    },
                    "runs": runs,
                }, indent=2) + "\n",
                encoding="utf-8",
            )
            output = temp / "folds.json"
            frozen = freeze(manifest, output)
            if set(frozen["folds"]) != {"31", "63", "127"}:
                raise RuntimeError("Self-test no produjo los tres folds.")
        finally:
            os.chdir(old_cwd)

    print("C1_BRATU_FOLD_FREEZE_SELF_TEST=PASS")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("run_manifest_json", nargs="?")
    parser.add_argument("output_json", nargs="?")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()

    if args.self_test:
        return self_test()
    if not args.run_manifest_json or not args.output_json:
        raise SystemExit("Uso: freeze_bratu_folds.py RUN_MANIFEST OUTPUT_JSON")

    frozen = freeze(
        Path(args.run_manifest_json).resolve(),
        Path(args.output_json).resolve(),
    )
    print("C1_BRATU_FOLD_FREEZE=PASS")
    for n in SIZES:
        item = frozen["folds"][str(n)]
        print(
            f"n={n} lambda_fold_h={item['lambda_fold_h']:.10f} "
            f"rel_diff={item['relative_difference']:.6g} "
            f"R2={item['fine_r2']:.6g}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
