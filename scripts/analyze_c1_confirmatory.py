#!/usr/bin/env python3
"""Análisis confirmatorio C1 fijado antes de observar resultados."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import random
import statistics
import tempfile
from collections import defaultdict
from pathlib import Path

BOOTSTRAP_REPS = 10000
BOOTSTRAP_LEVEL = 0.95
BOOTSTRAP_MASTER_SEED = 20261029

METHODS = ("nk", "cheb", "rec", "tf_cheb", "petsc")
METRICS = ("W_eq", "t_total")
COMPARISONS = (
    ("cheb", "nk"),
    ("rec", "nk"),
    ("tf_cheb", "nk"),
    ("petsc", "nk"),
    ("rec", "cheb"),
)
INSTANCE_FIELDS = (
    "problem", "family", "n", "N", "param", "instance_id",
    "start_kind", "start_id", "distance", "rho_x", "seed", "threads",
)
REQUIRED = set(INSTANCE_FIELDS) | {
    "phase", "method", "success", "repetition", "commit", "W_eq", "t_total",
}


def file_hash(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def key(row: dict[str, str]) -> tuple[str, ...]:
    return tuple(row[name] for name in INSTANCE_FIELDS)


def percentile(values: list[float], p: float) -> float:
    xs = sorted(values)
    if not xs:
        raise SystemExit("Bootstrap sin observaciones.")
    if len(xs) == 1:
        return xs[0]
    x = p * (len(xs) - 1)
    lo = int(math.floor(x))
    hi = int(math.ceil(x))
    if lo == hi:
        return xs[lo]
    w = x - lo
    return xs[lo] * (1.0 - w) + xs[hi] * w


def seed_for(label: str) -> int:
    payload = f"{BOOTSTRAP_MASTER_SEED}:{label}".encode()
    return int(hashlib.sha256(payload).hexdigest()[:16], 16)


def bootstrap_mean(values: list[float], label: str) -> tuple[float, float]:
    if len(values) < 2:
        raise SystemExit(f"Menos de dos unidades para bootstrap: {label}.")
    rng = random.Random(seed_for(label))
    n = len(values)
    estimates = []
    for _ in range(BOOTSTRAP_REPS):
        estimates.append(statistics.fmean(values[rng.randrange(n)] for _ in range(n)))
    alpha = 1.0 - BOOTSTRAP_LEVEL
    return percentile(estimates, alpha / 2), percentile(estimates, 1 - alpha / 2)


def load(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        missing = REQUIRED - set(reader.fieldnames or [])
        if missing:
            raise SystemExit("Faltan columnas: " + ",".join(sorted(missing)))
        rows = list(reader)
    if not rows:
        raise SystemExit("CSV confirmatorio vacío.")
    return rows


def validate(rows: list[dict[str, str]]) -> str:
    if any(r["phase"] != "confirmation" for r in rows):
        raise SystemExit("Sólo se acepta phase=confirmation.")
    commits = {r["commit"] for r in rows}
    if len(commits) != 1:
        raise SystemExit("La confirmación mezcla commits.")
    commit = next(iter(commits))
    if "dirty" in commit.lower():
        raise SystemExit("Commit dirty no admisible.")

    methods = {r["method"] for r in rows}
    if methods != set(METHODS):
        missing = set(METHODS) - methods
        extra = methods - set(METHODS)
        raise SystemExit(
            "Métodos incompatibles. "
            f"missing={','.join(sorted(missing))} extra={','.join(sorted(extra))}"
        )
    return commit


def collapse(rows: list[dict[str, str]]):
    grouped = defaultdict(list)
    for row in rows:
        grouped[(row["method"], key(row))].append(row)

    result = {method: {} for method in METHODS}
    for (method, instance), reps in grouped.items():
        ids = [int(r["repetition"]) for r in reps]
        if len(ids) != len(set(ids)):
            raise SystemExit(f"Repetición duplicada: {method} {instance}.")
        success_values = {int(r["success"]) for r in reps}
        if not success_values <= {0, 1} or len(success_values) != 1:
            raise SystemExit(f"Success inconsistente: {method} {instance}.")
        success = bool(next(iter(success_values)))

        record = {"success": success, "repetitions": tuple(sorted(ids))}
        for metric in METRICS:
            values = [float(r[metric]) for r in reps]
            if success:
                if not all(math.isfinite(v) and v > 0 for v in values):
                    raise SystemExit(f"{metric} inválido en éxito: {method} {instance}.")
                record[metric] = statistics.median(values)
            else:
                record[metric] = math.inf
        result[method][instance] = record

    reference = set(result[METHODS[0]])
    if not reference:
        raise SystemExit("No hay instancias.")
    for method in METHODS:
        if set(result[method]) != reference:
            raise SystemExit(f"Cobertura desigual para {method}.")
    for instance in reference:
        reps = result[METHODS[0]][instance]["repetitions"]
        for method in METHODS[1:]:
            if result[method][instance]["repetitions"] != reps:
                raise SystemExit(f"Repeticiones desiguales en {instance}.")
    return result


def write_summary(path: Path, data) -> None:
    total = len(data[METHODS[0]])
    with path.open("w", newline="", encoding="utf-8") as f:
        fields = ("method", "instances", "successes", "success_rate")
        w = csv.DictWriter(f, fieldnames=fields, lineterminator="\n")
        w.writeheader()
        for method in METHODS:
            successes = sum(bool(x["success"]) for x in data[method].values())
            w.writerow({
                "method": method,
                "instances": total,
                "successes": successes,
                "success_rate": f"{successes / total:.17g}",
            })


def write_pairs(effects_path: Path, pairs_path: Path, data) -> None:
    effect_fields = (
        "candidate", "baseline", "metric", "instances_total", "common_successes",
        "geomean_ratio", "ci95_low", "ci95_high",
        "candidate_success_rate", "baseline_success_rate",
        "success_diff_pp", "success_diff_ci95_low_pp", "success_diff_ci95_high_pp",
    )
    pair_fields = (
        "candidate", "baseline", *INSTANCE_FIELDS,
        "candidate_success", "baseline_success",
        "candidate_W_eq", "baseline_W_eq", "log_ratio_W_eq",
        "candidate_t_total", "baseline_t_total", "log_ratio_t_total",
    )
    instances = sorted(data[METHODS[0]])

    with effects_path.open("w", newline="", encoding="utf-8") as ef, \
         pairs_path.open("w", newline="", encoding="utf-8") as pf:
        ew = csv.DictWriter(ef, fieldnames=effect_fields, lineterminator="\n")
        pw = csv.DictWriter(pf, fieldnames=pair_fields, lineterminator="\n")
        ew.writeheader()
        pw.writeheader()

        for candidate, baseline in COMPARISONS:
            logs = {metric: [] for metric in METRICS}
            success_diffs = []

            for instance in instances:
                c = data[candidate][instance]
                b = data[baseline][instance]
                cs, bs = bool(c["success"]), bool(b["success"])
                success_diffs.append(float(cs) - float(bs))
                row = {
                    "candidate": candidate,
                    "baseline": baseline,
                    **dict(zip(INSTANCE_FIELDS, instance)),
                    "candidate_success": int(cs),
                    "baseline_success": int(bs),
                }
                for metric in METRICS:
                    cv, bv = float(c[metric]), float(b[metric])
                    lr = ""
                    if cs and bs:
                        value = math.log(cv / bv)
                        logs[metric].append(value)
                        lr = f"{value:.17g}"
                    row[f"candidate_{metric}"] = f"{cv:.17g}" if math.isfinite(cv) else "inf"
                    row[f"baseline_{metric}"] = f"{bv:.17g}" if math.isfinite(bv) else "inf"
                    row[f"log_ratio_{metric}"] = lr
                pw.writerow(row)

            csr = statistics.fmean(float(data[candidate][i]["success"]) for i in instances)
            bsr = statistics.fmean(float(data[baseline][i]["success"]) for i in instances)
            diff_pp = 100.0 * statistics.fmean(success_diffs)
            success_ci = bootstrap_mean(
                [100.0 * x for x in success_diffs],
                f"{candidate}/{baseline}:success",
            )

            for metric in METRICS:
                if len(logs[metric]) < 2:
                    raise SystemExit(
                        f"Pares exitosos insuficientes: {candidate}/{baseline} {metric}."
                    )
                estimate = statistics.fmean(logs[metric])
                low, high = bootstrap_mean(
                    logs[metric], f"{candidate}/{baseline}:{metric}"
                )
                ew.writerow({
                    "candidate": candidate,
                    "baseline": baseline,
                    "metric": metric,
                    "instances_total": len(instances),
                    "common_successes": len(logs[metric]),
                    "geomean_ratio": f"{math.exp(estimate):.17g}",
                    "ci95_low": f"{math.exp(low):.17g}",
                    "ci95_high": f"{math.exp(high):.17g}",
                    "candidate_success_rate": f"{csr:.17g}",
                    "baseline_success_rate": f"{bsr:.17g}",
                    "success_diff_pp": f"{diff_pp:.17g}",
                    "success_diff_ci95_low_pp": f"{success_ci[0]:.17g}",
                    "success_diff_ci95_high_pp": f"{success_ci[1]:.17g}",
                })


def write_profile(path: Path, data) -> None:
    instances = sorted(data[METHODS[0]])
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(
            f, fieldnames=("metric", "method", "tau", "fraction"),
            lineterminator="\n",
        )
        w.writeheader()
        for metric in METRICS:
            ratios = {method: [] for method in METHODS}
            taus = {1.0}
            for instance in instances:
                finite = [
                    float(data[m][instance][metric])
                    for m in METHODS if bool(data[m][instance]["success"])
                ]
                best = min(finite) if finite else math.inf
                for method in METHODS:
                    record = data[method][instance]
                    if bool(record["success"]) and math.isfinite(best):
                        ratio = float(record[metric]) / best
                        if 1.0 - 1e-12 < ratio < 1.0:
                            ratio = 1.0
                        ratios[method].append(ratio)
                        taus.add(ratio)
                    else:
                        ratios[method].append(math.inf)

            for method in METHODS:
                for tau in sorted(taus):
                    fraction = sum(r <= tau for r in ratios[method]) / len(instances)
                    w.writerow({
                        "metric": metric,
                        "method": method,
                        "tau": f"{tau:.17g}",
                        "fraction": f"{fraction:.17g}",
                    })


def analyze(input_path: Path, output_dir: Path) -> None:
    rows = load(input_path)
    commit = validate(rows)
    data = collapse(rows)
    output_dir.mkdir(parents=True, exist_ok=True)

    summary = output_dir / "summary_methods.csv"
    effects = output_dir / "paired_effects.csv"
    pairs = output_dir / "paired_instances.csv"
    profile = output_dir / "performance_profile.csv"
    manifest = output_dir / "analysis_manifest.json"

    write_summary(summary, data)
    write_pairs(effects, pairs, data)
    write_profile(profile, data)

    outputs = (summary, effects, pairs, profile)
    payload = {
        "schema_version": 1,
        "campaign": "C1",
        "phase": "confirmation",
        "status": "PREDEFINED_BEFORE_C1",
        "source_commit": commit,
        "input_sha256": file_hash(input_path),
        "methods": list(METHODS),
        "comparisons": [
            {"candidate": a, "baseline": b} for a, b in COMPARISONS
        ],
        "metrics": list(METRICS),
        "estimator": "exp(mean(paired_log_ratio))",
        "bootstrap": {
            "type": "paired_percentile",
            "reps": BOOTSTRAP_REPS,
            "level": BOOTSTRAP_LEVEL,
            "master_seed": BOOTSTRAP_MASTER_SEED,
            "unit": "scientific_instance",
        },
        "success_effect": "candidate_minus_baseline_percentage_points",
        "performance_profile": {
            "failure_ratio": "inf",
            "denominator": "all_paired_instances",
        },
        "instances": len(data[METHODS[0]]),
        "outputs_sha256": {p.name: file_hash(p) for p in outputs},
    }
    manifest.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(
        "C1_STATISTICAL_ANALYSIS=PASS "
        f"instances={payload['instances']} bootstrap={BOOTSTRAP_REPS}"
    )


def self_test() -> None:
    fields = (
        "phase", "problem", "family", "n", "N", "param", "instance_id",
        "method", "precond", "start_kind", "start_id", "distance", "rho_x",
        "seed", "success", "t_total", "threads", "repetition", "commit", "W_eq",
    )
    with tempfile.TemporaryDirectory(prefix="c1-analysis-") as temp:
        root = Path(temp)
        source = root / "confirmation.csv"
        out = root / "out"
        with source.open("w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=fields, lineterminator="\n")
            w.writeheader()
            for instance in range(4):
                for mi, method in enumerate(METHODS):
                    success = 0 if instance == 3 and method == "petsc" else 1
                    for rep in range(1, 4):
                        w.writerow({
                            "phase": "confirmation",
                            "problem": "contract-test",
                            "family": "Btest",
                            "n": 3,
                            "N": 27,
                            "param": instance + 1,
                            "instance_id": f"case={instance}",
                            "method": method,
                            "precond": "frozen",
                            "start_kind": "natural",
                            "start_id": "natural",
                            "distance": "natural",
                            "rho_x": 0,
                            "seed": instance,
                            "success": success,
                            "t_total": (1 + 0.1 * mi) * (instance + 1),
                            "threads": 1,
                            "repetition": rep,
                            "commit": "0123456789abcdef",
                            "W_eq": (2 + 0.2 * mi) * (instance + 1),
                        })
        analyze(source, out)
        with (out / "summary_methods.csv").open(newline="", encoding="utf-8") as f:
            rows = {r["method"]: r for r in csv.DictReader(f)}
        assert rows["nk"]["successes"] == "4"
        assert rows["petsc"]["successes"] == "3"
        manifest = json.loads((out / "analysis_manifest.json").read_text(encoding="utf-8"))
        assert manifest["bootstrap"]["reps"] == 10000
        assert manifest["performance_profile"]["failure_ratio"] == "inf"
        print("C1_STATISTICAL_ANALYSIS_SELF_TEST=PASS")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input_csv", nargs="?")
    parser.add_argument("output_dir", nargs="?")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return
    if not args.input_csv or not args.output_dir:
        parser.error("input_csv y output_dir son obligatorios.")
    analyze(Path(args.input_csv).resolve(), Path(args.output_dir).resolve())


if __name__ == "__main__":
    main()
