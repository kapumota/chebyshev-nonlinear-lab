#!/usr/bin/env python3
"""Selecciona y congela PETSc usando sólo calibración de C1.

La entrada de resultados debe contener W_eq producido con la tabla de costos
congelada. La selección replica la semántica de select_tuning.py:
success -> W_eq -> t_total.
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
import hashlib
import json
import math
from pathlib import Path
import statistics


FAMILIES = ("B1", "B5", "B3a", "Bext")
SELECTION_RULE = (
    "maximize_success_rate",
    "minimize_median_W_eq_on_success",
    "minimize_median_t_total_on_success",
)
INSTANCE_FIELDS = (
    "problem", "family", "n", "N", "param", "precond",
    "start_kind", "start_id", "distance", "rho_x", "seed", "threads",
    "calibration_instance_id",
)
REQUIRED_FIELDS = set(INSTANCE_FIELDS) | {
    "phase", "method", "success", "t_total", "W_eq", "repetition", "commit",
    "petsc_candidate_id", "petsc_ksp_rtol", "petsc_restart",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def finite_positive(value: float) -> bool:
    return math.isfinite(value) and value > 0.0


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def tuple_from(row: dict[str, str], fields: tuple[str, ...]) -> tuple[str, ...]:
    return tuple(row[field] for field in fields)


def validate_plan(plan: dict) -> None:
    if plan.get("schema_version") != 1:
        raise SystemExit("Plan PETSc con schema_version inválido.")
    if plan.get("campaign") != "C1" or plan.get("phase") != "calibration":
        raise SystemExit("El freeze sólo acepta el plan C1 de calibration.")
    fixed = plan.get("fixed", {})
    expected = {
        "snes_type": "newtonls",
        "linesearch_type": "bt",
        "ksp_type": "gmres",
        "matrix_mode": "MatShell analytic Jv",
    }
    for key, value in expected.items():
        if fixed.get(key) != value:
            raise SystemExit(f"fixed.{key} no coincide con el mecanismo PETSc congelable.")
    if int(fixed.get("threads", 0)) != 1:
        raise SystemExit("PETSc C1 debe estar congelado a un hilo.")
    if tuple(sorted(plan.get("families", {}))) != tuple(sorted(FAMILIES)):
        raise SystemExit("El plan no contiene exactamente las cuatro familias PETSc.")


def validate_run_manifest(run: dict, plan_path: Path, plan: dict, raw_path: Path) -> None:
    if run.get("artifact") != "C1_PETSC_CALIBRATION_RUN":
        raise SystemExit("Manifest de ejecución PETSc incorrecto.")
    if run.get("phase") != "calibration":
        raise SystemExit("El manifest de PETSc no corresponde a calibration.")
    if run.get("source_commit") != plan.get("source_commit"):
        raise SystemExit("source_commit difiere entre plan y ejecución PETSc.")
    if run.get("plan", {}).get("sha256") != sha256(plan_path):
        raise SystemExit("El plan PETSc cambió después de la ejecución.")
    if run.get("output_csv", {}).get("sha256") != sha256(raw_path):
        raise SystemExit("El CSV crudo PETSc cambió después de la ejecución.")


def summarize(rows: list[dict[str, str]], family: str, candidate_id: str) -> dict:
    instances: dict[tuple[str, ...], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        instances[tuple_from(row, INSTANCE_FIELDS)].append(row)

    coverage: dict[tuple[str, ...], frozenset[int]] = {}
    summarized = []
    for instance, reps in instances.items():
        repetition_ids = [int(row["repetition"]) for row in reps]
        if len(set(repetition_ids)) != len(repetition_ids):
            raise SystemExit(
                f"Repetición duplicada en {family}/{candidate_id}/{instance}."
            )
        coverage[instance] = frozenset(repetition_ids)

        success_values = {int(row["success"]) for row in reps}
        if not success_values <= {0, 1} or len(success_values) != 1:
            raise SystemExit(
                f"Éxito inconsistente entre repeticiones en {family}/{candidate_id}."
            )
        success = success_values.pop()
        if success:
            weq_values = [float(row["W_eq"]) for row in reps]
            total_values = [float(row["t_total"]) for row in reps]
            if not all(finite_positive(value) for value in weq_values):
                raise SystemExit(f"W_eq inválido en {family}/{candidate_id}.")
            if not all(finite_positive(value) for value in total_values):
                raise SystemExit(f"t_total inválido en {family}/{candidate_id}.")
            weq = statistics.median(weq_values)
            total = statistics.median(total_values)
        else:
            weq = math.inf
            total = math.inf
        summarized.append((success, weq, total))

    successes = sum(item[0] for item in summarized)
    successful = [item for item in summarized if item[0]]
    return {
        "coverage": coverage,
        "instances": len(summarized),
        "successes": successes,
        "success_rate": successes / len(summarized),
        "median_W_eq": (
            statistics.median(item[1] for item in successful)
            if successful else math.inf
        ),
        "median_t_total": (
            statistics.median(item[2] for item in successful)
            if successful else math.inf
        ),
    }


def ranking_key(candidate: dict) -> tuple:
    return (
        -candidate["successes"],
        candidate["median_W_eq"],
        candidate["median_t_total"],
        candidate["candidate_id"],
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("plan_json")
    parser.add_argument("run_manifest_json")
    parser.add_argument("results_weq_csv")
    parser.add_argument("primitive_costs_csv")
    parser.add_argument("ranking_csv")
    parser.add_argument("frozen_config_json")
    parser.add_argument("freeze_manifest_json")
    parser.add_argument(
        "--provenance",
        action="append",
        default=[],
        help="Archivo adicional cuya huella debe quedar incluida en el freeze.",
    )
    args = parser.parse_args()

    plan_path = Path(args.plan_json).resolve()
    run_manifest_path = Path(args.run_manifest_json).resolve()
    results_path = Path(args.results_weq_csv).resolve()
    costs_path = Path(args.primitive_costs_csv).resolve()
    ranking_path = Path(args.ranking_csv).resolve()
    config_path = Path(args.frozen_config_json).resolve()
    freeze_path = Path(args.freeze_manifest_json).resolve()

    for path in (plan_path, run_manifest_path, results_path, costs_path):
        if not path.is_file():
            raise SystemExit(f"No existe archivo requerido: {path}")

    plan = load_json(plan_path)
    validate_plan(plan)
    run_manifest = load_json(run_manifest_path)
    raw_path = Path(run_manifest.get("output_csv", {}).get("path", "")).resolve()
    if not raw_path.is_file():
        raise SystemExit("No se encuentra el CSV crudo declarado por el run manifest.")
    validate_run_manifest(run_manifest, plan_path, plan, raw_path)

    with results_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        missing = REQUIRED_FIELDS - set(reader.fieldnames or [])
        if missing:
            raise SystemExit(
                "Faltan columnas requeridas en PETSc+W_eq: " + ",".join(sorted(missing))
            )
        rows = list(reader)
    if not rows:
        raise SystemExit("CSV PETSc+W_eq vacío.")

    source_commit = plan["source_commit"]
    for row in rows:
        if row["phase"] != "calibration":
            raise SystemExit("El freeze PETSc rechaza cualquier fila no calibration.")
        if row["method"] != "petsc":
            raise SystemExit("El freeze PETSc sólo acepta method=petsc.")
        if row["family"] not in FAMILIES:
            raise SystemExit("Familia PETSc inesperada.")
        if row["commit"] != source_commit or "dirty" in row["commit"].lower():
            raise SystemExit("El freeze PETSc exige un único commit limpio del plan.")
        if int(row["threads"]) != 1:
            raise SystemExit("El freeze PETSc exige threads=1.")

    expected_count = int(plan["expected_candidate_count"])
    by_family_candidate: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        by_family_candidate[(row["family"], row["petsc_candidate_id"])].append(row)

    ranking_rows = []
    frozen_families = {}
    input_result_hash = sha256(results_path)

    for family in FAMILIES:
        block = plan["families"][family]
        plan_candidates = {str(item["id"]): item for item in block["candidates"]}
        if len(plan_candidates) != expected_count:
            raise SystemExit(f"El plan cambió el presupuesto de candidatos en {family}.")

        observed_ids = {
            candidate_id
            for current_family, candidate_id in by_family_candidate
            if current_family == family
        }
        if observed_ids != set(plan_candidates):
            raise SystemExit(
                f"Candidatos observados en {family} no coinciden con el plan."
            )

        candidates = []
        reference_coverage = None
        for candidate_id, candidate in plan_candidates.items():
            subset = by_family_candidate[(family, candidate_id)]
            expected_precond = str(candidate["precond"])
            expected_ksp_rtol = float(candidate["ksp_rtol"])
            expected_restart = int(candidate["restart"])
            for row in subset:
                if row["precond"] != expected_precond:
                    raise SystemExit(f"precond no coincide con plan en {family}/{candidate_id}.")
                if not math.isclose(
                    float(row["petsc_ksp_rtol"]), expected_ksp_rtol,
                    rel_tol=0.0, abs_tol=0.0,
                ):
                    raise SystemExit(f"ksp_rtol no coincide con plan en {family}/{candidate_id}.")
                if int(row["petsc_restart"]) != expected_restart:
                    raise SystemExit(f"restart no coincide con plan en {family}/{candidate_id}.")

            stats = summarize(subset, family, candidate_id)
            if reference_coverage is None:
                reference_coverage = stats["coverage"]
            elif stats["coverage"] != reference_coverage:
                raise SystemExit(
                    f"Instancias o repeticiones desiguales entre candidatos de {family}."
                )
            candidates.append({
                "candidate_id": candidate_id,
                "precond": expected_precond,
                "ksp_rtol": expected_ksp_rtol,
                "restart": expected_restart,
                **stats,
            })

        ranked = sorted(candidates, key=ranking_key)
        for rank, candidate in enumerate(ranked, start=1):
            ranking_rows.append({
                "family": family,
                "rank": rank,
                "selected": 1 if rank == 1 else 0,
                "candidate_id": candidate["candidate_id"],
                "precond": candidate["precond"],
                "ksp_rtol": f"{candidate['ksp_rtol']:.17g}",
                "restart": candidate["restart"],
                "instances": candidate["instances"],
                "successes": candidate["successes"],
                "success_rate": f"{candidate['success_rate']:.17g}",
                "median_W_eq": (
                    f"{candidate['median_W_eq']:.17g}"
                    if math.isfinite(candidate["median_W_eq"]) else "inf"
                ),
                "median_t_total": (
                    f"{candidate['median_t_total']:.17g}"
                    if math.isfinite(candidate["median_t_total"]) else "inf"
                ),
                "source_commit": source_commit,
            })

        selected = ranked[0]
        frozen_families[family] = {
            "selected_candidate_id": selected["candidate_id"],
            "precond": selected["precond"],
            "ksp_rtol": selected["ksp_rtol"],
            "gmres_restart": selected["restart"],
            "success_rate_calibration": selected["success_rate"],
            "median_W_eq_calibration": (
                selected["median_W_eq"]
                if math.isfinite(selected["median_W_eq"]) else None
            ),
            "median_t_total_calibration": (
                selected["median_t_total"]
                if math.isfinite(selected["median_t_total"]) else None
            ),
        }

    ranking_path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "family", "rank", "selected", "candidate_id", "precond", "ksp_rtol",
        "restart", "instances", "successes", "success_rate", "median_W_eq",
        "median_t_total", "source_commit",
    ]
    with ranking_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(ranking_rows)

    frozen_config = {
        "schema_version": 1,
        "artifact": "C1_PETSC_FROZEN_CONFIG",
        "status": "FROZEN_CALIBRATION_ONLY",
        "campaign": "C1",
        "method": "petsc",
        "source_commit": source_commit,
        "selection_rule": list(SELECTION_RULE),
        "fixed": plan["fixed"],
        "families": frozen_families,
    }
    config_path.parent.mkdir(parents=True, exist_ok=True)
    config_path.write_text(
        json.dumps(frozen_config, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    provenance = []
    for item in args.provenance:
        path = Path(item).resolve()
        if not path.is_file():
            raise SystemExit(f"No existe provenance: {path}")
        provenance.append({"path": item, "sha256": sha256(path)})

    generator_path = Path(__file__).resolve()
    freeze_manifest = {
        "schema_version": 1,
        "artifact": "C1_PETSC_CALIBRATION_FREEZE",
        "status": "FROZEN_CALIBRATION_ONLY",
        "campaign": "C1",
        "source_commit": source_commit,
        "selection_rule": list(SELECTION_RULE),
        "inputs": {
            "plan": {"path": args.plan_json, "sha256": sha256(plan_path)},
            "run_manifest": {
                "path": args.run_manifest_json,
                "sha256": sha256(run_manifest_path),
            },
            "raw_results": {
                "path": run_manifest["output_csv"]["path"],
                "sha256": run_manifest["output_csv"]["sha256"],
            },
            "results_with_W_eq": {
                "path": args.results_weq_csv,
                "sha256": input_result_hash,
            },
            "primitive_costs": {
                "path": args.primitive_costs_csv,
                "sha256": sha256(costs_path),
            },
            "provenance": provenance,
        },
        "outputs": {
            "ranking": {"path": args.ranking_csv, "sha256": sha256(ranking_path)},
            "frozen_config": {
                "path": args.frozen_config_json,
                "sha256": sha256(config_path),
            },
        },
        "generator": {
            "path": "scripts/freeze_petsc_calibration.py",
            "sha256": sha256(generator_path),
        },
    }
    freeze_path.parent.mkdir(parents=True, exist_ok=True)
    freeze_path.write_text(
        json.dumps(freeze_manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    print("PETSC_CALIBRATION_FREEZE=PASS")
    for family in FAMILIES:
        selected = frozen_families[family]
        print(
            family,
            selected["selected_candidate_id"],
            f"success={selected['success_rate_calibration']:.6g}",
            f"W_eq={selected['median_W_eq_calibration']}",
            f"t_total={selected['median_t_total_calibration']}",
        )
    print(f"config_sha256={sha256(config_path)}")
    print(f"ranking_sha256={sha256(ranking_path)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
