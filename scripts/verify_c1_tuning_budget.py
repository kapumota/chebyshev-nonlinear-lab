#!/usr/bin/env python3
"""Verifica el presupuesto común de tuning de C1 contra la malla PETSc."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

FAMILIES = ("B1", "B5", "B3a", "Bext")


def load(path: Path) -> dict:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise SystemExit(f"JSON inválido: {path}: {exc}") from exc


def same_float(a, b) -> bool:
    return math.isclose(float(a), float(b), rel_tol=0.0, abs_tol=0.0)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--budget",
        default="experiments/c1/C1_TUNING_BUDGET_V1.json",
    )
    parser.add_argument(
        "--grid",
        default="experiments/c1/C1_INSTANCE_GRID_V1.json",
    )
    args = parser.parse_args()

    budget = load(Path(args.budget))
    grid = load(Path(args.grid))

    if budget.get("schema_version") != 1:
        raise SystemExit("schema_version inválido en presupuesto.")
    if budget.get("artifact") != "C1_TUNING_BUDGET":
        raise SystemExit("artifact inválido en presupuesto.")
    if budget.get("status") != "FROZEN_BEFORE_CALIBRATION":
        raise SystemExit("El presupuesto no está congelado antes de calibración.")
    if budget.get("campaign") != "C1":
        raise SystemExit("El presupuesto no pertenece a C1.")

    count = int(budget.get("common_candidate_count", 0))
    if count != 4:
        raise SystemExit("C1 fija exactamente cuatro candidatos competitivos.")

    if budget.get("start_policy") != "natural_only":
        raise SystemExit("La slice de tuning debe usar solamente start natural.")

    fixed = budget.get("fixed_budget", {})
    expected_fixed = {
        "rtol": 1e-8,
        "max_nl": 50,
        "maxlin": 200,
        "warmup": 1,
        "reps": 7,
        "threads": 1,
    }
    for key, expected in expected_fixed.items():
        observed = fixed.get(key)
        if isinstance(expected, float):
            if observed is None or not same_float(observed, expected):
                raise SystemExit(f"{key} no coincide con presupuesto común.")
        elif observed != expected:
            raise SystemExit(f"{key} no coincide con presupuesto común.")

    internal = budget.get("internal_candidates", [])
    if len(internal) != count:
        raise SystemExit("El tuning interno no tiene cuatro candidatos.")

    ids = [str(x.get("id", "")) for x in internal]
    if len(set(ids)) != count or any(not x for x in ids):
        raise SystemExit("IDs internos vacíos o duplicados.")

    observed_internal = {
        (float(x["eta_init"]), int(x["restart"]))
        for x in internal
    }
    expected_internal = {
        (0.25, 30), (0.25, 50), (0.5, 30), (0.5, 50)
    }
    if observed_internal != expected_internal:
        raise SystemExit("El grid interno 2x2 cambió.")

    family_pc = budget.get("family_preconditioner", {})
    if set(family_pc) != set(FAMILIES):
        raise SystemExit("Faltan precondicionadores de alguna familia.")

    petsc = grid.get("petsc_tuning", {})
    if petsc.get("start_policy") != "natural_only":
        raise SystemExit("PETSc no usa la misma slice natural_only.")

    petsc_candidates = petsc.get("candidates", [])
    if len(petsc_candidates) != count:
        raise SystemExit("PETSc no tiene el presupuesto común de cuatro candidatos.")

    observed_petsc = {
        (float(x["ksp_rtol"]), int(x["restart"]))
        for x in petsc_candidates
    }
    expected_petsc = {
        (1e-4, 30), (1e-4, 50), (1e-5, 30), (1e-5, 50)
    }
    if observed_petsc != expected_petsc:
        raise SystemExit("El grid PETSc 2x2 cambió.")

    if petsc.get("fixed_preconditioner") != family_pc:
        raise SystemExit("Los precondicionadores PETSc no coinciden con el freeze común.")

    petsc_fixed = petsc.get("fixed", {})
    for key, expected in expected_fixed.items():
        observed = petsc_fixed.get(key)
        if isinstance(expected, float):
            if observed is None or not same_float(observed, expected):
                raise SystemExit(f"PETSc {key} no coincide con presupuesto común.")
        elif observed != expected:
            raise SystemExit(f"PETSc {key} no coincide con presupuesto común.")

    methods = budget.get("internal_methods", [])
    if methods != ["nk", "cheb", "rec"]:
        raise SystemExit("Los métodos internos competitivos cambiaron.")

    if budget.get("selection") != ["success", "W_eq", "t_total"]:
        raise SystemExit("La regla de selección cambió.")

    print("C1_TUNING_BUDGET=PASS candidates=4 start=natural_only")
    print("internal_grid=2x2 eta_init,restart")
    print("petsc_grid=2x2 ksp_rtol,restart")
    print("preconditioners=B1:jacobi,B5:none,B3a:mg,Bext:jacobi")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
