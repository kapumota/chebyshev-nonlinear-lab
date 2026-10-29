#!/usr/bin/env python3
"""Selecciona hiperparámetros de C1 con prioridad success -> W_eq -> t_total.

Entrada: CSV de calibración con W_eq ya calculado.
Las repeticiones de timing se condensan primero por instancia científica.
Las tolerancias y presupuestos de parada se mantienen fijos, no se tunean.
Todos los métodos reciben el mismo conjunto de instancias, repeticiones y número
de configuraciones candidatas.
"""

import argparse
import csv
import math
import statistics
from collections import defaultdict


CONFIG_FIELDS = (
    "eta_init", "eta_max", "gamma", "alpha",
    "eta2_factor", "eta2_max", "restart",
)

BUDGET_FIELDS = (
    "rtol", "max_nl", "maxlin",
)

INSTANCE_FIELDS = (
    "problem", "family", "n", "N", "param", "instance_id", "precond",
    "start_kind", "start_id", "distance", "rho_x", "seed", "threads",
)

REQUIRED_FIELDS = set(CONFIG_FIELDS) | set(BUDGET_FIELDS) | set(INSTANCE_FIELDS) | {
    "phase", "method", "success", "t_total", "W_eq", "repetition", "commit",
}


def tuple_from(row, fields):
    return tuple(row[field] for field in fields)


def finite_positive(value):
    return math.isfinite(value) and value > 0.0


def ranking_key(candidate):
    return (
        -candidate["successes"],
        candidate["median_W_eq"],
        candidate["median_t_total"],
        candidate["config"],
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_csv")
    parser.add_argument("output_csv")
    parser.add_argument(
        "--allow-dirty",
        action="store_true",
        help="Permite commits marcados como dirty sólo para pruebas de infraestructura.",
    )
    args = parser.parse_args()

    with open(args.input_csv, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        missing = REQUIRED_FIELDS - set(reader.fieldnames or [])
        if missing:
            raise SystemExit(
                "Faltan columnas requeridas: " + ",".join(sorted(missing))
            )
        rows = list(reader)

    if not rows:
        raise SystemExit("CSV vacío.")
    if any(row["phase"] != "calibration" for row in rows):
        raise SystemExit("select_tuning.py sólo acepta phase=calibration.")

    commits = {row["commit"] for row in rows}
    if len(commits) != 1:
        raise SystemExit("El tuning no puede mezclar commits distintos.")
    commit = next(iter(commits))
    if "dirty" in commit.lower() and not args.allow_dirty:
        raise SystemExit(
            "El tuning final requiere un commit limpio. Use --allow-dirty sólo para smoke tests."
        )

    budgets = {tuple_from(row, BUDGET_FIELDS) for row in rows}
    if len(budgets) != 1:
        raise SystemExit(
            "rtol, max_nl y maxlin son presupuesto fijo y no pueden mezclarse en tuning."
        )
    budget = next(iter(budgets))

    grouped = defaultdict(lambda: defaultdict(list))
    for row in rows:
        method = row["method"]
        config = tuple_from(row, CONFIG_FIELDS)
        instance = tuple_from(row, INSTANCE_FIELDS)
        grouped[(method, config)][instance].append(row)

    by_method = defaultdict(list)
    for (method, config), instances in grouped.items():
        summarized = []
        coverage = {}

        for instance, reps in instances.items():
            repetition_ids = [int(row["repetition"]) for row in reps]
            if len(set(repetition_ids)) != len(repetition_ids):
                raise SystemExit(
                    f"Repetición duplicada: {method} {config} {instance}"
                )
            coverage[instance] = frozenset(repetition_ids)

            success_values = {int(row["success"]) for row in reps}
            if not success_values <= {0, 1} or len(success_values) != 1:
                raise SystemExit(
                    f"Éxito inconsistente entre repeticiones: "
                    f"{method} {config} {instance}"
                )
            success = success_values.pop()

            if success:
                weq_values = [float(row["W_eq"]) for row in reps]
                total_values = [float(row["t_total"]) for row in reps]
                if not all(finite_positive(value) for value in weq_values):
                    raise SystemExit(
                        f"W_eq inválido en {method} {config} {instance}."
                    )
                if not all(finite_positive(value) for value in total_values):
                    raise SystemExit(
                        f"t_total inválido en {method} {config} {instance}."
                    )
                weq = statistics.median(weq_values)
                total = statistics.median(total_values)
            else:
                weq = float("inf")
                total = float("inf")

            summarized.append((instance, success, weq, total))

        success_count = sum(item[1] for item in summarized)
        success_rate = success_count / len(summarized)
        successful = [item for item in summarized if item[1]]
        median_weq = (
            statistics.median(item[2] for item in successful)
            if successful else float("inf")
        )
        median_total = (
            statistics.median(item[3] for item in successful)
            if successful else float("inf")
        )

        by_method[method].append({
            "config": config,
            "coverage": coverage,
            "instances": len(instances),
            "successes": success_count,
            "success_rate": success_rate,
            "median_W_eq": median_weq,
            "median_t_total": median_total,
        })

    methods = sorted(by_method)
    if not methods:
        raise SystemExit("No hay métodos para seleccionar.")

    candidate_counts = {len(by_method[method]) for method in methods}
    if len(candidate_counts) != 1:
        raise SystemExit(
            "Los métodos no tienen el mismo número de configuraciones candidatas."
        )

    reference_coverage = None
    for method in methods:
        candidates = by_method[method]
        method_coverage = candidates[0]["coverage"]
        for candidate in candidates[1:]:
            if candidate["coverage"] != method_coverage:
                raise SystemExit(
                    f"Instancias o repeticiones desiguales entre configuraciones de {method}."
                )
        if reference_coverage is None:
            reference_coverage = method_coverage
        elif method_coverage != reference_coverage:
            raise SystemExit(
                "Los métodos no cubren exactamente las mismas instancias y repeticiones."
            )

    output_fields = [
        "method", "rank", "selected", *CONFIG_FIELDS, *BUDGET_FIELDS,
        "instances", "successes", "success_rate", "median_W_eq",
        "median_t_total", "commit",
    ]

    with open(args.output_csv, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=output_fields, lineterminator="\n"
        )
        writer.writeheader()

        for method in methods:
            ranked = sorted(by_method[method], key=ranking_key)
            for rank, candidate in enumerate(ranked, start=1):
                row = {
                    "method": method,
                    "rank": rank,
                    "selected": 1 if rank == 1 else 0,
                    **dict(zip(CONFIG_FIELDS, candidate["config"])),
                    **dict(zip(BUDGET_FIELDS, budget)),
                    "instances": candidate["instances"],
                    "successes": candidate["successes"],
                    "success_rate": f"{candidate['success_rate']:.17g}",
                    "median_W_eq": f"{candidate['median_W_eq']:.17g}",
                    "median_t_total": f"{candidate['median_t_total']:.17g}",
                    "commit": commit,
                }
                writer.writerow(row)

    print(
        "Ranking de calibración: success -> W_eq -> t_total. "
        "La fila selected=1 de cada método se versiona antes de confirmación."
    )


if __name__ == "__main__":
    main()
