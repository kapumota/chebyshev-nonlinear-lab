#!/usr/bin/env python3
"""Congela costos primitivos a partir de CSV de calibración de C1.

Sólo acepta filas phase=calibration. Agrupa por familia, tamaño,
precondicionador y número de hilos, sin usar param como clave, porque las
mallas de calibración y confirmación son deliberadamente diferentes.
"""

import argparse
import csv
import math
import statistics
from collections import defaultdict


KEY_FIELDS = (
    "problem", "family", "n", "N", "precond", "threads"
)

REQUIRED_COLUMNS = set(KEY_FIELDS) | {
    "phase", "method", "n_F", "n_JVP", "n_D2", "n_prec",
    "t_F", "t_JVP", "t_D2", "t_prec",
}


def finite_positive(value):
    return math.isfinite(value) and value > 0.0


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_csv")
    parser.add_argument("output_csv")
    parser.add_argument("--min-samples", type=int, default=5)
    args = parser.parse_args()

    if args.min_samples < 1:
        raise SystemExit("--min-samples debe ser al menos 1.")

    groups = defaultdict(lambda: defaultdict(list))
    methods = defaultdict(set)
    observed_prec_calls = defaultdict(bool)

    with open(args.input_csv, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        missing = REQUIRED_COLUMNS - set(reader.fieldnames or [])
        if missing:
            raise SystemExit(
                "Faltan columnas requeridas: " + ",".join(sorted(missing))
            )
        rows = list(reader)

    if not rows:
        raise SystemExit("CSV de calibración vacío.")

    pairs = (
        ("c_F", "t_F", "n_F"),
        ("c_JVP", "t_JVP", "n_JVP"),
        ("c_D2", "t_D2", "n_D2"),
        ("c_prec", "t_prec", "n_prec"),
    )

    for row in rows:
        if row["phase"] != "calibration":
            raise SystemExit(
                "freeze_primitive_costs.py sólo acepta phase=calibration."
            )

        key = tuple(row[field] for field in KEY_FIELDS)
        methods[key].add(row["method"])

        for name, time_field, count_field in pairs:
            count = int(row[count_field])
            elapsed = float(row[time_field])
            if count < 0 or not math.isfinite(elapsed) or elapsed < 0.0:
                raise SystemExit(f"Conteo o tiempo inválido en {key} para {name}.")
            if name == "c_prec" and count > 0:
                observed_prec_calls[key] = True
            if count > 0 and finite_positive(elapsed):
                groups[key][name].append(elapsed / count)

    all_keys = sorted(set(groups) | set(methods))
    fieldnames = list(KEY_FIELDS) + [
        "c_F", "c_JVP", "c_D2", "c_prec",
        "samples_F", "samples_JVP", "samples_D2", "samples_prec",
        "methods_calibration",
    ]

    with open(args.output_csv, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()

        for key in all_keys:
            data = groups[key]
            precond = key[KEY_FIELDS.index("precond")]
            method_set = methods[key]

            for name in ("c_F", "c_JVP"):
                if len(data[name]) < args.min_samples:
                    raise SystemExit(
                        f"Muestras insuficientes para {name} en {key}: "
                        f"{len(data[name])} < {args.min_samples}"
                    )

            uses_d2 = bool(method_set & {"cheb", "rec"})
            if uses_d2 and len(data["c_D2"]) < args.min_samples:
                raise SystemExit(
                    f"Muestras insuficientes para c_D2 en {key}: "
                    f"{len(data['c_D2'])} < {args.min_samples}"
                )
            if not uses_d2 and data["c_D2"] and len(data["c_D2"]) < args.min_samples:
                raise SystemExit(
                    f"Muestras insuficientes para c_D2 en {key}: "
                    f"{len(data['c_D2'])} < {args.min_samples}"
                )

            requires_prec = precond != "none" or observed_prec_calls[key]
            if requires_prec and len(data["c_prec"]) < args.min_samples:
                raise SystemExit(
                    f"Muestras insuficientes para c_prec en {key}: "
                    f"{len(data['c_prec'])} < {args.min_samples}"
                )

            row = dict(zip(KEY_FIELDS, key))
            for name in ("c_F", "c_JVP", "c_D2", "c_prec"):
                values = data[name]
                row[name] = (
                    f"{statistics.median(values):.17g}"
                    if values else "0"
                )
                row["samples_" + name[2:]] = str(len(values))
            row["methods_calibration"] = ";".join(sorted(method_set))
            writer.writerow(row)

    print(
        "Costos congelados desde calibración. "
        "Revise, audite por param y versione el CSV antes de confirmación."
    )


if __name__ == "__main__":
    main()
