#!/usr/bin/env python3
"""Añade W_eq a un CSV de campaña usando una tabla de costos congelada.

La tabla debe haberse producido exclusivamente con calibración. La unión no
usa param porque calibración y confirmación emplean mallas distintas.
"""

import argparse
import csv
import math


KEY_FIELDS = (
    "problem", "family", "n", "N", "precond", "threads"
)


def key(row):
    return tuple(row[field] for field in KEY_FIELDS)


def finite_nonnegative(value):
    return math.isfinite(value) and value >= 0.0


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("costs_csv")
    parser.add_argument("campaign_csv")
    parser.add_argument("output_csv")
    args = parser.parse_args()

    costs = {}
    with open(args.costs_csv, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            row_key = key(row)
            if row_key in costs:
                raise SystemExit(
                    "Clave duplicada en tabla de costos: " + repr(row_key)
                )
            costs[row_key] = row

    if not costs:
        raise SystemExit("Tabla de costos vacía.")

    with open(args.campaign_csv, newline="", encoding="utf-8") as source:
        reader = csv.DictReader(source)
        rows = list(reader)
        if not rows:
            raise SystemExit("CSV de campaña vacío.")
        fieldnames = list(reader.fieldnames or [])
        if "W_eq" in fieldnames:
            raise SystemExit("El CSV de campaña ya contiene W_eq.")
        fieldnames.append("W_eq")

    with open(args.output_csv, "w", newline="", encoding="utf-8") as target:
        writer = csv.DictWriter(target, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()

        for row in rows:
            frozen = costs.get(key(row))
            if frozen is None:
                raise SystemExit(
                    "No existe costo congelado para la clave: "
                    + repr(key(row))
                )

            c_F = float(frozen["c_F"])
            c_JVP = float(frozen["c_JVP"])
            c_D2 = float(frozen["c_D2"])
            c_prec = float(frozen["c_prec"])

            if not finite_nonnegative(c_F):
                raise SystemExit("c_F inválido en tabla congelada.")
            if not (math.isfinite(c_JVP) and c_JVP > 0.0):
                raise SystemExit("c_JVP inválido en tabla congelada.")
            if not finite_nonnegative(c_D2):
                raise SystemExit("c_D2 inválido en tabla congelada.")
            if not finite_nonnegative(c_prec):
                raise SystemExit("c_prec inválido en tabla congelada.")

            n_F = int(row["n_F"])
            n_JVP = int(row["n_JVP"])
            n_D2 = int(row["n_D2"])
            n_prec = int(row["n_prec"])
            t_setup = float(row["t_setup"])
            if min(n_F, n_JVP, n_D2, n_prec) < 0:
                raise SystemExit("Conteo negativo en CSV de campaña.")
            if not finite_nonnegative(t_setup):
                raise SystemExit("t_setup inválido en CSV de campaña.")
            if n_D2 > 0 and c_D2 <= 0.0:
                raise SystemExit("La fila usa D2 pero c_D2 no está calibrado.")
            if n_prec > 0 and c_prec <= 0.0:
                raise SystemExit(
                    "La fila usa precondicionador pero c_prec no está calibrado."
                )

            numerator = (
                n_F * c_F
                + n_JVP * c_JVP
                + n_D2 * c_D2
                + n_prec * c_prec
                + t_setup
            )
            row["W_eq"] = f"{numerator / c_JVP:.17g}"
            writer.writerow(row)


if __name__ == "__main__":
    main()
