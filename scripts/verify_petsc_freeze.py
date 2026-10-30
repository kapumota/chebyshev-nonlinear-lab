#!/usr/bin/env python3
"""Verifica que el freeze PETSc de C1 no haya cambiado después de calibración."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path


FAMILIES = {"B1", "B5", "B3a", "Bext"}
SELECTION_RULE = [
    "maximize_success_rate",
    "minimize_median_W_eq_on_success",
    "minimize_median_t_total_on_success",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def check_file(record: dict, label: str) -> None:
    path = Path(record.get("path", "")).resolve()
    expected = record.get("sha256")
    if not path.is_file():
        raise SystemExit(f"Falta {label}: {path}")
    observed = sha256(path)
    if observed != expected:
        raise SystemExit(
            f"Hash inválido para {label}: esperado {expected}, observado {observed}."
        )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("freeze_manifest_json")
    parser.add_argument("--expected-source-commit")
    args = parser.parse_args()

    freeze_path = Path(args.freeze_manifest_json).resolve()
    freeze = json.loads(freeze_path.read_text(encoding="utf-8"))
    if freeze.get("artifact") != "C1_PETSC_CALIBRATION_FREEZE":
        raise SystemExit("Artefacto de freeze PETSc incorrecto.")
    if freeze.get("status") != "FROZEN_CALIBRATION_ONLY":
        raise SystemExit("El freeze PETSc no tiene estado FROZEN_CALIBRATION_ONLY.")
    if freeze.get("selection_rule") != SELECTION_RULE:
        raise SystemExit("La regla de selección PETSc cambió.")
    if args.expected_source_commit and freeze.get("source_commit") != args.expected_source_commit:
        raise SystemExit("source_commit no coincide con el esperado por el gate.")

    inputs = freeze.get("inputs", {})
    for label in (
        "plan", "run_manifest", "raw_results", "results_with_W_eq", "primitive_costs"
    ):
        check_file(inputs.get(label, {}), label)
    for index, record in enumerate(inputs.get("provenance", [])):
        check_file(record, f"provenance[{index}]")

    outputs = freeze.get("outputs", {})
    check_file(outputs.get("ranking", {}), "ranking")
    check_file(outputs.get("frozen_config", {}), "frozen_config")
    check_file(freeze.get("generator", {}), "generator")

    config_path = Path(outputs["frozen_config"]["path"]).resolve()
    config = json.loads(config_path.read_text(encoding="utf-8"))
    if config.get("artifact") != "C1_PETSC_FROZEN_CONFIG":
        raise SystemExit("Configuración PETSc congelada inválida.")
    if config.get("status") != "FROZEN_CALIBRATION_ONLY":
        raise SystemExit("Configuración PETSc no está congelada.")
    if config.get("selection_rule") != SELECTION_RULE:
        raise SystemExit("Configuración PETSc usa otra regla de selección.")
    if set(config.get("families", {})) != FAMILIES:
        raise SystemExit("La configuración PETSc no contiene las cuatro familias.")

    ranking_path = Path(outputs["ranking"]["path"]).resolve()
    with ranking_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    selected = {family: 0 for family in FAMILIES}
    for row in rows:
        family = row.get("family")
        if family not in FAMILIES:
            raise SystemExit("Ranking PETSc contiene una familia inesperada.")
        if row.get("selected") == "1":
            selected[family] += 1
    if any(value != 1 for value in selected.values()):
        raise SystemExit("El ranking debe tener exactamente un candidato seleccionado por familia.")

    print("PETSC_FREEZE_VERIFY=PASS")
    print(f"source_commit={freeze['source_commit']}")
    print(f"freeze_manifest_sha256={sha256(freeze_path)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
