#!/usr/bin/env python3
"""Valida el manifest de un baseline externo antes del freeze de C1."""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

REQUIRED = [
    "schema_version",
    "solver",
    "solver_version",
    "nonlinear_method",
    "linear_method",
    "preconditioner",
    "matrix_mode",
    "options",
    "repository_commit",
    "executable_sha256",
    "build",
    "hardware",
    "command",
    "problems",
    "created_utc",
]

PLACEHOLDERS = {"", "POR_FIJAR", "TODO", "TBD", "UNKNOWN"}


def fail(message: str) -> None:
    raise ValueError(message)


def validate_text(value: object, field: str) -> str:
    if not isinstance(value, str) or value.strip() in PLACEHOLDERS:
        fail(f"{field} debe contener un valor congelado, no un marcador.")
    return value.strip()


def validate_manifest(data: object) -> None:
    if not isinstance(data, dict):
        fail("El manifest debe ser un objeto JSON.")
    missing = [key for key in REQUIRED if key not in data]
    if missing:
        fail("Faltan campos obligatorios: " + ", ".join(missing))

    if data["schema_version"] != 1:
        fail("schema_version debe ser 1.")
    if validate_text(data["solver"], "solver") != "PETSc":
        fail("El baseline requerido en este gate debe ser PETSc.")
    validate_text(data["solver_version"], "solver_version")
    if validate_text(data["nonlinear_method"], "nonlinear_method") != "SNESNEWTONLS":
        fail("nonlinear_method debe ser SNESNEWTONLS.")
    validate_text(data["linear_method"], "linear_method")
    validate_text(data["preconditioner"], "preconditioner")
    validate_text(data["matrix_mode"], "matrix_mode")
    validate_text(data["repository_commit"], "repository_commit")
    validate_text(data["command"], "command")
    validate_text(data["created_utc"], "created_utc")

    digest = validate_text(data["executable_sha256"], "executable_sha256")
    if not re.fullmatch(r"[0-9a-fA-F]{64}", digest):
        fail("executable_sha256 debe tener 64 dígitos hexadecimales.")

    if not isinstance(data["options"], list) or not data["options"]:
        fail("options debe ser una lista no vacía con las opciones efectivas.")
    if not all(isinstance(item, str) and item.strip() for item in data["options"]):
        fail("Cada entrada de options debe ser texto no vacío.")

    build = data["build"]
    if not isinstance(build, dict):
        fail("build debe ser un objeto.")
    validate_text(build.get("compiler"), "build.compiler")
    validate_text(build.get("flags"), "build.flags")

    hardware = data["hardware"]
    if not isinstance(hardware, dict):
        fail("hardware debe ser un objeto.")
    validate_text(hardware.get("cpu"), "hardware.cpu")
    threads = hardware.get("threads")
    if not isinstance(threads, int) or threads <= 0:
        fail("hardware.threads debe ser un entero positivo.")

    problems = data["problems"]
    if not isinstance(problems, list) or not problems:
        fail("problems debe registrar al menos una familia evaluada.")
    if not all(isinstance(item, str) and item.strip() for item in problems):
        fail("Cada entrada de problems debe ser texto no vacío.")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("manifest", type=Path)
    args = parser.parse_args()

    try:
        with args.manifest.open(encoding="utf-8") as handle:
            data = json.load(handle)
        validate_manifest(data)
    except (OSError, json.JSONDecodeError, ValueError) as exc:
        print(f"MANIFEST_EXTERNAL_BASELINE=FAIL: {exc}", file=sys.stderr)
        return 1

    print("MANIFEST_EXTERNAL_BASELINE=PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
