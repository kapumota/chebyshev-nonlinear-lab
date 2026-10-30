#!/usr/bin/env python3
"""Ejecuta el smoke PETSc y genera procedencia verificable del baseline."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shlex
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

BEGIN = "PETSC_BASELINE_BEGIN"
END = "PETSC_BASELINE_END"


def run_text(command: list[str], cwd: Path | None = None) -> str:
    result = subprocess.run(
        command,
        cwd=cwd,
        check=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    return result.stdout.strip()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def git_value(repo: Path, *args: str) -> str:
    return run_text(["git", "-C", str(repo), *args])


def parse_metadata(stdout: str) -> dict[str, str]:
    inside = False
    values: dict[str, str] = {}
    for raw in stdout.splitlines():
        line = raw.strip()
        if line == BEGIN:
            inside = True
            continue
        if line == END:
            inside = False
            break
        if inside and "=" in line:
            key, value = line.split("=", 1)
            values[key.strip()] = value.strip()
    if not values:
        raise ValueError("El ejecutable no emitió metadatos PETSC_BASELINE_BEGIN/END.")
    required = {
        "petsc_version",
        "snes_type",
        "linesearch_type",
        "ksp_type",
        "pc_type",
        "matrix_mode",
        "snes_rtol",
        "snes_max_it",
        "ksp_rtol",
        "ksp_max_it",
        "gmres_restart",
    }
    missing = sorted(required - values.keys())
    if missing:
        raise ValueError("Faltan metadatos efectivos: " + ", ".join(missing))
    return values


def compile_entry(build_dir: Path, executable_source: str) -> tuple[str, str]:
    path = build_dir / "compile_commands.json"
    data = json.loads(path.read_text(encoding="utf-8"))
    for entry in data:
        source = str(entry.get("file", ""))
        if source.endswith(executable_source):
            command = entry.get("command")
            if command:
                parts = shlex.split(command)
                return parts[0], command
            arguments = entry.get("arguments")
            if arguments:
                return arguments[0], shlex.join(arguments)
    raise ValueError(f"No se encontró {executable_source} en {path}.")


def cpu_model() -> str:
    output = run_text(["lscpu"])
    for line in output.splitlines():
        if line.lower().startswith("model name:"):
            return line.split(":", 1)[1].strip()
    return "CPU_NO_IDENTIFICADA"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--petsc-dir", type=Path, required=True)
    parser.add_argument("--build-dir", type=Path, required=True)
    parser.add_argument("--petsc-arch", default=os.environ.get("PETSC_ARCH", ""))
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--raw-output", type=Path, required=True)
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()

    command = list(args.command)
    if command and command[0] == "--":
        command = command[1:]
    if not command:
        raise SystemExit("Debe indicar argumentos del baseline después de --.")
    if not args.petsc_arch:
        raise SystemExit("Falta --petsc-arch o la variable PETSC_ARCH.")

    executable = args.executable.resolve()
    repo = args.repo.resolve()
    petsc_dir = args.petsc_dir.resolve()
    build_dir = args.build_dir.resolve()

    if git_value(repo, "status", "--porcelain"):
        raise SystemExit("El repositorio del experimento debe estar limpio para capturar el manifest.")
    if git_value(petsc_dir, "status", "--porcelain"):
        raise SystemExit("El árbol fuente PETSc debe estar limpio para capturar el manifest.")

    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "1"
    env["OPENBLAS_NUM_THREADS"] = "1"
    env["MKL_NUM_THREADS"] = "1"
    env["NUMA_COMMIT"] = git_value(repo, "rev-parse", "HEAD")

    full_command = [str(executable), *command]
    result = subprocess.run(
        full_command,
        cwd=repo,
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    args.raw_output.parent.mkdir(parents=True, exist_ok=True)
    args.raw_output.write_text(
        result.stdout + ("\n### STDERR\n" + result.stderr if result.stderr else ""),
        encoding="utf-8",
    )
    if result.returncode != 0:
        print(result.stdout, end="")
        print(result.stderr, end="", file=sys.stderr)
        return result.returncode

    metadata = parse_metadata(result.stdout)
    compiler, compile_command = compile_entry(
        build_dir, "/apps/petsc_baseline.cpp"
    )

    petsc_commit = git_value(petsc_dir, "rev-parse", "HEAD")
    petsc_tag = git_value(petsc_dir, "describe", "--tags", "--always")
    petsc_library = (petsc_dir / args.petsc_arch / "lib" / "libpetsc.so").resolve()
    if not petsc_library.exists():
        raise SystemExit(f"No existe la biblioteca PETSc esperada: {petsc_library}")
    openblas_library = (petsc_dir / args.petsc_arch / "lib" / "libopenblas.so.0").resolve()
    if not openblas_library.exists():
        raise SystemExit(f"No existe OpenBLAS esperado: {openblas_library}")
    petsc_variables = petsc_dir / args.petsc_arch / "lib" / "petsc" / "conf" / "petscvariables"
    configure_options = ""
    for line in petsc_variables.read_text(encoding="utf-8").splitlines():
        if line.startswith("CONFIGURE_OPTIONS ="):
            configure_options = line.split("=", 1)[1].strip()
            break
    if not configure_options:
        raise SystemExit("No se encontraron CONFIGURE_OPTIONS de PETSc.")
    options = [
        f"-snes_type {metadata['snes_type']}",
        f"-snes_linesearch_type {metadata['linesearch_type']}",
        f"-snes_rtol {metadata['snes_rtol']}",
        f"-snes_max_it {metadata['snes_max_it']}",
        f"-ksp_type {metadata['ksp_type']}",
        f"-ksp_rtol {metadata['ksp_rtol']}",
        f"-ksp_max_it {metadata['ksp_max_it']}",
        f"-ksp_gmres_restart {metadata['gmres_restart']}",
        f"-pc_type {metadata['pc_type']}",
    ]

    manifest = {
        "schema_version": 1,
        "solver": "PETSc",
        "solver_version": metadata["petsc_version"],
        "nonlinear_method": "SNESNEWTONLS",
        "linear_method": metadata["ksp_type"].upper(),
        "preconditioner": (
            "BratuJacobi via PCSHELL"
            if metadata["pc_type"] == "shell"
            else "none"
        ),
        "matrix_mode": metadata["matrix_mode"],
        "options": options,
        "repository_commit": env["NUMA_COMMIT"],
        "executable_sha256": sha256(executable),
        "build": {
            "compiler": compiler,
            "flags": compile_command,
        },
        "hardware": {
            "cpu": cpu_model(),
            "threads": 1,
        },
        "command": shlex.join(full_command),
        "problems": ["B1"],
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "petsc_source_commit": petsc_commit,
        "petsc_source_tag": petsc_tag,
        "petsc_arch": args.petsc_arch,
        "petsc_configure_options": configure_options,
        "petsc_library_sha256": sha256(petsc_library),
        "openblas_library_sha256": sha256(openblas_library),
        "thread_environment": {
            "OMP_NUM_THREADS": "1",
            "OPENBLAS_NUM_THREADS": "1",
            "MKL_NUM_THREADS": "1",
        },
        "raw_output_sha256": sha256(args.raw_output),
    }

    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    args.manifest.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print("PETSC_BASELINE_CAPTURE=PASS")
    print(f"manifest={args.manifest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
