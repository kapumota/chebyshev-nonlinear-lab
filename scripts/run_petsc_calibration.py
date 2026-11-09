#!/usr/bin/env python3
"""Ejecuta el plan congelado de calibración PETSc para C1.

El plan define instancias y candidatos antes de observar resultados. Este runner
fuerza phase=calibration, un commit limpio, un hilo y el núcleo PETSc fijado.
No calcula W_eq ni selecciona configuraciones.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys
import tempfile


FAMILIES = ("B1", "B5", "B3a", "Bext")
FORBIDDEN_INSTANCE_OPTIONS = {
    "--phase", "--out", "--precond", "--rtol", "--max-nl", "--restart",
    "--maxlin", "--warmup", "--reps", "--commit",
    "-snes_type", "-snes_linesearch_type", "-ksp_type", "-ksp_rtol",
}
EXTRA_COLUMNS = (
    "petsc_candidate_id",
    "petsc_ksp_rtol",
    "petsc_restart",
    "calibration_instance_id",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def git_output(*args: str) -> str:
    result = subprocess.run(
        ["git", *args], check=True, text=True, capture_output=True
    )
    return result.stdout.strip()


def require_clean_commit(expected: str) -> None:
    head = git_output("rev-parse", "HEAD")
    if head != expected:
        raise SystemExit(
            f"El plan exige commit {expected}, pero HEAD es {head}."
        )

    # El plan y las salidas pueden ser artefactos aún no rastreados. Lo que se
    # prohíbe es ejecutar calibración con cambios de código ya rastreados que no
    # pertenezcan al source_commit declarado.
    worktree = subprocess.run(["git", "diff", "--quiet"], check=False)
    index = subprocess.run(["git", "diff", "--cached", "--quiet"], check=False)
    if worktree.returncode != 0 or index.returncode != 0:
        raise SystemExit(
            "La calibración PETSc exige código sin cambios rastreados respecto de HEAD."
        )


def load_plan(path: Path) -> dict:
    data = json.loads(path.read_text(encoding="utf-8"))
    if data.get("schema_version") != 1:
        raise SystemExit("schema_version de plan PETSc debe ser 1.")
    if data.get("campaign") != "C1" or data.get("phase") != "calibration":
        raise SystemExit("El plan debe declarar campaign=C1 y phase=calibration.")

    commit = str(data.get("source_commit", ""))
    if len(commit) != 40 or any(c not in "0123456789abcdef" for c in commit):
        raise SystemExit("source_commit debe ser un SHA-1 Git completo de 40 hexadecimales.")

    count = data.get("expected_candidate_count")
    if not isinstance(count, int) or count < 1:
        raise SystemExit("expected_candidate_count debe ser un entero positivo.")

    fixed = data.get("fixed")
    if not isinstance(fixed, dict):
        raise SystemExit("Falta el objeto fixed en el plan.")
    if fixed.get("snes_type") != "newtonls":
        raise SystemExit("C1 PETSc exige snes_type=newtonls.")
    if fixed.get("linesearch_type") != "bt":
        raise SystemExit("C1 PETSc exige linesearch_type=bt.")
    if fixed.get("ksp_type") != "gmres":
        raise SystemExit("C1 PETSc exige ksp_type=gmres.")
    if fixed.get("matrix_mode") != "MatShell analytic Jv":
        raise SystemExit("C1 PETSc exige MatShell analytic Jv.")
    if int(fixed.get("threads", 0)) != 1:
        raise SystemExit("La calibración PETSc de C1 exige threads=1.")
    for name in ("rtol", "max_nl", "maxlin", "warmup", "reps"):
        if name not in fixed:
            raise SystemExit(f"Falta fixed.{name} en el plan.")
    if float(fixed["rtol"]) <= 0.0:
        raise SystemExit("fixed.rtol debe ser positivo.")
    if int(fixed["max_nl"]) < 1 or int(fixed["maxlin"]) < 1:
        raise SystemExit("fixed.max_nl y fixed.maxlin deben ser positivos.")
    if int(fixed["warmup"]) < 0 or int(fixed["reps"]) < 1:
        raise SystemExit("fixed.warmup >= 0 y fixed.reps >= 1.")

    families = data.get("families")
    if not isinstance(families, dict) or tuple(sorted(families)) != tuple(sorted(FAMILIES)):
        raise SystemExit("El plan debe contener exactamente B1, B5, B3a y Bext.")

    for family in FAMILIES:
        block = families[family]
        instances = block.get("instances")
        candidates = block.get("candidates")
        if not isinstance(instances, list) or not instances:
            raise SystemExit(f"{family} no tiene instancias de calibración.")
        if not isinstance(candidates, list) or len(candidates) != count:
            raise SystemExit(
                f"{family} debe tener exactamente {count} candidatos PETSc."
            )

        instance_ids = set()
        for instance in instances:
            instance_id = str(instance.get("id", ""))
            args = instance.get("args")
            if not instance_id or instance_id in instance_ids:
                raise SystemExit(f"ID de instancia vacío o duplicado en {family}.")
            instance_ids.add(instance_id)
            if not isinstance(args, list) or not args or not all(isinstance(x, str) for x in args):
                raise SystemExit(f"args inválidos para {family}/{instance_id}.")
            forbidden = FORBIDDEN_INSTANCE_OPTIONS & set(args)
            if forbidden:
                raise SystemExit(
                    f"{family}/{instance_id} intenta fijar opciones reservadas: "
                    + ",".join(sorted(forbidden))
                )
            for required in ("--problem", "--family", "--n", "--param", "--start"):
                if required not in args:
                    raise SystemExit(f"{family}/{instance_id} carece de {required}.")
            family_index = args.index("--family") + 1
            if family_index >= len(args) or args[family_index] != family:
                raise SystemExit(f"{family}/{instance_id} declara otra familia en args.")

        candidate_ids = set()
        for candidate in candidates:
            candidate_id = str(candidate.get("id", ""))
            if not candidate_id or candidate_id in candidate_ids:
                raise SystemExit(f"ID de candidato vacío o duplicado en {family}.")
            candidate_ids.add(candidate_id)
            precond = str(candidate.get("precond", ""))
            if precond not in {"none", "jacobi", "mg"}:
                raise SystemExit(f"precond inválido en {family}/{candidate_id}.")
            ksp_rtol = float(candidate.get("ksp_rtol", 0.0))
            restart = int(candidate.get("restart", 0))
            if not (0.0 < ksp_rtol < 1.0) or restart < 1:
                raise SystemExit(f"ksp_rtol/restart inválidos en {family}/{candidate_id}.")

    return data


def command_for(
    executable: Path,
    fixed: dict,
    instance: dict,
    candidate: dict,
    output: Path,
) -> list[str]:
    return [
        str(executable),
        "--phase", "calibration",
        *instance["args"],
        "--precond", str(candidate["precond"]),
        "--rtol", str(fixed["rtol"]),
        "--max-nl", str(fixed["max_nl"]),
        "--restart", str(candidate["restart"]),
        "--maxlin", str(fixed["maxlin"]),
        "--warmup", str(fixed["warmup"]),
        "--reps", str(fixed["reps"]),
        "--out", str(output),
        "-snes_type", "newtonls",
        "-snes_linesearch_type", "bt",
        "-ksp_type", "gmres",
        "-ksp_rtol", str(candidate["ksp_rtol"]),
    ]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("plan_json")
    parser.add_argument("output_csv")
    parser.add_argument("run_manifest_json")
    args = parser.parse_args()

    plan_path = Path(args.plan_json).resolve()
    output_path = Path(args.output_csv).resolve()
    manifest_path = Path(args.run_manifest_json).resolve()
    plan = load_plan(plan_path)
    require_clean_commit(plan["source_commit"])

    executable = Path(plan.get("executable", "build-petsc/petsc_baseline")).resolve()
    if not executable.is_file() or not os.access(executable, os.X_OK):
        raise SystemExit(f"No existe ejecutable PETSc usable: {executable}")

    fixed = plan["fixed"]
    env = os.environ.copy()
    env["NUMA_COMMIT"] = plan["source_commit"]
    env["OMP_NUM_THREADS"] = "1"
    env["OPENBLAS_NUM_THREADS"] = "1"
    env["MKL_NUM_THREADS"] = "1"

    all_rows: list[dict[str, str]] = []
    fieldnames: list[str] | None = None
    commands: list[dict[str, str]] = []

    with tempfile.TemporaryDirectory(prefix="c1-petsc-cal-") as tmp:
        tmpdir = Path(tmp)
        for family in FAMILIES:
            block = plan["families"][family]
            for candidate in block["candidates"]:
                for instance in block["instances"]:
                    raw = tmpdir / f"{family}-{candidate['id']}-{instance['id']}.csv"
                    command = command_for(executable, fixed, instance, candidate, raw)
                    result = subprocess.run(
                        command, env=env, text=True, capture_output=True
                    )
                    if result.returncode != 0:
                        print(result.stdout, end="")
                        print(result.stderr, end="", file=sys.stderr)
                        raise SystemExit(
                            f"Falló calibración {family}/{candidate['id']}/{instance['id']}."
                        )
                    if not raw.is_file():
                        raise SystemExit(f"PETSc no produjo {raw}.")

                    with raw.open(newline="", encoding="utf-8") as handle:
                        reader = csv.DictReader(handle)
                        rows = list(reader)
                        current = list(reader.fieldnames or [])
                    if not rows:
                        raise SystemExit(
                            f"CSV vacío en {family}/{candidate['id']}/{instance['id']}."
                        )
                    if fieldnames is None:
                        fieldnames = current
                    elif current != fieldnames:
                        raise SystemExit("Los CSV PETSc no comparten el mismo esquema.")

                    for row in rows:
                        if row.get("phase") != "calibration":
                            raise SystemExit("El runner PETSc produjo una fila no calibration.")
                        if row.get("method") != "petsc" or row.get("family") != family:
                            raise SystemExit("Método o familia inesperados en salida PETSc.")
                        if row.get("commit") != plan["source_commit"]:
                            raise SystemExit("La fila PETSc no conserva source_commit.")
                        row["petsc_candidate_id"] = str(candidate["id"])
                        row["petsc_ksp_rtol"] = f"{float(candidate['ksp_rtol']):.17g}"
                        row["petsc_restart"] = str(int(candidate["restart"]))
                        row["calibration_instance_id"] = str(instance["id"])
                        all_rows.append(row)

                    commands.append({
                        "family": family,
                        "candidate_id": str(candidate["id"]),
                        "instance_id": str(instance["id"]),
                        "command": shlex.join(command),
                    })

    assert fieldnames is not None
    output_path.parent.mkdir(parents=True, exist_ok=True)
    final_fields = [*fieldnames, *EXTRA_COLUMNS]
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=final_fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(all_rows)

    run_manifest = {
        "schema_version": 1,
        "artifact": "C1_PETSC_CALIBRATION_RUN",
        "campaign": "C1",
        "phase": "calibration",
        "source_commit": plan["source_commit"],
        "plan": {"path": args.plan_json, "sha256": sha256(plan_path)},
        "executable": {"path": str(executable), "sha256": sha256(executable)},
        "output_csv": {"path": args.output_csv, "sha256": sha256(output_path)},
        "rows": len(all_rows),
        "commands": commands,
    }
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(
        json.dumps(run_manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    print(f"PETSC_CALIBRATION_RUN=PASS rows={len(all_rows)}")
    print(f"raw_sha256={sha256(output_path)}")
    print(f"plan_sha256={sha256(plan_path)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
