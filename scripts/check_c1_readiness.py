#!/usr/bin/env python3
"""Verifica el freeze completo antes de permitir C1 confirmatoria."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from pathlib import Path

REQUIRED_ARTIFACTS = (
    "preregistration",
    "instances",
    "starts",
    "primitive_costs",
    "internal_tuning",
    "petsc_freeze",
    "tf2018",
    "hardware",
    "statistical_protocol",
    "statistical_script",
    "campaign_manifest",
)

FROZEN_CODE_PREFIXES = (
    "apps/",
    "include/",
    "src/",
    "scripts/",
    "tests/",
    "cmake/",
    "bench/",
    "benchmarks/",
)

FROZEN_CODE_FILES = (
    "CMakeLists.txt",
)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def git(*args: str) -> str:
    result = subprocess.run(
        ["git", *args], text=True, capture_output=True, check=False
    )
    if result.returncode != 0:
        raise SystemExit(result.stderr.strip() or "Falló Git.")
    return result.stdout.strip()


def load_json(path: Path) -> dict:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise SystemExit(f"JSON inválido: {path}: {exc}") from exc


def resolve_artifact(root: Path, entry: dict, name: str) -> Path:
    if not isinstance(entry, dict):
        raise SystemExit(f"Artefacto inválido: {name}.")
    rel = entry.get("path")
    expected = entry.get("sha256")
    if not isinstance(rel, str) or not rel:
        raise SystemExit(f"Falta path de {name}.")
    if not isinstance(expected, str) or len(expected) != 64:
        raise SystemExit(f"SHA-256 inválido de {name}.")
    path = (root / rel).resolve()
    try:
        path.relative_to(root.resolve())
    except ValueError as exc:
        raise SystemExit(f"Ruta fuera del repositorio: {name}.") from exc
    if not path.is_file():
        raise SystemExit(f"No existe artefacto requerido: {name}: {path}.")
    actual = sha256(path)
    if actual != expected:
        raise SystemExit(
            f"Hash cambiado para {name}: esperado={expected} actual={actual}"
        )
    return path


def check_json_states(
    artifacts: dict[str, Path], source_commit: str
) -> None:
    prereg = load_json(artifacts["preregistration"])
    if prereg.get("status") != "FROZEN_BEFORE_C1":
        raise SystemExit("Preregistro no está FROZEN_BEFORE_C1.")
    if prereg.get("source_commit") != source_commit:
        raise SystemExit("Preregistro usa otro source_commit.")

    petsc = load_json(artifacts["petsc_freeze"])
    if petsc.get("status") != "FROZEN_CALIBRATION_ONLY":
        raise SystemExit("PETSc no está FROZEN_CALIBRATION_ONLY.")
    if petsc.get("source_commit") != source_commit:
        raise SystemExit("Freeze PETSc usa otro source_commit.")

    hardware = load_json(artifacts["hardware"])
    for field in ("cpu", "compiler", "threads"):
        if field not in hardware:
            raise SystemExit(f"Hardware manifest sin {field}.")

    campaign = load_json(artifacts["campaign_manifest"])
    if campaign.get("campaign") != "C1":
        raise SystemExit("Campaign manifest no corresponde a C1.")
    if campaign.get("source_commit") != source_commit:
        raise SystemExit("Campaign manifest usa otro source_commit.")


def git_is_ancestor(ancestor: str, descendant: str) -> bool:
    result = subprocess.run(
        ["git", "merge-base", "--is-ancestor", ancestor, descendant],
        text=True,
        capture_output=True,
        check=False,
    )
    if result.returncode == 0:
        return True
    if result.returncode == 1:
        return False
    raise SystemExit(result.stderr.strip() or "Falló git merge-base.")


def check_frozen_code_unchanged(source_commit: str, head: str) -> None:
    changed = [
        path
        for path in git(
            "diff", "--name-only", f"{source_commit}..{head}"
        ).splitlines()
        if path
    ]
    protected = [
        path
        for path in changed
        if path in FROZEN_CODE_FILES
        or path.startswith(FROZEN_CODE_PREFIXES)
    ]
    if protected:
        raise SystemExit(
            "Cambió código ejecutable después de source_commit: "
            + ", ".join(protected[:10])
        )


def check_git(source_commit: str, tag: str) -> None:
    head = git("rev-parse", "HEAD")
    if not git_is_ancestor(source_commit, head):
        raise SystemExit(
            f"source_commit={source_commit} no es ancestro de HEAD={head}."
        )

    status = git("status", "--porcelain")
    if status:
        raise SystemExit("El árbol Git no está limpio.")

    check_frozen_code_unchanged(source_commit, head)

    tag_commit = git("rev-list", "-n", "1", tag)
    if tag_commit != head:
        raise SystemExit(
            "El tag de preregistro no apunta al commit preregistrado actual."
        )


def check_confirmation_not_started(root: Path) -> None:
    confirmation = root / "experiments/c1/confirmation"
    if confirmation.exists():
        csvs = sorted(confirmation.rglob("*.csv"))
        if csvs:
            raise SystemExit(
                "Ya existen CSV confirmatorios antes del readiness: "
                + ", ".join(str(p.relative_to(root)) for p in csvs[:5])
            )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("freeze_manifest")
    parser.add_argument(
        "--skip-git",
        action="store_true",
        help="Sólo para smoke del checker, nunca concede readiness real.",
    )
    args = parser.parse_args()

    manifest_path = Path(args.freeze_manifest).resolve()
    root = Path.cwd().resolve()
    manifest = load_json(manifest_path)

    if manifest.get("schema_version") != 1:
        raise SystemExit("schema_version de readiness debe ser 1.")
    if manifest.get("campaign") != "C1":
        raise SystemExit("El manifest no corresponde a C1.")
    if manifest.get("status") != "FROZEN_BEFORE_C1":
        raise SystemExit("El manifest maestro no está FROZEN_BEFORE_C1.")
    if manifest.get("author_ratified") is not True:
        raise SystemExit("Falta ratificación explícita del autor.")
    if manifest.get("confirmation_started") is not False:
        raise SystemExit("C1 ya figura como iniciada.")

    source_commit = manifest.get("source_commit")
    tag = manifest.get("preregistration_tag")
    if not isinstance(source_commit, str) or len(source_commit) < 7:
        raise SystemExit("source_commit inválido.")
    if not isinstance(tag, str) or not tag:
        raise SystemExit("preregistration_tag inválido.")

    entries = manifest.get("artifacts")
    if not isinstance(entries, dict):
        raise SystemExit("Falta mapa artifacts.")

    missing = [name for name in REQUIRED_ARTIFACTS if name not in entries]
    if missing:
        raise SystemExit("Faltan artefactos: " + ",".join(missing))

    artifacts = {
        name: resolve_artifact(root, entries[name], name)
        for name in REQUIRED_ARTIFACTS
    }
    check_json_states(artifacts, source_commit)
    check_confirmation_not_started(root)

    if args.skip_git:
        print("C1_READINESS_SMOKE=PASS git=SKIPPED")
        return

    check_git(source_commit, tag)
    print("C1_READINESS=PASS")
    print("READY_FOR_C1")


if __name__ == "__main__":
    main()
