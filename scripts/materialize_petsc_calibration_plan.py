#!/usr/bin/env python3
"""Materializa el plan PETSc de calibration desde la malla C1 predefinida."""

from __future__ import annotations
import argparse, itertools, json, math, subprocess
from pathlib import Path

FAMILIES=("B1","B5","B3a","Bext")
PRECOND={"B1":"jacobi","B5":"none","B3a":"mg","Bext":"jacobi"}

def load_json(path):
    try: return json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError,json.JSONDecodeError) as exc: raise SystemExit(f"JSON inválido: {path}: {exc}")

def sha1(value):
    value=str(value).strip().lower()
    if len(value)!=40 or any(c not in "0123456789abcdef" for c in value):
        raise SystemExit("source_commit debe ser un SHA-1 completo de 40 hexadecimales.")
    return value

def head():
    return sha1(subprocess.run(["git","rev-parse","HEAD"],text=True,capture_output=True,check=True).stdout)

def fmt(x): return f"{float(x):.17g}"

def validate_grid(g):
    if g.get("schema_version")!=1 or g.get("artifact")!="C1_INSTANCE_GRID_V1":
        raise SystemExit("Malla C1 inválida.")
    if g.get("status")!="PREDEFINED_BEFORE_CALIBRATION" or g.get("campaign")!="C1":
        raise SystemExit("La malla no está predefinida para C1.")
    if set(g.get("families",{}))!=set(FAMILIES): raise SystemExit("Familias incompletas.")
    b1=g["families"]["B1"]
    if set(map(float,b1["calibration"])) & set(map(float,b1["confirmation"])): raise SystemExit("Solapamiento B1.")
    b5=g["families"]["B5"]
    if set(map(float,b5["calibration_q"])) & set(map(float,b5["confirmation_q"])): raise SystemExit("Solapamiento B5.")
    b3=g["families"]["B3a"]
    ca=set(itertools.product(map(float,b3["calibration"]["nu"]),map(float,b3["calibration"]["dt"])))
    co=set(itertools.product(map(float,b3["confirmation"]["nu"]),map(float,b3["confirmation"]["dt"])))
    if ca & co: raise SystemExit("Solapamiento B3a.")
    if (0.05,0.1) in ca or (0.05,0.1) in co: raise SystemExit("Smoke B3a reutilizado.")
    bx=g["families"]["Bext"]
    ca=set(map(float,bx["calibration_dt"])); co=set(map(float,bx["confirmation_dt"]))
    if ca & co: raise SystemExit("Solapamiento Bext.")
    if 0.01 in ca or 0.01 in co: raise SystemExit("Smoke Bext reutilizado.")
    t=g.get("petsc_tuning",{})
    if t.get("start_policy")!="natural_only" or len(t.get("candidates",[]))!=4:
        raise SystemExit("Slice PETSc inválida.")
    if t.get("fixed_preconditioner")!=PRECOND: raise SystemExit("Precondicionadores PETSc cambiados.")

def folds(path,sizes):
    d=load_json(path)
    if d.get("schema_version")!=1 or d.get("artifact")!="C1_BRATU_FOLDS":
        raise SystemExit("Artefacto de folds inválido.")
    if d.get("independent_of_compared_solver") is not True:
        raise SystemExit("Los folds deben ser independientes del solver comparado.")
    raw=d.get("folds",{}); out={}
    for n in sizes:
        e=raw.get(str(n)); v=e.get("lambda_fold_h") if isinstance(e,dict) else e
        try: v=float(v)
        except (TypeError,ValueError): raise SystemExit(f"Falta lambda_fold_h para n={n}.")
        if not math.isfinite(v) or v<=0: raise SystemExit(f"lambda_fold_h inválido para n={n}.")
        out[n]=v
    return out

def candidates(g,fam):
    p=PRECOND[fam]
    return [{"id":x["id"],"precond":p,"ksp_rtol":float(x["ksp_rtol"]),"restart":int(x["restart"])}
            for x in g["petsc_tuning"]["candidates"]]

def instances(g,fold):
    f=g["families"]; out={k:[] for k in FAMILIES}
    for n in f["B1"]["sizes_calibration"]:
        for rho in f["B1"]["calibration"]:
            lam=float(rho)*fold[int(n)]
            out["B1"].append({"id":f"B1-n{n}-rho{fmt(rho)}","args":
                ["--problem","bratu","--family","B1","--n",str(n),"--param",fmt(lam),"--start","natural"]})
    for n in f["B5"]["sizes_calibration"]:
        for q in f["B5"]["calibration_q"]:
            c=1.0-10.0**(-float(q))
            out["B5"].append({"id":f"B5-n{n}-q{fmt(q)}","args":
                ["--problem","chandra","--family","B5","--n",str(n),"--param",fmt(c),"--start","natural"]})
    b3=f["B3a"]; amp=float(b3["amp"])
    for n in b3["sizes_calibration"]:
        for nu,dt in itertools.product(b3["calibration"]["nu"],b3["calibration"]["dt"]):
            out["B3a"].append({"id":f"B3a-n{n}-nu{fmt(nu)}-dt{fmt(dt)}","args":
                ["--problem","burgers","--family","B3a","--n",str(n),"--param",fmt(nu),
                 "--dt",fmt(dt),"--amp",fmt(amp),"--start","natural"]})
    bx=f["Bext"]
    for n in bx["sizes_calibration"]:
        for dt in bx["calibration_dt"]:
            out["Bext"].append({"id":f"Bext-n{n}-dt{fmt(dt)}","args":
                ["--problem","brusselator","--family","Bext","--n",str(n),"--param",fmt(dt),"--start","natural"]})
    return out

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--grid",default="experiments/c1/C1_INSTANCE_GRID_V1.json")
    ap.add_argument("--bratu-folds",required=True)
    ap.add_argument("--output",required=True)
    ap.add_argument("--source-commit")
    a=ap.parse_args()
    g=load_json(a.grid); validate_grid(g)
    commit=sha1(a.source_commit) if a.source_commit else head()
    sizes=[int(x) for x in g["families"]["B1"]["sizes_calibration"]]
    ins=instances(g,folds(a.bratu_folds,sizes))
    plan={"schema_version":1,"campaign":"C1","phase":"calibration","source_commit":commit,
          "executable":"build-petsc/petsc_baseline","expected_candidate_count":4,
          "fixed":g["petsc_tuning"]["fixed"],
          "families":{fam:{"instances":ins[fam],"candidates":candidates(g,fam)} for fam in FAMILIES}}
    p=Path(a.output); p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(plan,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print("C1_PETSC_PLAN_MATERIALIZED=PASS "+" ".join(f"{f}={len(ins[f])}" for f in FAMILIES))

if __name__=="__main__": main()
