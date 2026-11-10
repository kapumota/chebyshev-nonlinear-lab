#!/usr/bin/env python3
"""Prepara binarios auxiliares, nunca ejecuta materializadores ni solvers."""
import argparse
import subprocess
from c1_ops import ROOT, TOOLS, sources, load

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build',action='store_true')
    parser.add_argument('--dry-run',action='store_true')
    args=parser.parse_args()
    if args.build and args.dry_run: parser.error('Modos incompatibles')
    sources()
    hardware=load(ROOT/'experiments/c1/frozen/hardware_manifest.json')
    compiler=hardware['compiler_path']
    original=ROOT/'experiments/c1/starts'
    files=[original/(name+'.cpp') for name in ['materialize_bratu_pair_arclength','materialize_coarse_start','materialize_perturbation','materialize_simple_reference']]
    files.append(TOOLS/'materialize_b3a_reference.cpp')
    for source in files:
        target=TOOLS/'bin'/source.stem
        command=[compiler,'-std=c++20','-O3','-fopenmp','-DNUMA_HAS_OPENMP=1','-I'+str(ROOT/'include'),str(source),'-o',str(target)]
        print(__import__('shlex').join(command))
        if args.build:
            target.parent.mkdir(exist_ok=True)
            if target.exists(): raise SystemExit('No se sobrescribe un binario existente: '+str(target))
            subprocess.run(command,check=True,cwd=ROOT)
if __name__=='__main__': main()
