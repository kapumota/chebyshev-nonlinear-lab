#!/usr/bin/env python3
"""Orquestación C1, sin decisiones de diseño ni cambios en fuentes congeladas."""
import argparse
import collections
import csv
import fcntl
import hashlib
import itertools
import json
import os
from pathlib import Path
import platform
import re
import subprocess
import sys
import tempfile

TOOLS = Path(__file__).resolve().parent
ROOT = TOOLS.parents[2]
PLAN = TOOLS / 'c1_execution_plan.json'
FROZEN = 'experiments/c1/frozen/'


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda: f.read(1048576), b''):
            h.update(chunk)
    return h.hexdigest()


def serial(obj):
    return json.dumps(obj, sort_keys=True, separators=(',', ':'), ensure_ascii=False)


def load(path):
    return json.loads(Path(path).read_text())


def save(path, obj):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + '.tmp')
    with tmp.open('w') as f:
        json.dump(obj, f, indent=2, sort_keys=True, ensure_ascii=False)
        f.write('\n')
        f.flush()
        os.fsync(f.fileno())
    os.replace(tmp, path)


def rows(path):
    with Path(path).open(newline='') as f:
        return list(csv.DictReader(f))


def resolve(path):
    p = Path(path)
    if p.is_absolute():
        return p
    p = (ROOT / p).resolve()
    if not p.is_relative_to(ROOT):
        raise ValueError('Ruta fuera de la copia local: ' + str(p))
    return p


def sources():
    """Recorre referencias congeladas y comprueba toda dependencia disponible."""
    refs, external = {}, {}
    def walk(obj):
        if isinstance(obj, dict):
            if isinstance(obj.get('path'), str) and re.fullmatch('[a-f0-9]{64}', str(obj.get('sha256', ''))):
                p, h = obj['path'], obj['sha256']
                if p in refs and refs[p] != h:
                    raise ValueError('Hashes congelados contradictorios: ' + p)
                refs[p] = h
            for v in obj.values():
                walk(v)
        elif isinstance(obj, list):
            for v in obj:
                walk(v)
    for p in sorted((ROOT / FROZEN).glob('*.json')):
        walk(load(p))
    required = {}
    for p, h in refs.items():
        target = resolve(p)
        if p.startswith('/') or p.startswith('build-'):
            external[p] = h
            if target.exists() and digest(target) != h:
                raise ValueError('SHA-256 externo incorrecto: ' + p)
        else:
            if not target.is_file() or digest(target) != h:
                raise ValueError('Artefacto congelado ausente o distinto: ' + p)
            required[p] = h
    # Incluye los propios manifests, aunque no tengan referencia entrante.
    for p in sorted((ROOT / FROZEN).glob('*')):
        if p.is_file():
            required[str(p.relative_to(ROOT))] = digest(p)
    return required, external


def build_plan():
    refs, external = sources()
    campaign = load(ROOT / FROZEN / 'campaign_manifest.json')
    starts = load(ROOT / FROZEN / 'starts_manifest.json')
    grid = load(resolve(campaign['instances']['path']))
    petsc = load(resolve(campaign['petsc']['config']['path']))
    tf = load(resolve(campaign['tf2018']['path']))
    policy = load(resolve(campaign['method_label_policy']['path']))
    pert = rows(resolve(starts['artifacts']['perturbed_starts']['path']))
    coarse = {r['instance_id']: r for r in rows(resolve(starts['artifacts']['coarse_starts']['path']))}
    b1 = {r['label']: r for r in rows(resolve(starts['artifacts']['b1_starts']['path']))}
    instances = {}
    for r in pert:
        label, family = r['instance_id'], r['family']
        if label in instances:
            continue
        n = int(re.search(r'_n(\d+)', label)[1])
        g = grid['families'][family]
        decode = lambda s: float(s.replace('p', '.'))
        if family == 'B1':
            z = b1[label]
            par = {'param': float(z['lambda_target']), 'rho_lambda': float(z['rho_lambda']), 'lambda_fold_h': float(z['lambda_fold_h'])}
        elif family == 'B3a':
            z = coarse[label]
            par = {'param': float(z['nu']), 'dt': float(z['dt']), 'amp': float(z['amp'])}
        elif family == 'B5':
            q = decode(label.split('_q')[1]); par = {'param': 1 - 10 ** (-q), 'q': q}
        else:
            par = {'param': decode(label.split('_dt')[1]), **g['fixed']}
        instances[label] = dict(instance=label, family=family, problem=g['problem'], n=n, parameters=par)
    artifacts = {}
    def add_artifact(key, h, argv, deps=()):
        path = 'experiments/c1/confirmation_tools/runtime/starts/' + key + '.txt'
        artifacts[key] = dict(path=path, sha256=h, argv=argv, dependencies=list(deps))
        return key
    def artpath(key):
        return artifacts[key]['path']
    def reference(i):
        label, f, n, p = i['instance'], i['family'], i['n'], i['parameters']
        key = label + '_reference'
        if key in artifacts:
            return key
        h = next(r['reference_sha256'] for r in pert if r['instance_id'] == label)
        out = 'experiments/c1/confirmation_tools/runtime/starts/' + key + '.txt'
        if f == 'B1':
            z = b1[label]
            argv = ['experiments/c1/confirmation_tools/bin/materialize_bratu_pair_arclength', str(n), z['lambda_fold_h'], z['rho_lambda'], out, out + '.continuation', label]
        elif f == 'B3a':
            argv = ['experiments/c1/confirmation_tools/bin/materialize_b3a_reference', str(n), str(p['param']), str(p['dt']), str(p['amp']), out]
        else:
            argv = ['experiments/c1/confirmation_tools/bin/materialize_simple_reference', f, str(n), str(p['param']), out]
        return add_artifact(key, h, argv)
    units = []
    methods = campaign['methods']['raw_internal'] + campaign['methods']['external']
    for label, i in sorted(instances.items()):
        f, n, p = i['family'], i['n'], i['parameters']
        configurations = []
        for kind in starts['confirmation']['families'][f]['starts']:
            if kind == 'natural':
                configurations.append((kind, 'natural', 0.0, [0], []))
            elif kind == 'continuation':
                key = reference(i)
                if artifacts[key]['sha256'] != b1[label]['continuation_sha256']:
                    raise ValueError('Referencia y continuation distintas, falta adaptador: ' + label)
                configurations.append((kind, kind, 0.0, [0], [key]))
            elif kind == 'coarse':
                z = coarse[label]; key = label + '_coarse'
                out = 'experiments/c1/confirmation_tools/runtime/starts/' + key + '.txt'
                args = [f, str(n), str(p['param'])]
                if f == 'B3a': args += [str(p['dt']), str(p['amp'])]
                add_artifact(key, z['coarse_sha256'], ['experiments/c1/confirmation_tools/bin/materialize_coarse_start', *args, out, label])
                configurations.append((kind, kind, 0.0, [0], [key]))
            else:
                distance = kind.removeprefix('perturbed_')
                ref = reference(i)
                selected = sorted((r for r in pert if r['instance_id'] == label and r['distance'] == distance), key=lambda r: int(r['seed']))
                keys = []
                for r in selected:
                    key = label + '_' + distance + '_' + r['seed']
                    out = 'experiments/c1/confirmation_tools/runtime/starts/' + key + '.txt'
                    add_artifact(key, r['perturb_sha256'], ['experiments/c1/confirmation_tools/bin/materialize_perturbation', artpath(ref), r['rho_x'], r['seed'], out], [ref])
                    keys.append(key)
                configurations.append(('perturbed', distance, starts['confirmation']['perturbation']['levels'][distance], [int(r['seed']) for r in selected], keys))
        for method in methods:
            for kind, distance, rho, seeds, keys in configurations:
                batches = [(seeds, keys)] if method != 'petsc' or kind != 'perturbed' else [([s], [k]) for s, k in zip(seeds, keys)]
                for ss, kk in batches:
                    fixed = petsc['fixed']
                    u = dict(i, method=method, preconditioner=tf['linear_configuration']['family_preconditioner'][f], start_kind=kind,
                             start_id='natural' if kind == 'natural' else label + '_' + kind,
                             distance=distance, rho_x=rho, seed=ss[0] if len(ss)==1 else ss,
                             warmup=campaign['execution']['warmup'], repetitions=campaign['execution']['measured_repetitions'],
                             threads=campaign['execution']['threads'], source_commit=campaign['source_commit'], start_artifact=kk,
                             expected_sha256=[artifacts[k]['sha256'] for k in kk])
                    if method == 'petsc': u['preconditioner'] = petsc['families'][f]['precond']
                    cid = hashlib.sha256(serial(u).encode()).hexdigest()
                    u['command_id'] = cid
                    u['output_path'] = campaign['outputs']['raw_directory'] + '/' + cid + '.csv'
                    argv = ['build-petsc/petsc_baseline' if method == 'petsc' else tf['evidence']['runner']['path'], '--phase', 'confirmation', '--family', f, '--problem', i['problem'], '--n', str(n), '--param', str(p['param']), '--precond', u['preconditioner'], '--start', kind, '--start-id', u['start_id'], '--warmup', str(u['warmup']), '--reps', str(u['repetitions']), '--commit', u['source_commit'], '--out', u['output_path']]
                    if f == 'B3a': argv += ['--dt', str(p['dt']), '--amp', str(p['amp'])]
                    if kk:
                        used = reference(i) if kind == 'perturbed' and method != 'petsc' else kk[0]
                        argv += ['--start-file', artpath(used)]
                    if method == 'petsc':
                        c = petsc['families'][f]
                        argv += ['--distance', distance, '--rho-x', str(rho), '--seed', str(ss[0]), '--rtol', str(fixed['rtol']), '--max-nl', str(fixed['max_nl']), '--maxlin', str(fixed['maxlin']), '--restart', str(c['gmres_restart']), '-snes_type', fixed['snes_type'], '-snes_linesearch_type', fixed['linesearch_type'], '-ksp_type', fixed['ksp_type'], '-ksp_rtol', str(c['ksp_rtol'])]
                    else:
                        c = campaign['internal_tuning']['selected'].get(method, {})
                        lc = tf['linear_configuration']
                        argv += ['--methods', method, '--rtol', str(lc['nonlinear_rtol']), '--max-nl', str(lc['max_nl']), '--maxlin', str(lc['maxlin']), '--restart', str(c.get('restart', lc['restart']))]
                        if 'eta_init' in c: argv += ['--eta-init', str(c['eta_init'])]
                        if kind == 'perturbed': argv += ['--dist', distance, '--seeds', str(max(ss))]
                    u['argv'] = argv
                    u['expected_rows'] = len(ss) * u['repetitions']
                    units.append(u)
    plan = dict(schema_version=2, source_commit=campaign['source_commit'], source_hashes=refs, external_dependencies=external,
                hardware=load(resolve(campaign['hardware']['path'])), label_policy=policy, instances=list(instances.values()),
                starts=artifacts, units=units, canonical_directory=campaign['outputs']['canonical_directory'],
                costs=campaign['primitive_costs']['path'], execution_started=False)
    import c1_amendment as amendment
    folder=ROOT/'experiments/c1/amendment01'
    if not (folder/'instance_id_policy.json').exists(): raise ValueError('Falta amendment01, aplique primero el patch de amendment')
    if load(folder/'instance_id_policy.json')!=amendment.policy(): raise ValueError('Política del amendment distinta')
    plan['canonicalization_policy']=load(folder/'instance_id_policy.json')
    plan['costs']='experiments/c1/amendment01/primitive_costs_v2.csv'
    state=load(folder/'c1_readiness_manifest_v2.json')['status']
    plan['amendment_status']=state
    for name in ['instance_id_policy.json','supplemental_calibration_plan.json','campaign_manifest_v2.json','preregistration_v2.json','c1_readiness_manifest_v2.json','primitive_costs_v2.csv','primitive_costs_provenance.json']:
        path=folder/name
        if path.exists(): plan['source_hashes'][str(path.relative_to(ROOT))]=digest(path)
    return plan


def audit(plan, reproduce=True):
    if reproduce and plan != build_plan():
        raise ValueError('El plan difiere de su reconstrucción determinista desde freeze')
    grid=load(ROOT/'experiments/c1/C1_INSTANCE_GRID_V1.json')
    for family,g in grid['families'].items():
        actual=[i for i in plan['instances'] if i['family']==family]
        if family=='B1':
            expected=set(itertools.product(g['sizes_confirmation'],g['confirmation']))
            observed={(i['n'],i['parameters']['rho_lambda']) for i in actual}
        elif family=='B3a':
            expected=set(itertools.product(g['sizes_confirmation'],g['confirmation']['nu'],g['confirmation']['dt']))
            observed={(i['n'],i['parameters']['param'],i['parameters']['dt']) for i in actual}
            if any(i['parameters']['amp']!=g['amp'] for i in actual): raise ValueError('Amplitud distinta del grid')
        elif family=='B5':
            expected=set(itertools.product(g['sizes_confirmation'],g['confirmation_q']))
            observed={(i['n'],i['parameters']['q']) for i in actual}
        else:
            expected=set(itertools.product(g['sizes_confirmation'],g['confirmation_dt']))
            observed={(i['n'],i['parameters']['param']) for i in actual}
        if expected!=observed: raise ValueError('Instancias distintas del grid: '+family)
    counts = collections.Counter(i['family'] for i in plan['instances'])
    if counts != dict(B1=24, B3a=27, B5=15, Bext=9): raise ValueError('Conteos de instancias incorrectos')
    if len({u['command_id'] for u in plan['units']}) != len(plan['units']): raise ValueError('command_id duplicado')
    logical = collections.Counter()
    for u in plan['units']:
        if (u['warmup'], u['repetitions'], u['threads']) != (1,7,1): raise ValueError('Parámetros de ejecución incorrectos')
        seeds = u['seed'] if isinstance(u['seed'], list) else [u['seed']]
        if u['start_kind'] == 'perturbed':
            if u['rho_x'] != dict(near=.05, mid=.20, far=.50)[u['distance']]: raise ValueError('rho_x incorrecto')
            if any(s not in range(1,13) for s in seeds): raise ValueError('Semilla incorrecta')
        for s in seeds:
            key = (u['instance'],u['method'],u['start_kind'],u['distance'],s)
            logical[key] += 1
    if any(v != 1 for v in logical.values()): raise ValueError('Combinación científica duplicada')
    methods = {k[1] for k in logical}
    if methods != {'nk','cheb','rec','tf2018','petsc'}: raise ValueError('Métodos incorrectos')
    for method in methods:
        c = collections.Counter(k[2] for k in logical if k[1] == method)
        if c != dict(natural=75, continuation=24, coarse=51, perturbed=2700): raise ValueError('Cobertura incorrecta: '+str(c))
    if len(logical) != 14250 or sum(u['expected_rows'] for u in plan['units']) != 99750: raise ValueError('Filas incorrectas')
    if plan['label_policy']['mapping'] != dict(nk='nk',cheb='cheb',rec='rec',tf2018='tf_cheb',petsc='petsc'): raise ValueError('Política canónica distinta')
    return dict(instances=dict(counts), commands=len(plan['units']), logical_combinations=len(logical), measured_rows=99750,
                warmup_solves=14250, total_solver_calls=114000)


def check_hashes(plan, external=False):
    for p,h in plan['source_hashes'].items():
        if not resolve(p).is_file() or digest(resolve(p)) != h: raise ValueError('Freeze distinto: '+p)
    if external:
        for p,h in plan['external_dependencies'].items():
            if not resolve(p).is_file() or digest(resolve(p)) != h: raise ValueError('Dependencia externa ausente o distinta: '+p)


def materialize(plan, key, verified, active=None):
    if key in verified: return
    a = plan['starts'][key]; path=resolve(a['path'])
    if not path.exists():
        for dep in a['dependencies']: materialize(plan, dep, verified)
        path.parent.mkdir(parents=True,exist_ok=True)
        argv=a['argv']; binary=resolve(argv[0])
        if not binary.is_file() or not os.access(binary,os.X_OK): raise ValueError('Materializador ausente: '+str(binary))
        log=path.with_suffix('.materialization.log')
        with log.open('ab') as out:
            result=subprocess.run([str(binary),*argv[1:]],cwd=ROOT,stdout=out,stderr=subprocess.STDOUT)
        if result.returncode: raise ValueError('Falló materializador: '+key)
    if digest(path)!=a['sha256']: raise ValueError('SHA-256 de start distinto: '+key)
    verified.add(key)


def release_perturbation(plan, key):
    a = plan['starts'][key]
    if a['dependencies']:
        path=resolve(a['path'])
        if path.exists() and digest(path)==a['sha256']: path.unlink()


def downstream_reasons(plan):
    import c1_amendment as amendment
    try:
        amendment.verify_amendment(plan)
        return []
    except Exception as exc:
        return ['Amendment pendiente o inválido: '+str(exc)]


def verify_starts(plan):
    verified=set()
    for key in plan['starts']:
        materialize(plan,key,verified)
        release_perturbation(plan,key)
    return len(verified)


def environment_reasons(plan):
    reasons=[]; h=plan['hardware']
    for k,v in {**h['thread_environment'],'NUMA_COMMIT':plan['source_commit']}.items():
        if os.environ.get(k)!=v: reasons.append('Variable '+k+' debe ser '+v)
    cpu=next((x.split(':',1)[1].strip() for x in Path('/proc/cpuinfo').read_text().splitlines() if x.startswith('model name')), '')
    if cpu!=h['cpu']: reasons.append('CPU distinta: '+cpu)
    if platform.machine()!=h['architecture']: reasons.append('Arquitectura distinta')
    if platform.release()!=h['kernel']: reasons.append('Kernel distinto')
    release=Path('/etc/os-release').read_text()
    if h['operating_system'] not in release: reasons.append('Sistema operativo distinto')
    mem=int(re.search(r'MemTotal:\s+(\d+)',Path('/proc/meminfo').read_text())[1])
    if mem!=h['memory_total_kib']: reasons.append('Memoria total distinta del manifest')
    compiler=Path(h['compiler_path'])
    if not compiler.is_file(): reasons.append('Compilador congelado ausente')
    else:
        version=subprocess.run([str(compiler),'--version'],capture_output=True,text=True).stdout.splitlines()[0]
        # El manifest puede registrar el alias c++, con el mismo GCC.
        if version[version.find('('):] != h['compiler'][h['compiler'].find('('):]: reasons.append('Versión del compilador distinta')
    for k,v in [('PETSC_DIR',h['petsc']['dir']),('PETSC_ARCH',h['petsc']['arch'])]:
        if os.environ.get(k)!=v: reasons.append(k+' distinto del manifest')
    version=Path(h['petsc']['dir'])/'include/petscversion.h'
    if not version.is_file(): reasons.append('PETSc congelado ausente')
    else:
        text=version.read_text(); parts=[]
        for k in ['MAJOR','MINOR','SUBMINOR']:
            match=re.search(r'#define\s+PETSC_VERSION_'+k+r'\s+(\d+)',text)
            parts.append(match[1] if match else '?')
        if '.'.join(parts)!=h['petsc']['version']: reasons.append('Versión PETSc distinta')
    for p in sorted({u['argv'][0] for u in plan['units']} | {a['argv'][0] for a in plan['starts'].values()}):
        if not resolve(p).is_file() or not os.access(resolve(p),os.X_OK): reasons.append('Ejecutable ausente: '+p)
    return reasons


def validate_effective_petsc(plan,u,path):
    if u['method']!='petsc': return
    text=path.with_suffix('.stdout').read_text()
    section=text.split('PETSC_BASELINE_BEGIN',1)[-1].split('PETSC_BASELINE_END',1)[0]
    metadata=dict(line.split('=',1) for line in section.splitlines() if '=' in line)
    argv=u['argv']
    expected={'petsc_version':plan['hardware']['petsc']['version'],
              'family':u['family'],'application_preconditioner':u['preconditioner'],
              'snes_type':argv[argv.index('-snes_type')+1],
              'linesearch_type':argv[argv.index('-snes_linesearch_type')+1],
              'ksp_type':argv[argv.index('-ksp_type')+1],
              'matrix_mode':'MatShell analytic Jv'}
    for k,v in expected.items():
        actual = metadata.get(k)
        if k == 'petsc_version':
            match = re.search(
                r'(?<!\d)(\d+)\.(\d+)\.(\d+)(?![\d.])',
                actual or ''
            )
            if match is None or match.group(0) != v:
                raise ValueError('Configuración PETSc efectiva distinta: '+k)
        elif actual != v:
            raise ValueError('Configuración PETSc efectiva distinta: '+k)
    for k,flag in [('ksp_rtol','-ksp_rtol'),('gmres_restart','--restart'),('snes_rtol','--rtol'),('snes_max_it','--max-nl'),('ksp_max_it','--maxlin')]:
        if float(metadata[k])!=float(argv[argv.index(flag)+1]): raise ValueError('Parámetro PETSc efectivo distinto: '+k)


def validate_raw(plan,u,path):
    data=rows(path); seeds=u['seed'] if isinstance(u['seed'],list) else [u['seed']]
    expected=set(itertools.product(seeds,range(1,u['repetitions']+1))); seen=set()
    for r in data:
        k=(int(r['seed']),int(r['repetition']))
        if k in seen or k not in expected: raise ValueError('Semilla/repetición duplicada o inesperada')
        seen.add(k)
        for field,value in [('method',u['method']),('family',u['family']),('problem',u['problem']),('n',u['n']),('precond',u['preconditioner']),('start_kind',u['start_kind']),('start_id',u['start_id']),('distance',u['distance']),('threads',u['threads']),('commit',u['source_commit']),('phase','confirmation')]:
            if r[field]!=str(value): raise ValueError('Metadato raw incorrecto: '+field)
        if float(r['rho_x'])!=u['rho_x']: raise ValueError('rho_x raw distinto')
        if float(r['param'])!=u['parameters']['param']: raise ValueError('Parámetro raw distinto')
        N=u['n']**3 if u['family']=='B1' else 2*u['n']**2 if u['family']=='B3a' else u['n'] if u['family']=='B5' else 2*u['n']
        if int(r['N'])!=N: raise ValueError('Dimensión raw distinta')
        for field,flag in [('rtol','--rtol'),('max_nl','--max-nl'),('restart','--restart'),('maxlin','--maxlin')]:
            if float(r[field])!=float(u['argv'][u['argv'].index(flag)+1]): raise ValueError('Configuración raw distinta: '+field)
        params=u['parameters'];fmt='.17g' if u['method']=='petsc' else '.6f'
        v=format(params['param'],fmt)
        if u['family']=='B1': iid='lambda='+v
        elif u['family']=='B5': iid='c='+v
        elif u['family']=='B3a': iid='nu='+v+';dt='+format(params['dt'],fmt)+';amp='+format(params['amp'],fmt)
        else: iid='dt='+v+';A=1;B=3;alpha=0.02'
        if r['instance_id']!=iid: raise ValueError('instance_id raw no corresponde a los parámetros del comando')
        if r['success'] not in ('0','1'): raise ValueError('success inválido')
    if seen!=expected: raise ValueError('Filas incompletas')
    return data


def preflight(plan,verify=True,initial=True):
    reasons=[]
    try: audit(plan); check_hashes(plan,external=True)
    except Exception as e: reasons.append(str(e))
    reasons+=environment_reasons(plan)
    reasons+=downstream_reasons(plan)
    try:
        __import__('c1_amendment').git_protected()
    except Exception as exc: reasons.append(str(exc))
    raw=ROOT/'experiments/c1/confirmation'
    if initial and raw.exists() and any(raw.rglob('*')): reasons.append('Ya existen archivos bajo confirmation')
    if initial and (TOOLS/'runtime/ledger.json').exists(): reasons.append('Ya existe ledger')
    if verify and not reasons:
        try: verify_starts(plan)
        except Exception as e: reasons.append(str(e))
    return reasons


def run(plan,execute=False,verify_only=False):
    audit(plan);check_hashes(plan)
    if not execute and not verify_only:
        print(json.dumps(audit(plan),indent=2)); return
    reasons=environment_reasons(plan)
    if reasons: raise ValueError('; '.join(reasons))
    check_hashes(plan,external=True)
    if verify_only:
        print('Starts verificados:',verify_starts(plan));return
    runtime=TOOLS/'runtime'; runtime.mkdir(exist_ok=True)
    with (runtime/'execution.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        ledger_path=runtime/'ledger.json'; marker=runtime/'confirmation_started.json'
        fingerprint=hashlib.sha256(serial(plan).encode()).hexdigest()
        if not marker.exists():
            reasons=preflight(plan)
            if reasons: raise ValueError('; '.join(reasons))
            with marker.open('x') as f:
                json.dump(dict(confirmation_started=True,plan_sha256=fingerprint),f); f.flush();os.fsync(f.fileno())
            save(ledger_path,dict(plan_sha256=fingerprint,units={u['command_id']:dict(state='not_started') for u in plan['units']}))
        if not ledger_path.exists(): raise ValueError('Inicio registrado sin ledger, requiere recuperación manual sin reejecución')
        ledger=load(ledger_path)
        if load(marker)['plan_sha256']!=fingerprint or ledger['plan_sha256']!=fingerprint: raise ValueError('Plan cambiado después del inicio')
        for u in plan['units']:
            cid=u['command_id']; entry=ledger['units'][cid]; path=resolve(u['output_path'])
            if entry['state'].startswith('completed_'):
                if digest(path)!=entry['raw_sha256']: raise ValueError('Raw completo alterado: '+cid)
                continue
            if entry['state'] in ('running','infrastructure_failure'):
                # Recupera sólo un CSV íntegro con prueba durable de salida normal.
                receipt=path.with_suffix('.receipt.json')
                if receipt.exists() and load(receipt).get('returncode')==0:
                    data=validate_raw(plan,u,path)
                    validate_effective_petsc(plan,u,path)
                    entry.update(state='completed_success' if all(r['success']=='1' for r in data) else 'completed_solver_failure',raw_sha256=digest(path))
                    save(ledger_path,ledger);continue
                # Una unidad incompleta se conserva como evidencia de infraestructura.
                # Se reanuda el comando entero, sin seleccionar semillas/repeticiones.
                attempt=int(entry.get('attempt',1))
                archive=runtime/'infrastructure_attempts'/cid/str(attempt)
                archive.mkdir(parents=True,exist_ok=False)
                for previous in [path,path.with_suffix('.stdout'),path.with_suffix('.stderr'),receipt]:
                    if previous.exists(): os.replace(previous,archive/previous.name)
                save(archive/'ledger_entry.json',entry)
                entry.clear();entry.update(state='not_started',attempt=attempt+1)
                save(ledger_path,ledger)
            if path.exists(): raise ValueError('Raw previo sin ledger completo: '+cid)
            try:
                for key in u['start_artifact']:
                    # Vuelve a calcular el hash antes de cada unidad, aun con caché.
                    materialize(plan,key,set())
                    if u['method']!='petsc': release_perturbation(plan,key)
                path.parent.mkdir(parents=True,exist_ok=True)
                entry.update(state='running');save(ledger_path,ledger)
                with path.with_suffix('.stdout').open('xb') as out,path.with_suffix('.stderr').open('xb') as err:
                    result=subprocess.run([str(resolve(u['argv'][0])),*u['argv'][1:]],cwd=ROOT,stdout=out,stderr=err)
                save(path.with_suffix('.receipt.json'),dict(returncode=result.returncode,command_id=cid,plan_sha256=fingerprint))
                if result.returncode: raise ValueError('Salida no cero del proceso: '+str(result.returncode))
                data=validate_raw(plan,u,path)
                validate_effective_petsc(plan,u,path)
                entry.update(state='completed_success' if all(r['success']=='1' for r in data) else 'completed_solver_failure',raw_sha256=digest(path))
                save(ledger_path,ledger)
                for key in u['start_artifact']: release_perturbation(plan,key)
            except BaseException as e:
                entry.update(state='infrastructure_failure',reason=str(e));save(ledger_path,ledger);raise


def postprocess(plan):
    audit(plan);check_hashes(plan)
    reasons=downstream_reasons(plan)
    if reasons: raise ValueError('; '.join(reasons))
    ledger=load(TOOLS/'runtime/ledger.json'); allrows=[]; canonical_instances=[]; hashes={}; fields=None
    if ledger['plan_sha256']!=hashlib.sha256(serial(plan).encode()).hexdigest(): raise ValueError('Ledger asociado a otro plan')
    for u in plan['units']:
        entry=ledger['units'][u['command_id']]
        if entry['state'] not in ('completed_success','completed_solver_failure'): raise ValueError('Unidad incompleta: '+u['command_id'])
        path=resolve(u['output_path'])
        if digest(path)!=entry['raw_sha256']: raise ValueError('Raw alterado')
        data=validate_raw(plan,u,path)
        validate_effective_petsc(plan,u,path)
        if fields is None: fields=list(data[0])
        if list(data[0])!=fields: raise ValueError('Esquemas CSV distintos')
        allrows+=data
        canonical_instances.extend([u]*len(data))
        hashes[u['output_path']]=digest(path)
        for suffix in ['.stdout','.stderr','.receipt.json']:
            companion=path.with_suffix(suffix)
            hashes[str(companion.relative_to(ROOT))]=digest(companion)
    if len(allrows)!=99750: raise ValueError('Completitud incorrecta')
    dest=resolve(plan['canonical_directory']);dest.mkdir(parents=True,exist_ok=True)
    canonical=dest/'c1_canonical.csv'
    # No sobrescribe análisis ni datasets previamente creados.
    with canonical.open('x',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
        for row,instance in zip(allrows,canonical_instances,strict=True):
            w.writerow(__import__('c1_amendment').canonicalize_row(row,instance))
    reread=rows(canonical)
    for raw,canon,instance in zip(allrows,reread,canonical_instances,strict=True):
        __import__('c1_amendment').audit_canonical_row(raw,canon,instance)
    save(dest/'raw_manifest.json',hashes)
    weq=dest/'c1_weq.csv'
    if weq.exists(): raise ValueError('Dataset W_eq ya existe')
    subprocess.run([sys.executable,str(ROOT/'scripts/compute_weq.py'),str(resolve(plan['costs'])),str(canonical),str(weq)],check=True,cwd=ROOT)
    save(dest/'canonical_audit.json',dict(rows=len(allrows),raw_canonical_equivalence=True,canonical_sha256=digest(canonical),weq_sha256=digest(weq),costs_sha256=digest(resolve(plan['costs'])),analysis_command=[sys.executable,'scripts/analyze_c1_confirmatory.py',str(weq),'experiments/c1/confirmation/analysis']))
    print('C1_CONFIRMATION_COMPLETENESS=PASS')


def main(mode=None):
    p=argparse.ArgumentParser(description=__doc__)
    if mode is None: p.add_argument('mode',choices=['generate','audit','run','preflight','postprocess'])
    p.add_argument('--plan',type=Path,default=PLAN)
    p.add_argument('--dry-run',action='store_true')
    p.add_argument('--execute',action='store_true')
    p.add_argument('--verify-starts-only',action='store_true')
    a=p.parse_args();mode=mode or a.mode
    try:
        if mode=='generate':
            plan=build_plan(); stats=audit(plan,False)
            if a.plan.exists() and load(a.plan)!=plan: raise ValueError('Plan existente distinto, no se sobrescribe')
            if not a.plan.exists(): save(a.plan,plan)
            print(json.dumps(stats,indent=2));return
        plan=load(a.plan)
        if mode=='audit': print(json.dumps(audit(plan),indent=2))
        elif mode=='preflight':
            reasons=preflight(plan)
            print('C1_EXECUTION_PREFLIGHT='+('BLOCKED' if reasons else 'PASS'))
            for reason in reasons: print(reason,file=sys.stderr)
            if reasons: sys.exit(1)
        elif mode=='run':
            if a.execute and a.dry_run: raise ValueError('--execute y --dry-run son incompatibles')
            run(plan,a.execute and not a.dry_run,a.verify_starts_only and not a.dry_run)
        elif mode=='postprocess':
            if a.dry_run: print('Preparado, no se leen ni generan resultados científicos')
            else: postprocess(plan)
    except Exception as e:
        if mode=='preflight': print('C1_EXECUTION_PREFLIGHT=BLOCKED')
        print(str(e),file=sys.stderr);sys.exit(1)

if __name__=='__main__': main()
