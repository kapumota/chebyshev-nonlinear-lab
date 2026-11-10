#!/usr/bin/env python3
"""Amendment C1: suplemento n=511 y congelación separada, sin tuning."""
import argparse
import collections
import copy
import csv
import fcntl
import hashlib
import itertools
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import c1_ops as ops

ROOT=ops.ROOT
TOOLS=ops.TOOLS
AMEND=ROOT/'experiments/c1/amendment01'
SUPPLEMENT=AMEND/'supplemental_calibration_plan.json'
COST_FIELDS=('problem','family','n','N','precond','threads')
PROTECTED=('apps','include','src','scripts','tests','cmake','bench','benchmarks','CMakeLists.txt')
SOURCE='70e565b726595697ccc09c7a1d58b3a0e2dc8eff'


def ref(path):
    p=Path(path)
    return {'path':str(p.relative_to(ROOT)), 'sha256':ops.digest(p)}


def assert_not_started():
    for path in [TOOLS/'runtime/confirmation_started.json',TOOLS/'runtime/ledger.json']:
        if path.exists(): raise ValueError('Ya existe registro de confirmation: '+str(path))
    directory=ROOT/'experiments/c1/confirmation'
    if directory.exists() and any(p.is_file() for p in directory.rglob('*')):
        raise ValueError('Ya existen archivos confirmatorios')
    for name in ['preregistration.json','campaign_manifest.json','c1_readiness_manifest.json']:
        if ops.load(ROOT/ops.FROZEN/name).get('confirmation_started') is not False:
            raise ValueError('El freeze original indica confirmation iniciada')


def git_protected():
    def git(*args):
        p=subprocess.run(['git','-C',str(ROOT),*args],capture_output=True,text=True)
        if p.returncode: raise ValueError('Verificación Git local: '+p.stderr.strip())
        return p.stdout.strip()
    git('cat-file','-e',SOURCE+'^{commit}')
    git('merge-base','--is-ancestor',SOURCE,'HEAD')
    head_changed=git('diff','--name-only',SOURCE,'HEAD','--',*PROTECTED)
    if head_changed: raise ValueError('HEAD contiene cambios en código protegido: '+head_changed)
    changed=git('diff','--name-only',SOURCE,'--',*PROTECTED)
    untracked=git('ls-files','--others','--exclude-standard','--',*PROTECTED)
    if changed or untracked:
        raise ValueError('Código protegido distinto de source_commit: '+changed+' '+untracked)


def supplement_plan():
    ops.sources()
    budget=ops.load(ROOT/'experiments/c1/C1_TUNING_BUDGET_V1.json')
    grid=ops.load(ROOT/'experiments/c1/C1_INSTANCE_GRID_V1.json')
    original=ops.load(ROOT/'experiments/c1/internal_calibration/internal_calibration_plan.json')
    original_run=ops.load(ROOT/'experiments/c1/internal_calibration/internal_calibration_run_manifest.json')
    campaign=ops.load(ROOT/ops.FROZEN/'campaign_manifest.json')
    if original_run['plan']['sha256']!=ops.digest(ROOT/'experiments/c1/internal_calibration/internal_calibration_plan.json'):
        raise ValueError('Plan de calibration original alterado')
    raw=ops.resolve(original_run['output_csv']['path'])
    if ops.digest(raw)!=original_run['output_csv']['sha256']: raise ValueError('Raw original alterado')
    if budget['internal_methods']!=['nk','cheb','rec']: raise ValueError('Métodos de calibration inesperados')
    if [c['id'] for c in budget['internal_candidates']]!=['t01','t02','t03','t04']: raise ValueError('Candidatos inesperados')
    units=[]
    for family in ['B3a','Bext']:
        g=grid['families'][family]
        if 511 not in g['sizes_confirmation']: raise ValueError('n=511 no está preregistrado')
        parameters=([{'param':nu,'dt':dt,'amp':g['amp']} for nu,dt in itertools.product(g['calibration']['nu'],g['calibration']['dt'])]
                    if family=='B3a' else [{'param':dt,**g['fixed']} for dt in g['calibration_dt']])
        for candidate in budget['internal_candidates']:
            for pars in parameters:
                def match(inv):
                    if inv['family']!=family or inv['candidate_id']!=candidate['id']: return False
                    cmd=inv['command']
                    if float(cmd[cmd.index('--param')+1])!=pars['param']: return False
                    if family=='B3a' and float(cmd[cmd.index('--dt')+1])!=pars['dt']: return False
                    return True
                template=next((inv for inv in original['invocations'] if match(inv)),None)
                if template is None: raise ValueError('Falta plantilla congelada para '+family+' '+str(pars))
                argv=list(template['command']);argv[0]='build-petsc/run_campaign'
                argv[argv.index('--n')+1]='511'
                for flag,value in [('--eta-init',candidate['eta_init']),('--restart',candidate['restart']),('--warmup',1),('--reps',7),('--maxlin',budget['fixed_budget']['maxlin']),('--max-nl',budget['fixed_budget']['max_nl'])]:
                    if float(argv[argv.index(flag)+1])!=value: raise ValueError('Plantilla y budget no concuerdan: '+flag)
                for flag,value in [('--phase','calibration'),('--start','natural'),('--methods',','.join(budget['internal_methods'])),('--commit',SOURCE),('--precond',budget['family_preconditioner'][family])]:
                    if argv[argv.index(flag)+1]!=value: raise ValueError('Plantilla congelada incompatible: '+flag)
                u=dict(family=family,problem=g['problem'],n=511,N=522242 if family=='B3a' else 1022,
                       parameters=pars,candidate=candidate,preconditioner=budget['family_preconditioner'][family],
                       methods=budget['internal_methods'],warmup=1,repetitions=7,threads=1,source_commit=SOURCE)
                cid=hashlib.sha256(ops.serial(u).encode()).hexdigest()
                u.update(command_id=cid,output_path='experiments/c1/amendment01/supplement/raw/'+cid+'.csv',expected_rows=21)
                argv[argv.index('--out')+1]=u['output_path'];u['argv']=argv;units.append(u)
    if collections.Counter(u['family'] for u in units)!=dict(B3a=36,Bext=12): raise ValueError('Conteos suplementarios distintos')
    inputs={str(p.relative_to(ROOT)):ops.digest(p) for p in [ROOT/'experiments/c1/C1_TUNING_BUDGET_V1.json',ROOT/'experiments/c1/C1_INSTANCE_GRID_V1.json',ROOT/'experiments/c1/internal_calibration/internal_calibration_plan.json',ROOT/'experiments/c1/internal_calibration/internal_calibration_run_manifest.json',raw,ROOT/'scripts/freeze_primitive_costs.py',ROOT/ops.FROZEN/'primitive_costs.csv']}
    return dict(artifact='C1_AMENDMENT01_SUPPLEMENTAL_CALIBRATION_PLAN',source_commit=SOURCE,phase='calibration',
                confirmation_started=False,no_confirmatory_data_observed=True,units=units,expected_rows=1008,
                input_hashes=inputs,executable=original['executable'],original_raw=str(raw.relative_to(ROOT)),
                hardware=ops.load(ROOT/ops.FROZEN/'hardware_manifest.json'),environment=campaign['execution'])


def validate_plan(plan):
    if plan!=supplement_plan(): raise ValueError('Plan suplementario distinto de los artefactos congelados')
    if len(plan['units'])!=48 or sum(u['expected_rows'] for u in plan['units'])!=1008: raise ValueError('Conteos suplementarios incorrectos')


def canonical_id(instance):
    """Identidad completa y sin redondeo numérico, desde parámetros del plan."""
    params={k:float(v).hex() for k,v in sorted(instance['parameters'].items())}
    # No usa comas ni comillas, es estable y reversible a binary64.
    return '|'.join(['C1v2',instance['family'],instance['problem'],'n='+str(instance['n']),*[k+'='+v for k,v in params.items()]])


def canonicalize_row(row,instance):
    canonical=dict(row)
    if canonical['method']=='tf2018': canonical['method']='tf_cheb'
    canonical['instance_id']=canonical_id(instance)
    audit_canonical_row(row,canonical,instance)
    return canonical


def audit_canonical_row(raw,canonical,instance):
    if set(raw)!=set(canonical): raise ValueError('Columnas distintas entre raw y canonical')
    for field in raw:
        expected=canonical_id(instance) if field=='instance_id' else 'tf_cheb' if field=='method' and raw[field]=='tf2018' else raw[field]
        if canonical[field]!=expected: raise ValueError('Celda alterada sin autorización: '+field)


def policy():
    return dict(artifact='C1_AMENDMENT01_INSTANCE_ID_POLICY',schema_version=1,source_commit=SOURCE,
                confirmation_started=False,no_confirmatory_data_observed=True,reason='pre-confirmation operational defect correction',
                plan_source='experiments/c1/confirmation_tools/c1_execution_plan.json',
                algorithm='C1v2|family|problem|n=integer|key=float.hex(value), parameter keys in lexicographic order; separators are vertical bars',
                precision='exact IEEE-754 binary64 round-trip, Python float.hex, no change to solver arguments',
                raw_preserved_byte_for_byte=True,allowed_changes=['method','instance_id'],
                method_mapping={'tf2018':'tf_cheb'},all_other_cells_unchanged=True,
                grouping='family, problem, n, all scientific parameters from the deterministic frozen instance plan',
                implementation='experiments/c1/confirmation_tools/c1_amendment.py:canonical_id')


def initial_manifests():
    """Plantillas pendientes, nunca declaran un hash ni costo inexistente."""
    common=dict(amendment='01',confirmation_started=False,no_confirmatory_data_observed=True,
                reason='pre-confirmation operational defect correction',status='PENDING_SUPPLEMENTAL_CALIBRATION',
                original_preregistration=ref(ROOT/ops.FROZEN/'preregistration.json'),
                original_preregistration_commit='26a83529be36d5f702821cccd5f616ab6cff5374',
                amendment_document=ref(ROOT/'docs/research/C1_PRECONFIRMATION_AMENDMENT_01.md'),
                canonical_changes=['method','instance_id'],
                unchanged=['H1-H4','families','confirmation_grids','starts','seeds','methods','tuning','PETSc','TF2018','W_eq_formula','bootstrap','comparisons','failure_treatment','hardware','source_commit'])
    result={}
    for name in ['campaign_manifest','preregistration','c1_readiness_manifest']:
        original=ROOT/ops.FROZEN/(name+'.json');d=ops.load(original)
        d.update(copy.deepcopy(common));d['original_manifest']=ref(original)
        d['artifact']=d['artifact']+'_AMENDMENT01_V2'
        d['instance_id_policy']=ref(AMEND/'instance_id_policy.json')
        costs={'path':'experiments/c1/amendment01/primitive_costs_v2.csv','sha256':None,'state':'PENDING_LOCAL_SUPPLEMENT'}
        if name=='c1_readiness_manifest': d['artifacts']['primitive_costs']=costs
        else: d['primitive_costs']=costs
        if name=='campaign_manifest':
            d['canonical_changes']=['method','instance_id']
        if name=='preregistration':
            d['canonical_changes']=['method','instance_id']
        result[name+'_v2.json']=d
    # Las referencias cruzadas se completan al finalizar localmente el suplemento.
    return result


def prepare():
    assert_not_started();AMEND.mkdir(parents=True,exist_ok=True)
    def create(path,obj):
        if path.exists():
            if ops.load(path)!=obj: raise ValueError('No se sobrescribe artefacto preparado distinto: '+str(path))
        else: ops.save(path,obj)
    create(AMEND/'instance_id_policy.json',policy())
    create(SUPPLEMENT,supplement_plan())
    for name,data in initial_manifests().items(): create(AMEND/name,data)
    print('C1_AMENDMENT_PREPARATION=READY\nC1_CONFIRMATION=NOT_STARTED')


def validate_unit(plan,u,path):
    data=ops.rows(path)
    expected=set(itertools.product(u['methods'],range(1,8)));seen=set()
    argv=u['argv']
    def flag(k): return argv[argv.index(k)+1]
    pars=u['parameters']
    iid=('nu='+format(pars['param'],'.6f')+';dt='+format(pars['dt'],'.6f')+';amp='+format(pars['amp'],'.6f') if u['family']=='B3a' else 'dt='+format(pars['param'],'.6f')+';A=1;B=3;alpha=0.02')
    for r in data:
        k=(r['method'],int(r['repetition']))
        if k not in expected or k in seen: raise ValueError('Método/repetición ausente, inesperado o duplicado')
        seen.add(k)
        for field,value in [('phase','calibration'),('family',u['family']),('problem',u['problem']),('n',511),('N',u['N']),('precond',u['preconditioner']),('threads',1),('commit',SOURCE),('start_kind','natural'),('start_id','natural'),('distance','natural'),('seed',0),('instance_id',iid)]:
            if r[field]!=str(value): raise ValueError('Fila suplementaria no permitida: '+field)
        if float(r['rho_x'])!=0 or float(r['param'])!=pars['param']: raise ValueError('Parámetro suplementario no permitido')
        for field,option in [('eta_init','--eta-init'),('restart','--restart'),('rtol','--rtol'),('max_nl','--max-nl'),('maxlin','--maxlin'),('eta_max','--eta-max'),('gamma','--gamma'),('alpha','--alpha'),('eta2_factor','--eta2-factor'),('eta2_max','--eta2-max')]:
            if float(r[field])!=float(flag(option)): raise ValueError('Configuración suplementaria distinta: '+field)
        if r['success'] not in ['0','1']: raise ValueError('success inválido')
    if seen!=expected or len(data)!=21: raise ValueError('Unidad suplementaria incompleta')
    return data


def supplement_gate(plan):
    assert_not_started();validate_plan(plan);git_protected()
    stub=dict(hardware=plan['hardware'],source_commit=SOURCE,units=plan['units'],starts={})
    reasons=ops.environment_reasons(stub)
    executable=ops.resolve(plan['units'][0]['argv'][0])
    if not executable.is_file() or ops.digest(executable)!=plan['executable']['sha256']:
        reasons.append('Binario run_campaign ausente o distinto del freeze')
    if reasons: raise ValueError('; '.join(reasons))


def run_supplement(execute=False):
    plan=ops.load(SUPPLEMENT);validate_plan(plan)
    if not execute:
        print(json.dumps({'invocations':48,'B3a':36,'Bext':12,'measured_rows':1008,'solver_executed':False},indent=2));return
    supplement_gate(plan)
    if (AMEND/'primitive_costs_v2.csv').exists(): raise ValueError('Suplemento ya congelado, no se ejecuta nuevamente')
    runtime=AMEND/'supplement';runtime.mkdir(exist_ok=True)
    with (runtime/'execution.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        ledger_path=runtime/'ledger.json';fingerprint=hashlib.sha256(ops.serial(plan).encode()).hexdigest()
        if not ledger_path.exists():
            if (runtime/'raw').exists() and any((runtime/'raw').iterdir()): raise ValueError('Raw suplementario previo sin ledger')
            ops.save(ledger_path,dict(plan_sha256=fingerprint,units={u['command_id']:{'state':'not_started','attempt':1} for u in plan['units']}))
        ledger=ops.load(ledger_path)
        if ledger['plan_sha256']!=fingerprint: raise ValueError('Plan del suplemento cambió')
        for u in plan['units']:
            cid=u['command_id'];e=ledger['units'][cid];path=ops.resolve(u['output_path']);receipt=path.with_suffix('.receipt.json')
            if e['state'].startswith('completed_'):
                if ops.digest(path)!=e['raw_sha256']: raise ValueError('Raw suplementario completo alterado')
                continue
            if e['state'] in ('running','infrastructure_failure'):
                if receipt.exists() and ops.load(receipt).get('returncode')==0:
                    data=validate_unit(plan,u,path)
                    e.update(state='completed_success' if all(r['success']=='1' for r in data) else 'completed_solver_failure',raw_sha256=ops.digest(path));ops.save(ledger_path,ledger);continue
                archive=runtime/'infrastructure_attempts'/cid/str(e['attempt']);archive.mkdir(parents=True,exist_ok=False)
                for p in [path,path.with_suffix('.stdout'),path.with_suffix('.stderr'),receipt]:
                    if p.exists(): os.replace(p,archive/p.name)
                ops.save(archive/'ledger_entry.json',e)
                e.update(state='not_started',attempt=e['attempt']+1);ops.save(ledger_path,ledger)
            if path.exists(): raise ValueError('Archivo raw previo inesperado')
            try:
                path.parent.mkdir(parents=True,exist_ok=True)
                e.update(state='running');ops.save(ledger_path,ledger)
                if ops.digest(ops.resolve(u['argv'][0]))!=plan['executable']['sha256']: raise ValueError('Binario congelado alterado antes de la unidad')
                with path.with_suffix('.stdout').open('xb') as out,path.with_suffix('.stderr').open('xb') as err:
                    result=subprocess.run([str(ops.resolve(u['argv'][0])),*u['argv'][1:]],cwd=ROOT,stdout=out,stderr=err)
                ops.save(receipt,dict(command_id=cid,plan_sha256=fingerprint,argv=u['argv'],warmup=u['warmup'],threads=u['threads'],environment={k:os.environ.get(k) for k in [*plan['hardware']['thread_environment'],'NUMA_COMMIT']},returncode=result.returncode))
                if result.returncode: raise ValueError('Fallo de infraestructura: código '+str(result.returncode))
                data=validate_unit(plan,u,path)
                e.update(state='completed_success' if all(r['success']=='1' for r in data) else 'completed_solver_failure',raw_sha256=ops.digest(path));ops.save(ledger_path,ledger)
            except BaseException as exc:
                e.update(state='infrastructure_failure',reason=str(exc));ops.save(ledger_path,ledger);raise


def collected_supplement(plan):
    validate_plan(plan)
    lp=AMEND/'supplement/ledger.json';ledger=ops.load(lp)
    fingerprint=hashlib.sha256(ops.serial(plan).encode()).hexdigest()
    if ledger['plan_sha256']!=fingerprint: raise ValueError('Ledger de otro suplemento')
    allrows=[];hashes={str(lp.relative_to(ROOT)):ops.digest(lp)}
    for u in plan['units']:
        path=ops.resolve(u['output_path']);entry=ledger['units'][u['command_id']]
        if entry['state'] not in ['completed_success','completed_solver_failure']: raise ValueError('Suplemento incompleto')
        if ops.digest(path)!=entry['raw_sha256']: raise ValueError('Raw suplementario alterado')
        receipt=path.with_suffix('.receipt.json');r=ops.load(receipt)
        if r['command_id']!=u['command_id'] or r['returncode']!=0 or r['argv']!=u['argv'] or r['plan_sha256']!=fingerprint or r['warmup']!=1 or r['threads']!=1: raise ValueError('Recibo o warmup del suplemento inválido')
        if r.get('environment')!={**plan['hardware']['thread_environment'],'NUMA_COMMIT':SOURCE}: raise ValueError('Entorno del suplemento inválido')
        for raw in validate_unit(plan,u,path):
            raw['internal_candidate_id']=u['candidate']['id'];raw['calibration_instance_id']=u['command_id'];allrows.append(raw)
        for p in [path,receipt,path.with_suffix('.stdout'),path.with_suffix('.stderr')]: hashes[str(p.relative_to(ROOT))]=ops.digest(p)
    if len(allrows)!=1008: raise ValueError('Se requieren exactamente 1008 filas suplementarias')
    return allrows,hashes


def write_csv(path,data,fields):
    with Path(path).open('x',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,lineterminator='\n');w.writeheader();w.writerows(data)


def cost_key(r): return tuple(r[k] for k in COST_FIELDS)


def expected_cost_keys():
    originals={cost_key(r) for r in ops.rows(ROOT/ops.FROZEN/'primitive_costs.csv')}
    return originals|{('burgers','B3a','511','522242','mg','1'),('brusselator','Bext','511','1022','jacobi','1')}


def verify_cost_table(path):
    data=ops.rows(path);keys=[cost_key(r) for r in data]
    if len(data)!=14 or len(set(keys))!=14 or set(keys)!=expected_cost_keys(): raise ValueError('Se requieren las doce claves originales más las dos n=511')
    original={cost_key(r):r for r in ops.rows(ROOT/ops.FROZEN/'primitive_costs.csv')}
    for r in data:
        if cost_key(r) in original and r!=original[cost_key(r)]: raise ValueError('Se alteraron costos originales')
    return data


def freeze_costs(input_csv,output_csv):
    subprocess.run([sys.executable,str(ROOT/'scripts/freeze_primitive_costs.py'),str(input_csv),str(output_csv)],check=True,cwd=ROOT,stdout=subprocess.PIPE,stderr=subprocess.PIPE)


def finalize():
    plan=ops.load(SUPPLEMENT);supplement_gate(plan)
    supplemental,hashes=collected_supplement(plan)
    if (AMEND/'primitive_costs_v2.csv').exists(): raise ValueError('La tabla v2 ya está congelada, no se sobrescribe')
    original=ops.rows(ops.resolve(plan['original_raw']))
    if len(original)!=6048 or any(r['phase']!='calibration' for r in original): raise ValueError('Calibration original inválida')
    fields=list(original[0]);combined=original+supplemental
    with tempfile.TemporaryDirectory(dir=AMEND) as tmp:
        tmp=Path(tmp);inp=tmp/'combined.csv';costs=tmp/'primitive_costs_v2.csv'
        write_csv(inp,combined,fields);freeze_costs(inp,costs);verify_cost_table(costs)
        # Se verifica primero todo el cálculo en temporal, después se publica v2.
        write_csv(AMEND/'supplemental_calibration.csv',supplemental,fields)
        write_csv(AMEND/'calibration_combined.csv',combined,fields)
        with (AMEND/'primitive_costs_v2.csv').open('xb') as f: f.write(costs.read_bytes())
    provenance=dict(artifact='C1_PRIMITIVE_COSTS_V2_PROVENANCE',source_commit=SOURCE,phase='calibration',confirmation_started=False,no_confirmatory_data_observed=True,
                    original=ref(ops.resolve(plan['original_raw'])),supplemental=ref(AMEND/'supplemental_calibration.csv'),combined=ref(AMEND/'calibration_combined.csv'),
                    costs=ref(AMEND/'primitive_costs_v2.csv'),generator=ref(ROOT/'scripts/freeze_primitive_costs.py'),plan=ref(SUPPLEMENT),
                    unit_artifacts=hashes,original_rows=6048,supplemental_rows=1008,total_rows=7056,retuning=False)
    ops.save(AMEND/'primitive_costs_provenance.json',provenance)
    result=initial_manifests()
    for name,d in result.items():
        d['status']='FROZEN_BEFORE_C1_AMENDED'
        d['primitive_costs_provenance']=ref(AMEND/'primitive_costs_provenance.json')
        if name=='c1_readiness_manifest_v2.json': d['artifacts']['primitive_costs']=ref(AMEND/'primitive_costs_v2.csv')
        else: d['primitive_costs']=ref(AMEND/'primitive_costs_v2.csv')
    campaign=result['campaign_manifest_v2.json'];ops.save(AMEND/'campaign_manifest_v2.json',campaign)
    result['preregistration_v2.json']['campaign_manifest']=ref(AMEND/'campaign_manifest_v2.json')
    ops.save(AMEND/'preregistration_v2.json',result['preregistration_v2.json'])
    readiness=result['c1_readiness_manifest_v2.json']
    readiness['artifacts']['campaign_manifest']=ref(AMEND/'campaign_manifest_v2.json')
    readiness['artifacts']['preregistration']=ref(AMEND/'preregistration_v2.json')
    ops.save(AMEND/'c1_readiness_manifest_v2.json',readiness)
    # Única transición autorizada del plan pendiente al plan con costos v2 reales.
    ops.save(ops.PLAN,ops.build_plan())
    check_readiness()
    print('C1_AMENDED_READINESS=PASS\nREADY_FOR_C1_AMENDED')


def verify_amendment(plan=None):
    """Verificación de datos y cadena v2, reutilizable también después de iniciar C1."""
    if ops.load(AMEND/'instance_id_policy.json')!=policy(): raise ValueError('Política instance_id no coincide con el amendment autorizado')
    manifests={name:ops.load(AMEND/(name+'_v2.json')) for name in ['campaign_manifest','preregistration','c1_readiness_manifest']}
    expected=initial_manifests()
    for name,d in manifests.items():
        e=expected[name+'_v2.json'];e['status']='FROZEN_BEFORE_C1_AMENDED'
        e['primitive_costs_provenance']=ref(AMEND/'primitive_costs_provenance.json')
        if name=='c1_readiness_manifest':
            e['artifacts']['primitive_costs']=ref(AMEND/'primitive_costs_v2.csv')
            e['artifacts']['campaign_manifest']=ref(AMEND/'campaign_manifest_v2.json')
            e['artifacts']['preregistration']=ref(AMEND/'preregistration_v2.json')
        else: e['primitive_costs']=ref(AMEND/'primitive_costs_v2.csv')
        if name=='preregistration': e['campaign_manifest']=ref(AMEND/'campaign_manifest_v2.json')
        if d!=e: raise ValueError('Manifest v2 incompleto, pendiente o alterado: '+name)
    verify_cost_table(AMEND/'primitive_costs_v2.csv')
    supplement=ops.load(SUPPLEMENT);supplemental,hashes=collected_supplement(supplement)
    p=ops.load(AMEND/'primitive_costs_provenance.json')
    if p['source_commit']!=SOURCE or p['confirmation_started'] is not False or p['no_confirmatory_data_observed'] is not True or p['retuning'] is not False: raise ValueError('Estado preconfirmatorio de procedencia inválido')
    if p['phase']!='calibration' or p['unit_artifacts']!=hashes or p['supplemental_rows']!=1008 or p['original_rows']!=6048 or p['total_rows']!=7056: raise ValueError('Procedencia de costos v2 incorrecta')
    for key in ['original','supplemental','combined','costs','generator','plan']:
        a=p[key]
        if ops.digest(ops.resolve(a['path']))!=a['sha256']: raise ValueError('Procedencia alterada: '+key)
    if p['generator']!=ref(ROOT/'scripts/freeze_primitive_costs.py') or p['plan']!=ref(SUPPLEMENT) or p['original']!=ref(ops.resolve(supplement['original_raw'])): raise ValueError('Fuentes de costos distintas')
    if ops.rows(AMEND/'supplemental_calibration.csv')!=supplemental: raise ValueError('Suplemento agregado distinto de raw')
    original=ops.rows(ops.resolve(supplement['original_raw']))
    combined=ops.rows(AMEND/'calibration_combined.csv')
    if combined!=original+supplemental or any(r['phase']!='calibration' for r in combined): raise ValueError('Combinación de calibration distinta')
    with tempfile.TemporaryDirectory() as d:
        output=Path(d)/'costs.csv';freeze_costs(AMEND/'calibration_combined.csv',output)
        if output.read_bytes()!=(AMEND/'primitive_costs_v2.csv').read_bytes(): raise ValueError('Costos v2 no reproducen freeze_primitive_costs.py')
    if plan is not None:
        if plan['costs']!='experiments/c1/amendment01/primitive_costs_v2.csv': raise ValueError('Plan no usa costos v2')
        keys={cost_key(r) for r in ops.rows(AMEND/'primitive_costs_v2.csv')}
        for u in plan['units']:
            N=u['n']**3 if u['family']=='B1' else 2*u['n']**2 if u['family']=='B3a' else u['n'] if u['family']=='B5' else 2*u['n']
            key=(u['problem'],u['family'],str(u['n']),str(N),u['preconditioner'],str(u['threads']))
            if key not in keys: raise ValueError('Clave confirmatoria sin costo: '+repr(key))


def check_readiness():
    assert_not_started();ops.sources()
    if not (AMEND/'primitive_costs_v2.csv').is_file(): raise ValueError('Pendiente: primitive_costs_v2.csv requiere las 1008 filas reales de calibration suplementaria en la workstation')
    git_protected()
    plan=ops.load(ops.PLAN);ops.audit(plan);verify_amendment(plan)
    reasons=ops.environment_reasons(dict(hardware=plan['hardware'],source_commit=SOURCE,units=[],starts={}))
    if reasons: raise ValueError('; '.join(reasons))


def main(mode=None):
    parser=argparse.ArgumentParser(description=__doc__)
    if mode is None: parser.add_argument('mode',choices=['prepare','supplement','finalize','readiness'])
    parser.add_argument('--dry-run',action='store_true');parser.add_argument('--execute',action='store_true')
    args=parser.parse_args();mode=mode or args.mode
    try:
        if args.execute and args.dry_run: raise ValueError('Modos incompatibles')
        if mode=='prepare': prepare()
        elif mode=='supplement': run_supplement(args.execute and not args.dry_run)
        elif mode=='finalize':
            if args.dry_run: print('Pendiente de 1008 filas reales del suplemento, no se escriben costos')
            elif not args.execute: raise ValueError('Para congelar localmente use --execute')
            else: finalize()
        else:
            check_readiness();print('C1_AMENDED_READINESS=PASS\nREADY_FOR_C1_AMENDED')
    except Exception as exc:
        if mode=='readiness': print('C1_AMENDED_READINESS=BLOCKED')
        print(str(exc),file=sys.stderr);sys.exit(1)

if __name__=='__main__': main()
