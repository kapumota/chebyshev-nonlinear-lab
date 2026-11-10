"""Contratos del amendment, con fixtures temporales que nunca son calibración."""
import collections
import contextlib
import copy
import csv
import hashlib
import io
from pathlib import Path
import shutil
import tempfile
import types
import unittest
from unittest.mock import patch
import c1_ops as ops
import c1_amendment as amendment

class AmendmentContracts(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.plan=amendment.supplement_plan()
        cls.confirm=ops.build_plan()
    def test_supplement_exact_design(self):
        self.assertEqual(self.plan,amendment.supplement_plan())
        self.assertEqual(collections.Counter(u['family'] for u in self.plan['units']),dict(B3a=36,Bext=12))
        self.assertEqual(sum(u['expected_rows'] for u in self.plan['units']),1008)
        for u in self.plan['units']:
            self.assertEqual(u['n'],511)
            argv=u['argv'];self.assertEqual(argv[argv.index('--phase')+1],'calibration')
            self.assertEqual(argv[argv.index('--start')+1],'natural')
    def test_changed_supplement_command_rejected(self):
        plan=copy.deepcopy(self.plan);plan['units'][0]['argv'][0]='otro'
        with self.assertRaises(ValueError): amendment.validate_plan(plan)
    def test_full_precision_identity(self):
        i=self.confirm['instances'][0]
        iid=amendment.canonical_id(i)
        for k,v in i['parameters'].items(): self.assertIn(k+'='+float(v).hex(),iid)
        ids={amendment.canonical_id(i) for i in self.confirm['instances']}
        self.assertEqual(len(ids),75)
    def test_canonical_only_two_columns_and_exact_cells(self):
        i=self.confirm['instances'][0]
        raw={'method':'tf2018','instance_id':'texto interno','param':'1.2345678901234567','success':'0','t_total':'1.000000000e+00','n_F':'7','n_JVP':'0','n_D2':'0','n_prec':'0','seed':'12','start_kind':'perturbed','distance':'far','rho_x':'0.5','repetition':'7'}
        c=amendment.canonicalize_row(raw,i)
        self.assertEqual(c['method'],'tf_cheb')
        self.assertEqual({k for k in raw if raw[k]!=c[k]},{'method','instance_id'})
        for method in ['nk','cheb','rec','tf2018','petsc']:
            self.assertEqual(amendment.canonicalize_row(dict(raw,method=method),i)['instance_id'],c['instance_id'])
        altered=dict(c,t_total='1')
        with self.assertRaises(ValueError): amendment.audit_canonical_row(raw,altered,i)
    def test_pending_does_not_contain_fictitious_costs(self):
        self.assertFalse((amendment.AMEND/'primitive_costs_v2.csv').exists())
        self.assertEqual(ops.load(amendment.AMEND/'campaign_manifest_v2.json')['primitive_costs']['sha256'],None)
    def test_original_costs_reproduce_without_calibration(self):
        with tempfile.TemporaryDirectory() as d:
            output=Path(d)/'original.csv'
            amendment.freeze_costs(ops.resolve(self.plan['original_raw']),output)
            self.assertEqual(output.read_bytes(),(ops.ROOT/ops.FROZEN/'primitive_costs.csv').read_bytes())

class IsolatedSupplementContracts(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory()
        self.root=Path(self.temp.name)/'repo'
        # Copia sólo las dependencias declaradas, nunca historiales, builds ni vectores masivos.
        required,_=ops.sources()
        paths=set(required)
        paths.update(str(p.relative_to(ops.ROOT)) for p in ops.TOOLS.glob('*') if p.is_file())
        paths.update(str(p.relative_to(ops.ROOT)) for p in amendment.AMEND.glob('*') if p.is_file())
        paths.update(['experiments/c1/internal_calibration/internal_calibration_raw.csv',
                      'experiments/c1/internal_calibration/internal_calibration_plan.json',
                      'experiments/c1/internal_calibration/internal_calibration_run_manifest.json',
                      'scripts/freeze_primitive_costs.py', 'docs/research/C1_PRECONFIRMATION_AMENDMENT_01.md'])
        for name in paths:
            target=self.root/name;target.parent.mkdir(parents=True,exist_ok=True)
            shutil.copy2(ops.ROOT/name,target)
        self.stack=contextlib.ExitStack()
        tools=self.root/'experiments/c1/confirmation_tools';amend=self.root/'experiments/c1/amendment01'
        for module,attr,value in [(ops,'ROOT',self.root),(ops,'TOOLS',tools),(ops,'PLAN',tools/'c1_execution_plan.json'),(amendment,'ROOT',self.root),(amendment,'TOOLS',tools),(amendment,'AMEND',amend),(amendment,'SUPPLEMENT',amend/'supplemental_calibration_plan.json')]:
            self.stack.enter_context(patch.object(module,attr,value))
        self.plan=amendment.supplement_plan()
        self.original=ops.rows(ops.resolve(self.plan['original_raw']))
        binary=self.root/'build-petsc/run_campaign';binary.parent.mkdir(exist_ok=True);binary.write_bytes(b'fixture, no ejecutable cientifico')
    def tearDown(self):
        self.stack.close();self.temp.cleanup()
    def fixture_rows(self,u):
        data=[]
        for method in u['methods']:
            sample=next(r for r in self.original if r['family']==u['family'] and r['method']==method and r['internal_candidate_id']==u['candidate']['id'])
            for rep in range(1,8):
                r={k:v for k,v in sample.items() if k not in ['internal_candidate_id','calibration_instance_id']}
                p=u['parameters'];r.update(n='511',N=str(u['N']),param=format(p['param'],'.17g'),repetition=str(rep))
                if u['family']=='B3a': r['instance_id']='nu='+format(p['param'],'.6f')+';dt='+format(p['dt'],'.6f')+';amp='+format(p['amp'],'.6f')
                else: r['instance_id']='dt='+format(p['param'],'.6f')+';A=1;B=3;alpha=0.02'
                data.append(r)
        return data
    def fixture_run(self,argv,**kwargs):
        path=ops.resolve(argv[argv.index('--out')+1])
        u=next(u for u in self.plan['units'] if ops.resolve(u['output_path'])==path)
        data=self.fixture_rows(u)
        with path.open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(data[0]),lineterminator='\n');w.writeheader();w.writerows(data)
        return types.SimpleNamespace(returncode=0)
    def execute_fixtures(self):
        original_digest=ops.digest
        def fixture_digest(path):
            if Path(path)==self.root/'build-petsc/run_campaign': return self.plan['executable']['sha256']
            return original_digest(path)
        env={**self.plan['hardware']['thread_environment'],'NUMA_COMMIT':amendment.SOURCE}
        self.stack.enter_context(patch.dict(__import__('os').environ,env))
        self.stack.enter_context(patch.object(ops,'digest',side_effect=fixture_digest))
        with patch.object(amendment,'supplement_gate'),patch('subprocess.run',side_effect=self.fixture_run) as solver:
            amendment.run_supplement(True)
            self.assertEqual(solver.call_count,48)
        with patch.object(amendment,'supplement_gate'),patch('subprocess.run') as solver:
            amendment.run_supplement(True);solver.assert_not_called()
    def test_complete_units_are_not_repeated(self): self.execute_fixtures()
    def test_wrong_phase_and_candidate_rejected(self):
        u=self.plan['units'][0];path=self.root/'test_fixture.csv';data=self.fixture_rows(u)
        def write():
            with path.open('w',newline='') as f:
                w=csv.DictWriter(f,fieldnames=list(data[0]));w.writeheader();w.writerows(data)
        data[0]['phase']='confirmation';write()
        with self.assertRaises(ValueError): amendment.validate_unit(self.plan,u,path)
        data[0]['phase']='calibration';data[0]['eta_init']='0.123';write()
        with self.assertRaises(ValueError): amendment.validate_unit(self.plan,u,path)
    def test_warmup_receipt_tampering_rejected(self):
        self.execute_fixtures()
        path=ops.resolve(self.plan['units'][0]['output_path']).with_suffix('.receipt.json')
        receipt=ops.load(path);receipt['warmup']=0;ops.save(path,receipt)
        with self.assertRaises(ValueError): amendment.collected_supplement(self.plan)
    def test_freeze_pipeline_in_temporary_fixture(self):
        self.execute_fixtures()
        with patch.object(amendment,'supplement_gate'),patch.object(amendment,'check_readiness'),contextlib.redirect_stdout(io.StringIO()):
            amendment.finalize()
        amendment.verify_amendment(ops.load(ops.PLAN))
        ops.audit(ops.load(ops.PLAN))
        costs=ops.rows(amendment.AMEND/'primitive_costs_v2.csv');self.assertEqual(len(costs),14)
        self.assertEqual({amendment.cost_key(r) for r in costs},amendment.expected_cost_keys())
        with patch.object(amendment,'supplement_gate'),self.assertRaises(ValueError): amendment.finalize()
        # Toda fixture y costo de esta prueba desaparece al salir del TemporaryDirectory.

if __name__=='__main__': unittest.main()
