"""Pruebas estructurales, sin solvers, starts ni datos científicos nuevos."""
import copy
import csv
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import c1_ops as ops

class OperationalTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls): cls.plan=ops.build_plan()
    def test_counts_and_determinism(self):
        self.assertEqual(self.plan,ops.build_plan())
        self.assertEqual(ops.audit(self.plan)['commands'],4350)
        self.assertEqual(sum(u['expected_rows'] for u in self.plan['units']),99750)
    def test_tampered_command_is_rejected(self):
        plan=copy.deepcopy(self.plan);plan['units'][0]['argv'][1]='alterado'
        with self.assertRaises(ValueError): ops.audit(plan)
    def test_duplicate_is_rejected(self):
        plan=copy.deepcopy(self.plan);plan['units'].append(plan['units'][0])
        with self.assertRaises(ValueError): ops.audit(plan,False)
    def test_seed_omission_is_rejected(self):
        plan=copy.deepcopy(self.plan)
        u=next(u for u in plan['units'] if isinstance(u['seed'],list));u['seed']=u['seed'][:-1]
        with self.assertRaises(ValueError): ops.audit(plan,False)
    def test_pending_amendment_blocks_execution(self):
        reasons=ops.downstream_reasons(self.plan)
        self.assertTrue(reasons)
        self.assertIn('Amendment pendiente',reasons[0])
    def test_wrong_start_hash_does_not_invoke_materializer(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'vector';p.write_text('prueba de integridad')
            plan={'starts':{'x':{'path':str(p),'sha256':'0'*64,'argv':['NO_EJECUTAR'],'dependencies':[]}}}
            with patch('subprocess.run') as run:
                with self.assertRaises(ValueError): ops.materialize(plan,'x',set())
                run.assert_not_called()
    def test_solver_failures_and_duplicate_repetitions(self):
        u=next(x for x in self.plan['units'] if x['start_kind']=='natural')
        data=[]
        for rep in range(1,8):
            r={k:str(v) for k,v in dict(method=u['method'],family=u['family'],problem=u['problem'],n=u['n'],precond=u['preconditioner'],start_kind=u['start_kind'],start_id=u['start_id'],distance=u['distance'],threads=1,commit=u['source_commit'],phase='confirmation',seed=0,repetition=rep,rho_x=0,param=u['parameters']['param'],success=0).items()}
            r.update(N=str(u['n']**3 if u['family']=='B1' else 2*u['n']**2 if u['family']=='B3a' else u['n'] if u['family']=='B5' else 2*u['n']),instance_id='lambda='+format(u['parameters']['param'],'.6f'))
            for field,flag in [('rtol','--rtol'),('max_nl','--max-nl'),('restart','--restart'),('maxlin','--maxlin')]: r[field]=u['argv'][u['argv'].index(flag)+1]
            data.append(r)
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'fixture.csv'
            def write(items):
                with p.open('w',newline='') as f:
                    w=csv.DictWriter(f,fieldnames=list(data[0]));w.writeheader();w.writerows(items)
            write(data);self.assertEqual(len(ops.validate_raw(self.plan,u,p)),7)
            write(data+[data[0]])
            with self.assertRaises(ValueError): ops.validate_raw(self.plan,u,p)

if __name__=='__main__': unittest.main()
