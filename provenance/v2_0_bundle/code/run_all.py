#!/usr/bin/env python3
"""Reproduce the active evidence; ordinary numerics, never a formal proof.
Run on a copy if the original SHA-256 manifest must remain unchanged.
--regenerate recomputes all 31 amplitude-held BVPs; otherwise their saved meshes
are reused, and all dynamical/adjoint/envelope checks are re-executed.
"""
from __future__ import annotations
import argparse,hashlib,json,os,platform,subprocess,sys,time
from pathlib import Path
from datetime import datetime,timezone
ROOT=Path(__file__).resolve().parents[1]

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def checked_json(path):
    def bad(value):raise ValueError(f'Non-standard JSON value: {value}')
    return json.loads(path.read_text(),parse_constant=bad)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--regenerate',action='store_true');args=ap.parse_args()
    (ROOT/'logs').mkdir(exist_ok=True)
    started=datetime.now(timezone.utc).isoformat()
    receipt={'schema_version':'1.0','started_utc':started,'python':sys.version,'platform':platform.platform(),'regenerated_cycles':args.regenerate,'code_sha256':{str(p.relative_to(ROOT)):digest(p) for p in sorted((ROOT/'code').glob('*.py'))},'steps':[],'passed':False}
    steps=[('exact_checks','exact_checks.py')]
    if args.regenerate:steps.append(('cycle_regeneration','cycles.py'))
    elif not (ROOT/'evidence/cycle_trials.json').exists():raise FileNotFoundError('Missing cycle trials: use --regenerate.')
    steps.extend([('numerical_audit','numerical_audit.py'),('additional_checks','additional_checks.py'),('manuscript_tables','make_tables.py')])
    env=os.environ.copy();env.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
    try:
        for name,script in steps:
            t=time.monotonic();command=[sys.executable,str(ROOT/'code'/script)]
            log=ROOT/'logs'/f'{name}_rerun.log'
            with log.open('w') as f:result=subprocess.run(command,cwd=ROOT,env=env,stdout=f,stderr=subprocess.STDOUT)
            row={'name':name,'script':script,'exit_status':result.returncode,'elapsed_seconds':time.monotonic()-t,'log':str(log.relative_to(ROOT))}
            receipt['steps'].append(row);print(name,'exit',result.returncode,flush=True)
            if result.returncode:raise RuntimeError(f'{name} failed; inspect {log}')
        trial=checked_json(ROOT/'evidence/cycle_trials.json')
        expected={f'mu_{m:.8g}_h_{h:.8g}' for m in [.1,.3,1.,10/3,10.] for h in [.05,.2,.5,1.,2.,3.]}|{'mu_0.3_h_0.01'}
        actual={r['id'] for r in trial['attempts']}
        if actual!=expected or len(trial['attempts'])!=31 or not all(r['accepted'] for r in trial['attempts']):raise RuntimeError('Expected 31 accepted target approximations not reproduced; failed attempts remain in the receipt.')
        for name in ['exact_checks','numerical_audit','cross_method','reciprocal_rates','envelope_checks','negative_controls','manuscript_summary']:
            data=checked_json(ROOT/'evidence'/f'{name}.json')
            if not data.get('passed',False):raise RuntimeError(f'Receipt {name} did not pass its scientific checks.')
        # The deliberately bad claim must fail, rather than be counted as evidence.
        wrong=checked_json(ROOT/'evidence/deliberately_wrong_exact_receipt.json')
        if wrong['passed']:raise RuntimeError('Negative control unexpectedly passed the wrong assertion.')
        receipt['passed']=True
        receipt['scope']='All specified finite checks reproduced. This is not a proof assistant or an interval certificate.'
    except Exception as exc:receipt['error']=repr(exc)
    receipt['completed_utc']=datetime.now(timezone.utc).isoformat()
    (ROOT/'evidence/run_receipt.json').write_text(json.dumps(receipt,indent=2,allow_nan=False)+'\n')
    return 0 if receipt['passed'] else 1
if __name__=='__main__':sys.exit(main())
