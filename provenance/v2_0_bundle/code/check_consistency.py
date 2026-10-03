#!/usr/bin/env python3
"""Check receipt counts, provenance, saved payloads and generated-table inputs.
This verifies bookkeeping and finite output predicates, not the mathematical proof.
"""
from __future__ import annotations
import hashlib,json,sys,zipfile
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):
    def bad(x):raise ValueError('Non-standard JSON constant '+x)
    return json.loads(p.read_text(),parse_constant=bad)

def run():
    checks=[]
    def ck(name,ok,detail=None):
        checks.append({'id':name,'passed':bool(ok),'detail':detail})
    try:
        exact=read(ROOT/'evidence/exact_checks.json')
        ck('exact_count',len(exact['checks'])==exact['number_checked']==50)
        ck('exact_actual_equalities',all(row['passed'] and row['residual'] in ('0','Matrix([[0, 0, 0], [0, 0, 0], [0, 0, 0]])','Matrix([[0], [0], [0]])','Matrix([[0, 0, 0]])','Matrix([[0], [0]])') for row in exact['checks']), 'Residual strings are checked as well as pass flags.')
        for name,key,n in [('numerical_audit','records',31),('cross_method','records',4),('reciprocal_rates','records',12),('envelope_checks','records',4),('negative_controls','controls',4)]:
            d=read(ROOT/'evidence'/f'{name}.json')
            ck(name+'_actual_count',len(d[key])==n and d.get('number_checked',d.get('number_attempted'))==n)
            ck(name+'_all_predicates',d['passed'] and all(x['passed'] for x in d[key]) and d['number_passed']==n)
        trials=read(ROOT/'evidence/cycle_trials.json');rows=trials['attempts']
        ck('BVP_count',len(rows)==31 and len({q['id'] for q in rows})==31 and trials['number_attempted']==31)
        ck('BVP_acceptance',trials['number_accepted']==31 and all(q['solver_success'] and q['accepted'] and all(q['acceptance_checks'].values()) and not q['guard_triggered'] for q in rows))
        for q in rows:
            p=ROOT/q['mesh']
            with np.load(p,allow_pickle=False) as a:
                ok=all(np.all(np.isfinite(a[k])) for k in a.files)
            ck('mesh_finite_'+q['id'],ok)
        audits=read(ROOT/'evidence/numerical_audit.json')['records']
        for q in audits:
            with np.load(ROOT/q['profile'],allow_pickle=False) as a:ok=all(np.all(np.isfinite(a[k])) for k in a.files)
            ck('profile_finite_'+q['id'],ok and all(q['checks'].values()))
        summary=read(ROOT/'evidence/manuscript_summary.json')
        for name,h in summary['source_sha256'].items():ck('table_source_'+name,sha(ROOT/'evidence'/name)==h)
        rr=read(ROOT/'evidence/run_receipt.json')
        ck('complete_active_run_record',rr['passed'] and all(x['exit_status']==0 for x in rr['steps']), {'regenerated_cycles':rr['regenerated_cycles']})
        for name,h in rr['code_sha256'].items():ck('executed_code_'+name,sha(ROOT/name)==h)
        inputs=read(ROOT/'provenance/INPUTS.json')
        for row in inputs['inputs']:ck('unchanged_input_'+row['snapshot'],sha(ROOT/row['snapshot'])==row['sha256'])
        with zipfile.ZipFile(ROOT/'provenance/earlier_review/v0_3_referee_checks.zip') as z:ck('historical_review_zip_crc',z.testzip() is None)
        old=read(ROOT/'evidence/exact_checks_initial_failure.json');ck('real_checker_failure_retained',not old['passed'] and old['number_passed']==49 and old['number_checked']==50)
        wrong=read(ROOT/'evidence/deliberately_wrong_exact_receipt.json');ck('deliberate_wrong_claim_rejected',not wrong['passed'] and wrong['number_passed']<wrong['number_checked'])
        for name,n in [('neutral_cycle_search.json',24),('fixed_rate_neutral_search.json',18)]:
            d=read(ROOT/'provenance/v1_2/evidence'/name);ck('historical_failed_search_'+name,len(d['rows'])==n and all(not q['accepted_neutral_candidate'] for q in d['rows']))
        status=read(ROOT/'STATUS.json');ck('claim_boundary',not status['original_all_parameter_conjecture_solved'] and not status['current_independent_referee_review'])
        ck('no_font_files',not any(p.suffix.lower() in ('.ttf','.otf','.woff','.woff2') for p in ROOT.rglob('*') if p.is_file()))
    except Exception as exc:ck('exception',False,repr(exc))
    result={'kind':'release_bookkeeping_and_finite_receipt_consistency','checks':checks,'number_checked':len(checks),'number_passed':sum(x['passed'] for x in checks),'passed':all(x['passed'] for x in checks),'scope':'Not proof verification and not a numerical enclosure.'}
    (ROOT/'evidence/release_checks.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(f"{result['number_passed']}/{result['number_checked']} consistency checks passed")
    for q in checks:
        if not q['passed']:print(q)
    return result['passed']
if __name__=='__main__':sys.exit(0 if run() else 1)
