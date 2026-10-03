#!/usr/bin/env python3
"""Generate manuscript numbers directly from strict JSON evidence receipts."""
from __future__ import annotations
import hashlib,json,sys
from pathlib import Path
import numpy as np
from cycles import ROOT,strict_write

def read(n): return json.loads((ROOT/'evidence'/n).read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def run():
    exact=read('exact_checks.json');trial=read('cycle_trials.json');num=read('numerical_audit.json');cross=read('cross_method.json');rec=read('reciprocal_rates.json');env=read('envelope_checks.json');neg=read('negative_controls.json')
    rows=num['records'];by={r['id']:r for r in rows}
    summary={'schema_version':'1.0','claim_boundary':'Finite numerical diagnostics, not interval certificates and not exhaustive cycle exclusion.',
       'counts':{'exact_equalities':[exact['number_passed'],exact['number_checked']], 'cycle_approximations':[trial['number_accepted'],trial['number_attempted']], 'cycle_geometry_and_adjoint_audits':[num['number_passed'],num['number_attempted']], 'cross_method':[cross['number_passed'],cross['number_checked']], 'reciprocal_rates':[rec['number_passed'],rec['number_checked']], 'envelope':[env['number_passed'],env['number_checked']], 'negative_controls':[neg['number_passed'],neg['number_checked']]},
       'maxima':{'exponent_cross_formulation':max(r['exponent_difference'] for r in rows),'orthogonality_scaled':max(r['orthogonality_max_scaled'] for r in rows),'state_return_closure':max(r['return_closure_error'] for r in rows),'moment_relative_residual':max(v for r in rows for v in r['moment_relative_residuals'].values()),'cross_method_matrix_difference':max(r['return_matrix_max_difference'] for r in cross['records']),'cross_method_exponent_difference':max(r['exponent_difference'] for r in cross['records']),'reciprocal_all_comparisons':max(v for r in rec['records'] for key,v in r.items() if key.startswith('max_difference') or key.startswith('difference_')),'envelope_hessian_scaled_error':max(r['relative_scaled_error'] for r in env['records']),'envelope_defect_residual':max(r['defect_identity_max_residual'] for r in env['records'])},
       'square_root_certificate_failures':sum(not all(v<1 for v in r['square_root_test_ratios']) for r in rows),'number_negative_computed_exponents':sum(r['lambda_c_projective']<0 for r in rows),
       'phase_example':by['mu_1_h_1']['phase_fractions'],'false_radial_example':{k:by['mu_0.3_h_0.01'][k] for k in ['mu','amplitude','square_root_test_ratios','weak_multiplier']}}
    src=['exact_checks.json','cycle_trials.json','numerical_audit.json','cross_method.json','reciprocal_rates.json','envelope_checks.json','negative_controls.json']
    summary['source_sha256']={n:sha(ROOT/'evidence'/n) for n in src}
    summary['passed']=all(exact['passed'] and v[0]==v[1] for v in summary['counts'].values())
    strict_write(ROOT/'evidence/manuscript_summary.json',summary)
    def sci(v):
        if v==0:return '0'
        s=f'{v:.3e}';a,b=s.split('e');return a+r'\times10^{'+str(int(b))+'}'
    macros=[]
    names={'exponent_cross_formulation':'ExponentAgreement','orthogonality_scaled':'OrthogonalityError','state_return_closure':'ReturnClosure','moment_relative_residual':'MomentError','cross_method_matrix_difference':'CrossMatrixError','cross_method_exponent_difference':'CrossExponentError','reciprocal_all_comparisons':'ReciprocalError','envelope_hessian_scaled_error':'EnvelopeHessianError','envelope_defect_residual':'EnvelopeDefectError'}
    for k,n in names.items():macros.append('\\newcommand{\\'+n+'}{'+sci(summary['maxima'][k])+'}')
    (ROOT/'paper/evidence_macros.tex').write_text('% Automatically generated; do not edit.\n'+'\n'.join(macros)+'\n')
    labels={'exact_equalities':'Exact symbolic equalities','cycle_approximations':'Fresh periodic boundary-value solves','cycle_geometry_and_adjoint_audits':'Cycle geometry and adjoint audits','cross_method':'DOP853--Radau repetitions','reciprocal_rates':'Reciprocal-rate coordinate comparisons','envelope':'Local envelope Hessian checks','negative_controls':'Deliberate-error controls'}
    lines=[r'\begin{tabular}{lrr}',r'\toprule',r'Check & Accepted/passed & Attempted \\',r'\midrule']
    for k,v in summary['counts'].items():lines.append(f'{labels[k]} & {v[0]} & {v[1]} '+r'\\')
    lines.extend([r'\bottomrule',r'\end{tabular}'])
    (ROOT/'paper/evidence_counts.tex').write_text('\n'.join(lines)+'\n')
    lines=[r'\begin{tabular}{rrrrr}',r'\toprule',r'$\mu$ & $h$ & $T$ (fast time) & $\lambda_c$ & $m_t$ \\',r'\midrule']
    for ident in ['mu_0.1_h_3','mu_0.3_h_0.01','mu_1_h_1','mu_10_h_3']:
        r=by[ident];lines.append(f"{r['mu']:.3g} & {r['amplitude']:.3g} & {r['period_fast']:.10f} & ${sci(r['lambda_c_projective'])}$ & {r['weak_multiplier']:.10f} "+r'\\')
    lines.extend([r'\bottomrule',r'\end{tabular}'])
    (ROOT/'paper/cycle_examples.tex').write_text('\n'.join(lines)+'\n')
    lines=[r'\begin{tabular}{lrrrrrr}',r'\toprule',r'Event & $X_{\max}$ & $Z_{\max}$ & $Y_{\max}$ & $X_{\min}$ & $Z_{\min}$ & $Y_{\min}$ \\',r'\midrule']
    lines.append('Fraction of period & '+' & '.join(f'{x:.6f}' for x in summary['phase_example'].values())+r'\\')
    lines.extend([r'\bottomrule',r'\end{tabular}'])
    (ROOT/'paper/phase_example.tex').write_text('\n'.join(lines)+'\n')
    print(json.dumps(summary,indent=2));return summary['passed']
if __name__=='__main__':sys.exit(0 if run() else 1)
