#!/usr/bin/env python3
"""High-precision recomputation of every numerical constant and scalar side condition of the large-gain proof.

Two independent scripts are executed and their individual comparisons collected:
  lib/lg_constants_author.py    the constants script written with the proof (mpmath, 40 significant digits)
  lib/lg_constants_referee.py   the adversarial referee's script, written from the text of the proof only
                                (mpmath, 50 significant digits); it also checks closed forms against quadrature,
                                the monotonicity claims that reduce every side condition to lambda = K_0 on grids,
                                and the bootstrap exclusion on a (mu, k, z) grid.
Each comparison is a strict or non-strict inequality between two numbers evaluated at >= 40 digits; the smallest
relative margin is reported.  mpmath is not interval arithmetic, so the receipt carries not_interval_certified;
the margins exceed the working precision by more than thirty orders of magnitude, except for the two equality
cases H = H_hat and L_0 = L_hat at mu = 1, k = K_0, which the referee script compares with a 1e-40 allowance.

--override-c1 / --override-cp replace one constant of the proof by a wrong value; used only by v3_controls, where
the run must then fail.
"""
from __future__ import annotations
import argparse, sys
from pathlib import Path
import mpmath as mp
from common import LIB, Clock, write_receipt, exec_script

HEADLINE = [  # (name in the author script's namespace, description)
    ('om0', 'omega_0(k_0): lower bound for the rotation frequency at the equilibrium'),
    ('T0', 'T_0(k_0) = 2 pi - pi/omega_0: lower bound for the length of the negative lobe'),
    ('Ri', 'Stage 1, case (i): lower bound for the rise time'),
    ('Rii', 'Stage 1, case (ii): lower bound for the rise time'),
    ('Ht', 'logarithmic majorant of mu H at k_0'),
    ('eps1', 'epsilon_1-bar at k_0'),
    ('d1', 'loss in the spike-mass lemma with p_k >= c_1 k, at k_0'),
    ('T1', 'T_1(k_0): first lower bound for T_-'),
    ('z10', 'z_1(k_0) = 0.3 T_1(k_0) - 3'),
    ('kap0', 'the constant of the bootstrap (kappa_0 in the research text)'),
    ('Hhat', 'H-hat at lambda = k_0'),
    ('Lhat', 'L-hat at lambda = k_0'),
    ('e1', 'epsilon_1-hat at lambda = k_0'),
    ('dp', "delta' at lambda = k_0"),
    ('ac', 'a_c at lambda = k_0'),
    ('A1', 'A_1 at lambda = k_0'),
    ('tau1', 'tau_1-bar at lambda = k_0'),
    ('W', 'W-hat at lambda = k_0'),
    ('vf', 'v_f at lambda = k_0'),
    ('tauc', 'crest-duration majorant at lambda = k_0'),
    ('Ctau', 'C_tau(k_0)'),
    ('cb_log', 'coefficient of log(lambda) in C_tau'),
    ('q', 'q(k_0)'),
    ('Phi_', 'final margin (c_T/4) k_0 - 2 Cbar_tau/(c_p sqrt(k_0)) - log(1 + 4 sqrt(k_0) + 16 k_0/3)'),
]
DERIVED = {   # further numbers of the constants table, evaluated in the namespace of the author script
    'k0_ci_a1': ('k_0 c_i(a_1)', lambda n: n['K0'] * n['ci'](n['A1W'])),
    'chi_ii_a2': ('k_0 (a_2 - 1 + exp(-a_2))/(k_0^(1/3) + 1)', lambda n: n['K0'] * (n['A2W'] - 1 + n['exp'](-n['A2W'])) / (n['K0'] ** (n['mpf'](1) / 3) + 1)),
    'k0_cb_a4': ('k_0 c_b(a_4)', lambda n: n['K0'] * n['cb'](n['A4W'])),
    'strip_ratio': ('ell/phi_ell', lambda n: n['ELL'] / n['PHI']),
    'phi_ell': ('phi_ell = 1 - exp(-2)', lambda n: n['PHI']),
    'c1_bound': ('phi_ell Psi_1(a_3) - ell/k_0, the bound for c_1', lambda n: n['PHI'] * n['Psi1'](n['A3W']) - n['ELL'] / n['K0']),
    'cp_bound': ('0.95 phi_ell (1 - exp(-c_sharp k_0)) - ell/k_0, the bound for c_p', lambda n: n['mpf']('0.95') * n['PHI'] * (1 - n['exp'](-n['CSHARP'] * n['K0'])) - n['ELL'] / n['K0']),
    'cT': ('c_T = c_sharp/0.3', lambda n: n['CT']),
    'window_sum': ('a_1 + a_2 + a_3 + a_4', lambda n: n['A1W'] + n['A2W'] + n['A3W'] + n['A4W']),
    'margin_term_1': ('(c_T/4) k_0', lambda n: n['CT'] * n['K0'] / 4),
    'margin_term_2': ('2 Cbar_tau(k_0)/(c_p sqrt(k_0))', lambda n: 2 * n['Cbar'] / (n['CP'] * n['sqrt'](n['K0']))),
    'margin_term_3': ('log(1 + 4 sqrt(k_0) + 16 k_0/3)', lambda n: n['logterm']),
}


def run(override=None, name='v3_largegain_constants', out=None, logs=True):
    clock = Clock()
    recs = []
    init = dict(override or {})
    a_ns, a_text, a_status = exec_script(LIB / 'lg_constants_author.py', init=init,
                                         log_name=(name + '__author.out') if logs else None)
    for r in a_ns.get('RECORDS', []):
        lhs, rhs = mp.mpf(r['lhs']), mp.mpf(r['rhs'])
        margin = (lhs - rhs) if r['op'] == '>=' else (rhs - lhs)
        scale = max(abs(lhs), abs(rhs), mp.mpf(1))
        recs.append(dict(source='author script', id=r['name'], lhs=r['lhs'], op=r['op'], rhs=r['rhs'],
                         margin=r['margin'], relative_margin=float(margin / scale), passed=bool(r['passed'])))
    r_ns, r_text, r_status = ({}, '', 0)
    if not override:
        r_ns, r_text, r_status = exec_script(LIB / 'lg_constants_referee.py', log_name=(name + '__referee.out') if logs else None)
        for r in r_ns.get('RECORDS', []):
            recs.append(dict(source='referee script', id=r['name'], detail=r['detail'][:300], passed=bool(r['passed'])))
    headline = {}
    for key, desc in HEADLINE:
        if key in a_ns:
            headline[key] = dict(value=mp.nstr(a_ns[key], 12), description=desc)
    for key, (desc, fn) in DERIVED.items():
        try:
            headline[key] = dict(value=mp.nstr(fn(a_ns), 12), description=desc)
        except Exception as exc:          # a missing name is recorded, never silently dropped
            headline[key] = dict(value=None, description=desc, error=repr(exc))
    rel = [r['relative_margin'] for r in recs if 'relative_margin' in r]
    worst = {'min_relative_margin_author_script': min(rel) if rel else None,
             'tightest_author_comparison': min((r for r in recs if 'relative_margin' in r), key=lambda r: r['relative_margin'])['id'] if rel else None,
             'author_script_exit_status': a_status, 'referee_script_exit_status': r_status}
    extra = {'software': {'mpmath': mp.__version__},
             'precision_digits': {'author_script': 40, 'referee_script': 50},
             'headline_constants_at_K0': headline,
             'number_author_comparisons': sum(r['source'] == 'author script' for r in recs),
             'number_referee_checks': sum(r['source'] == 'referee script' for r in recs),
             'override': {k: str(v) for k, v in (override or {}).items()}}
    passed = bool(recs) and all(r['passed'] for r in recs) and a_status == 0 and r_status == 0
    return write_receipt(name, 'mp', recs, worst=worst, clock=clock, extra=extra, out=out, passed=passed,
                         parameters={'K0': 1000, 'range': '0 < mu <= 1, lambda >= K0'},
                         scope='Scalar inequalities between explicitly computable numbers. They are the only numerical input '
                               'of the large-gain proof; the proof itself is in the paper.')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--override-c1'); ap.add_argument('--override-cp'); ap.add_argument('--out', type=Path)
    a = ap.parse_args()
    ov = {}
    if a.override_c1:
        ov['OVERRIDE_C1'] = a.override_c1
    if a.override_cp:
        ov['OVERRIDE_CP'] = a.override_cp
    ok = run(ov or None, name='v3_largegain_constants' if not ov else 'v3_largegain_constants_with_wrong_constant', out=a.out,
             logs=not ov)
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
