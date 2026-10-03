#!/usr/bin/env python3
"""Floating-point tests of the weak-adjoint theorems on 20 computed cycles (wraps lib/mono_floquet.py).

Cycles: mu in {1, 0.3, 0.02} x rho in {1.01, 1.5, 4, 12}, and (mu, rho) in {(3,1.01), (3,1.5), (3,4), (3,12), (50,1.5),
(50,4), (0.1,1.2), (7,2.5)}.  The wrapped script locates each cycle in logarithmic coordinates, integrates the
strong-stable ratios backward and the adjoint ratios forward, and decides every sign in signed-log arithmetic.

Per cycle:
  twelve_event_order     each of p, p', y', Psi_1, Psi_2, Psi_3 changes sign exactly twice per period, in the cyclic order
                         of the theorem; Psi > 0 componentwise at the crest and Psi < 0 at the trough.  Two events falling in
                         one double-precision time cell are admitted only for the pair (zero of p, zero of Psi_2), where the
                         order is decided exactly from Psi_2 = varpi y' at p = 0.
  supporting_plane_chords   det(z(t+s) - z(t), z'(t), e_s(t)) < 0 on a (base, chord) grid of about 2.9e5 pairs: no positive
                         value among pairs separated by more than 1e-2 in logarithmic coordinates, no resolved positive value
                         nearer the diagonal, and no sign mismatch against a plain determinant where that is computable
  supporting_plane_E     det(E - z(t), z'(t), e_s(t)) < 0 at every base point
  supporting_plane_curvature   det(z''(t), z'(t), e_s(t)) < 0 at every base point
  supporting_line        the tangency phase t_1 lies on the arc p < 0, p' > 0, y' < 0, with U(t_1) < 0; the supporting line of
                         the projected cycle at t_1 is positive away from t_1, positive at E, with positive second derivative
  rank_one_identity      the window identity for phi in {1, p', -U, cos(omega t), sin(3 omega t)}: relative error <= 1e-4,
                         with lambda_c taken from the independent forward (eta, zeta) system
  special_signs          the three special window integrals have the signs stated in the proposition
  exponent_three_ways    <eta> - mu, <r> - A and the logarithm of varpi over the window agree to 1e-7 (1 + |lambda_c|)
  strong_stable_ratios   the periodic ratios r, s are reached (defect <= 1e-7), s > 0, and the margin r - max(1, mu),
                         integrated in logarithmic variables, agrees with the direct value where both are resolved
The test function phi = y is evaluated and reported but not counted: on deep-trough cycles with mu < 1 its right-hand side
is many orders of magnitude smaller than the two integrals whose difference it is, so the comparison is rounding noise.
Nothing here is an interval enclosure.
"""
from __future__ import annotations
import math, sys, time, warnings
import numpy as np
from common import Clock, write_receipt, capture

GRID = [(mu, rho) for mu in (1.0, 0.3, 0.02) for rho in (1.01, 1.5, 4.0, 12.0)] + \
       [(3.0, 1.01), (3.0, 1.5), (3.0, 4.0), (3.0, 12.0), (50.0, 1.5), (50.0, 4.0), (0.1, 1.2), (7.0, 2.5)]
COUNTED = ['1', "p'", '-U', 'cos', 'sin3']
# cyclic order of the twelve events, in the names of the wrapped script: q1, q2, q3 are the components of the weak adjoint,
# pp = p', yp = y'; '+' is an up-crossing.  Starting from the trough (up-crossing of p).
CLAIMED_ORDER = ['p+', 'q2+', 'yp+', 'q1+', 'pp-', 'q3+', 'p-', 'q2-', 'yp-', 'q1-', 'pp+', 'q3-']


def order_matches(seq):
    """the predicate shared with v3_controls: the event sequence, rotated to start at the trough, is the claimed one"""
    if len(seq) != 12 or 'p+' not in seq:
        return False
    i = seq.index('p+')
    return list(seq[i:]) + list(seq[:i]) == CLAIMED_ORDER


def audit(mu, rho, first):
    import mono_floquet as mf
    o, _ = capture(mf.run, mu, rho, log_name='v3_weak_adjoint__floquet.out', append=not first)
    errs = {n: float(e[0]) for n, e in o['window_errs'].items()}
    lam = [o['lam_eta'], o['lam_r'], o['lam_window']]
    spread = max(lam) - min(lam)
    mc = o['margin_consistency']
    checks = {
        'twelve_event_order': all(v == 2 for v in o['counts'].values()) and bool(o['order_ok']) and order_matches(o['seq']) and bool(o['Psi_crest_pos'])
                              and bool(o['Psi_trough_neg']) and bool(o['ties_ok']),
        'supporting_plane_chords': o['chord_viol'] == 0 and o['neardiag_pos_resolved'] == 0 and o['det_crosscheck_mismatch'] == 0,
        'supporting_plane_E': o['E_max_ratio'] < 0,
        'supporting_plane_curvature': o['curv_max_ratio'] < 0,
        'supporting_line': bool(o['t1_in_F']) and o['U_t1'] < 0 and o['ell_neg_count'] == 0 and o['L_at_E'] > 0 and o['ell_dd_t1'] > 0,
        'rank_one_identity': max(errs[n] for n in COUNTED) <= 1e-4 and bool(o['Psi2_t1_neg']),
        'special_signs': all(tuple(o['window_signs'][n]) == (1.0, 1.0) for n in ('1', "p'", '-U')),
        'exponent_three_ways': spread <= 1e-7 * (1 + abs(o['lam_eta'])),
        'strong_stable_ratios': o['rs_defect'] <= 1e-7 and math.isfinite(o['min_ln_s']) and math.isfinite(o['margin_min_log10'])
                                and (not math.isfinite(mc) or mc <= 1e-2),
    }
    return dict(id=f'mu_{mu:g}_rho_{rho:g}', mu=mu, rho=rho, T=o['T'], lambda_c_from_eta=o['lam_eta'], lambda_c_from_r=o['lam_r'],
                lambda_c_from_window=o['lam_window'], lambda_c_from_return_map_fd=o['lam_fd'], exponent_spread=spread,
                sign_changes=o['counts'], event_sequence=o['seq'], events=o['events'], ties=o['ties'], min_event_gap=o['min_event_gap'],
                chord_pairs=o['chord_pairs'], chord_violations=o['chord_viol'], chord_max_ratio=o['chord_max_ratio'],
                chord_unresolved=o['chord_unres'], chord_small_pairs=o['chord_small_n'], chord_small_max_dt_over_T=o['chord_small_maxdt'],
                plain_determinant_crosschecks=o['det_crosscheck_n'], E_max_ratio=o['E_max_ratio'], curvature_max_ratio=o['curv_max_ratio'],
                t1_over_T=o['t1_over_T'], U_at_t1=o['U_t1'], supporting_line_min=o['min_ell_over_scale'],
                window_relative_errors=errs, window_error_counted=max(errs[n] for n in COUNTED), window_error_phi_y=errs['y'],
                window_signs={n: list(v) for n, v in o['window_signs'].items()},
                min_r_minus_max_1_mu_double=o['min_r_minus_max'], margin_min_log10=o['margin_min_log10'],
                margin_consistency=mc, min_log_s=o['min_ln_s'], ratio_defect=o['rs_defect'], adjoint_ratio_defect=o['ez_defect'],
                perron_slope_fd=o['perron_slope'], zeta_crest=o['zeta_crest'], return_map_positive_fd=o['Mpos'],
                attracting=bool(o['lam_eta'] < 0), checks=checks, attempted=len(checks), n_passed=sum(map(bool, checks.values())),
                passed=all(checks.values()), seconds=o['secs'])


def main():
    warnings.filterwarnings('ignore'); np.seterr(all='ignore')
    clock = Clock()
    recs = []
    for i, (mu, rho) in enumerate(GRID):
        try:
            r = audit(mu, rho, i == 0)
        except Exception as exc:
            r = dict(id=f'mu_{mu:g}_rho_{rho:g}', mu=mu, rho=rho, attempted=1, n_passed=0, passed=False, error=repr(exc))
        recs.append(r)
        print(r['id'], 'PASS' if r['passed'] else 'FAIL', {k_: v for k_, v in r.get('checks', {}).items() if not v},
              r.get('error', ''), f"{r.get('seconds', 0):.1f}s", flush=True)
    good = [r for r in recs if 'checks' in r]
    worst = {
        'chord_pairs_total': int(sum(r['chord_pairs'] for r in good)),
        'chord_violations_total': int(sum(r['chord_violations'] for r in good)),
        'max_chord_ratio': max(r['chord_max_ratio'] for r in good),
        'max_E_ratio': max(r['E_max_ratio'] for r in good),
        'max_curvature_ratio': max(r['curvature_max_ratio'] for r in good),
        'max_window_error_counted_test_functions': max(r['window_error_counted'] for r in good),
        'max_exponent_spread': max(r['exponent_spread'] for r in good),
        'cycles_with_tied_events': sum(bool(r['ties']) for r in good),
        'cycles_with_phi_y_error_below_1e-5': sum(r['window_error_phi_y'] <= 1e-5 for r in good),
        'max_window_error_phi_y_not_counted': max(r['window_error_phi_y'] for r in good),
        'min_margin_log10': min(r['margin_min_log10'] for r in good),
        'max_lambda_c': max(r['lambda_c_from_eta'] for r in good),
        'min_supporting_line_ratio': min(r['supporting_line_min'] for r in good),
    }
    extra = {'cycles': len(GRID), 'number_attracting': sum(r.get('attracting', False) for r in recs),
             'counting_rule': 'one count per (cycle, predicate)', 'test_functions_counted': COUNTED,
             'test_function_reported_not_counted': 'y'}
    ok = write_receipt('v3_weak_adjoint', 'float', recs, worst=worst, clock=clock, extra=extra,
                       attempted=sum(r['attempted'] for r in recs), n_passed=sum(r['n_passed'] for r in recs),
                       passed=all(r['passed'] for r in recs), parameters={'cycles_mu_rho': GRID},
                       scope='Sampled sign tests and quadratures on computed cycles; not an enclosure and not a statement about all cycles.')
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
