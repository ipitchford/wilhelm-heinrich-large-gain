#!/usr/bin/env python3
"""Deliberate-error controls.  A control passes only when the deliberately wrong assertion is rejected by the same
predicate, script and exit status that the real checks use.

C1  exact algebra        v3_exact_identities.py --inject-error adds the star-gauge identity with the wrong sign of p K.
                         Required: non-zero exit status, receipt not passed, exactly the injected test fails.
C2  star gauge off a cycle   Sigma = p^2 - U p' is evaluated, with the code used for cycles, on a non-periodic orbit started
                         at (U, p, y) = (0.5, 0, k/2) for mu = 1, rho = 2, where Sigma(0) = U (y - k) < 0.  Required: the sign
                         predicate of v3_gauges rejects "Sigma > 0", although it accepts it on the cycle at the same parameters.
C3  large-gain constant  v3_largegain_constants.py with c_1 = 0.38 instead of 0.374, and with c_p = 0.83 instead of 0.819.
                         Required: non-zero exit status and a failed comparison in each run.
C4  exponent formula     on the cycle at mu = 1, rho = 2 the predicate of v3_gauges for lambda_c = -<p K/Sigma> is applied to the
                         wrong-sign formula lambda_c = +<p K/Sigma>.  Required: rejected, while the correct formula is accepted.
C5  event order          the twelve-event predicate of v3_weak_adjoint is applied to the reversed event sequence of a computed
                         cycle.  Required: rejected, while the computed sequence is accepted.
The auxiliary receipts written by C1 and C3 are failing receipts by construction and are named *_control_*.
"""
from __future__ import annotations
import subprocess, sys, warnings
import numpy as np
from common import ROOT, V3, EVID, LOGS, Clock, write_receipt, read_json


def run_script(args, log):
    r = subprocess.run([sys.executable] + args, cwd=V3, text=True, capture_output=True)
    LOGS.mkdir(parents=True, exist_ok=True)
    (LOGS / log).write_text(r.stdout + r.stderr)
    return r.returncode


def main():
    warnings.filterwarnings('ignore'); np.seterr(all='ignore')
    clock = Clock()
    recs = []
    # ---- C1
    out = EVID / 'v3_control_exact_wrong_sign.json'
    status = run_script([str(V3 / 'v3_exact_identities.py'), '--inject-error', '--out', str(out)], 'v3_controls__exact_wrong_sign.log')
    rc = read_json(out)
    failed = [r['id'] for r in rc['records'] if not r['passed']]
    recs.append(dict(id='C1 wrong sign in the star-gauge identity is rejected', observed_exit_status=status, receipt_passed=rc['passed'],
                     failed_tests=failed, tests_in_receipt=rc['counts']['attempted'],
                     passed=bool(status != 0 and not rc['passed'] and failed == ['injected_wrong_sign_star_gauge'])))
    # ---- C2 and C4 share one cycle
    import cyclekit as ck
    import vlog
    from v3_gauges import sign_test, TOL_AVG
    from scipy.integrate import solve_ivp
    mu, rho = 1.0, 2.0
    k = rho * ck.kH(mu)
    cyc = ck.find(mu, rho)
    pr = ck.profile(mu, k, cyc)
    on_cycle = sign_test(pr['g']['Sg'])
    s0 = vlog.to_log(mu, k, 0.5, 0.0, 0.5 * k, 1.0, 1.0)
    sol = solve_ivp(vlog.log_rhs(mu, k), [0.0, 2 * cyc['T']], s0, method='DOP853', rtol=1e-12, atol=1e-13, max_step=cyc['T'] / 500)
    g = ck.gauges(mu, k, *sol.y[:5])
    off_cycle = sign_test(g['Sg'])
    recs.append(dict(id='C2 "Sigma > 0" is rejected on a non-periodic orbit where it is false', mu=mu, rho=rho,
                     initial_state={'U': 0.5, 'p': 0.0, 'y_over_k': 0.5}, samples=int(sol.t.size),
                     predicate_on_cycle=on_cycle, predicate_on_nonperiodic_orbit=off_cycle,
                     passed=bool(on_cycle['passed'] and not off_cycle['passed'] and off_cycle['negative_resolved'] > 0)))
    # ---- C3
    for key, flag, val, tag in (('c_1', '--override-c1', '0.38', 'c1'), ('c_p', '--override-cp', '0.83', 'cp')):
        out = EVID / f'v3_control_largegain_wrong_{tag}.json'
        status = run_script([str(V3 / 'v3_largegain_constants.py'), flag, val, '--out', str(out)], f'v3_controls__largegain_wrong_{tag}.log')
        rc = read_json(out)
        failed = [r['id'] for r in rc['records'] if not r['passed']]
        recs.append(dict(id=f'C3 wrong large-gain constant {key} = {val} is rejected', observed_exit_status=status,
                         receipt_passed=rc['passed'], failed_comparisons=failed,
                         passed=bool(status != 0 and not rc['passed'] and len(failed) >= 1)))
    # ---- C4
    I = pr['means']
    lam_eta, mK = float(I[0]), float(I[1])
    scale = float(I[4] + I[5] + I[6])
    agree = lambda a, b: abs(a - b) <= TOL_AVG * scale
    recs.append(dict(id='C4 exponent formula with the wrong sign is rejected', mu=mu, rho=rho, lambda_c_from_eta=lam_eta,
                     minus_mean_pK_over_Sigma=-mK, plus_mean_pK_over_Sigma=mK, correct_formula_accepted=bool(agree(lam_eta, -mK)),
                     wrong_formula_accepted=bool(agree(lam_eta, mK)), passed=bool(agree(lam_eta, -mK) and not agree(lam_eta, mK))))
    # ---- C5
    from v3_weak_adjoint import order_matches
    import mono_floquet as mf
    from common import capture
    o, _ = capture(mf.run, 1.0, 1.5, quick=True)
    seq = list(o['seq'])
    recs.append(dict(id='C5 reversed twelve-event order is rejected', mu=1.0, rho=1.5, computed_sequence=seq,
                     computed_order_accepted=bool(order_matches(seq)), reversed_order_accepted=bool(order_matches(seq[::-1])),
                     passed=bool(order_matches(seq) and not order_matches(seq[::-1]))))
    ok = write_receipt('v3_controls', 'control', recs, worst={}, clock=clock,
                       extra={'interpretation': 'A control passes only when the deliberately wrong assertion is rejected.',
                              'failing_receipts_by_construction': ['v3_control_exact_wrong_sign.json', 'v3_control_largegain_wrong_c1.json',
                                                                   'v3_control_largegain_wrong_cp.json']},
                       parameters={'cycle_for_C2_C4': {'mu': 1.0, 'rho': 2.0}, 'cycle_for_C5': {'mu': 1.0, 'rho': 1.5}},
                       scope='Regression tests of the failure semantics of the suite.')
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
