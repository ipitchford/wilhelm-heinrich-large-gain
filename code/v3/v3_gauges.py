#!/usr/bin/env python3
"""Floating-point check of the gauge theorem and its corollaries on computed cycles.

Cycles: the 31 saved amplitude-fixed cycles of the framework suite (re-polished by Newton iteration)
and an extended set  mu in {0.02, 0.1, 0.5, 1, 3, 10, 50} x rho in {1.001, 1.1, 2, 6, 20}.

Per cycle (predicates fixed in advance; TOL is a rounding allowance, not a fitted threshold):
  star_gauge_positive        min Sigma/scale > -TOL_SIGN  and  Sigma > 0 wherever |Sigma|/scale > TOL_SIGN
                             (scale = sum of the absolute values of the monomials of the gauge; signed-log arithmetic)
  convexity_gauge_positive   the same for Xi
  K_positive, W_positive     the same for K and W
  exponents_agree            |<eta>-mu + <pK/Sigma>| and |<eta>-mu + <pW/Xi>|  <=  TOL_AVG * (<|eta-mu|>+<|pK/Sigma|>+<|pW/Xi|>)
  eta0_mean                  |<eta_0> - mu| <= TOL_AVG * (mu + <|eta_0 - mu|>)
  polar_winding              |int_0^T Sigma/(U^2+p^2) dt - 2 pi| <= TOL_AVG
  monotone_integrals         the symmetric matrix [[<eta p>, <zeta p>],[<zeta p>, <zeta p'>]] is negative definite,
                             and |<zeta p> - <eta p'>| is at rounding level
  holder_sandwich            <x^{1/3}>^3 < gamma(gamma-1)(gamma-mu) < k,  gamma = 1 + <eta>
  attracting                 <eta> - mu < 0   (reported; it is the conjecture's content at this cycle, not a theorem)
  cycle_located              the crest return of the located point reproduces it to 1e-9 relative to the crest state
                             (|dU| <= 1e-9 (1+|U_max|), |dy| <= 1e-9 y_crest), and one period of the profile system closes
Nothing here is an interval enclosure, and sampled positivity is not positivity on a continuum.

Checker history (kept, not hidden).  The first version of `cycle_located` scaled the y-residual of the crest
return by k.  At mu = 0.02, rho = 20 the crest value of y is 775 k, the return is computed to relative accuracy
2.5e-12, and the k-scaled residual 1.8e-9 exceeded the 1e-9 gate; every mathematical predicate passed.  That
receipt is retained as evidence/v3/v3_gauges_initial_failure.json.  The residual is now measured relative to the
crest state.  In the same revision W was re-expressed, for mu > 1, in variables that do not cancel in the trough
(see cyclekit); before that the sign of W was unresolved at 11236 of 137303 samples, all on cycles with mu >= 3.
"""
from __future__ import annotations
import sys, time
import numpy as np
from common import ROOT, Clock, write_receipt, read_json
import cyclekit as ck

TOL_SIGN = 1e-9
TOL_AVG = 1e-7
MUS = [0.02, 0.1, 0.5, 1.0, 3.0, 10.0, 50.0]
RHOS = [1.001, 1.1, 2.0, 6.0, 20.0]


def saved_cycles():
    rows = read_json(ROOT / 'evidence/cycle_trials.json')['attempts']
    out = []
    for r in rows:
        if not r['accepted']:
            continue
        with np.load(ROOT / r['mesh'], allow_pickle=False) as q:
            mu = float(q['mu']); lam = float(q['lam']); st = q['states'][:, 0]
        k = mu * lam
        out.append(dict(id='saved_' + r['id'], mu=mu, k=k, U0=float(st[0]), y0=float(k * st[2]), source='saved v2.0 mesh'))
    return out


def sign_test(pair):
    ratio = pair[0]
    resolved = np.abs(ratio) > TOL_SIGN
    return dict(min_ratio=float(ratio.min()), unresolved=int((~resolved).sum()),
                negative_resolved=int(((ratio < 0) & resolved).sum()),
                passed=bool(ratio.min() > -TOL_SIGN and not ((ratio < 0) & resolved).any()))


def state_functions(mu, k, g):
    """signs of the explicit state functions eta_0 = V3/Sigma and eta_1 = (x Sigma + k psi(U) p')/Xi - 1 on the samples
    (reported, not a pass/fail predicate: the paper states that they need not be positive)"""
    U, q, y = g['U'], g['q'], g['y']
    (rS, mS), (rX, mX), (r3, m3) = g['Sg'], g['Xg'], g['V3']
    eta0_neg = bool(np.any(r3 < -TOL_SIGN))
    lx = np.log(k) + U
    with np.errstate(all='ignore'):
        kpsi = np.where(np.abs(U) < 1e-2, k * (U**2 / 2 + U**3 / 3 + U**4 / 8), k * (np.exp(U) * (U - 1.0) + 1.0))
        pd = mu * q - y
        L = np.stack([lx + mS + np.log(rS), np.log(kpsi) + np.log(np.abs(pd)), mX + np.log(rX)])
        S = np.stack([np.ones_like(U), np.sign(pd), -np.ones_like(U)])
        L = np.where(np.isfinite(L), L, -np.inf)
        M = np.max(L, axis=0)
        w = np.exp(L - M)
        ratio1 = np.sum(S * w, axis=0) / np.sum(w, axis=0)
        eta0 = (r3 / rS) * np.exp(np.clip(m3 - mS, -700, 700))
    return dict(eta0_negative_somewhere=eta0_neg, eta0_min=float(np.min(eta0)),
                eta0_fraction_of_samples_negative=float(np.mean(r3 < -TOL_SIGN)),
                eta1_negative_somewhere=bool(np.any(ratio1 < -1e-6)), eta1_min_normalised=float(np.nanmin(ratio1)),
                eta1_fraction_of_samples_negative=float(np.mean(ratio1 < -1e-6)))


def audit(case):
    t0 = time.time()
    mu, k = case['mu'], case['k']
    if 'U0' in case:
        cyc = ck.from_crest(mu, k, case['U0'], case['y0'])
    else:
        cyc = ck.find(mu, case['rho'])
    pr = ck.profile(mu, k, cyc)
    g, I, T = pr['g'], pr['means'], pr['T']
    lam_eta, mK, mW, e0 = I[0], I[1], I[2], I[3]
    scale = I[4] + I[5] + I[6]
    gam = 1.0 + lam_eta + mu
    chi = gam * (gam - 1.0) * (gam - mu)
    x13 = I[13] ** 3
    Mq = np.array([[I[9], I[11]], [I[11], I[10]]])
    ev = np.linalg.eigvalsh(Mq)
    s = {n: sign_test(g[n]) for n in ('Sg', 'Xg', 'K', 'W')}
    odef = float(np.max(np.abs(g['O'][0])))
    checks = {
        'cycle_located': cyc['res_rel'] < 1e-9 and pr['closure'] < 1e-7 * (1 + abs(cyc['U0'])),
        'star_gauge_positive': s['Sg']['passed'],
        'convexity_gauge_positive': s['Xg']['passed'],
        'K_positive': s['K']['passed'],
        'W_positive': s['W']['passed'],
        'orthogonality': odef < 1e-7,
        'exponents_agree': abs(lam_eta + mK) <= TOL_AVG * scale and abs(lam_eta + mW) <= TOL_AVG * scale,
        'eta0_mean': abs(e0) <= TOL_AVG * (mu + I[7]),
        'polar_winding': abs(I[8] * T - 2 * np.pi) <= TOL_AVG,
        'monotone_integrals': bool(ev.max() < 0) and abs(I[11] - I[12]) <= TOL_AVG * I[15],
        'holder_sandwich': bool(x13 < chi < k),
        'mean_x': abs(I[14]) <= TOL_AVG * k,
    }
    rec = dict(id=case['id'], source=case['source'], mu=mu, k=k, rho=k / ck.kH(mu), T=T, U_max=cyc['U0'],
               U_min=float(g['U'].min()), newton_residual=cyc['res'], location_residual_relative=cyc['res_rel'],
               profile_closure=pr['closure'], W_formula=pr['W_formula'], alt_variable_consistency=pr.get('alt_consistency'),
               samples=int(pr['t'].size),
               lambda_c_from_eta=lam_eta, lambda_c_from_K=-mK, lambda_c_from_W=-mW,
               exponent_gap_K=abs(lam_eta + mK), exponent_gap_W=abs(lam_eta + mW), exponent_scale=scale,
               exponent_gap_relative=max(abs(lam_eta + mK), abs(lam_eta + mW)) / scale,
               eta0_mean_minus_mu=e0, eta0_relative=abs(e0) / (mu + I[7]),
               polar_angle_error=abs(I[8] * T - 2 * np.pi),
               monotone_form_eigenvalues=ev.tolist(), cross_term_gap=abs(I[11] - I[12]),
               holder_lower=x13, holder_middle=chi, holder_upper=k,
               sign_tests=s, max_orthogonality_defect=odef, located_by=cyc.get('located_by', 'saved crest'),
               state_functions=state_functions(mu, k, g),
               attracting=bool(lam_eta < 0), checks=checks, passed=bool(all(checks.values())), seconds=time.time() - t0)
    return rec


def main():
    clock = Clock()
    if not ck.w_identity_check():
        raise RuntimeError('the cancellation-free form of W is not identical to its definition')
    cases = saved_cycles()
    for mu in MUS:
        for rho in RHOS:
            cases.append(dict(id=f'ext_mu_{mu:g}_rho_{rho:g}', mu=mu, rho=rho, k=rho * ck.kH(mu), source='extended set'))
    recs = []
    for c in cases:
        try:
            r = audit(c)
        except Exception as exc:       # a failure to compute is recorded as a failure, never dropped
            r = dict(id=c['id'], source=c['source'], mu=c['mu'], k=c['k'], passed=False, error=repr(exc))
        recs.append(r)
        print(c['id'], 'PASS' if r['passed'] else 'FAIL', f"lambda_c={r.get('lambda_c_from_eta')}",
              f"gap={r.get('exponent_gap_relative')}", f"{r.get('seconds', 0):.1f}s", flush=True)
        if not r['passed']:
            print('   ', {k_: v for k_, v in r.get('checks', {}).items() if not v}, r.get('error', ''), flush=True)
    good = [r for r in recs if 'checks' in r]
    worst = {
        'max_exponent_gap_relative': max(r['exponent_gap_relative'] for r in good),
        'max_exponent_gap_absolute': max(max(r['exponent_gap_K'], r['exponent_gap_W']) for r in good),
        'max_eta0_mean_relative': max(r['eta0_relative'] for r in good),
        'max_polar_angle_error': max(r['polar_angle_error'] for r in good),
        'min_star_gauge_ratio': min(r['sign_tests']['Sg']['min_ratio'] for r in good),
        'min_convexity_gauge_ratio': min(r['sign_tests']['Xg']['min_ratio'] for r in good),
        'min_K_ratio': min(r['sign_tests']['K']['min_ratio'] for r in good),
        'min_W_ratio': min(r['sign_tests']['W']['min_ratio'] for r in good),
        'max_monotone_form_eigenvalue': max(max(r['monotone_form_eigenvalues']) for r in good),
        'max_lambda_c': max(r['lambda_c_from_eta'] for r in good),
        'unresolved_sign_samples': int(sum(sum(v['unresolved'] for v in r['sign_tests'].values()) for r in good)),
        'unresolved_sign_samples_by_gauge': {n: int(sum(r['sign_tests'][n]['unresolved'] for r in good)) for n in ('Sg', 'Xg', 'K', 'W')},
        'total_samples': int(sum(r['samples'] for r in good)),
        'max_location_residual_relative': max(r['location_residual_relative'] for r in good),
        'max_orthogonality_defect': max(r['max_orthogonality_defect'] for r in good),
        'max_holder_lower_over_middle': max(r['holder_lower'] / r['holder_middle'] for r in good),
        'max_holder_middle_over_upper': max(r['holder_middle'] / r['holder_upper'] for r in good),
    }
    extra = {'tolerances': {'TOL_SIGN': TOL_SIGN, 'TOL_AVG': TOL_AVG},
             'checker_history': 'see module docstring; the receipt of the first version is v3_gauges_initial_failure.json',
             'ranges': {'mu': [min(r['mu'] for r in good), max(r['mu'] for r in good)],
                        'rho': [min(r['rho'] for r in good), max(r['rho'] for r in good)],
                        'U_min': min(r['U_min'] for r in good)},
             'number_attracting': sum(r.get('attracting', False) for r in recs),
             'cycles_with_eta0_negative_somewhere': sum(r['state_functions']['eta0_negative_somewhere'] for r in good),
             'cycles_with_eta1_negative_somewhere': sum(r['state_functions']['eta1_negative_somewhere'] for r in good),
             'ids_with_eta0_negative_somewhere': [r['id'] for r in good if r['state_functions']['eta0_negative_somewhere']],
             'number_saved_cycles': sum(r['source'] == 'saved v2.0 mesh' for r in recs),
             'number_extended_cycles': sum(r['source'] == 'extended set' for r in recs)}
    ok = write_receipt('v3_gauges', 'float', recs, worst=worst, clock=clock, extra=extra,
                       parameters={'saved_cycles': 31, 'extended_mu': MUS, 'extended_rho': RHOS},
                       scope='Sampled diagnostics on computed approximations of cycles; not an enclosure and not a statement about all cycles.')
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
