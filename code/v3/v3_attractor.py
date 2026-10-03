#!/usr/bin/env python3
"""Floating-point tests of the attractor theorem and of the crest-axis return formula.

For each parameter pair the referee battery lib/attractor_battery.py samples the attractor (orbits of the unstable
manifold of E at several phases, in logarithmic variables, plus the cycle) and evaluates, without underflow:

  supporting_planes    Delta_ab = zeta_a (U_b-U_a) - eta_a (p_b-p_a) - (y_b-y_a) < 0 for every ordered pair of sampled
                       attractor points.  Pass: no pair exceeds the rounding-amplification threshold of the battery
                       (10 (1e-11/delta + 1e-11), delta the relative separation of the pair), and the largest value among
                       pairs with relative separation above 1e-6 is negative.
  K_positive, W_positive   on every sample of the unstable-manifold spiral and of the cycle
  orthogonality        |O| / largest monomial <= 1e-8 on every sample
  equilibrium_cones    no sample in E + C_+ or E + C_-, and Lambda_E > 0 on every sample
  riccati_bounds       log eta, log zeta inside the a-priori bounds
  pullback             the periodic adjoint ratios are reached from the four extreme data (1e-8 or 1e8 in each slot)
  exponent_two_ways    <eta> - mu against the sum of the two leading exponents from a discrete QR computation that
                       never uses the Riccati equation: difference <= 1e-8 (1 + |lambda_c|)
  left_eigenvector     Phi(T)^T n = m_s n at the crest, relative error <= 1e-9
  crest_axis_value     h(u,0) - k = k (e^u - 1)/(1 + eta) on the crests of the spiral, relative error <= 1e-8
  crests_increasing    successive crests of the spiral increase to the cycle crest
  return_derivative    R_1'(u) = g(u) e^{iota(u)} / g(R_1 u) against the U-component of M (1, zeta(u)), M the derivative of
                       the three-dimensional crest return map from the variational equations, and M (1, zeta(u)) parallel
                       to (1, zeta(R_1 u)); rows are used where |prediction| > 1000 x the noise floor of the variational
                       computation; pass: at least 3 such rows, relative error <= 1e-6, parallel defect <= 1e-6
  return_derivative_fd the same derivative from finite differences between neighbouring orbits of the unstable manifold
                       (no variational equation, no zeta); pass: relative error <= 1e-3 on the resolvable rows
  gradient_fd          grad h = (zeta, -eta) by finite differences.  Three orbits of the unstable manifold at phases
                       theta, theta +- 1e-3 are sampled at common times; with the central difference d = z_+ - z_- across
                       orbits and the flow direction f, the plane through the three points gives a gradient g_fd, compared
                       with (zeta, -eta) on the middle orbit.  Samples are used where the largest logarithmic variation
                       across the stencil is below 0.01 and the sine of the angle between the projections of d and f
                       exceeds 0.05, so that the truncation error is at most about 0.01^2/0.05 = 2e-3; pass: at least
                       20 such samples and relative error <= 5e-3 on all of them.  Near the Hopf curve neighbouring orbits
                       differ almost by a time shift, no sample is admissible, and the test is recorded as not applicable.

For mu > 1 and deep troughs the pair (y, eta) cancels in Delta; lib/attractor_mu_gt1.py repeats the supporting-plane and
W tests in variables without cancellation at two further pairs.
Sampled inequalities on computed orbits; nothing is an enclosure.
"""
from __future__ import annotations
import math, sys, time, warnings
import numpy as np
from common import LIB, Clock, write_receipt, exec_script, capture

PAIRS = [(0.02, 1.5), (0.02, 4.0), (0.3, 1.01), (0.3, 1.5), (0.3, 4.0), (0.3, 12.0),
         (1.0, 1.1), (1.0, 1.5), (1.0, 4.0), (1.0, 12.0), (3.0, 1.5)]
MU_GT1 = [(3.0, 12.0), (10.0, 4.0)]
FD_VAR, FD_SIN, FD_TOL, FD_MIN = 0.01, 0.05, 5e-3, 20


def gradient_fd(mu, k, ctx, dth=1e-3, nt=2000):
    import vlog
    lam = k / mu
    te = np.linspace(ctx['t_skip'], ctx['tconv'], nt + 2)[1:-1]
    th0 = ctx['phases'][0]
    Sm, S0, Sp = [vlog.run_log(mu, k, vlog.spiral_start(mu, k, ctx['r0'], th0 + d), 0.0, ctx['tconv'], t_eval=te).y
                  for d in (-dth, 0.0, dth)]
    U0, w0, c0, a0, b0 = S0[:5]
    dU = Sp[0] - Sm[0]
    dp = -(np.exp(Sp[1]) - np.exp(Sm[1]))
    dy = np.exp(Sp[2]) - np.exp(Sm[2])
    eta, zeta = np.exp(a0), np.exp(b0)
    var = np.maximum.reduce([np.abs(dU) / (1 + np.abs(U0)), np.abs(Sp[1] - Sm[1]), np.abs(Sp[2] - Sm[2])])
    p = lam - np.exp(w0); y = np.exp(c0); x = k * np.exp(U0)
    pd = mu * np.exp(w0) - y; yd = x - y
    det = dU * pd - dp * p
    with np.errstate(divide='ignore', invalid='ignore'):
        g1 = (dy * pd - dp * yd) / det
        g2 = (dU * yd - dy * p) / det
        err = np.maximum(np.abs(g1 - zeta) / zeta, np.abs(g2 + eta) / eta)
        sin = np.abs(det) / (np.hypot(dU, dp) * np.hypot(p, pd))
    use = (var < FD_VAR) & (sin > FD_SIN) & np.isfinite(err)
    n = int(use.sum())
    if n < FD_MIN:
        return dict(applicable=False, admissible_samples=n, samples=int(nt))
    return dict(applicable=True, admissible_samples=n, samples=int(nt), median_relative_error=float(np.median(err[use])),
                max_relative_error=float(err[use].max()), U_range_used=[float(U0[use].min()), float(U0[use].max())])


def audit(mu, rho, first):
    import attractor_battery as ab
    t0 = time.time()
    (out, ctx), _ = capture(ab.battery, mu, rho, log_name='v3_attractor__battery.out', append=not first)
    k = out['k']
    res = [r for r in out['prop81'] if r['var'] == r['var'] and abs(r['pred']) > 1e3 * r['noise_floor']]
    fd = out.get('prop81_fd', [])
    g = gradient_fd(mu, k, ctx)
    lam_c = out['lamc_riccati']
    checks = {
        'supporting_planes': out['support_significant'] == 0 and out['support_worst_rel_far'] < 0,
        'K_positive': out['K_rel_min'] > 0,
        'W_positive': out['W_rel_min'] > 0,
        'orthogonality': out['O_rel_max'] <= 1e-8,
        'equilibrium_cones': out['n_in_Cplus'] + out['n_in_Cminus'] == 0 and out['LamE_rel_min'] > 0,
        'riccati_bounds': bool(out['bounds_ok']),
        'pullback': all(p_['final_gap'] < 1e-8 for p_ in out['pullback']),
        'exponent_two_ways': abs(lam_c - out['lamc_qr']) <= 1e-8 * (1 + abs(lam_c)),
        'left_eigenvector': out['MTn_relerr'] <= 1e-9,
        'crest_axis_value': out['gtilde_formula_relerr'] <= 1e-8,
        'crests_increasing': bool(out['crests_increasing']),
        'return_derivative': len(res) >= 3 and max(r['rel'] for r in res) <= 1e-6 and max(abs(r['parallel_defect']) for r in res) <= 1e-6,
    }
    not_applicable = []
    if fd:
        checks['return_derivative_fd'] = max(r['rel'] for r in fd) <= 1e-3
    else:
        not_applicable.append('return_derivative_fd')
    if g['applicable']:
        checks['gradient_fd'] = g['max_relative_error'] <= FD_TOL
    else:
        not_applicable.append('gradient_fd')
    rec = dict(id=f'mu_{mu:g}_rho_{rho:g}', mu=mu, rho=rho, k=k, T=out['T'], U_crest=out['U_crest'], U_min=out['Umin'],
               spiral_turns=out['spiral_turns'], attractor_points=out['N_points'], ordered_pairs=out['support_pairs'],
               support_nonnegative_pairs=out['support_nonneg'], support_beyond_rounding=out['support_significant'],
               support_worst_separated=out['support_worst_rel_far'], support_close_pairs=out['support_n_close'],
               support_histogram=out['support_hist'], K_min_ratio=out['K_rel_min'], W_min_ratio=out['W_rel_min'],
               orthogonality_max=out['O_rel_max'], LambdaE_min_ratio=out['LamE_rel_min'],
               lambda_c_riccati=lam_c, lambda_c_qr=out['lamc_qr'], lambda_c_from_return_map=out['lamc_M'],
               left_eigenvector_error=out['MTn_relerr'], crest_axis_value_error=out['gtilde_formula_relerr'],
               return_derivative_rows=len(res), return_derivative_max_error=max((r['rel'] for r in res), default=None),
               return_derivative_parallel_defect=max((abs(r['parallel_defect']) for r in res), default=None),
               return_derivative_fd_rows=len(fd), return_derivative_fd_max_error=max((r['rel'] for r in fd), default=None),
               gradient_fd=g, pullback_periods=[p_['periods'] for p_ in out['pullback']],
               pullback_max_gap=max(p_['final_gap'] for p_ in out['pullback']), attracting=bool(lam_c < 0),
               checks=checks, not_applicable=not_applicable, attempted=len(checks), n_passed=sum(map(bool, checks.values())),
               passed=all(checks.values()), seconds=time.time() - t0)
    return rec


def audit_mu_gt1(mu, rho, first):
    t0 = time.time()
    ns, text, st = exec_script(LIB / 'attractor_mu_gt1.py', argv=[repr(mu), repr(rho)])
    from common import LOGS
    LOGS.mkdir(parents=True, exist_ok=True)
    with (LOGS / 'v3_attractor__mu_gt1.out').open('w' if first else 'a') as f:
        f.write(text)
    r = ns['RESULT']
    checks = {'supporting_planes_cancellation_free': r['significant'] == 0 and r['far_worst'] < 0,
              'W_positive_cancellation_free': r['W_rel_min'] > 0}
    return dict(id=f'mu_{mu:g}_rho_{rho:g}_cancellation_free', mu=mu, rho=rho, attractor_points=r['N'], ordered_pairs=r['pairs'],
                U_min=r['U_min'], support_nonnegative_pairs=r['n_nonneg'], support_beyond_rounding=r['significant'],
                support_worst_separated=r['far_worst'], support_close_pairs=r['n_close'], support_histogram=r['hist'],
                W_min_ratio=r['W_rel_min'], checks=checks, not_applicable=[], attempted=len(checks),
                n_passed=sum(map(bool, checks.values())), passed=all(checks.values()), seconds=time.time() - t0)


def main():
    warnings.filterwarnings('ignore'); np.seterr(all='ignore')
    clock = Clock()
    recs = []
    for i, (mu, rho) in enumerate(PAIRS):
        try:
            r = audit(mu, rho, i == 0)
        except Exception as exc:
            r = dict(id=f'mu_{mu:g}_rho_{rho:g}', mu=mu, rho=rho, attempted=1, n_passed=0, passed=False, error=repr(exc))
        recs.append(r)
        print(r['id'], 'PASS' if r['passed'] else 'FAIL', {k_: v for k_, v in r.get('checks', {}).items() if not v},
              'n/a:', r.get('not_applicable'), r.get('error', ''), f"{r.get('seconds', 0):.1f}s", flush=True)
    for i, (mu, rho) in enumerate(MU_GT1):
        try:
            r = audit_mu_gt1(mu, rho, i == 0)
        except Exception as exc:
            r = dict(id=f'mu_{mu:g}_rho_{rho:g}_cancellation_free', mu=mu, rho=rho, attempted=1, n_passed=0, passed=False, error=repr(exc))
        recs.append(r)
        print(r['id'], 'PASS' if r['passed'] else 'FAIL', r.get('error', ''), f"{r.get('seconds', 0):.1f}s", flush=True)
    good = [r for r in recs if 'checks' in r]
    bat = [r for r in good if 'gradient_fd' in r]
    gfd = [r['gradient_fd'] for r in bat if r['gradient_fd']['applicable']]
    worst = {
        'ordered_pairs_total': int(sum(r['ordered_pairs'] for r in good)),
        'pairs_beyond_rounding_total': int(sum(r['support_beyond_rounding'] for r in good)),
        'pairs_nonnegative_total': int(sum(r['support_nonnegative_pairs'] for r in good)),
        'max_support_value_separated_pairs': max(r['support_worst_separated'] for r in good),
        'min_K_ratio': min(r['K_min_ratio'] for r in bat), 'min_W_ratio': min(r['W_min_ratio'] for r in good),
        'max_orthogonality': max(r['orthogonality_max'] for r in bat),
        'max_exponent_difference_riccati_vs_qr': max(abs(r['lambda_c_riccati'] - r['lambda_c_qr']) for r in bat),
        'max_left_eigenvector_error': max(r['left_eigenvector_error'] for r in bat),
        'max_crest_axis_value_error': max(r['crest_axis_value_error'] for r in bat),
        'max_return_derivative_error': max(r['return_derivative_max_error'] for r in bat if r['return_derivative_max_error'] is not None),
        'max_return_derivative_parallel_defect': max(r['return_derivative_parallel_defect'] for r in bat if r['return_derivative_parallel_defect'] is not None),
        'max_return_derivative_fd_error': max((r['return_derivative_fd_max_error'] for r in bat if r['return_derivative_fd_max_error'] is not None), default=None),
        'max_gradient_fd_error': max((g['max_relative_error'] for g in gfd), default=None),
        'max_gradient_fd_median_error': max((g['median_relative_error'] for g in gfd), default=None),
        'gradient_fd_pairs_applicable': len(gfd), 'gradient_fd_samples_total': int(sum(g['admissible_samples'] for g in gfd)),
        'max_lambda_c': max(r['lambda_c_riccati'] for r in bat), 'min_U': min(r['U_min'] for r in good),
    }
    extra = {'parameter_pairs': len(PAIRS), 'cancellation_free_pairs': len(MU_GT1),
             'not_applicable': {r['id']: r['not_applicable'] for r in good if r['not_applicable']},
             'counting_rule': 'one count per (parameter pair, predicate); predicates recorded as not applicable are not counted',
             'gradient_fd_filters': {'max_log_variation': FD_VAR, 'min_sine': FD_SIN, 'tolerance': FD_TOL, 'min_samples': FD_MIN}}
    ok = write_receipt('v3_attractor', 'float', recs, worst=worst, clock=clock, extra=extra,
                       attempted=sum(r['attempted'] for r in recs), n_passed=sum(r['n_passed'] for r in recs),
                       passed=all(r['passed'] for r in recs),
                       parameters={'battery_pairs_mu_rho': PAIRS, 'cancellation_free_pairs_mu_rho': MU_GT1},
                       scope='Sampled inequalities and finite-difference comparisons on computed orbits of the unstable manifold '
                             'and on computed cycles. Not an enclosure, and not a test of the topological part of the theorem.')
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
