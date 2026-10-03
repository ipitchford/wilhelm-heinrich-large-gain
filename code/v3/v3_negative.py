#!/usr/bin/env python3
"""Reproduction of the refutations cited in the section on negative results, at their stated parameters.

A record passes when the refutation is reproduced: the candidate property FAILS in the computation, as the paper states.

N1  crest return map not concave, not sub-homogeneous (mu = 0.0027, rho = 5.02).  On the unstable manifold of E the
    successive crests U_n give r(U_n) = log(U_{n+1}/U_n).  A concave return map phi with phi(0) = 0 has phi(U)/U
    non-increasing.  Reproduced when r has an interior minimum followed by a rise of at least 0.02 (stated: r = 0.2871 at E,
    0.2202 near U = 0.25 U*, 0.2452 near U = 0.54 U*).  Eight orbits, bundle integrator (logarithmic variables).
N2  Zhang-type monotonicity of (mu - zeta)/(k (e^U - 1)) along the section p = 0 of the attractor (mu = 1, rho = 5).
    The property requires the ratio to be increasing in U.  Reproduced when the fitted slope just inside the cycle,
    U in (U*(1 - 3e-6), U*(1 - 1e-8)), is below -0.1 (stated: -0.135) while the slope on 0.5 U* < U < 0.9 U* is positive.
    Twenty-four orbits of the unstable manifold.
N3  no certificate of the form M(t) (c1 x, c2 y, c3 (lambda - p)).  (a) the three quotient identities behind the averaging
    argument, exactly (sympy); (b) the Lotka-Volterra means <x/y> = 1, <y/(lambda-p)> = mu, <lambda-p> = lambda on 18 computed
    cycles (1e-8 relative); (c) the averaged inequalities c2 <= c3, c1 <= c2, c3/c1 <= 1/rho sum to 0 <= -log rho.
N4  pointwise candidates change sign on computed cycles (mu in {0.02, 0.1, 0.5, 1, 3, 10} x rho in {1.1, 2, 6}): each of
    x/zeta - A, zeta/eta + mu - A, k psi(U)/K - A, x p^2/W - A takes both signs on every cycle although its mean equals
    lambda_c < 0; zeta/eta takes values on both sides of 1; K - (zeta - mu) chi(U) takes both signs on the cycles with
    mu <= 1 and rho = 6.
N5  no amplitude selection at first order: on levels of p^2/2 + e^U - U - 1 the loop integrals of e^U p^2 and of (e^U - 1)^2
    both equal that of p^2, so the first-order Melnikov function is (1 - rho) times a positive integral.
Floating point (N3(a), N3(c) exact).  Nothing here is an enclosure.
"""
from __future__ import annotations
import math, sys, time, warnings
import numpy as np
from common import Clock, write_receipt

N4_MUS = [0.02, 0.1, 0.5, 1.0, 3.0, 10.0]
N4_RHOS = [1.1, 2.0, 6.0]


def n1(recs):
    import vlog
    t0 = time.time()
    mu, rho = 0.0027, 5.02
    k = rho * mu * (1 + mu)
    ed = vlog.eq_data(mu, k)
    r0 = 2 * math.pi * ed['lam'].real / ed['lam'].imag
    vals = []
    for th in np.linspace(0, 2 * np.pi, 9)[:-1]:
        _, cr = vlog.spiral_until_converged(mu, k, 1e-5, th, rtol=1e-12)
        Uc = cr[:, 1]
        vals.append(np.column_stack((Uc[2:-1], np.log(Uc[3:] / Uc[2:-1]))))
    V = np.vstack(vals)
    V = V[np.argsort(V[:, 0])]
    Us = V[-1, 0]
    i = int(np.argmin(np.where(V[:, 0] < 0.4 * Us, V[:, 1], np.inf)))
    j = int(np.argmax(np.where(V[:, 0] > V[i, 0], V[:, 1], -np.inf)))
    rise = float(V[j, 1] - V[i, 1])
    rec = dict(group='N1 return map', id='r(U) = log(phi(U)/U) is not monotone at mu=0.0027, rho=5.02', mu=mu, rho=rho,
               r_at_E=r0, U_star=float(Us), r_interior_minimum=float(V[i, 1]), minimum_at_U_over_Ustar=float(V[i, 0] / Us),
               r_later_maximum=float(V[j, 1]), maximum_at_U_over_Ustar=float(V[j, 0] / Us), rise=rise, crest_pairs=int(len(V)),
               stated={'r_at_E': 0.2871, 'r_interior_minimum': 0.2202, 'r_later_maximum': 0.2452},
               passed=bool(rise >= 0.02 and abs(r0 - 0.2871) <= 1e-4 and abs(V[i, 1] - 0.2202) <= 5e-3 and abs(V[j, 1] - 0.2452) <= 5e-3),
               seconds=time.time() - t0)
    recs.append(rec)
    return dict(return_ratio_rise=rise, return_ratio_at_E=r0)


def n2(recs):
    import vlog
    t0 = time.time()
    mu, rho = 1.0, 5.0
    k = rho * mu * (1 + mu)
    pts = []
    for th in np.linspace(0, 2 * np.pi, 25)[:-1]:
        _, cr = vlog.spiral_until_converged(mu, k, 1e-6, th)
        pts.append(np.column_stack((cr[4:, 1], np.exp(cr[4:, 4]))))
    P = np.vstack(pts)
    P = P[np.argsort(P[:, 0])]
    U, zeta = P[:, 0], P[:, 1]
    Us = U[-1]
    keep = np.concatenate(([True], np.diff(U) > 1e-9))
    U, zeta = U[keep], zeta[keep]
    R = (mu - zeta) / (k * np.expm1(U))
    near = (U > Us * (1 - 3e-6)) & (U < Us * (1 - 1e-8))
    mid = (U > 0.5 * Us) & (U < 0.9 * Us)
    s_near = float(np.polyfit(U[near], R[near], 1)[0]) if near.sum() >= 3 else float('nan')
    s_mid = float(np.polyfit(U[mid], R[mid], 1)[0]) if mid.sum() >= 3 else float('nan')
    recs.append(dict(group='N2 section ratio', id='(mu - zeta)/(k(e^U - 1)) decreases next to the cycle at mu=1, rho=5', mu=mu, rho=rho,
                     U_star=float(Us), slope_just_inside_cycle=s_near, points_just_inside=int(near.sum()), slope_mid_amplitude=s_mid,
                     points_mid=int(mid.sum()), stated={'slope_just_inside_cycle': -0.135},
                     passed=bool(near.sum() >= 3 and s_near <= -0.1 and abs(s_near + 0.135) <= 0.02 and s_mid > 0),
                     seconds=time.time() - t0))
    return dict(section_ratio_slope_near_cycle=s_near, section_ratio_slope_mid=s_mid)


def n3_exact(recs):
    import sympy as sp
    t = sp.Symbol('t')
    mu, k, c1, c2, c3 = sp.symbols('mu k c1 c2 c3', positive=True)
    lam, A = k / mu, 1 + mu
    U, p, y, M = [sp.Function(n)(t) for n in ('U', 'p', 'y', 'M')]
    x, q = k * sp.exp(U), lam - p
    flow = {sp.Derivative(U, t): p, sp.Derivative(p, t): k - mu * p - y, sp.Derivative(y, t): x - y}
    D = lambda e: sp.diff(e, t).subs(flow)
    m = sp.Derivative(M, t) / M
    G1, G2, G3 = M * c1 * x, M * c2 * y, M * c3 * q
    r1, r2 = x / y, y / q
    tests = [('quotient identity 3: (G3\' + mu G3 - G2)/G3 = M\'/M + (1 - c2/c3) y/q', (D(G3) + mu * G3 - G2) / G3 - (m + (1 - c2 / c3) * r2)),
             ('quotient identity 2: (G2\' + G2 - G1)/G2 = M\'/M + (1 - c1/c2) x/y', (D(G2) + G2 - G1) / G2 - (m + (1 - c1 / c2) * r1)),
             ('quotient identity 1: (G1\' + A G1 - x G3)/G1 = M\'/M + p + A - (c3/c1)(lambda - p)', (D(G1) + A * G1 - x * G3) / G1 - (m + p + A - (c3 / c1) * q))]
    for name, e in tests:
        r = sp.simplify(e)
        recs.append(dict(group='N3 monomial certificate', id=name, sub_kind='exact symbolic equality', residual=str(r), passed=bool(r == 0)))
    u1, u2, u3, lr = sp.symbols('u1 u2 u3 log_rho')
    lhs = (u2 - u3) + (u1 - u2) + (u3 - u1)
    recs.append(dict(group='N3 monomial certificate', id='averaged inequalities sum to 0 <= -log rho (Farkas multipliers 1, 1, 1)',
                     sub_kind='exact symbolic equality', residual=str(sp.expand(lhs)), passed=bool(sp.expand(lhs) == 0),
                     note='u_i = log c_i; the inequalities are u2 - u3 <= 0, u1 - u2 <= 0, u3 - u1 <= -log rho'))


def psi_k(k, U, x):
    """k psi(U) = x U - x + k without cancellation for small |U|"""
    if abs(U) < 1e-2:
        return k * (U**2 / 2 + U**3 / 3 + U**4 / 8 + U**5 / 30 + U**6 / 144)
    return x * U - x + k


def chi(U):
    """chi(U) = U - 2 (e^U - U - 1)/(e^U - 1) = U^2/6 - U^4/360 + ..."""
    U = np.asarray(U, float)
    small = np.abs(U) < 1e-2
    with np.errstate(all='ignore'):
        big = U - 2 * (np.expm1(U) - U) / np.expm1(U)
    return np.where(small, U**2 / 6.0 - U**4 / 360.0, big)


def n34_cycles(recs):
    import cyclekit as ck
    worst = dict(lv=0.0, mean=0.0)
    n_sign = 0
    for mu in N4_MUS:
        for rho in N4_RHOS:
            t0 = time.time()
            ident = f'cycle mu={mu:g} rho={rho:g}'
            k = rho * ck.kH(mu); lam = k / mu; A = 1 + mu
            try:
                cyc = ck.find(mu, rho)
                ext = lambda d: [d['x'] / d['zeta'] - A, d['zeta'] / d['eta'] + mu - A, psi_k(k, d['U'], d['x']) / d['K'] - A,
                                 d['x'] * d['p']**2 / d['W'] - A, d['x'] / d['y'] - 1.0, d['y'] / d['q'] - mu, d['q'] - lam,
                                 abs(d['x'] / d['zeta'] - A) + abs(d['zeta'] / d['eta'] + mu - A) + abs(psi_k(k, d['U'], d['x']) / d['K'] - A)
                                 + abs(d['x'] * d['p']**2 / d['W'] - A), abs(d['x'] / d['y'] - 1.0), abs(d['y'] / d['q'] - mu)]
                pr = ck.profile(mu, k, cyc, extras=ext, n_extras=10)
            except Exception as exc:
                recs.append(dict(group='N4 pointwise candidates', id=ident, passed=False, error=repr(exc)))
                continue
            g, em = pr['g'], pr['extra_means']
            lam_c = float(pr['means'][0])
            U, p, eta, zeta = g['U'], g['p'], g['eta'], g['zeta']
            x = k * np.exp(U)
            K = g['K'][0] * np.exp(g['K'][1]); W = g['W'][0] * np.exp(g['W'][1])
            kps = np.array([psi_k(k, u, xx) for u, xx in zip(U, x)])
            with np.errstate(all='ignore'):
                R = {'x/zeta - A': x / zeta - A, 'zeta/eta + mu - A': zeta / eta + mu - A, 'k psi(U)/K - A': kps / K - A,
                     'x p^2/W - A': x * p**2 / W - A}
                Ehat = (K - (zeta - mu) * chi(U)) / k
            scale = float(em[7])
            # N3(b): Lotka-Volterra means
            lv = max(abs(em[4]) / max(em[8], 1e-300), abs(em[5]) / max(em[9], 1e-300), abs(em[6]) / lam)
            worst['lv'] = max(worst['lv'], float(lv))
            recs.append(dict(group='N3 monomial certificate', id=ident + ': Lotka-Volterra means', sub_kind='floating point',
                             mean_x_over_y_minus_1=float(em[4]), mean_y_over_q_minus_mu=float(em[5]), mean_q_minus_lambda=float(em[6]),
                             relative_defect=float(lv), passed=bool(lv <= 1e-8)))
            # N4: sign changes and means
            for j, (name, v) in enumerate(R.items()):
                v = v[np.isfinite(v)]
                gap = abs(float(em[j]) - lam_c) / max(scale, 1e-300)
                worst['mean'] = max(worst['mean'], gap)
                both = bool(v.max() > 0 and v.min() < 0)
                n_sign += both
                recs.append(dict(group='N4 pointwise candidates', id=f'{ident}: {name} changes sign; mean = lambda_c', minimum=float(v.min()),
                                 maximum=float(v.max()), mean=float(em[j]), lambda_c=lam_c, mean_gap_relative=float(gap),
                                 passed=bool(both and gap <= 1e-7 and lam_c < 0)))
            ze = zeta / eta
            recs.append(dict(group='N4 pointwise candidates', id=f'{ident}: zeta/eta on both sides of 1', minimum=float(ze.min()),
                             maximum=float(ze.max()), passed=bool(ze.min() < 1 < ze.max())))
            if mu <= 1 and rho == 6.0:
                Eh = Ehat[np.isfinite(Ehat)]
                recs.append(dict(group='N4 pointwise candidates', id=f'{ident}: K - (zeta - mu) chi(U) takes both signs',
                                 minimum_over_k=float(Eh.min()), maximum_over_k=float(Eh.max()),
                                 fraction_of_samples_negative=float(np.mean(Eh < 0)), passed=bool(Eh.min() < 0 < Eh.max())))
            print(ident, f'lambda_c={lam_c:.6g}', f'{time.time()-t0:.1f}s', flush=True)
    return dict(max_lotka_volterra_mean_defect=worst['lv'], max_candidate_mean_gap=worst['mean'], candidate_sign_changes=n_sign)


def n5(recs):
    from toda import level_data
    worst = 0.0
    for h in (0.01, 0.1, 1.0, 10.0, 100.0):
        d = level_data(h)
        a, b = abs(d['IE'] / d['I'] - 1), abs(d['G2'] / d['I'] - 1)
        worst = max(worst, a, b)
        recs.append(dict(group='N5 first-order Melnikov function', id=f'level h={h:g}: loop integrals of e^U p^2 and (e^U-1)^2 equal that of p^2',
                         defect_exp_weight=float(a), defect_force_square=float(b), loop_integral_p2=float(d['I']), passed=bool(max(a, b) <= 1e-8)))
    return dict(max_melnikov_identity_defect=worst)


def main():
    warnings.filterwarnings('ignore'); np.seterr(all='ignore')
    clock = Clock()
    recs = []
    worst = {}
    for part in (n1, n2, n3_exact, n34_cycles, n5):
        try:
            worst.update(part(recs) or {})
        except Exception as exc:
            recs.append(dict(group=part.__name__, id='exception', passed=False, error=repr(exc)))
        print(part.__name__, 'done', f'{clock.seconds():.0f}s', flush=True)
    groups = {}
    for r in recs:
        g = groups.setdefault(r['group'], {'attempted': 0, 'passed': 0})
        g['attempted'] += 1; g['passed'] += int(r['passed'])
    extra = {'groups': groups, 'meaning_of_passed': 'the refutation stated in the paper is reproduced (the candidate property fails)',
             'counting_rule': 'one count per refutation, identity, or (cycle, candidate) pair'}
    ok = write_receipt('v3_negative', 'float', recs, worst=worst, clock=clock, extra=extra,
                       parameters={'N1': {'mu': 0.0027, 'rho': 5.02}, 'N2': {'mu': 1.0, 'rho': 5.0}, 'N4_mu': N4_MUS, 'N4_rho': N4_RHOS},
                       scope='Counterexamples to candidate certificates on computed orbits; the sign changes are far above rounding, '
                             'but nothing is an enclosure.')
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
