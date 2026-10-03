#!/usr/bin/env python3
"""Cycle-level tests on all grid cycles: convexity of the projected cycle (S4b), sign pattern of (1+eta)E (S5),
pointwise ratios (S6), plus brief identities (c),(f),(g).  Writes cycle_tests.csv."""
import sys, csv
sys.path.insert(0, '/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03/work/lanes/cartography')
import numpy as np
import cart as C
from stests import chi, g2

rows = [{k: (float(v) if v not in ('', 'nan') else np.nan) for k, v in r.items()} for r in csv.DictReader(open('grid_table.csv'))]
out = []
for r in rows:
    mu, k, T = r['mu'], r['k'], r['T']
    A = 1 + mu
    U0, y0 = r['U_max'], r['y0_over_k'] * k
    pr = C.profile(mu, k, U0, y0, T, nmin=6000)
    t, U, p, y, eta, zeta, K, W = (pr[n] for n in ('t', 'U', 'p', 'y', 'eta', 'zeta', 'K', 'W'))
    x = k * np.exp(U)
    pd = k - mu * p - y
    pdd = -mu * pd - (x - y)
    cross = p * pdd - pd * pd
    # scale-free curvature sign; ignore points where both are at round-off
    sc = np.abs(p * pdd) + pd * pd
    sg = np.sign(np.where(np.abs(cross) > 1e-9 * sc.max() * 1e-6, cross, 0.0))
    nz = sg[sg != 0]
    n_infl = int(np.sum(nz[1:] != nz[:-1]))
    frac_pos = float(np.trapezoid((cross > 0).astype(float), t) / T)
    # where is the curve non-convex (cross > 0 for a clockwise curve)
    bad = cross > 0
    lam = k / mu
    with np.errstate(divide='ignore', invalid='ignore'):
        ze = zeta / eta
        R1 = x / zeta
        R2 = ze + mu
        R3 = k * g2(U) / K
        R4 = x * p * p / W
        Eh = K - (zeta - mu) * chi(U)
        mean = lambda f: float(np.trapezoid(f, t) / T)
        ok = zeta > 1e-150
        row = dict(mu=mu, rho=r['rho'], lambda_c=r['lambda_c'], log_mt=r['log_mt'],
                   convex_projection=bool(n_infl == 0), inflections=n_infl, time_fraction_nonconvex=frac_pos,
                   nonconvex_U_min=float(U[bad].min()) if bad.any() else np.nan,
                   nonconvex_U_max=float(U[bad].max()) if bad.any() else np.nan,
                   nonconvex_p_sign=float(np.sign(p[bad]).mean()) if bad.any() else np.nan,
                   min_zeta_over_eta=float(np.min(ze[ok])), max_zeta_over_eta=float(np.max(ze[ok])),
                   time_fraction_zeta_over_eta_le_1=mean((ze <= 1).astype(float)),
                   R1_minus_A_min=float(np.min(R1[ok]) - A), R1_minus_A_max=float(np.max(R1[ok]) - A),
                   R2_minus_A_min=float(np.min(R2[ok]) - A), R2_minus_A_max=float(np.max(R2[ok]) - A),
                   R3_minus_A_min=float(np.min(R3) - A), R3_minus_A_max=float(np.max(R3) - A),
                   R4_minus_A_min=float(np.min(R4[ok]) - A), R4_minus_A_max=float(np.max(R4[ok]) - A),
                   Ehat_min_over_k=float(Eh.min() / k), Ehat_max_over_k=float(Eh.max() / k),
                   Ehat_time_fraction_negative=mean((Eh < 0).astype(float)),
                   min_K_over_k=float(K.min() / k), min_W=float(W.min()),
                   mean_x_over_y_minus_1=mean(x / y) - 1.0, mean_y_over_lam_minus_p_minus_mu=mean(y / (lam - p)) - mu,
                   LambdaE_min_over_k=float(np.min((y - k) - C.equilibrium_data(mu, k)['zetaE'] * U + C.equilibrium_data(mu, k)['etaE'] * p) / k),
                   max_p_over_lambda=float(p.max() / lam))
    out.append(row)
keys = list(out[0].keys())
with open('cycle_tests.csv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=keys)
    w.writeheader()
    for r in out:
        w.writerow({kk: (repr(float(v)) if isinstance(v, (float, np.floating)) else v) for kk, v in r.items()})
conv = [r for r in out if r['convex_projection']]
print(len(out), 'cycles; convex projection:', len(conv))
print('non-convex cycles (mu, rho-1, inflections, time frac, U-range, p sign):')
for r in out:
    if not r['convex_projection']:
        print(f"  mu={r['mu']:<6} rho-1={r['rho']-1:<10.4g} infl={r['inflections']} frac={r['time_fraction_nonconvex']:.4f} U in [{r['nonconvex_U_min']:.2f},{r['nonconvex_U_max']:.2f}] p-sign={r['nonconvex_p_sign']:+.2f}")
print('S6 on cycles: min zeta/eta <= 1 on every cycle?', all(r['min_zeta_over_eta'] <= 1 for r in out), ' range of min:', min(r['min_zeta_over_eta'] for r in out), max(r['min_zeta_over_eta'] for r in out))
print('   max zeta/eta > 1 on every cycle?', all(r['max_zeta_over_eta'] > 1 for r in out))
for n in ('R1', 'R2', 'R3', 'R4'):
    print(f'   {n}-A sign-changing on every cycle:', all(r[f'{n}_minus_A_min'] < 0 < r[f'{n}_minus_A_max'] for r in out))
print('S5 on cycles: Ehat>0 everywhere on', sum(r['Ehat_time_fraction_negative'] == 0 for r in out), 'of', len(out), '; first rho-1 with a negative part per mu:')
for mu in sorted(set(r['mu'] for r in out)):
    neg = [r['rho'] - 1 for r in out if r['mu'] == mu and r['Ehat_time_fraction_negative'] > 0]
    print(f'   mu={mu}: {min(neg) if neg else None}')
print('identity checks: max |<x/y>-1| =', max(abs(r['mean_x_over_y_minus_1']) for r in out), ' max |<y/(lam-p)>-mu|/mu =', max(abs(r['mean_y_over_lam_minus_p_minus_mu']) / r['mu'] for r in out))
print('min K/k', min(r['min_K_over_k'] for r in out), 'min W', min(r['min_W'] for r in out), 'min Lambda_E/k', min(r['LambdaE_min_over_k'] for r in out), 'max p/lambda', max(r['max_p_over_lambda'] for r in out))
