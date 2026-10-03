#!/usr/bin/env python3
"""Near-tangency scan at small mu: shape of r(U) = log(phi(U)/U) on the attractor inside the cycle.
A second cycle inside would need r to reach 0 before U*.  Reports min of r/r(0+) over 0.01 U* < U < 0.9 U*
(interior dip), its location, and whether r is monotone."""
import sys, json
sys.path.insert(0, '/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03/work/lanes/cartography')
import numpy as np
import cart as C
from hunt import primary_cycle

out = []
print('  mu       rho-1    U*       r(0+)       min r/r(0+) on [0.01,0.9]U*   at U/U*   r at 0.9U* / r(0+)   monotone   max rise')
for mu in (0.03, 0.01, 0.003, 0.001, 0.0003):
    for d in (0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0):
        rho = 1 + d
        k = rho * C.kH(mu)
        try:
            cyc, mono, hb, rb = primary_cycle(mu, rho)
        except Exception as e:
            print(mu, d, 'cycle failed', e)
            continue
        Us, ys = cyc['U0'], cyc['y0']
        pr = C.profile(mu, k, Us, ys, cyc['T'], record=False)
        ed = C.equilibrium_data(mu, k)
        rE = 2 * np.pi * ed['alpha'] / ed['omega']
        seeds = Us * np.concatenate((np.geomspace(3e-3, 0.3, 50), np.linspace(0.31, 0.97, 60)))
        Uv, rv = [], []
        for U0 in seeds:
            y0 = k + k * np.expm1(U0) / (1.0 + pr['eta0'])
            U, Y, tt, st = C.crest_iter(mu, k, U0, y0, 3, maxsteps=20_000_000)
            if st != 1:
                continue
            Uv.append(U[2])
            rv.append(np.log(U[3] / U[2]))
        Uv, rv = np.array(Uv), np.array(rv)
        o = np.argsort(Uv)
        Uv, rv = Uv[o], rv[o]
        sel = (Uv > 0.01 * Us) & (Uv < 0.9 * Us)
        i = np.argmin(rv[sel])
        dr = np.diff(rv)
        rise = float(dr.max())
        r09 = float(np.interp(0.9 * Us, Uv, rv))
        row = dict(mu=mu, rho=rho, U_star=Us, r_at_E=rE, min_ratio=float(rv[sel][i] / rE), at=float(Uv[sel][i] / Us),
                   r09_ratio=r09 / rE, monotone=bool(rise < 1e-9), max_rise=rise, log_mt=pr['log_mt'])
        out.append(row)
        print(f"  {mu:<8} {d:<7} {Us:8.4f} {rE:.4e}  {row['min_ratio']:10.4f}                 {row['at']:.3f}     {row['r09_ratio']:8.4f}           {row['monotone']}      {rise:+.3e}", flush=True)
json.dump(out, open('dipscan.json', 'w'), indent=1)
