#!/usr/bin/env python3
"""Numerical test of PROGRAM.md (L5): both branches of W^s(E) stay in the backward-invariant octants
{U>0,p<0,y<k} and {U<0,p>0,y>k}; in particular W^s(E) never meets {p=0} away from E."""
import sys
import numpy as np
from scipy.integrate import solve_ivp
sys.path.insert(0, '/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03/work')
from whlib import equilibrium_data, kH

for mu, rho in [(1.0, 1.2), (1.0, 5.0), (0.2, 1.5), (0.05, 3.0), (0.05, 1.02)]:
    k = rho * kH(mu)
    g = equilibrium_data(mu, k)["gammaE"]
    vs = np.array([1.0, -g, k / (1.0 - g)])
    for sgn in (+1, -1):
        z0 = np.array([0, 0, k]) + sgn * 1e-6 * vs / np.linalg.norm(vs)
        f = lambda t, z: [-z[1], -(k - mu * z[1] - z[2]), -(k * np.exp(z[0]) - z[2])]
        ev = lambda t, z: (abs(z[0]) - 30.0) * (z[2] - 1e-9)   # stop at |U| = 30 or y -> 0
        ev.terminal = True
        s = solve_ivp(f, [0, 200], z0, rtol=1e-10, atol=1e-13, events=ev, max_step=0.05)
        U, p, y = s.y
        ok = np.all(sgn * U > 0) and np.all(sgn * p < 0) and np.all(sgn * (y - k) < 0)
        print(f"mu={mu:<5g} rho={rho:<5g} branch {'+' if sgn>0 else '-'}: stays in octant: {bool(ok)}   "
              f"(final U={U[-1]:.3g}, p={p[-1]:.3g}, y-k={y[-1]-k:.3g}, backward time {s.t[-1]:.2f})")
