#!/usr/bin/env python3
"""First-order consequences of convexity+monotonicity of the backward time-T map F_T (coords abc=(U,-p,-(y-k))):
   (C1) (DF_T(z0) - I)(z0 - E)   >= 0      [0 = F_T(E) >= F_T(z0) + DF_T(z0)(E - z0)]
   (C2) (DF_T(z0) - I)(z(s) - z0) <= 0      [z(s) = F_T(z(s)) >= z0 + DF_T(z0)(z(s)-z0)]
   (C3) (DF_T(z0) - I) z''(0)     <= 0      [limit s->0 of (C2)]
 DF_T(z0) = S Phi(T)^{-1} S, S=diag(1,-1,-1), Phi = forward monodromy at the crest phase.
Also: DF_T entrywise positive.
"""
import numpy as np, warnings
from scipy.integrate import solve_ivp
from mono_lib import *
import whlib
warnings.filterwarnings('ignore')

S = np.diag([1.0, -1.0, -1.0])
for mu, rho in [(1.0, 1.05), (1.0, 1.3), (1.0, 2.0), (0.5, 1.1), (0.5, 2.0), (0.2, 1.05), (0.2, 1.5), (3.0, 1.5), (0.3, 1.2), (2.0, 1.1)]:
    k = rho * kH(mu)
    c = find_cycle(mu, k)
    T = c['T']
    s0 = np.concatenate(([c['U0'], 0.0, c['y0']], np.eye(3).ravel()))
    sol = solve_ivp(whlib._var_rhs(mu, k), [0, T], s0, method='DOP853', rtol=1e-12, atol=1e-14, dense_output=True)
    Phi = sol.y[3:12, -1].reshape(3, 3)
    DF = S @ np.linalg.inv(Phi) @ S
    z0 = np.array([c['U0'], 0.0, c['y0']])
    E = np.array([0.0, 0.0, k])
    c1 = (DF - np.eye(3)) @ (S @ (z0 - E))
    x0 = k * np.exp(z0[0]); pp = k - z0[2]; yp = x0 - z0[2]
    zpp = np.array([pp, -mu * pp - yp, -yp])           # (p', p'', y'') at crest (p=0)
    c3 = (DF - np.eye(3)) @ (S @ zpp)
    worst2 = -np.inf
    for sfrac in np.linspace(0.02, 0.98, 49):
        zs = sol.sol(sfrac * T)[:3]
        v = (DF - np.eye(3)) @ (S @ (zs - z0))
        worst2 = max(worst2, (v / np.abs(v).max()).max())
    print(f"mu={mu} rho={rho}: DF>0 {bool((DF > 0).all())}  C1 min rel comp={c1.min()/np.abs(c1).max():+.3e} (>=0)  "
          f"C2 max rel comp={worst2:+.3e} (<=0)  C3 max rel comp={c3.max()/np.abs(c3).max():+.3e} (<=0)")
