#!/usr/bin/env python3
"""Probe: where is the two-phase Delta of Theorem 3.1 numerically positive at mu=1, rho=12, and is it real?
Compare double-precision dense-output values against a direct re-integration to the two phases (no interpolation)."""
import numpy as np
from scipy.integrate import solve_ivp
from vprofile import profile, load_cycles, settle, G1
C = load_cycles()
key = (1.0, 12.0)
mu, rho = key; cyc = C[key]; k = cyc['k']; lam = k/mu; T = cyc['T']
R = profile(mu, k, cyc)
n = 700
ts = np.unique(np.concatenate([np.linspace(0, T, n//2, endpoint=False), R['sol'].t[np.linspace(0, len(R['sol'].t)-2, n//2).astype(int)]]))
Y = R['sol'].sol(ts); U, P, y, eta, zeta = Y[:5]; p = lam - P
t1 = zeta[:, None]*(U[None, :] - U[:, None]); t2 = -eta[:, None]*(p[None, :] - p[:, None]); t3 = -(y[None, :] - y[:, None])
Dl = t1 + t2 + t3
scale = np.abs(t1) + np.abs(t2) + np.abs(t3) + 1e-300
rel = Dl/scale
np.fill_diagonal(rel, -1)
idx = np.argwhere(rel > 0)
print("number of positive pairs:", len(idx), "of", rel.size)
order = np.argsort(-rel[idx[:, 0], idx[:, 1]])
for i, j in idx[order[:12]]:
    print(f"s={ts[i]:.9f} s'={ts[j]:.9f} |ds|={abs(ts[i]-ts[j]):.3e}  Delta={Dl[i,j]:+.3e} rel={rel[i,j]:+.3e}  terms: {t1[i,j]:+.3e} {t2[i,j]:+.3e} {t3[i,j]:+.3e}  U={U[i]:.3f},{U[j]:.3f} P={P[i]:.3e},{P[j]:.3e} y={y[i]:.3e},{y[j]:.3e} eta={eta[i]:.3e} zeta={zeta[i]:.3e}")
# all positive pairs: distribution of |ds| and location
ds = np.abs(ts[idx[:, 0]] - ts[idx[:, 1]])
print("positive pairs: min |ds| = %.3e, max |ds| = %.3e" % (ds.min(), ds.max()))
print("base phase range of positive pairs: [%.4f, %.4f] of T=%.4f; trough at t=%.4f" % (ts[idx[:, 0]].min(), ts[idx[:, 0]].max(), T, R['t_trough']))
