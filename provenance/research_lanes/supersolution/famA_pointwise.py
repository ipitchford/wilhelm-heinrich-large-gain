#!/usr/bin/env python3
"""Pointwise (on-cycle) supersolution in log-linear families:
   omega = phi*(1, eta~, zeta~), log phi = th0.F, log eta~ = th1.F, log zeta~ = th2.F
   g_i(t) = th0.DF + [mu - eta~, th1.DF + 1 - zeta~/eta~, th2.DF + A - x/zeta~]
maximise  m(theta) = min_{t,i} g_i(t)   (concave).  m>0  <=>  D omega >= C omega + m*omega pointwise on the cycle,
which forces lambda_c <= -m.  Necessary for any region-based (pointwise) proof within the family."""
import sys
import numpy as np
from scipy.optimize import minimize
from base import load
from cvx import make_blocks, gfun
D = load()
feats = tuple(sys.argv[1].split(','))
nf = 1 + len(feats)
def val_grad(q, B, s):
    t0 = q[:nf-1]; th = q[nf-1:]
    F, DF = B['F'], B['DF']
    G, (e1, e21, e3) = gfun(th, B)
    G = G + (t0 @ DF[1:])[None, :]
    m = G.min()
    E = np.exp(-(G - m)/s); Z = E.sum()
    val = m - s*np.log(Z) + s*np.log(G.size)
    if not np.isfinite(val): return 1e300, np.zeros_like(q)
    pi = E/Z
    g0 = pi.sum(axis=0) @ DF[1:].T
    g1 = -(pi[0]*e1) @ F.T + pi[1] @ DF.T + (pi[1]*e21) @ F.T
    g2 = -(pi[1]*e21) @ F.T + pi[2] @ DF.T + (pi[2]*e3) @ F.T
    return -val, -np.concatenate([g0, g1, g2])
print("features:", feats)
print(f"{'mu':>5} {'rho':>6} {'-lam_c':>10} {'m*=max min g':>13} {'m*/-lam':>8}")
for key in sorted(D):
    P = D[key]; mu = P['mu']
    B = make_blocks(P, feats, nsub=4096)
    q = np.zeros(nf-1 + 2*nf); q[nf-1] = np.log(mu); q[nf-1+nf] = np.log(mu)
    sc = max(abs(P['lam_M']), 1e-4)
    for s in (30, 10, 3, 1, 0.3, 0.1, 0.03, 0.01, 0.003):
        r = minimize(val_grad, q, args=(B, s*sc), jac=True, method='L-BFGS-B', options=dict(maxiter=4000, maxfun=8000, ftol=1e-15, gtol=1e-11))
        q = r.x
    t0 = q[:nf-1]; th = q[nf-1:]
    G, _ = gfun(th, B); G = G + (t0 @ B['DF'][1:])[None, :]
    m = G.min()
    print(f"{key[0]:5.2f} {key[1]:6.3f} {-P['lam_M']:10.3e} {m:+13.4e} {m/(-P['lam_M']):8.3f}", flush=True)
