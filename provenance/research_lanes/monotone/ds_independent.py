#!/usr/bin/env python3
"""Independent check of the dual-supporting-plane theorem (no Riccati, no mono_lib objects):
for base phases t0 on a cycle, e_s(t0) := eigenvector of the monodromy Phi_{t0}(T) for the smallest
multiplier (sign: delta U > 0);  D(t0,s) := det( z(s)-z(t0), z'(t0), e_s(t0) ).
Claims:  e_s in cone (+,-,-);  D(t0,s) < 0 for all s != t0;  det(E - z(t0), z'(t0), e_s(t0)) < 0;
         det(z''(t0), z'(t0), e_s(t0)) < 0.
Also checks return-map facts on R={p=0, x>y>k}: DP>0 entrywise and P(R) subset R at random points.
"""
import sys, numpy as np, warnings
from scipy.integrate import solve_ivp
sys.path.insert(0, '/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03/work')
import whlib
from whlib import find_cycle, kH, crest_return
warnings.filterwarnings('ignore')

rng = np.random.default_rng(7)
for mu, rho in [(1.0, 1.1), (1.0, 1.6), (0.4, 1.3), (0.4, 2.2), (0.15, 1.3), (2.5, 1.4), (6.0, 1.8), (0.8, 2.8)]:
    k = rho * kH(mu)
    c = find_cycle(mu, k)
    T = c['T']
    f = whlib.field(mu, k)
    base = solve_ivp(f, [0, T], [c['U0'], 0.0, c['y0']], method='DOP853', rtol=1e-12, atol=1e-14, dense_output=True)
    ss = np.linspace(0, T, 721)[:-1]
    Z = base.sol(ss)
    worst = -np.inf; worstE = -np.inf; worstC = -np.inf; cone_ok = True
    for j in range(0, 720, 30):
        t0 = ss[j]
        z0 = Z[:, j]
        s0 = np.concatenate((z0, np.eye(3).ravel()))
        sol = solve_ivp(whlib._var_rhs(mu, k), [0, T], s0, method='DOP853', rtol=1e-12, atol=1e-14)
        Phi = sol.y[3:12, -1].reshape(3, 3)
        w, V = np.linalg.eig(Phi)
        i = int(np.argmin(np.abs(w)))
        es = V[:, i].real
        es = es / es[0]
        cone_ok &= (es[1] < 0 and es[2] < 0)
        U, p, y = z0
        x = k * np.exp(U)
        zp = np.array([p, k - mu * p - y, x - y])
        zpp = np.array([zp[1], -mu * zp[1] - zp[2], x * p - zp[2]])
        scale = np.linalg.norm(zp) * np.linalg.norm(es)
        for jj in range(720):
            if min((jj - j) % 720, (j - jj) % 720) < 4:
                continue
            d = Z[:, jj] - z0
            worst = max(worst, np.linalg.det(np.array([d, zp, es]).T) / (scale * np.linalg.norm(d)))
        dE = np.array([0, 0, k]) - z0
        worstE = max(worstE, np.linalg.det(np.array([dE, zp, es]).T) / (scale * np.linalg.norm(dE)))
        worstC = max(worstC, np.linalg.det(np.array([zpp, zp, es]).T) / (scale * np.linalg.norm(zpp)))
    # return map on R
    okR = True; npts = 0; unres = 0
    for _ in range(25):
        U0 = c['U0'] * rng.uniform(0.2, 1.8)
        y0 = k + (k * np.exp(U0) - k) * rng.uniform(0.02, 0.98)       # k < y0 < x0
        try:
            U1, y1, _, M = crest_return(mu, k, U0, y0)
        except Exception:
            continue
        if np.abs(M).max() < 1e-6:      # super-contracting excursion: entries below roundoff, sign unresolved
            unres += 1
            okR &= bool(k * np.exp(U1) > y1 > k)
            continue
        npts += 1
        okR &= bool((M > 0).all() and k * np.exp(U1) > y1 > k)
    print(f"mu={mu} rho={rho}: e_s in cone {cone_ok}; max normalised det chord={worst:+.3e} E={worstE:+.3e} curvature={worstC:+.3e} (all must be <0); "
          f"return map on R: DP>0 & P(R) in R at {npts} resolved random points: {okR} ({unres} unresolved: |M|<1e-6, only P(R) in R tested)")
