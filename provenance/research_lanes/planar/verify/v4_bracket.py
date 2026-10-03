#!/usr/bin/env python3
"""Referee: sign of the bracket in NOTES 2(b) / claim 5:
      B(theta,q) = theta (H(q)-k) + k - H(theta q),    F(theta q)^F(q) = -v B(theta,q).
The notes say 'the bracket is >= 0 by convexity' with no range for theta.  Convexity of H with H(0)=k gives
B >= 0 for 0<=theta<=1 and B <= 0 for theta >= 1; more precisely (supporting plane at q, K = gauge of (A5))
      B(theta,q) <= (1-theta) K(q)   for all theta>=0,   so  B <= -(theta-1) K(q) < 0  for theta > 1, q != 0.
Numerical test on a dense sample of the attractor:
   H_lo(q') = max over sampled supporting planes   (<= H(q'))
   H_up(q') = piecewise-linear interpolation of the sampled graph over a Delaunay triangulation (>= H(q') inside the hull)
so  B_lo = theta(h(q)-k)+k-H_up(theta q) <= B <= B_hi = theta(h(q)-k)+k-H_lo(theta q).
usage: v4_bracket.py mu rho"""
import sys, math, json
import numpy as np
from scipy.spatial import Delaunay
from scipy.interpolate import LinearNDInterpolator
sys.path.insert(0, '/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03/work/lanes/planar/verify')
from vlog import *

mu = float(sys.argv[1]); rho = float(sys.argv[2])
k = rho * mu * (1 + mu); lam = k / mu
ed = eq_data(mu, k)
t_skip = 35.0 / ed['gam']
r0 = min(1e-3, 1e-3 * math.exp(-ed['lam'].real * t_skip))
tconv, cr = spiral_until_converged(mu, k, r0, 0.0)
nph = 2 if len(cr) > 400 else 16
rows = []
for j in range(nph):
    te = np.linspace(t_skip, tconv, 40000 // nph)
    sol = run_log(mu, k, spiral_start(mu, k, r0, 2 * math.pi * j / nph), 0.0, tconv, t_eval=te)
    rows.append(sol.y[:5].T)
S = np.vstack(rows)
U = S[:, 0]; p = lam - np.exp(S[:, 1]); y = np.exp(S[:, 2]); eta = np.exp(S[:, 3]); zeta = np.exp(S[:, 4])
keep = U > -12.0                                  # stay where y, zeta are well above underflow / mesh is meaningful
U, p, y, eta, zeta = U[keep], p[keep], y[keep], eta[keep], zeta[keep]
# add E itself
U = np.append(U, 0.0); p = np.append(p, 0.0); y = np.append(y, k); eta = np.append(eta, ed['eta']); zeta = np.append(zeta, ed['zeta'])
su, sp_ = np.ptp(U), np.ptp(p)
tri = Delaunay(np.column_stack((U / su, p / sp_)))
Hup = LinearNDInterpolator(tri, y)
c0 = y - zeta * U + eta * p                        # plane_a(u,v) = c0 + zeta u - eta v


def Hlo(u, v):
    out = np.full(len(u), -np.inf)
    for i0 in range(0, len(U), 4000):
        sl = slice(i0, i0 + 4000)
        out = np.maximum(out, np.max(c0[sl, None] + zeta[sl, None] * u[None, :] - eta[sl, None] * v[None, :], axis=0))
    return out

rng = np.random.default_rng(7)
idx = rng.choice(len(U) - 1, size=3000, replace=False)
qu, qv, qy, qK = U[idx], p[idx], y[idx], (zeta * U - eta * p - y + k)[idx]
res = []
print(f"mu={mu} rho={rho}: {len(U)} attractor samples, 3000 test points q")
print(" theta   min B_hi      max B_hi      min B_lo (in hull)   frac(B_hi>0)  frac(B_hi<0)   max[B_hi-(1-theta)K(q)]   sign(F(theta q)^F(q)) = -sign(v)?")
for th in (0.05, 0.3, 0.6, 0.9, 0.99, 1.01, 1.1, 1.5, 3.0):
    tu, tv = th * qu, th * qv
    hl = Hlo(tu, tv)
    hu = Hup(tu / su, tv / sp_)
    Bhi = th * (qy - k) + k - hl
    Blo = th * (qy - k) + k - hu
    inh = np.isfinite(hu)
    # wedge with the lower-bound H (sign test of the notes' statement)
    wedge = -qv * Bhi
    agree = float(np.mean(np.sign(wedge) == -np.sign(qv)))
    r = dict(theta=th, Bhi_min=float(Bhi.min()), Bhi_max=float(Bhi.max()),
             Blo_min=float(np.nanmin(Blo[inh])) if inh.any() else None, frac_pos=float(np.mean(Bhi > 0)),
             frac_neg=float(np.mean(Bhi < 0)), plane_bound_excess=float(np.max(Bhi - (1 - th) * qK)), in_hull=int(inh.sum()),
             wedge_sign_agree=agree)
    res.append(r)
    print(f" {th:5.2f}  {r['Bhi_min']:+.4e}  {r['Bhi_max']:+.4e}   {r['Blo_min'] if r['Blo_min'] is not None else float('nan'):+.4e} ({r['in_hull']:4d})"
          f"      {r['frac_pos']:.3f}        {r['frac_neg']:.3f}          {r['plane_bound_excess']:+.2e}              {agree:.3f}")
json.dump(dict(mu=mu, rho=rho, rows=res), open(
    f"/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03/work/lanes/planar/verify/out_v4_{mu:g}_{rho:g}.json", 'w'), indent=1)
