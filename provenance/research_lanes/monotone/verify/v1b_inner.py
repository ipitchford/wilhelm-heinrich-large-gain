#!/usr/bin/env python3
"""Claim 1, last sentence: det(z2 - z(t), z'(t), e_s(t)) < 0 for points z2 of entire solutions strictly INSIDE the cycle.
No second cycle exists at any tested parameter, so the test bodies are (a) the W^u(E) spiral (an entire solution),
(b) orbits started inside on the crest section, after 2 periods (strong contraction e^{-AT} per period onto the disc).
"Insideness" of a sample is delta = 1 - U_nextcrest/U_crest(cycle) (crest amplitude deficit of the inner orbit at its next
crest), and the result is binned by delta, so that samples numerically ON the cycle are not mis-reported."""
import sys, warnings
import numpy as np
from scipy.integrate import solve_ivp
sys.path.insert(0, '/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03/work/lanes/monotone/verify')
from vlib import *
from v1_floquet import SL, atoms, adapted_grid
warnings.filterwarnings('ignore'); np.seterr(all='ignore')

grid = [(mu, rho) for mu in (1.0, 0.3, 0.02) for rho in (1.01, 1.5, 4.0, 12.0)] + [(3.0, 4.0), (50.0, 4.0), (50.0, 1.5), (7.0, 2.5)]
if len(sys.argv) > 1:
    a = sys.argv[1:]; grid = [(float(a[i]), float(a[i + 1])) for i in range(0, len(a), 2)]
for mu, rho in grid:
    k = rho * kH(mu); lam = k / mu
    cyc = find_cycle(mu, k); T = cyc['T']; ed = eq_data(mu, k)

    def gb(tau, w):
        U = cyc['sol'](np.mod(T - tau, T))[0, 0]
        r, ls = w
        return [-(r * r - mu * r - np.exp(ls)), -((r - 1.0) - k * np.exp(U - ls))]
    w0 = np.array([ed['gE'], np.log(ed['gE'] * (ed['gE'] - mu))])
    for _ in range(3):
        sb = solve_ivp(gb, [0, T], w0, method='DOP853', rtol=RT, atol=AT, dense_output=True); w0 = sb.y[:, -1].copy()
    tb = np.unique(np.concatenate((np.linspace(0, T, 201)[:-1], adapted_grid(cyc, 200))))
    Ub, Vb, Qb = cyc['sol'](tb); rb, lsb = sb.sol(T - tb)
    B = atoms(mu, k, Ub[:, None], Vb[:, None], Qb[:, None], rb[:, None], lsb[:, None])
    v = np.real(ed['vec_u']); v = v / np.linalg.norm(v); gr = max(ed['lam_u'].real, 1e-6)
    L = np.log(lam)

    def crest(t, w):
        return w[2] - L
    crest.direction = 1
    fam = []; dl = []

    def add(sol_obj, ts):
        tc_ = sol_obj.t_events[0]; Uc_ = sol_obj.y_events[0][:, 0] if len(tc_) else np.array([])
        idx = np.searchsorted(tc_, ts, side='left')          # next crest at or after ts
        ok = idx < len(tc_)
        fam.append(sol_obj.sol(ts[ok])); dl.append(1.0 - Uc_[idx[ok]] / cyc['U0'])
    eps = 1e-6
    tend = min((np.log(1 / eps) + 10) / gr, 400 * T)
    su = solve_ivp(rhs(mu, k), [0, tend], [eps * v[0], np.log(k + eps * v[2]), np.log(lam - eps * v[1])],
                   method='DOP853', rtol=1e-12, atol=1e-13, dense_output=True, events=crest)
    add(su, np.concatenate((np.linspace(0, tend, 2000), np.linspace(0.7 * tend, tend, 3000))))
    for frac in (0.02, 0.1, 0.3, 0.5, 0.7, 0.9):
        a_ = frac * cyc['U0']
        sv_ = solve_ivp(rhs(mu, k), [0, 14 * T], [a_, np.log(k + ed['zetaE'] * a_), np.log(lam)],
                        method='DOP853', rtol=1e-12, atol=1e-13, dense_output=True, events=crest)
        add(sv_, np.linspace(2 * T, 14 * T, 1500))
    W = np.concatenate(fam, axis=1); dist = np.concatenate(dl)
    Uu, Vu, Qu = W
    d_U = SL.real(Uu[None, :] - Ub[:, None])
    d_p = SL.pos(Qb[:, None] + 0 * Qu[None, :]) - SL.pos(Qu[None, :] + 0 * Qb[:, None])
    d_y = SL.pos(Vu[None, :] + 0 * Vb[:, None]) - SL.pos(Vb[:, None] + 0 * Vu[None, :])
    ratio = (d_U * B['q1'] + d_p * B['q2'] + d_y * B['q3']).ev()[0]
    colmax = ratio.max(axis=0)
    line = f"mu={mu} rho={rho}: inner points {W.shape[1]}, bases {len(tb)};"
    for lo, hi in [(1e-1, np.inf), (1e-3, 1e-1), (1e-5, 1e-3), (1e-7, 1e-5), (-np.inf, 1e-7)]:
        m = (dist >= lo) & (dist < hi)
        if m.any():
            line += f"  delta[{lo:g},{hi:g}): n={m.sum()} viol={(colmax[m] > 0).sum()} max={colmax[m].max():+.2e};"
    print(line, flush=True)
