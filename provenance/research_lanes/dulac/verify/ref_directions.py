#!/usr/bin/env python3
"""Pointwise step in the proof of Theorem D(d) for GENERAL directions (alpha,beta), not only the coordinate ones:
for l = alpha*U + beta*p, at every intermediate level the cycle has exactly two points A (l'>0) and B (l'<0),
q_A - q_B = s*(beta,-alpha) with s<0, and (alpha*eta+beta*zeta)_A < (alpha*eta+beta*zeta)_B.
Double precision (whlib.find_cycle + whlib.cycle_profile, dense profile, linear interpolation in time), so only
moderate rho.  Also evaluates the mean <(alpha p+beta p')(alpha eta+beta zeta)> directly for each direction."""
import sys
import numpy as np
sys.path.insert(0, '/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03/work')
import whlib

rng = np.random.default_rng(20261003)
allok = True
for mu, rho in ((0.3, 1.5), (1.0, 1.1), (1.0, 2.5), (0.1, 2.0), (3.0, 2.0)):
    k = rho * mu * (1 + mu)
    c = whlib.find_cycle(mu, k)
    pr = whlib.cycle_profile(mu, k, c, n=40001)
    t, U, p, y, eta, zeta = (pr[s] for s in ('t', 'U', 'p', 'y', 'eta', 'zeta'))
    pd = k - mu * p - y
    T = c['T']
    worst_s, worst_g, worst_mean, ncross_bad = -np.inf, -np.inf, -np.inf, 0
    for trial in range(40):
        th = rng.uniform(0, 2 * np.pi)
        al, be = np.cos(th), np.sin(th)
        l = al * U + be * p
        ld = al * p + be * pd
        g = al * eta + be * zeta
        mean = np.trapezoid(ld * g, t) / T
        worst_mean = max(worst_mean, mean / (np.trapezoid(np.abs(ld * g), t) / T))
        nsign = int(np.sum(np.diff(np.sign(ld[:-1])) != 0))          # critical points of l on the cycle
        if nsign != 2:
            ncross_bad += 1
        lo, hi = l.min(), l.max()
        for lev in lo + (hi - lo) * np.linspace(0.02, 0.98, 49):
            d = l[:-1] - lev
            idx = np.where(d * (l[1:] - lev) < 0)[0]
            if len(idx) != 2:
                ncross_bad += 1
                continue
            pts = []
            for i in idx:
                w = (lev - l[i]) / (l[i + 1] - l[i])
                pts.append((ld[i], U[i] + w * (U[i + 1] - U[i]), p[i] + w * (p[i + 1] - p[i]), g[i] + w * (g[i + 1] - g[i])))
            A_ = pts[0] if pts[0][0] > 0 else pts[1]
            B_ = pts[1] if pts[0][0] > 0 else pts[0]
            dq = np.array([A_[1] - B_[1], A_[2] - B_[2]])
            s = dq @ np.array([be, -al])                               # q_A - q_B = s (beta,-alpha), |(beta,-alpha)|=1
            worst_s = max(worst_s, s)
            worst_g = max(worst_g, A_[3] - B_[3])
    ok = worst_s < 0 and worst_g < 0 and worst_mean < 0 and ncross_bad == 0
    allok &= ok
    print(f'mu={mu:<4} rho={rho:<4} 40 random directions x 49 levels: max s = {worst_s:+.3e} (<0), max [(a eta+b zeta)_A-(..)_B] = {worst_g:+.3e} (<0), '
          f'max normalised mean = {worst_mean:+.3e} (<0), levels/directions with wrong crossing count = {ncross_bad}  {"OK" if ok else "CHECK"}')
print('ALL OK' if allok else 'FLAGGED')
