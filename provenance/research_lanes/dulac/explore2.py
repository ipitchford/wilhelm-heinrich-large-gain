#!/usr/bin/env python3
"""Test: is M = K/V_rad (and W/Q_tan) on the upper half p>0 >= its values on the lower half p<0 at the same U?
If sup_{p<0} M(U,.) <= inf_{p>0} M(U,.) for every U then B = exp(G(U))/V_rad is a Dulac function."""
import numpy as np
from lpkit import *

for mu, rho in SETS:
    u = load(mu, rho)
    U, p, K, W, pd, eta, zeta = (u[s] for s in ('U', 'p', 'K', 'W', 'pd', 'eta', 'zeta'))
    Vrad = p * p - U * pd
    M = K / Vrad
    # restrict to outer region (r > 1e-3) and bin by U quantiles
    r2 = U * U + p * p
    for label, sel in (('outer r>0.05', r2 > 2.5e-3), ('inner 1e-6<r<1e-3', (r2 < 1e-6) & (r2 > 1e-12))):
        Us, ps, Ms = U[sel], p[sel], M[sel]
        if len(Us) < 100:
            continue
        edges = np.quantile(Us, np.linspace(0, 1, 61))
        worst = -np.inf; nb = 0; nviol = 0
        for a, b in zip(edges[:-1], edges[1:]):
            m = (Us >= a) & (Us <= b)
            up = Ms[m & (ps > 0)]; lo = Ms[m & (ps < 0)]
            if len(up) and len(lo):
                gap = lo.max() - up.min()      # want <= 0
                nb += 1; nviol += gap > 0
                worst = max(worst, gap / max(up.min(), 1e-300))
        print(f'mu={mu:<5} rho={rho:<5} {label:<20} bins={nb} violating={nviol} worst rel gap (sup_lower - inf_upper)/inf_upper = {worst:+.3f}')
