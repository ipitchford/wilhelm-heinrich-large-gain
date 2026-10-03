#!/usr/bin/env python3
"""Averaged score <min_i g_i> (must be >0 for a certificate; maximum possible = -lambda_c) of explicit,
parameter-free trial directions on the constraint line  zeta p - eta p' = y':
   T_K(c):  K <- c*Lam_E :  eta~ = (S2 - c p Lam_E)/S1,  zeta~ = (S3 - c p' Lam_E)/S1   (first-order exact at E for c=1)
and of the chain trials T_Y=(y/P,x/P), T_c=(mu,mu).  Dlog computed spectrally (rho<=3 only)."""
import numpy as np
from base import load, ddt
D = load()
def score_dir(et, zt, P):
    mu = P['mu']; A = 1+mu; T = P['T']
    if et.min() <= 0 or zt.min() <= 0: return np.nan
    g1 = mu - et
    g2 = ddt(np.log(et), T) + 1 - zt/et
    g3 = ddt(np.log(zt), T) + A - P['x']/zt
    return np.mean(np.minimum(np.minimum(g1, g2), g3))
print(f"{'mu':>5} {'rho':>6} {'-lam_c':>10} {'T_K(1)':>11} {'ratio':>8} {'T_K(.5)':>11} {'T_Y':>10} {'T_c':>10}  {'<eta~_K>-mu':>12}")
for key in sorted(D):
    if key[1] > 3.01: continue
    P = D[key]; mu = P['mu']; k = P['k']; lam = k/mu
    U, p, y, x = P['U'], P['p'], P['y'], P['x']
    ed = P['ed']; H = y-k; pd = k-mu*p-y; yd = x-y
    S1 = p*p-U*pd; S2 = yd*U-p*H; S3 = yd*p-pd*H
    L = H - ed['zetaE']*U + ed['etaE']*p
    out = []
    for c in (1.0, 0.5):
        out.append(score_dir((S2-c*p*L)/S1, (S3-c*pd*L)/S1, P))
    sY = score_dir(y/(lam-p), x/(lam-p), P)
    sc = np.mean(np.minimum(0, (mu*(1+mu)-x)/mu))
    print(f"{key[0]:5.2f} {key[1]:6.3f} {-P['lam_M']:10.3e} {out[0]:+11.3e} {out[0]/(-P['lam_M']):8.2f} {out[1]:+11.3e} {sY:+10.2e} {sc:+10.2e}  {np.mean((S2-p*L)/S1)-mu:+12.3e}")
