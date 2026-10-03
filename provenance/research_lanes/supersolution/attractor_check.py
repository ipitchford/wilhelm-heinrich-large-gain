#!/usr/bin/env python3
"""(1) Explicit identities with the E-invariant quadratic form Q_E = p^2 - 2 alpha U p + zeta_E U^2 (2 alpha = eta_E - mu):
        S1 = Q_E + U Lam_E,      D Q_E = 2 alpha Q_E - 2 (p - alpha U) Lam_E   =>  <(p - alpha U) Lam_E / Q_E> = alpha.
   (2) Star-shape S1>0 and local convexity S4>0 on the WHOLE attractor disc: test along orbits on W^u(E)
       (they spiral from E to the innermost cycle), using the pull-back eta,zeta carried by whlib."""
import sys
import numpy as np
from base import load, ROOT
sys.path.insert(0, ROOT + '/work')
from whlib import unstable_manifold_orbit, kH

D = load()
print("(1) cycle identities")
print(f"{'mu':>5} {'rho':>6} {'max|S1-QE-U LamE|/max S1':>26} {'<(p-aU)LamE/QE>':>16} {'alpha':>12} {'<pLamE/S1>':>12} {'-lam_c=<pK/S1>':>14}")
for key in sorted(D):
    P = D[key]; mu = P['mu']; k = P['k']
    U, p, y, x, eta, zeta = P['U'], P['p'], P['y'], P['x'], P['eta'], P['zeta']
    ed = P['ed']; al = 0.5*(ed['etaE'] - mu)
    H = y - k; pd = k - mu*p - y
    S1 = p*p - U*pd
    QE = p*p - 2*al*U*p + ed['zetaE']*U*U
    L = H - ed['zetaE']*U + ed['etaE']*p
    K = zeta*U - eta*p - H
    print(f"{key[0]:5.2f} {key[1]:6.3f} {np.max(np.abs(S1-QE-U*L))/S1.max():26.2e} {np.mean((p-al*U)*L/QE):16.8e} {al:12.8e} {np.mean(p*L/S1):12.4e} {np.mean(p*K/S1):14.4e}")

print("\n(2) W^u(E) spirals: min over the orbit of S1/r^2, S4/|tau|^2, K, W (t-window after transient of the pull-back)")
print(f"{'mu':>5} {'rho':>6} {'theta0':>6} {'min S1/r^2':>11} {'min S4/|tau|^2':>14} {'min K/r^2':>10} {'min W/|tau|^2':>13} {'r range':>22}")
for mu in (0.05, 0.2, 0.5, 1.0):
    for rho in (1.03, 1.3, 2.0, 5.0):
        k = rho*kH(mu)
        for th0 in (0.0, 1.5):
            o = unstable_manifold_orbit(mu, k, r0=1e-6, theta=th0, n=40001)
            U, p, y, x, eta, zeta, K, W = (o[s] for s in ('U', 'p', 'y', 'x', 'eta', 'zeta', 'K', 'W'))
            H = y - k; pd = k - mu*p - y; yd = x - y; pdd = -mu*pd - yd
            r2 = U*U + p*p; t2 = p*p + pd*pd
            S1 = p*p - U*pd; S4 = pd*pd - p*pdd
            # the K,W carried by whlib start at 0 (transient); recompute K from the definition instead
            Kd = zeta*U - eta*p - H
            etad = zeta + (mu-1)*eta - eta*eta; zetad = x - (1+eta)*zeta
            Wd = zetad*p - etad*pd
            m = r2 > 1e-8   # skip the first part where r ~ 1e-6 (round-off in differences)
            print(f"{mu:5.2f} {rho:6.2f} {th0:6.1f} {np.min(S1[m]/r2[m]):11.3e} {np.min(S4[m]/t2[m]):14.3e} {np.min(Kd[m]/r2[m]):10.3e} {np.min(Wd[m]/t2[m]):13.3e} {np.sqrt(r2[m].min()):10.2e} {np.sqrt(r2[m].max()):10.3f}")
