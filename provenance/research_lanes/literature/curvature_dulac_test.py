#!/usr/bin/env python3
"""Lane 'literature': test of the Gasull-Giacomini 'orbit-curvature' Dulac-Cherkas function (arXiv:2101.03874, Thm 2.3)
on the reduced planar Wilhelm-Heinrich flow  u' = v, v' = Q := k - mu v - h(u,v).

Their function  V := Q*Pdot - P*Qdot  becomes  V = p'^2 - p p''  and their  M_1 = D  becomes (derived in NOTES.md)
        D = p * W,      W = zeta' p - eta' p' > 0   (manuscript Lemma 3.3),
so D changes sign with p and Theorem 2.3 does NOT apply.  By-product (exact identity on any cycle with V > 0):
        (log V)' = (eta - mu) + p W / V     =>     lambda_c = <eta> - mu = - < p W / V >.
Checks: (i) residual of the pointwise identity  V' - (eta-mu) V - p W = 0 ;  (ii) sign of V on the cycle
(V > 0 <=> projected cycle in the (U,p) plane is strictly convex);  (iii) lambda_c vs -<pW/V>.   Floating point only.
"""
import sys, numpy as np
sys.path.insert(0, '/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03/work')
from whlib import find_cycle, cycle_profile, kH

print('  mu    rho    lambda_c(M)   <eta>-mu      -<pW/V>       min V/max V   V>0   max|ident. resid|/scale')
for mu in [0.05, 0.2, 1.0]:
    for rho in [1.05, 1.2, 1.5, 2.0, 3.0, 5.0, 8.0]:
        k = rho*kH(mu)
        try:
            c = find_cycle(mu, k); pr = cycle_profile(mu, k, c, n=20001)
        except Exception as e:
            print('  %-5g %-6g FAILED %r' % (mu, rho, e)); continue
        t, U, p, y, eta, zeta, W = (pr[s] for s in ('t', 'U', 'p', 'y', 'eta', 'zeta', 'W'))
        x = k*np.exp(U)
        p1 = k - mu*p - y
        p2 = (eta - mu)*p1 - zeta*p            # tangent equation (2.8) of the manuscript
        V = p1**2 - p*p2
        # analytic derivative of V along the flow: V' = p' p'' - p p'''
        etad = zeta + (mu - 1)*eta - eta**2; zetad = x - (1 + eta)*zeta
        p3 = etad*p1 + (eta - mu)*p2 - zetad*p - zeta*p1
        Vd = p1*p2 - p*p3
        resid = np.abs(Vd - (eta - mu)*V - p*W).max()/np.abs(Vd).max()
        T = t[-1] - t[0]
        if (V > 0).all():
            est = -np.trapezoid(p*W/V, t)/T
        else:
            est = np.nan
        print('  %-5g %-6g %-13.6e %-13.6e %-13.6e %-13.3e %-5s %.2e' % (mu, rho, c['lambda_c'], pr['lambda_c_riccati'], est, V.min()/V.max(), bool((V > 0).all()), resid))
