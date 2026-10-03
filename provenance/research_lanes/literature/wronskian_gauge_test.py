#!/usr/bin/env python3
"""Lane 'literature': numerical check of the Wronskian-gauge identities (proved in NOTES.md, section 3.1).
For a periodic C^2 function phi along a cycle put  S_phi = p phi' - p' phi  and  L2 = D^2 + (mu - eta) D + zeta.  Then
      S_phi' = (eta - mu) S_phi + p L2[phi].
 phi = U   : S* = p^2 - U p'   (star-shapedness about E in the (U,p) plane),   L2[U]   = K > 0
 phi = -p' : V  = p'^2 - p p'' (strict convexity of the projected cycle),      L2[-p'] = W > 0
Claims checked on each computed cycle: identity residuals; S* > 0; V > 0; lambda_c = -<pK/S*> = -<pW/V>
and the mixed family  lambda_c = -< p (a K + b W) / (a S* + b V) >  for a few (a,b).   Floating point only."""
import sys, numpy as np
sys.path.insert(0, '/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03/work')
from whlib import find_cycle, cycle_profile, kH
print('  mu    rho   lambda_c      -<pK/S*>      -<pW/V>       mixed(1,1)    mixed(1,10)   minS*/max  minV/max   resid S*   resid V')
worst = 0.0
for mu in [0.05, 0.3, 1.0, 4.0]:
    for rho in [1.05, 1.3, 2.0, 4.0, 8.0]:
        k = rho*kH(mu); c = find_cycle(mu, k); pr = cycle_profile(mu, k, c, n=20001)
        t, U, p, y, eta, zeta, K, W = (pr[s] for s in ('t', 'U', 'p', 'y', 'eta', 'zeta', 'K', 'W'))
        x = k*np.exp(U); T = t[-1]-t[0]
        p1 = k - mu*p - y; p2 = (eta-mu)*p1 - zeta*p
        etad = zeta + (mu-1)*eta - eta**2; zetad = x - (1+eta)*zeta
        p3 = etad*p1 + (eta-mu)*p2 - zetad*p - zeta*p1
        S = p**2 - U*p1; Sd = 2*p*p1 - p*p1 - U*p2          # S*' = p p' - U p''
        V = p1**2 - p*p2; Vd = p1*p2 - p*p3
        rS = np.abs(Sd - (eta-mu)*S - p*K).max()/np.abs(Sd).max(); rV = np.abs(Vd - (eta-mu)*V - p*W).max()/np.abs(Vd).max()
        avg = lambda f: np.trapezoid(f, t)/T
        e1 = -avg(p*K/S); e2 = -avg(p*W/V); e3 = -avg(p*(K+W)/(S+V)); e4 = -avg(p*(K+10*W)/(S+10*V))
        worst = max(worst, abs(e1-c['lambda_c']), abs(e2-c['lambda_c']), abs(e3-c['lambda_c']), abs(e4-c['lambda_c']))
        print('  %-5g %-5g %-13.6e %-13.6e %-13.6e %-13.6e %-13.6e %-10.2e %-10.2e %-10.1e %.1e' % (mu, rho, c['lambda_c'], e1, e2, e3, e4, S.min()/S.max(), V.min()/V.max(), rS, rV))
print('worst |formula - lambda_c| over all cycles and formulas: %.2e' % worst)
