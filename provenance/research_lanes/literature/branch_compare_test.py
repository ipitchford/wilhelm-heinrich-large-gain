#!/usr/bin/env python3
"""Lane 'literature' follow-up to curvature_dulac_test.py.
With V = p'^2 - p p'' > 0 (proved in NOTES.md) and dU = p dt:
     lambda_c * T = - oint (W/V) dU = - int_{Umin}^{Umax} [ (W/V)_rising(U) - (W/V)_falling(U) ] dU .
Test A: is the branch difference  d(U) = (W/V)_rising - (W/V)_falling  pointwise positive in U?  (would be a sufficient
        'pointwise' route to stability).  Test B: half-cycle exponents I_+ = int_rising (eta-mu) dt, I_- = int_falling,
        and the proved bounds  I_+ < 2 log r,  I_- > -2 log r,  r = (y_crest-k)/(k-y_trough).
Floating point only."""
import sys, numpy as np
sys.path.insert(0, '/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03/work')
from whlib import find_cycle, cycle_profile, kH
print('  mu    rho    lam_c*T      I_+         I_-         2log r      bounds ok  frac of U-range with d(U)<0   min d/max|d|')
for mu in [0.05, 0.2, 1.0]:
    for rho in [1.05, 1.2, 1.5, 2.0, 3.0, 5.0]:
        k = rho*kH(mu)
        c = find_cycle(mu, k); pr = cycle_profile(mu, k, c, n=40001)
        t, U, p, y, eta, zeta, W = (pr[s] for s in ('t', 'U', 'p', 'y', 'eta', 'zeta', 'W'))
        p1 = k - mu*p - y; p2 = (eta - mu)*p1 - zeta*p; V = p1**2 - p*p2
        itr = int(np.argmin(U))                       # trough index; profile starts at crest, falling first
        fall = slice(0, itr+1); rise = slice(itr, len(t))
        Ip = np.trapezoid((eta-mu)[rise], t[rise]); Im = np.trapezoid((eta-mu)[fall], t[fall])
        r = (y[0]-k)/(k-y[itr]); B = 2*np.log(r)
        Ug = np.linspace(U[itr], U[0], 4003)[1:-1]
        rf = np.interp(Ug, U[fall][::-1], (W/V)[fall][::-1]); rr = np.interp(Ug, U[rise], (W/V)[rise])
        d = rr - rf
        print('  %-5g %-6g %-12.5e %-11.4e %-11.4e %-11.4e %-10s %-29.4f %.3e' % (mu, rho, c['lambda_c']*c['T'], Ip, Im, B, bool(Ip < B and Im > -B), float((d < 0).mean()), d.min()/np.abs(d).max()))
