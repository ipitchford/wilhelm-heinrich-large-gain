#!/usr/bin/env python3
"""Referee check 2b (Theorem R, harder corners): shifted centres close to the admissible limit (c = 0.99|U_min|, 0.999|U_min|),
mu > 1 (the theorem is stated for every cycle, not only mu <= 1), and closure of the referee's own integration."""
import sys, time
import numpy as np
ROOT = '/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03'
sys.path.insert(0, ROOT + '/work')
from whlib import find_cycle, kH
from vprofile import profile, load_cycles
C = load_cycles()
cases = [(mu, rho, C[(mu, rho)]) for mu in (0.02, 0.3, 1.0) for rho in (1.01, 1.5, 4.0, 12.0)]
for mu in (5.0, 20.0):
    guess = None
    for rho in (1.5, 1.1, 1.01):
        c = find_cycle(mu, rho*kH(mu), guess=guess); guess = (c['U0'], c['y0']); cases.append((mu, rho, c))
    guess = None
    for rho in (1.5, 2.5, 4.0):
        c = find_cycle(mu, rho*kH(mu), guess=guess); guess = (c['U0'], c['y0'])
        if rho > 1.5: cases.append((mu, rho, c))
print(f"{'mu':>5} {'rho':>5} {'T':>9} {'closure':>8} {'lam_M':>13} {'<eta>-mu':>13} {'-<pK/S1>':>13} {'-<pW/S4>':>13} {'c=.99|Umin|':>13} {'c=.999|Umin|':>13} {'minS1/max':>9} {'minS4/max':>9} {'min sc/max':>10} {'minK':>9} {'minW':>10}")
worst = 0
for mu, rho, cyc in cases:
    k = cyc['k']
    R = profile(mu, k, cyc, cfracs=(0.99, 0.999))
    minsc = min(float(np.min(R[('sc', c)])/np.max(R[('sc', c)])) for c in R['cs'])
    devs = [abs(v/R['lam_eta'] - 1) for v in (R['lam_K'], R['lam_W'], *R['lam_c'])]
    worst = max(worst, *devs)
    print(f"{mu:5.2f} {rho:5.2f} {R['T']:9.3f} {R['closure']:8.1e} {cyc['lambda_c']:+13.6e} {R['lam_eta']:+13.6e} {R['lam_K']:+13.6e} {R['lam_W']:+13.6e} {R['lam_c'][0]:+13.6e} {R['lam_c'][1]:+13.6e} "
          f"{R['S1'].min()/R['S1'].max():9.2e} {R['S4'].min()/R['S4'].max():9.2e} {minsc:10.2e} {R['K'].min():9.2e} {R['W'].min():10.2e}", flush=True)
print("worst relative deviation from <eta>-mu:", worst)
