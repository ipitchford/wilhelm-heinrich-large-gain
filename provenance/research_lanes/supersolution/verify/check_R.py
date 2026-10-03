#!/usr/bin/env python3
"""Referee check 2 (Theorem R, claim 6 identities): extreme parameters mu in {0.02,0.3,1}, rho in {1.01,1.5,4,12}.
Attempts to falsify S1>0, S4>0, s_c>0 and the exponent formulas; also checks the imported facts K>0, W>0, O=0."""
import sys, time
import numpy as np
from vprofile import profile, load_cycles
C = load_cycles()
keys = [(mu, rho) for mu in (0.02, 0.3, 1.0) for rho in (1.01, 1.5, 4.0, 12.0)]
if len(sys.argv) > 1:
    keys = [k for k in sorted(C)]
print(f"{'mu':>5} {'rho':>5} {'T':>8} {'Umin':>10} {'lam_M(whlib)':>13} {'<eta>-mu':>13} {'-<pK/S1>':>13} {'-<pW/S4>':>13} {'c=.5|Umin|':>13} {'c=.9|Umin|':>13} "
      f"{'minS1/max':>9} {'minS4/max':>9} {'min sc/max':>10} {'minK/max':>9} {'minW/max':>9} {'max|O|rel':>9} {'K ode/def':>9} {'W ode/def':>9} {'<S2/S1>-mu':>10} {'QE-id':>9} {'perdef':>8} {'steps':>6} {'s':>5}")
worst = 0.0
for key in keys:
    mu, rho = key; cyc = C[key]; k = cyc['k']
    t0 = time.time()
    R = profile(mu, k, cyc)
    lamM = cyc['lambda_c']
    minsc = min(float(np.min(R[('sc', c)])/np.max(R[('sc', c)])) for c in R['cs'])
    Orel = float(np.max(np.abs(R['O'])/(np.abs(R['zeta']*R['p']) + np.abs(R['eta']*R['pd']) + np.abs(R['yd']) + 1e-300)))
    Kc = float(np.max(np.abs(R['K_def']/R['K'] - 1)))
    Wc = float(np.max(np.abs(R['W_def']/R['W'] - 1)))
    devs = [abs(v/R['lam_eta'] - 1) for v in (R['lam_K'], R['lam_W'], *R['lam_c'])]
    worst = max(worst, *devs)
    print(f"{mu:5.2f} {rho:5.2f} {R['T']:8.3f} {R['Umin']:10.3f} {lamM:+13.6e} {R['lam_eta']:+13.6e} {R['lam_K']:+13.6e} {R['lam_W']:+13.6e} {R['lam_c'][0]:+13.6e} {R['lam_c'][1]:+13.6e} "
          f"{R['S1'].min()/R['S1'].max():9.2e} {R['S4'].min()/R['S4'].max():9.2e} {minsc:10.2e} {R['K'].min()/R['K'].max():9.2e} {R['W'].min()/R['W'].max():9.2e} {Orel:9.1e} {Kc:9.1e} {Wc:9.1e} "
          f"{R['mean_S2S1']-mu:+10.1e} {R['mean_QE']/R['alpha']-1:+9.1e} {R['per_defect']:8.1e} {R['nsteps']:6d} {time.time()-t0:5.1f}", flush=True)
print("worst relative deviation of the four radial formulas from <eta>-mu:", worst)
