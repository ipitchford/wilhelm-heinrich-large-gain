#!/usr/bin/env python3
"""Referee: v6 (v) at (0.02,1.01) left a crest deviation 2.3e-4 after t = 60T + 40/|lambda_c|.  Is the orbit still
converging to the same cycle (slow, lambda_c = -1.8e-4) or stuck on something else?  Print crest gaps vs time."""
import sys, math, json
import numpy as np
sys.path.insert(0, '/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03/work/lanes/planar/verify')
from vlog import *
mu, rho = 0.02, 1.01
k = rho * mu * (1 + mu); lam = k / mu
o = json.load(open(f"out_v1b_{mu:g}_{rho:g}.json"))
Uc, lamc, T = o['U_crest'], o['lamc_riccati'], o['T']
rng = np.random.default_rng(11)
crest = make_crest(mu, k)
for trial in range(8):
    X0 = k * math.exp(rng.uniform(-6, 6)); y0 = k * math.exp(rng.uniform(-6, 6)); q0 = lam * math.exp(rng.uniform(-6, 3))
    s0 = np.array([math.log(X0 / k), math.log(q0), math.log(y0), 0.0, 0.0, 0.0])
    sol = run_log(mu, k, s0, 0.0, 2.6e5, events=[crest], rtol=1e-11)
    te = sol.t_events[0]; ue = sol.y_events[0][:, 0]
    pick = [np.searchsorted(te, tt) for tt in (5e4, 1e5, 1.5e5, 2e5, 2.5e5)]
    gaps = [ue[min(j, len(ue) - 1)] - Uc for j in pick]
    rate = math.log(abs(gaps[-1] / gaps[-2])) / 5e4 if gaps[-1] != 0 and gaps[-2] != 0 else float('nan')
    print(f"trial {trial}: start (U,q,y)=({s0[0]:+.2f},{q0:.3g},{y0:.3g}); crest gap at t=5e4..2.5e5: " + ' '.join(f"{g:+.2e}" for g in gaps)
          + f";  decay rate {rate:.3e} (lambda_c = {lamc:.3e})")
