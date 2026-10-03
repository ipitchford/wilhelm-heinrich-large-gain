#!/usr/bin/env python3
"""Referee side-finding: does the shared whlib (original variables, absolute tolerance 1e-13) keep the Riccati pair
positive through deep troughs?  (Not a statement of the lane; relevant to any lane using whlib at rho >= 4.)"""
import sys, json, warnings
import numpy as np
warnings.filterwarnings('ignore')
sys.path.insert(0, '/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03/work')
import whlib
for mu, rho in ((1, 4), (0.3, 4), (0.3, 12), (0.02, 4)):
    k = rho * mu * (1 + mu)
    o = json.load(open(f"out_v1b_{mu:g}_{rho:g}.json"))
    try:
        c = whlib.find_cycle(mu, k, guess=(o['U_crest'], o['y_crest']))
        pr = whlib.cycle_profile(mu, k, c)
        print(f"mu={mu} rho={rho}: whlib lambda_c(M)={c['lambda_c']:.9f} lambda_c(Riccati)={pr['lambda_c_riccati']:.9f}"
              f" [referee log-variables: {o['lamc_riccati']:.9f}];  min eta={pr['eta'].min():.3e} min zeta={pr['zeta'].min():.3e}"
              f" min K={pr['K'].min():.3e} min W={pr['W'].min():.3e}; true min eta = e^{o['bounds_log']['a_min']:.1f},"
              f" true min zeta = e^{o['bounds_log']['b_min']:.1f}")
    except Exception as e:
        print(f"mu={mu} rho={rho}: whlib failed: {type(e).__name__}: {e}")
