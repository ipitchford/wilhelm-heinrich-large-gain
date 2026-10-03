#!/usr/bin/env python3
"""Referee: cross-check of the referee integrators (vlog.py) against the v2.0 bundle receipts quoted in
whlib.py's self-test and against whlib.find_cycle itself (whlib is used ONLY here)."""
import sys, math
import numpy as np
sys.path.insert(0, '/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03/work/lanes/planar/verify')
sys.path.insert(0, '/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03/work')
from vlog import *
import whlib
ref = [(0.1, 0.196052776, 38.028889, 0.43642544581), (1.0, 2.201497, 6.502509, 0.77861840302),
       (1.0, 5.841219, 11.121221, 0.030745762), (10.0, 196.052776, 3.802889, 0.43642544581)]
for mu, k, Tref, mref in ref:
    c = whlib.find_cycle(mu, k)
    cyc = newton_cycle(mu, k, c['U0'] * 1.01, c['y0'] * 0.99)      # perturbed seed, own Newton
    prof = cycle_riccati(mu, k, cyc)
    mt = math.exp(prof['S'][5, -1])
    qr = qr_exponents(mu, k, cyc)
    print(f"mu={mu:<5} k={k:<12} T={cyc['T']:.6f} (bundle {Tref:.6f}, whlib {c['T']:.6f})   m_t: Riccati {mt:.9f},"
          f" QR {math.exp(np.sum(qr[-1])):.9f}, bundle {mref:.9f}, whlib {c['m_t']:.9f}")
