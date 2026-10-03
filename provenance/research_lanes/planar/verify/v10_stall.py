#!/usr/bin/env python3
"""Referee: (A2) rate remark.  Along the cycle at (0.3,12), log-gap between the pull-back solution started from
(eta,zeta)=(1e-8,1e-8) at the crest and the periodic reference, in LOG variables, at times through the trough
and across the next peak.  Shows that the Hilbert distance stalls while x ~ 0 and collapses at the peak."""
import sys, math, json
import numpy as np
sys.path.insert(0, '/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03/work/lanes/planar/verify')
from vlog import *
for mu, rho in ((0.3, 12.0), (1.0, 4.0)):
    k = rho * mu * (1 + mu)
    o = json.load(open(f"out_v1b_{mu:g}_{rho:g}.json"))
    cyc = newton_cycle(mu, k, o['U_crest'], o['y_crest'])
    T = cyc['T']
    prof = cycle_riccati(mu, k, cyc)
    base = to_log(mu, k, cyc['U0'], 0.0, cyc['y0'], 1.0, 1.0)
    ref = base.copy(); ref[3], ref[4] = prof['S'][3, 0], prof['S'][4, 0]
    tt = np.array([1, 2, 5, 10, 20, 0.5 * T, 0.9 * T, 0.99 * T, 1.0 * T, 1.02 * T, 1.1 * T])
    tt = np.unique(np.clip(tt, 0, 1.1 * T))
    for data in ((1e-8, 1e-8), (1e8, 1e-8)):
        s = base.copy(); s[3], s[4] = math.log(data[0]), math.log(data[1])
        a = run_log(mu, k, s, 0.0, tt[-1], t_eval=tt); b = run_log(mu, k, ref, 0.0, tt[-1], t_eval=tt)
        gap = np.maximum(np.abs(a.y[3] - b.y[3]), np.abs(a.y[4] - b.y[4]))
        print(f"mu={mu} rho={rho} T={T:.2f} data={data}: " + '  '.join(f"t={t:.1f}:{g:.2e}" for t, g in zip(tt, gap)))
