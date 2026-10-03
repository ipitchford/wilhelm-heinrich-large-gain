#!/usr/bin/env python3
"""Referee: locate cycles at the extreme parameters (whlib.find_cycle with continuation) and cache crest data.
Only the crest point (U0, y0), T and the return-map multiplier are taken from whlib; every profile quantity
is recomputed by the referee's own integrations in vprofile.py."""
import sys, time, pickle
import numpy as np
ROOT = '/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03'
sys.path.insert(0, ROOT + '/work')
from whlib import find_cycle, kH
HERE = ROOT + '/work/lanes/supersolution/verify'
MUS = [0.02, 0.3, 1.0]
RHOS = [1.01, 1.5, 4.0, 12.0]

def build():
    out = {}
    for mu in MUS:
        # continuation ladder: start at 1.5, go down to 1.01, then up to 12
        ladders = [[1.5, 1.2, 1.08, 1.03, 1.01], [1.5, 2.2, 3.0, 4.0, 6.0, 8.0, 10.0, 12.0]]
        for lad in ladders:
            guess = None
            for rho in lad:
                k = rho*kH(mu)
                t0 = time.time()
                try:
                    c = find_cycle(mu, k, guess=guess)
                    if c['U0'] < 1e-4:
                        raise RuntimeError('collapsed to E')
                    guess = (c['U0'], c['y0'])
                    out[(mu, rho)] = c
                    print(f"mu={mu:<5} rho={rho:<5} T={c['T']:10.4f} U0={c['U0']:9.5f} y0/k={c['y0']/k:9.5f} m_t={c['m_t']:.6e} lam_M={c['lambda_c']:+.6e} res={c['residual']:.1e} [{time.time()-t0:.1f}s]", flush=True)
                except Exception as e:
                    print(f"mu={mu} rho={rho} FAILED {type(e).__name__}: {e} [{time.time()-t0:.1f}s]", flush=True)
                    guess = None
    with open(HERE + '/vcycles.pkl', 'wb') as f:
        pickle.dump(out, f)

if __name__ == '__main__':
    build()
