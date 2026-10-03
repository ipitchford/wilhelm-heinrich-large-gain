#!/usr/bin/env python3
"""Per-cycle concave optimum of the monomial families (Target A)."""
import sys, pickle
import numpy as np
from base import load
from cvx import *
D = load()
feats = tuple(sys.argv[1].split(',')) if len(sys.argv) > 1 else ('U', 'y', 'P')
print("features:", feats)
quiet = len(sys.argv) > 2
print(f"{'mu':>5} {'rho':>5} {'-lam_c':>10} {'S(omY)':>10} {'S(omc)':>10} {'S*':>11} {'S*/-lam':>8}  theta*")
res = {}
for key in sorted(D):
    P = D[key]; mu = P['mu']
    B = make_blocks(P, feats)
    tc = theta_c(mu, feats)
    tY = theta_Y(mu, feats) if all(q in feats for q in ('U', 'y', 'P')) else tc
    best = (-np.inf, None)
    for t0 in (tc, tY):
        v, th = maximise(B, t0, scale=max(abs(P['lam_M']), 1e-3) * 10)
        if v > best[0]: best = (v, th)
    v, th = best
    if th is None:
        print(key, 'optimiser failed'); continue
    res[key] = (v, th)
    print(f"{key[0]:5.2f} {key[1]:5.2f} {-P['lam_M']:10.3e} {score(tY,B):+10.2e} {score(tc,B):+10.2e} {v:+11.3e} {v/(-P['lam_M']):8.3f}  "
          + ('' if quiet else np.array2string(th, precision=3, suppress_small=True, max_line_width=200)), flush=True)
pickle.dump(res, open('famA2_%s.pkl' % '_'.join(feats), 'wb'))
