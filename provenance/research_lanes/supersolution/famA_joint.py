#!/usr/bin/env python3
"""Joint (common exponents over all rho at fixed mu) concave programme for monomial families.
maximise over theta:  min_rho  S_rho(theta)/(-lambda_c(rho)).   Alive iff optimum > 0."""
import sys, pickle
import numpy as np
from base import load
from cvx import *
D = load()
feats = tuple(sys.argv[1].split(','))
rmin = float(sys.argv[2]) if len(sys.argv) > 2 else 0.0
rmax = float(sys.argv[3]) if len(sys.argv) > 3 else 1e9
print("features:", feats, " rho range", rmin, rmax)
for mu in sorted(set(k[0] for k in D)):
    keys = [k for k in sorted(D) if k[0] == mu and rmin <= k[1] <= rmax]
    Bs = [make_blocks(D[k], feats, nsub=4096) for k in keys]
    norms = [abs(D[k]['lam_M']) for k in keys]
    best = None
    for t0 in (theta_c(mu, feats),):
        vals, th = joint_maximise(Bs, t0, norms)
        if best is None or vals.min() > best[0].min(): best = (vals, th)
    vals, th = best
    print(f"mu={mu:5.2f}  min ratio = {vals.min():+8.4f}   ratios: " + ' '.join(f"{v:+.3f}" for v in vals))
    print("      theta:", np.array2string(th, precision=3, suppress_small=True, max_line_width=200))
