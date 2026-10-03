#!/usr/bin/env python3
"""Compact table of per-cycle averaged optima S*/(-lambda_c) for all feature sets saved by famA2.py."""
import glob, pickle
import numpy as np
from base import load
D = load()
files = sorted(glob.glob('famA2_*.pkl'))
names = [f[len('famA2_'):-4] for f in files]
res = [pickle.load(open(f, 'rb')) for f in files]
print("per-cycle concave-programme optimum  S*/(-lambda_c)   (1 = exact; <=0 = family dead on that cycle)")
print(f"{'mu':>5} {'rho':>6} {'-lam_c':>10} | " + ' | '.join(f"{n:>18}" for n in names))
for key in sorted(D):
    row = []
    for r in res:
        row.append(f"{r[key][0]/(-D[key]['lam_M']):18.3f}" if key in r else f"{'-':>18}")
    print(f"{key[0]:5.2f} {key[1]:6.3f} {-D[key]['lam_M']:10.3e} | " + ' | '.join(row))
print("min over cycles:" + ' ' * 9 + ' | '.join(f"{min(r[k][0]/(-D[k]['lam_M']) for k in r):18.3f}" for r in res))
