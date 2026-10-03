#!/usr/bin/env python3
"""Which explicit state functions are positive on all computed cycles? (candidate features)"""
import numpy as np
from base import load
D = load()
names = ['S1=p^2-U p\'', 'S2=y\'U-p(y-k)', 'S3=y\'p-p\'(y-k)', 'LamE', 'P*=P+P\'', 'T1=p^2+U(y-k)', 'p^2-p\'U+..']
worst = {}
for key in sorted(D):
    P = D[key]; mu = P['mu']; k = P['k']; lam = k/mu
    U, p, y, x = P['U'], P['p'], P['y'], P['x']
    ed = P['ed']
    pd = k - mu*p - y; yd = x - y; H = y - k
    vals = {
        "S1=p^2-Up'": p*p - U*pd,
        "S2=y'U-p(y-k)": yd*U - p*H,
        "S3=y'p-p'(y-k)": yd*p - pd*H,
        "LamE": H - ed['zetaE']*U + ed['etaE']*p,
        "P*": (1-mu)*(lam-p) + y,
        "UH=U*(y-k)+p^2": U*H + p*p,
        "(x-k)U": (x-k)*U,
        "x-k-zE*U... ": 0*U,
    }
    for nme, v in vals.items():
        r = v.min() / (np.abs(v).max() + 1e-300)
        worst[nme] = min(worst.get(nme, 1), r)
for nme, r in worst.items():
    print(f"{nme:20s} min over cycles of (min/max|.|) = {r:+.3e}")
