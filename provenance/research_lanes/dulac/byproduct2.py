#!/usr/bin/env python3
"""By-product (iii): Cherkas-type inverse-polynomial Dulac function B = 1/Q(U,p) on the annulus r >= rmin of the attractor disc.
   need  Q > 0  and  d := DQ - (eta-mu) Q >= 0.   LP: maximise t s.t. d >= t r^2, Q >= r^2 (rows normalised by r^2)."""
import sys
import numpy as np
from scipy.optimize import linprog
from lpkit import *

def run(u, degs, rmin):
    U, p, pd = u['U'], u['p'], u['pd']; q = u['eta'] - u['mu']
    r2 = U * U + p * p
    sel = r2 >= rmin**2
    U, p, pd, q, r2 = U[sel], p[sel], pd[sel], q[sel], r2[sel]
    colsQ, colsD = [], []
    for d in degs:
        for a in range(d + 1):
            b = d - a
            Qm = U**a * p**b
            Dm = np.zeros_like(U)
            if a > 0: Dm += a * U**(a - 1) * p**b * p
            if b > 0: Dm += b * U**a * p**(b - 1) * pd
            colsQ.append(Qm / r2); colsD.append((Dm - q * Qm) / r2)
    Qm = np.column_stack(colsQ); Dm = np.column_stack(colsD)
    n, m = Qm.shape
    sc = np.maximum(np.max(np.abs(Qm), axis=0), np.max(np.abs(Dm), axis=0))
    step = max(1, n // 12000)
    idx = np.arange(0, n, step)
    cobj = np.zeros(m + 1); cobj[-1] = -1
    A_ub = np.vstack([np.column_stack([-Dm[idx] / sc, np.ones(len(idx))]), np.column_stack([-Qm[idx] / sc, np.zeros(len(idx))])])
    b_ub = np.concatenate([np.zeros(len(idx)), -np.ones(len(idx))])
    res = linprog(cobj, A_ub=A_ub, b_ub=b_ub, bounds=[(-1e5, 1e5)] * m + [(None, 1e3)], method='highs')
    if res.status != 0:
        return None
    c = res.x[:-1] / sc
    return float(np.min(Dm @ c)), float(np.min(Qm @ c))

if __name__ == '__main__':
    for s in SETS:
        if s[1] > 2.0: continue
        u = load(*s)
        outs = []
        for degs in ((2,), (2, 3), (2, 3, 4), (2, 3, 4, 5, 6)):
            for rmin in (1e-2,):
                o = run(u, degs, rmin)
                outs.append(f'deg{degs[0]}-{degs[-1]}: ' + ('LP fail' if o is None else f't={o[0]/u["sig"]:+.3f}sig minQ/r2={o[1]:.2f}'))
        print(f'mu={s[0]:<5} rho={s[1]:<5} annulus r>=0.01:  ' + ' | '.join(outs), flush=True)
