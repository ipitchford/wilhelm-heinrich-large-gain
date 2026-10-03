#!/usr/bin/env python3
"""Farkas/dual certificates of infeasibility for small Dulac families.
For a family {Lambda_i} the LP  max t : sum_i c_i f_i(s) + (eta-mu)(s) <= -t  (all samples s)
has dual:  probability measure nu on the samples with  int f_i d nu = 0 (all i),  value  -t* = int (eta-mu) d nu.
If t* < 0 the measure nu is an explicit obstruction: every member of the family has  int (D Lambda + eta - mu) d nu = -t* > 0,
so D Lambda + eta - mu cannot be <= 0 on supp(nu).  We print supp(nu) (radius, quadrant, weight)."""
import sys
import numpy as np
from scipy.optimize import linprog
from lpkit import *

FAMS = {
    'G: gauge logs (zeta,eta,K,W,L,U,y,lamp,1+eta)': ['logzeta', 'logeta', 'logK', 'logW', 'logL', 'U', 'logy', 'loglamp', 'log1eta'],
    'a: logVrad,U,p,H': ['logVrad', 'U', 'p', 'H'],
    'e: 14 logs/linear': ['logVrad', 'logQtan', 'logQ0', 'logzeta', 'logeta', 'logy', 'loglamp', 'log1eta', 'U', 'p', 'H', 'logK', 'logW', 'logL'],
}

def dual(u, names):
    f = features(u)
    M = np.column_stack([f[n] for n in names]); b = u['eta'] - u['mu']
    n, m = M.shape
    sc = np.max(np.abs(M), axis=0); Ms = M / sc
    # dual LP directly: min_nu  -int b dnu ... i.e. maximise int b d nu s.t. Ms^T nu = 0, sum nu = 1, nu >= 0  (value = -t*)
    res = linprog(-b, A_eq=np.vstack([Ms.T, np.ones(n)]), b_eq=np.concatenate([np.zeros(m), [1.0]]), bounds=[(0, None)] * n, method='highs')
    if res.status != 0:
        return None
    nu = res.x
    return -res.fun, nu

if __name__ == '__main__':
    for s in [(0.05, 1.02), (0.05, 1.5), (1.0, 1.5)]:
        u = sub(load(*s), 4)
        U, p, pd = u['U'], u['p'], u['pd']
        r = np.sqrt(U**2 + p**2)
        quad = np.where(p > 0, np.where(pd > 0, 1, 2), np.where(pd <= 0, 3, 4))
        for fam, names in FAMS.items():
            out = dual(u, names)
            if out is None:
                print(f'mu={s[0]} rho={s[1]} {fam}: dual LP failed'); continue
            val, nu = out
            idx = np.flatnonzero(nu > 1e-9)
            print(f'mu={s[0]} rho={s[1]}  family {fam}:  -t*/sigma = int(eta-mu)dnu / sigma = {val/u["sig"]:+.3f}   support size {len(idx)}')
            if val > 0:
                for i in idx[np.argsort(-nu[idx])][:16]:
                    print(f'      weight={nu[i]:.4f}  r={r[i]:.3e}  U={U[i]:+.3e} p={p[i]:+.3e}  Q{quad[i]}  eta-mu={u["eta"][i]-u["mu"]:+.4f}')
