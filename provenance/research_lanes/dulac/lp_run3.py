#!/usr/bin/env python3
"""Disc-wide LP with the two-point 'monotone cone' terms allowed (convex-Dulac criterion):
     eta - mu + D(Lambda) - sum_j d_j N_j <= -t   on all disc samples,   d_j >= 0.
Compares feasibility of small Lambda families with and without the cone terms."""
import time, sys
import numpy as np
from scipy.optimize import linprog
from lpkit import *

def cone_feats(u):
    lam = u['lam']; k = u['k']
    p, pd, eta, zeta, y, U, x = (u[s] for s in ('p', 'pd', 'eta', 'zeta', 'y', 'U', 'x'))
    return {'p*eta': p * eta, 'pd*zeta': pd * zeta, '(p+pd)(eta+zeta)': (p + pd) * (eta + zeta),
            '(p-pd)(eta-zeta)': (p - pd) * (eta - zeta), 'p*y': p * y, 'pd*y': pd * y, 'pd*U': pd * U, 'pd*x': pd * x,
            'p*eta*x': p * eta * x / k, 'p*eta/x': p * eta * k / x, 'pd*zeta*lamp': pd * zeta * (lam - p) / lam}

def cplp_cone(u, names, cone=True, bound=1e4, n0=4000, nadd=3000, maxit=40, seed=0):
    f = features(u)
    cols = [f[n] for n in names]
    cn = []
    if cone:
        cf = cone_feats(u); cn = list(cf); cols += [-cf[n] for n in cn]
    M = np.column_stack(cols); b = u['eta'] - u['mu']
    n, m = M.shape
    sc = np.max(np.abs(M), axis=0); sc[sc == 0] = 1; Ms = M / sc
    rng = np.random.default_rng(seed)
    act = np.zeros(n, dtype=bool); act[rng.choice(n, size=min(n0, n), replace=False)] = True
    cobj = np.zeros(m + 1); cobj[-1] = -1
    bnds = [(-bound, bound)] * len(names) + [(0, bound)] * len(cn) + [(None, 10)]
    for it in range(maxit):
        idx = np.flatnonzero(act)
        res = linprog(cobj, A_ub=np.column_stack([Ms[idx], np.ones(len(idx))]), b_ub=-b[idx], bounds=bnds, method='highs')
        if res.status != 0:
            return None, None
        t = res.x[-1]; cs = res.x[:-1]
        R = Ms @ cs + b
        viol = R + t
        if viol.max() <= 1e-9 * max(1, abs(t)):
            break
        order = np.argsort(-viol)[:nadd]; act[order[viol[order] > 0]] = True
    return float(np.min(-R)), dict(zip(list(names) + cn, cs / sc))

FAMS = {
    'L: U,p,H': ['U', 'p', 'H'],
    'a: Vrad,U,p,H': ['logVrad', 'U', 'p', 'H'],
    'b: Vrad,Qtan,U,p,H': ['logVrad', 'logQtan', 'U', 'p', 'H'],
    'c: b+zeta,eta': ['logVrad', 'logQtan', 'logzeta', 'logeta', 'U', 'p', 'H'],
    'd: c+y,lamp,1eta,Q0': ['logVrad', 'logQtan', 'logQ0', 'logzeta', 'logeta', 'logy', 'loglamp', 'log1eta', 'U', 'p', 'H'],
}
if __name__ == '__main__':
    data = {s: load(*s) for s in SETS}
    for fam, names in FAMS.items():
        print(f'=== {fam} ===', flush=True)
        for s in SETS:
            u = data[s]
            t0, _ = cplp_cone(u, names, cone=False)
            t1, co = cplp_cone(u, names, cone=True)
            cs = ' '.join(f'{n}={v:+.3f}' for n, v in co.items() if abs(v) > 1e-3) if co else ''
            print(f'  mu={s[0]:<5} rho={s[1]:<5} t*/sigma no-cone={t0/u["sig"]:+.3f}  with-cone={t1/u["sig"]:+.3f}   {cs}', flush=True)
