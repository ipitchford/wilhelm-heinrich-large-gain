#!/usr/bin/env python3
"""Cutting-plane LP on the full sample:  maximise t  s.t.  sign*(M c + b) >= t*w  (sign=+1)
                                                       or  M c + b <= -t*w        (sign=-1)."""
from __future__ import annotations
import numpy as np
from scipy.optimize import linprog
from lpkit import features, weight, load, SETS


def build(us, names, polydeg=0, wkind='one', extra=None):
    """Stack feature matrices of several unpacked datasets."""
    Ms, bs, ws = [], [], []
    for u in us:
        f = features(u, polydeg=polydeg)
        if extra is not None:
            f.update(extra(u))
        Ms.append(np.column_stack([f[n] for n in names]))
        bs.append(u['eta'] - u['mu'])
        ws.append(weight(u, wkind))
    return np.vstack(Ms), np.concatenate(bs), np.concatenate(ws)


def cplp(M, b, w, sign=-1, bound=50.0, n0=4000, nadd=3000, maxit=40, tol=1e-9, seed=0, eq=None, tcap=1e3):
    """Returns (t, c).  bound applies to scaled coefficients (|c_i| * max|M_i| <= bound).
    eq = (Aeq, beq) optional equality constraints on c (unscaled)."""
    n, m = M.shape
    scale = np.max(np.abs(M), axis=0); scale[scale == 0] = 1.0
    Ms = M / scale
    rng = np.random.default_rng(seed)
    act = np.zeros(n, dtype=bool)
    act[rng.choice(n, size=min(n0, n), replace=False)] = True
    cobj = np.zeros(m + 1); cobj[-1] = -1.0
    bnds = [(-bound, bound)] * m + [(None, tcap)]
    Aeq = beq = None
    if eq is not None:
        Aeq = np.column_stack([eq[0] / scale, np.zeros(len(eq[1]))]); beq = eq[1]
    t = None; c = None
    for it in range(maxit):
        idx = np.flatnonzero(act)
        if sign < 0:
            A_ub = np.column_stack([Ms[idx], w[idx]]); b_ub = -b[idx]
        else:
            A_ub = np.column_stack([-Ms[idx], w[idx]]); b_ub = b[idx]
        res = linprog(cobj, A_ub=A_ub, b_ub=b_ub, A_eq=Aeq, b_eq=beq, bounds=bnds, method='highs')
        if res.status != 0:
            return None, None
        t = res.x[-1]; cs = res.x[:-1]
        R = Ms @ cs + b
        viol = (R + t * w) if sign < 0 else (t * w - R)       # want <= 0
        vmax = viol.max()
        if vmax <= tol * max(1.0, abs(t)):
            break
        order = np.argsort(-viol)
        new = order[:nadd]
        new = new[viol[new] > 0]
        act[new] = True
    # exact margin on the full sample for the final coefficient vector
    if sign < 0:
        tfull = np.min(-R / w)
    else:
        tfull = np.min(R / w)
    return tfull, cs / scale


if __name__ == '__main__':
    import sys, time
    names = ['logzeta', 'logeta', 'logK', 'logW', 'logL', 'U', 'logy', 'loglamp', 'log1eta', 'logVrad', 'logQtan', 'logQ0']
    for s in SETS:
        u = load(*s)
        M, b, w = build([u], names)
        t0 = time.time()
        t, c = cplp(M, b, w, sign=-1)
        print(s, f't*/sigma={t/u["sig"]:+.4f}', f'[{time.time()-t0:.1f}s]', {n: round(float(v), 3) for n, v in zip(names, c)}, flush=True)
