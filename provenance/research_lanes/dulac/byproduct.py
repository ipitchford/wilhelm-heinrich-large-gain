#!/usr/bin/env python3
"""By-product: Dulac functions of the simple forms
  (i)   B = exp(beta(U))                       [beta piecewise linear on a 40-node grid]
  (ii)  B = exp(P(U,p)),  P polynomial deg<=3   (regular at E)
  (iii) B = 1/Q(U,p),     Q polynomial, 2<=deg<=4, Q>0 off E   (Cherkas form;  need  DQ - (eta-mu) Q >= 0)
for the reduced planar field u'=v, v'=k-mu v-h(u,v), using the sampled h (H = h-k along W^u(E))."""
import itertools, sys
import numpy as np
from scipy.optimize import linprog
from lpkit import *
from cplp import cplp
from cdkit import hat_basis

def run_i(u, m=40):
    U, p = u['U'], u['p']
    nodes = np.unique(np.quantile(U, np.linspace(0, 1, m)))
    M = hat_basis(U, nodes) * p[:, None]            # D beta(U) = beta'(U) p ; unknown = beta' at nodes
    out = []
    for sign in (-1, 1):
        t, c = cplp(M, u['eta'] - u['mu'], np.ones(len(U)), sign=sign, bound=1e4)
        out.append(t)
    return out

def run_ii(u, deg=3):
    U, p, pd = u['U'], u['p'], u['pd']
    cols = []
    for a in range(deg + 1):
        for b in range(deg + 1 - a):
            if a + b == 0: continue
            col = np.zeros_like(U)
            if a > 0: col += a * U**(a - 1) * p**b * p
            if b > 0: col += b * U**a * p**(b - 1) * pd
            cols.append(col)
    M = np.column_stack(cols)
    return [cplp(M, u['eta'] - u['mu'], np.ones(len(U)), sign=s, bound=1e4)[0] for s in (-1, 1)]

def run_iii(u, degs=(2, 3, 4), n0=6000):
    """maximise t s.t.  DQ - (eta-mu)Q >= t * r^4/(1+r^2),  Q >= r^2,  |coeff| bounded (scaled)."""
    U, p, pd = u['U'], u['p'], u['pd']; q = u['eta'] - u['mu']
    r2 = U * U + p * p
    colsQ, colsD = [], []
    for d in degs:
        for a in range(d + 1):
            b = d - a
            Qm = U**a * p**b
            Dm = np.zeros_like(U)
            if a > 0: Dm += a * U**(a - 1) * p**b * p
            if b > 0: Dm += b * U**a * p**(b - 1) * pd
            colsQ.append(Qm); colsD.append(Dm - q * Qm)
    Qm = np.column_stack(colsQ); Dm = np.column_stack(colsD)
    n, m = Qm.shape
    sc = np.maximum(np.max(np.abs(Qm), axis=0), 1e-300)
    w = r2 * r2 / (1 + r2)
    # constraint generation
    rng = np.random.default_rng(0)
    act = np.zeros(n, dtype=bool); act[rng.choice(n, size=min(n0, n), replace=False)] = True
    cobj = np.zeros(m + 1); cobj[-1] = -1
    for it in range(40):
        idx = np.flatnonzero(act)
        A_ub = np.vstack([np.column_stack([-Dm[idx] / sc, w[idx]]), np.column_stack([-Qm[idx] / sc, np.zeros(len(idx))])])
        b_ub = np.concatenate([np.zeros(len(idx)), -r2[idx]])
        res = linprog(cobj, A_ub=A_ub, b_ub=b_ub, bounds=[(-1e6, 1e6)] * m + [(None, 1e3)], method='highs')
        if res.status != 0:
            return None
        c = res.x[:-1] / sc; t = res.x[-1]
        v1 = t * w - Dm @ c; v2 = r2 - Qm @ c
        viol = np.maximum(v1 / (w + 1e-300) * 0 + v1, v2)
        if viol.max() <= 1e-9 * (1 + abs(t)):
            break
        order = np.argsort(-viol)[:3000]; act[order[viol[order] > 0]] = True
    d = Dm @ c; Qv = Qm @ c
    return float(np.min(d / w)), float(np.min(Qv / r2))

if __name__ == '__main__':
    for s in SETS:
        u = load(*s)
        ti = run_i(u); tii = run_ii(u)
        tiii = run_iii(u) if s[1] <= 2.0 else None
        f = lambda t: 'fail' if t is None else f'{t/u["sig"]:+.3f}'
        print(f'mu={s[0]:<5} rho={s[1]:<5} (i) exp(beta(U)): t*/sig neg={f(ti[0])} pos={f(ti[1])} | (ii) exp(poly3(U,p)): neg={f(tii[0])} pos={f(tii[1])} | '
              f'(iii) 1/Q deg2-4: ' + ('skipped (deep trough)' if s[1] > 2.0 else ('LP fail' if tiii is None else f'min d/w={tiii[0]:+.3e}  min Q/r2={tiii[1]:.3f}')), flush=True)
