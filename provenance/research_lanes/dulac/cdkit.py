#!/usr/bin/env python3
"""'Convex-Dulac' certificates on cycles.

Criterion (CD): if F(U,p,y,eta,zeta) satisfies
  (1) F_U nondecreasing in eta and y, nonincreasing in p, independent of zeta,
  (2) F_p nondecreasing in U, y, zeta, independent of eta,
  (3) R' := eta - mu + F_y (x-y) + F_eta eta' + F_zeta zeta' <= 0 on the cycle (not identically 0),
then lambda_c < 0, because  oint F_U dU <= 0 and oint F_p dp <= 0 by the two-point supporting inequalities.

Ansatz A':  F = N(U) eta + M(p) zeta + (n1(U) + m1(p)) y + S0(y),   N, M, n1, m1 nondecreasing.
  R' = eta - mu + N(U) eta' + M(p) zeta' + (n1(U) + m1(p) + s0(y)) (x - y),   s0 = S0'.
All unknown functions are piecewise linear on quantile grids; LP maximises t with R' <= -t.
"""
import sys
import numpy as np
from scipy.optimize import linprog
sys.path.insert(0, '/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03/work')
import whlib


def cycle(mu, rho, n=4001):
    k = rho * mu * (1 + mu)
    c = whlib.find_cycle(mu, k)
    pr = whlib.cycle_profile(mu, k, c, n=n)
    U, p, y, x, eta, zeta = (pr[s][:-1] for s in ('U', 'p', 'y', 'x', 'eta', 'zeta'))
    d = dict(mu=mu, k=k, rho=rho, T=c['T'], lam_c=pr['lambda_c_riccati'], t=pr['t'][:-1], U=U, p=p, y=y, x=x, eta=eta, zeta=zeta,
             K=pr['K'][:-1], W=pr['W'][:-1])
    d['pd'] = k - mu * p - y
    d['etad'] = zeta + (mu - 1) * eta - eta**2
    d['zetad'] = x - (1 + eta) * zeta
    d['yd'] = x - y
    return d


def hat_basis(v, nodes):
    """Piecewise-linear interpolation weights: matrix B with f(v) = B @ f(nodes)."""
    v = np.clip(v, nodes[0], nodes[-1])
    j = np.clip(np.searchsorted(nodes, v, side='right') - 1, 0, len(nodes) - 2)
    w = (v - nodes[j]) / (nodes[j + 1] - nodes[j])
    B = np.zeros((len(v), len(nodes)))
    B[np.arange(len(v)), j] = 1 - w
    B[np.arange(len(v)), j + 1] = w
    return B


def nodes_for(v, m):
    q = np.quantile(v, np.linspace(0, 1, m))
    q = np.unique(q)
    # also make sure that the grid is not degenerate
    return q


def lp_ansatzA(c, m=16, use=('N', 'M', 'n1', 'm1', 's0'), bound=1e3, monotone=True, verbose=False):
    U, p, y = c['U'], c['p'], c['y']
    blocks, names, mono = [], [], []
    def add(name, v, factor, is_mono):
        nd = nodes_for(v, m)
        B = hat_basis(v, nd) * factor[:, None]
        blocks.append(B); names.append((name, nd)); mono.append(is_mono)
    if 'N' in use: add('N', U, c['etad'], True)
    if 'M' in use: add('M', p, c['zetad'], True)
    if 'n1' in use: add('n1', U, c['yd'], True)
    if 'm1' in use: add('m1', p, c['yd'], True)
    if 's0' in use: add('s0', y, c['yd'], False)
    A = np.hstack(blocks)
    nv = A.shape[1]
    # constraints: A c + (eta - mu) + t <= 0
    A_ub = [np.column_stack([A, np.ones(len(U))])]
    b_ub = [-(c['eta'] - c['mu'])]
    # monotonicity:  f_j - f_{j+1} <= 0
    off = 0
    for (name, nd), is_mono in zip(names, mono):
        L = len(nd)
        if is_mono and monotone:
            D = np.zeros((L - 1, nv + 1))
            for j in range(L - 1):
                D[j, off + j] = 1; D[j, off + j + 1] = -1
            A_ub.append(D); b_ub.append(np.zeros(L - 1))
        off += L
    A_ub = np.vstack(A_ub); b_ub = np.concatenate(b_ub)
    cobj = np.zeros(nv + 1); cobj[-1] = -1
    res = linprog(cobj, A_ub=A_ub, b_ub=b_ub, bounds=[(-bound, bound)] * nv + [(None, 10.0)], method='highs')
    if res.status != 0:
        return None, None
    sol = {}
    off = 0
    for (name, nd) in names:
        sol[name] = (nd, res.x[off:off + len(nd)]); off += len(nd)
    return res.x[-1], sol


if __name__ == '__main__':
    from gen_data import SETS
    for mu, rho in SETS:
        c = cycle(mu, rho)
        outs = []
        for use in (('N',), ('M',), ('N', 'M'), ('N', 'M', 's0'), ('N', 'M', 'n1', 'm1', 's0')):
            t, sol = lp_ansatzA(c, use=use)
            outs.append(f'{"+".join(use)}: ' + ('fail' if t is None else f'{t/abs(c["lam_c"]):+.3f}'))
        print(f'mu={mu:<5} rho={rho:<5} lam_c={c["lam_c"]:+.5f}  t*/|lam_c|:  ' + '   '.join(outs), flush=True)
