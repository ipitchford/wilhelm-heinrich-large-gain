#!/usr/bin/env python3
"""Cycle-only LP:  eta - mu + D(Lambda) - sum_j d_j N_j <= -t,  d_j >= 0, N_j in the monotone cone
(N_j has negative cycle mean by the two-point supporting inequalities):
   type-U:  p*eta*g(U),  p*y*g(U)        g >= 0
   type-p:  p'*zeta*g(p), p'*y*g(p), p'*U*g(p), p'*x*g(p)   g >= 0
Also prints the monotone integrals themselves."""
import sys, itertools
import numpy as np
from scipy.optimize import linprog
from cdkit import cycle
from gen_data import SETS

def feats_L(c, kind):
    mu, k = c['mu'], c['k']; lam = k / mu
    U, p, y, x, eta, zeta, K, W, pd = (c[s] for s in ('U', 'p', 'y', 'x', 'eta', 'zeta', 'K', 'W', 'pd'))
    H = y - k; Hd = x - y
    f = {'U': p, 'p': pd, 'H': Hd}
    if kind >= 1:
        f.update({'x': x * p, 'U2': 2 * U * p, 'Up': p * p + U * pd, 'p2': 2 * p * pd, 'UH': p * H + U * Hd,
                  'pH': pd * H + p * Hd, 'H2': 2 * H * Hd})
    if kind >= 2:
        Vrad = p * p - U * pd
        Qtan = pd * pd - (eta - mu) * p * pd + zeta * p * p
        f.update({'logzeta': x / zeta - 1 - eta, 'logeta': zeta / eta + mu - 1 - eta, 'logy': x / y - 1,
                  'loglamp': y / (lam - p) - mu, 'logVrad': (eta - mu) + p * K / Vrad, 'logQtan': (eta - mu) + p * W / Qtan,
                  'eta': c['etad'], 'zeta': c['zetad']})
    return f

def feats_N(c, kind):
    mu, k = c['mu'], c['k']; lam = k / mu
    U, p, y, x, eta, zeta, pd = (c[s] for s in ('U', 'p', 'y', 'x', 'eta', 'zeta', 'pd'))
    f = {'p*eta': p * eta, 'pd*zeta': pd * zeta}
    if kind >= 1:
        f.update({'p*y': p * y, 'pd*y': pd * y, 'pd*U': pd * U, 'pd*x': pd * x, 'p*eta*x': p * eta * x, 'p*eta/x': p * eta * k / x,
                  'pd*zeta*lamp': pd * zeta * (lam - p), 'pd*zeta/lamp': pd * zeta / (lam - p), 'p*eta2': p * eta**2})
    return f

def run(c, kL, kN, bound=1e4):
    fL = feats_L(c, kL); fN = feats_N(c, kN)
    A = np.column_stack([fL[n] for n in fL] + [-fN[n] for n in fN])
    sc = np.max(np.abs(A), axis=0); A = A / sc
    nL, nN = len(fL), len(fN)
    A_ub = np.column_stack([A, np.ones(len(c['U']))]); b_ub = -(c['eta'] - c['mu'])
    cobj = np.zeros(nL + nN + 1); cobj[-1] = -1
    bnds = [(-bound, bound)] * nL + [(0, bound)] * nN + [(None, 10)]
    res = linprog(cobj, A_ub=A_ub, b_ub=b_ub, bounds=bnds, method='highs')
    if res.status != 0:
        return None, None
    co = dict(zip(list(fL) + list(fN), res.x[:-1] / sc))
    return res.x[-1], co

if __name__ == '__main__':
    for mu, rho in SETS:
        c = cycle(mu, rho)
        T = c['T']; t = c['t']
        mean = lambda f: np.mean(f)            # uniform grid, periodic
        print(f'mu={mu:<5} rho={rho:<5} lam_c={c["lam_c"]:+.5f}  <eta p>={mean(c["eta"]*c["p"]):+.5f}  <zeta pd>={mean(c["zeta"]*c["pd"]):+.5f}', end='  ')
        outs = []
        for kL, kN in ((0, 0), (1, 0), (1, 1), (2, 0), (2, 1)):
            tt, co = run(c, kL, kN)
            outs.append(f'L{kL}N{kN}:' + ('fail' if tt is None else f'{tt/abs(c["lam_c"]):+.3f}'))
        print(' '.join(outs), flush=True)
