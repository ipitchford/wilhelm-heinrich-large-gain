#!/usr/bin/env python3
"""Feature library and LP driver for Dulac-function discovery on the sampled attractor disc.

For Lambda = sum_i c_i Lambda_i the Dulac quantity is
        R = D Lambda + eta - mu = sum_i c_i f_i + (eta - mu),     f_i = D Lambda_i  (exact formulas).
LP:  maximise t  s.t.  R(s) + t*w(s) <= 0  for every sample s   (sign=-1),
 or  R(s) - t*w(s) >= 0                                          (sign=+1).
"""
from __future__ import annotations
import itertools
import numpy as np
from scipy.optimize import linprog
from gen_data import SETS, fname, unpack


def load(mu, rho, tag='data'):
    return unpack(np.load(fname(mu, rho, tag)))


def features(u, names=None, polydeg=0, polyvars=('U', 'p', 'H')):
    """dict name -> D(Lambda_i) evaluated on the samples."""
    mu, k, lam, A, gE, etaE, zetaE, sig = (u[s] for s in ('mu', 'k', 'lam', 'A', 'gE', 'etaE', 'zetaE', 'sig'))
    U, p, H, x, y, eta, zeta, K, W, L, pd = (u[s] for s in ('U', 'p', 'H', 'x', 'y', 'eta', 'zeta', 'K', 'W', 'L', 'pd'))
    etad, zetad = u['etad'], u['zetad']
    Hd = x - y
    f = {}
    f['logzeta'] = x / zeta - (1 + eta)
    f['logeta'] = zeta / eta + mu - 1 - eta
    f['logK'] = k * u['psi'] / K - (1 + eta)
    f['logW'] = x * p * p / W - (1 + eta)
    f['logL'] = k * u['phi'] / L - gE
    f['U'] = p
    f['logy'] = x / y - 1
    f['loglamp'] = y / u['lamp'] - mu
    f['log1eta'] = etad / (1 + eta)
    Vrad = p * p - U * pd
    Qtan = pd * (-H) + p * Hd             # = p'^2 - p p'' = p'(k-y) + p(x-y)  (state formula)
    f['logVrad'] = (eta - mu) + p * K / Vrad
    f['logQtan'] = (eta - mu) + p * W / Qtan
    Q0 = p * p - 2 * sig * U * p + zetaE * U * U
    f['logQ0'] = (2 * p * pd - 2 * sig * (p * p + U * pd) + 2 * zetaE * U * p) / Q0
    f['p'] = pd
    f['H'] = Hd
    f['x'] = x * p
    f['eta'] = etad
    f['zeta'] = zetad
    f['p/lamp'] = lam * pd / u['lamp']**2
    Ld = k * u['phi'] - gE * L
    f['H/L'] = (Hd * L - H * Ld) / L**2
    f['K/L'] = (K / L) * (f['logK'] - f['logL'])
    f['W/K'] = (W / K) * (f['logW'] - f['logK'])
    f['K/Vrad'] = (K / Vrad) * (f['logK'] - f['logVrad'])
    f['W/Qtan'] = (W / Qtan) * (f['logW'] - f['logQtan'])
    f['W/Vrad'] = (W / Vrad) * (f['logW'] - f['logVrad'])
    f['Qtan/Vrad'] = (Qtan / Vrad) * (f['logQtan'] - f['logVrad'])
    f['L/Vrad'] = (L / Vrad) * (f['logL'] - f['logVrad'])
    # polynomial state functions
    if polydeg > 0:
        vals = dict(U=U, p=p, H=H)
        ders = dict(U=p, p=pd, H=Hd)
        for deg in range(1, polydeg + 1):
            for combo in itertools.combinations_with_replacement(polyvars, deg):
                name = '*'.join(combo)
                if name in f:
                    continue
                tot = np.zeros_like(U)
                for i in range(deg):
                    term = ders[combo[i]].copy()
                    for j in range(deg):
                        if j != i:
                            term = term * vals[combo[j]]
                    tot += term
                f[name] = tot
    if names is not None:
        f = {n: f[n] for n in names}
    return f


def sub(u, step):
    """Every step-th sample of an unpacked dict."""
    n = len(u['U'])
    return {k_: (v[::step] if isinstance(v, np.ndarray) and v.shape == (n,) else v) for k_, v in u.items()}


def subn(u, nmax=8000):
    return sub(u, max(1, len(u['U']) // nmax))


def weight(u, kind='r2'):
    r2 = u['U']**2 + u['p']**2
    if kind == 'one':
        return np.ones_like(r2)
    if kind == 'r2floor':
        return (r2 + 1e-2) / (1 + r2)
    return r2 / (1 + r2)


def solve_lp(datasets, names, sign=-1, polydeg=0, bound=None, wkind='r2', shared=True, tcap=10.0,
             l1=None, tfix=None, extra=None, verbose=False):
    """datasets: list of unpacked dicts.  If shared, one coefficient vector for all sets.
    Returns (t_opt, coeff dict, per-set max of R/w)."""
    rows, rhs, wts = [], [], []
    featnames = None
    for u in datasets:
        f = features(u, names=None, polydeg=polydeg)
        if names is not None:
            f = {n: f[n] for n in names}
        if extra is not None:
            f.update(extra(u))
        if featnames is None:
            featnames = list(f)
        M = np.column_stack([f[n] for n in featnames])
        rows.append(M); rhs.append(u['eta'] - u['mu']); wts.append(weight(u, wkind))
    M = np.vstack(rows); b = np.concatenate(rhs); w = np.concatenate(wts)
    m = M.shape[1]
    scale = np.max(np.abs(M), axis=0)
    scale[scale == 0] = 1.0
    Ms = M / scale
    # variables: c_scaled (m), t
    if sign < 0:
        A_ub = np.column_stack([Ms, w]); b_ub = -b
    else:
        A_ub = np.column_stack([-Ms, w]); b_ub = b
    cobj = np.zeros(m + 1); cobj[-1] = -1.0
    bnds = [(None, None)] * m + [(None, tcap)]
    if bound is not None:
        bnds = [(-bound * s, bound * s) for s in scale] + [(None, tcap)]
    res = linprog(cobj, A_ub=A_ub, b_ub=b_ub, bounds=bnds, method='highs')
    if res.status != 0:
        return None, None, res.message
    t = res.x[-1]
    c = res.x[:-1] / scale
    if l1 is not None and t > 0:
        # second stage: minimise weighted L1 norm subject to margin >= l1 * t
        tf = l1 * t if tfix is None else tfix
        # variables c+ , c- >= 0
        A2 = np.column_stack([A_ub[:, :m], -A_ub[:, :m]])
        b2 = b_ub - tf * w
        res2 = linprog(np.ones(2 * m), A_ub=A2, b_ub=b2, bounds=[(0, None)] * (2 * m), method='highs')
        if res2.status == 0:
            c = (res2.x[:m] - res2.x[m:]) / scale
            t = tf
    return t, dict(zip(featnames, c)), None


def residual(u, coeffs, polydeg=3):
    f = features(u, polydeg=polydeg)
    R = u['eta'] - u['mu']
    for n, c in coeffs.items():
        R = R + c * f[n]
    return R
