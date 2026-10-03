#!/usr/bin/env python3
"""Concave-programming engine for log-linear certificate families (Target A).

Family:  log eta~ = th1 . F(z),  log zeta~ = th2 . F(z)   (F_0 = 1), features F_j explicit functions of the state
with flow derivatives DF_j.  With omega1 = 1:
   g1 = mu - exp(th1.F),  g2 = th1.DF + 1 - exp((th2-th1).F),  g3 = th2.DF + A - x exp(-th2.F).
Each g_i is CONCAVE in theta, so  S(theta) = <min_i g_i>  is concave: the per-cycle optimum and the joint
(common-theta) max-min over several cycles are convex programs -> no local optima; a negative optimum is a
genuine refutation of the family on that cycle (up to solver tolerance).
"""
import numpy as np
from scipy.optimize import minimize

def gfun(th, B):
    nf = B['F'].shape[0]
    t1, t2 = th[:nf], th[nf:]
    F, DF = B['F'], B['DF']
    le = t1 @ F; lz = t2 @ F
    e1 = np.exp(np.clip(le, -300, 300)); e21 = np.exp(np.clip(lz - le, -300, 300)); e3 = np.exp(np.clip(B['logx'] - lz, -300, 300))
    G = np.vstack([B['mu'] - e1, t1 @ DF + 1 - e21, t2 @ DF + B['A'] - e3])
    return G, (e1, e21, e3)

def score(th, B):
    G, _ = gfun(th, B)
    return float(np.sum(B['wq'] * G.min(axis=0)))

def soft_score_grad(th, B, s):
    nf = B['F'].shape[0]
    F, DF = B['F'], B['DF']
    G, (e1, e21, e3) = gfun(th, B)
    m = G.min(axis=0)
    E = np.exp(-(G - m) / s)
    Z = E.sum(axis=0)
    val = np.sum(B['wq'] * (m - s * np.log(Z)))
    if not np.isfinite(val):
        return -1e300, np.zeros_like(th)
    pi = E / Z * B['wq']
    g_t1 = -(pi[0] * e1) @ F.T + pi[1] @ DF.T + (pi[1] * e21) @ F.T
    g_t2 = -(pi[1] * e21) @ F.T + pi[2] @ DF.T + (pi[2] * e3) @ F.T
    return val, np.concatenate([g_t1, g_t2])

def maximise(B, th0, softs=(0.3, 0.1, 0.03, 0.01, 0.003, 0.001, 3e-4, 1e-4), scale=None):
    th = np.array(th0, float)
    sc = scale if scale is not None else 1.0
    for s in softs:
        ss = s * sc
        f = lambda q: tuple(-v for v in soft_score_grad(q, B, ss))
        r = minimize(lambda q: (lambda v, g: (-v, -g))(*soft_score_grad(q, B, ss)), th, jac=True, method='L-BFGS-B',
                     options=dict(maxiter=3000, maxfun=6000, ftol=1e-15, gtol=1e-10))
        th = r.x
    return score(th, B), th

def joint_maximise(Bs, th0, norms, softs=(0.3, 0.1, 0.03, 0.01, 0.003, 0.001, 3e-4)):
    """maximise  min_c  score_c(theta)/norm_c   (soft-min over cycles as well)"""
    th = np.array(th0, float)
    for s in softs:
        def f(q):
            vals = []; grads = []
            for B, nm in zip(Bs, norms):
                v, g = soft_score_grad(q, B, s * nm)
                vals.append(v / nm); grads.append(g / nm)
            vals = np.array(vals); grads = np.array(grads)
            m = vals.min(); E = np.exp(-(vals - m) / s); Z = E.sum()
            return -(m - s * np.log(Z)), -(E / Z) @ grads
        r = minimize(f, th, jac=True, method='L-BFGS-B', options=dict(maxiter=4000, maxfun=8000, ftol=1e-15, gtol=1e-10))
        th = r.x
    return np.array([score(th, B) / nm for B, nm in zip(Bs, norms)]), th

def make_blocks(P, feats=('U', 'y', 'P'), nsub=None):
    """feature blocks on a cycle profile. feats subset of: U, y, P, LamE, star(=p^2-U p'), X2 ..."""
    mu = P['mu']; k = P['k']; lam = k / mu; A = 1 + mu
    U, p, y, x = P['U'], P['p'], P['y'], P['x']
    ed = P['ed']; gE, eE, zE = ed['gammaE'], ed['etaE'], ed['zetaE']
    Pp = lam - p; pd = k - mu * p - y; yd = x - y
    F = [np.ones_like(U)]; DF = [np.zeros_like(U)]
    for f in feats:
        if f == 'U':
            F.append(U); DF.append(p)
        elif f == 'y':
            F.append(P['logZ']); DF.append(np.expm1(U - P['logZ']))
        elif f == 'P':
            F.append(P['logY']); DF.append(mu * np.expm1(P['logZ'] - P['logY']))
        elif f == 'LamE':
            L = (y - k) - zE * U + eE * p
            F.append(np.log(L / k)); DF.append(k * (np.expm1(U) - U) / L - gE)
        elif f == 'p':      # non-monomial: exp(c p)
            F.append(p); DF.append(pd)
        elif f == 'ylin':   # exp(c (y-k)/k)
            F.append((y - k) / k); DF.append(yd / k)
        elif f == 'xlin':   # exp(c (x-k)/k)
            F.append((x - k) / k); DF.append(x * p / k)
        elif f == 'Pstar':  # log P*, P* = P + P' = (1-mu)P + y ;  P*' = -mu P* + x
            Ps = (1 - mu) * Pp + y
            F.append(np.log(Ps / lam)); DF.append(x / Ps - mu)
        elif f == 'S3':     # log S3, S3 = y'p - p'(y-k) > 0 numerically
            H = y - k
            S3 = yd * p - pd * H
            pdd = -mu * pd - yd; ydd = x * p - yd
            F.append(np.log(S3)); DF.append((ydd * p + yd * pd - pdd * H - pd * yd) / S3)
        elif f == 'QE':     # log Q_E, Q_E = p^2 - 2 al U p + zE U^2 (positive definite); D log Q_E = 2 al - 2 (p - al U) LamE / Q_E
            al = 0.5 * (eE - mu)
            QE = p * p - 2 * al * U * p + zE * U * U
            L = (y - k) - zE * U + eE * p
            F.append(np.log(QE)); DF.append(2 * al - 2 * (p - al * U) * L / QE)
        elif f == 'S4':     # log S4, S4 = p'^2 - p p'' > 0 (proved)
            pdd = -mu * pd - yd
            pddd = -mu * pdd - (x * p - yd)
            S4 = pd * pd - p * pdd
            F.append(np.log(S4)); DF.append((pd * pdd - p * pddd) / S4)
        elif f == 'star':   # log(p^2 - U p')  (needs positivity)
            S = p * p - U * pd
            pdd = -mu * pd - yd
            F.append(np.log(S)); DF.append((2 * p * pd - p * pd - U * pdd) / S)
        else:
            raise ValueError(f)
    F = np.vstack(F); DF = np.vstack(DF)
    n = len(U)
    idx = np.arange(n) if nsub is None else np.arange(0, n, max(1, n // nsub))
    wq = np.full(len(idx), 1.0 / len(idx))
    return dict(F=F[:, idx], DF=DF[:, idx], x=x[idx], logx=np.log(k) + U[idx], mu=mu, A=A, k=k, lam=lam, wq=wq)

def theta_Y(mu, feats):
    nf = 1 + len(feats); t = np.zeros((2, nf)); t[:, 0] = np.log(mu)
    fl = list(feats)
    t[0, 1 + fl.index('y')] = 1; t[0, 1 + fl.index('P')] = -1
    t[1, 1 + fl.index('U')] = 1; t[1, 1 + fl.index('P')] = -1
    return t.ravel()

def theta_c(mu, feats):
    nf = 1 + len(feats); t = np.zeros((2, nf)); t[:, 0] = np.log(mu)
    return t.ravel()
