#!/usr/bin/env python3
"""Target A machinery.  A certificate is a positive direction field (eta~, zeta~) (= omega2/omega1, omega3/omega1);
the common scalar factor of omega drops out of <min g>.  With omega1 = 1:
   g1 = mu - eta~,   g2 = Dlog eta~ + 1 - zeta~/eta~,   g3 = Dlog zeta~ + A - x/zeta~
and  lambda_c <= -<min_i g_i>.  For the true (eta,zeta) all g_i = mu - eta and <min g> = -lambda_c (the maximum).
Monomial family:  log eta~ = l1 + a1 U + b1 log(y/k) + c1 log(P/lam) + d1 log(Lam_E/k),   P = lam - p,
                  log zeta~ = l2 + a2 U + b2 log(y/k) + c2 log(P/lam) + d2 log(Lam_E/k).
"""
import sys
import numpy as np
from scipy.optimize import minimize
from base import load, avg

def blocks(P):
    mu = P['mu']; k = P['k']; lam = k / mu; A = 1 + mu
    U, p, y, x = P['U'], P['p'], P['y'], P['x']
    ed = P['ed']; gE, eE, zE = ed['gammaE'], ed['etaE'], ed['zetaE']
    Pp = lam - p
    LamE = (y - k) - zE * U + eE * p
    gU = np.expm1(U) - U
    # log-features and their flow derivatives
    F = np.vstack([np.ones_like(U), U, np.log(y / k), np.log(Pp / lam), np.log(LamE / k)])
    DF = np.vstack([np.zeros_like(U), p, x / y - 1, y / Pp - mu, k * gU / LamE - gE])
    return dict(F=F, DF=DF, x=x, mu=mu, A=A, k=k, lam=lam)

def score(theta, B, soft=None, use_d=True):
    th = np.asarray(theta, float).reshape(2, -1)
    F, DF = B['F'], B['DF']
    if not use_d:
        F = F[:4]; DF = DF[:4]
    le = th[0] @ F; lz = th[1] @ F
    et = np.exp(le); zt = np.exp(lz)
    g1 = B['mu'] - et
    g2 = th[0] @ DF + 1 - np.exp(lz - le)
    g3 = th[1] @ DF + B['A'] - B['x'] / zt
    G = np.vstack([g1, g2, g3])
    if soft is None:
        return float(np.mean(G.min(axis=0)))
    m = G.min(axis=0)
    return float(np.mean(m - soft * np.log(np.sum(np.exp(-(G - m) / soft), axis=0))))

def optimise(B, use_d=True, starts=None, verbose=False):
    mu = B['mu']; nf = 5 if use_d else 4
    if starts is None:
        s_Y = np.zeros((2, nf)); s_Y[0, 0] = np.log(mu); s_Y[0, 2] = 1; s_Y[0, 3] = -1
        s_Y[1, 0] = np.log(mu); s_Y[1, 1] = 1; s_Y[1, 3] = -1
        s_c = np.zeros((2, nf)); s_c[:, 0] = np.log(mu)
        starts = [s_Y.ravel(), s_c.ravel(), 0.5 * (s_Y + s_c).ravel()]
    best = (-np.inf, None)
    for s in starts:
        th = np.array(s, float)
        for soft in (0.05, 0.01, 0.002, None):
            f = lambda q: -score(q, B, soft=soft, use_d=use_d)
            r = minimize(f, th, method='Nelder-Mead', options=dict(maxiter=6000, maxfev=6000, xatol=1e-9, fatol=1e-13, adaptive=True))
            th = r.x
            r = minimize(f, th, method='Powell', options=dict(maxiter=4000, xtol=1e-10, ftol=1e-14))
            th = r.x
        val = score(th, B, use_d=use_d)
        if verbose: print('   start ->', val)
        if val > best[0]:
            best = (val, th.copy())
    return best

if __name__ == '__main__':
    D = load()
    use_d = (len(sys.argv) > 1 and sys.argv[1] == 'd')
    sel = sorted(D.keys())
    print(f"monomial family, use Lam_E exponent d: {use_d}")
    print(f"{'mu':>5} {'rho':>5} {'-lam_c':>10} {'<min g> Y':>11} {'<min g> c':>11} {'best':>11} {'best/-lam':>9}  theta(eta~: l,a,b,c[,d] | zeta~)")
    res = {}
    for key in sel:
        P = D[key]; B = blocks(P); mu = P['mu']
        nf = 5 if use_d else 4
        s_Y = np.zeros((2, nf)); s_Y[0, 0] = np.log(mu); s_Y[0, 2] = 1; s_Y[0, 3] = -1
        s_Y[1, 0] = np.log(mu); s_Y[1, 1] = 1; s_Y[1, 3] = -1
        s_c = np.zeros((2, nf)); s_c[:, 0] = np.log(mu)
        val, th = optimise(B, use_d=use_d)
        res[key] = (val, th)
        print(f"{key[0]:5.2f} {key[1]:5.2f} {-P['lam_M']:10.3e} {score(s_Y.ravel(), B, use_d=use_d):+11.3e} {score(s_c.ravel(), B, use_d=use_d):+11.3e} {val:+11.3e} {val/(-P['lam_M']):9.3f}  "
              + np.array2string(th.reshape(2, -1)[0], precision=4, suppress_small=True) + ' | ' + np.array2string(th.reshape(2, -1)[1], precision=4, suppress_small=True), flush=True)
    import pickle
    pickle.dump(res, open('famA_monomial_%s.pkl' % ('d' if use_d else 'nod'), 'wb'))
