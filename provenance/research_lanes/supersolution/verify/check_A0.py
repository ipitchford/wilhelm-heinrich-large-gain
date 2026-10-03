#!/usr/bin/env python3
"""Referee check 4 (Lemma A0): try to BREAK  lambda_c <= -S[eta~,zeta~]  with random and optimised positive periodic
trial directions, and test concavity of S in the exponents of log-linear families.
S = < min( mu - eta~, Dlog eta~ + 1 - zeta~/eta~, Dlog zeta~ + A - x/zeta~ ) >   (non-uniform trapezoid on a refined grid)."""
import sys, time
import numpy as np
from scipy.optimize import minimize
from vprofile import profile, load_cycles
C = load_cycles()
rng = np.random.default_rng(20261003)
keys = [(0.3, 1.01), (0.3, 1.5), (1.0, 1.01), (1.0, 4.0), (0.02, 1.5), (0.02, 1.01), (0.3, 12.0)]

def grid(R, sub=6):
    ts = R['sol'].t
    tt = np.unique(np.concatenate([ts] + [ts[:-1] + j/sub*np.diff(ts) for j in range(1, sub)]))
    return tt

def setup(key):
    mu, rho = key; cyc = C[key]; k = cyc['k']; lam = k/mu; A = 1+mu; T = cyc['T']
    R = profile(mu, k, cyc)
    tt = grid(R)
    Y = R['sol'].sol(tt)
    U, P, y, eta, zeta = Y[:5]
    p = lam - P; x = k*np.exp(U); pd = mu*P - y; yd = x - y
    Dle = (zeta + (mu-1)*eta - eta*eta)/eta
    Dlz = (x - (1+eta)*zeta)/zeta
    w = 2*np.pi/T
    # features along the cycle and their time derivatives (state functions and pure phase functions)
    # bounded features (so that stiff cycles with U ~ -1e3 do not overflow): 1, tanh U, p/(lam+|p|), (y-k)/(y+k), phases
    q = lam + np.abs(p)
    F = np.array([np.ones_like(tt), np.tanh(U), p/q, (y-k)/(y+k), np.sin(w*tt), np.cos(w*tt), np.sin(2*w*tt)])
    DF = np.array([np.zeros_like(tt), p/np.cosh(U)**2, pd*lam/q**2, 2*k*yd/(y+k)**2, w*np.cos(w*tt), -w*np.sin(w*tt), 2*w*np.cos(2*w*tt)])
    # absolute (state-only) family: features 1, U, log y, log P
    Fa = np.array([np.ones_like(tt), U, np.log(y), np.log(P)])
    DFa = np.array([np.zeros_like(tt), p, yd/y, -pd/P])
    return dict(mu=mu, A=A, T=T, tt=tt, x=x, eta=eta, zeta=zeta, Dle=Dle, Dlz=Dlz, F=F, DF=DF, Fa=Fa, DFa=DFa, lam_c=R['lam_eta'])

def S_rel(D, th):     # trial = true * exp(th.F)
    n = D['F'].shape[0]
    a, b = th[:n] @ D['F'], th[n:] @ D['F']
    e = D['eta']*np.exp(a); z = D['zeta']*np.exp(b)
    g = np.minimum.reduce([D['mu'] - e, D['Dle'] + th[:n] @ D['DF'] + 1 - z/e, D['Dlz'] + th[n:] @ D['DF'] + D['A'] - D['x']/z])
    return np.trapezoid(g, D['tt'])/D['T']

def S_abs(D, th):     # trial = exp(th.Fa)
    n = D['Fa'].shape[0]
    le, lz = th[:n] @ D['Fa'], th[n:] @ D['Fa']
    g = np.minimum.reduce([D['mu'] - np.exp(le), th[:n] @ D['DFa'] + 1 - np.exp(lz - le), th[n:] @ D['DFa'] + D['A'] - D['x']*np.exp(-lz)])
    return np.trapezoid(g, D['tt'])/D['T']

print(f"{'mu':>5} {'rho':>5} {'-lam_c':>12} {'S(true)+lam_c':>13} {'max_rand(S+lam_c)':>17} {'max_opt(S-Strue)':>16} {'max_abs_opt S/(-lam)':>20} {'concav viol rel':>15} {'concav viol abs':>15} {'nrand':>6} {'ncvab':>5} {'s':>5}")
for key in keys:
    t0 = time.time()
    D = setup(key)
    n = D['F'].shape[0]; na = D['Fa'].shape[0]
    L = -D['lam_c']
    s_true = S_rel(D, np.zeros(2*n))
    # (1) random trials at several scales: S + lam_c must be <= 0
    worst = -np.inf; nvalid = 0
    for sc in (1e-3, 1e-2, 1e-1, 0.5, 2.0):
        for _ in range(400):
            th = sc*rng.standard_normal(2*n)
            v = S_rel(D, th)
            assert np.isfinite(v), 'non-finite trial score'
            nvalid += 1
            worst = max(worst, v - L)
    # (2) adversarial: maximise S from random starts (Nelder-Mead then Powell) -- must not exceed -lam_c
    best = -np.inf
    for sc in (0.0, 0.05, 0.5):
        th0 = sc*rng.standard_normal(2*n)
        r = minimize(lambda th: -S_rel(D, th), th0, method='Powell', options=dict(xtol=1e-10, ftol=1e-14, maxiter=4000))
        assert np.isfinite(r.fun)
        best = max(best, -r.fun - s_true)
    # (3) absolute monomial family, optimised
    babs = -np.inf
    for sc in (0.1, 1.0):
        th0 = np.concatenate(([np.log(D['mu']), 0, 0, 0], [np.log(D['mu']), 0, 0, 0])) + sc*rng.standard_normal(2*na)*0.3
        if not np.isfinite(S_abs(D, th0)):
            continue
        r = minimize(lambda th: -S_abs(D, th), th0, method='Powell', options=dict(xtol=1e-9, ftol=1e-13, maxiter=6000))
        if np.isfinite(r.fun):
            babs = max(babs, -r.fun)
    # (4) midpoint concavity
    cv_rel = cv_abs = -np.inf; ncv = 0
    for _ in range(300):
        a, b = rng.standard_normal(2*n), rng.standard_normal(2*n)
        cv_rel = max(cv_rel, 0.5*(S_rel(D, a) + S_rel(D, b)) - S_rel(D, 0.5*(a+b)))
        a = np.array([np.log(D['mu']), 0, 0, 0]*2) + 0.7*rng.standard_normal(2*na)
        b = np.array([np.log(D['mu']), 0, 0, 0]*2) + 0.7*rng.standard_normal(2*na)
        va, vb, vm = S_abs(D, a), S_abs(D, b), S_abs(D, 0.5*(a+b))
        if np.isfinite(va) and np.isfinite(vb) and np.isfinite(vm):
            cv_abs = max(cv_abs, 0.5*(va + vb) - vm); ncv += 1
    print(f"{key[0]:5.2f} {key[1]:5.2f} {L:12.5e} {s_true-L:+13.2e} {worst:+17.3e} {best:+16.3e} {babs/L:20.4f} {cv_rel:+15.2e} {cv_abs:+15.2e} {nvalid:6d} {ncv:5d} {time.time()-t0:5.1f}", flush=True)
print("S(true)+lam_c is the quadrature error of the grid. Lemma A0 is falsified iff max_rand(S+lam_c) or max_opt(S-Strue) is positive beyond that error, or 'max_abs_opt' > 1, or a concavity violation is positive.")
