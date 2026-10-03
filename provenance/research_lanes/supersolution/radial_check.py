#!/usr/bin/env python3
"""Numerical check of the star-shape / convexity theorem and the exact exponent formulas
   S1 = p^2 - U p' > 0,   S1' = (eta-mu) S1 + p K        =>  lambda_c = -<p K / S1>
   S4 = p'^2 - p p'' > 0, S4' = (eta-mu) S4 + p W        =>  lambda_c = -<p W / S4>
   shifted centre q=(-c,0,k), 0<=c<|U_min|:  s_c = p^2 - p'(U+c) > 0,  lambda_c = -<p (K + c zeta)/s_c>
   beta > 1 and kappa > (1-mu) beta (mu<=1) for the strong-stable ratios.
"""
import numpy as np
from base import load, ddt
D = load()
mx = dict(e1=0, e4=0, ec=0, r1=0, r4=0)
print(f"{'mu':>5} {'rho':>6} {'lam_c':>12} {'-<pK/S1>':>12} {'-<pW/S4>':>12} {'c=.5|Umin|':>12} {'c=.9|Umin|':>12} {'minS1/max':>9} {'minS4/max':>9} {'ODE res S1':>10} {'ODE res S4':>10} {'min(beta-1)':>11} {'min(kap-(1-mu)beta)':>19}")
for key in sorted(D):
    P = D[key]; mu = P['mu']; k = P['k']; T = P['T']
    U, p, y, x, eta, zeta = P['U'], P['p'], P['y'], P['x'], P['eta'], P['zeta']
    H = y - k; pd = k - mu*p - y; yd = x - y; pdd = -mu*pd - yd
    S1 = p*p - U*pd; S4 = pd*pd - p*pdd
    K = zeta*U - eta*p - H
    etad = zeta + (mu-1)*eta - eta*eta; zetad = x - (1+eta)*zeta
    W = zetad*p - etad*pd
    l1 = -np.mean(p*K/S1); l4 = -np.mean(p*W/S4)
    outs = []
    for frac in (0.5, 0.9):
        c = frac*abs(U.min()); sc = p*p - pd*(U+c)
        outs.append(-np.mean(p*(K + c*zeta)/sc) if sc.min() > 0 else np.nan)
    # ODE residuals with exact derivative expressions
    S1d = p*pd + mu*U*pd + U*yd
    r1 = np.max(np.abs(S1d - (eta-mu)*S1 - p*K))/np.max(np.abs(S1d))
    pddd = -mu*pdd - (x*p - yd)
    S4d = pd*pdd - p*pddd
    r4 = np.max(np.abs(S4d - (eta-mu)*S4 - p*W))/np.max(np.abs(S4d))
    lam = P['lam_eta']
    print(f"{key[0]:5.2f} {key[1]:6.3f} {lam:+12.5e} {l1:+12.5e} {l4:+12.5e} {outs[0]:+12.5e} {outs[1]:+12.5e} {S1.min()/S1.max():9.2e} {S4.min()/S4.max():9.2e} {r1:10.1e} {r4:10.1e} {P['beta'].min()-1:11.2e} {np.min(P['kappa']-(1-mu)*P['beta']):19.2e}")
