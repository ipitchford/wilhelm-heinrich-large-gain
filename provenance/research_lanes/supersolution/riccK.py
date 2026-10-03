#!/usr/bin/env python3
"""Scalar reduction on the constraint line, parametrised by K = zeta U - eta p - (y-k) > 0:
   eta = (S2 - p K)/S1,  zeta = (S3 - p' K)/S1,   S1 = p^2 - U p', S2 = y'U - p(y-k), S3 = y'p - p'(y-k)
   (R_K)  K' = k G1(U) - (1 + S2/S1) K + (p/S1) K^2,   G1 = U e^U - e^U + 1
   lambda_c = <S2/S1> - mu - <p K / S1>.
Checks: identities, range of K/Lam_E (the two Bregman divergences of h between z and E), sign pattern."""
import numpy as np
from base import load, ddt
D = load()
print(f"{'mu':>5} {'rho':>6} {'lam_c':>11} {'formula':>11} {'resid R_K':>9} {'minS1/max':>9} {'K/LamE min':>10} {'max':>9} {'<pK/S1>':>10} {'<p LamE/S1>':>11} {'<S2/S1>-mu':>10} {'lam(K=LamE)':>11}")
for key in sorted(D):
    P = D[key]; mu = P['mu']; k = P['k']; T = P['T']
    U, p, y, x, eta, zeta = P['U'], P['p'], P['y'], P['x'], P['eta'], P['zeta']
    ed = P['ed']
    pd = k - mu*p - y; yd = x - y; H = y - k
    S1 = p*p - U*pd; S2 = yd*U - p*H; S3 = yd*p - pd*H
    K = zeta*U - eta*p - H
    LamE = H - ed['zetaE']*U + ed['etaE']*p
    G1 = U*np.exp(U) - np.expm1(U)
    lam_formula = np.mean(S2/S1) - mu - np.mean(p*K/S1)
    # residual of R_K using K' = kG1 - (1+eta)K (exact) vs RHS
    rhs = k*G1 - (1 + S2/S1)*K + (p/S1)*K*K
    lhs = k*G1 - (1 + eta)*K
    res = np.max(np.abs(rhs - lhs))/np.max(np.abs(lhs))
    r = K/LamE
    print(f"{key[0]:5.2f} {key[1]:6.3f} {P['lam_eta']:+11.4e} {lam_formula:+11.4e} {res:9.1e} {S1.min()/S1.max():9.2e} {r.min():10.4f} {r.max():9.4f} {np.mean(p*K/S1):+10.3e} {np.mean(p*LamE/S1):+11.3e} {np.mean(S2/S1)-mu:+10.3e} {np.mean(S2/S1)-mu-np.mean(p*LamE/S1):+11.3e}")
