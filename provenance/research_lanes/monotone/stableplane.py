#!/usr/bin/env python3
"""Test: does the stable plane span(z'(t), e_s(t)) through z(t) support the whole cycle (and E)?
beta_s(t) = Psi^(t).(z(t+s)-z(t)) ; sign for all s in (0,T), all t.
Also the family of identities  lam <beta_s> = - <Psi3 x c(U(.+s)-U)>,  c(d)=e^d-1-d.
"""
import numpy as np, warnings
from mono_lib import *
warnings.filterwarnings('ignore')

cases = [(1.0, 1.05), (1.0, 1.3), (1.0, 2.0), (0.5, 1.1), (0.5, 2.0), (0.2, 1.05), (0.2, 1.5),
         (0.05, 1.2), (3.0, 1.5), (5.0, 2.5), (1.0, 3.0), (0.2, 2.5)]
for mu, rho in cases:
    k = rho * kH(mu)
    c = find_cycle(mu, k)
    n = 2001
    d = floquet_data(mu, k, c, n=n)
    t = d['t']
    P = d['Psi']
    Z = np.array([d['U'], d['p'], d['y']])[:, :-1]      # drop duplicate endpoint
    Pp = P[:, :-1]
    m = n - 1
    worst_pos = -np.inf; worst_neg = np.inf
    nshift = 200
    ident_err = 0.0
    signs = []
    for j in range(1, nshift):
        sh = j * m // nshift
        Zs = np.roll(Z, -sh, axis=1)
        beta = np.sum(Pp * (Zs - Z), axis=0)
        worst_pos = max(worst_pos, beta.max()); worst_neg = min(worst_neg, beta.min())
        dU = Zs[0] - Z[0]
        w = d['x'][:-1] * (np.exp(dU) - 1 - dU)
        lhs = d['lam'] * beta.mean(); rhs = -(Pp[2] * w).mean()
        ident_err = max(ident_err, abs(lhs - rhs) / (abs(lhs) + 1e-30))
        signs.append(np.sign(rhs))
    print(f"mu={mu} rho={rho}: beta_s(t) range over all s,t: [{worst_neg:.4g}, {worst_pos:.4g}]  "
          f"identity rel.err {ident_err:.1e}  sign(-<Psi3 x c>) all {set(signs)}")
