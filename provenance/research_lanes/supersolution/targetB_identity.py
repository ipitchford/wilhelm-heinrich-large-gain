#!/usr/bin/env python3
"""Target B: verify  gamma(gamma-1)(gamma-mu) = <x ahat>/<ahat>  and the sign identity
   lambda_c * <nu.omegaY> = <ahat x (lambda - A - 2p)>   (nu = adjoint Perron of compound system, nu3 ~ ahat).
We check the *quantitative* version: lambda_c = <ahat x (lam-A-2p)> / <nuhat . omegaY> with
nuhat = (chat', bhat', ahat) up to the right combination; here we simply test sign and the ratio
against candidates built from (ahat,bhat,chat)."""
import numpy as np
from base import load, avg
D = load()
print(f"{'mu':>5} {'rho':>5} {'lam_c':>11} {'chi(g)':>10} {'<xa>/<a>':>10} {'kH':>7} {'S=<a x(lam-A-2p)>/<a>':>22} {'S/lam_c':>10} {'N1':>10} {'S/(lam N)':>10}")
for (mu, rho), P in sorted(D.items()):
    k = P['k']; A = 1 + mu; lam = k / mu
    a = np.exp(P['log_ahat']); x = P['x']; p = P['p']; y = P['y']
    g = P['gamma_beta']
    chi = g * (g - 1) * (g - mu)
    b = P['beta'] * a; c = P['kappa'] * a
    S = avg(a * x * (lam - A - 2 * p)) / avg(a)
    # adjoint Perron of compound: nu.w const for w' = Cw.  Guess nu = (c, b, a) * e^{-lam_c t}? test: d/dt(nu.omegaY) = nu3 x g3 ...
    # candidate normaliser N = <c*(lam-p) + b*y + a*x>/<a>
    N = avg(c * (lam - p) + b * y + a * x) / avg(a)
    print(f"{mu:5.2f} {rho:5.2f} {P['lam_M']:+11.4e} {chi:10.6f} {avg(x*a)/avg(a):10.6f} {mu*A:7.4f} {S:+22.6e} {S/P['lam_M']:10.4f} {N:10.4f} {S/(P['lam_M']*N):10.6f}")
