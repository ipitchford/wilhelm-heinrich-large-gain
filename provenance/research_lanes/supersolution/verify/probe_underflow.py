#!/usr/bin/env python3
"""Probe: the only negative values seen (min W = -3e-196 in check_R, min(beta-1) = -3e-245 in check_B, both at mu=0.02, rho=12)
sit at the absolute-tolerance floor of the integrator.  Re-integrate the W equation and the (delta, nu) system for that cycle
with the floor moved (atol 1e-200 -> 1e-300 resp. 1e-250 -> 1e-320): a genuine sign change would not move with atol."""
import numpy as np
from scipy.integrate import solve_ivp
from vprofile import profile, load_cycles, G1
C = load_cycles()
mu, rho = 0.02, 12.0; cyc = C[(mu, rho)]; k = cyc['k']; lam = k/mu; T = cyc['T']
R = profile(mu, k, cyc)
fwd = R['sol']
for atol in (1e-200, 1e-250, 1e-300):
    def f(t, s):
        U, P, y, eta = fwd.sol(t)[:4]
        return [-(1 + eta)*s[0] + k*np.exp(U)*(lam - P)**2]
    w0 = [R['W'][0]]
    r = solve_ivp(f, [0, T], w0, method='DOP853', rtol=1e-11, atol=atol, dense_output=True)
    tt = np.unique(np.concatenate([r.t] + [r.t[:-1] + j/4*np.diff(r.t) for j in range(1, 4)]))
    W = r.sol(tt)[0]
    print(f"W equation, atol={atol:.0e}: min W = {W.min():+.3e}  (W true order in the trough ~ exp(-800) ~ 1e-347, below double range)")
gE = R['ed']['gammaE']
for atol in (1e-250, 1e-290, 1e-320):
    def f(t, s):
        d, n = s
        U = fwd.sol(t)[0]
        return [d*(1+d) - n, n*(d + 1 - mu) - k*np.exp(U)]
    s = np.array([gE - 1, gE*(gE - 1)])
    for _ in range(3):
        r = solve_ivp(f, [T, 0], s, method='DOP853', rtol=1e-12, atol=atol, dense_output=True); s = r.y[:, -1]
    tt = np.unique(np.concatenate([r.t] + [r.t[:-1] + j/4*np.diff(r.t) for j in range(1, 4)]))
    d, n = r.sol(tt)
    print(f"(delta,nu) backward, atol={atol:.0e}: min(beta-1) = {d.min():+.3e}, min(kappa-(1-mu)beta) = {n.min():+.3e}")
