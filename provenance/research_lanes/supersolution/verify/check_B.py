#!/usr/bin/env python3
"""Referee check 3 (Proposition B): strong-stable Floquet ratios by backward integration in the variables
   delta = beta - 1,  nu = kappa - (1-mu) beta,  L = log a      (beta = b/a, kappa = c/a)
   delta' = delta(1+delta) - nu,   nu' = nu(delta + 1 - mu) - x,   L' = -(1+delta)
(these are the lane's v/a and (b-a)/a; integrating them directly keeps relative accuracy when they are tiny).
Tests: lambda_c = <beta> - A;  chi(gamma) = <x a^>/<a^>;  lambda_c = <a^ x(lam-A-2p)>/<c^P + b^y + a^x>;
       lambda_c = <a^(x-kH)>/<c^ + mu b^ + mu a^>;  k = <chi(beta)>;  beta > 1, kappa > (1-mu)beta."""
import sys, time
import numpy as np
from scipy.integrate import solve_ivp
from vprofile import profile, load_cycles
C = load_cycles()
keys = [(mu, rho) for mu in (0.02, 0.3, 1.0) for rho in (1.01, 1.5, 4.0, 12.0)]
print(f"{'mu':>5} {'rho':>5} {'<eta>-mu':>13} {'<beta>-A':>13} {'F1/lam':>10} {'F2/lam':>10} {'chi(g)':>12} {'<xa>/<a>':>12} {'kH':>8} {'<chi(beta)>/k':>13} {'min(beta-1)':>11} {'min(kap-(1-mu)beta)':>19} {'amax/amin':>9} {'perdef':>8} {'s':>5}")
for key in keys:
    mu, rho = key; cyc = C[key]; k = cyc['k']; lam = k/mu; A = 1+mu; T = cyc['T']
    t0 = time.time()
    R = profile(mu, k, cyc)
    fwd = R['sol']; gE = R['ed']['gammaE']
    def st(t):
        U, P, y = fwd.sol(t)[:3]
        return U, P, y
    def f(t, s):
        d, n = s
        U, P, y = st(t)
        x = k*np.exp(U)
        return [d*(1+d) - n, n*(d + 1 - mu) - x]
    s = np.array([gE - 1, gE*(gE - 1)])
    hist = []
    for _ in range(5):
        r = solve_ivp(f, [T, 0], s, method='DOP853', rtol=1e-12, atol=1e-250)
        new = r.y[:, -1]; hist.append(float(np.max(np.abs(new/s - 1)))); s = new
    # gamma
    def f2(t, s):
        d, n, Ib = s
        U, P, y = st(t); x = k*np.exp(U)
        return [d*(1+d) - n, n*(d + 1 - mu) - x, 1 + d]
    r = solve_ivp(f2, [T, 0], [s[0], s[1], 0.0], method='DOP853', rtol=1e-12, atol=[1e-250, 1e-250, 1e-13])
    gamma = -r.y[2, -1]/T
    def f3(t, s):
        d, n, L = s[:3]
        U, P, y = st(t); x = k*np.exp(U); p = lam - P
        be = 1 + d; ka = n + (1-mu)*be
        ah = np.exp(L + gamma*(t - T))
        return [d*(1+d) - n, n*(d + 1 - mu) - x, -be,
                ah*x, ah, ah*x*(lam - A - 2*p), ah*(ka*P + be*y + x), ah*(x - mu*A), ah*(ka + mu*be + mu), be*(be-1)*(be-mu)]
    r = solve_ivp(f3, [T, 0], [s[0], s[1], 0.0] + [0.0]*7, method='DOP853', rtol=1e-12,
                  atol=[1e-250, 1e-250, 1e-13] + [1e-14]*7, dense_output=True)
    J = -r.y[3:, -1]/T
    tt = np.unique(np.concatenate([r.t] + [r.t[:-1] + fr*np.diff(r.t) for fr in (0.25, 0.5, 0.75)]))
    Y = r.sol(tt)
    d, n, L = Y[:3]
    ah = np.exp(L + gamma*(tt - T))
    perdef = float(np.max(np.abs(r.y[:2, -1]/r.y[:2, 0] - 1)))
    lam_eta = R['lam_eta']
    F1 = J[2]/J[3]; F2 = J[4]/J[5]
    chi = gamma*(gamma-1)*(gamma-mu)
    print(f"{mu:5.2f} {rho:5.2f} {lam_eta:+13.6e} {gamma-A:+13.6e} {F1/lam_eta:10.7f} {F2/lam_eta:10.7f} {chi:12.8f} {J[0]/J[1]:12.8f} {mu*A:8.5f} {J[6]/k:13.9f} {d.min():11.3e} {n.min():19.3e} {ah.max()/ah.min():9.2e} {perdef:8.1e} {time.time()-t0:5.1f}", flush=True)
