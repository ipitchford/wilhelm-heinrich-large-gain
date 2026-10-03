#!/usr/bin/env python3
# --- v3.0 bundle copy of research file `lanes/largegain/lgcycle.py`.
# --- Changes from the original: none.
"""Large-gain cycle finder + event extraction (floating point; diagnostics only).
State (U, p, w) with w = y - k:  U' = p,  p' = -mu p - w,  w' = -w + k(exp(U)-1).
The cycle is super-attracting for large k, so we just integrate forward and take the last period
between successive down-crossings of U = 0."""
import numpy as np
METHOD = 'DOP853'
from scipy.integrate import solve_ivp

def rhs(mu, k):
    def f(t, s):
        U, p, w = s[:3]
        return [p, -mu*p - w, -w + k*np.expm1(min(U, 60.0))]
    return f

def rhs_full(mu, k):
    # state: U,p,w, log a, eta, zeta   (eta = b/a, zeta = c/a;  (log a)' = eta - mu)
    def f(t, s):
        U, p, w, la, eta, zeta = s
        x = k*np.exp(min(U, 60.0))
        return [p, -mu*p - w, -w + (x - k), eta - mu, zeta + (mu-1)*eta - eta*eta, x - (1+eta)*zeta]
    return f

def down(t, s): return s[0]
down.direction = -1

def find_cycle_lg(mu, k, nper=4, rtol=1e-10, atol=1e-12, tmax=None):
    lam = k/mu
    s0 = np.array([-5.0, 0.5*lam, -0.9*k])
    if tmax is None:
        tmax = nper*4.0*lam/mu + 200/np.sqrt(mu)
    sol = solve_ivp(rhs(mu, k), [0, tmax], s0, method=METHOD, rtol=rtol, atol=atol, events=down)
    te = sol.t_events[0]
    if len(te) < 3:
        raise RuntimeError(f'too few crossings {len(te)}')
    ye = sol.y_events[0]
    T = te[-1] - te[-2]
    Tprev = te[-2] - te[-3]
    return dict(mu=mu, k=k, T=T, Tprev=Tprev, s_start=ye[-2].copy(), s_end=ye[-1].copy())

def cycle_dense(mu, k, cyc, rtol=1e-11, atol=1e-13):
    """Integrate one period from the down-crossing with dense output, including Riccati (after warmup)."""
    T = cyc['T']
    s0 = np.concatenate([cyc['s_start'], [0.0, 1.0, 1.0]])
    s0[0] = 0.0
    f = rhs_full(mu, k)
    # warm up the Riccati variables over two periods
    for _ in range(2):
        sol = solve_ivp(f, [0, T], s0, method=METHOD, rtol=rtol, atol=atol)
        e = sol.y[:, -1]
        s0 = np.concatenate([cyc['s_start'], [0.0, e[4], e[5]]]); s0[0] = 0.0
    sol = solve_ivp(f, [0, T], s0, method=METHOD, rtol=rtol, atol=atol, dense_output=True)
    return sol

if __name__ == '__main__':
    for mu, k in [(1.0, 16.0), (1.0, 160.0), (1.0, 1000.0), (0.3, 15.6), (0.1, 100.0)]:
        c = find_cycle_lg(mu, k)
        sol = cycle_dense(mu, k, c)
        e = sol.y[:, -1]
        print(mu, k, 'T', c['T'], 'Tprev', c['Tprev'], 'closure', np.abs(c['s_end']-c['s_start']).max(),
              'lam_c*T', e[3], 'lam_c', e[3]/c['T'])
