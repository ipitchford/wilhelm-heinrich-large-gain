#!/usr/bin/env python3
"""Referee's own double-precision cycle finder (independent of whlib.find_cycle and of the lane's gen_data.rhs).

Log variables:  U,  a = log(y/k),  b = log(1 - p/lambda)   (lambda = k/mu), so
    U' = lambda (1 - e^b),   a' = e^{U-a} - 1,   b' = mu (e^{a-b} - 1).
Crest section: b = 0 crossed upward (p goes from + to -).  Return map (U0,a0) -> (U1,a1); Newton with
central finite differences.  Nothing here is rigorous; it only seeds the high-precision Taylor integrator.
"""
import math
import numpy as np
from scipy.integrate import solve_ivp

RT, AT = 1e-12, 1e-14


def rhs_log(mu, k):
    lam = k / mu

    def ex(v):          # trial steps of the integrator may overshoot; clip instead of raising
        return math.exp(v) if v < 700.0 else math.exp(700.0)

    def f(t, s):
        U, a, b = s
        return [-lam * math.expm1(min(b, 700.0)), ex(U - a) - 1.0, mu * (ex(a - b) - 1.0)]
    return f


def _up(t, s):
    return s[2]
_up.direction = 1
_up.terminal = True


def ret(mu, k, U0, a0, tmax=2.0e4):
    """first return to the crest section; returns (U1, a1, T, Umin)"""
    f = rhs_log(mu, k)
    t_off = 1e-4
    s1 = solve_ivp(f, [0, t_off], [U0, a0, 0.0], method='DOP853', rtol=RT, atol=AT).y[:, -1]
    if s1[2] <= 0:
        raise RuntimeError('start point is not a crest (b does not increase)')
    sol = solve_ivp(f, [t_off, tmax], s1, method='DOP853', rtol=RT, atol=AT, events=_up)
    if sol.status != 1:
        raise RuntimeError('no return')
    se = sol.y_events[0][0]
    return se[0], se[1], float(sol.t_events[0][0]), float(sol.y[0].min())


def find(mu, k, U0=1.0, a0=None, warm=60, verbose=False):
    if a0 is None:
        a0 = 0.3 * U0
    T = None
    for i in range(warm):
        U1, a1, T, Um = ret(mu, k, U0, a0)
        d = abs(U1 - U0) + abs(a1 - a0)
        U0, a0 = U1, a1
        if verbose:
            print('  warm', i, U0, a0, T, d)
        if d < 2e-3 * max(abs(U0), 1e-3) and i >= 3:
            break
    res = np.inf
    J = None
    for it in range(40):
        U1, a1, T, Um = ret(mu, k, U0, a0)
        r = np.array([U1 - U0, a1 - a0])
        res = float(np.abs(r).max())
        if verbose:
            print('  newton', it, U0, a0, T, res)
        hU = 1e-6 * max(1.0, abs(U0)); ha = 1e-6 * max(1.0, abs(a0))
        Up, ap, _, _ = ret(mu, k, U0 + hU, a0); Um_, am_, _, _ = ret(mu, k, U0 - hU, a0)
        Vp, bp, _, _ = ret(mu, k, U0, a0 + ha); Vm, bm, _, _ = ret(mu, k, U0, a0 - ha)
        J = np.array([[(Up - Um_) / (2 * hU), (Vp - Vm) / (2 * ha)], [(ap - am_) / (2 * hU), (bp - bm) / (2 * ha)]])
        if res < 5e-12:
            break
        dlt = np.linalg.solve(J - np.eye(2), -r)
        # damping: never let U0 collapse to the equilibrium
        sc = 1.0
        while U0 + sc * dlt[0] < 0.3 * U0:
            sc *= 0.5
        U0, a0 = U0 + sc * dlt[0], a0 + sc * dlt[1]
    U1, a1, T, Umin = ret(mu, k, U0, a0)
    ev = np.linalg.eigvals(J)
    m_t = float(np.max(ev.real))
    return dict(mu=mu, k=k, U0=float(U0), a0=float(a0), y0=k * math.exp(a0), T=T, Umin=Umin, J=J, m_t=m_t,
                lam_M=(math.log(m_t) / T if m_t > 0 else float('nan')), res=res)


if __name__ == '__main__':
    import sys, time
    for mu in (0.02, 0.3, 1.0, 5.0):
        guess = 1.0
        for rho in (1.01, 1.5, 4.0, 12.0):
            k = rho * mu * (1 + mu)
            t0 = time.time()
            try:
                c = find(mu, k, U0=guess)
                print(f'mu={mu:<5} rho={rho:<5} k={k:.6g} U0={c["U0"]:.10f} y0/k={math.exp(c["a0"]):.10f} T={c["T"]:.8f} Umin={c["Umin"]:.3f} '
                      f'm_t={c["m_t"]:.6e} lam_M={c["lam_M"]:+.6e} res={c["res"]:.1e} [{time.time()-t0:.1f}s]', flush=True)
                guess = max(c['U0'], 0.5)
            except Exception as ex:
                print(f'mu={mu} rho={rho}: FAILED {ex!r}', flush=True)
