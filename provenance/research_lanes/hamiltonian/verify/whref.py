#!/usr/bin/env python3
"""Referee's own integrator / cycle finder for the Wilhelm-Heinrich system (no code shared with whlib or the lane).

    U' = p,  p' = k - mu p - y,  y' = k e^U - y            (t-time)

Crest section: p = 0, p' < 0.  The return map is computed in two legs (crest -> trough -> crest) so that no
step-off from the section is needed.  The variational equations are integrated alongside (12 states).
"""
import numpy as np
from scipy.integrate import solve_ivp

RT, AT = 1e-12, 1e-14


def rhs_var(mu, k):
    def f(t, s):
        U, p, y = s[0], s[1], s[2]
        x = k * np.exp(U)
        Phi = s[3:12].reshape(3, 3)
        J = np.array([[0.0, 1.0, 0.0], [0.0, -mu, -1.0], [x, 0.0, -1.0]])
        return np.concatenate(([p, k - mu * p - y, x - y], (J @ Phi).ravel()))
    return f


def rhs(mu, k):
    def f(t, s):
        U, p, y = s
        return [p, k - mu * p - y, k * np.exp(U) - y]
    return f


def _ev_up(t, s):
    return s[1]
_ev_up.terminal = True
_ev_up.direction = 1


def _ev_dn(t, s):
    return s[1]
_ev_dn.terminal = True
_ev_dn.direction = -1


def crest_map(mu, k, U0, y0, var=True, tmax=None):
    """One return crest -> crest.  Returns (U1, y1, T, M2x2 or None, Umin)."""
    if tmax is None:
        tmax = 2000.0 / min(mu, 1.0) + 2000.0
    if var:
        f = rhs_var(mu, k)
        s0 = np.concatenate(([U0, 0.0, y0], np.eye(3).ravel()))
    else:
        f = rhs(mu, k)
        s0 = np.array([U0, 0.0, y0])
    a = solve_ivp(f, [0.0, tmax], s0, method="DOP853", rtol=RT, atol=AT, events=_ev_up)
    if a.status != 1:
        raise RuntimeError("no trough")
    tt = a.t_events[0][0]
    st = a.y_events[0][0]
    b = solve_ivp(f, [tt, tt + tmax], st, method="DOP853", rtol=RT, atol=AT, events=_ev_dn)
    if b.status != 1:
        raise RuntimeError("no crest")
    T = b.t_events[0][0]
    s = b.y_events[0][0]
    U1, p1, y1 = s[:3]
    if not var:
        return U1, y1, T, None, st[0]
    Phi = s[3:12].reshape(3, 3)
    fz = np.array([p1, k - mu * p1 - y1, k * np.exp(U1) - y1])
    M = np.zeros((2, 2))
    for j, v in enumerate((np.array([1.0, 0.0, 0.0]), np.array([0.0, 0.0, 1.0]))):
        d = Phi @ v
        d = d - fz * (d[1] / fz[1])
        M[0, j], M[1, j] = d[0], d[2]
    return U1, y1, T, M, st[0]


def cycle(mu, rho, guess=None, warm=60, newton=40, tol=1e-12, verbose=False):
    """Find a cycle at (mu, rho).  guess = (U0, y0) on the crest section."""
    k = rho * mu * (1 + mu)
    if guess is None:
        # generic start: a point with moderate amplitude on the approximate slow manifold
        a0 = 0.5
        guess = (a0, k + mu * rho * np.expm1(a0))
    U0, y0 = guess
    for i in range(warm):
        U1, y1, T, _, _ = crest_map(mu, k, U0, y0, var=False)
        d = abs(U1 - U0) + abs(y1 - y0) / k
        U0, y0 = U1, y1
        if verbose:
            print("  warm", i, U0, y0, d)
        if d < 1e-13 * (1 + abs(U0)):
            break
    res = np.inf
    for i in range(newton):
        U1, y1, T, M, Umin = crest_map(mu, k, U0, y0)
        r = np.array([U1 - U0, y1 - y0])
        res = abs(r[0]) / (1 + abs(U0)) + abs(r[1]) / k
        if verbose:
            print("  newton", i, U0, y0, res)
        if res < tol:
            break
        try:
            d = np.linalg.solve(M - np.eye(2), -r)
        except np.linalg.LinAlgError:
            d = r
        # damp absurd steps
        sc = max(1.0, abs(d[0]) / 0.5)
        U0, y0 = U0 + d[0] / sc, y0 + d[1] / sc
    U1, y1, T, M, Umin = crest_map(mu, k, U0, y0)
    ev = np.linalg.eigvals(M)
    m_t = float(np.max(ev.real))
    return dict(mu=mu, rho=rho, k=k, U0=float(U0), y0=float(y0), T=float(T), M=M, m_t=m_t,
                logmt=float(np.log(m_t)) if m_t > 0 else float('nan'), res=float(res), Umin=float(Umin))


def averages(mu, k, c, extra=None):
    """One period from the crest with running integrals.  Returns dict of period-averages in tau-variables:
       P = p/sqrt(mu), Q = p'/mu, R = p''/mu^{3/2} (= Q')."""
    T = c["T"]
    sm = np.sqrt(mu)

    def f(t, s):
        U, p, y = s[0], s[1], s[2]
        E = np.exp(U)
        x = k * E
        dp = k - mu * p - y
        ddp = -mu * dp - (x - y)
        P, Q, R = p / sm, dp / mu, ddp / mu ** 1.5
        return [p, dp, x - y,
                E, P * P, Q * Q, E * P * P, P ** 3, Q * Q * P, E * P ** 3, R * R, P ** 4, Q * Q * P * P, Q ** 3, P, Q, R,
                p * p, p ** 4, dp * dp * p * p, dp ** 3, ddp * ddp, dp * dp]
    s0 = np.zeros(23)
    s0[0], s0[2] = c["U0"], c["y0"]
    sol = solve_ivp(f, [0, T], s0, method="DOP853", rtol=1e-13, atol=1e-16)
    v = sol.y[:, -1]
    names = ["E", "P2", "Q2", "EP2", "P3", "Q2P", "EP3", "R2", "P4", "Q2P2", "Q3", "P", "Q", "R",
             "p2", "p4", "dp2p2", "dp3", "ddp2", "dp2"]
    out = {n: v[3 + i] / T for i, n in enumerate(names)}
    out["closure"] = float(abs(v[0] - c["U0"]) + abs(v[1]) + abs(v[2] - c["y0"]) / k)
    out["nsteps"] = sol.t.size
    return out
