#!/usr/bin/env python3
# --- v3.0 bundle copy of research file `lanes/planar/verify/vlog.py`.
# --- Changes from the original: none.
"""Referee library, log-variable formulation (keeps RELATIVE accuracy of y, lambda-p, eta, zeta through
the deep troughs at large rho, where x = k e^U underflows and the original-variable Riccati integration
falls below any absolute tolerance).  Independent of the lane's scripts and of whlib.

state s = (U, w, c, a, b, I):  p = lam - e^w,  y = e^c,  eta = e^a,  zeta = e^b,  I' = eta - mu
   U' = lam - e^w
   w' = e^{c-w} - mu                      (q = lam - p,  q' = y - mu q)
   c' = k e^{U-c} - 1                     (y' = x - y)
   a' = e^{b-a} + (mu-1) - e^a            (eta' = zeta + (mu-1) eta - eta^2)
   b' = k e^{U-b} - 1 - e^a               (zeta' = x - (1+eta) zeta)
"""
import math
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

RT = 1e-13
AT = np.array([1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-13])


def ex(z):
    return math.exp(z) if z < 700.0 else math.exp(700.0)


def eq_data(mu, k):
    f = lambda g: g * (g - 1.0) * (g - mu) - k
    lo = max(1.0, mu)
    hi = lo + 1.0
    while f(hi) <= 0:
        hi *= 2.0
    gam = brentq(f, lo, hi, xtol=1e-15, rtol=1e-15, maxiter=500)
    eta = gam - 1.0
    zeta = eta * (eta + 1.0 - mu)
    re = (eta - mu) / 2.0
    lam = complex(re, math.sqrt(zeta - re * re))
    v = np.array([1.0, lam, -(lam * lam + mu * lam)])
    return dict(gam=gam, eta=eta, zeta=zeta, lam=lam, v=v)


def log_rhs(mu, k):
    lam = k / mu

    def f(t, s):
        U, w, c, a, b, I = s
        ea = ex(a)
        return [lam - ex(w), ex(c - w) - mu, k * ex(U - c) - 1.0, ex(b - a) + (mu - 1.0) - ea,
                k * ex(U - b) - 1.0 - ea, ea - mu]
    return f


def to_log(mu, k, U, p, y, eta, zeta, I=0.0):
    return np.array([U, math.log(k / mu - p), math.log(y), math.log(eta), math.log(zeta), I])


def spiral_start(mu, k, r0, theta):
    ed = eq_data(mu, k)
    v = ed['v'] / np.linalg.norm(ed['v'])
    d = np.real(np.exp(1j * theta) * v)
    return to_log(mu, k, r0 * d[0], r0 * d[1], k + r0 * d[2], ed['eta'], ed['zeta'])


def make_crest(mu, k):
    ll = math.log(k / mu)

    def g(t, s):
        return ll - s[1]            # sign of p
    g.direction = -1
    return g


def run_log(mu, k, s0, t0, t1, t_eval=None, events=None, rtol=RT, atol=AT):
    return solve_ivp(log_rhs(mu, k), [t0, t1], s0, method='DOP853', rtol=rtol, atol=atol,
                     t_eval=t_eval, events=events)


def spiral_until_converged(mu, k, r0, theta, tol=1e-9, max_turns=6000, chunk_turns=40, rtol=RT, atol=AT):
    """crest rows: (t, U, c, a, b, I)"""
    ed = eq_data(mu, k)
    per = 2 * math.pi / ed['lam'].imag
    s = spiral_start(mu, k, r0, theta)
    t = 0.0
    crests = []
    ev = make_crest(mu, k)
    while len(crests) < max_turns:
        t1 = t + chunk_turns * per
        sol = run_log(mu, k, s, t, t1, events=[ev], rtol=rtol, atol=atol)
        if sol.status < 0:
            raise RuntimeError(sol.message)
        for te, ye in zip(sol.t_events[0], sol.y_events[0]):
            crests.append([te, ye[0], ye[2], ye[3], ye[4], ye[5]])
        s = sol.y[:, -1]
        t = t1
        if len(crests) >= 3:
            c = np.array(crests[-3:])
            if abs(c[-1, 1] - c[-2, 1]) < tol * abs(c[-1, 1]) and abs(c[-2, 1] - c[-3, 1]) < tol * abs(c[-1, 1]):
                break
    return t, np.array(crests)


# ---------------- variational system in original variables (no Riccati involved) ----------------
def var_rhs(mu, k, m):
    def f(t, s):
        U, p, y = s[0], s[1], s[2]
        x = k * ex(U)
        out = [p, k - mu * p - y, x - y]
        for j in range(m):
            a, b, c = s[3 + 3 * j], s[4 + 3 * j], s[5 + 3 * j]
            out += [b, -mu * b - c, x * a - c]
        return out
    return f


def _crest3(t, s):
    return s[1]
_crest3.direction = -1
_crest3.terminal = True


def crest_map(mu, k, U0, y0, tmax=None):
    s0 = np.concatenate(([U0, 0.0, y0], np.eye(3).ravel()))
    f = var_rhs(mu, k, 3)
    lam = k / mu
    at = np.concatenate(([1e-14, 1e-14 * max(1.0, lam), 1e-14 * k], 1e-14 * np.ones(9)))
    if tmax is None:
        tmax = 4000.0 * max(1.0, 1.0 / math.sqrt(mu))
    toff = 1e-3 * min(1.0, 1.0 / max(1.0, abs(k - y0)))
    a = solve_ivp(f, [0, toff], s0, method='DOP853', rtol=1e-12, atol=at)
    b = solve_ivp(f, [toff, tmax], a.y[:, -1], method='DOP853', rtol=1e-12, atol=at, events=_crest3)
    if b.status != 1:
        raise RuntimeError('no return')
    tau = float(b.t_events[0][0])
    s = b.y_events[0][0]
    U1, p1, y1 = s[:3]
    Phi = s[3:].reshape(3, 3).T
    fz = np.array([p1, k - mu * p1 - y1, k * ex(U1) - y1])
    M = np.zeros((2, 2))
    for j, e in enumerate((np.array([1.0, 0, 0]), np.array([0, 0, 1.0]))):
        d = Phi @ e
        d = d - fz * (d[1] / fz[1])
        M[:, j] = [d[0], d[2]]
    return U1, y1, tau, M, Phi


def newton_cycle(mu, k, U0, y0, iters=40, tol=1e-12):
    res = np.inf
    for _ in range(iters):
        U1, y1, T, M, Phi = crest_map(mu, k, U0, y0)
        r = np.array([U1 - U0, y1 - y0])
        res = float(np.hypot(r[0], r[1] / k))
        if res < tol:
            break
        d = np.linalg.solve(M - np.eye(2), -r)
        U0, y0 = U0 + d[0], y0 + d[1]
    U1, y1, T, M, Phi = crest_map(mu, k, U0, y0)
    return dict(U0=float(U0), y0=float(y0), T=T, M=M, Phi=Phi, res=res)


def cycle_riccati(mu, k, cyc, n=2401, max_periods=400, tol=1e-11, a0=None, b0=None):
    ed = eq_data(mu, k)
    T = cyc['T']
    a = math.log(ed['eta']) if a0 is None else a0
    b = math.log(ed['zeta']) if b0 is None else b0
    base = to_log(mu, k, cyc['U0'], 0.0, cyc['y0'], 1.0, 1.0)
    last = None
    hist = []
    for it in range(max_periods):
        s = base.copy(); s[3] = a; s[4] = b
        sol = run_log(mu, k, s, 0.0, T)
        a, b = sol.y[3, -1], sol.y[4, -1]
        if last is not None:
            hist.append(max(abs(a - last[0]), abs(b - last[1])))
            if hist[-1] < tol:
                break
        last = (a, b)
    s = base.copy(); s[3] = a; s[4] = b
    t = np.linspace(0, T, n)
    sol = run_log(mu, k, s, 0.0, T, t_eval=t)
    return dict(t=t, S=sol.y.copy(), T=T, periods_used=it + 1, hist=hist,
                defect=float(np.max(np.abs(sol.y[:5, -1] - sol.y[:5, 0]))))


def qr_exponents(mu, k, cyc, nsub=300, periods=3):
    T = cyc['T']
    f = var_rhs(mu, k, 2)
    lam = k / mu
    at = np.concatenate(([1e-14, 1e-14 * max(1.0, lam), 1e-14 * k], 1e-14 * np.ones(6)))
    Q = np.linalg.qr(np.random.default_rng(1).standard_normal((3, 2)))[0]
    out = []
    for per in range(periods):
        z = np.array([cyc['U0'], 0.0, cyc['y0']])
        acc = np.zeros(2)
        ts = np.linspace(0, T, nsub + 1)
        for i in range(nsub):
            s0 = np.concatenate((z, Q.T.ravel()))
            sol = solve_ivp(f, [ts[i], ts[i + 1]], s0, method='DOP853', rtol=1e-12, atol=at)
            z = sol.y[:3, -1]
            V = sol.y[3:, -1].reshape(2, 3).T
            Q, R = np.linalg.qr(V)
            acc += np.log(np.abs(np.diag(R)))
        out.append(acc.copy())
    return out


def nsum(terms):
    """terms: list of (sign array, log-magnitude array).  Returns sum / (largest |term|), elementwise,
    computed without underflow."""
    L = np.stack([t[1] for t in terms])
    S = np.stack([np.broadcast_to(t[0], L.shape[1:]) for t in terms])
    M = np.max(L, axis=0)
    with np.errstate(invalid='ignore', over='ignore'):
        R = np.sum(S * np.exp(L - M), axis=0)
    return R, M
