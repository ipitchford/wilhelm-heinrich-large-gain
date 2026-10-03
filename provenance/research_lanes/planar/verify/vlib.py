#!/usr/bin/env python3
"""Referee library for lane `planar` (written independently; shares no code with the lane's
scripts and does not import whlib).  Ordinary floating point; nothing here is a certificate.

System:  U' = p,  p' = k - mu p - y,  y' = k e^U - y,   Riccati: eta' = zeta + (mu-1) eta - eta^2,
zeta' = x - (1+eta) zeta.  Augmented state s = (U, p, y, eta, zeta, I),  I' = eta - mu.
"""
import math
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

RT = 1e-12


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
    lam = complex(re, math.sqrt(zeta - re * re))      # roots of z^2 + (mu-eta) z + zeta
    v = np.array([1.0, lam, -(lam * lam + mu * lam)])  # (dU, dp, dy)
    return dict(gam=gam, eta=eta, zeta=zeta, lam=lam, v=v)


def aug_rhs(mu, k):
    ex = math.exp

    def f(t, s):
        U, p, y, eta, zeta, I = s
        x = k * ex(U)
        return [p, k - mu * p - y, x - y, zeta + (mu - 1.0) * eta - eta * eta, x - (1.0 + eta) * zeta, eta - mu]
    return f


def atol_vec(mu, k):
    lam = k / mu
    return np.array([1e-14, 1e-14 * max(1.0, lam), 1e-14 * k, 1e-14, 1e-14 * k, 1e-13])


def _crest(t, s):
    return s[1]
_crest.direction = -1


def _trough(t, s):
    return s[1]
_trough.direction = +1


def spiral_start(mu, k, r0, theta):
    ed = eq_data(mu, k)
    v = ed['v'] / np.linalg.norm(ed['v'])
    d = np.real(np.exp(1j * theta) * v)
    return np.array([r0 * d[0], r0 * d[1], k + r0 * d[2], ed['eta'], ed['zeta'], 0.0])


def run_aug(mu, k, s0, t0, t1, t_eval=None, events=None, dense=False, method='DOP853', rtol=RT):
    return solve_ivp(aug_rhs(mu, k), [t0, t1], s0, method=method, rtol=rtol, atol=atol_vec(mu, k),
                     t_eval=t_eval, events=events, dense_output=dense)


def spiral_until_converged(mu, k, r0, theta, tol=1e-9, max_turns=6000, chunk_turns=40):
    """Pass 1: follow the W^u(E) spiral, collecting crest events, until the crest increments are < tol
    (relative).  Returns (t_conv, crest array rows = (t, U, y, eta, zeta, I))."""
    ed = eq_data(mu, k)
    per = 2 * math.pi / ed['lam'].imag
    s = spiral_start(mu, k, r0, theta)
    t = 0.0
    crests = []
    turns = 0
    while turns < max_turns:
        t1 = t + chunk_turns * per
        sol = run_aug(mu, k, s, t, t1, events=[_crest])
        for te, ye in zip(sol.t_events[0], sol.y_events[0]):
            crests.append([te, ye[0], ye[2], ye[3], ye[4], ye[5]])
        s = sol.y[:, -1]
        t = t1
        turns = len(crests)
        if len(crests) >= 3:
            c = np.array(crests)
            if abs(c[-1, 1] - c[-2, 1]) < tol * abs(c[-1, 1]) and abs(c[-2, 1] - c[-3, 1]) < tol * abs(c[-1, 1]):
                break
    return t, np.array(crests)


def var_rhs(mu, k, m):
    """state (3) + m tangent vectors (3 each)."""
    ex = math.exp

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


def crest_map(mu, k, U0, y0, tmax=None, deriv=True):
    """First return to {p=0, p'<0} from (U0,0,y0).  Returns U1,y1,tau,M (M = 2x2 section derivative in (U,y))."""
    m = 3 if deriv else 0
    s0 = np.concatenate(([U0, 0.0, y0], np.eye(3).T.ravel() if deriv else []))
    f = var_rhs(mu, k, m)
    at = np.concatenate((atol_vec(mu, k)[:3], 1e-13 * np.ones(3 * m)))
    if tmax is None:
        tmax = 2000.0 * max(1.0, 1.0 / math.sqrt(mu))
    # leave the section: p' = k - y0 < 0 at a crest, integrate a short time without events
    toff = 1e-3 * min(1.0, 1.0 / max(1.0, abs(k - y0)))
    a = solve_ivp(f, [0, toff], s0, method='DOP853', rtol=RT, atol=at)
    b = solve_ivp(f, [toff, tmax], a.y[:, -1], method='DOP853', rtol=RT, atol=at, events=_crest3)
    if b.status != 1:
        raise RuntimeError('no return')
    tau = float(b.t_events[0][0])
    s = b.y_events[0][0]
    U1, p1, y1 = s[:3]
    if not deriv:
        return U1, y1, tau, None
    Phi = s[3:].reshape(3, 3).T          # columns = tangent vectors
    fz = np.array([p1, k - mu * p1 - y1, k * math.exp(U1) - y1])
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


def cycle_riccati(mu, k, cyc, n=2001, max_periods=60, tol=1e-12):
    """Periodic Riccati pair on the cycle by forward (pull-back) iteration from the equilibrium values.
    Returns dict with arrays over one period starting at the crest."""
    ed = eq_data(mu, k)
    T = cyc['T']
    s = np.array([cyc['U0'], 0.0, cyc['y0'], ed['eta'], ed['zeta'], 0.0])
    last = None
    for it in range(max_periods):
        sol = run_aug(mu, k, s, 0.0, T)
        e = sol.y[3:5, -1].copy()
        s = np.array([cyc['U0'], 0.0, cyc['y0'], e[0], e[1], 0.0])
        if last is not None and np.max(np.abs(e - last) / np.abs(e)) < tol:
            break
        last = e
    t = np.linspace(0, T, n)
    sol = run_aug(mu, k, s, 0.0, T, t_eval=t)
    U, p, y, eta, zeta, I = sol.y
    return dict(t=t, U=U, p=p, y=y, eta=eta, zeta=zeta, I=I, T=T, periods_used=it + 1,
                defect=float(np.max(np.abs(sol.y[:5, -1] - sol.y[:5, 0]) / (np.abs(sol.y[:5, 0]) + 1e-300))))


def qr_exponents(mu, k, cyc, nsub=400, periods=3):
    """Independent of the Riccati equation: top-two Lyapunov sum over one period by discrete QR of the
    variational flow (2 tangent vectors).  Returns (sum of log diag R over last period) = lambda_c * T."""
    T = cyc['T']
    f = var_rhs(mu, k, 2)
    at = np.concatenate((atol_vec(mu, k)[:3], 1e-13 * np.ones(6)))
    Q = np.linalg.qr(np.random.default_rng(1).standard_normal((3, 2)))[0]
    out = []
    for per in range(periods):
        z = np.array([cyc['U0'], 0.0, cyc['y0']])
        acc = np.zeros(2)
        ts = np.linspace(0, T, nsub + 1)
        for i in range(nsub):
            s0 = np.concatenate((z, Q.T.ravel()))
            sol = solve_ivp(f, [ts[i], ts[i + 1]], s0, method='DOP853', rtol=RT, atol=at)
            z = sol.y[:3, -1]
            V = sol.y[3:, -1].reshape(2, 3).T
            Q, R = np.linalg.qr(V)
            acc += np.log(np.abs(np.diag(R)))
        out.append(acc.copy())
    return out
