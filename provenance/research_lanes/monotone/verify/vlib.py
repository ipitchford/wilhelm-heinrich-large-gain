#!/usr/bin/env python3
"""Referee library (independent of the lane's mono_lib.py and of whlib's integrator).

Log-coordinates (U, V, Q) with V = ln y, Q = ln(lam - p), lam = k/mu:
    U' = lam - e^Q,   V' = k e^{U-V} - 1,   Q' = e^{V-Q} - mu.
These keep RELATIVE accuracy of y and of lam-p in astronomically deep troughs.
Crest section: p = 0, p' < 0  <=>  Q = ln(lam), Q increasing.
"""
import numpy as np
from scipy.integrate import solve_ivp

RT, AT = 1e-12, 1e-13


def kH(mu):
    return mu * (1 + mu)


def rhs(mu, k):
    lam = k / mu

    def f(t, w):
        U, V, Q = w
        return [lam - np.exp(Q), k * np.exp(U - V) - 1.0, np.exp(V - Q) - mu]
    return f


def crest_event(mu, k):
    L = np.log(k / mu)

    def g(t, w):
        return w[2] - L
    g.direction = 1
    g.terminal = True
    return g


def ret(mu, k, U0, y0, tmax=None, dense=False, method='DOP853'):
    """first return to the crest section from (U0, p=0, y0).  returns U1, y1, T (and sol)."""
    lam = k / mu
    f = rhs(mu, k)
    w0 = np.array([U0, np.log(y0), np.log(lam)])
    toff = 1e-4
    a = solve_ivp(f, [0, toff], w0, method=method, rtol=RT, atol=AT, dense_output=dense)
    if tmax is None:
        tmax = 4000.0 * max(1.0, 1.0 / mu)
    b = solve_ivp(f, [toff, tmax], a.y[:, -1], method=method, rtol=RT, atol=AT,
                  events=crest_event(mu, k), dense_output=dense)
    if b.status != 1:
        raise RuntimeError('no crest return')
    T = float(b.t_events[0][0])
    U1, V1, Q1 = b.y_events[0][0]
    if dense:
        def sol(t):
            t = np.atleast_1d(t)
            out = np.empty((3, t.size))
            m = t <= toff
            if m.any():
                out[:, m] = a.sol(t[m])
            if (~m).any():
                out[:, ~m] = b.sol(t[~m])
            return out
        return U1, float(np.exp(V1)), T, sol
    return U1, float(np.exp(V1)), T


def jac_fd(mu, k, U0, y0, h=1e-6):
    """2x2 derivative of the crest return map in (U, y) by central differences."""
    M = np.empty((2, 2))
    hU = h * max(1.0, abs(U0)); hy = h * y0
    a = ret(mu, k, U0 + hU, y0); b = ret(mu, k, U0 - hU, y0)
    M[:, 0] = [(a[0] - b[0]) / (2 * hU), (a[1] - b[1]) / (2 * hU)]
    a = ret(mu, k, U0, y0 + hy); b = ret(mu, k, U0, y0 - hy)
    M[:, 1] = [(a[0] - b[0]) / (2 * hy), (a[1] - b[1]) / (2 * hy)]
    return M


def eq_data(mu, k):
    roots = np.roots([1.0, 1.0 + mu, mu, k])
    gE = -min(r.real for r in roots if abs(r.imag) < 1e-10)
    J = np.array([[0, 1, 0], [0, -mu, -1], [k, 0, -1.0]])
    w, V = np.linalg.eig(J)
    i = int(np.argmax(w.imag))
    return dict(gE=gE, etaE=gE - 1, zetaE=(gE - 1) * (gE - mu), lam_u=w[i], vec_u=V[:, i])


def find_cycle(mu, k, verbose=False):
    """own cycle finder: slow-manifold bisection on the crest amplitude, then Newton (FD Jacobian)."""
    ed = eq_data(mu, k)
    # crest section points E + a*(1, zetaE) (direction inside the centre-unstable plane of E, p=0)
    def pt(a):
        return a, k + ed['zetaE'] * a

    def growth(a):
        U0, y0 = pt(a)
        U1, y1, _ = ret(mu, k, U0, y0)
        U2, y2, _ = ret(mu, k, U1, y1)
        return U2 - U1, (U2, y2)
    lo, hi = 1e-3, 1e-3
    g, _ = growth(lo)
    if g <= 0:
        raise RuntimeError('no outward growth near E')
    for _ in range(60):
        hi *= 2.0
        g, last = growth(hi)
        if g < 0:
            break
        lo = hi
    else:
        raise RuntimeError('no bracket')
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        g, last_mid = growth(mid)
        if g > 0:
            lo = mid
        else:
            hi = mid; last = last_mid
        if hi - lo < 1e-7 * hi:
            break
    U0, y0 = last
    res = np.inf
    for it in range(25):
        U1, y1, T = ret(mu, k, U0, y0)
        rvec = np.array([U1 - U0, y1 - y0])
        res = abs(rvec[0]) + abs(rvec[1]) / k
        if verbose:
            print('  newton', it, U0, y0, res)
        if res < 3e-11:
            break
        M = jac_fd(mu, k, U0, y0)
        d = np.linalg.solve(M - np.eye(2), -rvec)
        U0, y0 = U0 + d[0], y0 + d[1]
    U1, y1, T, sol = ret(mu, k, U0, y0, dense=True)
    return dict(mu=mu, k=k, U0=U0, y0=y0, T=T, sol=sol, res=res)


def state(mu, k, cyc, t):
    """all state-derived quantities at times t (mod T) on the cycle (crest at t=0)."""
    lam = k / mu
    T = cyc['T']
    tt = np.mod(np.atleast_1d(t), T)
    U, V, Q = cyc['sol'](tt)
    y = np.exp(V); q = np.exp(Q); x = k * np.exp(U)
    p = lam - q
    pp = mu * q - y          # p' = k - mu p - y
    yp = x - y
    ppp = -mu * pp - yp
    ypp = x * p - yp
    return dict(t=tt, U=U, p=p, y=y, x=x, q=q, pp=pp, yp=yp, ppp=ppp, ypp=ypp)


def strong_stable(mu, k, cyc, t_eval, nper=3):
    """r=b/a, s=c/a of the strong-stable Floquet vector e_s=(a,-b,-c) by BACKWARD integration of
       r' = r^2 - mu r - s,  s' = s(r-1) - x   in variables (r, ln s).  Returns r, s on t_eval (in [0,T])."""
    T = cyc['T']
    ed = eq_data(mu, k)

    def g(tau, w):      # tau = T - t (backward)
        U = cyc['sol'](np.mod(T - tau, T))[0, 0]
        r, ls = w
        s = np.exp(ls)
        return [-(r * r - mu * r - s), -((r - 1.0) - k * np.exp(U - ls))]
    w0 = np.array([ed['gE'], np.log(ed['gE'] * (ed['gE'] - mu))])
    for _ in range(nper):
        sb = solve_ivp(g, [0, T], w0, method='DOP853', rtol=RT, atol=AT)
        w1 = sb.y[:, -1].copy()
        defect = np.abs(w1 - w0).max()
        w0 = w1
    te = np.asarray(t_eval)
    order = np.argsort(T - te)
    sb = solve_ivp(g, [0, T], w0, method='DOP853', rtol=RT, atol=AT, t_eval=(T - te)[order], dense_output=True)
    r = np.empty(te.size); s = np.empty(te.size)
    r[order] = sb.y[0]; s[order] = np.exp(sb.y[1])
    return r, s, defect, sb.sol


def adjoint_ratios(mu, k, cyc, t_eval, nper=4):
    """eta, zeta (forward Riccati; forward-stable), in variables (eta, ln zeta)."""
    T = cyc['T']
    ed = eq_data(mu, k)

    def f(t, w):
        U = cyc['sol'](np.mod(t, T))[0, 0]
        eta, lz = w
        z = np.exp(lz)
        return [z + (mu - 1) * eta - eta * eta, k * np.exp(U - lz) - (1 + eta)]
    w0 = np.array([ed['etaE'], np.log(ed['zetaE'])])
    for _ in range(nper):
        sf = solve_ivp(f, [0, T], w0, method='DOP853', rtol=RT, atol=AT)
        w1 = sf.y[:, -1].copy(); defect = np.abs(w1 - w0).max(); w0 = w1
    sf = solve_ivp(f, [0, T], w0, method='DOP853', rtol=RT, atol=AT, t_eval=np.asarray(t_eval))
    return sf.y[0], np.exp(sf.y[1]), defect
