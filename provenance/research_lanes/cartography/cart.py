#!/usr/bin/env python3
"""Cartography lane: compiled (numba) DOP853 toolkit for the Wilhelm--Heinrich oscillator.

Coordinates as in whlib / paper v2.0 eq. (1.5):
    U' = p,  p' = k - mu p - y,  y' = k e^U - y,     x = k e^U,  E = (0,0,k),  k_H = mu(1+mu).
Adjoint ratios (pull-back attracting):  eta' = zeta + (mu-1) eta - eta^2,  zeta' = x - (1+eta) zeta.
Everything is ordinary floating point.  Nothing here is an interval certificate.

Integrator: explicit DOP853 (coefficients imported from scipy), per-component rtol/atol,
events located to round-off by Newton on the step length (full RK step, no interpolant),
integration restarted at every event so return maps are deterministic functions of the section point.
"""
from __future__ import annotations
import numpy as np
from numba import njit
from scipy.integrate._ivp import dop853_coefficients as _dc

_A = np.ascontiguousarray(_dc.A[:12, :12])
_B = np.ascontiguousarray(_dc.B)
_E3 = np.ascontiguousarray(_dc.E3)
_E5 = np.ascontiguousarray(_dc.E5)

# ----------------------------------------------------------------------------- right-hand sides
# mode 0: (U,p,y)
# mode 1: (U,p,y,eta,zeta,K,W, I_eta, I_U, I_q2, I_q3, I_x)          q = eta - mu
# mode 2: (U,p,y, Phi[9] row-major, dz/dk[3])
# mode 3: (U,p,y,eta,zeta,K,W, dU,dp,dy,deta,dzeta, I_eta)            one transversal variation
NDIM = (3, 12, 15, 13)


@njit(cache=True)
def g2(U):
    """U e^U - e^U + 1 without cancellation near U = 0."""
    if abs(U) < 0.5:
        term = U * U / 2.0          # n = 2: U^2 (n-1)/n!
        acc = term
        un = U * U
        fact = 2.0
        for n in range(3, 30):
            un *= U
            fact *= n
            acc += un * (n - 1) / fact
        return acc
    return U * np.exp(U) - np.exp(U) + 1.0


@njit(cache=True)
def rhs(mode, s, par, out):
    mu = par[0]
    k = par[1]
    U = s[0]
    p = s[1]
    y = s[2]
    x = k * np.exp(U)
    out[0] = p
    out[1] = k - mu * p - y
    out[2] = x - y
    if mode == 0:
        return
    if mode == 1 or mode == 3:
        eta = s[3]
        zeta = s[4]
        out[3] = zeta + (mu - 1.0) * eta - eta * eta
        out[4] = x - (1.0 + eta) * zeta
        out[5] = -(1.0 + eta) * s[5] + k * g2(U)
        out[6] = -(1.0 + eta) * s[6] + x * p * p
        if mode == 1:
            q = eta - mu
            out[7] = eta
            out[8] = U
            out[9] = q * q
            out[10] = q * q * q
            out[11] = x
        else:
            dU = s[7]
            dp = s[8]
            dy = s[9]
            de = s[10]
            dz = s[11]
            out[7] = dp
            out[8] = -mu * dp - dy
            out[9] = x * dU - dy
            out[10] = dz + (mu - 1.0) * de - 2.0 * eta * de
            out[11] = x * dU - (1.0 + eta) * dz - zeta * de
            out[12] = eta
        return
    if mode == 2:
        for c in range(3):
            a = s[3 + c]        # Phi[0,c]
            b = s[6 + c]        # Phi[1,c]
            d = s[9 + c]        # Phi[2,c]
            out[3 + c] = b
            out[6 + c] = -mu * b - d
            out[9 + c] = x * a - d
        a = s[12]
        b = s[13]
        d = s[14]
        out[12] = b
        out[13] = -mu * b - d + 1.0
        out[14] = x * a - d + x / k
        return


@njit(cache=True)
def _step(mode, par, y, h, K, ytmp, ynew):
    """One DOP853 step of size h from y.  K[0] must already hold f(y).  K[12] = f(ynew) on exit."""
    n = y.shape[0]
    for s in range(1, 12):
        for i in range(n):
            acc = 0.0
            for j in range(s):
                acc += _A[s, j] * K[j, i]
            ytmp[i] = y[i] + h * acc
        rhs(mode, ytmp, par, K[s])
    for i in range(n):
        acc = 0.0
        for j in range(12):
            acc += _B[j] * K[j, i]
        ynew[i] = y[i] + h * acc
    rhs(mode, ynew, par, K[12])


@njit(cache=True)
def _errnorm(K, h, y, ynew, rtol, atol):
    n = y.shape[0]
    e5 = 0.0
    e3 = 0.0
    for i in range(n):
        sc = atol[i] + rtol[i] * max(abs(y[i]), abs(ynew[i]))
        a5 = 0.0
        a3 = 0.0
        for j in range(13):
            a5 += K[j, i] * _E5[j]
            a3 += K[j, i] * _E3[j]
        a5 /= sc
        a3 /= sc
        e5 += a5 * a5
        e3 += a3 * a3
    if e5 == 0.0 and e3 == 0.0:
        return 0.0
    denom = e5 + 0.01 * e3
    return abs(h) * e5 / np.sqrt(denom * n)


@njit(cache=True)
def integrate(mode, par, y, t0, tmax, rtol, atol, ev_idx, ev_dir, max_ev, ev_t, ev_y,
              rec_t, rec_y, hmax, maxsteps):
    """Integrate in place.  Returns (status, t, nsteps, nev, nrec).
    status 0: reached tmax;  1: stopped after max_ev events;  -1: step size underflow / maxsteps;  -2: non-finite.
    Event: component ev_idx crossing zero; ev_dir -1 (downward), +1 (upward), 0 (both); ev_idx<0 disables.
    ev_dir -2 / +2 / 3: same as -1 / +1 / 0 but a downward crossing is accepted only at U>0 and an upward one
    only at U<0 (on the attractor crests have U>0, troughs U<0; this rejects spurious re-crossings caused by
    stability-limited wiggles of p around its quasi-static value in stiff troughs).
    The event component is set to exactly 0 in the stored / restarted state."""
    n = y.shape[0]
    K = np.zeros((13, n))
    ytmp = np.zeros(n)
    ynew = np.zeros(n)
    yev = np.zeros(n)
    t = t0
    nsteps = 0
    nev = 0
    nrec = 0
    nrec_max = rec_t.shape[0]
    rhs(mode, y, par, K[0])
    # initial step
    d0 = 0.0
    d1 = 0.0
    for i in range(n):
        sc = atol[i] + rtol[i] * abs(y[i])
        d0 += (y[i] / sc) ** 2
        d1 += (K[0, i] / sc) ** 2
    d0 = np.sqrt(d0 / n)
    d1 = np.sqrt(d1 / n)
    if d0 < 1e-5 or d1 < 1e-5 or not np.isfinite(d0) or not np.isfinite(d1):
        h = 1e-6
    else:
        h = 0.01 * d0 / d1
    if not (h > 1e-9):
        h = 1e-9
    h = min(h, hmax, tmax - t0)
    if nrec_max > 0:
        rec_t[0] = t
        for i in range(n):
            rec_y[0, i] = y[i]
        nrec = 1
    rejected = False
    while t < tmax:
        if nsteps >= maxsteps:
            return -1, t, nsteps, nev, nrec
        if h < 1e-14 * max(1.0, abs(t)):
            return -1, t, nsteps, nev, nrec
        last = False
        if t + h >= tmax:
            h = tmax - t
            last = True
        _step(mode, par, y, h, K, ytmp, ynew)
        err = _errnorm(K, h, y, ynew, rtol, atol)
        if not np.isfinite(err):
            h *= 0.2
            rejected = True
            if h < 1e-300:
                return -2, t, nsteps, nev, nrec
            continue
        if err >= 1.0:
            h *= max(0.2, 0.9 * err ** (-0.125))
            rejected = True
            continue
        # accepted
        nsteps += 1
        if err == 0.0:
            fac = 10.0
        else:
            fac = min(10.0, 0.9 * err ** (-0.125))
        if rejected:
            fac = min(1.0, fac)
        rejected = False
        hit = False
        if ev_idx >= 0:
            g0 = y[ev_idx]
            g1 = ynew[ev_idx]
            down_ok = ev_dir == -1 or ev_dir == 0 or ev_dir == -2 or ev_dir == 3
            up_ok = ev_dir == 1 or ev_dir == 0 or ev_dir == 2 or ev_dir == 3
            filt = ev_dir <= -2 or ev_dir >= 2
            if down_ok and g0 > 0.0 and g1 <= 0.0:
                hit = (not filt) or ynew[0] > 0.0
            if up_ok and g0 < 0.0 and g1 >= 0.0:
                hit = (not filt) or ynew[0] < 0.0
        if hit:
            a = 0.0
            ga = g0
            b = h
            gb = g1
            hs = h * g0 / (g0 - g1)
            for it in range(60):
                _step(mode, par, y, hs, K, ytmp, yev)
                gs = yev[ev_idx]
                dg = K[12, ev_idx]
                if abs(gs) <= 4e-16 * abs(dg) * h or (b - a) <= 4e-16 * h:
                    break
                if (gs > 0.0) == (ga > 0.0):
                    a = hs
                    ga = gs
                else:
                    b = hs
                    gb = gs
                hn = hs - gs / dg
                if not (hn > a and hn < b):
                    hn = a + (b - a) * ga / (ga - gb)
                    if not (hn > a and hn < b):
                        hn = 0.5 * (a + b)
                hs = hn
            t = t + hs
            for i in range(n):
                y[i] = yev[i]
            y[ev_idx] = 0.0
            rhs(mode, y, par, K[0])
            ev_t[nev] = t
            for i in range(n):
                ev_y[nev, i] = y[i]
            nev += 1
            if nrec < nrec_max:
                rec_t[nrec] = t
                for i in range(n):
                    rec_y[nrec, i] = y[i]
                nrec += 1
            if nev >= max_ev:
                return 1, t, nsteps, nev, nrec
            h = min(h * fac, hmax)
            continue
        t = t + h
        for i in range(n):
            y[i] = ynew[i]
            K[0, i] = K[12, i]
        if nrec < nrec_max:
            rec_t[nrec] = t
            for i in range(n):
                rec_y[nrec, i] = y[i]
            nrec += 1
        if last:
            break
        h = min(h * fac, hmax)
    return 0, t, nsteps, nev, nrec


# ----------------------------------------------------------------------------- python wrappers
RT = 1e-12
_EMPTY_T = np.zeros(0)
_EMPTY_Y = np.zeros((0, 1))


def kH(mu):
    return mu * (1.0 + mu)


def tolerances(mode, mu, k, rt=RT):
    """Per-component tolerances.  U, p, y-k cross zero -> absolute floors scaled to the problem.
    eta, zeta, K, W are positive -> essentially pure relative control."""
    n = NDIM[mode]
    rtol = np.full(n, rt)
    atol = np.full(n, rt * 1e-1)
    sp = np.sqrt(mu)
    rtol[0] = 0.0          # U may be astronomically negative in troughs: absolute control only
    atol[0] = rt
    atol[1] = rt * sp
    atol[2] = rt * k
    if mode == 1 or mode == 3:
        atol[3] = 1e-200
        atol[4] = 1e-200
        atol[5] = 1e-200 if mode == 1 else rt * 1e-10 * k
        atol[6] = 1e-200 if mode == 1 else rt * 1e-10 * k * mu
        if mode == 1:
            atol[7:] = rt * 1e-1
            atol[11] = rt * k
        else:
            rtol[7:] = rt
            atol[7:] = 1e-200
            atol[12] = rt * 1e-1
    if mode == 2:
        rtol[3:] = rt
        atol[3:] = 1e-200
    return rtol, atol


def equilibrium_data(mu, k):
    roots = np.roots([1.0, 1.0 + mu, mu, k])
    gammaE = -min(r.real for r in roots if abs(r.imag) < 1e-9 * max(1.0, abs(r)))
    # polish the real root by Newton
    g = gammaE
    for _ in range(50):
        F = g * (g - 1) * (g - mu) - k
        dF = (g - 1) * (g - mu) + g * (g - mu) + g * (g - 1)
        g2 = g - F / dF
        if abs(g2 - g) < 1e-16 * abs(g):
            g = g2
            break
        g = g2
    gammaE = g
    etaE = gammaE - 1.0
    zetaE = etaE * (etaE + 1.0 - mu)
    alpha = 0.5 * (gammaE - (1.0 + mu))            # real part of unstable pair
    omega2 = k / gammaE - alpha * alpha            # |lam|^2 = k/gammaE
    omega = np.sqrt(max(omega2, 0.0))
    lam = complex(alpha, omega)
    # eigenvector of J = [[0,1,0],[0,-mu,-1],[k,0,-1]] : (1, lam, k/(lam+1))
    vec = np.array([1.0, lam, k / (lam + 1.0)])
    return dict(gammaE=gammaE, etaE=etaE, zetaE=zetaE, lam_u=lam, vec_u=vec, alpha=alpha, omega=omega)


def tmax_default(mu, k):
    return 1e7 * max(1.0, 1.0 / np.sqrt(mu))


def crest_return(mu, k, U0, y0, rt=RT, tmax=None):
    """First return to {p=0, p'<0}.  Returns (U1, y1, T) or raises."""
    par = np.array([mu, k])
    s = np.array([U0, 0.0, y0])
    rtol, atol = tolerances(0, mu, k, rt)
    ev_t = np.zeros(1)
    ev_y = np.zeros((1, 3))
    st, t, ns, nev, nrec = integrate(0, par, s, 0.0, tmax or tmax_default(mu, k), rtol, atol, 1, -2, 1,
                                     ev_t, ev_y, _EMPTY_T, _EMPTY_Y, np.inf, 50_000_000)
    if st != 1:
        raise RuntimeError(f"no crest return (status {st}) mu={mu} k={k} U0={U0} y0={y0}")
    return s[0], s[2], t


def crest_iter(mu, k, U0, y0, n, rt=RT, tmax=None, maxsteps=500_000_000):
    """n successive crest returns.  Returns arrays U[n+1], y[n+1], t[n+1] (index 0 = start)."""
    par = np.array([mu, k])
    s = np.array([U0, 0.0, y0])
    rtol, atol = tolerances(0, mu, k, rt)
    ev_t = np.zeros(n)
    ev_y = np.zeros((n, 3))
    st, t, ns, nev, nrec = integrate(0, par, s, 0.0, tmax or tmax_default(mu, k) * n, rtol, atol, 1, -2, n,
                                     ev_t, ev_y, _EMPTY_T, _EMPTY_Y, np.inf, maxsteps)
    U = np.concatenate(([U0], ev_y[:nev, 0]))
    Y = np.concatenate(([y0], ev_y[:nev, 2]))
    tt = np.concatenate(([0.0], ev_t[:nev]))
    return U, Y, tt, st


def crest_return_var(mu, k, U0, y0, rt=RT, tmax=None, maxsteps=50_000_000):
    """Return (U1, y1, T, M, dk) with M the 2x2 section derivative in (U,y) and
    dk = (dU1/dk, dy1/dk) at fixed start point (section-projected)."""
    par = np.array([mu, k])
    s = np.zeros(15)
    s[0] = U0
    s[2] = y0
    s[3] = 1.0
    s[7] = 1.0
    s[11] = 1.0
    rtol, atol = tolerances(2, mu, k, rt)
    ev_t = np.zeros(1)
    ev_y = np.zeros((1, 15))
    st, t, ns, nev, nrec = integrate(2, par, s, 0.0, tmax or tmax_default(mu, k), rtol, atol, 1, -2, 1,
                                     ev_t, ev_y, _EMPTY_T, _EMPTY_Y, np.inf, maxsteps)
    if st != 1:
        raise RuntimeError(f"no crest return (status {st}) mu={mu} k={k} U0={U0} y0={y0}")
    U1, p1, y1 = s[0], s[1], s[2]
    Phi = s[3:12].reshape(3, 3)
    fz = np.array([p1, k - mu * p1 - y1, k * np.exp(U1) - y1])
    cols = []
    for v in (np.array([1.0, 0, 0]), np.array([0, 0, 1.0])):
        d = Phi @ v
        d = d - fz * (d[1] / fz[1])
        cols.append([d[0], d[2]])
    M = np.array(cols).T
    zk = s[12:15]
    zk = zk - fz * (zk[1] / fz[1])
    return U1, y1, t, M, np.array([zk[0], zk[2]])


def multipliers(M, mu, T):
    tr = M[0, 0] + M[1, 1]
    det = M[0, 0] * M[1, 1] - M[0, 1] * M[1, 0]
    disc = tr * tr - 4 * det
    m_t = 0.5 * (tr + np.sqrt(max(disc, 0.0)))
    return m_t, np.exp(-(1 + mu) * T) / m_t if m_t > 0 else np.nan


def newton_fixed_k(mu, k, U0, y0, maxit=40, tol=1e-12, rt=RT, maxsteps=50_000_000):
    """Newton on the crest map at fixed k.  Returns dict or raises."""
    res = np.inf
    for it in range(maxit):
        U1, y1, T, M, _ = crest_return_var(mu, k, U0, y0, rt, maxsteps=maxsteps)
        r = np.array([U1 - U0, y1 - y0])
        prev_res = res
        res = float(np.hypot(r[0], r[1] / k))
        if res < tol or (res < 1e-9 and res > 0.3 * prev_res):
            break
        d = np.linalg.solve(M - np.eye(2), -r)
        # damp wild steps
        sc = max(abs(d[0]) / max(0.5, 0.5 * abs(U0)), abs(d[1]) / (0.5 * abs(y0)))
        if sc > 1:
            d = d / sc
        U0, y0 = U0 + d[0], y0 + d[1]
        if not (y0 > 0 and np.isfinite(U0)):
            raise RuntimeError("newton left the domain")
    U1, y1, T, M, _ = crest_return_var(mu, k, U0, y0, rt)
    m_t, m_s = multipliers(M, mu, T)
    return dict(mu=mu, k=k, rho=k / kH(mu), U0=float(U0), y0=float(y0), T=float(T), M=M, m_t=float(m_t),
                residual=float(np.hypot(U1 - U0, (y1 - y0) / k)), iters=it + 1)


def newton_fixed_h(mu, h, y0, k, maxit=40, tol=1e-12, rt=RT):
    """Amplitude-fixed shooting: unknowns (y0, k) with crest U0 = h.  Well conditioned through Hopf."""
    res = np.inf
    for it in range(maxit):
        U1, y1, T, M, dk = crest_return_var(mu, k, h, y0, rt)
        r = np.array([U1 - h, y1 - y0])
        prev_res = res
        res = float(np.hypot(r[0] / max(1.0, abs(h)), r[1] / k))
        if res < tol or (res < 1e-9 and res > 0.3 * prev_res):
            break
        Jm = np.array([[M[0, 1], dk[0]], [M[1, 1] - 1.0, dk[1]]])
        d = np.linalg.solve(Jm, -r)
        sc = max(abs(d[0]) / (0.3 * abs(y0)), abs(d[1]) / (0.3 * k))
        if sc > 1:
            d = d / sc
        y0, k = y0 + d[0], k + d[1]
        if not (y0 > 0 and k > 0):
            raise RuntimeError("newton_fixed_h left the domain")
    U1, y1, T, M, dk = crest_return_var(mu, k, h, y0, rt)
    m_t, m_s = multipliers(M, mu, T)
    return dict(mu=mu, k=float(k), rho=float(k / kH(mu)), U0=float(h), y0=float(y0), T=float(T), M=M,
                m_t=float(m_t), residual=float(np.hypot((U1 - h) / max(1.0, abs(h)), (y1 - y0) / k)), iters=it + 1)


def profile(mu, k, U0, y0, T, nmin=4000, settle_max=60, settle_tol=1e-10, rt=RT, record=True):
    """One period from the crest with pulled-back (eta, zeta, K, W) and exact quadratures.
    Returns dict with arrays (if record) and scalars.
    The pull-back is iterated period by period (state reset to the crest each time) until the crest values of
    (eta, zeta) repeat to settle_tol (relative)."""
    par = np.array([mu, k])
    ed = equilibrium_data(mu, k)
    rtol, atol = tolerances(1, mu, k, rt)
    s = np.zeros(12)
    s[0], s[1], s[2], s[3], s[4] = U0, 0.0, y0, ed["etaE"], ed["zetaE"]
    s[5] = 1e-3 * k
    s[6] = 1e-3 * k * mu
    ev_t = np.zeros(4)
    ev_y = np.zeros((4, 12))
    defect = np.inf
    nset = 0
    for nset in range(1, settle_max + 1):
        prev = s[3:7].copy()
        s[0], s[1], s[2] = U0, 0.0, y0
        s[7:] = 0.0
        st, t, ns, nev, nrec = integrate(1, par, s, 0.0, T, rtol, atol, -1, 0, 1, ev_t, ev_y,
                                         _EMPTY_T, _EMPTY_Y, np.inf, 50_000_000)
        if st != 0:
            raise RuntimeError(f"profile settle failed status {st}")
        newdefect = float(np.max(np.abs(s[3:5] - prev[:2]) / np.abs(s[3:5])))
        stalled = nset >= 4 and newdefect > 0.5 * defect and newdefect < 1e-7
        defect = newdefect
        if (defect < settle_tol and nset >= 2) or stalled:
            break
    start = s.copy()
    start[0], start[1], start[2] = U0, 0.0, y0
    start[7:] = 0.0
    out = dict(settle_periods=nset, settle_defect=defect)
    nrec_max = 0
    if record:
        nrec_max = 40 * nmin + 200000
    rec_t = np.zeros(nrec_max)
    rec_y = np.zeros((nrec_max, 12))
    s = start.copy()
    # both-direction event on p to capture the trough exactly (crest at t=0 and t=T are endpoints)
    st, t, ns, nev, nrec = integrate(1, par, s, 0.0, T, rtol, atol, 1, +1, 4, ev_t, ev_y,
                                     rec_t, rec_y, T / nmin if record else np.inf, 50_000_000)
    if st != 0:
        raise RuntimeError(f"profile failed status {st} nev {nev}")
    out.update(
        closure=float(max(abs(s[0] - U0), abs(s[2] - y0) / k)),
        ratio_defect=float(np.max(np.abs(s[3:7] - start[3:7]) / np.abs(start[3:7]))),
        mean_eta=float(s[7] / T), lambda_c=float(s[7] / T - mu), log_mt=float(s[7] - mu * T),
        mean_U=float(s[8] / T), q2=float(s[9] / T), q3=float(s[10] / T), mean_x=float(s[11] / T),
        U_max=float(U0), U_min=float(ev_y[0, 0]) if nev >= 1 else np.nan,
        t_trough=float(ev_t[0]) if nev >= 1 else np.nan, n_up_events=int(nev),
        eta0=float(start[3]), zeta0=float(start[4]), K0=float(start[5]), W0=float(start[6]))
    if record:
        r = rec_y[:nrec]
        out.update(t=rec_t[:nrec].copy(), U=r[:, 0].copy(), p=r[:, 1].copy(), y=r[:, 2].copy(),
                   eta=r[:, 3].copy(), zeta=r[:, 4].copy(), K=r[:, 5].copy(), W=r[:, 6].copy())
        out["x"] = k * np.exp(out["U"])
    return out


def hopf_guess(mu, h):
    """Leading-order Hopf branch at crest amplitude h: returns (y0, k)."""
    sig2 = (1 + mu) ** 2 / mu
    k = kH(mu) * (1.0 + h * h / (2.0 * (sig2 + 4.0)))
    y0 = k * (1.0 + (np.exp(h) - 1.0) / (1.0 + mu) * 1.0)
    # crest identity x-k = (1+eta)(y-k) with eta ~ mu  ->  y-k ~ k(e^h-1)/(1+mu)
    return y0, k


def unstable_seed(mu, k, r0, theta=0.0):
    ed = equilibrium_data(mu, k)
    v = np.real(np.exp(1j * theta) * ed["vec_u"])
    scale = np.array([1.0, np.sqrt(mu) if mu < 1 else mu, k])
    v = v / np.linalg.norm(v / scale)
    return np.array([0.0, 0.0, k]) + r0 * v, ed
