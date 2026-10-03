#!/usr/bin/env python3
"""Shared floating-point toolkit for the Wilhelm--Heinrich oscillator.

Coordinates (paper v2.0, eq. (1.5)):   U' = p,  p' = k - mu*p - y,  y' = k*exp(U) - y
with x = k*exp(U), equilibrium E = (0, 0, k), Hopf at k_H = mu*(1+mu), rho = k/k_H > 1.
Adjoint ratios:  eta' = zeta + (mu-1)*eta - eta**2,   zeta' = x - (1+eta)*zeta.
Weak transverse exponent: lambda_c = <eta> - mu = log(m_t)/T.  Attracting  <=>  lambda_c < 0.

Everything here is ordinary floating point (DOP853); nothing is an interval certificate.
"""
from __future__ import annotations
import numpy as np
from scipy.integrate import solve_ivp

RTOL, ATOL = 1e-11, 1e-13


def kH(mu):
    return mu * (1.0 + mu)


def field(mu, k):
    def f(t, z):
        U, p, y = z
        return [p, k - mu * p - y, k * np.exp(U) - y]
    return f


def equilibrium_data(mu, k):
    """gamma_E, eta_E, zeta_E and the unstable complex eigenpair of E."""
    roots = np.roots([1.0, 1.0 + mu, mu, k])
    real = [r.real for r in roots if abs(r.imag) < 1e-12]
    gammaE = -min(real)
    etaE = gammaE - 1.0
    zetaE = etaE * (etaE + 1.0 - mu)
    J = np.array([[0, 1, 0], [0, -mu, -1], [k, 0, -1.0]])
    w, V = np.linalg.eig(J)
    i = int(np.argmax(w.imag))
    return dict(gammaE=gammaE, etaE=etaE, zetaE=zetaE, lam_u=w[i], vec_u=V[:, i])


def _var_rhs(mu, k):
    def f(t, s):
        U, p, y = s[:3]
        x = k * np.exp(U)
        Phi = s[3:12].reshape(3, 3)
        J = np.array([[0, 1, 0], [0, -mu, -1], [x, 0, -1.0]])
        return np.concatenate(([p, k - mu * p - y, x - y], (J @ Phi).ravel()))
    return f


def _down(t, s):
    return s[1]
_down.direction = -1
_down.terminal = True


def crest_return(mu, k, U0, y0, tmax=None, with_derivative=True):
    """First return to the crest section {p=0, p'<0} from (U0, 0, y0).
    Returns (U1, y1, T, M) with M the 2x2 section derivative in (U,y)."""
    s0 = np.concatenate(([U0, 0.0, y0], np.eye(3).ravel()))
    f = _var_rhs(mu, k)
    # step off the section first so the terminal event does not fire at t=0
    t_off = 1e-6
    a = solve_ivp(f, [0, t_off], s0, method="DOP853", rtol=RTOL, atol=ATOL)
    if tmax is None:
        tmax = 400.0 * max(1.0, 1.0 / np.sqrt(mu))
    b = solve_ivp(f, [t_off, tmax], a.y[:, -1], method="DOP853", rtol=RTOL, atol=ATOL, events=_down)
    if b.status != 1:
        raise RuntimeError("no crest return found")
    T = float(b.t_events[0][0])
    s = b.y_events[0][0]
    U1, p1, y1 = s[:3]
    if not with_derivative:
        return U1, y1, T, None
    Phi = s[3:12].reshape(3, 3)
    fz = np.array([p1, k - mu * p1 - y1, k * np.exp(U1) - y1])
    cols = []
    for v in (np.array([1.0, 0, 0]), np.array([0, 0, 1.0])):
        d = Phi @ v
        d = d - fz * (d[1] / fz[1])
        cols.append([d[0], d[2]])
    M = np.array(cols).T
    return U1, y1, T, M


def find_cycle(mu, k, guess=None, warm=40, newton=30, tol=1e-11):
    """Locate a periodic orbit through a crest (U0, 0, y0) by warm-up returns then Newton.
    Returns dict(U0, y0, T, M, m_t, m_s, lambda_c, residual)."""
    if guess is None:
        ed = equilibrium_data(mu, k)
        v = np.real(ed["vec_u"])
        z = np.array([0.0, 0.0, k]) + 1e-3 * v / np.linalg.norm(v)
        sol = solve_ivp(field(mu, k), [0, 50.0 / max(ed["lam_u"].real, 1e-3)], z, method="DOP853",
                        rtol=1e-9, atol=1e-11, events=_down3)
        if len(sol.t_events[0]) == 0:
            raise RuntimeError("warm-up found no crest")
        zc = sol.y_events[0][-1]
        U0, y0 = zc[0], zc[2]
    else:
        U0, y0 = guess
    for _ in range(warm):
        U1, y1, T, _M = crest_return(mu, k, U0, y0, with_derivative=False)
        done = abs(U1 - U0) + abs(y1 - y0) / k < 1e-4
        U0, y0 = U1, y1
        if done:
            break
    res = np.inf
    for _ in range(newton):
        U1, y1, T, M = crest_return(mu, k, U0, y0)
        r = np.array([U1 - U0, y1 - y0])
        res = float(np.linalg.norm(r / np.array([1.0, k])))
        if res < tol:
            break
        d = np.linalg.solve(M - np.eye(2), -r)
        U0, y0 = U0 + d[0], y0 + d[1]
    U1, y1, T, M = crest_return(mu, k, U0, y0)
    ev = np.linalg.eigvals(M)
    m_t = float(np.max(ev.real))
    return dict(mu=mu, k=k, rho=k / kH(mu), U0=float(U0), y0=float(y0), T=T, M=M, m_t=m_t,
                m_s=float(np.exp(-(1 + mu) * T) / m_t), lambda_c=float(np.log(m_t) / T), residual=res)


def _down3(t, z):
    return z[1]
_down3.direction = -1


def cycle_profile(mu, k, cyc, n=4001, settle_periods=4):
    """Dense periodic profile on one period starting at the crest.
    Returns dict of arrays t,U,p,y,x,eta,zeta,K,W plus lambda_c_riccati = <eta>-mu.
    K = zeta*U - eta*p - y + k > 0 and W = zeta'*p - eta'*p' > 0 (paper eqs (3.5),(3.8))."""
    T = cyc["T"]
    ed = equilibrium_data(mu, k)

    def f(t, s):
        U, p, y, eta, zeta, K, W = s
        x = k * np.exp(U)
        return [p, k - mu * p - y, x - y,
                zeta + (mu - 1) * eta - eta * eta,
                x - (1 + eta) * zeta,
                -(1 + eta) * K + (x * U - x + k),
                -(1 + eta) * W + x * p * p]
    s = np.array([cyc["U0"], 0.0, cyc["y0"], ed["etaE"], ed["zetaE"], 0.0, 0.0])
    for _ in range(settle_periods):
        sol = solve_ivp(f, [0, T], s, method="DOP853", rtol=RTOL, atol=ATOL)
        s = sol.y[:, -1].copy()
        s[0], s[1], s[2] = cyc["U0"], 0.0, cyc["y0"]
    t = np.linspace(0, T, n)
    sol = solve_ivp(f, [0, T], s, method="DOP853", rtol=RTOL, atol=ATOL, t_eval=t)
    U, p, y, eta, zeta, K, W = sol.y
    lam = float(np.trapezoid(eta, t) / T - mu)
    return dict(t=t, U=U, p=p, y=y, x=k * np.exp(U), eta=eta, zeta=zeta, K=K, W=W,
                lambda_c_riccati=lam, periodicity_defect=float(np.max(np.abs(sol.y[:, -1] - sol.y[:, 0]))))


def unstable_manifold_orbit(mu, k, r0=1e-7, theta=0.0, t_end=None, n=20001):
    """One orbit on W^u(E): the spiral that fills the attractor disc inside the innermost cycle.
    Carries the pull-back adjoint ratios eta,zeta and the convexity gauges K, W."""
    ed = equilibrium_data(mu, k)
    v = np.real(np.exp(1j * theta) * ed["vec_u"])
    z0 = np.array([0.0, 0.0, k]) + r0 * v / np.linalg.norm(v)
    if t_end is None:
        t_end = (np.log(1.0 / r0) + 12.0) / ed["lam_u"].real

    def f(t, s):
        U, p, y, eta, zeta, K, W = s
        x = k * np.exp(U)
        return [p, k - mu * p - y, x - y,
                zeta + (mu - 1) * eta - eta * eta,
                x - (1 + eta) * zeta,
                -(1 + eta) * K + (x * U - x + k),
                -(1 + eta) * W + x * p * p]
    s0 = np.concatenate((z0, [ed["etaE"], ed["zetaE"], 0.0, 0.0]))
    t = np.linspace(0, t_end, n)
    sol = solve_ivp(f, [0, t_end], s0, method="DOP853", rtol=RTOL, atol=1e-16, t_eval=t)
    U, p, y, eta, zeta, K, W = sol.y
    return dict(t=sol.t, U=U, p=p, y=y, x=k * np.exp(U), eta=eta, zeta=zeta, K=K, W=W, eq=ed)


def continue_in_rho(mu, rhos, guess=None):
    """Cycles along increasing rho at fixed mu, each seeded from the previous crest."""
    out = []
    for rho in rhos:
        c = find_cycle(mu, rho * kH(mu), guess=guess)
        out.append(c)
        guess = (c["U0"], c["y0"] * 1.0)
    return out


if __name__ == "__main__":
    # self-test against the v2.0 bundle receipts
    ref = [(0.1, 0.196052776, 38.028889, 0.43642544581),
           (1.0, 2.201497, 6.502509, 0.77861840302),
           (1.0, 5.841219, 11.121221, 0.030745762),
           (10.0, 196.052776, 3.802889, 0.43642544581)]
    for mu, k, Tref, mref in ref:
        c = find_cycle(mu, k)
        pr = cycle_profile(mu, k, c)
        print(f"mu={mu:<5} k={k:<12} T={c['T']:.6f} (ref {Tref:.6f})  m_t={c['m_t']:.9f} (ref {mref:.9f})"
              f"  lam_M={c['lambda_c']:.9f} lam_Riccati={pr['lambda_c_riccati']:.9f}"
              f"  minK={pr['K'].min():.3e} minW={pr['W'].min():.3e} res={c['residual']:.1e}")
