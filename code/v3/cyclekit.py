#!/usr/bin/env python3
"""Cycle location and gauge profiles in logarithmic variables (floating point; nothing is an enclosure).

State (U, w, c, a, b) with  p = lam - e^w,  y = e^c,  eta = e^a,  zeta = e^b,  lam = k/mu.  These keep the
RELATIVE accuracy of y, lam - p, eta and zeta through deep troughs, where x = k e^U underflows.

    find(mu, rho)             forward iteration of the crest return map from a moderate crest; if that has
                              not settled after 60 returns (slow convergence near the Hopf curve), the
                              referee's bisection on the crest amplitude; then a variational Newton polish
    from_crest(mu, k, U0, y0) the Newton polish from a given crest
    profile(mu, k, cyc)       one period with running integrals of every average used by the checks, and the
                              gauges K, W, Sigma, Xi, O sampled on the solver's own accepted steps
                              (max_step = T/2000; no interpolated points)

Every gauge is a sum of signed monomials in (q, y, eta, zeta, x, p, U).  It is evaluated in signed-log
arithmetic and returned as  ratio = sum / (sum of absolute values)  in [-1, 1], together with the logarithm
of the scale.  A sign test is "ratio against a rounding tolerance"; no monomial is allowed to underflow.

For mu > 1 the trough has eta -> mu - 1 and y -> (mu - 1)(lam - p), and W = zeta' p - eta' p' cancels to
leading order in the variables above.  Two further logarithmic variables are then integrated,
    eps = eta - (mu - 1) > 0   (eps' = zeta - eta eps),      r = y - (mu - 1)(lam - p) > 0   (r' = x - mu r),
and W is evaluated from the identity (exact; checked by sympy in `w_identity_check`)
    W = x lam - x q - zeta lam - eta zeta lam + eta zeta q + zeta r + eta eps q - eta eps r,   q = lam - p.
"""
from __future__ import annotations
import math
import numpy as np
from scipy.integrate import solve_ivp
import common  # noqa: F401  (sets sys.path for lib)
import vlog
import mono_vlib

NINT = 16      # number of running integrals


def kH(mu):
    return mu * (1.0 + mu)


def from_crest(mu, k, U0, y0):
    """Newton iteration on the crest return map (variational equations, original variables) from a given crest."""
    cyc = vlog.newton_cycle(mu, k, U0, y0, iters=40, tol=1e-12)
    if not (abs(cyc['U0']) > 1e-7 and np.isfinite(cyc['T'])):
        raise RuntimeError('Newton iteration returned the equilibrium or a non-finite state')
    # location residual of the returned crest, relative to the crest state itself: the k-scaled residual of
    # vlog.newton_cycle is not meaningful when y at the crest is hundreds of times k (deep-trough cycles)
    U1, y1, _T, _M, _Phi = vlog.crest_map(mu, k, cyc['U0'], cyc['y0'])
    cyc['res_rel'] = float(max(abs(U1 - cyc['U0']) / (1.0 + abs(cyc['U0'])), abs(y1 - cyc['y0']) / cyc['y0']))
    cyc.update(mu=mu, k=k, rho=k / kH(mu))
    return cyc


def find(mu, rho):
    k = rho * kH(mu)
    ed = mono_vlib.eq_data(mu, k)
    U0, y0 = 0.5, k + ed['zetaE'] * 0.5
    settled = False
    try:
        for _ in range(60):
            U1, y1, _T = mono_vlib.ret(mu, k, U0, y0)
            d = abs(U1 - U0) + abs(y1 - y0) / k
            U0, y0 = U1, y1
            if d < 1e-6 * (1 + abs(U0)):
                settled = True
                break
    except RuntimeError:
        settled = False
    if not settled:
        c0 = mono_vlib.find_cycle(mu, k)
        U0, y0 = c0['U0'], c0['y0']
    cyc = from_crest(mu, k, U0, y0)
    cyc['located_by'] = 'forward iteration' if settled else 'amplitude bisection'
    return cyc


# ---------------------------------------------------------------- signed-log sums
def _sl_scalar(terms):
    """terms: list of (sign, log|.|) floats.  Returns (sum/exp(M), M) with M the largest log-magnitude."""
    M = max(t[1] for t in terms)
    if M == -math.inf:
        return 0.0, M
    return sum(s * math.exp(L - M) for s, L in terms), M


def _sl_vec(terms):
    """vectorised: returns ratio = sum/sum|.| and log(sum|.|)"""
    L = np.stack([np.asarray(t[1], float) for t in terms])
    S = np.stack([np.broadcast_to(np.asarray(t[0], float), L.shape[1:]) for t in terms])
    L = np.where(S == 0, -np.inf, L)
    M = np.max(L, axis=0)
    M = np.where(np.isfinite(M), M, 0.0)
    w = np.exp(L - M)
    ab = np.sum(w, axis=0)
    return np.sum(S * w, axis=0) / np.where(ab > 0, ab, 1.0), M + np.log(np.where(ab > 0, ab, 1.0))


def _terms(mu, k, U, w, c, a, b, p, lg, sg):
    """signed monomials of the gauges; lg = log, sg = sign (scalar or vector versions)"""
    lam = k / mu
    A = 1.0 + mu
    lmu, lk, llam, lA = math.log(mu), math.log(k), math.log(lam), math.log(A)
    lx = lk + U
    lp, sp_ = lg(abs(p)), sg(p)
    lU, sU = lg(abs(U)), sg(U)
    one = sp_ * 0 + 1.0
    Sg = [(one, 2 * lp), (-sU, lU + lmu + w), (sU, lU + c)]
    Xg = [(one, 2 * lmu + 2 * w), (-one, math.log(2.0) + lmu + w + c), (one, 2 * c),
          (sp_, lp + 2 * lmu + w), (-sp_, lp + lA + c), (sp_, lp + lx)]
    K = [(sU, b + lU), (-one, a + llam), (one, a + w), (-one, c), (one, lk + 0 * U)]
    W = [(one, lx + llam), (-one, lx + w), (-one, b + llam), (one, b + w), (-one, a + b + llam), (one, a + b + w),
         (-one, lmu + b + w), (one, b + c), (one, lmu + 2 * a + w), (-one, 2 * a + c)]
    if mu != 1.0:
        s1 = math.copysign(1.0, mu - 1.0)
        lm1 = math.log(abs(mu - 1.0))
        W += [(-s1 * one, lm1 + lmu + a + w), (s1 * one, lm1 + a + c)]
    O = [(sp_, b + lp), (-one, a + lmu + w), (one, a + c), (-one, lx), (one, c)]
    V3 = [(sU, lU + lx), (-sU, lU + c), (-sp_, lp + c), (sp_, lp + lk)]
    return dict(Sg=Sg, Xg=Xg, K=K, W=W, O=O, V3=V3)


def _w_terms_alt(mu, k, U, w, a, b, le, lr, one):
    """monomials of W in the cancellation-free variables (mu > 1); see the module docstring"""
    lam = k / mu
    lk, llam = math.log(k), math.log(lam)
    lx = lk + U
    return [(one, lx + llam), (-one, lx + w), (-one, b + llam), (-one, a + b + llam), (one, a + b + w),
            (one, b + lr), (one, a + le + w), (-one, a + le + lr)]


def w_identity_check():
    """sympy: the 8-monomial form equals zeta' p - eta' p' identically (used by v3_gauges once per run)"""
    import sympy as sp
    x, q, zeta, eta, mu, lam, y = sp.symbols('x q zeta eta mu lam y')
    eps, r = eta - (mu - 1), y - (mu - 1) * q
    p, pd = lam - q, mu * q - y
    W = (x - (1 + eta) * zeta) * p - (zeta + (mu - 1) * eta - eta**2) * pd
    alt = x * lam - x * q - zeta * lam - eta * zeta * lam + eta * zeta * q + zeta * r + eta * eps * q - eta * eps * r
    return sp.expand(W - alt) == 0


def _slog(v):
    return math.log(v) if v > 0 else -math.inf


def _ssign(v):
    return 1.0 if v > 0 else (-1.0 if v < 0 else 0.0)


def gauges(mu, k, U, w, c, a, b, le=None, lr=None):
    """vectorised gauges: dict name -> (ratio, logscale), plus plain state arrays"""
    lam = k / mu
    p = lam - np.exp(w)
    with np.errstate(divide='ignore'):
        T = _terms(mu, k, U, w, c, a, b, p, lambda v: np.log(v), np.sign)
        if le is not None:
            T['W'] = _w_terms_alt(mu, k, U, w, a, b, le, lr, np.sign(p) * 0 + 1.0)
        out = {n: _sl_vec(t) for n, t in T.items()}
    out.update(U=U, p=p, eta=np.exp(a), zeta=np.exp(b), y=np.exp(c), q=np.exp(w))
    return out


def _rhs(mu, k, alt=False, extras=None):
    """alt=True (mu > 1): state is (U, w, c, a, b, le, lr, integrals); otherwise (U, w, c, a, b, integrals).
    extras: optional callable(dict of scalars) -> list of further integrands, appended after the NINT standard ones."""
    lam = k / mu
    lk = math.log(k)
    ex = vlog.ex
    n0 = 7 if alt else 5

    def f(t, s):
        U, w, c, a, b = s[0], s[1], s[2], s[3], s[4]
        q = ex(w); y = ex(c); eta = ex(a); zeta = ex(b); x = k * ex(U)
        p = lam - q
        pd = mu * q - y
        T = _terms(mu, k, U, w, c, a, b, p, _slog, _ssign)
        extra = []
        if alt:
            le, lr = s[5], s[6]
            T['W'] = _w_terms_alt(mu, k, U, w, a, b, le, lr, 1.0)
            extra = [ex(b - le) - eta, k * ex(U - lr) - mu]
        rS, mS = _sl_scalar(T['Sg']); rX, mX = _sl_scalar(T['Xg']); rK, mK = _sl_scalar(T['K'])
        rW, mW = _sl_scalar(T['W']); r3, m3 = _sl_scalar(T['V3'])
        e0 = (r3 / rS) * ex(m3 - mS)
        iK = p * (rK / rS) * ex(mK - mS)
        iW = p * (rW / rX) * ex(mW - mX)
        Sg = rS * ex(mS)
        more = [] if extras is None else list(extras(dict(U=U, p=p, pd=pd, q=q, y=y, x=x, eta=eta, zeta=zeta,
                                                             K=rK * ex(mK), W=rW * ex(mW), Sg=Sg, Xg=rX * ex(mX))))
        return [p, ex(c - w) - mu, k * ex(U - c) - 1.0, ex(b - a) + (mu - 1.0) - eta, k * ex(U - b) - 1.0 - eta] + extra + [
                eta - mu,                 # 0  <eta> - mu
                iK,                       # 1  <p K / Sigma>
                iW,                       # 2  <p W / Xi>
                e0 - mu,                  # 3  <eta_0> - mu
                abs(eta - mu),            # 4
                abs(iK),                  # 5
                abs(iW),                  # 6
                abs(e0 - mu),             # 7
                Sg / (U * U + p * p),     # 8  total polar angle
                eta * p,                  # 9
                zeta * pd,                # 10
                zeta * p,                 # 11
                eta * pd,                 # 12
                ex((lk + U) / 3.0),       # 13 x^(1/3)
                x - k,                    # 14
                abs(eta * p) + abs(zeta * pd)] + more   # 15, then the optional extra integrands
    return f


def profile(mu, k, cyc, rtol=1e-12, min_steps=2000, extras=None, n_extras=0):
    T = cyc['T']
    ric = vlog.cycle_riccati(mu, k, cyc, n=3)
    s5 = ric['S'][:5, 0]
    alt = mu > 1.0
    n0 = 5
    head = s5
    if alt:
        # crest values of eps = eta - (mu-1) and r = y - (mu-1)(lam - p); both are O(1) differences at the crest
        eps0 = math.exp(s5[3]) - (mu - 1.0)
        r0 = math.exp(s5[2]) - (mu - 1.0) * math.exp(s5[1])
        if not (eps0 > 0 and r0 > 0):
            raise RuntimeError('non-positive eps or r at the crest')
        head = np.concatenate((s5, [math.log(eps0), math.log(r0)]))
        n0 = 7
    s0 = np.concatenate((head, np.zeros(NINT + n_extras)))
    atol = np.concatenate((1e-13 * np.ones(n0), 1e-14 * np.ones(NINT + n_extras)))
    sol = solve_ivp(_rhs(mu, k, alt, extras), [0.0, T], s0, method='DOP853', rtol=rtol, atol=atol, max_step=T / min_steps)
    if sol.status != 0:
        raise RuntimeError('profile integration failed: ' + sol.message)
    g = gauges(mu, k, *sol.y[:n0])
    out = dict(T=T, t=sol.t, steps=int(sol.t.size), g=g, means=sol.y[n0:n0 + NINT, -1] / T,
               extra_means=sol.y[n0 + NINT:, -1] / T, S=sol.y[:5],
               closure=float(np.max(np.abs(sol.y[:5, -1] - sol.y[:5, 0]))),
               riccati_periods=ric['periods_used'], riccati_defect=ric['defect'],
               W_formula='cancellation-free variables (mu > 1)' if alt else 'state monomials')
    if alt:
        # consistency of the two extra variables with their definitions, where double precision resolves them
        eps_def = np.exp(sol.y[3]) - (mu - 1.0)
        r_def = np.exp(sol.y[2]) - (mu - 1.0) * np.exp(sol.y[1])
        m1 = eps_def > 1e-6 * np.exp(sol.y[3])
        m2 = r_def > 1e-6 * np.exp(sol.y[2])
        out['alt_consistency'] = float(max(np.max(np.abs(np.exp(sol.y[5][m1]) / eps_def[m1] - 1.0)) if m1.any() else 0.0,
                                           np.max(np.abs(np.exp(sol.y[6][m2]) / r_def[m2] - 1.0)) if m2.any() else 0.0))
        out['alt_closure'] = float(np.max(np.abs(sol.y[5:7, -1] - sol.y[5:7, 0])))
    return out
