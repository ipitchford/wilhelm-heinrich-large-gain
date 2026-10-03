#!/usr/bin/env python3
# --- v3.0 bundle copy of research file `lanes/planar/verify/v5_mu_gt1.py`.
# --- Changes from the original: absolute import path removed; the result dictionary is bound to RESULT instead of being written to a JSON file.
"""Referee: supporting planes (A3) for mu > 1 in deep troughs.
In the battery v1b the case (mu,rho)=(3,12) gave ~1e5 pairs with Delta/maxterm in [0, 5e-10]: for mu>1 the trough
has eta -> mu-1 and y -> (mu-1)(lam-p), so  -eta dp - dy  cancels to leading order and the sign of Delta is below
double precision in the variables (y, eta).  Cancellation-free variables (all positive on the attractor):
    q = lam - p,   r = y - (mu-1) q  (r' = x - mu r),   eps = eta - (mu-1)  (eps' = zeta - eta eps),   zeta.
Then   Delta_ab = zeta_a (U_b-U_a) + eps_a (q_b-q_a) - (r_b-r_a)   exactly.
state (U, w=log q, l=log r, e=log eps, b=log zeta):
    U' = lam - e^w,  w' = e^{l-w} - 1,  l' = k e^{U-l} - mu,  e' = e^{b-e} - (mu-1) - e^e,  b' = k e^{U-b} - mu - e^e.
usage: v5_mu_gt1.py mu rho     (mu > 1)"""
import sys, math, json
import numpy as np
from scipy.integrate import solve_ivp
from vlog import eq_data, ex, newton_cycle

mu = float(sys.argv[1]); rho = float(sys.argv[2])
assert mu > 1
k = rho * mu * (1 + mu); lam = k / mu; ll = math.log(lam)
ed = eq_data(mu, k)


def rhs(t, s):
    U, w, l, e, b = s
    ee = ex(e)
    return [lam - ex(w), ex(l - w) - 1.0, k * ex(U - l) - mu, ex(b - e) - (mu - 1.0) - ee, k * ex(U - b) - mu - ee]


def tolog(U, p, y, eta, zeta):
    q = lam - p
    return np.array([U, math.log(q), math.log(y - (mu - 1) * q), math.log(eta - (mu - 1)), math.log(zeta)])


def crest(t, s):
    return ll - s[1]
crest.direction = -1

RT, AT = 1e-13, 1e-14
t_skip = 35.0 / ed['gam']
r0 = min(1e-3, 1e-3 * math.exp(-ed['lam'].real * t_skip))
v = ed['v'] / np.linalg.norm(ed['v'])
pts = []
nph = 6
tconv = None
for j in range(nph):
    d = np.real(np.exp(2j * math.pi * j / nph) * v)
    s0 = tolog(r0 * d[0], r0 * d[1], k + r0 * d[2], ed['eta'], ed['zeta'])
    if tconv is None:                      # find convergence time from crest increments
        t = 0.0; s = s0.copy(); cr = []
        per = 2 * math.pi / ed['lam'].imag
        while True:
            sol = solve_ivp(rhs, [t, t + 20 * per], s, method='DOP853', rtol=RT, atol=AT, events=[crest])
            cr += [ye[0] for ye in sol.y_events[0]]
            s = sol.y[:, -1]; t += 20 * per
            if len(cr) > 3 and abs(cr[-1] - cr[-2]) < 1e-9 * cr[-1] and abs(cr[-2] - cr[-3]) < 1e-9 * cr[-1]:
                break
        tconv = t
        last = sol.y_events[0][-1]
    te = np.linspace(t_skip, tconv, 2400)
    sol = solve_ivp(rhs, [0, tconv], s0, method='DOP853', rtol=RT, atol=AT, t_eval=te)
    S = sol.y.T
    pS = lam - np.exp(S[:, 1])

    def arc_idx(Uc, n):
        dd = np.hypot(np.diff(Uc) / (np.ptp(Uc) + 1e-300), np.diff(pS) / (np.ptp(pS) + 1e-300))
        arc = np.concatenate(([0], np.cumsum(dd)))
        return np.clip(np.searchsorted(arc, np.linspace(0, arc[-1], n)), 0, len(arc) - 1)
    idx = np.unique(np.concatenate((np.arange(0, len(te), 16), arc_idx(S[:, 0], 150), arc_idx(np.maximum(S[:, 0], -8.0), 150))))
    pts.append(S[idx])
P = np.vstack(pts)
# cycle points: continue the converged spiral for one more period (it is on the cycle to 1e-9 relative) -- plus
# an independent Newton cycle to get T
qL = math.exp(last[1]); yL = math.exp(last[2]) + (mu - 1) * qL
cyc = newton_cycle(mu, k, last[0], yL)
solc = solve_ivp(rhs, [0, cyc['T']], s, method='DOP853', rtol=RT, atol=AT, t_eval=np.linspace(0, cyc['T'], 900))
P = np.vstack((P, solc.y.T))
U, w, l, e, b = P.T
N = len(P)
print(f"mu={mu} rho={rho}: N={N} points; Umin={U.min():.1f}; min eps = eta-(mu-1): e^{e.min():.1f}; min r: e^{l.min():.1f}")
worst = -np.inf; n_pos = 0; n_sig = 0; far_worst = -np.inf; n_close = 0
edges = np.array([-3.1, -1e-3, -1e-6, -1e-9, -1e-12, 0.0, 1e-12, 1e-9, 1e-6, 3.1]); hist = np.zeros(9, dtype=np.int64)
for i0 in range(0, N, 300):
    sl = slice(i0, min(N, i0 + 300))
    dU = U[None, :] - U[sl, None]; dw = w[None, :] - w[sl, None]; dl = l[None, :] - l[sl, None]
    with np.errstate(divide='ignore', invalid='ignore'):
        L1 = b[sl, None] + np.log(np.abs(dU))
        L2 = e[sl, None] + np.maximum(w[None, :], w[sl, None]) + np.log(-np.expm1(-np.abs(dw)))
        L3 = np.maximum(l[None, :], l[sl, None]) + np.log(-np.expm1(-np.abs(dl)))
    M = np.maximum(np.maximum(L1, L2), L3)
    with np.errstate(invalid='ignore'):
        R = np.sign(dU) * np.exp(L1 - M) + np.sign(dw) * np.exp(L2 - M) - np.sign(dl) * np.exp(L3 - M)
    ii = np.arange(sl.start, sl.stop)
    R[ii - sl.start, ii] = -1.0
    R[~np.isfinite(R)] = 9.0
    delta = np.maximum(np.maximum(np.abs(dU) / (1 + np.abs(U[sl, None])), np.abs(dw) / (1 + np.abs(w[sl, None]))),
                       np.abs(dl) / (1 + np.abs(l[sl, None])))
    delta[ii - sl.start, ii] = 1.0
    with np.errstate(divide='ignore'):
        thr = 10 * (1e-11 / delta + 1e-11)
    n_sig += int((R > thr).sum()); n_pos += int((R >= 0).sum())
    close = delta < 1e-6
    n_close += int(close.sum())
    far_worst = max(far_worst, float(np.max(np.where(close, -9.0, R))))
    worst = max(worst, float(R.max()))
    hist += np.histogram(R, bins=edges)[0]
print(f"   A3 in cancellation-free variables: {N*(N-1)} ordered pairs; worst Delta/maxterm = {worst:.3e}; #(>=0) = {n_pos};"
      f" significant = {n_sig}; worst among pairs with rel. separation > 1e-6: {far_worst:.3e} ({n_close} closer pairs)")
print(f"   histogram over {edges.tolist()}: {hist.tolist()}")
# W = zeta' p - eta' p' in the same variables
one = np.ones(N)
lx = math.log(k) + U
lm, lm1 = math.log(mu), math.log(mu - 1)
terms = [(one, lx + ll), (-one, lx + w), (-one, lm + b + ll), (one, lm + b + w), (-one, e + b + ll), (one, e + b + w),
         (-one, b + w), (one, b + l), (one, lm1 + e + w), (-one, lm1 + e + l), (one, 2 * e + w), (-one, 2 * e + l)]
L = np.stack([t[1] for t in terms]); Sg = np.stack([t[0] for t in terms]); Mx = L.max(axis=0)
Wn = np.sum(Sg * np.exp(L - Mx), axis=0)
print(f"   W/maxterm: min = {Wn.min():.3e}   (W>0 claimed in (A5))")
RESULT = dict(mu=mu, rho=rho, N=N, worst=worst, n_nonneg=n_pos, significant=n_sig, far_worst=far_worst, hist=hist.tolist(),
              W_rel_min=float(Wn.min()), U_min=float(U.min()), n_close=n_close, pairs=int(N * (N - 1)))
