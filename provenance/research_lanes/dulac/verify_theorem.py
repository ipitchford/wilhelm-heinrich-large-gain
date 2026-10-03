#!/usr/bin/env python3
"""Numerical test of every statement of Theorem D (NOTES.md) on 20 fresh parameter sets
(mu in {0.1,0.3,0.7,1.0} x rho in {1.001,1.1,2.5,8,20}), cycles computed independently of the LP data.

Checks per cycle:
  A  V = p^2 - U p' > 0,  Q = p'^2 - p p'' > 0                       (min of scale-free ratios)
  B  pointwise identities  DV = (eta-mu)V + pK,  DQ = (eta-mu)Q + pW,  V_UH = pK + eta V,
     W = x p - (1+eta-mu) y' + zeta (y-k)                              (max relative residual)
  C  lambda_c = <eta>-mu = -<pK/V> = -<pW/Q> = log(m_t)/T ;  <eta0> = <eta1> = mu
  D  sign pattern  (eta-eta0)*p <= 0, (eta-eta1)*p <= 0
  E  <eta p> < 0, <zeta p'> < 0, <eta p><zeta p'> > <zeta p>^2 ; eta_low > eta_up, y_low > y_up at equal U;
     zeta_right > zeta_left at equal p
"""
import sys, time
import numpy as np
from scipy.integrate import solve_ivp
sys.path.insert(0, '/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03/work')
import whlib
from gen_data import rhs, eqdata, psi, phi


def cycle_log(mu, k, n=12001, settle=6, guess=None):
    c = whlib.find_cycle(mu, k, guess=guess)
    ed = eqdata(mu, k); lam = k / mu
    T = c['T']
    f = rhs(mu, k, ed['gammaE'])
    s = np.array([c['U0'], np.log(c['y0'] / k), 0.0, np.log(ed['etaE']), np.log(ed['zetaE']), 0.0, 0.0, 0.0])
    for _ in range(settle):
        sol = solve_ivp(f, [0, T], s, method='DOP853', rtol=1e-11, atol=1e-13)
        s = sol.y[:, -1].copy()
        s[0], s[1], s[2] = c['U0'], np.log(c['y0'] / k), 0.0
    t = np.linspace(0, T, n)
    sol = solve_ivp(f, [0, T], s, method='DOP853', rtol=1e-11, atol=1e-13, t_eval=t)
    U, d, e, le, lz, lK, lW, lL = sol.y
    p = -lam * np.expm1(e); H = k * np.expm1(d); y = k * np.exp(d); x = k * np.exp(U)
    eta = np.exp(le); zeta = np.exp(lz); K = np.exp(lK); W = np.exp(lW)
    return dict(mu=mu, k=k, lam=lam, T=T, m_t=c['m_t'], lam_M=c['lambda_c'], t=t, U=U, p=p, H=H, y=y, x=x, eta=eta, zeta=zeta,
                K=K, W=W, pd=k * (np.exp(e) - np.exp(d)), crest=(c['U0'], c['y0']), closure=float(np.max(np.abs(sol.y[:5, -1] - sol.y[:5, 0]))))


def checks(c):
    mu, k = c['mu'], c['k']; A = 1 + mu
    t, U, p, H, y, x, eta, zeta, K, W, pd = (c[s] for s in ('t', 'U', 'p', 'H', 'y', 'x', 'eta', 'zeta', 'K', 'W', 'pd'))
    T = c['T']
    Hd = x - y
    V = p * p - U * pd
    Q = pd * (-H) + p * Hd                 # p'^2 - p p'' = p'(k-y) + p(x-y)
    VUH = U * Hd - p * H
    mean = lambda f: np.trapezoid(f, t) / T
    out = {}
    out['minV'] = np.min(V / (U * U + p * p)); out['minQ'] = np.min(Q / (p * p + pd * pd))
    # identities (state formulas for DV, DQ)
    DV = -mu * V + VUH
    DQ = x * V - A * Q + k * psi(U) * pd
    r1 = np.max(np.abs(DV - ((eta - mu) * V + p * K)) / (np.abs(DV) + np.abs(V) + 1e-300))
    r2 = np.max(np.abs(DQ - ((eta - mu) * Q + p * W)) / (np.abs(DQ) + np.abs(Q) + 1e-300))
    r3 = np.max(np.abs(VUH - (p * K + eta * V)) / (np.abs(VUH) + eta * V))
    Wst = x * p - (1 + eta - mu) * Hd + zeta * H
    r4 = np.max(np.abs(Wst - W) / (np.abs(x * p) + np.abs((1 + eta - mu) * Hd) + np.abs(zeta * H)))
    out['id'] = max(r1, r2, r3, r4)
    lam_eta = mean(eta) - mu
    out['lam_eta'] = lam_eta; out['lam_M'] = c['lam_M']
    out['lam_V'] = -mean(p * K / V); out['lam_Q'] = -mean(p * W / Q)
    eta0 = VUH / V; eta1 = mu + DQ / Q
    out['eta0'] = mean(eta0) - mu; out['eta1'] = mean(eta1) - mu
    out['sign0'] = np.max((eta - eta0) * p / (np.abs(p) + 1e-300) / (1 + eta)); out['sign1'] = np.max((eta - eta1) * np.sign(p) / (1 + np.abs(eta1)))
    out['etap'] = mean(eta * p); out['zetapd'] = mean(zeta * pd); out['zetap'] = mean(zeta * p)
    out['det'] = out['etap'] * out['zetapd'] - out['zetap']**2
    # branch comparisons at equal U
    up = p > 0; lo = p < 0
    Ug = np.linspace(U.min(), U.max(), 203)[1:-1]
    def itp(q, sel, v, f):
        o = np.argsort(v[sel]); return np.interp(q, v[sel][o], f[sel][o])
    out['d_eta'] = np.min(itp(Ug, lo, U, eta) - itp(Ug, up, U, eta))
    out['d_y'] = np.min((itp(Ug, lo, U, y) - itp(Ug, up, U, y)) / k)
    ri = pd < 0; le = pd > 0
    pg = np.linspace(p.min(), p.max(), 203)[1:-1]
    out['d_zeta'] = np.min(itp(pg, ri, p, zeta) - itp(pg, le, p, zeta))
    return out


if __name__ == '__main__':
    ok_all = True
    for mu in (0.1, 0.3, 0.7, 1.0):
        guess = None
        for rho in (1.001, 1.1, 2.5, 8.0, 20.0):
            k = rho * mu * (1 + mu)
            t0 = time.time()
            if rho < 1.01:
                # Hopf asymptotics for the crest:  eps^2 = 4 sigma (A^2+mu) D / (A mu^2),  H = mu*U at the crest
                ed = eqdata(mu, k); A = 1 + mu; D = mu * mu + 6 * mu + 1
                eps = np.sqrt(4 * ed['sigma'] * (A * A + mu) * D / (A * mu * mu))
                guess = (eps, k + mu * eps)
            try:
                c = cycle_log(mu, k, guess=guess)
            except Exception as ex:
                print(f'mu={mu} rho={rho}: cycle FAILED {ex}'); guess = None; ok_all = False; continue
            guess = c['crest']
            o = checks(c)
            lam = o['lam_eta']
            ok = (o['minV'] > 0 and o['minQ'] > 0 and o['id'] < 1e-6 and abs(o['lam_V'] - lam) < 2e-3 * abs(lam) + 1e-7
                  and abs(o['lam_Q'] - lam) < 2e-3 * abs(lam) + 1e-7 and abs(o['eta0']) < 1e-6 + 1e-3 * abs(lam) and abs(o['eta1']) < 1e-6 + 1e-3 * abs(lam)
                  and o['sign0'] <= 1e-9 and o['sign1'] <= 1e-9 and o['etap'] < 0 and o['zetapd'] < 0 and o['det'] > 0
                  and o['d_eta'] > 0 and o['d_y'] > 0 and o['d_zeta'] > 0 and lam < 0)
            ok_all &= ok
            print(f'mu={mu:<4} rho={rho:<6} T={c["T"]:8.3f} Umin={c["U"].min():9.2f} clos={c["closure"]:.0e} | minV={o["minV"]:.2e} minQ={o["minQ"]:.2e} id={o["id"]:.0e} | '
                  f'lam: eta={lam:+.6e} M={o["lam_M"]:+.6e} V={o["lam_V"]:+.6e} Q={o["lam_Q"]:+.6e} | <eta0>-mu={o["eta0"]:+.0e} <eta1>-mu={o["eta1"]:+.0e} '
                  f'sign={max(o["sign0"], o["sign1"]):+.0e} | <eta p>={o["etap"]:+.3e} <zeta pd>={o["zetapd"]:+.3e} det={o["det"]:+.2e} '
                  f'd_eta={o["d_eta"]:+.1e} d_y={o["d_y"]:+.1e} d_zeta={o["d_zeta"]:+.1e} {"OK" if ok else "CHECK"} [{time.time()-t0:.0f}s]', flush=True)
    print('ALL OK' if ok_all else 'SOME CHECKS FLAGGED')
