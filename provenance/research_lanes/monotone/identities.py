#!/usr/bin/env python3
"""Numerical verification (relative residuals) of the averaged / window identities recorded in NOTES.md."""
import numpy as np, warnings
from mono_lib import *
warnings.filterwarnings('ignore')


def rel(a, b):
    return abs(a - b) / (abs(a) + abs(b) + 1e-300)


for mu, rho in [(1.0, 1.3), (0.5, 2.0), (0.2, 1.5), (3.0, 1.5), (1.0, 1.05), (0.3, 1.2)]:
    k = rho * kH(mu); A = 1 + mu; lam_p = k / mu
    c = find_cycle(mu, k)
    d = floquet_data(mu, k, c, n=20001)
    t = d['t']; T = d['T']
    U, p, y, x, pp, yp = d['U'], d['p'], d['y'], d['x'], d['pp'], d['yp']
    eta, zeta, r, s, ah, bh, ch = d['eta'], d['zeta'], d['r'], d['s'], d['ah'], d['bh'], d['ch']
    lam, gam = d['lam'], d['gam']
    chi = lambda g: g * (g - 1) * (g - mu)
    out = {}
    out['<r>=1+<eta>'] = rel(avg(r, t), 1 + avg(eta, t))
    out['<x a^>=chi(gamma)<a^>'] = rel(avg(x * ah, t), chi(gam) * avg(ah, t))
    out['<b^>=gamma<a^>'] = rel(avg(bh, t), gam * avg(ah, t))
    out['<c^>=gamma(gamma-mu)<a^>'] = rel(avg(ch, t), gam * (gam - mu) * avg(ah, t))
    out['<chi(r)>=k'] = rel(avg(chi(r), t), k)
    q = lam_p - p
    out['(f) lam<c^q+b^y+a^x>=<a^x(lambda-A-2p)>'] = rel(lam * avg(ch * q + bh * y + ah * x, t), avg(ah * x * (lam_p - A - 2 * p), t))
    # unnormalised-sign weak adjoint Psi0 = z' x e_s^ (periodic)
    P = np.array([bh * yp - ch * pp, ah * yp + ch * p, -(bh * p + ah * pp)])
    out['<P2>=-(1+lam)<P3>'] = rel(avg(P[1], t), -(1 + lam) * avg(P[2], t))
    out['<P1>=(mu+lam)<P2>'] = rel(avg(P[0], t), (mu + lam) * avg(P[1], t))
    out['lam<P1>=<x P3>'] = rel(lam * avg(P[0], t), avg(x * P[2], t))
    out['-<P3>=<a^(2r-gamma)p>'] = rel(-avg(P[2], t), avg(ah * (2 * r - gam) * p, t))
    out['<x P3>=-<a^ x p(2r-gamma-p)>'] = rel(avg(x * P[2], t), -avg(ah * x * p * (2 * r - gam - p), t))
    sig = P[0] * pp + P[1] * d['ppp'] + P[2] * d['ypp']
    g = 1 - np.exp(U) + U * np.exp(U)
    kap = -(P[0] * U + P[1] * p + P[2] * (y - k))
    out['lam<sig>=-<P3 x p^2>'] = rel(lam * avg(sig, t), -avg(P[2] * x * p ** 2, t))
    out['lam<kap>=-k<P3 g>'] = rel(lam * avg(kap, t), -k * avg(P[2] * g, t))
    # compound Perron vector e_t ^ z' proportional to (1, eta, zeta)
    et = d['et']
    w = np.array([et[1] * p - et[0] * pp, et[0] * yp - et[2] * p, et[1] * yp - et[2] * pp])
    out["(e_t^z')_2/(e_t^z')_1 = eta"] = np.max(np.abs(w[1] / w[0] - eta)) / np.max(eta)
    out["(e_t^z')_3/(e_t^z')_1 = zeta"] = np.max(np.abs(w[2] / w[0] - zeta)) / np.max(zeta)
    # Perron eigenvector of crest return derivative = (1, zeta(0))
    M = c['M']; v = np.array([1.0, zeta[0]])
    out['M(1,zeta0)=m_t(1,zeta0)'] = np.linalg.norm(M @ v - c['m_t'] * v) / np.linalg.norm(v)
    # window identity at t1: int_{t1}^{t1+T} e^{-lam t} P3 x (p^2 + alpha*) dt = 0, alpha* = |p|(r-mu)+|y'|/r at t1
    psi = pp + r * p
    zz = [a for a, b in zeros_of(psi, t) if b > 0]
    i1 = int(np.argmin(np.abs(t / T - zz[0])))
    alpha = abs(p[i1]) * (r[i1] - mu) + abs(yp[i1]) / r[i1]
    idx = np.r_[i1:len(t) - 1, 0:i1 + 1]
    tt = np.r_[t[i1:len(t) - 1], t[0:i1 + 1] + T]
    integrand = np.exp(-lam * tt) * (P[2] * x * (p ** 2 + alpha))[idx]
    scale = np.trapezoid(np.abs(integrand), tt)
    out['window orthogonality (alpha*)'] = abs(np.trapezoid(integrand, tt)) / scale
    # e^{-lam T}-1 from the three window identities (x, x p^2, k g): must all equal m-1
    m = np.exp(-lam * T)
    wgt = np.exp(-lam * tt)
    I_x = np.trapezoid(wgt * (P[2] * x)[idx], tt); I_p2 = np.trapezoid(wgt * (P[2] * x * p ** 2)[idx], tt); I_g = np.trapezoid(wgt * (P[2] * k * g)[idx], tt)
    e1 = np.exp(-lam * t[i1])
    out['(m-1)P1(t1) = -int P3 x'] = rel((m - 1) * e1 * P[0][i1], -I_x)
    out['(m-1)sig(t1) = int P3 x p^2'] = rel((m - 1) * e1 * sig[i1], I_p2)
    out['(m-1)kap(t1) = int P3 k g'] = rel((m - 1) * e1 * kap[i1], I_g)
    worst = max(out.values())
    print(f"mu={mu} rho={rho} lam={lam:+.6f}: worst residual {worst:.2e}")
    for kname, v_ in out.items():
        if v_ > 1e-4:
            print(f"     {kname}: {v_:.2e}")
