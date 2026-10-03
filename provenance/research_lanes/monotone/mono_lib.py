#!/usr/bin/env python3
"""Lane 'monotone': Floquet data along a Wilhelm-Heinrich cycle.

For a cycle (crest phase at t=0) this computes on a uniform grid over one period
  z=(U,p,y), z'=(p,p',y'), z''
  eta, zeta                     positive adjoint ratios (strong adjoint n=(zeta,-eta,-1))
  r=b/a, s=c/a                  strong-stable ratios: e_s = (a,-b,-c), a,b,c>0, a'=-b, b'=-mu b-c, c'=-x a-c
  ah,bh,ch                      periodic normalisation e^{gamma t}(a,b,c), gamma = <r> = 1+<eta>
  nu,nup,dy_t                   weak Floquet vector e_t=(nu,nu',dy), multiplied by e^{-lambda_c t} (periodic)
  Psi (3 x n)                   periodic weak adjoint  Psi^ = z' x e_s^  (annihilates z' and e_s)
  sig                           sigma^ = Psi^ . z''
All floating point, single process.
"""
import sys
import numpy as np
from scipy.integrate import solve_ivp
sys.path.insert(0, '/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03/work')
from whlib import find_cycle, kH, equilibrium_data, crest_return  # noqa

RT, AT = 1e-11, 1e-13


def cumtrapz0(f, t):
    out = np.zeros_like(f)
    out[1:] = np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(t))
    return out


def floquet_data(mu, k, cyc, n=6001):
    T = cyc['T']
    ed = equilibrium_data(mu, k)
    A = 1 + mu
    t = np.linspace(0, T, n)

    # forward: state + (eta,zeta) (forward-stable Riccati)
    def f(tt, s):
        U, p, y, eta, zeta = s
        x = k * np.exp(U)
        return [p, k - mu * p - y, x - y, zeta + (mu - 1) * eta - eta * eta, x - (1 + eta) * zeta]
    s0 = np.array([cyc['U0'], 0.0, cyc['y0'], ed['etaE'], ed['zetaE']])
    for _ in range(5):
        sol = solve_ivp(f, [0, T], s0, method='DOP853', rtol=RT, atol=AT)
        s0 = sol.y[:, -1].copy(); s0[:3] = [cyc['U0'], 0.0, cyc['y0']]
    sol = solve_ivp(f, [0, T], s0, method='DOP853', rtol=RT, atol=AT, t_eval=t, dense_output=True)
    U, p, y, eta, zeta = sol.y
    x = k * np.exp(U)
    pp = k - mu * p - y
    yp = x - y
    ppp = -mu * pp - yp
    ypp = x * p - yp
    lam = np.trapezoid(eta, t) / T - mu
    gam = 1 + mu + lam
    dense = sol.sol

    # backward: strong-stable ratios r=b/a, s=c/a :  r' = r^2 - mu r - s ,  s' = s r - s - x
    def g(tau, w):  # tau = T - t
        tt = T - tau
        Ux = dense(tt)[0]
        xx = k * np.exp(Ux)
        r, s = w
        return [-(r * r - mu * r - s), -(s * r - s - xx)]
    gE = ed['gammaE']
    w0 = np.array([gE, gE * (gE - mu)])  # equilibrium values r=gamma, s=gamma(gamma-mu)
    for _ in range(5):
        sb = solve_ivp(g, [0, T], w0, method='DOP853', rtol=RT, atol=AT)
        w0 = sb.y[:, -1].copy()
    sb = solve_ivp(g, [0, T], w0, method='DOP853', rtol=RT, atol=AT, t_eval=T - t[::-1])
    r = sb.y[0][::-1].copy(); s = sb.y[1][::-1].copy()
    gam_r = np.trapezoid(r, t) / T
    ah = np.exp(cumtrapz0(gam_r - r, t))
    ah /= np.exp(np.trapezoid(np.log(ah), t) / T)  # geometric mean 1
    bh = r * ah; ch = s * ah

    # weak Floquet: nu'' = (eta-mu) nu' - zeta nu
    def h(tt, v):
        st = dense(tt)
        e, zt = st[3], st[4]
        return [v[1], (e - mu) * v[1] - zt * v[0], v[3], (e - mu) * v[3] - zt * v[2]]
    sv = solve_ivp(h, [0, T], [1, 0, 0, 1.0], method='DOP853', rtol=RT, atol=AT, t_eval=t)
    Mon = np.array([[sv.y[0, -1], sv.y[2, -1]], [sv.y[1, -1], sv.y[3, -1]]])
    ev, V = np.linalg.eig(Mon)
    # eigenvalue 1 has eigenvector (p(0),p'(0)) = (0,p'(0)); the other is m_t
    i = int(np.argmax(np.abs(V[0, :] / np.linalg.norm(V, axis=0))))  # the one with nonzero nu(0)
    m_t = ev[i].real
    v0 = V[:, i].real; v0 = v0 / v0[0]
    nu = v0[0] * sv.y[0] + v0[1] * sv.y[2]
    nup = v0[0] * sv.y[1] + v0[1] * sv.y[3]
    damp = np.exp(-lam * t)
    nu_h = nu * damp; nup_h = nup * damp          # e^{-lam t} * (nu, nu')
    nupp_h = ((eta - mu) * nup - zeta * nu) * damp
    dy_h = -(nupp_h + mu * nup_h)                 # delta y = -(delta p' + mu delta p)

    # weak adjoint (periodic): Psi^ = z' x e_s^ , e_s^ = (ah,-bh,-ch)
    Psi = np.array([bh * yp - ch * pp, ah * yp + ch * p, -(bh * p + ah * pp)])
    # normalise so that Psi . e_t^ = 1 at t=0 (constant in t up to numerical error)
    et = np.array([nu_h, nup_h, dy_h])
    pair = np.sum(Psi * et, axis=0)
    sgn_scale = pair[0]
    Psi = Psi / sgn_scale
    sig = Psi[0] * pp + Psi[1] * ppp + Psi[2] * ypp
    return dict(mu=mu, k=k, T=T, t=t, U=U, p=p, y=y, x=x, pp=pp, yp=yp, ppp=ppp, ypp=ypp,
                eta=eta, zeta=zeta, lam=lam, gam=gam, gam_r=gam_r, r=r, s=s, ah=ah, bh=bh, ch=ch,
                nu=nu_h, nup=nup_h, dy=dy_h, m_t_2nd=m_t, Psi=Psi, sig=sig, pair=pair / sgn_scale,
                et=et)


def zeros_of(f, t):
    """phases (fraction of period) of sign changes of f, with direction (+1 up, -1 down)."""
    out = []
    for i in range(len(f) - 1):
        if f[i] == 0:
            continue
        if f[i] * f[i + 1] < 0:
            tz = t[i] - f[i] * (t[i + 1] - t[i]) / (f[i + 1] - f[i])
            out.append((tz / t[-1], +1 if f[i + 1] > f[i] else -1))
    return out


def avg(f, t):
    return np.trapezoid(f, t) / (t[-1] - t[0])


if __name__ == '__main__':
    for mu, rho in [(1.0, 1.3), (0.2, 2.0)]:
        k = rho * kH(mu)
        c = find_cycle(mu, k)
        d = floquet_data(mu, k, c)
        t = d['t']
        print(f"mu={mu} rho={rho} T={d['T']:.5f} lam(Riccati)={d['lam']:+.8f} lam(M)={c['lambda_c']:+.8f} "
              f"gamma: 1+mu+lam={d['gam']:.8f} <r>={d['gam_r']:.8f}  m_t(M)={c['m_t']:.8f} m_t(2nd)={d['m_t_2nd']:.8f}")
        print("   pairing Psi.e_t range", d['pair'].min(), d['pair'].max(),
              " Psi.z' max", np.abs(d['Psi'][0] * d['p'] + d['Psi'][1] * d['pp'] + d['Psi'][2] * d['yp']).max(),
              " periodicity r,s,a:", abs(d['r'][0] - d['r'][-1]), abs(d['s'][0] - d['s'][-1]), abs(d['ah'][0] - d['ah'][-1]))
        # identity: sigma^' = lam sigma^ + Psi3 x p^2  -> lam <sig> = -<Psi3 x p^2>
        print("   lam*<sig> =", d['lam'] * avg(d['sig'], t), "  -<Psi3 x p^2> =", -avg(d['Psi'][2] * d['x'] * d['p'] ** 2, t))
