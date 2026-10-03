#!/usr/bin/env python3
"""Lane 'supersolution': shared data builder.

Builds a cache of cycles on a (mu, rho) grid with dense periodic profiles of
  state (U,p,y), adjoint ratios (eta,zeta), strong-stable Riccati ratios (beta,kappa)=(b/a,c/a),
  and the periodic part ahat of the strong-stable U-component.
Floating point only.
"""
import sys, os, time, pickle
import numpy as np
from scipy.integrate import solve_ivp
ROOT = '/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03'
sys.path.insert(0, ROOT + '/work')
from whlib import find_cycle, equilibrium_data, kH, field

HERE = ROOT + '/work/lanes/supersolution'
CACHE = HERE + '/cycles.pkl'

MUS = [0.05, 0.1, 0.2, 0.5, 1.0]
RHOS = [1.001, 1.003, 1.01, 1.03, 1.1, 1.3, 1.6, 2.0, 3.0, 5.0, 10.0]


def profile(mu, k, cyc, n=None):
    T = cyc['T']
    rho = k / kH(mu)
    if n is None:
        n = 4096 if rho <= 3.5 else 32768
    ed = equilibrium_data(mu, k)
    # integrate in log-chemical variables (U, v=log Y, w=log Z): relative accuracy in the deep trough
    lam = k / mu
    def flog(t, s):
        U, v, w = s
        return [lam * (-np.expm1(v)), mu * np.expm1(w - v), np.expm1(U + np.log(k) - (w + np.log(k)))]
    z0 = [cyc['U0'], 0.0, np.log(cyc['y0'] / k)]
    sol = solve_ivp(flog, [0, T], z0, method='DOP853', rtol=1e-12, atol=1e-13, dense_output=True)
    closure = float(np.max(np.abs(sol.y[:, -1] - np.array(z0))))

    def Uof(t):
        return sol.sol(np.mod(t, T))[0]

    # forward Riccati (eta,zeta): attracting forward
    def fr(t, s):
        eta, zeta = s
        x = k * np.exp(Uof(t))
        return [zeta + (mu - 1) * eta - eta * eta, x - (1 + eta) * zeta]
    s = np.array([ed['etaE'], ed['zetaE']])
    for _ in range(6):
        r = solve_ivp(fr, [0, T], s, method='DOP853', rtol=1e-11, atol=1e-13)
        s_new = r.y[:, -1]
        done = np.max(np.abs(s_new - s)) < 1e-10 * (1 + np.max(np.abs(s)))
        s = s_new
        if done:
            break
    t = np.linspace(0, T, n, endpoint=False)
    tt = np.append(t, T)
    r = solve_ivp(fr, [0, T], s, method='DOP853', rtol=1e-11, atol=1e-13, t_eval=tt)
    eta, zeta = r.y[0, :-1], r.y[1, :-1]
    # backward Riccati (beta,kappa): attracting backward. beta'=beta^2-mu beta-kappa, kappa'=kappa(beta-1)-x
    gE = ed['gammaE']

    def br(sg, s):  # sg = -t
        beta, kappa = s
        x = k * np.exp(Uof(-sg))
        return [-(beta * beta - mu * beta - kappa), -(kappa * (beta - 1) - x)]
    s2 = np.array([gE, gE * (gE - mu)])
    for _ in range(6):
        r2 = solve_ivp(br, [0, T], s2, method='DOP853', rtol=1e-11, atol=1e-13)
        s_new = r2.y[:, -1]
        done = np.max(np.abs(s_new - s2)) < 1e-10 * (1 + np.max(np.abs(s2)))
        s2 = s_new
        if done:
            break
    # state at t=0 (== -sg=0); go backward one period, evaluate at sg = T - t  -> t in (0,T]
    sg_eval = np.linspace(0, T, n + 1)
    r2 = solve_ivp(br, [0, T], s2, method='DOP853', rtol=1e-11, atol=1e-13, t_eval=sg_eval)
    # value at time t_j = T - sg_j (mod T). reverse to increasing t
    beta = r2.y[0, ::-1][:-1]   # index j -> t = T - sg_{n-j} = t_j  for j=0..n-1 (t_0 = 0 == T)
    kappa = r2.y[1, ::-1][:-1]
    st = sol.sol(t)
    U, v, wz = st
    logY = v.copy(); logZ = wz.copy()
    p = -lam * np.expm1(v); y = k * np.exp(wz)
    x = k * np.exp(U)
    gamma = float(np.mean(beta))
    # ahat = exp(gamma t - int_0^t beta): spectral antiderivative of the mean-free part
    bt = beta - gamma
    Bh = np.fft.rfft(bt)
    w = 2 * np.pi * np.fft.rfftfreq(n, d=T / n)
    Ih = np.zeros_like(Bh)
    Ih[1:] = Bh[1:] / (1j * w[1:])
    Ib = np.fft.irfft(Ih, n)
    loga = -(Ib - Ib.mean())
    return dict(mu=mu, k=k, rho=rho, T=T, t=t, U=U, p=p, y=y, x=x, logY=logY, logZ=logZ, eta=eta, zeta=zeta, beta=beta, kappa=kappa,
                log_ahat=loga, gamma_beta=gamma, lam_eta=float(np.mean(eta) - mu), lam_beta=float(gamma - 1 - mu),
                lam_M=cyc['lambda_c'], closure=closure, ed=ed, cyc=cyc)


def build(mus=MUS, rhos=RHOS, verbose=True):
    out = {}
    for mu in mus:
        guess = None
        # go from moderate rho down to Hopf, then up, seeding with neighbours
        order_down = [r for r in sorted(rhos, reverse=True) if r <= 1.3]
        order_up = [r for r in sorted(rhos) if r > 1.3]
        cyc_by = {}
        for seq in (order_down, order_up):
            guess = None if seq is order_down else cyc_by.get(1.3, None)
            if guess is not None:
                guess = (guess['U0'], guess['y0'])
            for rho in seq:
                k = rho * kH(mu)
                t0 = time.time()
                try:
                    if guess is not None and seq is order_down:
                        # scale amplitude ~ sqrt(rho-1)
                        pass
                    c = find_cycle(mu, k, guess=guess)
                    if c['U0'] < 1e-4:
                        raise RuntimeError('collapsed to E')
                    cyc_by[rho] = c
                    guess = (c['U0'], c['y0'])
                    pr = profile(mu, k, c)
                    out[(mu, rho)] = pr
                    if verbose:
                        print(f"mu={mu:<5} rho={rho:<5} T={c['T']:9.4f} Umax={c['U0']:8.4f} Umin={pr['U'].min():10.3f} "
                              f"lam_M={c['lambda_c']:+.6e} lam_eta={pr['lam_eta']:+.6e} lam_beta={pr['lam_beta']:+.6e} "
                              f"clos={pr['closure']:.1e} res={c['residual']:.0e} [{time.time()-t0:.1f}s]", flush=True)
                except Exception as e:
                    print(f"mu={mu} rho={rho} FAILED {type(e).__name__}: {e} [{time.time()-t0:.1f}s]", flush=True)
                    guess = None
    with open(CACHE, 'wb') as f:
        pickle.dump(out, f)
    return out


def load():
    with open(CACHE, 'rb') as f:
        return pickle.load(f)


def avg(f):
    """periodic trapezoid mean on the uniform grid"""
    return float(np.mean(f))


def ddt(f, T):
    """spectral derivative of a periodic sample (use only for smooth, well-resolved profiles)"""
    n = len(f)
    w = 2 * np.pi * np.fft.rfftfreq(n, d=T / n)
    return np.fft.irfft(1j * w * np.fft.rfft(f), n)


if __name__ == '__main__':
    build()
