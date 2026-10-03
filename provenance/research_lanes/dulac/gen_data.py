#!/usr/bin/env python3
"""Dense sampling of the attractor disc D \\ {E} of the Wilhelm-Heinrich oscillator.

Everything is integrated in logarithmic variables so that relative accuracy is kept both next to E
and in the astronomically deep troughs at large rho:

    U = log X,  d = log Z = log(y/k),  e = log Y = log(1 - p/lambda)
    U' = -lambda*expm1(e),   d' = expm1(U-d),   e' = mu*expm1(d-e)
    p = -lambda*expm1(e),  H = y-k = k*expm1(d),  x = k e^U

and the pull-back adjoint ratios / convexity gauges in log form

    (log eta)'  = zeta/eta + mu - 1 - eta
    (log zeta)' = x/zeta - (1+eta)
    (log K)'    = k psi(U)/K - (1+eta)          psi(U) = U e^U - e^U + 1
    (log W)'    = x p^2/W - (1+eta)
    (log L)'    = k phi(U)/L - gamma_E           L = Lambda_E,  phi(U) = e^U - 1 - U

Several orbits of W^u(E) (start radii geometrically spaced inside one linear growth factor) are used
so the spirals fill the disc even when the growth per revolution is large.
"""
from __future__ import annotations
import sys, time
import numpy as np
from scipy.integrate import solve_ivp

ROOT = '/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03'
LANE = ROOT + '/work/lanes/dulac'


def phi(U):
    """e^U - 1 - U, accurate for small U."""
    U = np.asarray(U, dtype=float)
    small = np.abs(U) < 1e-3
    out = np.empty_like(U)
    Us = U[small]
    out[small] = Us**2 / 2 * (1 + Us / 3 * (1 + Us / 4 * (1 + Us / 5)))
    Ub = U[~small]
    out[~small] = np.expm1(Ub) - Ub
    return out


def psi(U):
    """U e^U - e^U + 1, accurate for small U."""
    U = np.asarray(U, dtype=float)
    small = np.abs(U) < 1e-3
    out = np.empty_like(U)
    Us = U[small]
    out[small] = Us**2 / 2 * (1 + 2 * Us / 3 * (1 + 3 * Us / 8 * (1 + 4 * Us / 15)))
    Ub = U[~small]
    out[~small] = Ub * np.expm1(Ub) - (np.expm1(Ub) - Ub)
    return out


def _phi1(U):
    if abs(U) < 1e-3:
        return U * U / 2 * (1 + U / 3 * (1 + U / 4 * (1 + U / 5)))
    return np.expm1(U) - U


def _psi1(U):
    if abs(U) < 1e-3:
        return U * U / 2 * (1 + 2 * U / 3 * (1 + 3 * U / 8 * (1 + 4 * U / 15)))
    em1 = np.expm1(U)
    return U * em1 - (em1 - U)


def eqdata(mu, k):
    roots = np.roots([1.0, 1.0 + mu, mu, k])
    gammaE = -min(r.real for r in roots if abs(r.imag) < 1e-10)
    # polish the real root
    for _ in range(3):
        g = gammaE
        gammaE = g - (g * (g - 1) * (g - mu) - k) / ((g - 1) * (g - mu) + g * (g - mu) + g * (g - 1))
    etaE = gammaE - 1.0
    zetaE = etaE * (etaE + 1.0 - mu)
    sigma = 0.5 * (gammaE - 1 - mu)
    omega = np.sqrt(zetaE - sigma**2)
    return dict(gammaE=gammaE, etaE=etaE, zetaE=zetaE, sigma=sigma, omega=omega)


def rhs(mu, k, gammaE):
    lam = k / mu

    def f(t, s):
        U, d, e, le, lz, lK, lW, lL = s
        p = -lam * np.expm1(e)
        eta = np.exp(le)
        x = k * np.exp(U)
        return [p, np.expm1(U - d), mu * np.expm1(d - e),
                np.exp(lz - le) + mu - 1 - eta,
                k * np.exp(U - lz) - 1 - eta,
                k * _psi1(U) * np.exp(-lK) - 1 - eta,
                x * p * p * np.exp(-lW) - 1 - eta,
                k * _phi1(U) * np.exp(-lL) - gammaE]
    return f


def sample(mu, k, r_keep=1e-7, pts_per_rev=360, max_orbits=16, extra_efolds=14.0, thin=0.12, thin_outer=0.03,
           verbose=True, n_orb=None):
    ed = eqdata(mu, k)
    sig, om, gE = ed['sigma'], ed['omega'], ed['gammaE']
    lam = k / mu
    TE = 2 * np.pi / om
    growth = sig * TE                       # log-growth per revolution in the linear regime
    if n_orb is None:
        n_orb = int(min(max_orbits, max(1, np.ceil(growth / 0.12))))
    t_trans = 45.0 / gE + 2 * TE             # transient for the pull-backs
    keys = ('U', 'd', 'e', 'leta', 'lzeta', 'lK', 'lW', 'lL')
    out = {key: [] for key in ('t', 'orb') + keys}
    f = rhs(mu, k, gE)
    sys.path.insert(0, ROOT + '/work')
    import whlib
    T_cyc = whlib.find_cycle(mu, k)['T']
    for j in range(n_orb):
        r_start = r_keep * np.exp(growth * j / n_orb) * np.exp(-sig * t_trans)
        U0 = r_start
        p0 = sig * r_start
        H0 = ed['zetaE'] * U0 - ed['etaE'] * p0          # unstable plane: Lambda_E = 0 to first order
        s0 = np.array([U0, np.log1p(H0 / k), np.log1p(-p0 / lam), np.log(ed['etaE']), np.log(ed['zetaE']),
                       2 * np.log(r_start), 2 * np.log(r_start), 2 * np.log(r_start)])
        t_end = t_trans + np.log(1.0 / r_keep) / sig + max(extra_efolds / sig, 6.0 * T_cyc)
        n_rev = t_end / TE
        nt = int(n_rev * pts_per_rev) + 2
        tt = np.linspace(0.0, t_end, nt)
        a = 1e-13 * r_start
        atol = np.array([a, a, a, 1e-12, 1e-12, 1e-12, 1e-12, 1e-12])
        sol = solve_ivp(f, [0, t_end], s0, method='DOP853', rtol=1e-11, atol=atol, t_eval=tt)
        if sol.status != 0:
            raise RuntimeError('integration failed: ' + sol.message)
        keep = sol.t >= t_trans
        t = sol.t[keep]
        Y = sol.y[:, keep]
        pp = -lam * np.expm1(Y[2])
        # revolutions are delimited by crest crossings (p: + -> -)
        cr = np.flatnonzero((pp[:-1] > 0) & (pp[1:] <= 0)) + 1
        rev = np.searchsorted(cr, np.arange(len(t)), side='right')     # revolution index of each sample
        nrev = rev.max() + 1
        kept = np.zeros(nrev, dtype=bool)
        last_logr = -np.inf
        for r in range(nrev):
            idx = cr[r - 1] if r > 0 else 0
            amp = abs(Y[0, idx]) + 1e-300
            lr = np.log(amp)
            thr = thin if amp < 0.02 else thin_outer
            if abs(lr - last_logr) >= thr or (j == 0 and r >= nrev - 3):
                kept[r] = True
                last_logr = lr
        kept[0] = False                      # partial first revolution
        m = kept[rev]
        out['t'].append(t[m]); out['orb'].append(np.full(m.sum(), j))
        for i, key in enumerate(keys):
            out[key].append(Y[i, m])
        if verbose:
            print(f'  orbit {j+1}/{n_orb}: t_end={t_end:.1f} revs={n_rev:.0f} kept_pts={m.sum()} nfev={sol.nfev}', flush=True)
    res = {key: np.concatenate(v) for key, v in out.items()}
    res.update(mu=mu, k=k, rho=k / (mu * (1 + mu)), T_cyc=T_cyc, **ed)
    return res


def unpack(d):
    """Derived quantities from a saved sample dict/npz."""
    mu = float(d['mu']); k = float(d['k']); lam = k / mu
    gE = float(d['gammaE']); etaE = float(d['etaE']); zetaE = float(d['zetaE'])
    sig = float(d['sigma']); om = float(d['omega'])
    U = d['U']; dd = d['d']; e = d['e']
    p = -lam * np.expm1(e)
    H = k * np.expm1(dd)
    y = k * np.exp(dd)
    x = k * np.exp(U)
    eta = np.exp(d['leta']); zeta = np.exp(d['lzeta'])
    K = np.exp(d['lK']); W = np.exp(d['lW']); L = np.exp(d['lL'])
    pd = k * (np.exp(e) - np.exp(dd))     # p' = k(Y - Z); avoids the cancellation in -mu*p - H when p ~ lambda
    return dict(mu=mu, k=k, lam=lam, A=1 + mu, gE=gE, etaE=etaE, zetaE=zetaE, sig=sig, om=om,
                rho=k / (mu * (1 + mu)), t=d['t'], orb=d['orb'],
                U=U, p=p, H=H, y=y, x=x, eta=eta, zeta=zeta, K=K, W=W, L=L, pd=pd,
                lamp=lam * np.exp(e), phi=phi(U), psi=psi(U),
                etad=zeta + (mu - 1) * eta - eta**2, zetad=x - (1 + eta) * zeta)


SETS = [(0.05, 1.02), (0.05, 1.5), (0.05, 5.0), (0.2, 1.05), (0.2, 2.0), (0.5, 1.2), (0.5, 3.0),
        (1.0, 1.02), (1.0, 1.5), (1.0, 5.0)]


def fname(mu, rho, tag='data'):
    return f'{LANE}/{tag}_mu{mu}_rho{rho}.npz'


if __name__ == '__main__':
    sets = SETS
    if len(sys.argv) > 1:
        sets = [(float(sys.argv[1]), float(sys.argv[2]))]
    for mu, rho in sets:
        k = rho * mu * (1 + mu)
        t0 = time.time()
        print(f'mu={mu} rho={rho} k={k:.6g}', flush=True)
        d = sample(mu, k, verbose=False)
        np.savez_compressed(fname(mu, rho), **d)
        u = unpack(d)
        # consistency checks: K = zeta U - eta p - H ; L = H - zetaE U + etaE p ; orthogonality
        Kdir = u['zeta'] * u['U'] - u['eta'] * u['p'] - u['H']
        Ldir = u['H'] - u['zetaE'] * u['U'] + u['etaE'] * u['p']
        Wdir = u['zetad'] * u['p'] - u['etad'] * u['pd']
        orth = u['zeta'] * u['p'] - u['eta'] * u['pd'] - (u['x'] - u['y'])
        r2 = u['U']**2 + u['p']**2
        big = r2 > 1e-6
        print(f'  n={len(u["t"])} Umin={u["U"].min():.3f} Umax={u["U"].max():.3f} '
              f'relerr(K)={np.max(np.abs(Kdir-u["K"])[big]/u["K"][big]):.1e} '
              f'relerr(L)={np.max(np.abs(Ldir-u["L"])[big]/u["L"][big]):.1e} '
              f'relerr(W)={np.max(np.abs(Wdir-u["W"])[big]/u["W"][big]):.1e} '
              f'orth={np.max(np.abs(orth)[big]/np.sqrt(r2[big])):.1e} [{time.time()-t0:.1f}s]', flush=True)
