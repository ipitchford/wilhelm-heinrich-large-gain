#!/usr/bin/env python3
"""Scaling survey of the computed cycle for large rho (floating point).
For each (mu, rho): T, T_plus (U>0), tau_c = |{x>c0}| for several c0, U_min, U_max, p_max, p_min,
int eta over crest / trough, lambda_c."""
import sys, numpy as np
sys.path.insert(0, '/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03/work')
from whlib import *
from scipy.integrate import solve_ivp

def profile(mu, k, cyc, n=200001):
    T = cyc['T']
    ed = equilibrium_data(mu, k)
    def f(t, s):
        U, p, y, eta, zeta = s
        x = k*np.exp(U)
        return [p, k-mu*p-y, x-y, zeta+(mu-1)*eta-eta*eta, x-(1+eta)*zeta]
    s = np.array([cyc['U0'], 0.0, cyc['y0'], ed['etaE'], ed['zetaE']])
    for _ in range(3):
        sol = solve_ivp(f, [0, T], s, method='LSODA', rtol=1e-10, atol=1e-13)
        s = sol.y[:, -1].copy(); s[0], s[1], s[2] = cyc['U0'], 0.0, cyc['y0']
    sol = solve_ivp(f, [0, T], s, method='LSODA', rtol=1e-10, atol=1e-13, dense_output=True)
    return sol

def survey(mu, rho, guess=None):
    k = rho*kH(mu)
    c = find_cycle(mu, k, guess=guess)
    T = c['T']
    sol = profile(mu, k, c)
    # nonuniform grid: use solver's own steps refined
    t = np.unique(np.concatenate([sol.t, np.linspace(0, T, 20001)]))
    tt = []
    for a, b in zip(t[:-1], t[1:]):
        tt.extend(np.linspace(a, b, 5)[:-1])
    tt.append(T); t = np.array(tt)
    U, p, y, eta, zeta = sol.sol(t)
    x = k*np.exp(U)
    dt = np.diff(t)
    def integ(f, mask=None):
        g = 0.5*(f[1:]+f[:-1])
        if mask is not None:
            m = mask[1:] & mask[:-1]
            return float(np.sum(g[m]*dt[m]))
        return float(np.sum(g*dt))
    out = dict(mu=mu, rho=rho, k=k, lam=k/mu, T=T, Umin=U.min(), Umax=U.max(), pmax=p.max(), pmin=p.min(),
               ymax=y.max(), etamax=eta.max(), mean_eta=integ(eta)/T, lam_c=integ(eta)/T-mu, m_t=c['m_t'],
               Tplus=integ(np.ones_like(t), U > 0))
    for c0 in (kH(mu), 0.5*mu, 0.125*mu*mu):
        mask = x > c0
        out[f'tau_{c0:.3g}'] = integ(np.ones_like(t), mask)
        out[f'Ieta_crest_{c0:.3g}'] = integ(eta, mask)
        out[f'Ieta_trough_{c0:.3g}'] = integ(eta, ~mask)
        out[f'Ix13_{c0:.3g}'] = integ(x**(1/3), mask)
    out['guess'] = (c['U0'], c['y0'])
    return out

if __name__ == '__main__':
    for mu in (1.0, 0.3, 0.1):
        g = None
        for rho in (1.5, 2, 3, 5, 8, 12, 20, 40, 80):
            try:
                o = survey(mu, rho, g)
            except Exception as e:
                print('FAIL', mu, rho, e); break
            g = o.pop('guess')
            print(' '.join(f'{a}={b:.4g}' for a, b in o.items()), flush=True)
