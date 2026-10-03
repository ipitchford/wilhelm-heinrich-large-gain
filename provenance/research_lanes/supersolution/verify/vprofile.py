#!/usr/bin/env python3
"""Referee's own profile integrator (independent of the lane's base.py).
State in (U, P = lambda - p, y) so that P and y keep relative accuracy in deep troughs; adjoint ratios by forward Riccati,
K and W by their own linear ODEs (and, as a cross-check, from their definitions); all period averages are integrated
as extra ODE components (no grid quadrature)."""
import sys, pickle
import numpy as np
from scipy.integrate import solve_ivp
ROOT = '/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03'
sys.path.insert(0, ROOT + '/work')
from whlib import equilibrium_data, kH
HERE = ROOT + '/work/lanes/supersolution/verify'

def G1(U):
    return U*np.exp(U) - np.expm1(U)

def derived(mu, k, U, P, y, eta, zeta, K, W, cs=()):
    lam = k/mu
    p = lam - P
    x = k*np.exp(U)
    pd = mu*P - y                    # p' = k - mu p - y
    yd = x - y
    pdd = -mu*pd - yd
    H = y - k
    S1 = p*p - U*pd
    S2 = U*yd - p*H
    S3 = p*yd - pd*H
    S4 = pd*pd - p*pdd
    out = dict(p=p, x=x, pd=pd, yd=yd, pdd=pdd, H=H, S1=S1, S2=S2, S3=S3, S4=S4)
    for c in cs:
        out[('sc', c)] = p*p - pd*(U + c)
    return out

def settle(mu, k, cyc, periods=5, rtol=1e-12):
    """converge (eta, zeta, K, W) to the periodic solution along the cycle; returns state at the crest."""
    lam = k/mu; T = cyc['T']; ed = equilibrium_data(mu, k)
    def f(t, s):
        U, P, y, eta, zeta, K, W = s
        x = k*np.exp(U); p = lam - P
        return [p, y - mu*P, x - y, zeta + (mu-1)*eta - eta*eta, x - (1+eta)*zeta,
                -(1+eta)*K + k*G1(U), -(1+eta)*W + x*p*p]
    s = np.array([cyc['U0'], lam, cyc['y0'], ed['etaE'], ed['zetaE'], 1e-3, 1e-3])
    atol = np.array([1e-13, 1e-200, 1e-200, 1e-200, 1e-200, 1e-200, 1e-200])
    hist = []
    for _ in range(periods):
        sol = solve_ivp(f, [0, T], s, method='DOP853', rtol=rtol, atol=atol)
        new = sol.y[:, -1].copy()
        hist.append(float(np.max(np.abs(new[3:]/s[3:] - 1))))
        closure = float(max(abs(new[0]-cyc['U0']), abs(new[1]-lam)/lam, abs(new[2]-cyc['y0'])/cyc['y0']))
        s = new
        s[0], s[1], s[2] = cyc['U0'], lam, cyc['y0']
    return s, hist, closure

def profile(mu, k, cyc, cfracs=(0.5, 0.9), rtol=1e-12, periods=5):
    lam = k/mu; T = cyc['T']; A = 1 + mu; ed = equilibrium_data(mu, k)
    al = 0.5*(ed['etaE'] - mu)
    s0, hist, closure = settle(mu, k, cyc, periods=periods, rtol=rtol)
    # first pass for U_min (needed for the shifted centres)
    def f0(t, s):
        U, P, y = s
        return [lam - P, y - mu*P, k*np.exp(U) - y]
    def ev(t, s):
        return lam - s[1]
    ev.direction = 1
    a0 = solve_ivp(f0, [0, T], s0[:3], method='DOP853', rtol=rtol, atol=[1e-13, 1e-200, 1e-200], events=ev)
    te = [t for t in a0.t_events[0] if 1e-6*T < t < T*(1-1e-6)]
    Umin = float(a0.y_events[0][[i for i, t in enumerate(a0.t_events[0]) if 1e-6*T < t < T*(1-1e-6)][0]][0])
    t_trough = float(te[0])
    cs = tuple(fr*abs(Umin) for fr in cfracs)
    nI = 5 + len(cs)
    def f(t, s):
        U, P, y, eta, zeta, K, W = s[:7]
        x = k*np.exp(U); p = lam - P
        d = derived(mu, k, U, P, y, eta, zeta, K, W, cs)
        LamE = d['H'] - ed['zetaE']*U + ed['etaE']*p
        QE = p*p - 2*al*U*p + ed['zetaE']*U*U
        out = [p, y - mu*P, x - y, zeta + (mu-1)*eta - eta*eta, x - (1+eta)*zeta,
               -(1+eta)*K + k*G1(U), -(1+eta)*W + x*p*p,
               eta, p*K/d['S1'], p*W/d['S4'], d['S2']/d['S1'], (p - al*U)*LamE/QE]
        for c in cs:
            out.append(p*(K + c*zeta)/d[('sc', c)])
        return out
    s = np.concatenate((s0, np.zeros(nI)))
    atol = np.concatenate(([1e-13], np.full(6, 1e-200), np.full(nI, 1e-13)))
    sol = solve_ivp(f, [0, T], s, method='DOP853', rtol=rtol, atol=atol, dense_output=True)
    # sample: all solver steps + refinement between them
    ts = sol.t
    tt = np.unique(np.concatenate([ts] + [ts[:-1] + fr*np.diff(ts) for fr in (0.25, 0.5, 0.75)]))
    Y = sol.sol(tt)
    U, P, y, eta, zeta, K, W = Y[:7]
    d = derived(mu, k, U, P, y, eta, zeta, K, W, cs)
    I = sol.y[7:, -1]
    res = dict(mu=mu, k=k, T=T, t=tt, U=U, P=P, y=y, eta=eta, zeta=zeta, K=K, W=W, cs=cs, Umin=Umin, t_trough=t_trough,
               lam_eta=I[0]/T - mu, lam_K=-I[1]/T, lam_W=-I[2]/T, mean_S2S1=I[3]/T, mean_QE=I[4]/T, alpha=al,
               lam_c=[-I[5+i]/T for i in range(len(cs))], settle_hist=hist, closure=closure,
               per_defect=float(np.max(np.abs(sol.y[3:7, -1]/sol.y[3:7, 0] - 1))), nsteps=len(ts), ed=ed, sol=sol)
    res.update(d)
    # definitions (independent of the K, W ODEs)
    res['K_def'] = zeta*U - eta*d['p'] - d['H']
    etad = zeta + (mu-1)*eta - eta*eta; zetad = d['x'] - (1+eta)*zeta
    res['W_def'] = zetad*d['p'] - etad*d['pd']
    res['O'] = zeta*d['p'] - eta*d['pd'] - d['yd']       # orthogonality defect
    return res

def load_cycles():
    with open(HERE + '/vcycles.pkl', 'rb') as f:
        return pickle.load(f)
