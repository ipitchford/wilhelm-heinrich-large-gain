#!/usr/bin/env python3
"""Sampling of the attractor disc W^u(E) with first AND second order data of h:
state (U,p,y), pull-back (eta,zeta), gauges K,W, and a transversal Jacobi field (dU,dp,dy) with the linearised
Riccati pair (deta,dzeta), from which D^2h = [GF Gw][F w]^{-1}, F=(p,p'), GF=(zeta',-eta'), Gw=(dzeta,-deta)."""
import sys
import numpy as np
from scipy.integrate import solve_ivp
sys.path.insert(0, '/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03/work')
from whlib import equilibrium_data, kH, find_cycle

def rhs(mu, k):
    def f(t, s):
        U, p, y, eta, zeta, K, W, dU, dp, dy, de, dz = s
        x = k*np.exp(U)
        return [p, k-mu*p-y, x-y,
                zeta+(mu-1)*eta-eta*eta, x-(1+eta)*zeta,
                -(1+eta)*K+(x*U-x+k), -(1+eta)*W+x*p*p,
                dp, -mu*dp-dy, x*dU-dy,
                dz+(mu-1-2*eta)*de, x*dU-(1+eta)*dz-zeta*de]
    return f

def sample(mu, k, thetas, r0=1e-6, n=3000, rmin=2e-3, extra_time=10.0):
    ed = equilibrium_data(mu, k)
    tend = (np.log(1/r0)+extra_time)/ed['lam_u'].real
    out = []
    for j, th in enumerate(thetas):
        vc = np.exp(1j*th)*ed['vec_u']; nv = np.linalg.norm(np.real(vc))
        z0 = np.array([0, 0, k])+r0*np.real(vc)/nv
        w0 = r0*np.real(1j*vc)/nv
        s0 = np.concatenate((z0, [ed['etaE'], ed['zetaE'], 0, 0], w0, [0, 0]))
        t = np.linspace(0, tend, n)
        sol = solve_ivp(rhs(mu, k), [0, tend], s0, method='DOP853', rtol=1e-11, atol=1e-18, t_eval=t)
        U, p, y, eta, zeta, K, W, dU, dp, dy, de, dz = sol.y
        x = k*np.exp(U); pp = k-mu*p-y
        etap = zeta+(mu-1)*eta-eta**2; zetap = x-(1+eta)*zeta
        det = p*dp-pp*dU
        # G = [GF Gw] [F w]^{-1},  [F w] = [[p,dU],[pp,dp]]
        inv = np.array([[dp, -dU], [-pp, p]])/det
        GF = np.array([zetap, -etap]); Gw = np.array([dz, -de])
        G11 = GF[0]*inv[0, 0]+Gw[0]*inv[1, 0]; G12 = GF[0]*inv[0, 1]+Gw[0]*inv[1, 1]
        G21 = GF[1]*inv[0, 0]+Gw[1]*inv[1, 0]; G22 = GF[1]*inv[0, 1]+Gw[1]*inv[1, 1]
        r = np.sqrt(U**2+p**2)
        m = r > rmin
        d = dict(orbit=np.full(m.sum(), j), t=sol.t[m], U=U[m], p=p[m], y=y[m], x=x[m], pp=pp[m], eta=eta[m], zeta=zeta[m],
                 etap=etap[m], zetap=zetap[m], K=K[m], W=W[m], huu=G11[m], huv=0.5*(G12+G21)[m], hvv=G22[m],
                 asym=(G12-G21)[m], r=r[m])
        out.append(d)
    keys = out[0].keys()
    S = {kk: np.concatenate([o[kk] for o in out]) for kk in keys}
    S['ed'] = ed; S['mu'] = mu; S['k'] = k
    return S
