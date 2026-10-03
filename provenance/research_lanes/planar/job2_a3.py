#!/usr/bin/env python3
"""Job 2(a) supplement: in which coordinates on the crest half-axis is the return map concave / sublinear?
Coordinates c(u): u,  e=G~(u) (energy),  yh=h(u,0)-k,  xh=k(e^u-1),  w=sqrt(2(e^u-1-u)).
For coordinate c, derivative of the return map is  R_c' = [c'(u1)/c'(u)] * [g~(u)/g~(u1)] * exp(iota(u)), u1=R1(u),
with g~(u)=h(u,0)-k, e'(u)=g~(u), yh'(u)=zeta(u,0), xh' = k e^u, w' = (e^u-1)/w.
Report: R_c' decreasing in u (concavity)?  R_c(c)/c decreasing (sublinear)?"""
import sys, json
import numpy as np
from scipy.integrate import solve_ivp
sys.path.insert(0, '/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03/work')
from whlib import equilibrium_data, kH

def ext(mu, k):
    def f(t, s):
        U, p, y, eta, zeta, I = s
        x = k*np.exp(U)
        return [p, k-mu*p-y, x-y, zeta+(mu-1)*eta-eta*eta, x-(1+eta)*zeta, eta-mu]
    return f
def ev_c(t, s): return s[1]
ev_c.direction = -1

def run(mu, rho, nth=12):
    k = rho*kH(mu); ed = equilibrium_data(mu, k)
    tend = (np.log(1e6)+16)/ed['lam_u'].real
    rows = []
    for th in np.linspace(0, 2*np.pi, nth, endpoint=False):
        v = np.real(np.exp(1j*th)*ed['vec_u']); z0 = np.array([0, 0, k])+1e-6*v/np.linalg.norm(v)
        sol = solve_ivp(ext(mu, k), [0, tend], list(z0)+[ed['etaE'], ed['zetaE'], 0.0], method='DOP853', rtol=1e-12, atol=1e-15, events=ev_c)
        Y = sol.y_events[0]
        for a, b in zip(Y[:-1], Y[1:]):
            if a[0] > 2e-3:
                rows.append((a[0], b[0], b[5]-a[5], a[2]-k, b[2]-k, a[4], b[4]))
    R = np.array(rows); R = R[np.argsort(R[:, 0])]
    keep = [0]
    for i in range(1, len(R)):
        if R[i, 0]-R[keep[-1], 0] > 1e-5: keep.append(i)
    R = R[keep]
    u, u1, iota, gt, gt1, z, z1 = R.T
    base = gt/gt1*np.exp(iota)                       # derivative in coordinate u
    w = lambda s: np.sqrt(2*(np.exp(s)-1-s))
    cder = {'u': (np.ones_like(u), np.ones_like(u)), 'energy': (gt, gt1), 'y-k': (z, z1),
            'x-k': (k*np.exp(u), k*np.exp(u1)), 'w': ((np.exp(u)-1)/w(u), (np.exp(u1)-1)/w(u1))}
    cval = {'u': (u, u1), 'y-k': (gt, gt1), 'x-k': (k*(np.exp(u)-1), k*(np.exp(u1)-1)), 'w': (w(u), w(u1))}
    res = dict(mu=mu, rho=rho, n=int(len(u)))
    for name, (c0, c1) in cder.items():
        d = base*c1/c0
        res[name+'_concave'] = bool(np.all(np.diff(d) < 0)); res[name+'_max_increase_of_derivative'] = float(np.diff(d).max())
    for name, (c0, c1) in cval.items():
        q = c1/c0
        res[name+'_sublinear'] = bool(np.all(np.diff(q) < 0))
    return res

if __name__ == '__main__':
    out = []
    for mu, rho in [(1.0, 1.02), (1.0, 1.1), (1.0, 1.5), (1.0, 2.0), (1.0, 3.0), (0.3, 1.1), (0.3, 1.5), (0.3, 2.5), (0.1, 1.2), (0.1, 2.0),
                    (3.0, 1.3), (3.0, 2.0), (10.0, 1.5), (0.03, 1.5)]:
        r = run(mu, rho); out.append(r)
        print(f"mu={mu} rho={rho} n={r['n']}: concave in " + ' '.join(f"{nm}={r[nm+'_concave']}" for nm in ('u', 'energy', 'y-k', 'x-k', 'w'))
              + ' | sublinear in ' + ' '.join(f"{nm}={r[nm+'_sublinear']}" for nm in ('u', 'y-k', 'x-k', 'w')), flush=True)
    json.dump(out, open('job2_a3.json', 'w'), indent=1)
