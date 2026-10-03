#!/usr/bin/env python3
"""Sparse (L1-minimal) state-explicit compound certificates on the attractor sample, and their structure.
   minimise sum |c_j|  s.t.  omega_l >= |f||r|,  conj(f^kappa)_l >= eps*|f|^2  (l=1,2,3) on sampled points r>=rmin."""
import json, sys
import numpy as np
from scipy.optimize import linprog
from job2_lib import sample
from job2_c3 import basis3, conj_wedge
from whlib import kH, find_cycle

def build(mu, rho, rfrac, nsub=4000, seed=0):
    k = rho*kH(mu)
    S = sample(mu, k, np.linspace(0, 2*np.pi, 8, endpoint=False), n=2500)
    cyc = find_cycle(mu, k)
    U, p, y, x = S['U'], S['p'], S['y'], S['x']
    yh = y-k; E = np.exp(U)-1
    f = [p, k-mu*p-y, x-y]
    B = basis3(U, p, yh, E)
    zero = np.zeros_like(U)
    om_cols = [[], [], []]; ka_cols = [[], [], []]; names = []
    for nm, (b, gb) in B.items():
        Db_f = gb[0]*f[0]+gb[1]*f[1]+gb[2]*f[2]
        for i in range(3):
            nu = [zero, zero, zero]; nu[i] = b
            Dnu_f = [zero, zero, zero]; Dnu_f[i] = Db_f
            Jnu = [nu[1], -mu*nu[1]-nu[2], x*nu[0]-nu[2]]
            kappa = [Dnu_f[l]-Jnu[l] for l in range(3)]
            om = conj_wedge(f, nu); ka = conj_wedge(f, kappa)
            for l in range(3):
                om_cols[l].append(om[l]); ka_cols[l].append(ka[l])
            names.append(f'{nm}*e{i+1}')
    r = np.sqrt(U*U+p*p+(yh/k)**2); sp = np.sqrt(f[0]**2+f[1]**2+f[2]**2)
    m = np.sqrt(U*U+p*p) > rfrac*cyc['U0']
    Om = [np.column_stack(c) for c in om_cols]; Ka = [np.column_stack(c) for c in ka_cols]
    return dict(names=names, Om=Om, Ka=Ka, r=r, sp=sp, m=m, k=k)

def sparse(mu, rho, rfrac, eps=0.0, nsub=4000):
    Bd = build(mu, rho, rfrac)
    names = Bd['names']; nv = len(names)
    idx = np.nonzero(Bd['m'])[0]
    rng = np.random.default_rng(1)
    if len(idx) > nsub: idx = rng.choice(idx, nsub, replace=False)
    # variables c = cp - cm, cp,cm >= 0
    A = []; bb = []
    for l in range(3):
        Om = Bd['Om'][l][idx]; Ka = Bd['Ka'][l][idx]
        A.append(np.hstack([-Om, Om])); bb.append(-(Bd['sp']*Bd['r'])[idx])
        A.append(np.hstack([-Ka, Ka])); bb.append(-eps*(Bd['sp']**2)[idx])
    lp = linprog(np.ones(2*nv), A_ub=np.vstack(A), b_ub=np.concatenate(bb), bounds=[(0, None)]*(2*nv), method='highs')
    if lp.status != 0:
        return dict(mu=mu, rho=rho, status=int(lp.status))
    c = lp.x[:nv]-lp.x[nv:]
    viol = 0.0
    for l in range(3):
        viol = max(viol, float(np.mean(Bd['Ka'][l][Bd['m']]@c < -1e-9)), float(np.mean(Bd['Om'][l][Bd['m']]@c <= 0)))
    co = {nm: round(float(v), 4) for nm, v in zip(names, c) if abs(v) > 1e-6}
    return dict(mu=mu, rho=rho, k=Bd['k'], status=0, l1=float(np.abs(c).sum()), nnz=len(co), validation_frac_viol=viol, coeffs=co)

if __name__ == '__main__':
    out = []
    for mu, rho in [(1.0, 1.1), (1.0, 1.5), (1.0, 2.0), (0.3, 1.5), (3.0, 1.3)]:
        r = sparse(mu, rho, 0.05); out.append(r)
        print(json.dumps(r), flush=True)
    json.dump(out, open('job2_c4.json', 'w'), indent=1)
