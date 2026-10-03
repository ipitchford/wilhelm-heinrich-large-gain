#!/usr/bin/env python3
"""Job 2(c), LP search for a transversal vector field nu on the attractor disc with one-signed rotation of the flow:
   X ^ nu >= |X||q|   and   X ^ [nu,X] <= tau |X|^2   on the sampled disc (r >= rmin),  [nu,X] = DX nu - Dnu X.
Then B = 1/(X ^ nu) is a Dulac function: div(B X) = B * X^[nu,X]/(X^nu).  tau* < 0 would be a numerical certificate
candidate.  nu ranges over the linear span of fields (b_j,0),(0,b_j), b_j monomials of degree <= 2 in
(u, v, yh=h-k, E=e^u-1) vanishing at the equilibrium.  X = (v, k - mu v - h), DX = [[0,1],[-zeta, eta-mu]]."""
import json
import numpy as np
from scipy.optimize import linprog
from job2_lib import sample
from whlib import kH, find_cycle

def basis(u, v, yh, E, zeta, eta):
    """list of (b, b_u, b_v) for scalar functions on the disc (total derivatives along the graph y=h(u,v))"""
    one = np.ones_like(u); zero = np.zeros_like(u)
    prim = {'u': (u, one, zero), 'v': (v, zero, one), 'yh': (yh, zeta, -eta), 'E': (E, E+1, zero)}
    out = dict(prim)
    names = list(prim)
    for i in range(4):
        for j in range(i, 4):
            a, b = names[i], names[j]
            fa, fau, fav = prim[a]; fb, fbu, fbv = prim[b]
            out[a+'*'+b] = (fa*fb, fau*fb+fa*fbu, fav*fb+fa*fbv)
    return out

def run(mu, rho, rfrac):
    k = rho*kH(mu)
    S = sample(mu, k, np.linspace(0, 2*np.pi, 8, endpoint=False), n=2500)
    cyc = find_cycle(mu, k)
    u, v, y, eta, zeta = S['U'], S['p'], S['y'], S['eta'], S['zeta']
    yh = y-k; E = np.exp(u)-1
    X1 = v; X2 = k-mu*v-y
    B = basis(u, v, yh, E, zeta, eta)
    cols_w = []; cols_c = []; names = []
    for nm, (b, bu, bv) in B.items():
        # nu = (b,0):  X^nu = X1*0 - X2*b = -X2 b ;  [nu,X] = DX nu - Dnu X = (0*b + 1*0, -zeta*b) - (bu X1 + bv X2, 0)
        w1 = -X2*b
        c1_1 = -(bu*X1+bv*X2); c1_2 = -zeta*b
        cols_w.append(w1); cols_c.append(X1*c1_2-X2*c1_1); names.append('('+nm+',0)')
        # nu = (0,b):  X^nu = X1 b ;  [nu,X] = DX nu - Dnu X = (b, (eta-mu) b) - (0, bu X1 + bv X2)
        w2 = X1*b
        c2_1 = b; c2_2 = (eta-mu)*b-(bu*X1+bv*X2)
        cols_w.append(w2); cols_c.append(X1*c2_2-X2*c2_1); names.append('(0,'+nm+')')
    Wm = np.column_stack(cols_w); Cm = np.column_stack(cols_c)
    r = np.sqrt(u*u+v*v); sp2 = X1*X1+X2*X2
    m = r > rfrac*cyc['U0']
    idx = np.nonzero(m)[0]
    rng = np.random.default_rng(0)
    if len(idx) > 5000: idx = rng.choice(idx, 5000, replace=False)
    nv = Wm.shape[1]
    # variables: c (nv), tau.   -W c <= -|X||q| ;  C c - tau |X|^2 <= 0
    A1 = np.hstack([-Wm[idx], np.zeros((len(idx), 1))]); b1 = -np.sqrt(sp2[idx])*r[idx]
    A2 = np.hstack([Cm[idx], -sp2[idx][:, None]]); b2 = np.zeros(len(idx))
    cobj = np.zeros(nv+1); cobj[-1] = 1
    lp = linprog(cobj, A_ub=np.vstack([A1, A2]), b_ub=np.concatenate([b1, b2]), bounds=[(-200, 200)]*nv+[(-10, 10)], method='highs')
    res = dict(mu=mu, rho=rho, rfrac=rfrac, n=int(len(idx)), status=int(lp.status))
    if lp.status == 0:
        c = lp.x[:-1]; res['tau'] = float(lp.x[-1])
        beta = (Cm@c)/(Wm@c)
        res['validation_min_wedge_ratio'] = float(np.min((Wm@c)[m]/(np.sqrt(sp2[m])*r[m])))
        res['validation_max_beta'] = float(beta[m].max()); res['validation_frac_beta_pos'] = float(np.mean(beta[m] > 0))
        res['coeffs'] = {nm: float(x) for nm, x in zip(names, c) if abs(x) > 1e-6}
    return res

if __name__ == '__main__':
    out = []
    for mu, rho in [(1.0, 1.1), (1.0, 2.0), (0.3, 1.5), (3.0, 1.3)]:
        for rfrac in (0.05, 0.3):
            r = run(mu, rho, rfrac); out.append(r)
            print(f"mu={mu} rho={rho} r>{rfrac}crest: status={r['status']} tau*={r.get('tau')} valid max beta={r.get('validation_max_beta')} frac beta>0={r.get('validation_frac_beta_pos')} min wedge ratio={r.get('validation_min_wedge_ratio')}", flush=True)
    json.dump(out, open('job2_c2.json', 'w'), indent=1)
