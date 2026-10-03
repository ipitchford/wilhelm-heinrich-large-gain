#!/usr/bin/env python3
"""Job 2(c), state-explicit version (no eta, zeta in the inequality).  For a 3D vector field nu(U,p,yh) put
   omega = conj(f ^ nu),  conj(W) = (W12, -W13, -W23),   kappa = Dnu f - J nu.
Then  D omega - C omega = conj(f ^ kappa)  with C the compound matrix of v2.0 (4.6).  If omega > 0 and
conj(f ^ kappa) >= 0 (not identically 0) along every cycle, then every cycle has lambda_c < 0  (brief item (e)).
LP: nu in the span of (b_j e_i), b_j monomials of degree <= 2 in (U, p, yh, E=e^U-1); maximise t subject to
   omega_i >= |f||r|,   conj(f ^ kappa)_i >= t |f|^2     at sampled attractor points with r >= rmin.
t* >= 0 is a numerical certificate candidate on the sample; t* < 0 means no certificate in this span."""
import json
import numpy as np
from scipy.optimize import linprog
from job2_lib import sample
from whlib import kH, find_cycle

def basis3(U, p, yh, E):
    one = np.ones_like(U); zero = np.zeros_like(U)
    prim = {'U': (U, (one, zero, zero)), 'p': (p, (zero, one, zero)), 'yh': (yh, (zero, zero, one)), 'E': (E, (E+1, zero, zero))}
    out = dict(prim); names = list(prim)
    for i in range(4):
        for j in range(i, 4):
            a, b = names[i], names[j]; fa, ga = prim[a]; fb, gb = prim[b]
            out[a+'*'+b] = (fa*fb, tuple(ga[l]*fb+fa*gb[l] for l in range(3)))
    return out

def conj_wedge(f, g):
    W12 = f[0]*g[1]-f[1]*g[0]; W13 = f[0]*g[2]-f[2]*g[0]; W23 = f[1]*g[2]-f[2]*g[1]
    return [W12, -W13, -W23]

def run(mu, rho, rfrac, deg2=True):
    k = rho*kH(mu)
    S = sample(mu, k, np.linspace(0, 2*np.pi, 8, endpoint=False), n=2500)
    cyc = find_cycle(mu, k)
    U, p, y, x = S['U'], S['p'], S['y'], S['x']
    yh = y-k; E = np.exp(U)-1
    f = [p, k-mu*p-y, x-y]
    B = basis3(U, p, yh, E)
    if not deg2: B = {nm: B[nm] for nm in ('U', 'p', 'yh', 'E')}
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
    idx = np.nonzero(m)[0]
    rng = np.random.default_rng(0)
    if len(idx) > 4000: idx = rng.choice(idx, 4000, replace=False)
    nv = len(names)
    A = []; bb = []
    for l in range(3):
        Om = np.column_stack(om_cols[l])[idx]; Ka = np.column_stack(ka_cols[l])[idx]
        A.append(np.hstack([-Om, np.zeros((len(idx), 1))])); bb.append(-(sp*r)[idx])
        A.append(np.hstack([-Ka, (sp**2)[idx][:, None]])); bb.append(np.zeros(len(idx)))
    cobj = np.zeros(nv+1); cobj[-1] = -1
    lp = linprog(cobj, A_ub=np.vstack(A), b_ub=np.concatenate(bb), bounds=[(-500, 500)]*nv+[(-10, 10)], method='highs')
    res = dict(mu=mu, rho=rho, rfrac=rfrac, deg2=deg2, n=int(len(idx)), status=int(lp.status))
    if lp.status == 0:
        res['t'] = float(lp.x[-1])
        c = lp.x[:-1]
        viol = 0.0
        for l in range(3):
            Ka = np.column_stack(ka_cols[l])[m]@c; Om = np.column_stack(om_cols[l])[m]@c
            viol = max(viol, float(np.mean(Ka < 0)), float(np.mean(Om <= 0)))
        res['validation_max_frac_violation'] = viol
        res['coeffs'] = {nm: float(v) for nm, v in zip(names, c) if abs(v) > 1e-7}
    return res

if __name__ == '__main__':
    out = []
    for mu, rho in [(1.0, 1.1), (1.0, 2.0), (0.3, 1.5), (3.0, 1.3)]:
        for rfrac, deg2 in ((0.3, False), (0.05, True), (0.3, True)):
            r = run(mu, rho, rfrac, deg2); out.append(r)
            print(f"mu={mu} rho={rho} r>{rfrac}crest deg2={deg2}: status={r['status']} t*={r.get('t')} validation frac viol={r.get('validation_max_frac_violation')}", flush=True)
    json.dump(out, open('job2_c3.json', 'w'), indent=1)
