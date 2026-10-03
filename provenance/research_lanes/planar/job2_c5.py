#!/usr/bin/env python3
"""Confirmation run for the state-explicit compound certificate LP (job2_c3) on the FULL attractor sample
(all points with planar radius > rfrac*crest, no subsampling), more parameter pairs, plus an out-of-sample test:
certificate fitted on spirals theta_j, validated on the cycle itself and on independent spirals."""
import json, sys
import numpy as np
from scipy.optimize import linprog
from job2_lib import sample
from job2_c3 import basis3, conj_wedge
from whlib import kH, find_cycle, cycle_profile

def design(mu, k, U, p, y):
    x = k*np.exp(U); yh = y-k; E = np.exp(U)-1
    f = [p, k-mu*p-y, x-y]
    B = basis3(U, p, yh, E)
    zero = np.zeros_like(U)
    om_cols = [[], [], []]; ka_cols = [[], [], []]
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
    r = np.sqrt(U*U+p*p+(yh/k)**2); sp = np.sqrt(f[0]**2+f[1]**2+f[2]**2)
    return [np.column_stack(c) for c in om_cols], [np.column_stack(c) for c in ka_cols], r, sp

def run(mu, rho, rfrac=0.05):
    k = rho*kH(mu); cyc = find_cycle(mu, k)
    S = sample(mu, k, np.linspace(0, 2*np.pi, 6, endpoint=False), n=2000)
    m = np.sqrt(S['U']**2+S['p']**2) > rfrac*cyc['U0']
    Om, Ka, r, sp = design(mu, k, S['U'][m], S['p'][m], S['y'][m])
    nv = Om[0].shape[1]
    A = []; bb = []
    for l in range(3):
        A.append(np.hstack([-Om[l], np.zeros((len(r), 1))])); bb.append(-(sp*r))
        A.append(np.hstack([-Ka[l], (sp**2)[:, None]])); bb.append(np.zeros(len(r)))
    cobj = np.zeros(nv+1); cobj[-1] = -1
    lp = linprog(cobj, A_ub=np.vstack(A), b_ub=np.concatenate(bb), bounds=[(-2000, 2000)]*nv+[(-10, 10)], method='highs')
    res = dict(mu=mu, rho=rho, n=int(len(r)), status=int(lp.status), lam_c=cyc['lambda_c'])
    if lp.status != 0: return res
    c = lp.x[:-1]; res['t'] = float(lp.x[-1])
    # out of sample: the cycle itself, and other spirals
    pr = cycle_profile(mu, k, cyc, n=4001)
    Om2, Ka2, r2, sp2 = design(mu, k, pr['U'], pr['p'], pr['y'])
    res['cycle_min_omega_ratio'] = float(min(np.min((Om2[l]@c)/(sp2*r2)) for l in range(3)))
    res['cycle_min_kappa_ratio'] = float(min(np.min((Ka2[l]@c)/(sp2**2)) for l in range(3)))
    S3 = sample(mu, k, np.array([0.5, 2.3, 4.4]), n=2000)
    m3 = np.sqrt(S3['U']**2+S3['p']**2) > rfrac*cyc['U0']
    Om3, Ka3, r3, sp3 = design(mu, k, S3['U'][m3], S3['p'][m3], S3['y'][m3])
    res['heldout_min_omega_ratio'] = float(min(np.min((Om3[l]@c)/(sp3*r3)) for l in range(3)))
    res['heldout_min_kappa_ratio'] = float(min(np.min((Ka3[l]@c)/(sp3**2)) for l in range(3)))
    return res

if __name__ == '__main__':
    out = []
    grid = [(1.0, 1.02), (1.0, 1.1), (1.0, 1.5), (1.0, 2.0), (1.0, 3.0), (1.0, 5.0), (0.3, 1.1), (0.3, 1.5), (0.3, 3.0),
            (0.1, 1.3), (0.1, 2.0), (3.0, 1.3), (3.0, 2.5), (10.0, 1.5)]
    for mu, rho in grid:
        try:
            r = run(mu, rho)
        except Exception as e:
            print(mu, rho, 'FAILED', repr(e), flush=True); continue
        out.append(r)
        print(json.dumps(r), flush=True)
    json.dump(out, open('job2_c5.json', 'w'), indent=1)
