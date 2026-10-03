#!/usr/bin/env python3
"""Scale-invariant version of the state-explicit compound certificate LP.
 nu in span{ b_j e_i }, b_j monomials of degree <= 2 in (U, p, yh, E=e^U-1)   (42 unknowns, unbounded)
 maximise t  s.t.   omega_l >= omega^rad_l   and   conj(f^kappa)_l >= t sqrt(mu) omega^rad_l   (l=1,2,3)
 on attractor points with planar radius > rfrac*crest: dense spirals on W^u(E) AND the cycle itself,
 where omega^rad = conj(f ^ (U,p,yh)) (>0 on the attractor, checked).  Validation on held-out spirals.
 t* > 0 : numerical certificate candidate (lambda_c<0 for any cycle through the sampled region);
 t* <= 0: no certificate in this span with this normalisation."""
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
    rad = conj_wedge(f, [U, p, yh])
    return [np.column_stack(c) for c in om_cols], [np.column_stack(c) for c in ka_cols], rad

def run(mu, rho, rfrac=0.05, nspir=5, n=5000, nsub=5000):
    k = rho*kH(mu); cyc = find_cycle(mu, k)
    S = sample(mu, k, np.linspace(0, 2*np.pi, nspir, endpoint=False), n=n)
    pr = cycle_profile(mu, k, cyc, n=3001)
    m = np.sqrt(S['U']**2+S['p']**2) > rfrac*cyc['U0']
    U = np.concatenate([S['U'][m], pr['U']]); p = np.concatenate([S['p'][m], pr['p']]); y = np.concatenate([S['y'][m], pr['y']])
    Om, Ka, rad = design(mu, k, U, p, y)
    res = dict(mu=mu, rho=rho, n=int(len(U)), lam_c=cyc['lambda_c'], min_rad=[float(a.min()) for a in rad])
    if len(U) > nsub:
        # keep all cycle points (last 3001) and a random subset of spiral points
        rng = np.random.default_rng(0); ns = len(U)-len(pr['U'])
        sel = np.concatenate([rng.choice(ns, nsub-1000, replace=False), ns+np.arange(0, len(pr['U']), 3)])
        U, p, y = U[sel], p[sel], y[sel]; Om = [a[sel] for a in Om]; Ka = [a[sel] for a in Ka]; rad = [a[sel] for a in rad]
        res['n_lp'] = int(len(U))
    if min(res['min_rad']) <= 0:
        res['note'] = 'radial wedge not positive'; return res
    nv = Om[0].shape[1]; s = np.sqrt(mu)
    A = []; bb = []
    for l in range(3):
        A.append(np.hstack([-Om[l]/rad[l][:, None], np.zeros((len(U), 1))])); bb.append(-np.ones(len(U)))
        A.append(np.hstack([-Ka[l]/rad[l][:, None], s*np.ones((len(U), 1))])); bb.append(np.zeros(len(U)))
    cobj = np.zeros(nv+1); cobj[-1] = -1
    lp = linprog(cobj, A_ub=np.vstack(A), b_ub=np.concatenate(bb), bounds=[(None, None)]*nv+[(-5, 5)], method='highs')
    res['status'] = int(lp.status)
    if lp.status != 0: return res
    c = lp.x[:-1]; res['t'] = float(lp.x[-1]); res['max_abs_coeff'] = float(np.abs(c).max())
    S3 = sample(mu, k, np.array([0.5, 2.9]), n=n)
    m3 = np.sqrt(S3['U']**2+S3['p']**2) > rfrac*cyc['U0']
    Om3, Ka3, rad3 = design(mu, k, S3['U'][m3], S3['p'][m3], S3['y'][m3])
    res['heldout_min_omega_over_rad'] = float(min(np.min((Om3[l]@c)/rad3[l]) for l in range(3)))
    res['heldout_min_kappa_over_s_rad'] = float(min(np.min((Ka3[l]@c)/(s*rad3[l])) for l in range(3)))
    res['heldout_frac_kappa_negative'] = float(max(np.mean((Ka3[l]@c) < 0) for l in range(3)))
    return res

if __name__ == '__main__':
    out = []
    grid = [(1.0, 1.02), (1.0, 1.1), (1.0, 1.5), (1.0, 2.0), (1.0, 3.0), (0.3, 1.1), (0.3, 1.5), (0.3, 3.0),
            (0.1, 1.3), (0.1, 2.0), (3.0, 1.3), (3.0, 2.5), (10.0, 1.5)]
    for mu, rho in grid:
        try:
            r = run(mu, rho)
        except Exception as e:
            print(mu, rho, 'FAILED', repr(e), flush=True); continue
        out.append(r)
        print(json.dumps(r), flush=True)
        json.dump(out, open('job2_c6.json', 'w'), indent=1)
    json.dump(out, open('job2_c6.json', 'w'), indent=1)
