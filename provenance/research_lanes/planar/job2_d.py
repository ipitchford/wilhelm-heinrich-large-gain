#!/usr/bin/env python3
"""Job 2(d): rotated-vector-field structure in k (numerical probes).
 D1  Are cycles at different k (fixed mu) disjoint/nested in a fixed plane?  planes: (U,p), (V=U+log k,p), (U,s/k)
     with s = p-(y-k)+mu U the Lienard variable (s' = -k(e^U-1)).
 D2  In the plane (U, st=s/k):  U' = Pt_k(U,st), st' = -(e^U-1).  F ^ dF/dk = (e^U-1) dPt_k/dk.  Rotated family
     iff (e^U-1) dPt/dk has one sign.  Estimated by scattered interpolation of p over (U,st) on the attractor
     discs at k and 1.03k.
 D3  cross-parameter ordering (brief item (h)): in V coordinates h_{k2} > H_{k1} on the common domain (checked on
     the section p=0 with Hermite curves) and the size of d h_k/dk there (is it < 1 ?).
"""
import sys, json
import numpy as np
from scipy.interpolate import LinearNDInterpolator, CubicHermiteSpline
from scipy.integrate import solve_ivp
sys.path.insert(0, '/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03/work')
from whlib import find_cycle, cycle_profile, equilibrium_data, kH, unstable_manifold_orbit

def n_cross(A, B):
    a0, a1 = A[:-1], A[1:]; b0, b1 = B[:-1], B[1:]
    da = a1-a0; db = b1-b0; cnt = 0
    for i in range(len(a0)):
        r = da[i]; qp = b0-a0[i]
        den = r[0]*db[:, 1]-r[1]*db[:, 0]
        with np.errstate(divide='ignore', invalid='ignore'):
            t = (qp[:, 0]*db[:, 1]-qp[:, 1]*db[:, 0])/den
            u = (qp[:, 0]*r[1]-qp[:, 1]*r[0])/den
        cnt += int(np.sum((den != 0) & (t >= 0) & (t < 1) & (u >= 0) & (u < 1)))
    return cnt

def D1(mu, rhos):
    cyc = None; curves = []
    for rho in rhos:
        k = rho*kH(mu)
        c = find_cycle(mu, k, guess=None if cyc is None else (cyc['U0'], cyc['y0']*k/cyc['k']))
        cyc = c
        pr = cycle_profile(mu, k, c, n=1601)
        s = pr['p']-(pr['y']-k)+mu*pr['U']
        curves.append(dict(rho=rho, Up=np.c_[pr['U'], pr['p']], Vp=np.c_[pr['U']+np.log(k), pr['p']],
                           Us=np.c_[pr['U'], s/k], Upl=np.c_[pr['U'], pr['p']/(k/mu)]))
    res = {}
    for key in ('Up', 'Vp', 'Us', 'Upl'):
        res[key] = [n_cross(curves[i][key], curves[i+1][key]) for i in range(len(curves)-1)]
    return res

def cloud(mu, k, nth=10, n=2500):
    P = []
    for th in np.linspace(0, 2*np.pi, nth, endpoint=False):
        o = unstable_manifold_orbit(mu, k, r0=1e-6, theta=th, n=n)
        m = np.sqrt(o['U']**2+o['p']**2) > 1e-3
        P.append(np.c_[o['U'][m], o['p'][m], o['y'][m]])
    return np.vstack(P)

def D2(mu, rho, dk=0.03):
    k1 = rho*kH(mu); k2 = k1*(1+dk)
    out = {}
    C1 = cloud(mu, k1); C2 = cloud(mu, k2)
    st1 = (C1[:, 1]-(C1[:, 2]-k1)+mu*C1[:, 0])/k1; st2 = (C2[:, 1]-(C2[:, 2]-k2)+mu*C2[:, 0])/k2
    I1 = LinearNDInterpolator(np.c_[C1[:, 0], st1], C1[:, 1]); I2 = LinearNDInterpolator(np.c_[C2[:, 0], st2], C2[:, 1])
    # evaluate on the points of cloud 1 (subsample) -> p2 - p1 at the same (U,st)
    idx = np.random.default_rng(0).choice(len(C1), 4000, replace=False)
    Q = np.c_[C1[idx, 0], st1[idx]]
    dP = (I2(Q)-C1[idx, 1])/(k2-k1)
    ok = np.isfinite(dP) & (np.abs(Q[:, 0]) > 0.02)
    rot = (np.exp(Q[ok, 0])-1)*dP[ok]
    out['n'] = int(ok.sum()); out['frac_rot_positive'] = float(np.mean(rot > 0))
    out['rot_quantiles'] = [float(np.quantile(rot, a)) for a in (0.01, 0.25, 0.5, 0.75, 0.99)]
    # interpolation noise estimate: self-consistency of cloud 1 against its own interpolant built from half the spirals
    return out

if __name__ == '__main__':
    res = {}
    for mu in (0.3, 1.0, 3.0):
        rhos = [1.05, 1.15, 1.3, 1.5, 1.75, 2.0, 2.4, 3.0]
        r = D1(mu, rhos); res[f'D1_mu{mu}'] = r
        print(f"D1 mu={mu}: crossings between consecutive cycles (rho={rhos}):")
        for key, v in r.items(): print(f"     plane {key}: {v}")
    for mu, rho in [(1.0, 1.3), (1.0, 2.0), (0.3, 1.5), (3.0, 1.5)]:
        r = D2(mu, rho); res[f'D2_mu{mu}_rho{rho}'] = r
        print(f"D2 mu={mu} rho={rho}: n={r['n']} frac[(e^U-1) dPt/dk > 0] = {r['frac_rot_positive']:.3f}  quantiles {['%.3g' % a for a in r['rot_quantiles']]}", flush=True)
    json.dump(res, open('job2_d.json', 'w'), indent=1)
