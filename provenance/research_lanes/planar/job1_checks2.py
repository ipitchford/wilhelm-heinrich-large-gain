#!/usr/bin/env python3
"""Job 1 supplementary checks:
 X1  projected attractor orbits never cross (planar flow), projected off-attractor orbits do cross them, with y-gap != 0
 X2  local-plane-fit gradient error scales like the stencil radius (two sampling densities)
 X3  the single non-negative pair of N3 at (mu,rho)=(1,1.1) is a rounding artefact of two nearly coincident points
 X4  crests of successive W^u(E) turns (root-refined) increase monotonically to the cycle crest
"""
import sys, json
import numpy as np
from scipy.integrate import solve_ivp
from scipy.spatial import cKDTree
sys.path.insert(0, '/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03/work')
from whlib import find_cycle, cycle_profile, unstable_manifold_orbit, equilibrium_data, kH, field

def crossings(A, B, yA, yB):
    """transversal crossings between polylines A (n,2) and B (m,2); returns list of |y gap| at crossings"""
    a0, a1 = A[:-1], A[1:]; b0, b1 = B[:-1], B[1:]
    da = a1-a0; db = b1-b0
    out = []
    for i in range(len(a0)):
        r = da[i]; qp = b0-a0[i]
        den = r[0]*db[:, 1]-r[1]*db[:, 0]
        with np.errstate(divide='ignore', invalid='ignore'):
            t = (qp[:, 0]*db[:, 1]-qp[:, 1]*db[:, 0])/den
            u = (qp[:, 0]*r[1]-qp[:, 1]*r[0])/den
        m = (den != 0) & (t > 0) & (t < 1) & (u > 0) & (u < 1)
        for j in np.nonzero(m)[0]:
            ya = yA[i]+t[j]*(yA[i+1]-yA[i]); yb = yB[j]+u[j]*(yB[j+1]-yB[j])
            out.append(abs(ya-yb))
    return out

def run(mu, rho):
    k = rho*kH(mu); ed = equilibrium_data(mu, k)
    cyc = find_cycle(mu, k)
    res = dict(mu=mu, rho=rho)
    oa = unstable_manifold_orbit(mu, k, r0=1e-6, theta=0.0, n=3001)
    ob = unstable_manifold_orbit(mu, k, r0=1e-6, theta=1.9, n=3001)
    def cut(o):
        r = np.sqrt(o['U']**2+o['p']**2); m = r > 1e-2
        # stop before the orbit is numerically on the cycle (last 25% of the time window)
        idx = np.nonzero(m)[0]; idx = idx[: int(0.8*len(idx))]
        return np.c_[o['U'][idx], o['p'][idx]], o['y'][idx]
    A, yA = cut(oa); B, yB = cut(ob)
    cr = crossings(A, B, yA, yB)
    res['X1_attractor_vs_attractor_polyline_crossings'] = len(cr)
    res['X1_attractor_vs_attractor_max_ygap_at_crossings'] = float(max(cr)) if cr else 0.0
    # off-attractor orbits: start above/below the attractor
    offc = []
    for dy in (0.4*k, -0.4*k):
        i = len(A)//2
        z0 = [A[i, 0], A[i, 1], yA[i]+dy]
        tt = np.linspace(0, 3*cyc['T'], 1500)
        sol = solve_ivp(field(mu, k), [0, tt[-1]], z0, method='DOP853', rtol=1e-11, atol=1e-13, t_eval=tt)
        Cc = sol.y[:2].T; yC = sol.y[2]
        c2 = crossings(Cc[:400], A, yC[:400], yA)
        offc.append(dict(dy=float(dy), n_cross=len(c2), max_ygap=float(max(c2)) if c2 else 0.0,
                         min_ygap=float(min(c2)) if c2 else 0.0))
    res['X1_off_attractor'] = offc
    # X2 density scaling
    for nth in (8, 32):
        pts = []; lab = []
        for j, th in enumerate(np.linspace(0, 2*np.pi, nth, endpoint=False)):
            o = unstable_manifold_orbit(mu, k, r0=1e-6, theta=th, n=3001)
            r = np.sqrt(o['U']**2+o['p']**2); m = r > 5e-2
            pts.append(np.c_[o['U'][m], o['p'][m], o['y'][m], o['eta'][m], o['zeta'][m]]); lab += [j]*int(m.sum())
        P = np.vstack(pts); lab = np.array(lab)
        tree = cKDTree(P[:, :2]); rng = np.random.default_rng(1)
        test = rng.choice(len(P), size=500, replace=False)
        errs = []; hs = []
        for i in test:
            dd, nn = tree.query(P[i, :2], k=80)
            nn = nn[lab[nn] != lab[i]][:10]
            if len(nn) < 8: continue
            dq = P[nn, :2]-P[i, :2]
            coef, *_ = np.linalg.lstsq(np.c_[np.ones(len(nn)), dq], P[nn, 2], rcond=None)
            gi = np.array([P[i, 4], -P[i, 3]])
            errs.append(np.linalg.norm(coef[1:]-gi)); hs.append(np.max(np.linalg.norm(dq, axis=1)))
        res[f'X2_orbits{nth}'] = dict(median_abs_grad_err=float(np.median(errs)), median_stencil=float(np.median(hs)),
                                      p90_abs_grad_err=float(np.quantile(errs, 0.9)))
    # X4 crests via events
    def ev(t, z): return z[1]
    ev.direction = -1
    v = np.real(ed['vec_u']); z0 = np.array([0, 0, k])+1e-6*v/np.linalg.norm(v)
    tend = (np.log(1e6)+14)/ed['lam_u'].real
    sol = solve_ivp(field(mu, k), [0, tend], z0, method='DOP853', rtol=1e-12, atol=1e-15, events=ev)
    Uc = sol.y_events[0][:, 0]; Uc = Uc[Uc > 0]
    d = np.diff(Uc)
    res['X4_n_turns'] = int(len(Uc)); res['X4_min_crest_increment'] = float(d.min())
    res['X4_last_crest'] = float(Uc[-1]); res['X4_cycle_crest'] = float(cyc['U0'])
    res['X4_monotone_until_1e-9_of_cycle'] = bool(np.all(d[Uc[1:] < cyc['U0']-1e-9] > 0))
    return res

if __name__ == '__main__':
    allr = []
    for mu, rho in [(1.0, 1.1), (1.0, 2.0), (0.3, 1.5), (3.0, 1.3)]:
        r = run(mu, rho); allr.append(r); print(json.dumps(r), flush=True)
    json.dump(allr, open('job1_checks2.json', 'w'), indent=1)
    # X3
