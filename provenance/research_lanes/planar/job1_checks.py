#!/usr/bin/env python3
"""Numerical checks of the attractor-graph theorem (Job 1). Floating point only.

N1 pull-back uniqueness of the positive Riccati solution (different positive initial data collapse)
N2 generalised tangent orthogonality O = zeta p - eta p' - (x-y) = 0 on W^u(E)
N3 strict supporting planes between ALL sampled pairs of attractor points (W^u(E) orbits, cycle, E)
N4 C^1 matching: grad h = (zeta,-eta) by local least-squares plane fits across neighbouring orbits/turns
N5 off-attractor orbits: y - H(U,p) != 0 and decays; on-attractor orbit: y - H = 0 up to sampling error
N6 uniform bounds eta_-,eta_+,zeta_-,zeta_+ (Lemma) and K>0, W>0, cone exclusion
N8 spiral turns of W^u(E) are nested and accumulate on the cycle (crest sequence monotone)
"""
import sys, json, time
import numpy as np
from scipy.integrate import solve_ivp
from scipy.spatial import cKDTree
sys.path.insert(0, '/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03/work')
from whlib import find_cycle, cycle_profile, unstable_manifold_orbit, equilibrium_data, kH, field

def riccati_along(mu, k, z0, e0, t_end, n=4001):
    def f(t, s):
        U, p, y, eta, zeta = s
        x = k*np.exp(U)
        return [p, k-mu*p-y, x-y, zeta+(mu-1)*eta-eta*eta, x-(1+eta)*zeta]
    t = np.linspace(0, t_end, n)
    sol = solve_ivp(f, [0, t_end], list(z0)+list(e0), method='DOP853', rtol=1e-11, atol=1e-14, t_eval=t)
    return sol.t, sol.y

def run(mu, rho, thetas=(0.0, 0.7, 1.4, 2.1, 2.8, 3.5, 4.2, 4.9, 5.6), verbose=True):
    k = rho*kH(mu)
    out = dict(mu=mu, rho=rho, k=k)
    ed = equilibrium_data(mu, k)
    cyc = find_cycle(mu, k)
    pr = cycle_profile(mu, k, cyc, n=1501)
    out.update(T=cyc['T'], m_t=cyc['m_t'], lam_c=cyc['lambda_c'], lam_c_riccati=pr['lambda_c_riccati'])
    orbs = [unstable_manifold_orbit(mu, k, r0=1e-6, theta=th, n=6001) for th in thetas]
    # ---- N1: pull-back uniqueness along first orbit, from very different positive data
    o = orbs[0]
    i0 = len(o['t'])//3
    z0 = [o['U'][i0], o['p'][i0], o['y'][i0]]
    tA, sA = riccati_along(mu, k, z0, [1e-3, 1e-3], o['t'][-1]-o['t'][i0], n=len(o['t'])-i0)
    tB, sB = riccati_along(mu, k, z0, [60.0, 300.0], o['t'][-1]-o['t'][i0], n=len(o['t'])-i0)
    ref_eta, ref_zeta = o['eta'][i0:], o['zeta'][i0:]
    dA = np.abs(sA[3]-ref_eta)+np.abs(sA[4]-ref_zeta)
    dB = np.abs(sB[3]-ref_eta)+np.abs(sB[4]-ref_zeta)
    jj = np.searchsorted(tA, 25.0)
    out['N1_initial_gap'] = float(max(dA[0], dB[0]))
    out['N1_gap_after_t25'] = float(max(dA[jj:].max(), dB[jj:].max()))
    # ---- collect attractor sample
    pts = []
    for o in orbs:
        r = np.sqrt(o['U']**2 + o['p']**2 + ((o['y']-k)/k)**2)
        m = r > 2e-3
        pts.append(np.c_[o['U'][m], o['p'][m], o['y'][m], o['eta'][m], o['zeta'][m], o['K'][m], o['W'][m]])
        # N2
        O = o['zeta']*o['p'] - o['eta']*(k-mu*o['p']-o['y']) - (o['x']-o['y'])
        sc = np.abs(o['zeta']*o['p'])+np.abs(o['eta']*(k-mu*o['p']-o['y']))+np.abs(o['x']-o['y'])
        out['N2_max_rel_O'] = max(out.get('N2_max_rel_O', 0.0), float(np.max(np.abs(O[m])/sc[m])))
    P = np.vstack(pts)
    C = np.c_[pr['U'], pr['p'], pr['y'], pr['eta'], pr['zeta'], pr['K'], pr['W']][:-1]
    Ept = np.array([[0.0, 0.0, k, ed['etaE'], ed['zetaE'], 0.0, 0.0]])
    out['n_Wu_points'] = int(len(P)); out['n_cycle_points'] = int(len(C))
    # ---- N6 bounds
    allp = np.vstack([P, C, Ept])
    xm, xM = k*np.exp(allp[:, 0].min()), k*np.exp(allp[:, 0].max())
    eta_p = (mu-1+np.sqrt((mu-1)**2+4*xM))/2; zeta_m = xm/(1+eta_p)
    eta_m = (mu-1+np.sqrt((mu-1)**2+4*zeta_m))/2; zeta_p = xM/(1+eta_m)
    out['N6_bounds_ok'] = bool((allp[:, 3] >= eta_m-1e-9).all() and (allp[:, 3] <= eta_p+1e-9).all()
                               and (allp[:, 4] >= zeta_m-1e-9).all() and (allp[:, 4] <= zeta_p+1e-9).all())
    out['N6_eta_range'] = [float(allp[:, 3].min()), float(allp[:, 3].max()), float(eta_m), float(eta_p)]
    out['N6_eta_minus_max0mu1'] = float(allp[:, 3].min()-max(0, mu-1))
    out['N6_minK_Wu'] = float(P[:, 5].min()); out['N6_minW_Wu'] = float(P[:, 6].min())
    # cone exclusion: no attractor point (other than E) in C+ or C-
    rU, rp, ry = P[:, 0], P[:, 1], P[:, 2]-k
    out['N6_points_in_cones'] = int(np.sum((rU >= 0) & (rp <= 0) & (ry <= 0)) + np.sum((rU <= 0) & (rp >= 0) & (ry >= 0)))
    LamE = ry - ed['zetaE']*rU + ed['etaE']*rp
    out['N6_min_LambdaE_Wu'] = float(LamE.min())
    # ---- N3: all pairs supporting plane, subsample to keep memory small
    rng = np.random.default_rng(0)
    idx = rng.choice(len(P), size=min(2500, len(P)), replace=False)
    S = np.vstack([P[idx], C[::3], Ept])
    q = S[:, :2]; yv = S[:, 2]; g = np.c_[S[:, 4], -S[:, 3]]
    # Delta_ab = g_a.(q_b-q_a) - (y_b-y_a)  (must be <0 for a != b)
    D = (g @ q.T) - np.sum(g*q, axis=1)[:, None] - (yv[None, :]-yv[:, None])
    d2 = ((q[:, None, :]-q[None, :, :])**2).sum(-1) + (yv[None, :]-yv[:, None])**2
    np.fill_diagonal(D, -1.0); np.fill_diagonal(d2, 1.0)
    out['N3_pairs'] = int(D.size-len(S))
    out['N3_max_Delta'] = float(D.max())
    out['N3_num_nonneg'] = int((D >= 0).sum())
    out['N3_max_Delta_over_dist2'] = float((D/d2).max())
    # ---- N4: C^1 matching by local plane fits using points from OTHER orbits only
    tree = cKDTree(P[:, :2])
    sizes = np.cumsum([0]+[len(a) for a in pts])
    lab = np.zeros(len(P), dtype=int)
    for j in range(len(pts)):
        lab[sizes[j]:sizes[j+1]] = j
    test = rng.choice(len(P), size=600, replace=False)
    errs = []; hs = []; resid = []
    for i in test:
        dd, nn = tree.query(P[i, :2], k=60)
        nn = nn[lab[nn] != lab[i]][:14]
        if len(nn) < 8:
            continue
        dq = P[nn, :2]-P[i, :2]
        Amat = np.c_[np.ones(len(nn)), dq]
        coef, *_ = np.linalg.lstsq(Amat, P[nn, 2], rcond=None)
        hloc = np.max(np.linalg.norm(dq, axis=1))
        gi = np.array([P[i, 4], -P[i, 3]])
        errs.append(np.linalg.norm(coef[1:]-gi)/np.linalg.norm(gi)); hs.append(hloc)
        resid.append(abs(coef[0]-P[i, 2])/max(1e-300, hloc))   # value mismatch / h  (C^0 injectivity)
    errs = np.array(errs); hs = np.array(hs); resid = np.array(resid)
    out['N4_n'] = int(len(errs))
    out['N4_median_rel_grad_err'] = float(np.median(errs)); out['N4_max_rel_grad_err'] = float(errs.max())
    out['N4_median_stencil_radius'] = float(np.median(hs))
    out['N4_max_value_mismatch_over_h'] = float(resid.max())
    # first-order Taylor remainder across different orbits: r_ab/|dq| -> 0 and >= 0
    dd, nn = tree.query(P[test, :2], k=40)
    rem = []; mono = []
    for row, i in enumerate(test):
        cand = nn[row][lab[nn[row]] != lab[i]][:5]
        for j in cand:
            dq = P[j, :2]-P[i, :2]; n2 = np.linalg.norm(dq)
            if n2 == 0: continue
            gi = np.array([P[i, 4], -P[i, 3]]); gj = np.array([P[j, 4], -P[j, 3]])
            rem.append((P[j, 2]-P[i, 2]-gi@dq)/n2)
            mono.append(((gj-gi)@dq)/n2**2)
    rem = np.array(rem); mono = np.array(mono)
    out['N4_remainder_over_dq_min'] = float(rem.min()); out['N4_remainder_over_dq_max'] = float(rem.max())
    out['N4_monotone_gradient_min'] = float(mono.min())
    # continuity of (eta,zeta) up to the boundary cycle: nearest W^u point to each cycle point
    dd, nn = tree.query(C[:, :2], k=1)
    out['N4_boundary_max_dist'] = float(dd.max())
    out['N4_boundary_max_grad_jump'] = float(np.max(np.abs(P[nn, 3:5]-C[:, 3:5])))
    out['N4_boundary_max_y_jump'] = float(np.max(np.abs(P[nn, 2]-C[:, 2])))
    # ---- N5: off-attractor orbits; H(q) = max over sampled planes
    Sall = np.vstack([P[::2], C, Ept])
    qS = Sall[:, :2]; yS = Sall[:, 2]; gS = np.c_[Sall[:, 4], -Sall[:, 3]]
    cS = yS - np.sum(gS*qS, axis=1)
    def Hfun(Q):
        return np.max(gS @ Q.T + cS[:, None], axis=0)
    # on-attractor orbit not used in the sample? use held-out half of P
    held = P[1::2][::7]
    gap_on = held[:, 2]-Hfun(held[:, :2])
    out['N5_on_attractor_gap_max_abs'] = float(np.max(np.abs(gap_on)))
    res5 = []
    for (i, dy) in [(len(P)//2, 0.3*k), (len(P)//3, -0.3*k), (2*len(P)//3, 0.05*k)]:
        z0 = P[i, :3].copy(); z0[2] += dy
        tt = np.linspace(0, 12.0, 241)
        sol = solve_ivp(field(mu, k), [0, 12.0], z0, method='DOP853', rtol=1e-11, atol=1e-13, t_eval=tt)
        gap = sol.y[2]-Hfun(sol.y[:2].T)
        res5.append(dict(dy=float(dy), gap0=float(gap[0]), gap_t3=float(gap[60]), gap_t6=float(gap[120]), gap_t12=float(gap[-1]),
                         min_abs_gap_first_2=float(np.min(np.abs(gap[:40])))))
    out['N5_off'] = res5
    # ---- N8: crests of spiral turns monotone, converge to cycle crest
    o = orbs[0]
    pp = o['p']; UU = o['U']
    cr = [UU[i] for i in range(1, len(pp)) if pp[i-1] > 0 >= pp[i] and UU[i] > 0]
    cr = np.array(cr)
    out['N8_turn_crests_monotone'] = bool(np.all(np.diff(cr) > -1e-9))
    out['N8_last_crest_vs_cycle'] = [float(cr[-1]), float(cyc['U0'])]
    if verbose:
        print(json.dumps(out, indent=1))
    return out

if __name__ == '__main__':
    t0 = time.time()
    res = []
    for mu, rho in [(1.0, 1.1), (1.0, 2.0), (0.3, 1.5), (3.0, 1.3)]:
        res.append(run(mu, rho, verbose=False))
        r = res[-1]
        print(f"mu={mu} rho={rho}: N1 gap {r['N1_initial_gap']:.2e}->{r['N1_gap_after_t25']:.2e} | N2 relO {r['N2_max_rel_O']:.1e} | "
              f"N3 maxDelta {r['N3_max_Delta']:.2e} nonneg {r['N3_num_nonneg']}/{r['N3_pairs']} | "
              f"N4 grad err med {r['N4_median_rel_grad_err']:.1e} max {r['N4_max_rel_grad_err']:.1e} (stencil {r['N4_median_stencil_radius']:.1e}) "
              f"rem/dq in [{r['N4_remainder_over_dq_min']:.1e},{r['N4_remainder_over_dq_max']:.1e}] | bdry jump {r['N4_boundary_max_grad_jump']:.1e} at dist {r['N4_boundary_max_dist']:.1e} | "
              f"N5 on {r['N5_on_attractor_gap_max_abs']:.1e} off {[ (round(a['gap0'],3), float('%.1e'%a['gap_t6']), float('%.1e'%a['gap_t12'])) for a in r['N5_off']]} | "
              f"N6 ok {r['N6_bounds_ok']} cones {r['N6_points_in_cones']} minK {r['N6_minK_Wu']:.1e} minW {r['N6_minW_Wu']:.1e} minLamE {r['N6_min_LambdaE_Wu']:.1e} | N8 {r['N8_turn_crests_monotone']} {r['N8_last_crest_vs_cycle']} [{time.time()-t0:.0f}s]", flush=True)
    json.dump(res, open('/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03/work/lanes/planar/job1_checks.json', 'w'), indent=1)
