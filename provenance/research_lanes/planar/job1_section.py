#!/usr/bin/env python3
"""Section test of the attractor-graph theorem (sharp, event-refined).
On the half-line section {p=0, U>0} every attractor orbit must hit ONE convex increasing curve y=h(U,0)
with dy/dU = zeta.  We merge the crest events of several W^u(E) orbits (different spirals), sort by U and test
   (a) y strictly increasing in U (graph property / injectivity),
   (b) zeta_i <= (y_{i+1}-y_i)/(U_{i+1}-U_i) <= zeta_{i+1}   (convexity + grad h = (zeta,-eta), across DIFFERENT spirals),
   (c) same on the section {U=0, p>0} with slope -eta (in p).
Off-attractor orbits started at the same (U,p) with shifted y: their section points are off the curve.
"""
import sys, json
import numpy as np
from scipy.integrate import solve_ivp
sys.path.insert(0, '/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03/work')
from whlib import find_cycle, cycle_profile, equilibrium_data, kH

def ext(mu, k):
    def f(t, s):
        U, p, y, eta, zeta = s
        x = k*np.exp(U)
        return [p, k-mu*p-y, x-y, zeta+(mu-1)*eta-eta*eta, x-(1+eta)*zeta]
    return f

def sections(mu, k, z0, e0, tend):
    def ev1(t, s): return s[1]
    ev1.direction = -1
    def ev2(t, s): return s[0]
    ev2.direction = 1
    sol = solve_ivp(ext(mu, k), [0, tend], list(z0)+list(e0), method='DOP853', rtol=1e-12, atol=1e-15, events=[ev1, ev2])
    return sol.y_events[0], sol.y_events[1]

def run(mu, rho, nth=7):
    k = rho*kH(mu); ed = equilibrium_data(mu, k); cyc = find_cycle(mu, k)
    tend = (np.log(1e6)+10)/ed['lam_u'].real
    S1 = []; S2 = []
    for j, th in enumerate(np.linspace(0, 2*np.pi, nth, endpoint=False)):
        v = np.real(np.exp(1j*th)*ed['vec_u']); z0 = np.array([0, 0, k])+1e-6*v/np.linalg.norm(v)
        a, b = sections(mu, k, z0, [ed['etaE'], ed['zetaE']], tend)
        a = a[(a[:, 0] > 1e-3) & (a[:, 0] < cyc['U0']-1e-5)]; b = b[(b[:, 1] > 1e-3)]
        S1.append(np.c_[a, np.full(len(a), j)]); S2.append(np.c_[b, np.full(len(b), j)])
    S1 = np.vstack(S1); S1 = S1[np.argsort(S1[:, 0])]
    res = dict(mu=mu, rho=rho, n_crest_points=int(len(S1)))
    dU = np.diff(S1[:, 0]); dy = np.diff(S1[:, 2]); sl = dy/dU
    mixed = S1[1:, 5] != S1[:-1, 5]          # neighbours from different spirals
    res['a_neighbour_pairs_from_different_spirals'] = int(mixed.sum())
    res['a_y_increasing'] = bool(np.all(dy > 0))
    lo = (sl-S1[:-1, 4])/S1[:-1, 4]; hi = (S1[1:, 4]-sl)/S1[1:, 4]
    res['b_min_rel_(slope-zeta_i)'] = float(lo.min()); res['b_min_rel_(zeta_i+1-slope)'] = float(hi.min())
    res['b_max_rel_bracket_width'] = float(((S1[1:, 4]-S1[:-1, 4])/S1[:-1, 4]).max())
    # violations expressed in y-units (to be compared with the integration tolerance ~1e-11..1e-10)
    res['b_max_abs_violation_y_units'] = float(max(np.max(S1[:-1, 4]*dU-dy), np.max(dy-S1[1:, 4]*dU), 0.0))
    big = dU > 1e-4
    res['b_pairs_with_dU>1e-4'] = int(big.sum())
    res['b_strict_bracket_on_those'] = bool(np.all(lo[big] > 0) and np.all(hi[big] > 0))
    # section U=0, p>0 : coordinate p, slope dh/dp = -eta
    cyc_p = None
    S2 = np.vstack(S2); S2 = S2[np.argsort(S2[:, 1])]
    pr = cycle_profile(mu, k, cyc, n=4001)
    iu = np.nonzero((pr['U'][:-1] < 0) & (pr['U'][1:] >= 0))[0][0]
    pcyc = pr['p'][iu]
    S2 = S2[S2[:, 1] < pcyc-1e-4]
    dp = np.diff(S2[:, 1]); dy2 = np.diff(S2[:, 2]); sl2 = dy2/dp
    # convex in p with derivative -eta: -eta_i <= slope <= -eta_{i+1}
    lo2 = (sl2+S2[:-1, 3])/S2[:-1, 3]; hi2 = (-S2[1:, 3]-sl2)/S2[1:, 3]
    res['c_n'] = int(len(S2)); res['c_y_decreasing_in_p'] = bool(np.all(dy2 < 0))
    res['c_min_rel_(slope+eta_i)'] = float(lo2.min()); res['c_min_rel_(-eta_i+1-slope)'] = float(hi2.min())
    res['c_max_abs_violation_y_units'] = float(max(np.max(-S2[:-1, 3]*dp-dy2), np.max(dy2+S2[1:, 3]*dp), 0.0))
    big2 = dp > 1e-4
    res['c_pairs_with_dp>1e-4'] = int(big2.sum())
    res['c_strict_bracket_on_those'] = bool(np.all(lo2[big2] > 0) and np.all(hi2[big2] > 0))
    # off-attractor: same starting (U,p) as an attractor crest, y shifted
    off = []
    for dyv in (0.2*k, -0.2*k, 0.02*k):
        i = len(S1)//2
        z0 = [S1[i, 0], 0.0, S1[i, 2]+dyv]
        a, b = sections(mu, k, z0, [S1[i, 3], S1[i, 4]], 4*cyc['T'])
        a = a[1:] if abs(a[0, 0]-z0[0]) < 1e-9 else a
        # distance of off-attractor crest points from the curve (linear interpolation in U), within range
        gaps = []
        for row in a[:4]:
            if S1[0, 0] < row[0] < S1[-1, 0]:
                gaps.append(float(row[2]-np.interp(row[0], S1[:, 0], S1[:, 2])))
            else:
                gaps.append(float('nan'))   # projects outside pi(W^u(E)) on this section
        off.append(dict(dy=float(dyv), crest_gaps=gaps, crest_U=[float(r[0]) for r in a[:4]], cycle_U0=float(cyc['U0'])))
    res['off_attractor'] = off
    return res

if __name__ == '__main__':
    allr = []
    for mu, rho in [(1.0, 1.1), (1.0, 2.0), (0.3, 1.5), (3.0, 1.3), (0.1, 1.3), (1.0, 3.0)]:
        r = run(mu, rho); allr.append(r)
        print(json.dumps({kk: vv for kk, vv in r.items() if kk != 'off_attractor'}))
        for o in r['off_attractor']:
            print('   off', json.dumps(o), flush=True)
    json.dump(allr, open('job1_section.json', 'w'), indent=1)
