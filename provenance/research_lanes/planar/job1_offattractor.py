#!/usr/bin/env python3
"""Off-attractor test (event-refined).  The attractor meets each half-axis section in a C^1 convex curve
(Hermite-interpolated from W^u(E) section points with the exact slopes zeta, -eta).  An orbit started at an
attractor crest with y shifted by dy has equal (U,p) and distinct y at t=0; at its next section crossings the
gap y - h(U,p) is predicted, to first order in dy, by  gap ~ dy * exp(-int_0^t (1+eta) ds) / (correction),
more precisely Delta(t) = n1(t).(z2(t)-z1(t)) obeys Delta' + (1+eta1) Delta = -x1[e^d-1-d] exactly.
We report gap/dy at successive crossings and the linear prediction exp(-int (1+eta))."""
import sys, json
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicHermiteSpline
sys.path.insert(0, '/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03/work')
from whlib import find_cycle, equilibrium_data, kH

def ext(mu, k):
    def f(t, s):
        U, p, y, eta, zeta, I = s
        x = k*np.exp(U)
        return [p, k-mu*p-y, x-y, zeta+(mu-1)*eta-eta*eta, x-(1+eta)*zeta, 1+eta]
    return f
def evp_dn(t, s): return s[1]
evp_dn.direction = -1
def evp_up(t, s): return s[1]
evp_up.direction = 1
def evU_up(t, s): return s[0]
evU_up.direction = 1
def evU_dn(t, s): return s[0]
evU_dn.direction = -1
EVS = [evp_dn, evU_dn, evp_up, evU_up]   # crest, U=0 (p<0), trough, U=0 (p>0)

def run(mu, rho, nth=9):
    k = rho*kH(mu); ed = equilibrium_data(mu, k); cyc = find_cycle(mu, k)
    tend = (np.log(1e6)+10)/ed['lam_u'].real
    S = [[], [], [], []]
    for th in np.linspace(0, 2*np.pi, nth, endpoint=False):
        v = np.real(np.exp(1j*th)*ed['vec_u']); z0 = np.array([0, 0, k])+1e-6*v/np.linalg.norm(v)
        sol = solve_ivp(ext(mu, k), [0, tend], list(z0)+[ed['etaE'], ed['zetaE'], 0.0], method='DOP853', rtol=1e-12, atol=1e-15, events=EVS)
        for j in range(4): S[j].append(sol.y_events[j])
    curves = []
    for j in range(4):
        A = np.vstack(S[j]); c = 0 if j in (0, 2) else 1       # coordinate along the section
        A = A[np.abs(A[:, c]) > 1e-3]; A = A[np.argsort(A[:, c])]
        keep = np.r_[True, np.diff(A[:, c]) > 1e-7]; A = A[keep]
        slope = A[:, 4] if c == 0 else -A[:, 3]
        curves.append((CubicHermiteSpline(A[:, c], A[:, 2], slope), A[0, c], A[-1, c], c))
    res = dict(mu=mu, rho=rho, tests=[])
    A0 = np.vstack(S[0]); A0 = A0[(A0[:, 0] > 0.3*cyc['U0']) & (A0[:, 0] < 0.6*cyc['U0'])]
    base = A0[len(A0)//2]
    for dy in (1e-2*k, -1e-2*k, 0.3*k):
        z0 = [base[0], 0.0, base[2]+dy, base[3], base[4], 0.0]
        sol = solve_ivp(ext(mu, k), [0, 1.2*cyc['T']], z0, method='DOP853', rtol=1e-12, atol=1e-15, events=EVS)
        solb = solve_ivp(ext(mu, k), [0, 1.2*cyc['T']], [base[0], 0.0, base[2], base[3], base[4], 0.0], method='DOP853', rtol=1e-12, atol=1e-15, dense_output=True)
        rows = []
        for j in (1, 2, 3, 0):
            for te, ye in zip(sol.t_events[j], sol.y_events[j]):
                if te < 1e-9: continue
                spl, lo, hi, c = curves[j]
                if not (lo < ye[c] < hi):
                    rows.append((float(te), j, None, None)); continue
                gap = ye[2]-float(spl(ye[c]))
                pred = float(np.exp(-solb.sol(te)[5]))
                rows.append((float(te), j, gap/dy, pred))
        rows.sort()
        res['tests'].append(dict(dy_over_k=dy/k, rows=rows[:5]))
    # control: an attractor orbit (held-out spiral) evaluated against the Hermite curves
    v = np.real(np.exp(1j*0.37)*ed['vec_u']); z0 = np.array([0, 0, k])+1e-6*v/np.linalg.norm(v)
    sol = solve_ivp(ext(mu, k), [0, tend], list(z0)+[ed['etaE'], ed['zetaE'], 0.0], method='DOP853', rtol=1e-12, atol=1e-15, events=EVS)
    mx = 0.0
    for j in range(4):
        spl, lo, hi, c = curves[j]
        for ye in sol.y_events[j]:
            if lo < ye[c] < hi:
                mx = max(mx, abs(ye[2]-float(spl(ye[c])))/k)
    res['control_heldout_attractor_orbit_max_gap_over_k'] = mx
    return res

if __name__ == '__main__':
    out = []
    for mu, rho in [(1.0, 1.1), (1.0, 2.0), (0.3, 1.5), (3.0, 1.3)]:
        r = run(mu, rho); out.append(r)
        print(f"mu={mu} rho={rho}: held-out attractor orbit max |gap|/k = {r['control_heldout_attractor_orbit_max_gap_over_k']:.2e}")
        for t in r['tests']:
            print('   dy/k=%+.2g :' % t['dy_over_k'], ' '.join(('[t=%.2f sec%d gap/dy=%s pred=%s]' % (a, j, ('%.3e' % g) if g is not None else 'outside', ('%.3e' % p) if p is not None else '-')) for a, j, g, p in t['rows']))
    json.dump(out, open('job1_offattractor.json', 'w'), indent=1)
