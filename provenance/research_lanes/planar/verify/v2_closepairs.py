#!/usr/bin/env python3
"""Referee: are the isolated positive values of the supporting-plane functional for NEAR-COINCIDENT attractor
points (lane: +3e-10 at (1,1.1); referee v1: +5e-14 at (1,1.01)) rounding, or a genuine failure of (A3)/(A4)
at small scales near the cycle?

Sharp test on the crest section {p=0, p'<0} of one orbit of W^u(E): consecutive crests a_j=(u_j,0,y_j) with
adjoint ratio zeta_j.  (A3) both ways is     zeta_j <= (y_{j+1}-y_j)/(u_{j+1}-u_j) <= zeta_{j+1},
i.e. c_fwd = -2 Delta(a_j -> a_{j+1})/d^2 > 0 and c_bwd = -2 Delta(a_{j+1} -> a_j)/d^2 > 0, d = u_{j+1}-u_j;
for a C^2 convex h both tend to h_uu(u_j,0).  Section points are hit EXACTLY (Henon's trick: integrate in the
variable p from the last step before the crossing to p=0), so no dense-output interpolation error enters.
Three tolerance levels.  If the theorem holds: the numbers agree where 1e-12/d^2 is small and any sign
anomaly moves/vanishes as the tolerance is tightened.   usage: v2_closepairs.py mu rho
"""
import sys, math, json
import numpy as np
from scipy.integrate import DOP853, solve_ivp
sys.path.insert(0, '/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03/work/lanes/planar/verify')
from vlog import eq_data

mu = float(sys.argv[1]); rho = float(sys.argv[2])
k = rho * mu * (1 + mu)
ed = eq_data(mu, k)


def rhs(t, s):
    U, p, y, eta, zeta = s
    x = k * math.exp(U)
    return np.array([p, k - mu * p - y, x - y, zeta + (mu - 1) * eta - eta * eta, x - (1 + eta) * zeta])


def henon(pvar, s):                    # s = (t, U, y, eta, zeta), independent variable p
    t, U, y, eta, zeta = s
    x = k * math.exp(U)
    pp = k - mu * pvar - y
    return np.array([1.0, pvar, x - y, zeta + (mu - 1) * eta - eta * eta, x - (1 + eta) * zeta]) / pp


def crests(rtol, atol, r0=1e-3, theta=0.0, max_turns=5000, tol=1e-10):
    v = ed['v'] / np.linalg.norm(ed['v'])
    d = np.real(np.exp(1j * theta) * v)
    s0 = np.array([r0 * d[0], r0 * d[1], k + r0 * d[2], ed['eta'], ed['zeta']])
    sol = DOP853(rhs, 0.0, s0, np.inf, rtol=rtol, atol=atol)
    out = []
    while len(out) < max_turns:
        t_old, y_old = sol.t, sol.y.copy()
        sol.step()
        if y_old[1] > 0 and sol.y[1] <= 0:
            h0 = np.array([t_old, y_old[0], y_old[2], y_old[3], y_old[4]])
            r = solve_ivp(henon, [y_old[1], 0.0], h0, method='DOP853', rtol=rtol, atol=atol)
            te, U, y, eta, zeta = r.y[:, -1]
            out.append([te, U, y, eta, zeta])
            if len(out) > 3 and abs(out[-1][1] - out[-2][1]) < tol * out[-1][1]:
                break
    return np.array(out)

levels = [(1e-10, 1e-12), (1e-12, 1e-14), (3e-14, 1e-16)]
tabs = [crests(rt, at) for rt, at in levels]
n = min(len(t) for t in tabs)
skip = int(np.searchsorted(tabs[-1][:, 0], 35.0 / ed['gam'])) + 1
res = []
print(f"mu={mu} rho={rho} k={k:.5g}: {n} crests; columns per tolerance level (rtol): c_fwd, c_bwd")
print("   j        u_j          d=u_{j+1}-u_j " + ''.join(f"|   rtol={rt:.0e}: c_fwd      c_bwd    " for rt, _ in levels))
# pick crests so that d spans the decades on the way out and on the approach to the cycle
dref = np.diff(tabs[-1][:n, 1])
jmax = int(np.argmax(dref))
targets = [1e-4, 1e-5, 1e-6, 3e-7, 1e-7, 3e-8, 1e-8]
sel = {skip, (skip + jmax) // 2, jmax}
for tg in targets:
    cand = np.where(dref[jmax:] < tg)[0]
    if len(cand):
        sel.add(int(jmax + cand[0]))
anom = [0, 0, 0]
allc = []
for li, tb in enumerate(tabs):
    d = np.diff(tb[:n, 1]); dy = np.diff(tb[:n, 2])
    cf = -2 * (tb[:n - 1, 4] * d - dy) / d ** 2
    cb = -2 * (-tb[1:n, 4] * d + dy) / d ** 2
    allc.append((d, cf, cb))
    anom[li] = int(((cf[skip:] <= 0) | (cb[skip:] <= 0)).sum())
for j in sorted(sel):
    if j >= n - 1:
        continue
    line = f"{j:5d} {tabs[-1][j,1]:14.10f} {dref[j]:12.3e}   "
    row = dict(j=int(j), u=float(tabs[-1][j, 1]), d=float(dref[j]))
    for li in range(3):
        d, cf, cb = allc[li]
        line += f"|        {cf[j]:+11.4e} {cb[j]:+11.4e}   "
        row[f"cf{li}"] = float(cf[j]); row[f"cb{li}"] = float(cb[j])
    res.append(row)
    print(line)
print("number of crest pairs (after the start-up transient) with c_fwd<=0 or c_bwd<=0, per tolerance level:", anom)
# smallest separation at which a sign anomaly occurs, per level
for li in range(3):
    d, cf, cb = allc[li]
    bad = np.where((cf[skip:] <= 0) | (cb[skip:] <= 0))[0] + skip
    if len(bad):
        print(f"   rtol={levels[li][0]:.0e}: anomalies only at separations d <= {np.abs(d[bad]).max():.2e}"
              f"  (largest |Delta| among them {np.max(np.abs(np.minimum(cf[bad], cb[bad]) * d[bad]**2 / 2)):.2e})")
    else:
        print(f"   rtol={levels[li][0]:.0e}: no anomaly")
json.dump(dict(mu=mu, rho=rho, n=n, anomalies=anom, rows=res),
          open(f"/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03/work/lanes/planar/verify/out_v2_{mu:g}_{rho:g}.json", 'w'), indent=1)
