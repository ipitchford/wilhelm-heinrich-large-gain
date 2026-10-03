#!/usr/bin/env python3
"""Referee: the single spiral point at (mu,rho)=(0.02,12) still flagged outside the 160001-vertex cycle polygon.
Locate it and compare with the cycle AT THE SAME U (root-find on the cycle), in the cancellation-free quantity
q = lam - p (log scale)."""
import sys, math, json
import numpy as np
from matplotlib.path import Path
from scipy.interpolate import CubicSpline
sys.path.insert(0, '/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03/work/lanes/planar/verify')
from vlog import *
import warnings; warnings.filterwarnings('ignore')
mu, rho = 0.02, 12.0
k = rho * mu * (1 + mu); lam = k / mu
ed = eq_data(mu, k)
o = json.load(open(f"out_v1b_{mu:g}_{rho:g}.json"))
cyc = newton_cycle(mu, k, o['U_crest'], o['y_crest'])
t_skip = 35.0 / ed['gam']
r0 = min(1e-3, 1e-3 * math.exp(-ed['lam'].real * t_skip))
tconv, cr = spiral_until_converged(mu, k, r0, 0.3)
crU = cr[:, 1]
prof = cycle_riccati(mu, k, cyc, n=160001)
S = prof['S']
poly = Path(np.column_stack((S[0], lam - np.exp(S[1]))))
jstop = int(np.argmax(crU > crU[-1] * (1 - 1e-6)))
tstop = cr[max(jstop - 1, 0), 0]
te = np.linspace(t_skip, tstop, 6000)
sol = run_log(mu, k, spiral_start(mu, k, r0, 0.3), 0.0, tstop, t_eval=te)
Pts = np.column_stack((sol.y[0], lam - np.exp(sol.y[1])))
inside = poly.contains_points(Pts)
bad = np.where(~inside)[0]
for i in bad:
    Ub, wb = sol.y[0, i], sol.y[1, i]
    print(f"flagged spiral point: t={te[i]:.3f}  U={Ub:.6f}  p={Pts[i,1]:.17g}  log q={wb:.6f}  (lam={lam})")
    # cycle branches at the same U: indices where S[0] crosses Ub
    sgn = np.sign(S[0] - Ub)
    cross = np.where(sgn[:-1] * sgn[1:] < 0)[0]
    for c in cross:
        f = (Ub - S[0, c]) / (S[0, c + 1] - S[0, c])
        wc = S[1, c] + f * (S[1, c + 1] - S[1, c])
        pc = lam - math.exp(wc)
        print(f"   cycle at the same U: branch {'rising (p>0)' if pc > 0 else 'falling (p<0)'}: p={pc:.17g}  log q={wc:.6f};"
              f"  spiral minus cycle: dp={Pts[i,1]-pc:.3e},  d log q={wb-wc:.3e}")
    # neighbours of the flagged point on the spiral
    for j in (i - 1, i + 1):
        if 0 <= j < len(te):
            print(f"   neighbour sample t={te[j]:.3f}: U={sol.y[0,j]:.4f} p={Pts[j,1]:.17g} inside={inside[j]}")

# exact comparison: integrate the cycle BACKWARD from its crest to the same U (event-refined), original variables
from scipy.integrate import solve_ivp
def f3(t, z):
    U, p, y = z
    return [p, k - mu * p - y, k * ex(U) - y]
for i in bad:
    Ub = sol.y[0, i]
    ev = lambda t, z: z[0] - Ub
    ev.terminal = True
    r = solve_ivp(f3, [0, -5.0], [cyc['U0'], 0.0, cyc['y0']], method='DOP853', rtol=1e-12, atol=1e-14, events=ev)
    pc = r.y_events[0][0][1]
    print(f"   EXACT cycle point at U={Ub:.6f} on the rising branch (backward from crest, t={r.t_events[0][0]:.5f}):"
          f" p_cycle={pc:.12f}; spiral p={Pts[i,1]:.12f}; spiral - cycle = {Pts[i,1]-pc:.3e}"
          f"  -> spiral point is {'INSIDE (below the upper arc)' if Pts[i,1] < pc else 'OUTSIDE'}")
