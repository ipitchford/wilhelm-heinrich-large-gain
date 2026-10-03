#!/usr/bin/env python3
"""Referee: follow-up to v6 (iv).  A few spiral points were reported outside the 4001-vertex polygon of the projected
cycle.  Chord artefact or genuine?  Re-test with (a) the same polygon, recording at which spiral turn the flagged
points sit and how far outside (relative), (b) a polygon with 40x more vertices, (c) a margin: spiral restricted to
turns whose crest is at least 1e-3 (relative) below the cycle crest.   usage: v6b_polygon.py mu rho"""
import sys, math, json
import numpy as np
from matplotlib.path import Path
sys.path.insert(0, '/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03/work/lanes/planar/verify')
from vlog import *
mu = float(sys.argv[1]); rho = float(sys.argv[2])
k = rho * mu * (1 + mu); lam = k / mu
ed = eq_data(mu, k)
o = json.load(open(f"/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03/work/lanes/planar/verify/out_v1b_{mu:g}_{rho:g}.json"))
cyc = newton_cycle(mu, k, o['U_crest'], o['y_crest'])
t_skip = 35.0 / ed['gam']
r0 = min(1e-3, 1e-3 * math.exp(-ed['lam'].real * t_skip))
tconv, cr = spiral_until_converged(mu, k, r0, 0.3)
crU = cr[:, 1]
for label, relgap, nvert in (("as v6 (gap 1e-6, 4001 vertices)", 1e-6, 4001), ("gap 1e-6, 160001 vertices", 1e-6, 160001),
                             ("gap 1e-3, 160001 vertices", 1e-3, 160001)):
    prof = cycle_riccati(mu, k, cyc, n=nvert)
    S = prof['S']
    poly = Path(np.column_stack((S[0], lam - np.exp(S[1]))))
    jstop = int(np.argmax(crU > crU[-1] * (1 - relgap)))
    tstop = cr[max(jstop - 1, 0), 0]
    te = np.linspace(t_skip, tstop, 6000)
    sol = run_log(mu, k, spiral_start(mu, k, r0, 0.3), 0.0, tstop, t_eval=te)
    Pts = np.column_stack((sol.y[0], lam - np.exp(sol.y[1])))
    inside = poly.contains_points(Pts)
    bad = np.where(~inside)[0]
    msg = ''
    if len(bad):
        turn = np.searchsorted(cr[:, 0], te[bad])
        msg = f"; flagged points sit on spiral turns {turn.min()}..{turn.max()} of {jstop} (crest gaps to cycle: {crU[-1]-crU[min(turn.min(), len(crU)-1)-1]:.2e})"
    print(f"mu={mu} rho={rho}: {label}: {inside.sum()}/{len(Pts)} inside{msg}")
