#!/usr/bin/env python3
"""Greedy forward selection: smallest feature family (per-set coefficients) that is Dulac-feasible (negative sign)
on the full disc sample of all 10 parameter sets.  Score of a family = min over sets of t*/sigma."""
import time, sys
import numpy as np
from lpkit import *
from cplp import build, cplp

POOL = ['logVrad', 'logQtan', 'logQ0', 'logzeta', 'logeta', 'logy', 'loglamp', 'log1eta', 'U', 'p', 'H', 'logK', 'logW', 'logL',
        'p/lamp', 'H/L', 'K/L', 'W/K', 'K/Vrad', 'W/Qtan', 'W/Vrad', 'Qtan/Vrad', 'L/Vrad', 'eta', 'zeta', 'x']
BASE = ['logVrad', 'U', 'p', 'H']
data = {s: sub(load(*s), 2) for s in SETS}
feat = {s: features(data[s]) for s in SETS}

def score(names):
    worst = np.inf; per = []
    for s in SETS:
        u = data[s]
        M = np.column_stack([feat[s][n] for n in names]); b = u['eta'] - u['mu']; w = np.ones(len(b))
        t, c = cplp(M, b, w, sign=-1, bound=1e4)
        v = -np.inf if t is None else t / u['sig']
        per.append(v); worst = min(worst, v)
        if worst < -5:
            break
    return worst, per

cur = list(BASE)
w0, per = score(cur)
print('base', cur, f'score={w0:+.3f}', flush=True)
for rnd in range(9):
    best = None
    for cand in POOL:
        if cand in cur:
            continue
        sc, per = score(cur + [cand])
        if best is None or sc > best[0]:
            best = (sc, cand, per)
    cur.append(best[1])
    print(f'round {rnd+1}: add {best[1]:<10} score={best[0]:+.3f}  per-set: ' + ' '.join(f'{v:+.2f}' for v in best[2]), flush=True)
    if best[0] > 0.05:
        break
print('final family:', cur)
