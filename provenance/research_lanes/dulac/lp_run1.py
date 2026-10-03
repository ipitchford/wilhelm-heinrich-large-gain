#!/usr/bin/env python3
"""Step 1: per-parameter-set LP feasibility for several structured families, both signs.
LP is solved on a ~8000-row subsample and the optimum is re-checked on the full sample."""
import sys, time
import numpy as np
from lpkit import *

GAUGE = ['logzeta', 'logeta', 'logK', 'logW', 'logL', 'U', 'logy', 'loglamp', 'log1eta']
WRON = ['logVrad', 'logQtan', 'logQ0']
RAT = ['p/lamp', 'H/L', 'K/L', 'W/K', 'K/Vrad', 'W/Qtan', 'W/Vrad', 'Qtan/Vrad', 'L/Vrad', 'eta', 'zeta', 'x']
FAMS = {
    'G  (gauge logs)': (GAUGE, 0),
    'GW (gauge+wronskian logs)': (GAUGE + WRON, 0),
    'GW+poly1': (GAUGE + WRON + ['p', 'H'], 0),
    'GW+poly2': (GAUGE + WRON, 2),
    'GW+rat': (GAUGE + WRON + RAT + ['p', 'H'], 0),
    'GW+rat+poly3': (GAUGE + WRON + RAT, 3),
}

if __name__ == '__main__':
    data = {s: load(*s) for s in SETS}
    for fam, (names, pdg) in FAMS.items():
        print(f'=== family {fam} ===', flush=True)
        for s in SETS:
            u = data[s]; us = subn(u)
            out = []
            for sign in (-1, +1):
                nm = list(names)
                if pdg > 0:
                    f = features(us, polydeg=pdg)
                    nm = list(names) + [n for n in f if ('*' in n or n in ('p', 'H')) and n not in names]
                t, c, msg = solve_lp([us], nm, sign=sign, polydeg=pdg, wkind='one', bound=200.0)
                if t is None:
                    out.append('fail')
                else:
                    R = residual(u, c)
                    tfull = np.min(-sign * R)
                    out.append(f'{t/u["sig"]:+.3f} (full {tfull/u["sig"]:+.3f})')
            print(f'  mu={s[0]:<5} rho={s[1]:<5}  t*/sigma(neg)={out[0]}   t*/sigma(pos)={out[1]}', flush=True)
