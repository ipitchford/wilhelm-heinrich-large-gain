#!/usr/bin/env python3
"""One coefficient vector for several parameter sets (family f, 26 features), negative sign."""
import numpy as np
from lpkit import *
from cplp import build, cplp
from lp_run2 import FAMS
names = FAMS['f: e+ratios']
data = {s: sub(load(*s), 3) for s in SETS}
groups = {'all 10 sets': SETS, 'mu=1 (rho 1.02,1.5,5)': [s for s in SETS if s[0] == 1.0], 'mu=0.05': [s for s in SETS if s[0] == 0.05],
          'near Hopf (rho<=1.05)': [s for s in SETS if s[1] <= 1.05], 'mu=1: rho 1.02,1.5': [(1.0, 1.02), (1.0, 1.5)]}
for gname, sets in groups.items():
    us = [data[s] for s in sets]
    M, b, w = build(us, names)
    # margin in units of each set's sigma
    w = np.concatenate([np.full(len(u['U']), u['sig']) for u in us])
    t, c = cplp(M, b, w, sign=-1, bound=1e4)
    print(f'{gname:<28} shared-coefficient t*/sigma = ' + ('LP fail' if t is None else f'{t:+.3f}'), flush=True)
