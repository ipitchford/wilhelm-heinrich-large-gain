#!/usr/bin/env python3
"""Small structured families, per-set cutting-plane LP on full samples (negative sign), margin relative to sigma."""
import time, sys
import numpy as np
from lpkit import *
from cplp import build, cplp

FAMS = {
    'a: Vrad,U,p,H': ['logVrad', 'U', 'p', 'H'],
    'b: Vrad,Qtan,U,p,H': ['logVrad', 'logQtan', 'U', 'p', 'H'],
    'c: b+zeta,eta': ['logVrad', 'logQtan', 'logzeta', 'logeta', 'U', 'p', 'H'],
    'd: c+y,lamp,1eta,Q0': ['logVrad', 'logQtan', 'logQ0', 'logzeta', 'logeta', 'logy', 'loglamp', 'log1eta', 'U', 'p', 'H'],
    'e: d+K,W,L': ['logVrad', 'logQtan', 'logQ0', 'logzeta', 'logeta', 'logy', 'loglamp', 'log1eta', 'U', 'p', 'H', 'logK', 'logW', 'logL'],
    'f: e+ratios': ['logVrad', 'logQtan', 'logQ0', 'logzeta', 'logeta', 'logy', 'loglamp', 'log1eta', 'U', 'p', 'H', 'logK', 'logW', 'logL',
                    'p/lamp', 'H/L', 'K/L', 'W/K', 'K/Vrad', 'W/Qtan', 'W/Vrad', 'Qtan/Vrad', 'L/Vrad', 'eta', 'zeta', 'x'],
}
data = {s: load(*s) for s in SETS}
for fam, names in FAMS.items():
    print(f'=== {fam} ===', flush=True)
    for s in SETS:
        u = data[s]
        M, b, w = build([u], names)
        t0 = time.time()
        t, c = cplp(M, b, w, sign=-1, bound=1e4)
        if t is None:
            print(f'  {s}: LP failed'); continue
        cs = ' '.join(f'{n}={v:+.3f}' for n, v in zip(names, c) if abs(v) > 1e-3)
        print(f'  mu={s[0]:<5} rho={s[1]:<5} t*/sigma={t/u["sig"]:+.3f} [{time.time()-t0:.1f}s]  {cs}', flush=True)
