#!/usr/bin/env python3
"""Pure gauge family G (brief 3(d) ansatz without Wronskians) on the full disc samples, both signs; plus the two
one-line a-priori bounds  zeta < y  and (mu <= 1/2)  eta < lambda - p."""
import numpy as np
from lpkit import *
from cplp import build, cplp
G = ['logzeta', 'logeta', 'logK', 'logW', 'logL', 'U', 'logy', 'loglamp', 'log1eta']
for s in SETS:
    u = load(*s)
    M, b, w = build([u], G)
    tn, cn = cplp(M, b, w, sign=-1, bound=1e4)
    tp, cp = cplp(M, b, w, sign=+1, bound=1e4)
    f = lambda t: 'fail' if t is None else f'{t/u["sig"]:+.4f}'
    extra = f'max(zeta-y)/k={np.max(u["zeta"]-u["y"])/u["k"]:+.2e}'
    if s[0] <= 0.5:
        extra += f'  max(eta-(lam-p))={np.max(u["eta"]-u["lamp"]):+.2e}'
    print(f'mu={s[0]:<5} rho={s[1]:<5} family G: t*/sigma  (R<=0): {f(tn)}   (R>=0): {f(tp)}   | {extra}', flush=True)
