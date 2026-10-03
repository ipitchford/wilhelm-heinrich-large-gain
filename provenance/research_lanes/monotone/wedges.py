#!/usr/bin/env python3
"""Componentwise signs of compound (wedge) vectors along cycles.
w = v ^ z' := (v2 p - v1 p', v1 y' - v3 p, v2 y' - v3 p')   (compound coords, w' = C w for solutions v)
  v = z - E      -> wE   (kappa~ = l.w,  l=(c,b,a))
  v = -z''       -> wC   (z' ^ z'')
  v = e_t        -> wT   (right Perron vector of compound; must be >0)
Also checks l.w identities:  l^.wE  = Psi^.(z-E),  l^.(z'' ^ z') = sigma^.
"""
import numpy as np, warnings
from mono_lib import *
warnings.filterwarnings('ignore')

cases = [(1.0, 1.05), (1.0, 1.3), (1.0, 2.0), (1.0, 4.0), (0.5, 1.1), (0.5, 2.0), (0.5, 3.5),
         (0.2, 1.05), (0.2, 1.5), (0.2, 3.0), (0.05, 1.2), (0.05, 2.0), (3.0, 1.5), (5.0, 2.5)]


def wedge(v, zp):
    v1, v2, v3 = v
    p, pp, yp = zp
    return np.array([v2 * p - v1 * pp, v1 * yp - v3 * p, v2 * yp - v3 * pp])


def frac_neg(f):
    return float(np.mean(f < 0))


for mu, rho in cases:
    k = rho * kH(mu)
    c = find_cycle(mu, k)
    d = floquet_data(mu, k, c, n=8001)
    t = d['t']
    zp = (d['p'], d['pp'], d['yp'])
    wE = wedge((d['U'], d['p'], d['y'] - k), zp)
    wC = wedge((-d['pp'], -d['ppp'], -d['ypp']), zp)
    wT = wedge(tuple(d['et']), zp)
    print(f"mu={mu} rho={rho}: frac(neg) wE={[round(frac_neg(w),3) for w in wE]}  "
          f"z'^z''={[round(frac_neg(w),3) for w in wC]}  e_t^z'={[round(frac_neg(w),3) for w in wT]}")
