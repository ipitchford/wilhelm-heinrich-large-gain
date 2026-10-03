#!/usr/bin/env python3
"""Twelve-event order: zeros of the weak adjoint Psi=(Psi1,Psi2,Psi3) interlace the six extrema.
Claimed cyclic order (phases, crest X_max at 0):
  X_max < [Psi2 down] < Z_max < [Psi1 down] < Y_max < [Psi3 down = t1] < X_min
        < [Psi2 up]   < Z_min < [Psi1 up]   < Y_min < [Psi3 up = t2]   < X_max + T
Also: r > max(1,mu); zero pattern of e_t = (nu, nu', dy); convexity diagnostics
  q1 = p'^2 - p p''  (hodograph rotation / convexity of the (U,p) projection)
  q2 = 4 zeta - (eta-mu)^2 (frozen leading 2x2 block has complex eigenvalues)
  q3 = p^2 - U p'   (star-shapedness of the (U,p) projection)
"""
import numpy as np, warnings
from mono_lib import *
warnings.filterwarnings('ignore')

cases = [(1.0, 1.02), (1.0, 1.3), (1.0, 2.0), (1.0, 3.0), (0.5, 1.1), (0.5, 2.0), (0.5, 3.5),
         (0.2, 1.05), (0.2, 1.5), (0.2, 3.0), (0.1, 1.5), (0.3, 1.01), (0.7, 5.0), (2.0, 1.2), (4.0, 2.0)]
allok = True
for mu, rho in cases:
    k = rho * kH(mu)
    c = find_cycle(mu, k)
    d = floquet_data(mu, k, c, n=6001)
    t = d['t']
    P = d['Psi']

    def z(f, direction):
        zz = [a for a, b in zeros_of(f, t) if b == direction]
        return zz

    try:
        Xmax = 0.0
        Zmax = z(d['yp'], -1)[0]; Zmin = z(d['yp'], +1)[0]
        Ymax = z(d['pp'], +1)[0]; Ymin = z(d['pp'], -1)[0]     # Y_max = p-min
        Xmin = z(d['p'], +1)[0]
        P2d, P2u = z(P[1], -1), z(P[1], +1)
        P1d, P1u = z(P[0], -1), z(P[0], +1)
        P3d, P3u = z(P[2], -1), z(P[2], +1)
        counts = [len(v) for v in (P1d, P1u, P2d, P2u, P3d, P3u)]
        seq = [Xmax, P2d[0], Zmax, P1d[0], Ymax, P3d[0], Xmin, P2u[0], Zmin, P1u[0], Ymin, P3u[0], 1.0]
        ordered = all(seq[i] < seq[i + 1] for i in range(len(seq) - 1)) and counts == [1] * 6
    except Exception as e:
        ordered = False; seq = []; counts = str(e)
    # e_t zeros
    nd, nu_ = z(d['nu'], -1), z(d['nu'], +1)
    npd, npu = z(d['nup'], -1), z(d['nup'], +1)
    dyd, dyu = z(d['dy'], -1), z(d['dy'], +1)
    try:
        seq_t = [Xmax, nd[0], Zmax, dyd[0], Ymax, npu[0], Xmin, nu_[0], Zmin, dyu[0], Ymin, npd[0], 1.0]
        ord_t = all(seq_t[i] < seq_t[i + 1] for i in range(len(seq_t) - 1)) and [len(v) for v in (nd, nu_, npd, npu, dyd, dyu)] == [1] * 6
    except Exception:
        ord_t = False
    q1 = d['pp'] ** 2 - d['p'] * d['ppp']
    q2 = 4 * d['zeta'] - (d['eta'] - mu) ** 2
    q3 = d['p'] ** 2 - d['U'] * d['pp']
    rmin = d['r'].min()
    allok &= ordered
    print(f"mu={mu} rho={rho}: Psi twelve-event order={ordered} counts={counts} | e_t analogue order={ord_t} | "
          f"min r - max(1,mu)={rmin - max(1, mu):+.3g} | min q1={q1.min():.3g} min q2={q2.min():.3g} min q3={q3.min():.3g}")
    if not ordered:
        print("   seq:", np.round(seq, 4))
print("ALL Psi orders OK" if allok else "SOME Psi order FAILED")
