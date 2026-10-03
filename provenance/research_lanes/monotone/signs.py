#!/usr/bin/env python3
"""Item 1: sign patterns / zero phases around the cycle for the Floquet objects."""
import numpy as np, warnings
from mono_lib import *
warnings.filterwarnings('ignore')

cases = [(1.0, 1.05), (1.0, 1.3), (1.0, 2.0), (1.0, 4.0),
         (0.5, 1.1), (0.5, 2.0), (0.5, 3.5),
         (0.2, 1.05), (0.2, 1.5), (0.2, 3.0),
         (0.05, 1.2), (0.05, 2.0), (0.05, 4.0),
         (3.0, 1.5), (5.0, 2.5)]


def fz(f, t):
    return ' '.join(f"{ph:.4f}{'^' if d > 0 else 'v'}" for ph, d in zeros_of(f, t))


for mu, rho in cases:
    k = rho * kH(mu)
    c = find_cycle(mu, k)
    d = floquet_data(mu, k, c, n=8001)
    t = d['t']
    P = d['Psi']
    print(f"=== mu={mu} rho={rho} T={d['T']:.4f} lam_c={d['lam']:+.6f} m_t={c['m_t']:.6f} Umax={d['U'].max():.3f} Umin={d['U'].min():.3f}")
    print("  zeros p   :", fz(d['p'], t), " (crest at 0)")
    print("  zeros p'  :", fz(d['pp'], t))
    print("  zeros y'  :", fz(d['yp'], t))
    print("  zeros U   :", fz(d['U'], t))
    print("  zeros y-k :", fz(d['y'] - k, t))
    print("  zeros Psi1:", fz(P[0], t))
    print("  zeros Psi2:", fz(P[1], t))
    print("  zeros Psi3:", fz(P[2], t))
    print("  zeros sig :", fz(d['sig'], t), "  min/max sig", d['sig'].min(), d['sig'].max())
    print("  zeros nu  :", fz(d['nu'], t))
    print("  zeros nu' :", fz(d['nup'], t))
    print("  zeros dy_t:", fz(d['dy'], t))
    kap = P[0] * (0 - d['U']) + P[1] * (0 - d['p']) + P[2] * (k - d['y'])   # Psi.(E-z)
    print("  zeros Psi.(E-z):", fz(kap, t), " min/max", kap.min(), kap.max())
    g = 1 - np.exp(d['U']) + d['U'] * np.exp(d['U'])
    print(f"  <sig>={avg(d['sig'],t):+.5e}  <Psi3 x p^2>={avg(P[2]*d['x']*d['p']**2,t):+.5e}  "
          f"<Psi.(E-z)>={avg(kap,t):+.5e}  k<Psi3 g>={k*avg(P[2]*g,t):+.5e}  check lam<kap>+k<Psi3 g>={d['lam']*avg(kap,t)+k*avg(P[2]*g,t):+.1e}")
    print(f"  <Psi1>={avg(P[0],t):+.5e} <Psi2>={avg(P[1],t):+.5e} <Psi3>={avg(P[2],t):+.5e}  <Psi3 x>={avg(P[2]*d['x'],t):+.5e}")
