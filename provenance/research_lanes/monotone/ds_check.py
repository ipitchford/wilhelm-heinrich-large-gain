#!/usr/bin/env python3
"""Numerical verification of every step of the dual-supporting-plane (DS) proof, on many cycles.

 S1  r > mu, s > 0 (strong-stable ratios)                              [min(r)-mu, min(s)]
 S2  psi_p = p' + r p has exactly two zeros: t1 in F=(Ymax,Xmin), t2 in C=(Ymin,Xmax)
 S3  on F and C: y'/p > 0                                               [six-extrema order]
 S4  theta=-p'/p crosses the level theta1=r(t1) exactly once in F and once in C
 S5  L(t) = p'(t1)(U-U(t1)) - p(t1)(p-p(t1)) > 0 for t != t1 ; L''(t1) = p^2 r (r-mu) + p y' > 0
 S6  U(t1) < 0
 S7  B_s(t) = Psi(t).(z(t+s)-z(t)) < 0 (all s,t), sigma<0, Psi.(z-E)>0
 Also: J0 = <a^ psi_p> (>0?),  I2-analogue <Psi^3 x nu^2>, and sign formulas with weights.
"""
import numpy as np, warnings, sys
from mono_lib import *
warnings.filterwarnings('ignore')

cases = [(1.0, 1.02), (1.0, 1.05), (1.0, 1.3), (1.0, 2.0), (1.0, 3.0), (0.5, 1.1), (0.5, 2.0), (0.5, 3.5),
         (0.2, 1.05), (0.2, 1.5), (0.2, 3.0), (0.05, 1.2), (0.1, 1.5), (0.3, 1.01), (0.7, 5.0)]
if len(sys.argv) > 1:
    cases = [(float(sys.argv[1]), float(sys.argv[2]))]

allok = True
for mu, rho in cases:
    k = rho * kH(mu)
    c = find_cycle(mu, k)
    n = 4001
    d = floquet_data(mu, k, c, n=n)
    t = d['t']; T = d['T']
    p, pp, yp, U, y, x = d['p'], d['pp'], d['yp'], d['U'], d['y'], d['x']
    r, s, ah = d['r'], d['s'], d['ah']
    psi = pp + r * p
    zpsi = zeros_of(psi, t)
    zp = zeros_of(p, t); zpp = zeros_of(pp, t); zyp = zeros_of(yp, t)
    # phases: crest at 0 (=1), trough = up-zero of p ; p-min = up-zero of p', p-max = down-zero of p'
    ph_tr = [a for a, b in zp if b > 0][0]
    ph_pmin = [a for a, b in zpp if b > 0][0]
    ph_pmax = [a for a, b in zpp if b < 0][0]
    S1 = (r.min() - mu, s.min())
    S2 = (len(zpsi) == 2 and ph_pmin < zpsi[0][0] < ph_tr and ph_pmax < zpsi[1][0] < 1.0)
    ph = t / T
    F = (ph > ph_pmin) & (ph < ph_tr); C = (ph > ph_pmax) & (ph < 1.0)
    S3 = bool(np.all((yp / p)[F][2:-2] > 0) and np.all((yp / p)[C][2:-2] > 0))
    i1 = int(np.argmin(np.abs(ph - zpsi[0][0])))
    th = -pp / p
    th1 = r[i1]
    nF = len(zeros_of((th - th1)[F], t[F])); nC = len(zeros_of((th - th1)[C], t[C]))
    S4 = (nF == 1 and nC == 1)
    L = pp[i1] * (U - U[i1]) - p[i1] * (p - p[i1])
    mask = np.abs(np.arange(n) - i1) > 3
    S5 = (L[mask][:-1].min() > 0, p[i1] ** 2 * r[i1] * (r[i1] - mu) + p[i1] * yp[i1])
    S6 = U[i1]
    # S7
    P = d['Psi']
    Z = np.array([U, p, y])[:, :-1]; Pq = P[:, :-1]; m = n - 1
    bmax = -np.inf
    for j in range(1, 100):
        sh = j * m // 100
        bmax = max(bmax, np.sum(Pq * (np.roll(Z, -sh, axis=1) - Z), axis=0).max())
    kapt = P[0] * U + P[1] * p + P[2] * (y - k)
    S7 = (bmax, d['sig'].max(), kapt.min())
    J0 = avg(ah * psi, t)
    I2 = avg(P[2] * x * d['nu'] ** 2, t)
    g = 1 - np.exp(U) + U * np.exp(U)
    sg = [np.sign(avg(P[2] * w, t)) for w in (-x, x * p ** 2, k * g)]
    ok = (S1[0] > 0 and S1[1] > 0 and S2 and S3 and S4 and S5[0] and S5[1] > 0 and S6 < 0 and bmax < 0 and S7[1] < 0 and S7[2] > 0)
    allok &= ok
    print(f"mu={mu} rho={rho} lam={d['lam']:+.5f} | S1 r-mu,s min=({S1[0]:.3g},{S1[1]:.3g}) S2={S2} t1={zpsi[0][0]:.4f} t2={zpsi[1][0]:.4f} "
          f"S3={S3} S4={S4} S5 minL={S5[0]},L''={S5[1]:.3g} S6 U(t1)={S6:.3g} S7 maxB={S7[0]:.3g} max sig={S7[1]:.3g} min kap~={S7[2]:.3g} "
          f"| J0={J0:+.4g} I2~={I2:+.4g} signs<Psi3 W> (W=-x,xp^2,kg)={sg} OK={ok}")
print("ALL OK" if allok else "SOME FAILED")
