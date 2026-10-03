#!/usr/bin/env python3
"""Explicit pointwise bounds on eta from the constraint line  zeta p - eta p' = y'  plus positivity and K>0.
  p>0:            eta < S2/S1               (from K>0, S1=p^2-Up'>0),  S2 = y'U - p(y-k)
  p p' < 0:       eta < |y'/p'|             (from zeta>0)
  p<0:            eta > S2/S1 ;  p'<0 => zeta > S3/S1 ; p'>0 => zeta < S3/S1
Report: violations (should be none), and the phase-averages  int eta dt  vs  int bound dt  on each phase."""
import numpy as np
from base import load
D = load()
print(f"{'mu':>5} {'rho':>6} {'lam_c':>10} | {'frac p>0':>8} {'<eta>_p>0':>10} {'<B1>_p>0':>10} | {'frac Q3':>7} {'<eta>_Q3':>9} {'<B2>_Q3':>9} | {'frac Q2':>7} {'<eta>_Q2':>9} | {'need':>9} {'viol':>5}")
for key in sorted(D):
    P = D[key]; mu = P['mu']; k = P['k']
    U, p, y, x, eta, zeta = P['U'], P['p'], P['y'], P['x'], P['eta'], P['zeta']
    pd = k - mu*p - y; yd = x - y; H = y - k
    S1 = p*p - U*pd; S2 = yd*U - p*H; S3 = yd*p - pd*H
    Kk = zeta*U - eta*p - H
    B1 = S2/S1
    with np.errstate(divide='ignore', invalid='ignore'):
        B2 = np.abs(yd/pd)
    m1 = p > 0
    mQ3 = (p < 0) & (pd > 0)       # after min of p, before trough
    mQ2 = (p < 0) & (pd < 0)       # after crest, before min of p
    mQ1 = (p > 0) & (pd < 0)       # after max p, before crest
    viol = 0
    viol += int(np.sum(eta[m1] > B1[m1]*(1+1e-7) + 1e-9))
    viol += int(np.sum(eta[mQ3] > B2[mQ3]*(1+1e-7) + 1e-9))
    viol += int(np.sum(eta[mQ1] > B2[mQ1]*(1+1e-7) + 1e-9))
    viol += int(np.sum(eta[~m1] < B1[~m1]*(1-1e-7) - 1e-9))
    Bp = np.where(mQ1, np.minimum(B1, B2), B1)
    f1 = m1.mean(); f3 = mQ3.mean(); f2 = mQ2.mean()
    # contributions to <eta> (already weighted by fraction)
    c1 = np.sum(eta[m1])/len(eta); c1b = np.sum(Bp[m1])/len(eta)
    c3 = np.sum(eta[mQ3])/len(eta); c3b = np.sum(B2[mQ3])/len(eta)
    c2 = np.sum(eta[mQ2])/len(eta)
    print(f"{key[0]:5.2f} {key[1]:6.3f} {P['lam_M']:+10.3e} | {f1:8.3f} {c1:10.4f} {c1b:10.4f} | {f3:7.3f} {c3:9.4f} {c3b:9.4f} | {f2:7.3f} {c2:9.4f} | mu={mu:5.2f} tot_bound(excl Q2 slack)={c1b+c3b+c2:8.4f} {viol:5d}")
