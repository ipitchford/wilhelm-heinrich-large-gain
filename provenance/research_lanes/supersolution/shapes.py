#!/usr/bin/env python3
"""Shape diagnostics for Targets B and C on cached cycles."""
import numpy as np
from base import load, avg, ddt
D = load()
def phase(t, T): return t / T
def nturn(f):
    """number of sign changes of f' around the period (periodic)"""
    d = np.diff(np.append(f, f[0])); s = np.sign(d); s = s[s != 0]
    return int(np.sum(s != np.roll(s, 1)))
def nzero(f):
    s = np.sign(f); s = s[s != 0]
    return int(np.sum(s != np.roll(s, 1)))
print("Target B shapes: phases (t/T, crest at 0) of extrema; trough phase = argmin U")
print(f"{'mu':>5} {'rho':>5} {'ph_trough':>9} {'ph_amax':>8} {'ph_amin':>8} {'a_max/min':>10} {'turn(a)':>7} {'turn(beta)':>10} {'turn(kap)':>9} {'beta_min':>9} {'beta_max':>9} {'gamma':>8} {'A':>5} {'(lam+A)/2':>9} {'<w beta>/<w>':>12} {'N-2<w>':>10}")
for (mu, rho), P in sorted(D.items()):
    k = P['k']; A = 1 + mu; lam = k / mu; T = P['T']; t = P['t']
    a = np.exp(P['log_ahat']); x = P['x']; p = P['p']; y = P['y']; beta = P['beta']; kap = P['kappa']
    w = a * x
    b = beta * a; c = kap * a
    N = avg(c * (lam - p) + b * y + a * x)
    print(f"{mu:5.2f} {rho:5.2f} {t[np.argmin(P['U'])]/T:9.3f} {t[np.argmax(a)]/T:8.3f} {t[np.argmin(a)]/T:8.3f} {a.max()/a.min():10.3e} {nturn(a):7d} {nturn(beta):10d} {nturn(kap):9d} {beta.min():9.4f} {beta.max():9.4f} {P['gamma_beta']:8.4f} {A:5.2f} {(lam+A)/2:9.4f} {avg(w*beta)/avg(w):12.4f} {(N-2*avg(w))/avg(a):+10.3e}")
print()
print("Target C shapes: d = eta - y/(lam-p), e = zeta - x/(lam-p)")
print(f"{'mu':>5} {'rho':>5} {'<d>':>11} {'lam_c':>11} {'d_min':>10} {'d_max':>10} {'frac(d>0)':>9} {'e_min':>10} {'e_max':>10} {'frac(e>0)':>9} {'zeros(d)':>8} {'zeros(e)':>8} {'<e>':>10} {'<d^2>':>10}")
for (mu, rho), P in sorted(D.items()):
    k = P['k']; A = 1 + mu; lam = k / mu
    et = P['y'] / (lam - P['p']); zt = P['x'] / (lam - P['p'])
    d = P['eta'] - et; e = P['zeta'] - zt
    print(f"{mu:5.2f} {rho:5.2f} {avg(d):+11.4e} {P['lam_M']:+11.4e} {d.min():+10.3e} {d.max():+10.3e} {np.mean(d>0):9.3f} {e.min():+10.3e} {e.max():+10.3e} {np.mean(e>0):9.3f} {nzero(d):8d} {nzero(e):8d} {avg(e):+10.3e} {avg(d*d):10.3e}")
