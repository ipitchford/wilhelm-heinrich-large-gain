#!/usr/bin/env python3
"""Item 3: the crest return map P on the section {p=0, p'<0} = {p=0, y>k}, coordinates (U,y).
 (a) is DP entrywise positive away from cycles (global monotonicity on {x>y>k}, and on y>k)?
 (b) sublinearity relative to E:  P(E+theta(z-E)) >= E+theta(P(z)-E)?  (componentwise), theta in (0,1)
 (c) second differences: sign of D^2P[h,h] components for h = e1, e2, e1+e2, radial.
"""
import numpy as np, warnings, sys
from mono_lib import *
warnings.filterwarnings('ignore')


def P(mu, k, U0, y0, deriv=False):
    U1, y1, T, M = crest_return(mu, k, U0, y0, with_derivative=deriv)
    return np.array([U1, y1]), T, M


def study(mu, rho, ngrid=7):
    k = rho * kH(mu)
    c = find_cycle(mu, k)
    Us, ys = c['U0'], c['y0']
    print(f"=== mu={mu} rho={rho} k={k:.4f} crest U*={Us:.4f} y*={ys:.4f} (x*={k*np.exp(Us):.4f}) m_t={c['m_t']:.5f}")
    # grid of section points: U in (0.15..1.6)*U*, y-k in (0.15..1.6)*(y*-k), flag x>y
    fr = np.linspace(0.15, 1.75, ngrid)
    negM = 0; tot = 0; worst = None
    sub_viol = 0; sub_tot = 0
    conc = {('e1', 0): [], ('e1', 1): [], ('e2', 0): [], ('e2', 1): [], ('rad', 0): [], ('rad', 1): [], ('e1+e2', 0): [], ('e1+e2', 1): []}
    for fu in fr:
        for fy in fr:
            U0 = fu * Us; y0 = k + fy * (ys - k)
            inreg = k * np.exp(U0) > y0
            try:
                z1, T, M = P(mu, k, U0, y0, deriv=True)
            except Exception as e:
                continue
            tot += 1
            if (M < 0).any():
                negM += 1
                if worst is None or M.min() < worst[0]:
                    worst = (M.min(), fu, fy, inreg, M.copy())
            # sublinearity relative to E=(0,k): theta=0.9
            for th in (0.5, 0.9):
                try:
                    zt, _, _ = P(mu, k, th * U0, k + th * (y0 - k))
                except Exception:
                    continue
                lhs = zt - np.array([0, k]); rhs = th * (z1 - np.array([0, k]))
                sub_tot += 1
                if (lhs < rhs - 1e-9).any():
                    sub_viol += 1
            # second differences
            eps = 1e-3
            for name, hvec in (('e1', np.array([Us, 0.0])), ('e2', np.array([0.0, ys - k])),
                               ('rad', np.array([U0, y0 - k])), ('e1+e2', np.array([Us, ys - k]))):
                try:
                    zp, _, _ = P(mu, k, U0 + eps * hvec[0], y0 + eps * hvec[1])
                    zm, _, _ = P(mu, k, U0 - eps * hvec[0], y0 - eps * hvec[1])
                except Exception:
                    continue
                d2 = (zp - 2 * z1 + zm) / eps ** 2
                conc[(name, 0)].append(d2[0]); conc[(name, 1)].append(d2[1])
    print(f"  DP entrywise>0 fails at {negM}/{tot} grid points; worst: {None if worst is None else (round(worst[0],4), round(worst[1],2), round(worst[2],2), 'x>y' if worst[3] else 'x<y')}")
    if worst is not None:
        print("   worst M =", np.round(worst[4], 4).tolist())
    print(f"  radial sublinearity (rel. E) violated at {sub_viol}/{sub_tot}")
    for key, v in conc.items():
        v = np.array(v)
        print(f"  D2P[{key[0]},{key[0]}]_{key[1]+1}: neg {np.mean(v<0):.2f} pos {np.mean(v>0):.2f}  range [{v.min():.3g},{v.max():.3g}]")


if __name__ == '__main__':
    for mu, rho in [(1.0, 1.3), (1.0, 2.5), (0.2, 1.5), (0.2, 1.1), (3.0, 1.5)]:
        study(mu, rho)
