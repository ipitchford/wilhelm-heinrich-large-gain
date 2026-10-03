#!/usr/bin/env python3
"""Numerical test of the exact formulas of NOTES 4.4 on true cycles (floating point):
 (X1) (1+eta) p' - eta' p - (eta-mu)(1+eta) p + x - k = 0
 (X2) E1 = V^2/2 + rho*Phi(U),  V = A P/(1+eta) = P/s:   dE1/dtau = sqrt(mu) qh V^2 + rho G V (s^3-1)/s^2
 and the size of the variance <qh^2>, qh = (eta-mu)/mu, against the non-attracting ceiling r/(A(2+mu))."""
import sys
import numpy as np
sys.path.insert(0, '/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03/work')
from whlib import cycle_profile, kH, find_cycle
from ampcycle import amp_cycle
from toda import Phi

print("# mu      rho       r        max|X1|/scale  max|X2 defect|/max|E1'|  <qh^2>    ceiling r/(A(2+mu))  <qh^2>/<G^2>  E1max/E1min")
cases = [(1e-3, None, 1.0), (1e-3, None, 2.5), (1e-2, None, 2.0), (0.1, 1.5, None), (1.0, 2.0, None), (0.05, 3.0, None)]
for mu, rho, a in cases:
    if a is not None:
        c = amp_cycle(mu, a)
        k = c["k"]; rho = c["rho"]
    else:
        k = rho * kH(mu); c = find_cycle(mu, k, warm=80)
    A = 1 + mu
    pr = cycle_profile(mu, k, c, n=200001)
    t, U, p, y, x, eta, zeta = (pr[s] for s in ("t", "U", "p", "y", "x", "eta", "zeta"))
    dp = k - mu * p - y
    deta = zeta + (mu - 1) * eta - eta ** 2
    X1 = (1 + eta) * dp - deta * p - (eta - mu) * (1 + eta) * p + x - k
    sc = np.max(np.abs(x - k))
    eps = np.sqrt(mu) / A
    r = (rho - 1) / eps ** 2
    qh = (eta - mu) / mu
    s = (1 + eta) / A
    P = p / np.sqrt(mu)
    V = P / s
    G = np.expm1(U)
    E1 = V ** 2 / 2 + rho * Phi(U)
    tau = np.sqrt(mu) * t
    dE1_num = np.gradient(E1, tau, edge_order=2)
    dE1_formula = np.sqrt(mu) * qh * V ** 2 + rho * G * V * (s ** 3 - 1) / s ** 2
    av = lambda f: np.trapezoid(f, t) / (t[-1] - t[0])
    print(f"{mu:<8g} {rho:<9.6f} {r:<8.4f} {np.max(np.abs(X1))/sc:<14.2e} {np.max(np.abs(dE1_num-dE1_formula))/np.max(np.abs(dE1_formula)):<24.2e} "
          f"{av(qh**2):<9.4f} {r/(A*(2+mu)):<20.4f} {av(qh**2)/av(G**2):<13.4f} {E1.max()/E1.min():.4f}")
