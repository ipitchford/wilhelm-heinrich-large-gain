#!/usr/bin/env python3
"""Refutation test of the first-order (brief) premise: at FIXED rho > 1 the cycle energy does not
converge to a finite root of M(h) as mu -> 0; it diverges like 1/eps^2 ~ 1/mu.
Also checks the exact identity <(rho e^U - 1) p^2> = 0 on true cycles."""
import sys
import numpy as np
sys.path.insert(0, '/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03/work')
from whlib import find_cycle, cycle_profile, kH
from toda import Phi

print("# rho   mu      eps      U_max     U_min       h=Phi(U_max)  h*eps^2    mu*U_min   T*mu      m_t          <(rho e^U-1)p^2>/<p^2>")
for rho in [1.5, 3.0]:
    guess = None
    for mu in [0.3, 0.1, 0.03, 0.01, 0.003, 0.001]:
        eps = np.sqrt(mu) / (1 + mu)
        k = rho * kH(mu)
        try:
            c = find_cycle(mu, k, warm=80)
            pr = cycle_profile(mu, k, c, n=20001)
        except Exception as e:
            print(rho, mu, "FAILED", e); continue
        U, p, t = pr["U"], pr["p"], pr["t"]
        ident = np.trapezoid((rho * np.exp(U) - 1) * p * p, t) / np.trapezoid(p * p, t)
        h = Phi(c["U0"])
        print(f"{rho:<5g} {mu:<7g} {eps:<8.5f} {c['U0']:<9.4f} {U.min():<11.3f} {h:<13.4f} {h*eps**2:<10.5f} "
              f"{mu*U.min():<10.5f} {c['T']*mu:<9.4f} {c['m_t']:<12.5e} {ident:+.2e}  res={c['residual']:.1e}", flush=True)
