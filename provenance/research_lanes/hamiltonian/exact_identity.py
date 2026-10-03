#!/usr/bin/env python3
"""Numerical test of the exact 'third-order selection identity' (NOTES.md, Lemma 1.1(d)) on true WH cycles:
   (rho-1)<P^2> = eps^2 [ (1/rho)(<P^4> - 3<Q^2 P^2>) - <Q^3> + 2(1-1/rho)(<Q'^2> - <Q^2>) ],
with P = U_tau, Q = U_tautau, Q' = U_tautautau, tau = sqrt(mu) t, eps = sqrt(mu)/(1+mu).
Also <Q^2> = <P^2> and <(rho e^U - 1)P^2> = 0."""
import sys
import numpy as np
sys.path.insert(0, '/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03/work')
from whlib import find_cycle, cycle_profile, kH

print("# mu    rho    lhs=(rho-1)<P^2>   rhs            rel.err     <Q^2>/<P^2>-1   <(rho e^U-1)P^2>/<P^2>")
for mu, rho in [(1.0, 1.2), (1.0, 2.0), (1.0, 5.0), (0.3, 1.5), (0.1, 1.1), (0.1, 3.0), (0.05, 2.0), (3.0, 2.0)]:
    k = rho * kH(mu)
    c = find_cycle(mu, k, warm=80)
    pr = cycle_profile(mu, k, c, n=200001)
    t, U, p, y, x = pr["t"], pr["U"], pr["p"], pr["y"], pr["x"]
    eps = np.sqrt(mu) / (1 + mu)
    dp = k - mu * p - y                    # p'
    ddp = -mu * dp - (x - y)               # p''
    P = p / np.sqrt(mu); Q = dp / mu; Qd = ddp / mu ** 1.5
    av = lambda f: np.trapezoid(f, t) / (t[-1] - t[0])
    lhs = (rho - 1) * av(P * P)
    rhs = eps ** 2 * ((av(P ** 4) - 3 * av(Q * Q * P * P)) / rho - av(Q ** 3) + 2 * (1 - 1 / rho) * (av(Qd * Qd) - av(Q * Q)))
    print(f"{mu:<6g} {rho:<6g} {lhs:<18.10g} {rhs:<14.10g} {abs(lhs-rhs)/abs(lhs):<11.2e} {av(Q*Q)/av(P*P)-1:<15.2e} {av((rho*np.exp(U)-1)*P*P)/av(P*P):+.2e}")
