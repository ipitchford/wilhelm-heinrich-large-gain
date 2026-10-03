#!/usr/bin/env python3
"""Referee check of Lemma 1 (NOTES 1): exact identities on every cycle, at the extreme grid
mu in {0.02, 0.3, 1} x rho in {1.01, 1.5, 4, 12}, plus one mu > 1 case (the lemma is stated for all parameters).
Own integrator (whref.py); period-averages by running integrals inside the ODE (robust for spiky cycles).

Each identity is written  LHS - RHS = 0  and the defect is normalised by the largest term entering it.
 (a) <e^U> = 1
 (b) <Q^2> = <P^2>
 (c) <Q^2> = rho <e^U P^2>
 (d) (rho-1)<P^2> = eps<P^3> - 2 eps <Q^2 P>
 (e) <e^U P^3> = -(2 eps/rho)(<Q'^2> - <Q^2>)
 (f) (rho-1)<P^2> = eps^2[(1/rho)(<P^4> - 3<Q^2P^2>) - <Q^3> + 2(1-1/rho)(<Q'^2> - <Q^2>)]
 (g) the same in the brief's variables:
     (k-k_H) A <p^2> = (1/rho)(mu<p^4> - 3<p'^2 p^2>) - <p'^3> + 2(1-1/rho)(<p''^2> - mu<p'^2>)
"""
import sys, time
import numpy as np
sys.path.insert(0, '/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03/work/lanes/hamiltonian/verify')
from whref import cycle, averages

grid = [(mu, rho) for mu in (0.02, 0.3, 1.0) for rho in (1.01, 1.5, 4.0, 12.0)] + [(3.0, 2.0), (0.02, 1.003)]
print("# mu    rho    U_max    U_min       T          m_t         res      |  (a)       (b)       (c)       (d)       (e)       (f)       (g)      closure")
worst = 0.0
for mu, rho in grid:
    t0 = time.time()
    eps = np.sqrt(mu) / (1 + mu)
    A = 1 + mu
    kH = mu * A
    k = rho * kH
    r = (rho - 1) / eps ** 2
    # near Hopf start from the small-amplitude prediction a^2 ~ 2r; otherwise generic
    guess = None
    if rho < 1.2:
        a0 = min(np.sqrt(2 * r), 1.5)
        guess = (a0, k + mu * rho * np.expm1(a0))
    try:
        c = cycle(mu, rho, guess=guess, warm=40 if rho < 1.2 else 80)
    except Exception as e:
        print(f"{mu:<6g} {rho:<6g} CYCLE FAILED: {e}")
        continue
    av = averages(mu, k, c)
    P2, Q2 = av["P2"], av["Q2"]
    d = {}
    d["a"] = abs(av["E"] - 1)
    d["b"] = abs(Q2 - P2) / P2
    d["c"] = abs(Q2 - rho * av["EP2"]) / P2
    t1, t2, t3 = (rho - 1) * P2, eps * av["P3"], -2 * eps * av["Q2P"]
    d["d"] = abs(t1 - t2 - t3) / max(abs(t1), abs(t2), abs(t3))
    t1, t2 = av["EP3"], -(2 * eps / rho) * (av["R2"] - Q2)
    d["e"] = abs(t1 - t2) / max(abs(t1), abs(t2), 1e-300)
    terms = [(rho - 1) * P2, eps ** 2 / rho * av["P4"], -3 * eps ** 2 / rho * av["Q2P2"], -eps ** 2 * av["Q3"],
             2 * eps ** 2 * (1 - 1 / rho) * av["R2"], -2 * eps ** 2 * (1 - 1 / rho) * Q2]
    d["f"] = abs(terms[0] - sum(terms[1:])) / max(abs(x) for x in terms)
    tb = [(k - kH) * A * av["p2"], mu * av["p4"] / rho, -3 * av["dp2p2"] / rho, -av["dp3"],
          2 * (1 - 1 / rho) * av["ddp2"], -2 * (1 - 1 / rho) * mu * av["dp2"]]
    d["g"] = abs(tb[0] - sum(tb[1:])) / max(abs(x) for x in tb)
    worst = max(worst, max(d.values()))
    print(f"{mu:<6g} {rho:<6g} {c['U0']:<8.4f} {c['Umin']:<11.4g} {c['T']:<10.4f} {c['m_t']:<11.4e} {c['res']:<8.1e} |  "
          + " ".join(f"{d[kk]:<9.1e}" for kk in "abcdefg") + f" {av['closure']:.1e}   [{time.time()-t0:.0f}s]", flush=True)
print("worst normalised defect over all identities and cycles:", f"{worst:.2e}")
