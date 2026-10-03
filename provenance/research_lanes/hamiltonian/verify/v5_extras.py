#!/usr/bin/env python3
"""Referee extras.
 (X1) Theorem A sharpness margin: 9/5 - kurt for asymptotically linear convex potentials at large h (must stay > 0).
 (X2) Lemma 1(f) evaluated on Toda orbits: <P^4> - 3<Q^2P^2> - <Q^3> = (2/3)<P^4>  (independent derivation of the
      selection rule r = Rsel(h) straight from the exact identity, bypassing the modified-energy computation).
 (X3) Theorem C's asymptotic formula on the brief's extreme grid mu in {0.02,0.3,1} x rho in {1.01,1.5,4,12}
      (recorded for information: the theorem is asymptotic in mu at fixed r and claims nothing at these points).
"""
import sys
import numpy as np
import mpmath as mp
from scipy.integrate import solve_ivp
sys.path.insert(0, '/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03/work/lanes/hamiltonian/verify')
from todaq import level, hstar, Phi
from whref import cycle

mp.mp.dps = 30


def turning(PhiF, h, side):
    a, b = mp.mpf(0), mp.mpf(side)
    while PhiF(b) < h:
        b *= 2
    for _ in range(200):
        m = (a + b) / 2
        if PhiF(m) < h:
            a = m
        else:
            b = m
    return (a + b) / 2


def mpkurt(PhiF, h):
    h = mp.mpf(h)
    um, up = turning(PhiF, h, -1), turning(PhiF, h, +1)
    P = lambda u: mp.sqrt(max(2 * (h - PhiF(u)), mp.mpf(0)))
    integ = lambda g: 2 * (mp.quad(g, mp.linspace(um, 0, 9)) + mp.quad(g, mp.linspace(0, up, 9)))
    T0 = integ(lambda u: 1 / P(u) if P(u) > 0 else mp.mpf(0))
    J0 = integ(P)
    J1 = integ(lambda u: P(u) ** 3)
    return (J1 / T0) / (J0 / T0) ** 2


print("(X1) 9/5 - kurt (30-digit quadrature)")
pots = {
    "sqrt(1+U^2)-1": lambda u: mp.sqrt(1 + u * u) - 1,
    "log cosh U": lambda u: mp.log(mp.cosh(u)) if abs(u) < 30 else abs(u) - mp.log(2),
    "e^U-U-1": lambda u: mp.expm1(u) - u,
}
for name, F in pots.items():
    print(f"  {name:<16}", "  ".join(f"h={h}: {mp.nstr(mp.mpf(9)/5 - mpkurt(F, h), 6)}" for h in ['1e2', '1e4', '1e6']))

print("\n(X2) Lemma 1(f) on Toda orbits (rho = 1): [<P^4> - 3<Q^2 P^2> - <Q^3>] / [(2/3)<P^4>] - 1")
for hh in [0.05, 1.0, 10.0, 100.0]:
    d = level(hh)

    def f(t, s):
        u, p = s[0], s[1]
        g = np.expm1(u)
        q = -g
        return [p, q, p ** 4, q * q * p * p, q ** 3]
    ev = lambda t, s: s[1]
    ev.terminal = True; ev.direction = 1
    sol = solve_ivp(f, [0, 10 * d["T0"]], [d["up"], 0, 0, 0, 0], method="DOP853", rtol=1e-13, atol=1e-15, events=ev)
    v = sol.y_events[0][0]
    print(f"  h={hh:<6g} {(v[2] - 3 * v[3] - v[4]) / (2.0 / 3.0 * v[2]) - 1:+.2e}")

print("\n(X3) Theorem C leading-order formula on the brief's grid (information only)")
print("# mu    rho    eps      r=(rho-1)/eps^2   h*(r)      h_crest    logmt(true)    -eps^3 J0 Rsel'(h*)   ratio    m_t<1")
for mu in (0.02, 0.3, 1.0):
    for rho in (1.01, 1.5, 4.0, 12.0):
        eps = np.sqrt(mu) / (1 + mu)
        r = (rho - 1) / eps ** 2
        k = rho * mu * (1 + mu)
        guess = None
        if rho < 1.2:
            a0 = min(np.sqrt(2 * r), 1.5)
            guess = (a0, k + mu * rho * np.expm1(a0))
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            c = cycle(mu, rho, guess=guess, warm=40 if rho < 1.2 else 80)
            hs = hstar(r)
            d = level(hs)
        pred = -eps ** 3 * d["J0"] * d["Rselp"]
        print(f"{mu:<6g} {rho:<6g} {eps:<8.5f} {r:<17.4f} {hs:<10.4f} {Phi(c['U0']):<10.4f} {c['logmt']:<14.6e} {pred:<21.6e} "
              f"{c['logmt']/pred:<8.4f} {c['m_t'] < 1}")
