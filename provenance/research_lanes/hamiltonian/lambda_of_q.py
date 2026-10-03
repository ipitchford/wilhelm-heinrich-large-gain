#!/usr/bin/env python3
"""Singular relaxation limit (mu -> 0, rho fixed): fixed points of F_lam = G_lam o S satisfy S(q) = lam*g(q/lam),
g(a) = s(a) - a with (e^s - 1)/s = 1/(1-a).  They are parametrised by the pulse amplitude q:
lam = Lambda(q) := q / a(q),  g(a)/a = S(q)/q.  Uniqueness for every lam > 1  <=>  Lambda strictly increasing.
Floating-point scan (not a proof)."""
import numpy as np
from scipy.optimize import brentq
from relax_map import pulse


def g(a):
    if a < 1e-8:
        return a + (2.0 / 3.0) * a * a
    f = lambda s: np.expm1(s) / s - 1.0 / (1.0 - a)
    hi = 1.0
    while f(hi) < 0:
        hi *= 2
    s = brentq(f, 1e-12, hi, xtol=1e-15, rtol=1e-14)
    return s - a


def a_of(ratio):
    """solve g(a)/a = ratio for a in (0,1)."""
    f = lambda a: g(a) / a - ratio
    return brentq(f, 1e-7, 1 - 1e-13, xtol=1e-15, rtol=1e-13)


if __name__ == "__main__":
    print("# check g(a)/a increasing:", bool(np.all(np.diff([g(a) / a for a in np.linspace(0.001, 0.999, 500)]) > 0)))
    qs = np.concatenate([np.geomspace(0.02, 1.0, 18)[:-1], np.geomspace(1.0, 12.0, 16)])
    prev = None
    mono = True
    print("# q        S(q)        S/q       a(q)=q/Lambda   Lambda(q)   dLambda>0")
    for q in qs:
        S = pulse(q)
        a = a_of(S / q)
        Lam = q / a
        ok = (prev is None) or (Lam > prev)
        mono &= ok
        print(f"{q:<9.5f} {S:<11.6f} {S/q:<9.5f} {a:<15.9f} {Lam:<11.6f} {ok}")
        prev = Lam
    print("# Lambda monotone increasing on grid:", mono)
