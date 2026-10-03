#!/usr/bin/env python3
"""Fast Toda level integrals by Gauss-Kronrod quadrature in U (endpoint substitution U = U_pm -/+ s^2).
Referee's own; cross-validated against the mpmath tanh-sinh values of v1_theoremA.py in v3."""
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq


def Phi(u):
    return np.expm1(u) - u


def turning(h):
    up = brentq(lambda u: Phi(u) - h, 0.0, np.log(h + 2.0) + 2.0 + np.sqrt(2 * h), xtol=1e-15, rtol=1e-15)
    um = brentq(lambda u: Phi(u) - h, -(h + 1.0), 0.0, xtol=1e-15, rtol=1e-15)
    return um, up


def level(h):
    um, up = turning(h)

    def P(u):
        return np.sqrt(max(2.0 * (h - Phi(u)), 0.0))

    def both(g):
        # int_{um}^{0} g(U) dU with U = um + s^2, and int_0^{up} with U = up - s^2
        a = quad(lambda s: 2 * s * g(um + s * s), 0.0, np.sqrt(-um), epsabs=1e-13, epsrel=1e-13, limit=400)[0]
        b = quad(lambda s: 2 * s * g(up - s * s), 0.0, np.sqrt(up), epsabs=1e-13, epsrel=1e-13, limit=400)[0]
        return 2 * (a + b)
    T0 = both(lambda u: 1.0 / P(u) if P(u) > 0 else 0.0)
    J0 = both(P)
    J1 = both(lambda u: P(u) ** 3)
    kurt = (J1 / T0) / (J0 / T0) ** 2
    return dict(h=h, um=um, up=up, T0=T0, J0=J0, J1=J1, Rsel=2.0 / 3.0 * J1 / J0, kurt=kurt,
                Rselp=2.0 / 3.0 * (3.0 - kurt))


def hstar(r):
    return brentq(lambda h: level(h)["Rsel"] - r, 1e-8, 1.3 * r + 1e-6, xtol=1e-13, rtol=1e-13)
