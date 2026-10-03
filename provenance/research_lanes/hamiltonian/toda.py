#!/usr/bin/env python3
"""Unperturbed exponential-potential (Toda / log-Lotka-Volterra) centre.

H0(U,p) = p^2/2 + Phi(U),  Phi(U) = e^U - U - 1,   U'' = -(e^U - 1)   (rho = 1 normalisation).
For general rho: H = p^2/2 + rho*Phi(U); with p = sqrt(rho) q, s = sqrt(rho) tau this is rho*H0(U,q).

level_data(h) returns period T0 and the time-averages / loop integrals needed by the averaging theory:
  I   = oint p^2 dtau = oint p dU,   J1 = oint p^4 dtau = oint p^3 dU,
  IE  = oint e^U p^2 dtau,  JE = oint e^U p^4 dtau,   G2 = oint (e^U-1)^2 dtau.
Floating point only.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq


def Phi(U):
    return np.expm1(U) - U


def turning_points(h):
    """U_- < 0 < U_+ with Phi(U_pm) = h."""
    up = brentq(lambda u: Phi(u) - h, 0.0, np.log(h + 2.0) + 2.0 + np.sqrt(2 * h), xtol=1e-15, rtol=1e-15)
    um = brentq(lambda u: Phi(u) - h, -(h + 1.0) - 1e-12, 0.0, xtol=1e-15, rtol=1e-15)
    return um, up


def level_data(h, rtol=1e-12, atol=1e-14):
    um, up = turning_points(h)

    def f(t, s):
        U, p = s[0], s[1]
        e = np.exp(U)
        return [p, -(e - 1.0), p * p, p ** 4, e * p * p, e * p ** 4, (e - 1.0) ** 2]

    def ev(t, s):
        return s[1]
    ev.terminal = True
    ev.direction = 1  # p goes negative first (starting at crest), returns to 0 from below at the trough
    s0 = [up, 0.0, 0, 0, 0, 0, 0]
    a = solve_ivp(f, [0, 1e-9], s0, method="DOP853", rtol=rtol, atol=atol)
    b = solve_ivp(f, [1e-9, 50 + 20 * np.sqrt(h + 1)], a.y[:, -1], method="DOP853", rtol=rtol, atol=atol, events=ev)
    assert b.status == 1, "half period not found"
    half = b.t_events[0][0]
    s = b.y_events[0][0]
    T0 = 2 * half
    I, J1, IE, JE, G2 = (2 * s[2], 2 * s[3], 2 * s[4], 2 * s[5], 2 * s[6])
    return dict(h=h, Um=um, Up=up, T0=T0, I=I, J1=J1, IE=IE, JE=JE, G2=G2,
                p2=I / T0, p4=J1 / T0, trough_err=abs(s[0] - um))


def R_of_h(h):
    d = level_data(h)
    return d["J1"] / d["I"]


def h_star(r):
    """Root of r = (2/3) R(h)."""
    f = lambda h: (2.0 / 3.0) * R_of_h(h) - r
    lo, hi = 1e-6, max(4 * r, 1.0)
    while f(hi) < 0:
        hi *= 2
    return brentq(f, lo, hi, xtol=1e-12, rtol=1e-12)


if __name__ == "__main__":
    print("# h      T0        I=oint p^2   IE/I (=1?)   G2/I (=1?)   JE/J1 (=1?)   R=J1/I    R/h     kurt=<p^4>/<p^2>^2")
    for h in [1e-3, 1e-2, 0.1, 0.5, 1, 2, 5, 10, 30, 100, 300, 1000]:
        d = level_data(h)
        kurt = d["p4"] / d["p2"] ** 2
        print(f"{h:<8g} {d['T0']:<9.5f} {d['I']:<12.6g} {d['IE']/d['I']:<12.9f} {d['G2']/d['I']:<12.9f} "
              f"{d['JE']/d['J1']:<12.9f} {d['J1']/d['I']:<9.5g} {d['J1']/d['I']/h:<7.4f} {kurt:.5f}")
