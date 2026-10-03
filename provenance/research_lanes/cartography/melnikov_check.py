#!/usr/bin/env python3
"""Check of the exact identity  oint U_ss^2 ds = oint U_s^2 ds  on every periodic orbit of U_ss = -(e^U - 1)
(first-order Melnikov function of the mu -> 0, rho - 1 = O(mu) Hamiltonian limit has no amplitude selection)."""
import numpy as np
from scipy.integrate import solve_ivp
for h in (0.1, 0.5, 1.0, 2.0, 4.0, 6.0):
    def f(s, z):
        return [z[1], -np.expm1(z[0]), z[1] ** 2, np.expm1(z[0]) ** 2]
    def ev(s, z):
        return z[1]
    ev.direction = -1
    ev.terminal = True
    sol = solve_ivp(f, [0, 1e5], [h, -1e-14, 0, 0], method='DOP853', rtol=1e-13, atol=1e-15, events=ev, first_step=1e-4)
    P = sol.t[-1]
    I1, I2 = sol.y[2, -1], sol.y[3, -1]
    print(f"h={h}: period={P:.8f}  oint U_s^2 = {I1:.12f}  oint U_ss^2 = {I2:.12f}  relative difference = {abs(I1 - I2) / I1:.2e}")
