#!/usr/bin/env python3
"""Exact (sympy) check of the third-order normal form used in NOTES.md, Section 2.

Fast-slow system (tau-time, eps = 1/sigma, rho = 1 + eps^2 r, G = e^U - 1, E = e^U):
    U' = P,  P' = Q,  eps Q' = -Q - eps P - rho G(U).
(1) slow-manifold series Q = w0 + eps w1 + eps^2 w2 + eps^3 w3 solves the invariance equation to O(eps^4);
(2) with H0 = P^2/2 + E - U - 1, h1 = P^3/3, h2 = r*Phi - Psi + s, s = sigma(U,H0),
    sigma(U,h) = int_{um(h)}^U 2(h - Phi(u)) du,  the reduced field X satisfies
    X(H0 + eps h1 + eps^2 h2) = eps^3 F3 + O(eps^4),  F3 = r P^2 - 4 G E P^2 + 2 G P^2 (U - um(H0));
(3) first-order exact identity for the FULL third-order equation (no reduction):
    U''*(eps U''' + U'' + eps U' + rho G) = d/dtau(...) + [U''^2 - rho E U'^2]  -> <U''^2> = rho <e^U U'^2>.
"""
import sympy as sp

U, P, eps, r, h = sp.symbols('U P epsilon r h', real=True)
E = sp.exp(U); G = E - 1
Phi = E - U - 1
rho = 1 + eps**2 * r

# (1) slow manifold
w0 = -G
w1 = G * P
w2 = -r * G - E * P**2 + G**2
w3 = r * E * P + E * P**3 - 4 * G * E * P - G**2 * P
w = w0 + eps * w1 + eps**2 * w2 + eps**3 * w3
inv = eps * (sp.diff(w, U) * P + sp.diff(w, P) * w) + w + eps * P + rho * G
inv = sp.expand(inv)
low = [sp.simplify(inv.coeff(eps, j)) for j in range(4)]
print("(1) invariance residual coefficients eps^0..eps^3:", low)
assert all(c == 0 for c in low)

# (2) modified energy.  um = U_-(h) abstract; sigma_h = 2(U - um), sigma_U = 2(h - Phi)
um = sp.Function('um')
H0 = P**2 / 2 + Phi
Psi_prime = G**2                      # Psi'(U) = G^2
# represent h2 through its partial derivatives only (that is all X needs):
sig_U = 2 * (h - Phi)                 # d sigma/dU at fixed h
sig_h = 2 * (U - um(h))               # d sigma/dh
dh2_dU = r * G - Psi_prime + (sig_U + sig_h * sp.diff(H0, U))
dh2_dP = sig_h * sp.diff(H0, P)
dh1_dU, dh1_dP = 0, P**2
dH0_dU, dH0_dP = sp.diff(H0, U), sp.diff(H0, P)
XU = P
XP = w
tot = (XU * (dH0_dU + eps * dh1_dU + eps**2 * dh2_dU) + XP * (dH0_dP + eps * dh1_dP + eps**2 * dh2_dP))
tot = sp.expand(tot.subs(h, H0))
c = [sp.simplify(tot.coeff(eps, j)) for j in range(4)]
F3 = r * P**2 - 4 * G * E * P**2 + 2 * G * P**2 * (U - um(H0))
print("(2) eps^0, eps^1, eps^2 coefficients of X(H~):", c[0], c[1], c[2])
print("    eps^3 coefficient minus F3:", sp.simplify(c[3] - F3))
assert c[0] == 0 and c[1] == 0 and c[2] == 0 and sp.simplify(c[3] - F3) == 0

# (3) exact identity on the full equation
t = sp.symbols('tau')
u = sp.Function('u')(t)
eq = eps * u.diff(t, 3) + u.diff(t, 2) + eps * u.diff(t) + (1 + eps**2 * r) * (sp.exp(u) - 1)
lhs = u.diff(t, 2) * eq
exact = sp.diff(eps * u.diff(t, 2)**2 / 2 + eps * u.diff(t)**2 / 2 + (1 + eps**2 * r) * (sp.exp(u) - 1) * u.diff(t), t)
rest = sp.simplify(sp.expand(lhs - exact - (u.diff(t, 2)**2 - (1 + eps**2 * r) * sp.exp(u) * u.diff(t)**2)))
print("(3) U''*eq - d/dtau[...] - (U''^2 - rho e^U U'^2) =", rest)
assert rest == 0
print("ALL EXACT CHECKS PASSED")
