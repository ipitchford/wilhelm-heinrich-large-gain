#!/usr/bin/env python3
"""Referee check of the algebra behind Theorem C (NOTES 3.1 and Step 0/Step 1 of 4.1).

 (A) slow-manifold coefficients w_j DERIVED by recursion from the invariance equation (not verified from the
     lane's formulas), then compared with the lane's w0..w3;
 (B) the exact equation for D = Q - w used in Step 1;
 (C) eps^0..eps^3 coefficients of X(H~), derived with sigma_h, sigma_U as abstract partials, compared with F3;
 (D) numerical average of F3 over Toda orbits (ODE integration) vs r J0 - (2/3) J1 (quadrature);
 (E) the three averaging identities used to reduce oint F3;
 (F) limit g3(0) = lim M3(Phi(a))/(a G(a)) = pi r  (Step 3);
 (G) divergence route to the multiplier: oint dw/dP dtau on the Toda orbit = r T0 - 2 J0 = M3'(h) (independent
     derivation of the formula in Theorem C(3)), and = -J0 Rsel' at r = Rsel(h).
"""
import sys
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
sys.path.insert(0, '/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03/work/lanes/hamiltonian/verify')
from todaq import level, Phi as PhiN, turning

FAIL = []
U, P, eps, r, h = sp.symbols('U P epsilon r h', real=True)
E = sp.exp(U)
G = E - 1
PhiS = E - U - 1

# (A) recursion: eps (w_U P + w_P w) = -w - eps P - (1+eps^2 r) G
N = 4
ws = [sp.Symbol(f'w{j}') for j in range(N)]
w = [None] * N
w[0] = -G
for j in range(1, N):
    # coefficient of eps^j:  w_j = -[ (w_U P + w_P w)_{j-1} + delta_{j1} P + delta_{j2} r G ]
    conv = sp.diff(w[j - 1], U) * P + sum(sp.diff(w[i], P) * w[j - 1 - i] for i in range(j))
    w[j] = sp.expand(-(conv + (P if j == 1 else 0) + (r * G if j == 2 else 0)))
lane = [-G, G * P, -r * G - E * P**2 + G**2, r * E * P + E * P**3 - 4 * G * E * P - G**2 * P]
dA = [sp.simplify(w[j] - lane[j]) for j in range(N)]
print("(A) derived w_j minus lane's w_j, j=0..3:", dA)
if any(x != 0 for x in dA):
    FAIL.append("A")

# (B) D-equation: eps D' = -(1 + eps w_P) D, with D = Q - w(U,P) for an exact invariant graph
Q = sp.Symbol('Q')
wf = sp.Function('w')(U, P)
rho = 1 + eps**2 * r
Qdot = (-Q - eps * P - rho * G) / eps
inv_rhs = (-wf - eps * P - rho * G) / eps          # = w_U P + w_P w   (invariance)
Ddot = Qdot - sp.diff(wf, U) * P - sp.diff(wf, P) * Q
Ddot = Ddot.subs(sp.diff(wf, U), (inv_rhs - sp.diff(wf, P) * wf) / P)
resB = sp.simplify(Ddot + (1 / eps + sp.diff(wf, P)) * (Q - wf))
print("(B) D' + (1/eps + w_P) D =", resB)
if resB != 0:
    FAIL.append("B")

# (C) modified energy
um = sp.Function('um')
H0 = P**2 / 2 + PhiS
wser = sum(eps**j * w[j] for j in range(N))
sig_U = 2 * (h - PhiS)
sig_h = 2 * (U - um(h))
dh2_dU = r * G - G**2 + sig_U + sig_h * sp.diff(H0, U)
dh2_dP = sig_h * sp.diff(H0, P)
XH = P * (sp.diff(H0, U) + eps**2 * dh2_dU) + wser * (sp.diff(H0, P) + eps * P**2 + eps**2 * dh2_dP)
XH = sp.expand(XH.subs(h, H0))
co = [sp.simplify(XH.coeff(eps, j)) for j in range(4)]
F3 = r * P**2 - 4 * G * E * P**2 + 2 * G * P**2 * (U - um(H0))
print("(C) coefficients eps^0,1,2 of X(H~):", co[:3], "; eps^3 minus F3:", sp.simplify(co[3] - F3))
if any(c != 0 for c in co[:3]) or sp.simplify(co[3] - F3) != 0:
    FAIL.append("C")
# check that sigma defined as int_{um(h)}^U 2(h-Phi) has the partials used
uu = sp.Symbol('u')
sig_h_check = sp.diff(sp.Integral(2 * (h - (sp.exp(uu) - uu - 1)), (uu, um(h), U)), h).doit()
sig_h_check = sp.simplify(sig_h_check.subs(sp.exp(um(h)), h + um(h) + 1))
print("    d sigma/dh from the definition (using Phi(um(h)) = h):", sp.simplify(sig_h_check), " expected 2(U-um(h))")

# (D),(E),(G) numerics on Toda orbits
print("\n(D) oint F3 dtau (ODE)  vs  r J0 - (2/3) J1 (quadrature);  (E) identities;  (G) oint w3_P = r T0 - 2 J0")
print("# h      r     ointF3        rJ0-(2/3)J1   rel.diff   ointGP2/J0   3ointGUP2/J1-1  3ointGEP2/J1-1  ointEP4/J1-1  ointw3P        rT0-2J0       -J0Rsel'(r=Rsel)")
for hh in [0.05, 0.7, 4.0, 25.0, 150.0]:
    d = level(hh)
    um_, up_ = d["um"], d["up"]
    for rr in [0.3, d["Rsel"], 3.0]:
        def f(t, s):
            u, p = s[0], s[1]
            e = np.exp(u); g = e - 1
            F3n = rr * p * p - 4 * g * e * p * p + 2 * g * p * p * (u - um_)
            w3P = rr * e + 3 * e * p * p - 4 * g * e - g * g
            return [p, -g, F3n, g * p * p, g * u * p * p, g * e * p * p, e * p ** 4, w3P, p ** 4, 1.0]
        ev = lambda t, s: s[1]
        ev.terminal = True; ev.direction = 1
        # crest -> trough (half period); all integrands are even under (tau,P) -> (-tau,-P), so double
        sol = solve_ivp(f, [0, 10 * d["T0"]], [up_, 0, 0, 0, 0, 0, 0, 0, 0, 0], method="DOP853", rtol=1e-13, atol=1e-15, events=ev)
        v = 2 * sol.y_events[0][0]
        M3 = rr * d["J0"] - 2.0 / 3.0 * d["J1"]
        scale = max(abs(rr * d["J0"]), abs(2.0 / 3.0 * d["J1"]))
        rel = abs(v[2] - M3) / scale
        e1 = v[3] / d["J0"]; e2 = 3 * v[4] / d["J1"] - 1; e3 = 3 * v[5] / d["J1"] - 1; e4 = v[6] / d["J1"] - 1
        gdiv = rr * d["T0"] - 2 * d["J0"]
        print(f"{hh:<7g} {rr:<5.3g} {v[2]:<13.6e} {M3:<13.6e} {rel:<10.1e} {e1:<12.1e} {e2:<15.1e} {e3:<15.1e} {e4:<13.1e} "
              f"{v[7]:<14.7e} {gdiv:<13.7e} {(-d['J0']*d['Rselp'] if abs(rr-d['Rsel'])<1e-14 else float('nan')):.7e}")
        if rel > 1e-9 or max(abs(e1), abs(e2), abs(e3), abs(e4)) > 1e-9 or abs(v[7] - gdiv) > 1e-8 * max(1, abs(gdiv)):
            FAIL.append(f"D/E/G h={hh} r={rr}")
        if abs(rr - d["Rsel"]) < 1e-14 and abs(gdiv + d["J0"] * d["Rselp"]) > 1e-9 * d["J0"]:
            FAIL.append(f"G-consistency h={hh}")

# (F)
print("\n(F) g3(a) = M3(Phi(a))/(a G(a)) as a -> 0, r = 1 (claim: -> pi r)")
for a in [0.3, 0.1, 0.03, 0.01]:
    d = level(PhiN(a))
    g3 = (1.0 * d["J0"] - 2.0 / 3.0 * d["J1"]) / (a * np.expm1(a))
    print(f"  a={a:<5g} g3={g3:.6f}   g3/pi={g3/np.pi:.6f}")

print("\nFAILURES:", FAIL if FAIL else "none")
