#!/usr/bin/env python3
"""Referee: symbolic verification (sympy) of every algebraic identity used in attractor_theorem.tex and in
NOTES.md 2(b).  Each line prints the simplified residual, which must be 0."""
import sympy as sp

U, p, y, eta, zeta, mu, k, th = sp.symbols('U p y eta zeta mu k theta', real=True)
x = k * sp.exp(U)
f = sp.Matrix([p, k - mu * p - y, x - y])
J = sp.Matrix([[0, 1, 0], [0, -mu, -1], [x, 0, -1]])
ric = {eta: zeta + (mu - 1) * eta - eta ** 2, zeta: x - (1 + eta) * zeta}


def D(F):
    """derivative along the flow of (U,p,y) and the Riccati system"""
    return (sp.diff(F, U) * f[0] + sp.diff(F, p) * f[1] + sp.diff(F, y) * f[2]
            + sp.diff(F, eta) * ric[eta] + sp.diff(F, zeta) * ric[zeta])


def show(name, expr):
    r = sp.simplify(sp.expand(expr))
    print(f"{name:72s} residual = {r}")
    return r == 0

ok = True
n = sp.Matrix([zeta, -eta, -1])
ok &= show("(1.3)  n' + J^T n + (1+eta) n = 0", (n.applyfunc(D) + J.T * n + (1 + eta) * n).norm() ** 2)
O = (n.T * f)[0]
ok &= show("Lemma 4.2(i)  O' + (1+eta) O = 0,  O = n.f", D(O) + (1 + eta) * O)
K = zeta * U - eta * p - y + k
ok &= show("Lemma 4.2(ii) K' + (1+eta)K - k(Ue^U-e^U+1) = 0", D(K) + (1 + eta) * K - k * (U * sp.exp(U) - sp.exp(U) + 1))
W = D(zeta) * p - D(eta) * f[1]
ok &= show("Lemma 4.2(iii) W' + (1+eta)W - x p^2 + eta' O = 0  (so =0 where O=0)", D(W) + (1 + eta) * W - x * p ** 2 + D(eta) * O)
ok &= show("W = n'.f", W - (n.applyfunc(D).T * f)[0])
# two-orbit comparison
U2, p2, y2 = sp.symbols('U2 p2 y2', real=True)
f2 = sp.Matrix([p2, k - mu * p2 - y2, k * sp.exp(U2) - y2])
Delta = zeta * (U2 - U) - eta * (p2 - p) - (y2 - y)
D2 = D(Delta) + sp.diff(Delta, U2) * f2[0] + sp.diff(Delta, p2) * f2[1] + sp.diff(Delta, y2) * f2[2]
r = x * (sp.exp(U2 - U) - 1 - (U2 - U))
ok &= show("Thm 4.1  Delta' + (1+eta)Delta + x1[e^{U2-U1}-1-(U2-U1)] = 0", D2 + (1 + eta) * Delta + r)
# linear system <-> Riccati
H1, H2, H3 = sp.symbols('H1 H2 H3', positive=True)
B = sp.Matrix([[0, 0, x], [1, mu, 0], [0, 1, 1]])
Hd = B * sp.Matrix([H1, H2, H3])
e_ = H2 / H3; z_ = H1 / H3
de = (Hd[1] * H3 - H2 * Hd[2]) / H3 ** 2; dz = (Hd[0] * H3 - H1 * Hd[2]) / H3 ** 2
ok &= show("(3.1) ratios of H' = B H solve Riccati (eta)", de - (z_ + (mu - 1) * e_ - e_ ** 2))
ok &= show("(3.1) ratios of H' = B H solve Riccati (zeta)", dz - (x - (1 + e_) * z_))
# equilibrium
g = sp.symbols('gamma', positive=True)
zz = sp.symbols('z')
etaE = g - 1; zetaE = (g - 1) * (g - mu)
ok &= show("char. poly  z^3+Az^2+mu z+k = (z+g)(z^2+(mu-etaE)z+zetaE)  [k=g(g-1)(g-mu)]",
           zz ** 3 + (1 + mu) * zz ** 2 + mu * zz + g * (g - 1) * (g - mu) - (zz + g) * (zz ** 2 + (mu - etaE) * zz + zetaE))
JE = J.subs(U, 0)
vs = sp.Matrix([1, -g, -g * (g - mu)])
ok &= show("J_E v_s = -gamma v_s", ((JE * vs + g * vs).subs(k, g * (g - 1) * (g - mu))).norm() ** 2)
nE = sp.Matrix([zetaE, -etaE, -1])
ok &= show("J_E^T n_E = -gamma n_E", ((JE.T * nE + g * nE).subs(k, g * (g - 1) * (g - mu))).norm() ** 2)
# Lemma 2.1 coordinates
q = k / mu - p
ok &= show("q' = y - mu q", -f[1] - (y - mu * q))
# Remark 5.1
h = sp.symbols('h')
sol = sp.solve(sp.Eq(-eta * (k - h), k * sp.exp(U) - h), h)[0]
ok &= show("Remark 5.1  h(u,0)-k = k(e^u-1)/(1+eta)", sol - k - k * (sp.exp(U) - 1) / (1 + eta))
# Remark 8.2 log-derivatives
ok &= show("D log(1+eta) = (zeta-mu)/(1+eta) + mu - eta", D(sp.log(1 + eta)) - ((zeta - mu) / (1 + eta) + mu - eta))
ok &= show("D log(zeta+mu eta) = x/(zeta+mu eta) - (1+eta-mu)", D(sp.log(zeta + mu * eta)) - (x / (zeta + mu * eta) - (1 + eta - mu)))
# Claim 5, planar wedge identity with a generic H
Hf = sp.Function('H')
u, v = sp.symbols('u v', real=True)
F = lambda a, b: sp.Matrix([b, k - mu * b - Hf(a, b)])
wedge = lambda a, b: a[0] * b[1] - a[1] * b[0]
lhs = wedge(F(th * u, th * v), F(u, v))
rhs = -v * (th * (Hf(u, v) - k) + k - Hf(th * u, th * v))
ok &= show("Claim 5  F(theta q)^F(q) = -v[theta(H(q)-k)+k-H(theta q)]", lhs - rhs)
# Claim 5, 3D compound identity
C = sp.Matrix([[-mu, 1, 0], [0, -1, 1], [x, 0, -(1 + mu)]])


def conj_wedge(a, b):
    W12 = a[0] * b[1] - a[1] * b[0]; W13 = a[0] * b[2] - a[2] * b[0]; W23 = a[1] * b[2] - a[2] * b[1]
    return sp.Matrix([W12, -W13, -W23])

nu1, nu2, nu3 = [sp.Function(s)(U, p, y) for s in ('nu1', 'nu2', 'nu3')]
nu = sp.Matrix([nu1, nu2, nu3])
om = conj_wedge(f, nu)
Dnu_f = nu.jacobian([U, p, y]) * f
gen = om.applyfunc(D) - C * om - conj_wedge(f, Dnu_f - J * nu)
ok &= show("general: D omega - C omega = conj(f ^ (Dnu f - J nu)), omega = conj(f^nu)", sum(sp.expand(e) ** 2 for e in gen))
nur = sp.Matrix([U, p, y - k])
omr = conj_wedge(f, nur)
spec = omr.applyfunc(D) - C * omr - k * (U * sp.exp(U) - sp.exp(U) + 1) * sp.Matrix([0, p, f[1]])
ok &= show("Claim 5 (3D): nu=(U,p,y-k):  D omega - C omega = k(Ue^U-e^U+1)(0,p,p')", sum(sp.expand(e) ** 2 for e in spec))
# planar wedge transport used in Prop 8.1(ii):  (F ^ Jv)' = div F * (F ^ Jv)
F1, F2 = [sp.Function(s)(u, v) for s in ('F1', 'F2')]
a1, a2 = sp.symbols('a1 a2')
DF = sp.Matrix([[sp.diff(F1, u), sp.diff(F1, v)], [sp.diff(F2, u), sp.diff(F2, v)]])
Fv = sp.Matrix([F1, F2]); Jv = sp.Matrix([a1, a2])
dw = wedge(DF * Fv, Jv) + wedge(Fv, DF * Jv)
ok &= show("Prop 8.1(ii)  d/dt (F^J) = (div F)(F^J)", dw - (DF[0, 0] + DF[1, 1]) * wedge(Fv, Jv))
print("ALL IDENTITIES OK" if ok else "SOME IDENTITY FAILED")
