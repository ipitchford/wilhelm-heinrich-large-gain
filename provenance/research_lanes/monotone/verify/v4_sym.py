#!/usr/bin/env python3
"""Referee re-derivation (sympy, exact) of every algebraic identity used in ds_theorem.tex and NOTES R1,R6,R7.
Written from the equations of motion only; does not import the lane's sym_checks.py."""
import sympy as sp
t = sp.symbols('t')
mu, k = sp.symbols('mu k', positive=True)
A = 1 + mu
U = sp.Function('U')(t); p = sp.Function('p')(t); y = sp.Function('y')(t)
a = sp.Function('a')(t); b = sp.Function('b')(t); c = sp.Function('c')(t)
x = k * sp.exp(U)
eqs = {U.diff(t): p, p.diff(t): k - mu * p - y, y.diff(t): x - y,
       a.diff(t): -b, b.diff(t): -mu * b - c, c.diff(t): -x * a - c}


def D(e, n=1):
    for _ in range(n):
        e = sp.diff(e, t).subs(eqs)
    return sp.expand(e)


def chk(name, e):
    e = sp.simplify(e)
    print(f"{'OK ' if e == 0 else 'FAIL'} {name}" + ('' if e == 0 else f'   residual: {e}'))
    return e == 0


ok = True
pp, yp = D(p), D(y)
J = sp.Matrix([[0, 1, 0], [0, -mu, -1], [x, 0, -1]])
# e_s=(a,-b,-c) solves v' = J v
es = sp.Matrix([a, -b, -c])
ok &= chk("e_s=(a,-b,-c) solves v'=Jv given a'=-b, b'=-mu b-c, c'=-x a-c", (es.applyfunc(D) - J * es).norm() ** 2)
r, s = b / a, c / a
ok &= chk("r' = r^2 - mu r - s", D(r) - (r ** 2 - mu * r - s))
ok &= chk("s' = s r - s - x", D(s) - (s * r - s - x))
ok &= chk("r'' + (A-3r) r' + r(r-1)(r-mu) = x", D(r, 2) + (A - 3 * r) * D(r) + r * (r - 1) * (r - mu) - x)
ct = c + (mu - 1) * b
ok &= chk("c~ = -(b'+b)", ct + D(b) + b)
ok &= chk("c~' = -x a - mu c~", D(ct) - (-x * a - mu * ct))
ok &= chk("r' = r^2 - r - c~/a", D(r) - (r ** 2 - r - ct / a))
ok &= chk("L3[a] = a''' + A a'' + mu a' + x a = 0", D(a, 3) + A * D(a, 2) + mu * D(a) + x * a)
# weak adjoint
zp = sp.Matrix([p, pp, yp])
Psi = sp.exp(A * t) * zp.cross(es)
dPsi = Psi.applyfunc(lambda e: sp.diff(e, t).subs(eqs))
ok &= chk("Psi' = -J^T Psi for Psi = e^{At} z' x e_s", sp.simplify((dPsi + J.T * Psi).norm() ** 2))
formula = sp.exp(A * t) * a * sp.Matrix([r * yp - s * pp, yp + s * p, -(pp + r * p)])
ok &= chk("Psi = e^{At} a (r y' - s p', y' + s p, -(p' + r p))", sp.simplify((Psi - formula).norm() ** 2))
ok &= chk("Psi . z' = 0", sp.simplify(Psi.dot(zp)))
ok &= chk("Psi . e_s = 0", sp.simplify(Psi.dot(es)))
# chord identity: second solution
U2 = sp.Function('U2')(t); p2 = sp.Function('p2')(t); y2 = sp.Function('y2')(t)
eqs2 = dict(eqs); eqs2.update({U2.diff(t): p2, p2.diff(t): k - mu * p2 - y2, y2.diff(t): k * sp.exp(U2) - y2})
d = sp.Matrix([U2 - U, p2 - p, y2 - y])
B = Psi.dot(d)
dB = sp.diff(B, t).subs(eqs2)
ok &= chk("(Psi.d)' = Psi_3 x (e^{dU} - 1 - dU)", sp.simplify(dB - Psi[2] * x * (sp.exp(U2 - U) - 1 - (U2 - U))))
zpp = zp.applyfunc(D)
ok &= chk("(Psi.z'')' = Psi_3 x p^2", sp.simplify(sp.diff(Psi.dot(zpp), t).subs(eqs) - Psi[2] * x * p ** 2))
# theta and the three crossing derivatives
theta = -pp / p
ok &= chk("theta' = theta^2 - mu theta + y'/p", D(theta) - (theta ** 2 - mu * theta + yp / p))
ok &= chk("(theta - r)' - [(theta-r)(theta+r-mu) + y'/p + s] = 0  (so = y'/p + s at theta=r)",
          D(theta - r) - ((theta - r) * (theta + r - mu) + yp / p + s))
psi_p = pp + r * p
g2 = yp / p + s
# g2' = s psi_p / p + g2 * (something regular): exhibit the cofactor
ok &= chk("g2' = s psi_p/p - g2 (1 + p'/p)   (so g2' = s psi_p/p on g2=0)", sp.simplify(D(g2) - s * psi_p / p + g2 * (1 + pp / p)))
g1 = yp / pp - s / r
ok &= chk("g1' = x psi_p/(r p') + g1 (mu - 1 + y'/p' + s/r)   (so g1' = x psi_p/(r p') on g1=0)",
          sp.simplify(D(g1) - x * psi_p / (r * pp) - g1 * (mu - 1 + yp / pp + s / r)))
# Lemma 5 curvature
ok &= chk("p'^2 - p p'' = p'^2 + mu p p' + p y'", pp ** 2 - p * D(p, 2) - (pp ** 2 + mu * p * pp + p * yp))
ok &= chk("at p'=-r p:  p'^2 - p p'' = p^2 r (r-mu) + p y'", sp.simplify((pp ** 2 + mu * p * pp + p * yp).subs(y, k - mu * p + r * p) - (p ** 2 * r * (r - mu) + p * (x - (k - mu * p + r * p)))))
# Prop 8
phi = sp.Function('phi')(t)
om = sp.Matrix([phi, phi.diff(t), -phi.diff(t, 2) - mu * phi.diff(t)])
L3 = phi.diff(t, 3) + A * phi.diff(t, 2) + mu * phi.diff(t) + x * phi
ok &= chk("omega' - J omega = -L3[phi] e3", sp.simplify((om.diff(t) - J * om + sp.Matrix([0, 0, L3])).norm() ** 2))
ok &= chk("L3[1] = x", L3.subs(phi, 1).doit() - x)
ok &= chk("L3[p] = 0", D(p, 3) + A * D(p, 2) + mu * D(p) + x * p)
ok &= chk("L3[p'] = -x p^2", D(pp, 3) + A * D(pp, 2) + mu * D(pp) + x * pp + x * p ** 2)
ok &= chk("L3[U] = k(1 - e^U + U e^U)", D(U, 3) + A * D(U, 2) + mu * D(U) + x * U - k * (1 - sp.exp(U) + U * sp.exp(U)))
ok &= chk("L3[y] = x(k + p^2)", D(y, 3) + A * D(y, 2) + mu * D(y) + x * y - x * (k + p ** 2))
ok &= chk("p'' + r p' = (r-mu) p' - y'", D(p, 2) + r * pp - ((r - mu) * pp - yp))
# compound system (R7 / v2.0 Prop 4.4)
v = sp.Matrix(sp.symbols('v1:4')); u = sp.Matrix(sp.symbols('u1:4'))
w = lambda v, u: sp.Matrix([v[1] * u[0] - v[0] * u[1], v[0] * u[2] - v[2] * u[0], v[1] * u[2] - v[2] * u[1]])
X = sp.symbols('X', positive=True)
JX = sp.Matrix([[0, 1, 0], [0, -mu, -1], [X, 0, -1]])
C = sp.Matrix([[-mu, 1, 0], [0, -1, 1], [X, 0, -A]])
ok &= chk("compound: d/dt (v^u) = C (v^u)", sp.simplify((w(JX * v, u) + w(v, JX * u) - C * w(v, u)).norm() ** 2))
dU, dy, Y, kk = sp.symbols('dU dy Y kk')
ok &= chk("wedge of (dU,0,dy) with z'=(0,k-y,x-y) is ((y-k)dU,(x-y)dU,(y-k)dy)",
          (w(sp.Matrix([dU, 0, dy]), sp.Matrix([0, kk - Y, X - Y])) - sp.Matrix([(Y - kk) * dU, (X - Y) * dU, (Y - kk) * dy])).norm() ** 2)
# cross-product identity for general M
Msym = sp.Matrix(3, 3, sp.symbols('m1:10'))
ok &= chk("(Mu) x v + u x (Mv) = (tr M - M^T)(u x v)", sp.simplify(((Msym * u).cross(v) + u.cross(Msym * v) - (Msym.trace() * sp.eye(3) - Msym.T) * u.cross(v)).norm() ** 2))
# adjoint g-equation mirror
eta = sp.Function('eta')(t); zeta = sp.Function('zeta')(t)
eq3 = {eta.diff(t): zeta + (mu - 1) * eta - eta ** 2, zeta.diff(t): x - (1 + eta) * zeta, U.diff(t): p, p.diff(t): k - mu * p - y, y.diff(t): x - y}
g = 1 + eta
g1_ = sp.diff(g, t).subs(eq3); g2_ = sp.diff(g1_, t).subs(eq3)
ok &= chk("g=1+eta: g'' + (3g-A) g' + g(g-1)(g-mu) = x", sp.simplify(g2_ + (3 * g - A) * g1_ + g * (g - 1) * (g - mu) - x))
print('ALL OK' if ok else 'SOME FAILED')
