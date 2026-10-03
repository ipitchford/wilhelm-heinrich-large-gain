#!/usr/bin/env python3
"""Exact (sympy) checks of the algebraic identities used in lane 'monotone'.
State: U,p,y ; x=k e^U.  Vector field F=(p, k-mu p-y, x-y).  D = derivative along the flow.
Strong-stable ratios r=b/a, s=c/a with   r' = r^2 - mu r - s,  s' = s r - s - x
(e_s=(a,-b,-c): a'=-b, b'=-mu b-c, c'=-x a-c).
"""
import sympy as sp

U, p, y, k, mu, r, s, a = sp.symbols('U p y k mu r s a', real=True)
A = 1 + mu
x = k * sp.exp(U)
F = sp.Matrix([p, k - mu * p - y, x - y])
J = F.jacobian([U, p, y])


def D(f):          # derivative along the flow, including r,s,a as extra state
    return sp.expand(sp.diff(f, U) * F[0] + sp.diff(f, p) * F[1] + sp.diff(f, y) * F[2]
                     + sp.diff(f, r) * (r ** 2 - mu * r - s) + sp.diff(f, s) * (s * r - s - x)
                     + sp.diff(f, a) * (-r * a))


def L3(f):
    return sp.simplify(D(D(D(f))) + A * D(D(f)) + mu * D(f) + x * f)


res = {}
pp, yp = F[1], F[2]
res['L3[1]=x'] = sp.simplify(L3(sp.Integer(1)) - x)
res['L3[U]=k g(U)'] = sp.simplify(L3(U) - k * (1 - sp.exp(U) + U * sp.exp(U)))
res['L3[p]=0'] = L3(p)
res["L3[p']=-x p^2"] = sp.simplify(L3(pp) + x * p ** 2)
res['L3[y]=x(k+p^2)'] = sp.simplify(L3(y) - x * (k + p ** 2))
res['L3[lambda-p]... Q[q]=x'] = sp.simplify(D(D(-p)) + A * D(-p) + mu * (k / mu - p) - x)
# strong-stable scalar: a''' + A a'' + mu a' + x a = 0  with a'=-r a
res['L3[a]=0'] = sp.simplify(L3(a))
# r'' + (A-3r) r' + r(r-1)(r-mu) = x
res["r''+(A-3r)r'+chi(r)=x"] = sp.simplify(D(D(r)) + (A - 3 * r) * D(r) + r * (r - 1) * (r - mu) - x)
# L2' = D^2 + (A-r) D + beta, beta = (r-1)(r-mu) - 2 r' ;  L3[phi] = L2'[phi' + r phi]
beta = (r - 1) * (r - mu) - 2 * D(r)
L2p = lambda psi: D(D(psi)) + (A - r) * D(psi) + beta * psi
for name, phi in (('1', sp.Integer(1)), ('U', U), ("p'", pp), ('y', y), ('p', p)):
    res[f"L3[{name}] = L2'[(D+r){name}]"] = sp.simplify(L3(phi) - L2p(D(phi) + r * phi))
res["L2'[r]=x"] = sp.simplify(L2p(r) - x)
# weak adjoint Psi^ = z' x e_s^/ a  (divide by a): components
psi_p = pp + r * p
P1 = r * yp - s * pp; P2 = yp + s * p; P3 = -psi_p
# adjoint equation for Psi = e^{At} a * (P1,P2,P3):   Psi' = -J^T Psi
Psi = sp.Matrix([P1, P2, P3])
lhs = sp.Matrix([D(q) for q in Psi]) + (A - r) * Psi          # derivative of e^{At} a Psi divided by e^{At} a
res["Psi' = -J^T Psi"] = sp.simplify(lhs + J.T * Psi)
res["Psi.z'=0"] = sp.simplify(Psi.dot(F))
# theta = -p'/p :  theta' = theta^2 - mu theta + y'/p
theta = -pp / p
res["theta'"] = sp.simplify(D(theta) - (theta ** 2 - mu * theta + yp / p))
# at a zero of psi_p (r = theta): (theta - r)' = y'/p + s
res["(theta-r)' at zero"] = sp.simplify((D(theta) - D(r)).subs(r, theta) - (yp / p + s))
# zeros of Psi2: y'/p = -s ; derivative of (y'/p + s) there equals s*psi_p/p
e2 = D(yp / p + s)
res["(y'/p+s)' at zero = s psi_p/p"] = sp.simplify(e2.subs(s, -yp / p) - ((-yp / p) * (pp + r * p) / p))
# zeros of Psi1: y'/p' = s/r ; derivative of (y'/p' - s/r) there equals x psi_p /(r p')
e1 = D(yp / pp - s / r)
res["(y'/p'-s/r)' at zero = x psi_p/(r p')"] = sp.simplify(e1.subs(s, r * yp / pp) - x * (pp + r * p) / (r * pp))
# L''(t1) = p'^2 - p p'' = p^2 r (r-mu) + p y'   when p' = -r p
ppp = D(pp)
res["L''(t1)"] = sp.simplify((pp ** 2 - p * ppp) - (pp ** 2 + mu * p * pp + p * yp))
# sigma at t1: Psi.z'' with Psi=c(p',-p,0): c(p'^2 - p p'')
# compound: w = v ^ u, w' = C w
v1, v2, v3, u1, u2, u3 = sp.symbols('v1 v2 v3 u1 u2 u3')
v = sp.Matrix([v1, v2, v3]); u = sp.Matrix([u1, u2, u3])
w = sp.Matrix([v2 * u1 - v1 * u2, v1 * u3 - v3 * u1, v2 * u3 - v3 * u2])
C = sp.Matrix([[-mu, 1, 0], [0, -1, 1], [x, 0, -A]])
dv = J * v; du = J * u
dw = sp.Matrix([sum(sp.diff(wi, q) * dq for q, dq in zip(list(v) + list(u), list(dv) + list(du))) for wi in w])
res["compound w'=C w"] = sp.simplify(dw - C * w)
cr = v.cross(u)
dcr = sp.Matrix([sum(sp.diff(ci, q) * dq for q, dq in zip(list(v) + list(u), list(dv) + list(du))) for ci in cr])
res["(v x u)' = (tr J - J^T)(v x u)"] = sp.simplify(dcr - (J.trace() * sp.eye(3) - J.T) * cr)
# Leibniz for L3 on products:  L3[f g] = g L3 f + f L3 g - x f g + 3 f'' g' + 3 f' g'' + 2A f' g'
# chord identity: F(z2)-F(z1)-J(z1)(z2-z1) = (0,0,x1 (e^{d}-1-d))
U2, p2, y2 = sp.symbols('U2 p2 y2')
F2 = F.subs({U: U2, p: p2, y: y2}, simultaneous=True)
dvec = sp.Matrix([U2 - U, p2 - p, y2 - y])
res['chord remainder'] = sp.simplify(F2 - F - J * dvec - sp.Matrix([0, 0, x * (sp.exp(U2 - U) - 1 - (U2 - U))]))
# scaled curve residual: E + th (z - E):  third component residual  k[th(e^U-1) - (e^{th U}-1)]
th = sp.symbols('theta', positive=True)
zs = sp.Matrix([th * U, th * p, k + th * (y - k)])
Fs = F.subs({U: zs[0], p: zs[1], y: zs[2]}, simultaneous=True)
res['scaled residual'] = sp.simplify(th * F - Fs - sp.Matrix([0, 0, k * (th * (sp.exp(U) - 1) - (sp.exp(th * U) - 1))]))

bad = 0
for name, val in res.items():
    z = val == 0 or (hasattr(val, 'is_zero_matrix') and val.is_zero_matrix) or sp.simplify(val) == 0 or (isinstance(val, sp.MatrixBase) and all(sp.simplify(e) == 0 for e in val))
    print(f"{'PASS' if z else 'FAIL'}  {name}")
    bad += (not z)
print(f"{len(res) - bad}/{len(res)} exact identities pass")
