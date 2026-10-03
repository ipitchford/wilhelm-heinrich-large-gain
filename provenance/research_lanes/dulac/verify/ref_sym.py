#!/usr/bin/env python3
"""Referee's own exact check of the algebra behind Theorem D / O1 / O2 (written from the statements in NOTES.md,
not from the lane's sym_check.py).  x is an independent symbol with D x = x p and U enters polynomially, so every
residual is a polynomial identity; exp never appears.  psi-term: k*psi(U) = x U - x + k.

Each line prints  residual  and, where the identity only holds on the tangent-orthogonality constraint O = 0,
the exact cofactor c with  residual = c * O  (so the identity holds on cycles, where O = 0)."""
import sympy as sp

U, p, y, x, eta, zeta, mu, k = sp.symbols('U p y x eta zeta mu k')
A = 1 + mu
vf = {U: p, p: k - mu * p - y, y: x - y, x: x * p, eta: zeta + (mu - 1) * eta - eta**2, zeta: x - (1 + eta) * zeta}


def D(f):
    return sp.expand(sum(sp.diff(f, v) * r for v, r in vf.items()))


pd = vf[p]; pdd = D(pd); H = y - k; Hd = x - y
O = zeta * p - eta * pd - Hd
K = zeta * U - eta * p - H
W = D(zeta) * p - D(eta) * pd
V = p**2 - U * pd
Q = pd**2 - p * pdd
V3 = U * Hd - p * H
kpsi = x * U - x + k
allok = True


def exact(name, expr):
    global allok
    r = sp.expand(expr)
    print(f'{name:<70} residual = {r}')
    allok &= (r == 0)


def modO(name, expr):
    """expr must be a polynomial multiple of O"""
    global allok
    num = sp.expand(expr)
    q, r = sp.div(sp.Poly(num, U, p, y, x, eta, zeta), sp.Poly(sp.expand(O), U, p, y, x, eta, zeta))
    ok = (r.as_expr() == 0)
    print(f'{name:<70} = ({sp.factor(q.as_expr())}) * O   remainder = {r.as_expr()}')
    allok &= ok


print('--- v2.0 inputs re-derived')
exact("D O + (1+eta) O", D(O) + (1 + eta) * O)
exact("K' + (1+eta) K - k psi", D(K) + (1 + eta) * K - kpsi)
modO("W' + (1+eta) W - x p^2", D(W) + (1 + eta) * W - x * p**2)
modO("p'' - [(eta-mu)p' - zeta p]", pdd - ((eta - mu) * pd - zeta * p))
print('--- Theorem D(b)')
modO("V' - [(eta-mu)V + pK]", D(V) - ((eta - mu) * V + p * K))
exact("V' - [-mu V + V3]", D(V) - (-mu * V + V3))
modO("V3 - [pK + eta V]", V3 - (p * K + eta * V))
exact("Q - [p'(k-y) + p(x-y)]", Q - (pd * (k - y) + p * (x - y)))
modO("Q' - [(eta-mu)Q + pW]", D(Q) - ((eta - mu) * Q + p * W))
exact("Q' - [xV - AQ + k psi p']", D(Q) - (x * V - A * Q + kpsi * pd))
print('--- Theorem D(c): linear system for (eta,zeta), kappa equation')
modO("eta V - (V3 - pK)", eta * V - (V3 - p * K))
modO("zeta V - (Q - p'K)", zeta * V - (Q - pd * K))
kap = K / V; eta0 = V3 / V; zeta0 = Q / V
num = sp.together(D(kap) - (kpsi / V - (1 + eta + eta0 - mu) * kap))
modO("V^2 * {kappa' - [k psi/V - (1+eta+eta0-mu) kappa]}", sp.numer(num) * (V**2 / sp.denom(num)))
eta1 = (x * V + kpsi * pd) / Q - 1
exact("(log V)' - (eta0 - mu)   [times V]", sp.simplify(D(V) - (eta0 - mu) * V))
exact("(log Q)' - (eta1 - mu)   [times Q]", sp.expand(sp.simplify(D(Q) - (eta1 - mu) * Q)))
print('--- Theorem D(e): compound system  omega\' = C omega + k psi (0,p,p\')')
exact("V3' - [-V3 + Q + k psi p]", D(V3) - (-V3 + Q + kpsi * p))
print('   wedge r ^ z\' with r=(U,p,H), z\'=(p,p\',x-y):  -(12), (13), (23) components vs (V, V3, Q):',
      sp.expand(-(U * pd - p * p) - V), sp.expand((U * Hd - H * p) - V3), sp.expand((p * Hd - H * pd) - Q))
r1 = sp.simplify(D(eta0) - (zeta0 + (mu - 1) * eta0 - eta0**2) - kpsi * p / V)
r2 = sp.simplify(D(zeta0) - (x - (1 + eta0) * zeta0) - kpsi * pd / V)
print(f'{"Riccati defects of (eta0,zeta0) minus k psi (p,p\')/V":<70} residual = {r1}, {r2}'); allok &= (r1 == 0 and r2 == 0)
print('--- general Wronskian principle')
f, g = sp.Function('f'), sp.Function('g')
L2 = lambda h: D(D(h)) + (mu - eta) * D(h) + zeta * h
exact("L2[1] - zeta", L2(sp.Integer(1)) - zeta)
exact("L2[U] - K", L2(U) - K)
modO("L2[p]", L2(p))
modO("L2[y-k] - W", L2(H) - W)
print('--- O1: W = W0(z) + kappa Q with W0 = x p - (1+eta0-mu)(x-y) + zeta0 (y-k);  on O=0')
modO("W - [x p - (1+eta-mu)(x-y) + zeta (y-k)]", W - (x * p - (1 + eta - mu) * Hd + zeta * H))
# substitute eta = eta0 - p kappa, zeta = zeta0 - p' kappa (kappa free symbol) into the state form of W:
kp = sp.Symbol('kappa')
Wst = (x * p - (1 + eta - mu) * Hd + zeta * H).subs({eta: eta0 - p * kp, zeta: zeta0 - pd * kp})
W0 = x * p - (1 + eta0 - mu) * Hd + zeta0 * H
exact("W(eta0 - p kappa, zeta0 - p' kappa) - [W0 + kappa Q]", sp.simplify(Wst - (W0 + kp * Q)))
print('--- O2: the kappa=0 plane through z contains z\' and E')
exact("V*(zeta0 p - eta0 p' - (x-y))", sp.simplify(V * (zeta0 * p - eta0 * pd - Hd)))
exact("V*(zeta0 U - eta0 p - (y-k))", sp.simplify(V * (zeta0 * U - eta0 * p - H)))
print('--- O2: crude dynamic upper bounds')
exact("(y-zeta)' + (y-zeta) - eta zeta", D(y - zeta) + (y - zeta) - eta * zeta)
lam = k / mu
exact("(lam-p-eta)' + mu(lam-p-eta) - [(y-zeta) + eta(eta+1-2mu)]", D(lam - p - eta) + mu * (lam - p - eta) - ((y - zeta) + eta * (eta + 1 - 2 * mu)))
print('--- B1(ii): at E, R(E) = eta_E - mu = gamma_E - A = 2 sigma')
L, g_, s_, w2 = sp.symbols('lambda gamma sigma w2')
c = sp.Poly(sp.expand((L + g_) * (L**2 - 2 * s_ * L + w2)) - (L**3 + A * L**2 + mu * L + k), L).all_coeffs()
print('   coefficient equations of (lam+gamma)(lam^2-2 sigma lam+|lam_u|^2) = lam^3+A lam^2+mu lam+k :', c)
print('   => gamma - 2 sigma = A (lambda^2 coefficient);  eta_E = gamma-1 fixed point of Riccati with x=k:',
      sp.expand((x - (1 + eta) * zeta).subs({x: g_ * (g_ - 1) * (g_ - mu), zeta: (g_ - 1) * (g_ - mu), eta: g_ - 1})),
      sp.expand((zeta + (mu - 1) * eta - eta**2).subs({zeta: (g_ - 1) * (g_ - mu), eta: g_ - 1})))
print('--- brief 3(d)/NOTES sec.3: Q0 = p^2 - 2 sigma U p + zeta_E U^2 has D Q0 = 2 sigma Q0 in the planar linearisation u\'=v, v\'=-zeta_E u+2 sigma v')
u_, v_, zE = sp.symbols('u v zeta_E')
Q0 = v_**2 - 2 * s_ * u_ * v_ + zE * u_**2
exact("D Q0 - 2 sigma Q0", sp.diff(Q0, u_) * v_ + sp.diff(Q0, v_) * (-zE * u_ + 2 * s_ * v_) - 2 * s_ * Q0)
print('ALL EXACT CHECKS PASSED' if allok else 'SOME CHECK FAILED')
