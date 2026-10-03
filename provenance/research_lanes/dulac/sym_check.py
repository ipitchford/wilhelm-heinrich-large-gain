#!/usr/bin/env python3
"""Exact (sympy) verification of the algebraic identities used in Theorem D.
U,p,y,eta,zeta are independent symbols; D is the derivation defined by the ODEs
   U'=p, p'=k-mu p-y, y'=x-y, eta'=zeta+(mu-1)eta-eta^2, zeta'=x-(1+eta)zeta,  x=k e^U.
O = zeta p - eta p' - (x-y) is the tangent-orthogonality constraint (O=0 on cycles; D O = -(1+eta) O identically).
Every identity is reported as  LHS - RHS = (explicit multiple of O)."""
import sympy as sp
U, p, y, eta, zeta, mu, k = sp.symbols('U p y eta zeta mu k')
x = k * sp.exp(U)
A = 1 + mu
vf = {U: p, p: k - mu * p - y, y: x - y, eta: zeta + (mu - 1) * eta - eta**2, zeta: x - (1 + eta) * zeta}
D = lambda f: sp.expand(sum(sp.diff(f, v) * r for v, r in vf.items()))
pd = vf[p]; H = y - k; Hd = x - y
pdd = D(pd)
O = zeta * p - eta * pd - (x - y)
K = zeta * U - eta * p - y + k
W = D(zeta) * p - D(eta) * pd
V = p**2 - U * pd
Q = pd**2 - p * pdd
VUH = U * Hd - p * H
psi = U * sp.exp(U) - sp.exp(U) + 1
ok = True
def chk(name, expr, mult=None):
    """expr should equal mult*O (mult=None means identically zero)."""
    global ok
    r = sp.simplify(expr - (mult * O if mult is not None else 0))
    print(f'{name:<62} residual = {r}' + ('' if mult is None else f'   [= ({mult})*O off-constraint]'))
    ok &= (r == 0)
chk("D O + (1+eta) O = 0", D(O) + (1 + eta) * O)
chk("K' + (1+eta)K = k psi(U)", D(K) + (1 + eta) * K - k * psi)
chk("W' + (1+eta)W = x p^2  (mod O)", D(W) + (1 + eta) * W - x * p**2, sp.simplify((D(W) + (1 + eta) * W - x * p**2) / O))
chk("V' = (eta-mu)V + pK   (mod O)", D(V) - ((eta - mu) * V + p * K), -U)
chk("V' = -mu V + V_UH      (state identity)", D(V) - (-mu * V + VUH))
chk("V_UH = pK + eta V      (mod O)", VUH - (p * K + eta * V), sp.simplify((VUH - (p * K + eta * V)) / O))
chk("Q = p'(k-y) + p(x-y)   (state identity)", Q - (pd * (k - y) + p * (x - y)))
chk("Q' = x V - A Q + k psi p'  (state identity)", D(Q) - (x * V - A * Q + k * psi * pd))
resQ = sp.simplify((D(Q) - ((eta - mu) * Q + p * W)) / O)
chk("Q' = (eta-mu)Q + pW    (mod O)", D(Q) - ((eta - mu) * Q + p * W), resQ)
chk("V_UH' = -V_UH + Q + k psi p (state identity)", D(VUH) - (-VUH + Q + k * psi * p))
chk("W = x p - (1+eta-mu)(x-y) + zeta (y-k)  (mod O)", W - (x * p - (1 + eta - mu) * Hd + zeta * H), sp.simplify((W - (x * p - (1 + eta - mu) * Hd + zeta * H)) / O))
# K-parametrisation: eta = eta0 - p K/V, zeta = zeta0 - p' K/V on O=0
chk("eta V - (V_UH - pK) = 0   (mod O)", eta * V - (VUH - p * K), sp.simplify((eta * V - (VUH - p * K)) / O))
chk("zeta V - (Q - p'K) = 0    (mod O)", zeta * V - (Q - pd * K), sp.simplify((zeta * V - (Q - pd * K)) / O))
# L2 building blocks: L2[f] = f'' + (mu-eta) f' + zeta f
L2 = lambda f: D(D(f)) + (mu - eta) * D(f) + zeta * f
chk("L2[1] = zeta", L2(sp.Integer(1)) - zeta)
chk("L2[U] = K", L2(U) - K)
chk("L2[p] = 0             (mod O)", L2(p), sp.simplify(L2(p) / O))
chk("L2[y-k] = W           (mod O)", L2(H) - W, sp.simplify((L2(H) - W) / O))
# Riccati defects of the explicit pair (eta0,zeta0) = (V_UH,Q)/V
eta0 = VUH / V; zeta0 = Q / V
chk("eta0' - [zeta0+(mu-1)eta0-eta0^2] = k psi p / V", sp.simplify(D(eta0) - (zeta0 + (mu - 1) * eta0 - eta0**2) - k * psi * p / V))
chk("zeta0' - [x-(1+eta0)zeta0] = k psi p' / V", sp.simplify(D(zeta0) - (x - (1 + eta0) * zeta0) - k * psi * pd / V))
print('ALL IDENTITIES VERIFIED' if ok else 'FAILURE')
