#!/usr/bin/env python3
"""Referee check 1: every algebraic identity used in the lane's proofs, verified symbolically (sympy).
Independent of the lane's scripts.  Along the flow:
   U' = p, p' = k - mu p - y, y' = k e^U - y,
   eta' = zeta + (mu-1) eta - eta^2, zeta' = x - (1+eta) zeta.
The tangent-orthogonality constraint O := zeta p - eta p' - (x-y) = 0 is imposed by eliminating zeta
(valid where p != 0; identities are polynomial after clearing p, so they extend by continuity)."""
import sympy as sp

U, p, y, eta, zeta, mu, k, c = sp.symbols('U p y eta zeta mu k c', real=True)
x = k*sp.exp(U)
A = 1 + mu
H = y - k
pd = k - mu*p - y
yd = x - y
etad = zeta + (mu-1)*eta - eta**2
zetad = x - (1+eta)*zeta
state = [U, p, y, eta, zeta]
rhs = [p, pd, yd, etad, zetad]

def D(F):
    return sum(sp.diff(F, s)*r for s, r in zip(state, rhs))

O = zeta*p - eta*pd - yd                      # orthogonality defect (=0 on a cycle)
zeta_on = sp.solve(sp.Eq(O, 0), zeta)[0]      # zeta on the constraint

def on(F):
    return sp.simplify(F.subs(zeta, zeta_on))

res = {}
pdd = D(pd)
S1 = p**2 - U*pd
S2 = U*yd - p*H
S3 = p*yd - pd*H
S4 = pd**2 - p*pdd
K = zeta*U - eta*p - H
W = zetad*p - etad*pd
G1 = U*sp.exp(U) - sp.exp(U) + 1

# --- invariance of the constraint: D(O) = -(1+eta) O (so imposing O=0 is consistent)
res['D(O) + (1+eta) O == 0 (identically)'] = sp.simplify(D(O) + (1+eta)*O)
# --- paper identities
res['K\' + (1+eta)K - k G1 (identically)'] = sp.simplify(D(K) + (1+eta)*K - k*G1)
res['W\' + (1+eta)W - x p^2  on O=0'] = on(D(W) + (1+eta)*W - x*p**2)
res['tangent eq p\'\' - (eta-mu)p\' + zeta p on O=0'] = on(pdd - (eta-mu)*pd + zeta*p)
# --- Theorem R
res['S1\' - (eta-mu)S1 - pK  on O=0'] = on(D(S1) - (eta-mu)*S1 - p*K)
res['S4\' - (eta-mu)S4 - pW  on O=0'] = on(D(S4) - (eta-mu)*S4 - p*W)
sc = p**2 - pd*(U + c)
res['s_c\' - (eta-mu)s_c - p(K + c zeta) on O=0'] = on(D(sc) - (eta-mu)*sc - p*(K + c*zeta))
res['S1 at p=0 equals U(y-k)'] = sp.simplify(S1.subs(p, 0) - U*(y-k))
res['S4 at p=0 equals (k-y)^2'] = sp.simplify(S4.subs(p, 0) - (k-y)**2)
res['s_c at p=0 equals (y-k)(U+c)'] = sp.simplify(sc.subs(p, 0) - (y-k)*(U+c))
# --- corollaries
res['S2 - eta S1 - pK on O=0'] = on(S2 - eta*S1 - p*K)
res['S3 - zeta S1 - p\'K on O=0'] = on(S3 - zeta*S1 - pd*K)
res['S1\' - (S2 - mu S1) (identically, no eta)'] = sp.simplify(D(S1) - (S2 - mu*S1))
res['(zeta - y)\' + (zeta - y) + eta zeta (identically)'] = sp.simplify(D(zeta - y) + (zeta - y) + eta*zeta)
# --- E quadratic form
gE = sp.symbols('gamma_E', positive=True)
kE = gE*(gE-1)*(gE-mu)
etaE = gE - 1; zetaE = etaE*(etaE + 1 - mu); al = (etaE - mu)/2
QE = p**2 - 2*al*U*p + zetaE*U**2
LamE = H - zetaE*U + etaE*p
res['S1 - QE - U LamE (identically)'] = sp.simplify((S1 - QE - U*LamE))
DQ = (D(QE) - 2*al*QE + 2*(p - al*U)*LamE).subs(k, kE)
res['D QE - 2 alpha QE + 2(p - alpha U) LamE  (k = chi(gamma_E))'] = sp.simplify(DQ)
DL = (D(LamE) + gE*LamE - k*(sp.exp(U) - 1 - U)).subs(k, kE)
res['LamE\' + gamma_E LamE - k(e^U-1-U)  (k = chi(gamma_E))'] = sp.simplify(DL)
# --- scalar Riccati for K, N on O=0 with eta = (S2 - pK)/S1
Ks = sp.symbols('Ks')
eta_from_K = (S2 - p*Ks)/S1
RK = k*G1 - (1 + S2/S1)*Ks + (p/S1)*Ks**2
res['(R_K): [kG1 - (1+eta)K] - RHS with eta=(S2-pK)/S1'] = sp.simplify(k*G1 - (1 + eta_from_K)*Ks - RK)
N = sp.symbols('N')
# N = K/S1, N' = K'/S1 - N S1'/S1 with S1' = S2 - mu S1
Np = (RK.subs(Ks, N*S1))/S1 - N*(S2 - mu*S1)/S1
res['N\' - [kG1/S1 - (1-mu+2 S2/S1) N + p N^2]'] = sp.simplify(Np - (k*G1/S1 - (1 - mu + 2*S2/S1)*N + p*N**2))
# --- Lemma A0: defects of omega = phi (1, e, z) in the compound system
e, z, phi, de, dz, dphi = sp.symbols('e z phi de dz dphi')
C = sp.Matrix([[-mu, 1, 0], [0, -1, 1], [x, 0, -A]])
om = sp.Matrix([phi, phi*e, phi*z])
omd = sp.Matrix([dphi, dphi*e + phi*de, dphi*z + phi*dz])
g = [dphi/phi + mu - e, dphi/phi + de/e + 1 - z/e, dphi/phi + dz/z + A - x/z]
dfc = omd - C*om
res['A0 defect rows'] = [sp.simplify(dfc[i] - g[i]*om[i]) for i in range(3)]
# at the true pair the three g_i coincide with mu - eta (phi=1)
gt = [mu - eta, etad/eta + 1 - zeta/eta, zetad/zeta + A - x/zeta]
res['A0 at true pair: g2-g1, g3-g1'] = [sp.simplify(gt[1]-gt[0]), sp.simplify(gt[2]-gt[0])]
# compound matrix: C = S J^[2] S
J = sp.Matrix([[0, 1, 0], [0, -mu, -1], [x, 0, -1]])
J2 = sp.Matrix([[J[0,0]+J[1,1], J[1,2], -J[0,2]], [J[2,1], J[0,0]+J[2,2], J[0,1]], [-J[2,0], J[1,0], J[1,1]+J[2,2]]])
Sg = sp.diag(1, -1, -1)
res['C - S J^[2] S'] = sp.simplify(Sg*J2*Sg - C)
# --- explicit almost-solutions
P = k/mu - p
omY = sp.Matrix([P, y, x]); res['omega^Y defect - (0,0,x(2p+A-lambda))'] = [sp.simplify(v) for v in (sp.Matrix([D(q) for q in omY]) - C*omY - sp.Matrix([0, 0, x*(2*p + A - k/mu)]))]
omc = sp.Matrix([1, mu, mu]); res['omega^c defect - (0,0,kH - x)'] = [sp.simplify(v) for v in (-C*omc - sp.Matrix([0, 0, mu*A - x]))]
omR = sp.Matrix([S1, S2, S3]); res['omega^R defect - kG1 (0,p,p\')'] = [sp.simplify(v) for v in (sp.Matrix([D(q) for q in omR]) - C*omR - k*G1*sp.Matrix([0, p, pd]))]
# --- Prop B: adjoint nu = e^{At}(c,b,a) ; a'=-b, b'=-mu b - c, c' = -c - x a
a_, b_, c_ = sp.symbols('a b c_')
nu = sp.Matrix([c_, b_, a_]); nud = sp.Matrix([A*c_ + (-c_ - x*a_), A*b_ + (-mu*b_ - c_), A*a_ + (-b_)])
res['nu\' + C^T nu'] = [sp.simplify(v) for v in (nud + C.T*nu)]
be, ka = sp.symbols('beta kappa')
bed = be**2 - mu*be - ka; kad = ka*(be - 1) - x
res['(kappa + beta^2/2 - beta)\' - (chi(beta) - x)'] = sp.simplify(kad + be*bed - bed - (be*(be-1)*(be-mu) - x))
# --- Prop T: polar equation and linear-fractional difference
th, r, m, m1, m2 = sp.symbols('theta r m m1 m2')
Up, pp = r*sp.cos(th), r*sp.sin(th)
ppd = -mu*pp - m*r
rdot = (Up*pp + pp*ppd)/r
thdot = (Up*ppd - pp*pp)/r**2
sig = sp.sin(th)**2 + mu*sp.sin(th)*sp.cos(th) + m*sp.cos(th)
Xi = sp.sin(th)*(m + mu*sp.sin(th) - sp.cos(th))/sig
res['theta\' + sigma'] = sp.simplify(thdot + sig)
res['(r\'/r)/theta\' - Xi'] = sp.simplify(rdot/r/thdot - Xi)
res['Xi(m2)-Xi(m1) - sin (m2-m1)/(sig1 sig2)'] = sp.simplify(Xi.subs(m, m2) - Xi.subs(m, m1) - sp.sin(th)*(m2-m1)/(sig.subs(m, m1)*sig.subs(m, m2)))
# --- Lemma A1: first-order coefficients solve the linearised Riccati system on the centre plane (U''=-mu U)
Dn = mu**2 + 6*mu + 1
al1, al2, be1, be2 = mu*A/Dn, -2*mu/Dn, mu*(mu**2 + 4*mu + 1)/Dn, -mu*A/Dn
U1, p1 = sp.symbols('U1 p1')
a1 = al1*U1 + al2*p1; b1 = be1*U1 + be2*p1
def D1(F):   # U1' = p1, p1' = -mu U1
    return sp.diff(F, U1)*p1 + sp.diff(F, p1)*(-mu*U1)
res['A1: a\' - (b - A a)'] = sp.simplify(D1(a1) - (b1 - A*a1))
res['A1: b\' - (mu A u - A b - mu a)'] = sp.simplify(D1(b1) - (mu*A*U1 - A*b1 - mu*a1))
res['A1: eigenvalues of [[-A,1],[-mu,-A]]'] = list(sp.Matrix([[-A, 1], [-mu, -A]]).eigenvals().keys())

bad = 0
for name, val in res.items():
    ok = (val == 0) or (isinstance(val, list) and all(v == 0 for v in val)) or name.startswith('A1: eigen')
    if hasattr(val, 'is_zero_matrix'):
        ok = val.is_zero_matrix
    bad += (not ok)
    print(('OK   ' if ok else 'FAIL ') + name + ('' if ok and not name.startswith('A1: eigen') else '   -> ' + str(val)))
print('\nfailures:', bad)
