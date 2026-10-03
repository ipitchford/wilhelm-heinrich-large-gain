#!/usr/bin/env python3
"""Checks for Section 5 (s05_gauges.tex) of the v3.0 paper, in the paper's global notation.

Part A: exact (sympy) verification of every identity stated in the section after the renaming
        V,S*,S1 -> Sigma;  Q,V,S4 -> Xi;  V3,S2 -> Sigma*eta0;  beta,kappa(lane) -> frak r, frak s;  q,P -> ell.
Part B: floating-point spot checks on cycles located with the shared library ROOT/work/whlib.py.

Location: ROOT/work/framework_checks/ (moved out of the paper's sections directory at the audit of
3 October 2026).  Block A7 and block (j) of Part B were added at that audit: remark formulas, intermediate
proof steps, the shifted centre, the dynamic upper bounds, supporting lines and planes on ordered pairs,
the branch orderings, gamma < gamma_E, and the compound monodromy.

Run:  /Users/admin/.venvs/analysis-py313/bin/python -W ignore s05_gauges.check.py > s05_gauges.check.out
Nothing here is an interval certificate.
"""
import sys
sys.path.insert(0, '/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03/work')
import numpy as np
import sympy as sp

FAIL = []


def ok(name, expr):
    if isinstance(expr, sp.MatrixBase):
        e = expr.applyfunc(sp.simplify)
        good = all(v == 0 for v in e)
    else:
        e = sp.simplify(expr)
        good = (e == 0)
        if not good:
            e = sp.simplify(sp.expand(sp.together(expr)))
            good = (e == 0)
    print(f"  [{'OK' if good else 'FAIL'}] {name}")
    if not good:
        FAIL.append(name)
        print("        residual:", e)


# ---------------------------------------------------------------------------------------------
print("PART A: exact identities (sympy)")
t = sp.symbols('t')
mu, k, cc = sp.symbols('mu k c', positive=True)
A = 1 + mu
lam = k / mu
U, p, y, eta, zeta = [sp.Function(n)(t) for n in ('U', 'p', 'y', 'eta', 'zeta')]
x = k * sp.exp(U)
pd = k - mu * p - y          # p'
yd = x - y                   # y'
flow = {sp.Derivative(U, t): p, sp.Derivative(p, t): pd, sp.Derivative(y, t): yd,
        sp.Derivative(eta, t): zeta + (mu - 1) * eta - eta ** 2,
        sp.Derivative(zeta, t): x - (1 + eta) * zeta}


def D(e, n=1):
    for _ in range(n):
        e = sp.diff(e, t).subs(flow)
    return e


zeta_on_O = (x - y + eta * pd) / p          # tangent orthogonality  zeta p - eta p' = x - y


def onO(e):
    return e.subs(zeta, zeta_on_O)


psi = U * sp.exp(U) - sp.exp(U) + 1
K = zeta * U - eta * p - y + k
W = D(zeta) * p - D(eta) * pd
Sg = p ** 2 - U * pd
Xg = pd ** 2 - p * D(pd)
N0 = U * (x - y) - p * (y - k)               # = Sigma * eta0
eta0 = N0 / Sg
zeta0 = Xg / Sg
kap = K / Sg
Oc = zeta * p - eta * pd - (x - y)

print(" A1. tangent equation and Wronskian gauge lemma")
ok("p'' = (eta-mu)p' - zeta p   [on O=0]", onO(D(pd) - ((eta - mu) * pd - zeta * p)))
L2 = lambda f: D(f, 2) + (mu - eta) * D(f) + zeta * f
for nm, phi in (("U", U), ("-p'", -pd), ("U+c", U + cc)):
    S = p * D(phi) - pd * phi
    ok(f"S_phi' = (eta-mu)S_phi + p L2[phi], phi={nm}   [on O=0]", onO(D(S) - ((eta - mu) * S + p * L2(phi))))
ok("S_U = Sigma", p * D(U) - pd * U - Sg)
ok("S_{-p'} = Xi", p * D(-pd) + pd * pd - Xg)
ok("L2[U] = K", L2(U) - K)
ok("L2[-p'] = W   [on O=0]", onO(L2(-pd) - W))
ok("L2[1] = zeta", L2(sp.Integer(1)) - zeta)
ok("L2[p] = 0   [on O=0]", onO(L2(p)))

print(" A2. gauge equations (eq:gaugeodes), state forms, eta0, eta1")
ok("Sigma' = (eta-mu)Sigma + pK   [on O=0]", onO(D(Sg) - ((eta - mu) * Sg + p * K)))
ok("Xi' = (eta-mu)Xi + pW   [on O=0]", onO(D(Xg) - ((eta - mu) * Xg + p * W)))
ok("Sigma' = -mu Sigma + [U(x-y) - p(y-k)]   (identity)", D(Sg) - (-mu * Sg + N0))
ok("Sigma eta0 = eta Sigma + pK   [on O=0]", onO(N0 - (eta * Sg + p * K)))
ok("Xi = p'(k-y) + p(x-y)   (identity)", Xg - (pd * (k - y) + p * (x - y)))
ok("Xi' = x Sigma - A Xi + k psi p'   (identity)", D(Xg) - (x * Sg - A * Xg + k * psi * pd))
eta1 = (x * Sg + k * psi * pd) / Xg - 1
ok("(log Xi)' = eta1 - mu   (identity)", D(Xg) / Xg - (eta1 - mu))
ok("eta - eta1 = -pW/Xi   [on O=0]", onO(eta - eta1 + p * W / Xg))
Sc = p ** 2 - (U + cc) * pd
ok("Sigma_c' = (eta-mu)Sigma_c + p(K + c zeta)   [on O=0]", onO(D(Sc) - ((eta - mu) * Sc + p * (K + cc * zeta))))
ok("K' + (1+eta)K = k psi(U)   (identity)", D(K) + (1 + eta) * K - k * psi)
ok("W' + (1+eta)W = x p^2   [on O=0]", onO(D(W) + (1 + eta) * W - x * p ** 2))

print(" A3. scalar reduction (prop:kappa)")
ok("eta = eta0 - p kappa   [on O=0]", onO(eta - (eta0 - p * kap)))
ok("zeta = zeta0 - p' kappa   [on O=0]", onO(zeta - (zeta0 - pd * kap)))
ok("kappa' = k psi/Sigma - (1+eta+eta0-mu) kappa   (identity)", D(kap) - (k * psi / Sg - (1 + eta + eta0 - mu) * kap))
ok("kappa' = k psi/Sigma - (1+2eta0-mu)kappa + p kappa^2   [on O=0]",
   onO(D(kap) - (k * psi / Sg - (1 + 2 * eta0 - mu) * kap + p * kap ** 2)))
kk = sp.Symbol('kappa_')
e_, z_ = sp.symbols('e_ z_')
Wfun = lambda e, z: p * (x - (1 + e) * z) - pd * (z + (mu - 1) * e - e ** 2)
ok("W as a function of (eta,zeta)", W - Wfun(eta, zeta))
ok("W(eta0-p kappa, zeta0-p' kappa) = W(eta0,zeta0) + kappa Xi   (identity in kappa)",
   Wfun(eta0 - p * kk, zeta0 - pd * kk) - (Wfun(eta0, zeta0) + kk * Xg))
ok("K(eta0-p kappa, zeta0-p' kappa) = kappa Sigma   (identity in kappa)",
   (zeta0 - pd * kk) * U - (eta0 - p * kk) * p - y + k - kk * Sg)
ok("kappa=0 plane contains the tangent: zeta0 p - eta0 p' = x-y", zeta0 * p - eta0 * pd - (x - y))
Cm = sp.Matrix([[-mu, 1, 0], [0, -1, 1], [x, 0, -A]])
om = sp.Matrix([Sg, N0, Xg])
ok("(Sigma, Sigma eta0, Xi)' = C(.) + k psi (0,p,p')   (identity)",
   sp.Matrix([D(om[i]) for i in range(3)]) - Cm * om - k * psi * sp.Matrix([0, p, pd]))
ok("eta0' - [zeta0 + (mu-1)eta0 - eta0^2] = k psi p/Sigma", D(eta0) - (zeta0 + (mu - 1) * eta0 - eta0 ** 2) - k * psi * p / Sg)
ok("zeta0' - [x - (1+eta0)zeta0] = k psi p'/Sigma", D(zeta0) - (x - (1 + eta0) * zeta0) - k * psi * pd / Sg)
# two-point identity (O2): base point 1 with slopes (zeta0 - p' kappa, -(eta0 - p kappa)), second point (U2,p2,y2)
U1, p1, y1, U2, p2, y2, k1 = sp.symbols('U1 p1 y1 U2 p2 y2 kappa1')
x1 = k * sp.exp(U1)
pd1 = k - mu * p1 - y1
Sg1 = p1 ** 2 - U1 * pd1
Xg1 = pd1 * (k - y1) + p1 * (x1 - y1)
e01 = (U1 * (x1 - y1) - p1 * (y1 - k)) / Sg1
z01 = Xg1 / Sg1
Delta = (z01 - pd1 * k1) * (U2 - U1) - (e01 - p1 * k1) * (p2 - p1) - (y2 - y1)
ok("Delta = [zeta0 dU - eta0 dp - dy] + kappa1 (p1 dp - p1' dU)",
   Delta - ((z01 * (U2 - U1) - e01 * (p2 - p1) - (y2 - y1)) + k1 * (p1 * (p2 - p1) - pd1 * (U2 - U1))))
ok("second point E: bracket vanishes, Delta = -kappa1 Sigma1",
   Delta.subs({U2: 0, p2: 0, y2: k}) + k1 * Sg1)

print(" A4. polar form (prop:polar, prop:nestedradii)")
r, th, s = sp.symbols('r theta s', real=True)
polar = {U: r * sp.cos(th), p: r * sp.sin(th), y: k + r * s}
pd_pol = pd.subs(polar)
ok("p' = -r(mu sin + s)", pd_pol + r * (mu * sp.sin(th) + s))
rdot = ((U * p + p * pd) / r).subs(polar)
thdot = ((U * pd - p ** 2) / r ** 2).subs(polar)
sig = sp.sin(th) ** 2 + mu * sp.sin(th) * sp.cos(th) + s * sp.cos(th)
ok("r'/r = sin(cos - mu sin - s)", rdot / r - sp.sin(th) * (sp.cos(th) - mu * sp.sin(th) - s))
ok("theta' = -(sin^2 + mu sin cos + s cos)", thdot + sig)
ok("Sigma = r^2 sigma", Sg.subs(polar) - r ** 2 * sig)
ydot_pol = (x - y).subs(polar)
sdot = ydot_pol / r - s * rdot / r
ok("s' = [k(e^{r cos}-1) - r s]/r - s sin(cos - mu sin - s)",
   sdot - ((k * (sp.exp(r * sp.cos(th)) - 1) - r * s) / r - s * sp.sin(th) * (sp.cos(th) - mu * sp.sin(th) - s)))
s1, s2 = sp.symbols('s1 s2', real=True)
g = lambda ss: sp.sin(th) * (ss + mu * sp.sin(th) - sp.cos(th)) / (sp.sin(th) ** 2 + mu * sp.sin(th) * sp.cos(th) + ss * sp.cos(th))
ok("dlog r/dtheta = g(s,theta)", (rdot / r) / thdot - g(s))
ok("g(s2)-g(s1) = sin (s2-s1)/(sigma1 sigma2)",
   g(s2) - g(s1) - sp.sin(th) * (s2 - s1) / (sig.subs(s, s1) * sig.subs(s, s2)))
ok("dg/ds = sin/sigma^2", sp.diff(g(s), s) - sp.sin(th) / sig ** 2)

print(" A5. Lotka-Volterra chain and second-order forms (prop:lvchain)")
ell = lam - p
r1, r2 = x / y, y / ell
ok("r1' = r1(lambda+1-r1-ell)", D(r1) - r1 * (lam + 1 - r1 - ell))
ok("r2' = r2(mu-1+r1-r2)", D(r2) - r2 * (mu - 1 + r1 - r2))
ok("ell' = ell(r2-mu)", D(ell) - ell * (r2 - mu))
fmu = lambda v: v * (v + 1) * (v + 1 - mu)
Lop = lambda w: D(w, 2) + (3 * w + 2 - mu) * D(w) + fmu(w)
ok("L[eta] = x", Lop(eta) - x)
ok("L[y/ell] = x(A+p)/ell", Lop(y / ell) - x * (A + p) / ell)
et, zt = y / ell, x / ell
ok("(y/ell)' = x/ell + (mu-1)y/ell - (y/ell)^2  (first Riccati equation, exact)", D(et) - (zt + (mu - 1) * et - et ** 2))
ok("(x/ell)' - [x - (1+y/ell)x/ell] = x(A+2p-lambda)/ell", D(zt) - (x - (1 + et) * zt) - x * (A + 2 * p - lam) / ell)
phi = sp.Function('phi')(t)
c0 = sp.Symbol('c0')
w = sp.diff(phi, t) / phi + c0 - 1
Lplain = sp.diff(w, t, 2) + (3 * w + 2 - mu) * sp.diff(w, t) + fmu(w)
Dc = lambda f, sh: sp.diff(f, t) + sh * f
ok("L[phi'/phi + c - 1] = (D+c)(D+c-1)(D+c-mu)phi/phi", Lplain - Dc(Dc(Dc(phi, c0 - mu), c0 - 1), c0) / phi)
omg = zeta + eta + eta ** 2 / 2
Gmu = lambda v: sp.Rational(3, 2) * v ** 2 + (2 - mu) * v
ok("eta' = omega - G_mu(eta)", D(eta) - (omg - Gmu(eta)))
ok("omega' = x - f_mu(eta)", D(omg) - (x - fmu(eta)))

print(" A6. strong-stable identities and the comparison lemma (prop:strongidentity, lem A0)")
a, b, c = [sp.Function(n)(t) for n in ('a', 'b', 'c')]
flow.update({sp.Derivative(a, t): -b, sp.Derivative(b, t): -mu * b - c, sp.Derivative(c, t): -c - x * a})
J = sp.Matrix([[0, 1, 0], [0, -mu, -1], [x, 0, -1]])
es = sp.Matrix([a, -b, -c])
ok("(a,-b,-c) solves delta' = J delta", sp.Matrix([D(es[i]) for i in range(3)]) - J * es)
ok("a''' + A a'' + mu a' + x a = 0", D(a, 3) + A * D(a, 2) + mu * D(a) + x * a)
nu = sp.exp(A * t) * sp.Matrix([c, b, a])
ok("nu = e^{At}(c,b,a) solves nu' = -C^T nu", sp.Matrix([D(nu[i]) for i in range(3)]) + Cm.T * nu)
oY = sp.Matrix([ell, y, x])
ok("omega^Y=(ell,y,x): omega' - C omega = (0,0,x(2p+A-lambda))",
   sp.Matrix([D(oY[i]) for i in range(3)]) - Cm * oY - sp.Matrix([0, 0, x * (2 * p + A - lam)]))
oc = sp.Matrix([1, mu, mu])
ok("omega^c=(1,mu,mu): omega' - C omega = (0,0,kH - x)", -Cm * oc - sp.Matrix([0, 0, mu * A - x]))
fr, fs = b / a, c / a
chi = lambda v: v * (v - 1) * (v - mu)
ok("frak r' = r^2 - mu r - s", D(fr) - (fr ** 2 - mu * fr - fs))
ok("frak s' = s r - s - x", D(fs) - (fs * fr - fs - x))
ok("(frak s + r^2/2 - r)' = chi(r) - x", D(fs + fr ** 2 / 2 - fr) - (chi(fr) - x))
ok("r'' + (A-3r) r' + chi(r) = x", D(fr, 2) + (A - 3 * fr) * D(fr) + chi(fr) - x)
v = c - (1 - mu) * b
ok("v = c-(1-mu)b: v' = -mu v - x a", D(v) + mu * v + x * a)
ok("(b-a)' = -v", D(b - a) + v)
ok("chi(s) = s^3 - A s^2 + mu s", sp.expand(chi(sp.Symbol('s_')) - (sp.Symbol('s_') ** 3 - A * sp.Symbol('s_') ** 2 + mu * sp.Symbol('s_'))))
# Lemma A0 defects
ph, et_, zt_ = [sp.Function(n)(t) for n in ('phi_', 'etat', 'zetat')]
om3 = sp.Matrix([ph, ph * et_, ph * zt_])
dl = lambda f: sp.diff(f, t) / f
g1 = dl(ph) + mu - et_
g2 = dl(ph) + dl(et_) + 1 - zt_ / et_
g3 = dl(ph) + dl(zt_) + A - x / zt_
res = sp.Matrix([sp.diff(om3[i], t) for i in range(3)]) - Cm * om3 - sp.diag(g1, g2, g3) * om3
ok("omega = phi(1,eta~,zeta~): omega' = C omega + diag(g) omega", res)
R1 = sp.diff(et_, t) - (zt_ + (mu - 1) * et_ - et_ ** 2)
R2 = sp.diff(zt_, t) - (x - (1 + et_) * zt_)
ok("g2 = g1 + R1/eta~", g2 - g1 - R1 / et_)
ok("g3 = g1 + R2/zeta~", g3 - g1 - R2 / zt_)
ok("true pair: g2 = mu - eta", (D(eta) / eta + 1 - zeta / eta) - (mu - eta))
ok("true pair: g3 = mu - eta", (D(zeta) / zeta + A - x / zeta) - (mu - eta))
# general wedge of two variational solutions solves w' = C w
d1 = sp.Matrix([sp.Function(f'd1{i}')(t) for i in range(3)])
d2 = sp.Matrix([sp.Function(f'd2{i}')(t) for i in range(3)])
sub = {}
for dv in (d1, d2):
    rhs = J * dv
    for i in range(3):
        sub[sp.Derivative(dv[i], t)] = rhs[i]
wd = sp.Matrix([d1[0] * d2[1] - d1[1] * d2[0], d1[2] * d2[0] - d1[0] * d2[2], d1[2] * d2[1] - d1[1] * d2[2]])
ok("wedge of two variational solutions: w' = C w", sp.Matrix([sp.diff(wd[i], t).subs(sub) for i in range(3)]) - Cm * wd)
fvec = sp.Matrix([p, pd, yd])
rvec = sp.Matrix([U, p, y - k])
wz = sp.Matrix([fvec[0] * rvec[1] - fvec[1] * rvec[0], fvec[2] * rvec[0] - fvec[0] * rvec[2], fvec[2] * rvec[1] - fvec[1] * rvec[2]])
ok("wedge(z', z-E) = (Sigma, Sigma eta0, Xi)", wz - om)
ok("gamma(gamma-1)(gamma-mu) = f_mu(gamma-1)", sp.expand(chi(sp.Symbol('g_')) - fmu(sp.Symbol('g_') - 1)))

print(" A7. formulas added at the audit (remarks, intermediate proof steps, shifted centre)")
ok("(p/U)' = -Sigma/U^2", D(p / U) + Sg / U ** 2)
ok("(p'/p)' = -Xi/p^2", D(pd / p) + Xg / p ** 2)
ok("(y-zeta)' + (y-zeta) = eta zeta", D(y - zeta) + (y - zeta) - eta * zeta)
ok("(ell-eta)' + mu(ell-eta) = (y-zeta) + eta(eta+1-2mu)", D(ell - eta) + mu * (ell - eta) - ((y - zeta) + eta * (eta + 1 - 2 * mu)))
ok("proof of (gau:Xistate): x Sigma - A Xi + k psi p' = x p^2 + A(p p'' - p'^2) - (x-k)p'",
   x * Sg - A * Xg + k * psi * pd - (x * p ** 2 + A * (p * D(pd) - pd ** 2) - (x - k) * pd))
ok("proof of (gau:Xistate): Xi' = (p'+mu p)p'' + x p^2 - p(x-y)", D(Xg) - ((pd + mu * p) * D(pd) + x * p ** 2 - p * (x - y)))
ok("p'' + A p' + x - k = -mu p   and   p'' + (x-y) = -mu p'", sp.Matrix([D(pd) + A * pd + x - k + mu * p, D(pd) + (x - y) + mu * pd]))
ok("(Sigma eta0)' = U[xp-(x-y)] - p'(y-k)", D(N0) - (U * (x * p - (x - y)) - pd * (y - k)))
ok("kappa-coefficient of W: p[(1+eta0)p' + p zeta0] + p'[p' + (mu-1)p - 2 eta0 p] = Xi",
   p * ((1 + eta0) * pd + p * zeta0) + pd * (pd + (mu - 1) * p - 2 * eta0 * p) - Xg)
ok("zeta0 U - eta0 p - y + k = 0", zeta0 * U - eta0 * p - y + k)
ok("f_mu'' = 6s+4-2mu,  G_mu' = 3s+2-mu,  chi'' = 6s-2A",
   sp.Matrix([sp.diff(fmu(sp.Symbol('s_')), sp.Symbol('s_'), 2) - (6 * sp.Symbol('s_') + 4 - 2 * mu),
              sp.diff(Gmu(sp.Symbol('s_')), sp.Symbol('s_')) - (3 * sp.Symbol('s_') + 2 - mu),
              sp.diff(chi(sp.Symbol('s_')), sp.Symbol('s_'), 2) - (6 * sp.Symbol('s_') - 2 * A)]))
ok("(D+mu)ell = y, (D+1)y = x, (D+A)x = x(p+A)", sp.Matrix([D(ell) + mu * ell - y, D(y) + y - x, D(x) + A * x - x * (p + A)]))
ok("phi=ell, c=A: phi'/phi + c - 1 = y/ell", D(ell) / ell + A - 1 - y / ell)
ok("a'' = mu b + c,  a''' = -mu^2 b - A c - x a", sp.Matrix([D(a, 2) - (mu * b + c), D(a, 3) - (-mu ** 2 * b - A * c - x * a)]))
polar_c = {U: -cc + r * sp.cos(th), p: r * sp.sin(th), y: k + r * s}
ok("shifted centre (-c,0): Sigma_c = r^2 sigma", Sc.subs(polar_c) - r ** 2 * sig)
ok("shifted centre: theta' = -sigma", (((U + cc) * pd - p ** 2) / r ** 2).subs(polar_c) + sig)
ok("shifted centre: (log r)' = sin(cos - mu sin - s)",
   ((((U + cc) * p + p * pd) / r ** 2).subs(polar_c)) - sp.sin(th) * (sp.cos(th) - mu * sp.sin(th) - s))
ok("shifted centre: K + c zeta = r(zeta cos - eta sin) - (y-k)",
   (K + cc * zeta - (r * (zeta * sp.cos(th) - eta * sp.sin(th)) - (y - k))).subs({U: -cc + r * sp.cos(th), p: r * sp.sin(th)}))
ok("polar exponent: pK/Sigma = K sin/(r sigma)  (so that dt = -dtheta/sigma gives K sin/(r sigma^2))",
   (p / Sg).subs(polar) - sp.sin(th) / (r * sig))
fgen, ggen, etg, ztg = [sp.Function(nm_)(t) for nm_ in ('f_gen', 'g_gen', 'eta_gen', 'zeta_gen')]
L2g = lambda h_: sp.diff(h_, t, 2) + (mu - etg) * sp.diff(h_, t) + ztg * h_
Vfg = fgen * sp.diff(ggen, t) - sp.diff(fgen, t) * ggen
ok("Lemma gau:wronskian for arbitrary f, g and arbitrary coefficients eta(t), zeta(t): V' = (eta-mu)V + f L2[g] - g L2[f]",
   sp.expand(sp.diff(Vfg, t) - ((etg - mu) * Vfg + fgen * L2g(ggen) - ggen * L2g(fgen))))
vth = sp.Symbol('vartheta')
Smix, Nmix = (1 - vth) * Sg + vth * Xg, (1 - vth) * K + vth * W
ok("mixed gauge: S' = (eta-mu)S + pN for S = (1-th)Sigma + th Xi, N = (1-th)K + th W   [on O=0]",
   onO(D(Smix) - ((eta - mu) * Smix + p * Nmix)))
p1s, pd1s, U1s = sp.symbols('p1s pd1s U1s')
ok("(p1,p1') wedge (E - q1) = -Sigma1", (p1s * (0 - p1s) - pd1s * (0 - U1s)) + (p1s ** 2 - U1s * pd1s))

# ---------------------------------------------------------------------------------------------
print()
print("PART B: floating-point spot checks (whlib cycles; DOP853; averages by periodic trapezoid rule)")
from scipy.integrate import solve_ivp
from whlib import find_cycle, kH, equilibrium_data


def profile(m, kv, cyc, n=24001):
    T = cyc['T']
    ed = equilibrium_data(m, kv)
    gE = ed['gammaE']

    def f(tt, st):
        Uv, pv, yv, e, z, Kv, Wv = st
        xv = kv * np.exp(Uv)
        return [pv, kv - m * pv - yv, xv - yv, z + (m - 1) * e - e * e, xv - (1 + e) * z,
                -(1 + e) * Kv + (xv * Uv - xv + kv), -(1 + e) * Wv + xv * pv * pv]
    st = np.array([cyc['U0'], 0.0, cyc['y0'], ed['etaE'], ed['zetaE'], 0.0, 0.0])
    for _ in range(6):
        sol = solve_ivp(f, [0, T], st, method='DOP853', rtol=1e-12, atol=1e-14)
        st = sol.y[:, -1].copy()
        st[0], st[1], st[2] = cyc['U0'], 0.0, cyc['y0']
    tg = np.linspace(0, T, n)
    sol = solve_ivp(f, [0, T], st, method='DOP853', rtol=1e-12, atol=1e-14, t_eval=tg, dense_output=False)
    Uv, pv, yv, e, z, Kv, Wv = sol.y
    # strong-stable ratios: integrate backward in time (attracting there)
    from scipy.interpolate import CubicSpline
    Us = CubicSpline(tg, Uv)

    def fb(sg, w):          # sg = -t
        tt = (-sg) % T
        xv = kv * np.exp(Us(tt))
        rr, ss = w
        return [-(rr * rr - m * rr - ss), -(ss * rr - ss - xv)]
    w0 = np.array([gE, gE * (gE - m)])
    for _ in range(8):
        sb = solve_ivp(fb, [0, T], w0, method='DOP853', rtol=1e-11, atol=1e-13)
        w0 = sb.y[:, -1].copy()
    sb = solve_ivp(fb, [0, T], w0, method='DOP853', rtol=1e-11, atol=1e-13, t_eval=np.linspace(0, T, n))
    frr = sb.y[0][::-1].copy()
    fss = sb.y[1][::-1].copy()
    return dict(t=tg, T=T, U=Uv, p=pv, y=yv, x=kv * np.exp(Uv), eta=e, zeta=z, K=Kv, W=Wv, fr=frr, fs=fss)


def avg(f, P):
    return np.trapezoid(f, P['t']) / P['T']


def ddt(f, P):
    return np.gradient(f, P['t'], edge_order=2)


cases = [(0.3, 1.5), (1.0, 1.2), (1.0, 2.0), (0.2, 2.0), (3.0, 1.5), (0.1, 1.1)]
rng = np.random.default_rng(20261003)
numfail = []


def chk(cond, msg):
    if not cond:
        numfail.append(msg)
        print("     NUMERICAL CHECK FAILED:", msg)


for (m, rho) in cases:
    kv = rho * kH(m)
    Av = 1 + m
    lamv = kv / m
    cyc = find_cycle(m, kv)
    P = profile(m, kv, cyc)
    U_, p_, y_, x_, e_, z_, K_, W_ = (P[q] for q in ('U', 'p', 'y', 'x', 'eta', 'zeta', 'K', 'W'))
    pd_ = kv - m * p_ - y_
    pdd_ = -m * pd_ - (x_ - y_)
    Sg_ = p_ ** 2 - U_ * pd_
    Xg_ = pd_ ** 2 - p_ * pdd_
    psi_ = U_ * np.exp(U_) - np.exp(U_) + 1
    lam_c = avg(e_, P) - m
    print(f"\n mu={m} rho={rho} k={kv:.6g} T={P['T']:.6f}  lambda_c=<eta>-mu={lam_c:+.10f}  (monodromy {cyc['lambda_c']:+.10f})")
    Kdef = z_ * U_ - e_ * p_ - y_ + kv
    Wdef = (x_ - (1 + e_) * z_) * p_ - (z_ + (m - 1) * e_ - e_ ** 2) * pd_
    print(f"   consistency: max|K-Kdef|={np.max(np.abs(K_-Kdef)):.1e}  max|W-Wdef|={np.max(np.abs(W_-Wdef)):.1e}  "
          f"max|O|={np.max(np.abs(z_*p_-e_*pd_-(x_-y_))):.1e}")
    print(f"   (a) min Sigma/(U^2+p^2)={np.min(Sg_/(U_**2+p_**2)):.3e}  min Xi/(p^2+p'^2)={np.min(Xg_/(p_**2+pd_**2)):.3e}  "
          f"min K={K_.min():.3e}  min W={W_.min():.3e}")
    chk(Sg_.min() > 0 and Xg_.min() > 0, "gauge positivity")
    r1_ = np.max(np.abs(ddt(Sg_, P) - ((e_ - m) * Sg_ + p_ * K_))) / np.max(np.abs(Sg_))
    r2_ = np.max(np.abs(ddt(Xg_, P) - ((e_ - m) * Xg_ + p_ * W_))) / np.max(np.abs(Xg_))
    print(f"       gauge-equation residuals (finite differences, relative): Sigma {r1_:.1e}, Xi {r2_:.1e}")
    f1 = -avg(p_ * K_ / Sg_, P)
    f2 = -avg(p_ * W_ / Xg_, P)
    f3 = -avg(p_ * (K_ + W_) / (Sg_ + Xg_), P)
    cs = 0.5 * abs(U_.min())
    Scs = p_ ** 2 - (U_ + cs) * pd_
    f4 = -avg(p_ * (K_ + cs * z_) / Scs, P)
    print(f"       -<pK/Sigma>={f1:+.10f}  -<pW/Xi>={f2:+.10f}  (a,b)=(1,1): {f3:+.10f}  shifted c=0.5|Umin|: {f4:+.10f} (min Sigma_c={Scs.min():.2e})")
    for fv in (f1, f2, f3, f4):
        chk(abs(fv - lam_c) < 2e-7 * max(1, abs(lam_c)), "exponent formula")
    eta0_ = (U_ * (x_ - y_) - p_ * (y_ - kv)) / Sg_
    eta1_ = (x_ * Sg_ + kv * psi_ * pd_) / Xg_ - 1
    viol = int(np.sum(((e_ - eta0_) * p_ > 1e-12)))
    print(f"       <eta0>-mu={avg(eta0_, P)-m:+.2e}  <eta1>-mu={avg(eta1_, P)-m:+.2e}  "
          f"max|eta-eta0+pK/Sigma|={np.max(np.abs(e_-eta0_+p_*K_/Sg_)):.1e}  sign violations of (eta-eta0)p<=0: {viol}  min eta0={eta0_.min():+.3f}")
    chk(abs(avg(eta0_, P) - m) < 1e-7 and abs(avg(eta1_, P) - m) < 1e-7 and viol == 0, "eta0/eta1")
    # half-cycle bounds: crest at t=0; trough = interior zero of p
    it = int(np.argmax(p_[1:-1] * p_[2:] < 0) + 1)       # first sign change (falling arc ends)
    I_minus = np.trapezoid((e_ - m)[:it + 1], P['t'][:it + 1])
    I_plus = np.trapezoid((e_ - m)[it:], P['t'][it:])
    rr_ = (y_[0] - kv) / (kv - y_[it])
    print(f"       half-cycle: I_+={I_plus:+.5f} < 2 log r={2*np.log(rr_):+.5f};  I_-={I_minus:+.5f} > {-2*np.log(rr_):+.5f}")
    chk(I_plus < 2 * np.log(rr_) + 1e-3 and I_minus > -2 * np.log(rr_) - 1e-3, "half-cycle bounds")
    # (b) kappa
    kap_ = K_ / Sg_
    rk = np.max(np.abs(ddt(kap_, P) - (kv * psi_ / Sg_ - (1 + 2 * eta0_ - m) * kap_ + p_ * kap_ ** 2))) / np.max(np.abs(kap_))
    zeta0_ = Xg_ / Sg_
    print(f"   (b) kappa-Riccati residual (relative) {rk:.1e};  max|eta-(eta0-p kappa)|={np.max(np.abs(e_-(eta0_-p_*kap_))):.1e}  "
          f"max|zeta-(zeta0-p' kappa)|={np.max(np.abs(z_-(zeta0_-pd_*kap_))):.1e}  max kappa={kap_.max():.4f}")
    # branch integral  lambda_c T = -int (kappa_up - kappa_down) dU
    fall = slice(0, it + 1)
    rise = slice(it, len(p_))
    Ug = np.linspace(U_.min() * 0.999999, U_.max() * 0.999999, 40001)
    k_dn = np.interp(Ug, U_[fall][::-1], kap_[fall][::-1])
    k_up = np.interp(Ug, U_[rise], kap_[rise])
    br = -np.trapezoid(k_up - k_dn, Ug) / P['T']
    print(f"       -(1/T) int (kappa_up - kappa_down) dU = {br:+.8f}   (lambda_c {lam_c:+.8f});  "
          f"fraction of U-range with kappa_up<kappa_down: {np.mean(k_up < k_dn):.2f}")
    chk(abs(br - lam_c) < 2e-4 * max(1, abs(lam_c)), "branch integral")
    # (c) polar
    r_ = np.hypot(U_, p_)
    th_ = np.unwrap(np.arctan2(p_, U_))
    s_ = (y_ - kv) / r_
    sg_ = np.sin(th_) ** 2 + m * np.sin(th_) * np.cos(th_) + s_ * np.cos(th_)
    res_r = np.max(np.abs(ddt(np.log(r_), P) - np.sin(th_) * (np.cos(th_) - m * np.sin(th_) - s_)))
    res_t = np.max(np.abs(ddt(th_, P) + sg_))
    tot = th_[-1] - th_[0]
    # integrate over the angle: |d theta| = sigma dt
    polar_int = -np.trapezoid((K_ / r_) * np.sin(th_) / sg_ ** 2, -th_) / P['T']
    polar_bad = -np.trapezoid((K_ / r_) * np.sin(th_), -th_) / P['T']
    print(f"   (c) polar residuals: (log r)' {res_r:.1e}, theta' {res_t:.1e};  total angle {tot/np.pi:+.6f} pi;  min sigma={sg_.min():.3e}")
    print(f"       -(1/T) oint (K/r) sin(theta) sigma^-2 |dtheta| = {polar_int:+.10f}   (lambda_c {lam_c:+.10f})")
    print(f"       same integral WITHOUT the factor sigma^-2      = {polar_bad:+.10f}   (not equal: the factor is needed)")
    chk(abs(polar_int - lam_c) < 2e-6 * max(1, abs(lam_c)) and abs(tot + 2 * np.pi) < 1e-6, "polar formula")
    # (d) monotone integrals
    q11, q12, q22 = avg(e_ * p_, P), avg(z_ * p_, P), avg(z_ * pd_, P)
    q12b = avg(e_ * pd_, P)
    worst = -np.inf
    for _ in range(200):
        al, be = rng.normal(size=2)
        worst = max(worst, avg((al * p_ + be * pd_) * (al * e_ + be * z_), P) / (al * al + be * be))
    print(f"   (d) <eta p>={q11:+.6e}  <zeta p'>={q22:+.6e}  <zeta p>={q12:+.6e}  <eta p'>={q12b:+.6e}  det={q11*q22-q12**2:+.3e}  "
          f"max over 200 random directions of the normalised form: {worst:+.3e}")
    chk(q11 < 0 and q22 < 0 and q11 * q22 - q12 ** 2 > 0 and worst < 0 and abs(q12 - q12b) < 1e-7, "monotone integrals")
    # (f) Hoelder
    gam = 1 + avg(e_, P)
    mm = gam * (gam - 1) * (gam - m)
    lhs = avg(x_ ** (1 / 3), P) ** 3
    print(f"   (f) <x^(1/3)>^3={lhs:.6f} < gamma(gamma-1)(gamma-mu)={mm:.6f} < k={kv:.6f};  kH={kH(m):.6f};  "
          f"<e^(U/3)>={avg(np.exp(U_/3), P):.6f} < rho^(-1/3)={rho**(-1/3):.6f};  <b>={avg(z_/e_, P):.6f} (gamma-mu={gam-m:.6f})  <c>={avg(x_/z_, P):.6f} (gamma={gam:.6f})")
    chk(lhs < mm < kv and avg(np.exp(U_ / 3), P) < rho ** (-1 / 3), "Hoelder sandwich")
    # (g) LV means
    ell_ = lamv - p_
    print(f"   (g) <x/y>={avg(x_/y_, P):.9f}  <y/ell>={avg(y_/ell_, P):.9f} (mu={m})  <ell>={avg(ell_, P):.9f} (lambda={lamv:.9f})")
    chk(abs(avg(x_ / y_, P) - 1) < 1e-7 and abs(avg(y_ / ell_, P) - m) < 1e-7 and abs(avg(ell_, P) - lamv) < 1e-7, "LV means")
    resL = np.max(np.abs(ddt(ddt(e_, P), P) + (3 * e_ + 2 - m) * ddt(e_, P) + e_ * (e_ + 1) * (e_ + 1 - m) - x_)) / np.max(x_)
    print(f"       second-order equation residual for eta (finite differences, relative to max x): {resL:.1e}")
    # (h) strong-stable identities
    fr_, fs_ = P['fr'], P['fs']
    gam_s = avg(fr_, P)
    loga = -np.concatenate(([0.0], np.cumsum(0.5 * (fr_[1:] + fr_[:-1]) * np.diff(P['t']))))
    ah = np.exp(loga + gam_s * P['t'])
    bh, ch = fr_ * ah, fs_ * ah
    lamB1 = avg(ah * x_ * (lamv - Av - 2 * p_), P) / avg(ch * ell_ + bh * y_ + ah * x_, P)
    lamB2 = avg(ah * (x_ - kH(m)), P) / avg(ch + m * bh + m * ah, P)
    chi_g = gam_s * (gam_s - 1) * (gam_s - m)
    print(f"   (h) <frak r>-A={gam_s-Av:+.10f}  omega^Y formula={lamB1:+.10f}  omega^c formula={lamB2:+.10f}  (lambda_c {lam_c:+.10f})")
    print(f"       chi(gamma)={chi_g:.8f}  <x a^>/<a^>={avg(x_*ah, P)/avg(ah, P):.8f};  <chi(frak r)>={avg(fr_*(fr_-1)*(fr_-m), P):.8f} (k={kv:.8f});  "
          f"min(frak r - 1)={np.min(fr_-1):.3e}  min(frak s-(1-mu)frak r)={np.min(fs_-(1-m)*fr_):.3e}  a^ periodicity defect={abs(ah[-1]/ah[0]-1):.1e}")
    for fv in (gam_s - Av, lamB1, lamB2):
        chk(abs(fv - lam_c) < 5e-6 * max(1, abs(lam_c)), "strong-stable formula")
    chk(abs(avg(fr_ * (fr_ - 1) * (fr_ - m), P) - kv) < 1e-5 * kv and np.min(fr_ - 1) > 0 and np.min(fs_ - (1 - m) * fr_) > 0, "strong-stable ratios")
    # Lemma A0: random positive trial directions never beat the true pair
    best = -np.inf
    Fs = [U_ / max(1, np.max(np.abs(U_))), p_ / np.max(np.abs(p_)), (y_ - kv) / kv,
          np.sin(2 * np.pi * P['t'] / P['T']), np.cos(2 * np.pi * P['t'] / P['T'])]
    dFs = [ddt(F, P) for F in Fs]
    for scale in (0.01, 0.1, 0.5):
        for _ in range(100):
            c1 = scale * rng.normal(size=5)
            c2 = scale * rng.normal(size=5)
            le = sum(ci * F for ci, F in zip(c1, Fs))
            lz = sum(ci * F for ci, F in zip(c2, Fs))
            et_n, zt_n = e_ * np.exp(le), z_ * np.exp(lz)
            dle = ddt(np.log(e_), P) + sum(ci * F for ci, F in zip(c1, dFs))
            dlz = ddt(np.log(z_), P) + sum(ci * F for ci, F in zip(c2, dFs))
            gmin = np.minimum(np.minimum(m - et_n, dle + 1 - zt_n / et_n), dlz + Av - x_ / zt_n)
            best = max(best, avg(gmin, P))
    g_true = avg(np.minimum(np.minimum(m - e_, ddt(np.log(e_), P) + 1 - z_ / e_), ddt(np.log(z_), P) + Av - x_ / z_), P)
    print(f"       comparison lemma: <min g> at the true pair = {g_true:+.8f} (= -lambda_c up to finite differences);  "
          f"largest <min g> over 300 random trial directions = {best:+.8f}  (must be <= {-lam_c:+.8f})")
    chk(best <= -lam_c + 1e-6 and abs(g_true + lam_c) < 1e-4 * max(1, abs(lam_c)), "comparison lemma")
    # (i) Gauss map: simple closed curve (eta, omega), clockwise
    om_ = z_ + e_ + e_ ** 2 / 2
    area = 0.5 * np.trapezoid(e_ * ddt(om_, P) - om_ * ddt(e_, P), P['t'])
    itop = int(np.argmax(om_))
    etad_top = (z_ + (m - 1) * e_ - e_ ** 2)[itop]
    print(f"   (i) signed area of the (eta,omega) curve = {area:+.5e} (negative = clockwise);  eta' at the top = {etad_top:+.4e};  "
          f"min |(eta',zeta')| = {np.min(np.hypot(z_+(m-1)*e_-e_**2, x_-(1+e_)*z_)):.2e}")
    chk(area < 0 and etad_top > 0, "Gauss map orientation")
    # ---- (j) checks added at the audit ------------------------------------------------------
    n_ = len(p_)
    # half-cycle bounds with S = Sigma (the script above tests S = Xi)
    lgS = np.log((U_[0] * (y_[0] - kv)) / (U_[it] * (y_[it] - kv)))
    chk(I_plus < lgS + 1e-3 and I_minus > -lgS - 1e-3, "half-cycle bounds with Sigma")
    # eta1 may be negative?  dynamic bounds of Remark gau:dynamicbounds
    ell_ = lamv - p_
    dyn_ok = (y_ - z_).min() > 0 and ((m > 0.5) or (ell_ - e_).min() > 0)
    tolk = 1e-9 * (1 + kap_.max())
    up = pd_ > 1e-9
    dn = pd_ < -1e-9
    fa = p_ < -1e-9
    kb1 = np.all(kap_[up] < zeta0_[up] / pd_[up] + tolk)
    kb2 = np.all(kap_[dn] < (y_[dn] - zeta0_[dn]) / np.abs(pd_[dn]) + tolk)
    kb3 = (m > 0.5) or np.all(kap_[fa] < (ell_[fa] - eta0_[fa]) / np.abs(p_[fa]) + tolk)
    print(f"   (j) half-cycle with Sigma: I_+={I_plus:+.5f} < {lgS:+.5f}, I_-={I_minus:+.5f} > {-lgS:+.5f};  min eta1={eta1_.min():+.4f};  "
          f"min(y-zeta)={(y_-z_).min():.3e}  min(ell-eta)={(ell_-e_).min():+.3e} (claimed >0 only for mu<=1/2);  kappa upper bounds hold: {bool(kb1 and kb2 and kb3)}")
    chk(dyn_ok and kb1 and kb2 and kb3, "dynamic upper bounds (Remark gau:dynamicbounds)")
    # supporting lines (gau:supportline) and supporting planes (Theorem thm:support) on 240 x n ordered pairs
    sel = np.arange(0, n_ - 1, max(1, (n_ - 1) // 240))
    worst_line, worst_plane = -np.inf, -np.inf
    tt_ = P['t']
    for i in sel:
        far = np.minimum(np.abs(tt_ - tt_[i]), P['T'] - np.abs(tt_ - tt_[i])) > 2e-3 * P['T']
        dq = np.hypot(U_ - U_[i], p_ - p_[i])[far]
        wl = (p_[i] * (p_ - p_[i]) - pd_[i] * (U_ - U_[i]))[far] / (np.hypot(p_[i], pd_[i]) * dq)
        wp = (z_[i] * (U_ - U_[i]) - e_[i] * (p_ - p_[i]) - (y_ - y_[i]))[far] / dq
        worst_line = max(worst_line, wl.max())
        worst_plane = max(worst_plane, wp.max())
    wE = (-(p_ ** 2 - U_ * pd_)).max()
    print(f"       supporting lines: max normalised (p1,p1')^(q-q1) = {worst_line:+.3e} (<0);  planes: max Delta/|dq| = {worst_plane:+.3e} (<0);  "
          f"origin: max(-Sigma) = {wE:+.3e}")
    chk(worst_line < 0 and worst_plane < 0 and wE < 0, "supporting lines / planes")
    # orderings of Proposition prop:monotoneintegrals at equal U and at equal p
    Uf, Ur = U_[:it + 1][::-1], U_[it:]                  # falling arc reversed (U increasing), rising arc
    lev = np.linspace(U_.min(), U_.max(), 402)[1:-1]
    f_at = lambda arr_f, arr_r: (np.interp(lev, Uf, arr_f[:it + 1][::-1]), np.interp(lev, Ur, arr_r[it:]))
    eF, eR = f_at(e_, e_)
    yF, yR = f_at(y_, y_)
    pF, pR = f_at(p_, p_)
    okU = np.all(eF > eR) and np.all(yF > yR) and np.all(eR * (pR - pF) < yF - yR + 1e-9) and np.all(yF - yR < eF * (pR - pF) + 1e-9)
    imax, imin = int(np.argmax(p_)), int(np.argmin(p_))   # crest(0) -> p min -> trough(it) -> p max -> crest
    Aarc = slice(imin, imax + 1)                           # p' > 0, p increasing
    idxB = np.concatenate((np.arange(imax, n_ - 1), np.arange(0, imin + 1)))   # p' < 0, p decreasing
    levp = np.linspace(p_.min(), p_.max(), 402)[1:-1]
    gA = lambda arr: np.interp(levp, p_[Aarc], arr[Aarc])
    gB = lambda arr: np.interp(levp, p_[idxB][::-1], arr[idxB][::-1])
    UA, UB, yA, yB, zA, zB = gA(U_), gB(U_), gA(y_), gB(y_), gA(z_), gB(z_)
    okP = np.all(UB > UA) and np.all(yB > yA) and np.all(zB > zA) and np.all(zA * (UB - UA) < yB - yA + 1e-9) and np.all(yB - yA < zB * (UB - UA) + 1e-9)
    print(f"       orderings at 400 equal-U levels (falling arc: larger eta, y): {bool(okU)};  at 400 equal-p levels (p'<0 point: larger U, y, zeta): {bool(okP)}")
    chk(okU and okP, "branch orderings")
    # Jensen consequence  gamma < gamma_E,  and the compound monodromy: nonnegative, spectral radius m_t
    gE_ = equilibrium_data(m, kv)['gammaE']
    from scipy.interpolate import CubicSpline as _CS
    Usp = _CS(tt_, U_)

    def fC(tq, wv):
        xq = kv * np.exp(Usp(tq % P['T']))
        Wm = wv.reshape(3, 3)
        Cq = np.array([[-m, 1.0, 0.0], [0.0, -1.0, 1.0], [xq, 0.0, -Av]])
        return (Cq @ Wm).ravel()
    sc = solve_ivp(fC, [0, P['T']], np.eye(3).ravel(), method='DOP853', rtol=1e-11, atol=1e-14)
    PhiC = sc.y[:, -1].reshape(3, 3)
    sprC = float(np.max(np.abs(np.linalg.eigvals(PhiC))))
    print(f"       gamma={gam:.6f} < gamma_E={gE_:.6f};  compound monodromy Phi_C(T,0): min entry {PhiC.min():.3e} (>=0), "
          f"log(spr)/T={np.log(sprC)/P['T']:+.9f} (lambda_c {lam_c:+.9f})")
    chk(gam < gE_ and PhiC.min() >= 0 and abs(np.log(sprC) / P['T'] - lam_c) < 1e-6 * max(1, abs(lam_c)), "gamma<gamma_E / compound monodromy")

print()
print("SUMMARY: exact identities failed:", len(FAIL), FAIL)
print("SUMMARY: numerical checks failed:", len(numfail), numfail)
