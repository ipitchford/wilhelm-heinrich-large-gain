#!/usr/bin/env python3
"""Re-check of every displayed formula in s02_adjoint, s03_support, s04_geometry and appA_hopf
after the v3.0 renaming (psi, Rem, Lambda_E, bold H, xi-variations, chi/calB/calC in App. A).

Part 1: exact symbolic identities (sympy).  Part 2: floating-point spot check with whlib.
Exit status is nonzero if any check fails."""
import sys
import sympy as sp

fails = []
def check(name, expr):
    if isinstance(expr, sp.MatrixBase):
        ok = all(sp.simplify(e) == 0 for e in expr)
    else:
        ok = sp.simplify(expr) == 0
    print(("PASS " if ok else "FAIL ") + name)
    if not ok:
        fails.append(name)

mu, k = sp.symbols('mu k', positive=True)
A = 1 + mu
U, p, y, eta, zeta = sp.symbols('U p y eta zeta', real=True)
x = k * sp.exp(U)
psi = lambda u: u * sp.exp(u) - sp.exp(u) + 1          # \psiU
Rem = lambda d: sp.exp(d) - 1 - d                      # \Rem

# flow derivatives of the state and of the adjoint ratios (eqs logsystem, projective)
dU, dp, dy = p, k - mu * p - y, x - y
deta = zeta + (mu - 1) * eta - eta**2
dzeta = x - (1 + eta) * zeta
def Dt(f):
    return (sp.diff(f, U) * dU + sp.diff(f, p) * dp + sp.diff(f, y) * dy
            + sp.diff(f, eta) * deta + sp.diff(f, zeta) * dzeta)

# ---------- s02 ----------
J = sp.Matrix([[0, 1, 0], [0, -mu, -1], [x, 0, -1]])
S = sp.diag(1, -1, -1)
B = sp.Matrix([[0, 0, x], [1, mu, 0], [0, 1, 1]])
check("s02 eq:adjointmatrices  B = -S J^T S", (B + S * J.T * S))
H1, H2, H3 = sp.symbols('H1 H2 H3', positive=True)
Hv = sp.Matrix([H1, H2, H3]); dH = B * Hv
e_, z_ = H2 / H3, H1 / H3
check("s02 eq:projective  eta' from H'=BH",
      (dH[1] * H3 - H2 * dH[2]) / H3**2 - (z_ + (mu - 1) * e_ - e_**2))
check("s02 eq:projective  zeta' from H'=BH",
      (dH[0] * H3 - H1 * dH[2]) / H3**2 - (x - (1 + e_) * z_))
check("s02 H3' = (1+eta) H3", dH[2] - (1 + e_) * H3)
n = sp.Matrix([zeta, -eta, -1]); dn = sp.Matrix([dzeta, -deta, 0])
check("s02 n' + J^T n + (1+eta) n = 0", dn + J.T * n + (1 + eta) * n)
O = zeta * p - eta * dp - (x - y)
check("s02 O' + (1+eta) O = 0 (O = left side of eq:orthogonality)", Dt(O) + (1 + eta) * O)
ddp = Dt(dp)
check("s02 p'' = -mu p' - x + y", ddp - (-mu * dp - x + y))
check("s02 eq:tangent modulo O:  p'' - [(eta-mu)p' - zeta p] = O", ddp - ((eta - mu) * dp - zeta * p) - O)
t = sp.symbols('t')
u = sp.Function('u')(t); et = sp.Function('eta')(t); ze = sp.Function('zeta')(t); xt = sp.Function('x')(t)
inner = u.diff(t, 2) + (mu - et) * u.diff(t) + ze * u
fact = inner.diff(t) + (1 + et) * inner
fact = fact.subs({et.diff(t): ze + (mu - 1) * et - et**2, ze.diff(t): xt - (1 + et) * ze})
check("s02 eq:factorisation", sp.expand(fact - (u.diff(t, 3) + A * u.diff(t, 2) + mu * u.diff(t) + xt * u)))
s = sp.symbols('s')
f_mu = s * (s + 1) * (s + 1 - mu)
check("s02 eq:gainprimitive", Dt(zeta + eta + eta**2 / 2) - (x - f_mu.subs(s, eta)))
xp, xm = sp.symbols('x_plus x_minus', positive=True)
eta_p = (mu - 1 + sp.sqrt((mu - 1)**2 + 4 * xp)) / 2
check("s02 lem:projectivebounds  eta_+ solves eta(eta+1-mu)=x_+", sp.expand(eta_p * (eta_p + 1 - mu) - xp))
zeta_m = xm / (1 + eta_p)
eta_m = (mu - 1 + sp.sqrt((mu - 1)**2 + 4 * zeta_m)) / 2
check("s02 lem:projectivebounds  eta_- solves eta(eta+1-mu)=zeta_-", sp.expand(eta_m * (eta_m + 1 - mu) - zeta_m))

# ---------- s03 ----------
U2, p2, y2 = sp.symbols('U2 p2 y2', real=True)
x2 = k * sp.exp(U2)
Delta = zeta * (U2 - U) - eta * (p2 - p) - (y2 - y)
dDelta = Dt(Delta) + sp.diff(Delta, U2) * p2 + sp.diff(Delta, p2) * (k - mu * p2 - y2) + sp.diff(Delta, y2) * (x2 - y2)
check("s03 eq:supportidentity  Delta'+(1+eta)Delta = -x_1 Rem(U_2-U_1)", dDelta + (1 + eta) * Delta + x * Rem(U2 - U))
check("s03 x_2 - x_1 - x_1(U_2-U_1) = x_1 Rem(U_2-U_1)", x2 - x - x * (U2 - U) - x * Rem(U2 - U))
K = zeta * U - eta * p - y + k
check("s03 eq:positiveK  K'+(1+eta)K = xU - x + k", Dt(K) + (1 + eta) * K - (x * U - x + k))
check("s03 eq:positiveK  xU - x + k = k psi(U)", x * U - x + k - k * psi(U))
check("s03 K = -Delta with second solution E", K + Delta.subs({U2: 0, p2: 0, y2: k}))
check("s03 x Rem(-U) = k psi(U)", x * Rem(-U) - k * psi(U))
check("s03 Lambda_E = -Delta with base E (ratios eta_E, zeta_E) and second solution (U,p,y)",
      ((y - k) - zeta * U + eta * p) + (zeta * (U - 0) - eta * (p - 0) - (y - k)))
g = sp.symbols('gamma_E', positive=True)
etaE, zetaE = g - 1, (g - 1) * (g - mu)
kE = g * (g - 1) * (g - mu)
BE = B.subs(sp.exp(U), 1).subs(k, kE)
check("s03 equilibrium adjoint eigenvector  B_E (zeta_E,eta_E,1)^T = gamma_E (...)",
      BE * sp.Matrix([zetaE, etaE, 1]) - g * sp.Matrix([zetaE, etaE, 1]))
check("s03 gamma_E at Hopf value: A(A-1)(A-mu) = mu(1+mu)", A * (A - 1) * (A - mu) - mu * (1 + mu))
LamE = (y - k) - zetaE * U + etaE * p
dLamE = sp.diff(LamE, U) * dU + sp.diff(LamE, p) * dp + sp.diff(LamE, y) * dy
check("s03 eq:LambdaE  Lambda_E' = -gamma_E Lambda_E + k Rem(U)",
      (dLamE - (-g * LamE + k * Rem(U))).subs(k, kE))
check("s03 Lambda_E = -(zeta_E U - eta_E p - y + k)", LamE + (zetaE * U - etaE * p - y + k))
W = dzeta * p - deta * dp
check("s03 lem:curvature  W'+(1+eta)W - x p^2 = -eta' O", Dt(W) + (1 + eta) * W - x * p**2 + deta * O)

# ---------- s04 ----------
O0 = O.subs(p, 0)      # -eta (k-y) - (x-y) = 0
xsol = sp.solve(O0, sp.exp(U))[0] * k   # x on the constraint at p=0
check("s04 eq:crestidentities  x-k = (1+eta)(y-k)", xsol - k - (1 + eta) * (y - k))
check("s04 eq:crestidentities  x-y = eta(y-k)", xsol - y - eta * (y - k))
check("s04 y'' = x p at y'=0", (Dt(dy) - x * p + dy))       # y'' = x p - y'
lam = k / mu
Ychem = 1 - p / lam
check("s04 Y'' - zeta p/lambda = -[(eta-mu) p' + O]/lambda  (so Y'' = zeta p/lambda where p'=0 on a cycle)",
      Dt(Dt(Ychem)) - zeta * p / lam + ((eta - mu) * dp + O) / lam)
check("s04 y' = zeta p - eta p' - O  (so y' = zeta p where p'=0 on a cycle)", dy - (zeta * p - eta * dp - O))
# compound matrix: w = (om12, -om13, -om23), om = f ^ v, v' = J v
v1, v2, v3 = sp.symbols('v1 v2 v3', real=True)
f = sp.Matrix([dU, dp, dy]); v = sp.Matrix([v1, v2, v3])
def wedge(a, b):
    return sp.Matrix([a[0] * b[1] - a[1] * b[0], -(a[0] * b[2] - a[2] * b[0]), -(a[1] * b[2] - a[2] * b[1])])
w = wedge(f, v)
dv = J * v; df = J * f
dw = wedge(df, v) + wedge(f, dv)
C = sp.Matrix([[-mu, 1, 0], [0, -1, 1], [x, 0, -A]])
check("s04 eq:compound  w' = C w", sp.expand(dw - C * w))
xiU, xiy = sp.symbols('xi_U xi_y', real=True)
w0 = wedge(f.subs(p, 0), sp.Matrix([xiU, 0, xiy]))
check("s04 initial wedge at a crest", w0 - sp.Matrix([(y - k) * xiU, (x - y) * xiU, (y - k) * xiy]))
m11, m12, m21, m22 = sp.symbols('m11 m12 m21 m22')
M = sp.Matrix([[m11, m12], [m21, m22]])
Mz = sp.diag(1, 1 / k) * M * sp.diag(1, k)
check("s04 similarity  entries of M_(U,z)", Mz - sp.Matrix([[m11, k * m12], [m21 / k, m22]]))
check("s04 similarity  det and trace invariant", sp.Matrix([Mz.det() - M.det(), Mz.trace() - M.trace()]))
check("s04 AT bound  2 pi (sqrt(mu)+1/sqrt(mu)) = A * 2 pi / sqrt(mu)",
      2 * sp.pi * (sp.sqrt(mu) + 1 / sp.sqrt(mu)) - A * 2 * sp.pi / sp.sqrt(mu))
c0, c1, th = sp.symbols('c0 c1 theta', real=True)
Uh = c0 + c1 * sp.cos(sp.sqrt(mu) * t + th)
check("s04 equality case: U'''+mu U' = 0 for the first harmonic", Uh.diff(t, 3) + mu * Uh.diff(t))

# ---------- appendix A ----------
sig, rho = sp.symbols('sigma rho', positive=True)
tau = sp.symbols('tau')
Uf = sp.Function('U')
lhs_t = (Uf(sp.sqrt(mu) * t).diff(t, 3) + A * Uf(sp.sqrt(mu) * t).diff(t, 2) + mu * Uf(sp.sqrt(mu) * t).diff(t)).doit()
# coefficients of U_ttt, U_tt, U_t after division by mu^(3/2)
check("appA sigma*rho = k / mu^(3/2)", ((1 + mu) / sp.sqrt(mu)) * (k / (mu * (1 + mu))) - k / mu**sp.Rational(3, 2))
check("appA coefficient of U_tautau", A * mu / mu**sp.Rational(3, 2) - (1 + mu) / sp.sqrt(mu))
nu = sp.symbols('nu')
P = nu**3 + sig * nu**2 + nu + sig * rho
check("appA characteristic polynomial at rho=1", sp.expand(P.subs(rho, 1) - (nu + sig) * (nu**2 + 1)))
dnu = (-sp.diff(P, rho) / sp.diff(P, nu)).subs({nu: sp.I, rho: 1})
check("appA d Re nu / d rho = sigma / (2 (1+sigma^2))", sp.re(sp.simplify(dnu)) - sig / (2 * (1 + sig**2)))
J0 = sp.Matrix([[0, 1, 0], [0, 0, 1], [-sig, -1, -sig]])
q = sp.Matrix([1, sp.I, -1]); qb = q.conjugate()
chi = sp.Matrix([[sp.I * sig, sig + sp.I, 1]]) / (-2 + 2 * sp.I * sig)
check("appA J0 q = i q", J0 * q - sp.I * q)
check("appA chi J0 = i chi", chi * J0 - sp.I * chi)
check("appA chi q = 1", (chi * q)[0] - 1)
e3 = sp.Matrix([0, 0, 1])
Bf = lambda a, b: -sig * a[0] * b[0] * e3
Cf = lambda a, b, c: -sig * a[0] * b[0] * c[0] * e3
h11 = J0.inv() * Bf(q, qb)
check("appA J0^{-1} B(q,qbar) = (1,0,0)^T", h11 - sp.Matrix([1, 0, 0]))
h20 = (2 * sp.I * sp.eye(3) - J0).inv() * Bf(q, q)
check("appA first component of (2iI-J0)^{-1} B(q,q) = sigma/[3(sigma+2i)]", h20[0] - sig / (3 * (sig + 2 * sp.I)))
G21 = (chi * (Cf(q, q, qb) - 2 * Bf(q, h11) + Bf(qb, h20)))[0]
check("appA G21", G21 - sig * (sig + 3 * sp.I) / (3 * (sig + 2 * sp.I) * (-1 + sp.I * sig)))
check("appA l1 = Re G21 / 2", sp.re(sp.simplify(G21)) / 2 + sig / ((sig**2 + 1) * (sig**2 + 4)))
check("appA sigma >= 2: sigma - 2 = (1-sqrt(mu))^2/sqrt(mu)", (1 + mu) / sp.sqrt(mu) - 2 - (1 - sp.sqrt(mu))**2 / sp.sqrt(mu))

# ---------- floating-point spot check (whlib) ----------
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import numpy as np
import whlib

def fcheck(name, ok, detail=""):
    print(("PASS " if ok else "FAIL ") + name + ("   " + detail if detail else ""))
    if not ok:
        fails.append(name)

for mu_, rho_ in [(0.1, 1.5), (1.0, 2.0), (3.0, 1.3)]:
    k_ = rho_ * whlib.kH(mu_)
    c = whlib.find_cycle(mu_, k_)
    pr = whlib.cycle_profile(mu_, k_, c)
    T, M = c["T"], c["M"]
    A_ = 1 + mu_
    tag = f"[mu={mu_}, rho={rho_}] "
    Mz_ = np.diag([1, 1 / k_]) @ M @ np.diag([1, k_])
    ev, evz = np.sort(np.linalg.eigvals(M).real), np.sort(np.linalg.eigvals(Mz_).real)
    avg_eta = float(np.trapezoid(pr["eta"], pr["t"]) / T)
    fcheck(tag + "M entrywise positive in (U,y) and in (U,z)", bool((M > 0).all() and (Mz_ > 0).all()),
           f"M={M.tolist()}")
    # the small eigenvalue m_s is lost to cancellation when exp(-AT) is below rounding; compare with an absolute floor
    fcheck(tag + "eig M_(U,z) = eig M", bool(np.allclose(ev, evz, rtol=1e-9, atol=1e-12 * np.abs(M).max())),
           f"eig M={ev.tolist()}, eig M_(U,z)={evz.tolist()}")
    fcheck(tag + "m_t = exp((<eta>-mu)T)", bool(abs(ev[1] / np.exp((avg_eta - mu_) * T) - 1) < 1e-6),
           f"m_t={ev[1]:.12g}, <eta>-mu={avg_eta - mu_:.9f}")
    if np.exp(-A_ * T) > 1e-9:
        fcheck(tag + "m_s = exp(-(1+<eta>)T)", bool(abs(ev[0] / np.exp(-(1 + avg_eta) * T) - 1) < 1e-4),
               f"m_s={ev[0]:.6e}")
    else:
        print("SKIP " + tag + "m_s = exp(-(1+<eta>)T): exp(-AT) is below the rounding level of M")
    # floating determinants this small carry no relative accuracy; test with an absolute tolerance
    fcheck(tag + "det M = exp(-AT) (absolute tolerance), exp(-AT) < exp(-4 pi), T > 2 pi/sqrt(mu)",
           bool(abs(np.linalg.det(M) - np.exp(-A_ * T)) < 1e-9 * np.abs(M).max() ** 2
                and np.exp(-A_ * T) < np.exp(-4 * np.pi) and T > 2 * np.pi / np.sqrt(mu_)),
           f"T={T:.6f}, det M={np.linalg.det(M):.3e}, exp(-AT)={np.exp(-A_ * T):.3e}")
    U_, p_, y_, x_, et_, ze_ = (pr[s_] for s_ in ("U", "p", "y", "x", "eta", "zeta"))
    dp_ = k_ - mu_ * p_ - y_
    ed = whlib.equilibrium_data(mu_, k_)
    LamE_ = (y_ - k_) - ed["zetaE"] * U_ + ed["etaE"] * p_
    scale = np.max(np.abs(x_ - y_)) + 1e-300
    fcheck(tag + "orthogonality zeta p - eta p' - (x-y) = 0", bool(np.max(np.abs(ze_ * p_ - et_ * dp_ - (x_ - y_))) / scale < 1e-6))
    fcheck(tag + "K > 0, W > 0, Lambda_E > 0 on the cycle",
           bool(pr["K"].min() > 0 and pr["W"].min() > 0 and LamE_.min() > 0),
           f"min K={pr['K'].min():.3e}, min W={pr['W'].min():.3e}, min Lambda_E={LamE_.min():.3e}")
    fcheck(tag + "<p'^2> = mu <p^2>", bool(abs(np.trapezoid(dp_**2, pr['t']) / (mu_ * np.trapezoid(p_**2, pr['t'])) - 1) < 1e-6))
    fcheck(tag + "<x> = k = <f_mu(eta)>",
           bool(abs(np.trapezoid(x_, pr['t']) / T / k_ - 1) < 1e-6
                and abs(np.trapezoid(et_ * (et_ + 1) * (et_ + 1 - mu_), pr['t']) / T / k_ - 1) < 1e-5))
    in_Cp = (U_ >= 0) & (p_ <= 0) & (y_ - k_ <= 0)
    in_Cm = (U_ <= 0) & (p_ >= 0) & (y_ - k_ >= 0)
    fcheck(tag + "no sample in C_+ or C_-", bool(not in_Cp.any() and not in_Cm.any()))

print()
print("ALL CHECKS PASSED" if not fails else f"{len(fails)} CHECK(S) FAILED: {fails}")
sys.exit(1 if fails else 0)
