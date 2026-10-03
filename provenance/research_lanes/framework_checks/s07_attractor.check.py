#!/usr/bin/env python3
"""Checks for Section 7 (s07_attractor.tex) of the v3.0 paper, in the paper's notation.

Part A: exact (sympy) re-derivation of every displayed identity of the section after the renaming
        (source \\cD -> \\Dom, D -> \\Dsc, q -> ell, psi_t -> chi_t, Delta_0,tau_0 -> d_0,nu_0, ...).
Part B: floating-point spot checks on cycles and on orbits of W^u(E), using ROOT/work/whlib.py.

Run:  /Users/admin/.venvs/analysis-py313/bin/python s07_attractor.check.py > s07_attractor.check.out
Single process, about five seconds.  Nothing here is an interval certificate, and nothing here is part of a proof.
"""
import sys
sys.path.insert(0, '/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03/work')
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp, dblquad
from scipy.linalg import expm
import whlib

FAIL = []


def ok(name, expr):
    if isinstance(expr, sp.MatrixBase):
        e = expr.applyfunc(lambda v: sp.simplify(sp.expand(v)))
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


def okb(name, cond, detail=""):
    print(f"  [{'OK' if cond else 'FAIL'}] {name} {detail}")
    if not cond:
        FAIL.append(name)


# ---------------------------------------------------------------------------------------------
print("PART A  exact identities (sympy)")
t = sp.symbols('t', real=True)
mu, k = sp.symbols('mu k', positive=True)
A = 1 + mu
lam = k / mu
U, p, y, eta, zeta = [sp.Function(n)(t) for n in ('U', 'p', 'y', 'eta', 'zeta')]
x = k * sp.exp(U)
flow = {U.diff(t): p, p.diff(t): k - mu * p - y, y.diff(t): x - y}
ric = {eta.diff(t): zeta + (mu - 1) * eta - eta**2, zeta.diff(t): x - (1 + eta) * zeta}


def Dt(e, extra=None):
    """derivative along the flow, with the Riccati system for (eta, zeta)"""
    d = sp.diff(e, t)
    sub = dict(flow)
    sub.update(ric)
    if extra:
        sub.update(extra)
    return sp.expand(d.subs(sub))


# (att:eq:ric) <=> (att:eq:normal)
xg = sp.Function('x')(t)
J = sp.Matrix([[0, 1, 0], [0, -mu, -1], [xg, 0, -1]])
n = sp.Matrix([zeta, -eta, -1])
res = n.diff(t) + J.T * n + (1 + eta) * n
ok("J^T n = (-x, zeta+mu*eta, eta+1)", J.T * n - sp.Matrix([-xg, zeta + mu * eta, eta + 1]))
ok("normal equation: components are zeta'-x+(1+eta)zeta, -eta'+zeta+(mu-1)eta-eta^2, 0",
   res - sp.Matrix([zeta.diff(t) - xg + (1 + eta) * zeta, -eta.diff(t) + zeta + (mu - 1) * eta - eta**2, 0]))

# (att:eq:xyl)
ell = lam - p
ok("x' = (lambda-ell) x", sp.diff(x, t).subs(flow) - (lam - ell) * x)
ok("ell' = y - mu*ell", sp.diff(ell, t).subs(flow) - (y - mu * ell))

# (att:eq:lin) <=> (att:eq:ric)
H1, H2, H3 = [sp.Function(n)(t) for n in ('H1', 'H2', 'H3')]
B = sp.Matrix([[0, 0, xg], [1, mu, 0], [0, 1, 1]])
Hv = sp.Matrix([H1, H2, H3])
lin = {H1.diff(t): (B * Hv)[0], H2.diff(t): (B * Hv)[1], H3.diff(t): (B * Hv)[2]}
e_ = H2 / H3
z_ = H1 / H3
ok("ratios of H'=BH solve the eta-equation", sp.diff(e_, t).subs(lin) - (z_ + (mu - 1) * e_ - e_**2))
ok("ratios of H'=BH solve the zeta-equation", sp.diff(z_, t).subs(lin) - (xg - (1 + e_) * z_))
ricg = {eta.diff(t): zeta + (mu - 1) * eta - eta**2, zeta.diff(t): xg - (1 + eta) * zeta}
Hc = sp.Matrix([zeta * H3, eta * H3, H3])
dHc = Hc.diff(t).subs({H3.diff(t): (1 + eta) * H3}).subs(ricg)
ok("converse: H=(zeta,eta,1)*H3 with H3'=(1+eta)H3 solves H'=BH", dHc - B * Hc)
ok("adjoint: (e^{int(1+eta)} n)' = -J^T (e^{...} n) given the normal equation",
   ((1 + eta) * n + n.diff(t)).subs(ricg) + J.T * n)

# Birkhoff constants
mm, MM = sp.symbols('m M', positive=True)
ok("nu_0 = tanh(d_0/4) = (M-m)/(M+m)",
   sp.tanh(2 * sp.log(MM / mm) / 4).rewrite(sp.exp) - (MM - mm) / (MM + mm))

# Theorem thm:twoorbit, (att:eq:deltaode)
U1, p1, y1, U2, p2, y2 = [sp.Function(n)(t) for n in ('U1', 'p1', 'y1', 'U2', 'p2', 'y2')]
x1 = k * sp.exp(U1)
two = {U1.diff(t): p1, p1.diff(t): k - mu * p1 - y1, y1.diff(t): x1 - y1,
       U2.diff(t): p2, p2.diff(t): k - mu * p2 - y2, y2.diff(t): k * sp.exp(U2) - y2,
       eta.diff(t): zeta + (mu - 1) * eta - eta**2, zeta.diff(t): x1 - (1 + eta) * zeta}
Delta = zeta * (U2 - U1) - eta * (p2 - p1) - (y2 - y1)
Rem = lambda d: sp.exp(d) - 1 - d
ok("Delta' + (1+eta) Delta = -x1*Rem(U2-U1)",
   sp.diff(Delta, t).subs(two) + (1 + eta) * Delta + x1 * Rem(U2 - U1))

# Lemma att:lem:gauges
pp = k - mu * p - y
O = zeta * p - eta * pp - (x - y)
ok("O' = -(1+eta) O", Dt(O) + (1 + eta) * O)
K = zeta * U - eta * p - y + k
psi = U * sp.exp(U) - sp.exp(U) + 1
ok("K' + (1+eta)K = k*psi(U)   (no constraint needed)", Dt(K) + (1 + eta) * K - k * psi)
ok("x*Rem(-U) = k*psi(U)", x * Rem(-U) - k * psi)
etap, zetap = ric[eta.diff(t)], ric[zeta.diff(t)]
W = zetap * p - etap * pp
ok("W' + (1+eta)W - x p^2 = -eta' * O", Dt(W) + (1 + eta) * W - x * p**2 + etap * O)
ok("zeta'' = x p - eta' zeta - (1+eta) zeta'", Dt(zetap) - (x * p - etap * zeta - (1 + eta) * zetap))
ok("eta'' = zeta' + (mu-1) eta' - 2 eta eta'", Dt(etap) - (zetap + (mu - 1) * etap - 2 * eta * etap))
ok("p'' = (eta-mu) p' - zeta p + O", Dt(pp) - ((eta - mu) * pp - zeta * p + O))

# Lemma att:lem:axis: (att:eq:pde) at v = 0
h, u, v, et, ze = sp.symbols('h u v eta zeta')
pde = ze * v - et * (k - mu * v - h) - (k * sp.exp(u) - h)
hsol = sp.solve(pde.subs(v, 0), h)[0]
ok("h(u,0)-k = k(e^u-1)/(1+eta)", hsol - k - k * (sp.exp(u) - 1) / (1 + et))

# Remark att:rem:convexgauges
ok("K = k - [h - u h_u - v h_v]", (ze * u - et * v - h + k) - (k - (h - u * ze - v * (-et))))

# Proposition att:prop:scaling
th, Hq, Hthq = sp.symbols('theta H_q H_thetaq')
Fq = sp.Matrix([v, k - mu * v - Hq])
Fth = sp.Matrix([th * v, k - mu * th * v - Hthq])
wedge = lambda a, b: a[0] * b[1] - a[1] * b[0]
bth = th * (Hq - k) + k - Hthq
ok("F(theta q) ^ F(q) = -v * b_theta(q)", wedge(Fth, Fq) + v * bth)
h1, K1, gq = sp.symbols('h1 K1 gq')   # (i): q1 = theta q, theta*g(q1).q = K1 - k + h1
ok("(i) lower bound algebra gives (1-theta) K(theta q)",
   th * (h1 + (1 - th) * (K1 - k + h1) / th) - th * k + k - h1 - (1 - th) * K1)
hq, Kq = sp.symbols('hq Kq')          # (ii): g(q).q = Kq - k + hq
ok("(ii) upper bound algebra gives (1-theta) K(q)",
   th * (hq - k) + k - (hq - (1 - th) * (Kq - k + hq)) - (1 - th) * Kq)

# Section att:sec:above: linearisation at E
g = sp.symbols('gamma', positive=True)
kE = g * (g - 1) * (g - mu)
JE = sp.Matrix([[0, 1, 0], [0, -mu, -1], [kE, 0, -1]])
zz = sp.symbols('z')
ok("char. polynomial of J_E is z^3 + A z^2 + mu z + k", (zz * sp.eye(3) - JE).det() - (zz**3 + A * zz**2 + mu * zz + kE))
ok("-gamma_E is a root", (zz**3 + A * zz**2 + mu * zz + kE).subs(zz, -g))
ok("cubic at gamma=A equals k_H", (g * (g - 1) * (g - mu)).subs(g, A) - mu * (1 + mu))
eE = sp.Matrix([1, -g, -g * (g - mu)])
ok("J_E e_E = -gamma_E e_E", JE * eE + g * eE)
nE = sp.Matrix([(g - 1) * (g - mu), -(g - 1), -1])
ok("J_E^T n_E = -gamma_E n_E", JE.T * nE + g * nE)
ok("constant Riccati solution at E", sp.Matrix([nE[0] + (mu - 1) * (g - 1) - (g - 1)**2, kE - g * nE[0]]))
ok("eta_E - mu = gamma_E - A", (g - 1) - mu - (g - A))
aa = sp.symbols('a')
ok("m_t exponent: -A + (1+<eta>) = <eta> - mu", -A + (1 + aa) - (aa - mu))

# Proposition prop:crestaxis(ii): wedge transport and the two boundary evaluations
F1, F2, F1u, F1v, F2u, F2v, j1, j2 = sp.symbols('F1 F2 F1u F1v F2u F2v j1 j2')
DF = sp.Matrix([[F1u, F1v], [F2u, F2v]])
Fv = sp.Matrix([F1, F2])
jv = sp.Matrix([j1, j2])
dF = DF * Fv
dj = DF * jv
wdot = dF[0] * jv[1] + Fv[0] * dj[1] - dF[1] * jv[0] - Fv[1] * dj[0]
ok("w=F^j satisfies w' = (div F) w", wdot - (F1u + F2v) * wedge(Fv, jv))
g0, g1, R1p, T1p = sp.symbols('g0 g1 R1p T1p')
ok("w(0) = (0,-g(u)) ^ (1,0) = g(u)", wedge(sp.Matrix([0, -g0]), sp.Matrix([1, 0])) - g0)
Fend = sp.Matrix([0, -g1])
ok("w(T_1) = F ^ ((R_1',0) - T_1' F) = g(R_1 u) R_1'", wedge(Fend, sp.Matrix([R1p, 0]) - T1p * Fend) - g1 * R1p)

# Remark att:rem:planes
ok("w_1' = k - y", Dt(p + mu * U) - (k - y))
ok("w_2' = -k(e^U-1)", Dt(p - (y - k) + mu * U) + k * (sp.exp(U) - 1))
hU, hp = sp.symbols('h_U h_p')   # h_U = zeta, h_p = -eta
jac = lambda rows: sp.Matrix(rows).det()
ok("Jacobian (U,p)->(U,y)", jac([[1, 0], [hU, hp]]).subs({hU: ze, hp: -et}) + et)
ok("Jacobian (U,p)->(p,y)", jac([[0, 1], [hU, hp]]).subs({hU: ze, hp: -et}) + ze)
ok("Jacobian (U,p)->(y,w1)", jac([[hU, hp], [mu, 1]]).subs({hU: ze, hp: -et}) - (ze + mu * et))
ok("Jacobian (U,p)->(U,w2)", jac([[1, 0], [mu - hU, 1 - hp]]).subs({hU: ze, hp: -et}) - (1 + et))
xs = sp.symbols('x', positive=True)
div1 = (-hU / hp - 1).subs({hU: ze, hp: -et})                       # d p_1/dU + d(ke^U-y)/dy
div2 = (-mu + xs / hU - 1).subs({hU: ze, hp: -et})                  # d(k-mu p-y)/dp + x dU/dy - 1
div3 = (xs / (hU - mu * hp) - 1).subs({hU: ze, hp: -et})            # x dU/dy|_{w1} - 1
div4 = ((hU - mu) / (1 - hp)).subs({hU: ze, hp: -et})               # d p_4/dU|_{w2}
ok("div in (U,y) plane = zeta/eta - 1", div1 - (ze / et - 1))
ok("div in (p,y) plane = x/zeta - A", div2 - (xs / ze - (1 + mu)))
ok("div in (y,w1) plane = x/(zeta+mu eta) - 1", div3 - (xs / (ze + mu * et) - 1))
ok("div in (U,w2) plane = (zeta-mu)/(1+eta)", div4 - (ze - mu) / (1 + et))
de = ze + (mu - 1) * et - et**2
dz = xs - (1 + et) * ze
Dlog = [de / et, dz / ze, (dz + mu * de) / (ze + mu * et), de / (1 + et)]
claimed = [ze / et + mu - 1 - et, xs / ze - 1 - et, xs / (ze + mu * et) - (1 + et - mu), (ze - mu) / (1 + et) + mu - et]
for nm, a_, b_, dv in zip(("eta", "zeta", "zeta+mu*eta", "1+eta"), Dlog, claimed, (div1, div2, div3, div4)):
    ok(f"D log({nm}) formula", a_ - b_)
    ok(f"divergence = (eta-mu) + D log({nm})", dv - (et - mu) - a_)

# Lemma lem:absorbing: the double integral and the convolution bound
s_, r_, a_, b_ = sp.symbols('s r a b', positive=True)
tt = sp.Rational(7, 5)
mv, lv = sp.Rational(7, 10), sp.Rational(13, 10)
I1 = sp.integrate(sp.exp(-mv * (tt - s_)) * sp.integrate(sp.exp(-(s_ - r_)) * sp.exp(-lv * (tt - r_)), (r_, 0, s_)), (s_, 0, tt))
I2 = sp.integrate(sp.integrate(sp.exp(-mv * a_ - b_ - lv * (a_ + b_)), (b_, 0, tt - a_)), (a_, 0, tt))
ok("beta(t): change of variables r1=t-s, r2=s-r (mu=7/10, lambda=13/10, t=7/5)", sp.simplify(I1 - I2))
m0 = sp.symbols('mu0', positive=True)
ok("int_0^s e^{-mu0(s-r)} e^{-mu0 r} dr = s e^{-mu0 s}",
   sp.integrate(sp.exp(-m0 * (s_ - r_)) * sp.exp(-m0 * r_), (r_, 0, s_)) - s_ * sp.exp(-m0 * s_))

# ---------------------------------------------------------------------------------------------
print("\nPART B  floating-point spot checks (whlib; DOP853)")


def spiral(mu, k, r0, theta=0.3, n_crests=400, t_max=None):
    """Orbit on W^u(E) with the pull-back Riccati pair and I = int (eta - mu) dt; returns crest events."""
    ed = whlib.equilibrium_data(mu, k)
    vec = np.real(np.exp(1j * theta) * ed["vec_u"])
    z0 = np.array([0.0, 0.0, k]) + r0 * vec / np.linalg.norm(vec)

    def f(t, s):
        U, p, y, eta, zeta, I = s
        x = k * np.exp(U)
        return [p, k - mu * p - y, x - y, zeta + (mu - 1) * eta - eta * eta, x - (1 + eta) * zeta, eta - mu]

    def ev(t, s):
        return s[1]
    ev.direction = -1
    if t_max is None:
        t_max = (np.log(1.0 / r0) + 14.0) / ed["lam_u"].real
    sol = solve_ivp(f, [0, t_max], np.concatenate((z0, [ed["etaE"], ed["zetaE"], 0.0])), method="DOP853",
                    rtol=1e-12, atol=1e-15, events=ev, dense_output=False)
    te = sol.t_events[0]
    ye = sol.y_events[0]
    return te, ye, ed, sol


for (mu_, rho_) in [(1.0, 1.5), (0.3, 1.5), (3.0, 1.3)]:
    k_ = rho_ * whlib.kH(mu_)
    lam_ = k_ / mu_
    print(f"\n mu={mu_} rho={rho_} k={k_:.6g}")
    cyc = whlib.find_cycle(mu_, k_)
    pr = whlib.cycle_profile(mu_, k_, cyc, n=4001)
    T = cyc["T"]
    Uc, pc, yc, ec, zc, xc = pr["U"], pr["p"], pr["y"], pr["eta"], pr["zeta"], pr["x"]
    ppc = k_ - mu_ * pc - yc
    ed = whlib.equilibrium_data(mu_, k_)
    # gauges on the cycle
    Oc = zc * pc - ec * ppc - (xc - yc)
    Kc = zc * Uc - ec * pc - yc + k_
    Wc = (xc - (1 + ec) * zc) * pc - (zc + (mu_ - 1) * ec - ec**2) * ppc
    Lc = (yc - k_) - ed["zetaE"] * Uc + ed["etaE"] * pc
    okb("O = 0 on the cycle", np.max(np.abs(Oc)) < 1e-7 * k_, f"max|O|/k={np.max(np.abs(Oc))/k_:.1e}")
    okb("K>0, W>0, Lambda_E>0 on the cycle", Kc.min() > 0 and Wc.min() > 0 and Lc.min() > 0,
        f"min K={Kc.min():.3e} min W={Wc.min():.3e} min Lambda_E={Lc.min():.3e}")
    inCp = np.sum((Uc >= 0) & (pc <= 0) & (yc - k_ <= 0))
    inCm = np.sum((Uc <= 0) & (pc >= 0) & (yc - k_ >= 0))
    okb("no cycle point in E+C_+ or E+C_-", inCp == 0 and inCm == 0)
    # Riccati bounds (att:eq:constants) with x_-+ from the cycle
    xm, xp_ = xc.min(), xc.max()
    eP = (mu_ - 1 + np.sqrt((mu_ - 1)**2 + 4 * xp_)) / 2
    zM = xm / (1 + eP)
    eM = (mu_ - 1 + np.sqrt((mu_ - 1)**2 + 4 * zM)) / 2
    zP = xp_ / (1 + eM)
    okb("bounds eta_-<=eta<=eta_+, zeta_-<=zeta<=zeta_+, eta_->max(0,mu-1)",
        eM <= ec.min() and ec.max() <= eP and zM <= zc.min() and zc.max() <= zP and eM > max(0, mu_ - 1),
        f"eta in [{ec.min():.4g},{ec.max():.4g}] within [{eM:.4g},{eP:.4g}]; zeta in [{zc.min():.4g},{zc.max():.4g}] within [{zM:.4g},{zP:.4g}]")
    # lambda_c five ways + monodromy
    tq = pr["t"]
    avg = lambda f_: np.trapezoid(f_, tq) / T
    l1 = avg(ec) - mu_
    alt = [avg(zc / ec) - 1, avg(xc / zc) - (1 + mu_), avg(xc / (zc + mu_ * ec)) - 1, avg((zc - mu_) / (1 + ec))]
    okb("lambda_c: five cycle averages agree with log(m_t)/T",
        max(abs(a - l1) for a in alt) < 2e-7 and abs(l1 - cyc["lambda_c"]) < 2e-7,
        f"<eta>-mu={l1:.9f} monodromy={cyc['lambda_c']:.9f} max dev of other averages={max(abs(a-l1) for a in alt):.1e}")
    # full 3x3 monodromy at the crest: M^T n = m_s n, n.f = 0, n.e_s != 0
    s0 = np.concatenate(([cyc["U0"], 0.0, cyc["y0"]], np.eye(3).ravel()))
    sv = solve_ivp(whlib._var_rhs(mu_, k_), [0, T], s0, method="DOP853", rtol=1e-12, atol=1e-14)
    Mono = sv.y[3:12, -1].reshape(3, 3)
    n0 = np.array([zc[0], -ec[0], -1.0])
    ms = np.exp(-(1 + avg(ec)) * T)
    relerr = np.linalg.norm(Mono.T @ n0 - ms * n0) / (ms * np.linalg.norm(n0))
    f0 = np.array([0.0, k_ - cyc["y0"], k_ * np.exp(cyc["U0"]) - cyc["y0"]])
    w_, V_ = np.linalg.eig(Mono)
    i_s = int(np.argmin(np.abs(w_ - ms)))
    es = np.real(V_[:, i_s])
    abserr = np.linalg.norm(Mono.T @ n0 - ms * n0) / (np.linalg.norm(Mono) * np.linalg.norm(n0))
    okb("M^T n(a_0) = m_s n(a_0), m_s = exp(-(1+<eta>)T)", abserr < 1e-9 and relerr < 1e-2,
        f"m_s={ms:.3e}; residual/(|M||n|)={abserr:.1e}; residual/(m_s|n|)={relerr:.1e}")
    okb("n.f = 0 and n.e_s != 0 at the crest",
        abs(n0 @ f0) < 1e-7 * np.linalg.norm(n0) * np.linalg.norm(f0) and abs(n0 @ es) > 1e-3 * np.linalg.norm(n0),
        f"n.f={n0@f0:.1e}  |n.e_s|/|n|={abs(n0@es)/np.linalg.norm(n0):.3f}")
    okb("multipliers {1, m_s, m_t} with m_s m_t = e^{-AT}",
        abs(np.linalg.det(Mono) - np.exp(-(1 + mu_) * T)) < 1e-9 and min(abs(w_ - 1)) < 1e-6
        and min(abs(w_ - cyc["m_t"])) < 1e-6 and min(abs(w_ - ms)) < 1e-9,
        f"eig={np.sort(w_.real)}  det={np.linalg.det(Mono):.6e}  e^(-AT)={np.exp(-(1+mu_)*T):.6e}")
    # axis identity at the cycle crest
    okb("crest: y-k = k(e^U-1)/(1+eta)", abs((yc[0] - k_) - k_ * (np.exp(Uc[0]) - 1) / (1 + ec[0])) < 1e-8 * k_)
    # absorbing box of Lemma lem:absorbing (upper constants) and the lower constant in logarithmic form
    b1 = dblquad(lambda b, a: np.exp(-mu_ * a - b - lam_ * (a + b)), 0, 1, lambda a: 0, lambda a: 1 - a)[0]
    xbar = 2 * lam_ / b1
    ybar = 2 * xbar
    lbar = 2 * ybar / mu_
    mu0 = min(1.0, mu_)
    sgrid = np.linspace(0, 400 / mu0, 400001)
    bad = np.where((lbar + ybar * sgrid) * np.exp(-mu0 * sgrid) > lam_ / 4)[0]
    L = sgrid[bad[-1] + 1]
    log_xlow = np.log(k_ / 4) - lbar * L
    okb("cycle inside the absorbing box (upper bounds; lower bound for x in logarithmic form)",
        xc.max() <= lam_ / b1 and yc.max() <= ybar and (lam_ - pc).max() <= lbar and np.log(xc.min()) >= log_xlow
        and xbar > 2 * k_,
        f"x_max={xc.max():.4g} <= lambda/beta_1={lam_/b1:.4g}; log x_min={np.log(xc.min()):.3f} >= {log_xlow:.3g}")
    # Birkhoff constants and the pull-back rate (att:eq:pullback), driving function = the cycle
    Bm = np.array([[0, 0, xm], [1, mu_, 0], [0, 1, 1.0]])
    Bp = np.array([[0, 0, xp_], [1, mu_, 0], [0, 1, 1.0]])
    mlow, mup = expm(Bm).min(), expm(Bp).max()
    d0 = 2 * np.log(mup / mlow)
    nu0 = np.tanh(d0 / 4)
    okb("underline m > 0 and nu_0 = (M-m)/(M+m) < 1", mlow > 0 and abs(nu0 - (mup - mlow) / (mup + mlow)) < 1e-12 and nu0 < 1,
        f"m={mlow:.4g} M={mup:.4g} d_0={d0:.4g} nu_0={nu0:.6f}")

    def fr(t, s):
        U, p, y = s[:3]
        x = k_ * np.exp(U)
        out = [p, k_ - mu_ * p - y, x - y]
        for i in range(3, len(s), 2):
            e, z = s[i], s[i + 1]
            out += [z + (mu_ - 1) * e - e * e, x - (1 + e) * z]
        return out
    data = [(ec[0], zc[0]), (0.0, 0.0), (50.0, 0.0), (0.0, 50.0), (50.0, 50.0)]
    s0 = [cyc["U0"], 0.0, cyc["y0"]] + [c for d in data for c in d]
    tgrid = [1.0, 2.0, 5.0, 10.0]
    sr = solve_ivp(fr, [0, 10.0], s0, method="DOP853", rtol=1e-12, atol=1e-14, t_eval=tgrid)
    worst = 0.0
    good = True
    for it, tv in enumerate(tgrid):
        ref_e, ref_z = sr.y[3, it], sr.y[4, it]
        bound = d0 * nu0 ** (np.floor(tv) - 1)
        for i in range(5, len(s0), 2):
            gap = max(abs(np.log(sr.y[i, it] / ref_e)), abs(np.log(sr.y[i + 1, it] / ref_z)))
            worst = max(worst, gap / bound)
            good &= gap <= bound
    okb("pull-back estimate (att:eq:pullback) from data (0,0),(50,0),(0,50),(50,50)", good, f"max gap/bound={worst:.3g}")

    # orbits on W^u(E): crest events
    te, ye, ed, sol = spiral(mu_, k_, 1e-6)
    uc, ycr, ecr, zcr, Icr = ye[:, 0], ye[:, 2], ye[:, 3], ye[:, 4], ye[:, 5]
    keep = uc > 1e-3 * cyc["U0"]
    te, uc, ycr, ecr, zcr, Icr = te[keep], uc[keep], ycr[keep], ecr[keep], zcr[keep], Icr[keep]
    gst = ycr - k_
    okb("W^u(E) crests: y-k = k(e^u-1)/(1+eta)  (Lemma att:lem:axis)",
        np.max(np.abs(gst - k_ * (np.exp(uc) - 1) / (1 + ecr)) / gst) < 1e-6,
        f"{len(uc)} crests, max rel. err {np.max(np.abs(gst - k_*(np.exp(uc)-1)/(1+ecr))/gst):.1e}")
    inc = np.diff(uc)
    growing = uc[:-1] < cyc["U0"] * (1 - 1e-7)
    okb("R_1(u) > u inside the cycle (att:eq:outward): crests increase", np.all(inc[growing] > 0),
        f"first crest {uc[0]:.4g}, last {uc[-1]:.9g}, cycle crest {cyc['U0']:.9g}")
    # supporting planes (att:eq:support) between crest points of different turns and the scaling inequalities
    Kcr = zcr * uc - ycr + k_
    sel = np.where(uc < cyc["U0"] * (1 - 1e-5))[0]
    viol_sup = 0
    viol_i = 0
    viol_ii = 0
    npairs = 0
    for a in sel:
        for b in sel:
            if a == b or abs(uc[a] - uc[b]) < 1e-3 * uc[b]:
                continue
            npairs += 1
            # supporting plane at a, evaluated at b (both have p = 0)
            if not (ycr[b] - ycr[a] > zcr[a] * (uc[b] - uc[a])):
                viol_sup += 1
            thv = uc[a] / uc[b]                      # theta q_b = q_a on the crest axis
            bval = thv * (ycr[b] - k_) + k_ - ycr[a]  # b_theta(q_b)
            if not (bval <= (1 - thv) * Kcr[b] + 1e-9 * k_):
                viol_ii += 1
            if thv < 1 and not (bval >= (1 - thv) * Kcr[a] - 1e-9 * k_ and bval > 0):
                viol_i += 1
            if thv > 1 and not (bval < 0):
                viol_ii += 1
    okb("supporting planes between crest points of different turns", viol_sup == 0, f"{npairs} ordered pairs")
    okb("scaling: b_theta >= (1-theta)K(theta q) > 0 for theta<1; b_theta <= (1-theta)K(q), so b_theta<0 for theta>1",
        viol_i == 0 and viol_ii == 0, f"{npairs} ordered pairs on the crest half-axis")
    # return-map derivative (att:eq:returnderivative) by central differences between neighbouring orbits of W^u(E)
    epsr = 2e-4
    tm, ym, _, _ = spiral(mu_, k_, 1e-6 * (1 - epsr))
    tp, yp, _, _ = spiral(mu_, k_, 1e-6 * (1 + epsr))
    um, up = ym[:, 0], yp[:, 0]
    full = ye[:, 0]
    nmin = min(len(um), len(up), len(full))
    errs = []
    for i in range(nmin - 1):
        u0_ = full[i]
        if u0_ < 0.02 * cyc["U0"] or u0_ > 0.999 * cyc["U0"]:
            continue
        fd = (up[i + 1] - um[i + 1]) / (up[i] - um[i])
        iota = ye[i + 1, 5] - ye[i, 5]
        formula = (ye[i, 2] - k_) * np.exp(iota) / (ye[i + 1, 2] - k_)
        errs.append(abs(fd / formula - 1))
    okb("R_1'(u) = g(u) e^{iota(u)} / g(R_1 u) against central differences", len(errs) > 3 and max(errs) < 2e-5,
        f"{len(errs)} crests, max rel. err {max(errs):.1e}")
    # Corollary att:cor:dulac(i): eta - mu at E
    okb("eta_E - mu = gamma_E - A > 0", abs(ed["etaE"] - mu_ - (ed["gammaE"] - 1 - mu_)) < 1e-12 and ed["etaE"] > mu_,
        f"eta_E-mu={ed['etaE']-mu_:.6f} = 2 Re(lambda_u)={2*ed['lam_u'].real:.6f}")

print("\nSUMMARY:", "ALL CHECKS PASSED" if not FAIL else f"{len(FAIL)} FAILED: {FAIL}")
