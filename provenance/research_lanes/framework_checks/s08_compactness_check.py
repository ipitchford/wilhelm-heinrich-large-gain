#!/usr/bin/env python3
"""Formula and inequality re-check for paper section s08_compactness.tex (v3.0), written at the audit of 3 October 2026.

Part A (exact, sympy): the identities used in Lemma lem:statebox, in the preliminaries of Lemma lem:eqseparation
  and in Lemma lem:speed.
Part B (floating point): the two-filter tail bound on a grid; on four cycles: every inequality of the state-box
  argument (ratio bounds, upper bounds, the excursion below X = 1/2, the lower bound beta_0), the cone exclusion and
  the angle bound, the quadratic form V and its derivative along the flow, the separation constants for the
  one-point parameter set Q = {(mu,k)}, the speed relations, the arclength bound and the period bounds.

Floating point is evidence only; nothing here is an interval certificate.
Run:  /Users/admin/.venvs/analysis-py313/bin/python -W ignore s08_compactness_check.py > s08_compactness_check.out
(single process, about one minute)
"""
import sys
import numpy as np
import sympy as sp
from scipy.linalg import expm, solve_continuous_lyapunov
from scipy.integrate import quad_vec

sys.path.insert(0, '/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03/work')
import whlib  # noqa: E402

FAIL = 0


def report(name, ok, extra=''):
    global FAIL
    FAIL += (not ok)
    print(f"{'PASS' if ok else 'FAIL'}  {name}  {extra}")


def Z(e):
    if isinstance(e, sp.MatrixBase):
        return all(sp.simplify(q) == 0 for q in e)
    return sp.simplify(e) == 0


# ----------------------------------------------------------------------------------------------
# Part A: exact identities
# ----------------------------------------------------------------------------------------------
def part_a():
    X, Y, Zc, mu, lam, k, s, g = sp.symbols('X Y Z mu lambda k s gamma', positive=True)
    A = 1 + mu
    # filters (eq:filters): X' = lam X (1-Y), Y' = mu (Z-Y), Z' = X - Z
    dX, dY, dZ = lam * X * (1 - Y), mu * (Zc - Y), X - Zc
    R, S = Zc / X, Y / X
    dR = sp.diff(R, X) * dX + sp.diff(R, Zc) * dZ
    dS = sp.diff(S, X) * dX + sp.diff(S, Y) * dY
    report("A01 R=Z/X: R' = 1 - R(1 + lam(1-Y))", Z(dR - (1 - R * (1 + lam * (1 - Y)))))
    report("A02 S=Y/X: S' = mu R - S(mu + lam(1-Y))", Z(dS - (mu * R - S * (mu + lam * (1 - Y)))))
    report("A03 at X'=0 (Y=1): X = 1/S", Z((1 / S).subs(Y, 1) - X))
    # two-filter tail: I(s) = int_0^s exp(-mu(s-sig)-sig) dsig solves I' = -mu I + e^{-s}, I(0)=0
    Ic = (sp.exp(-s) - sp.exp(-mu * s)) / (mu - 1)
    report("A04 I(s) = (e^{-s}-e^{-mu s})/(mu-1): I' = -mu I + e^{-s}, I(0) = 0  (mu != 1)",
           Z(sp.diff(Ic, s) + mu * Ic - sp.exp(-s)) and Z(Ic.subs(s, 0)))
    I1 = s * sp.exp(-s)
    report("A05 mu = 1: I(s) = s e^{-s} solves I' = -I + e^{-s}, I(0) = 0", Z(sp.diff(I1, s) + I1 - sp.exp(-s)) and Z(I1.subs(s, 0)))
    th1 = sp.exp(-mu * s) + mu * Ic
    th2 = sp.exp(-s) + Ic
    report("A06 theta_mu: e^{-mu s} + mu I = e^{-s} + I", Z(th1 - th2))
    report("A07 theta_mu is the second filter started at 1 with first filter e^{-s}: theta' = mu(e^{-s} - theta), theta(0) = 1",
           Z(sp.diff(th1, s) - mu * (sp.exp(-s) - th1)) and Z(th1.subs(s, 0) - 1))
    report("A08 theta_mu = survival function of Exp(1)+Exp(mu): -theta' = density mu(e^{-s}-e^{-mu s})/(mu-1)",
           Z(-sp.diff(th2, s) - mu * (sp.exp(-s) - sp.exp(-mu * s)) / (mu - 1)))
    v = sp.symbols('v', positive=True)
    report("A09 d/dv[(1+v)e^{-v/2}] = (1-v)e^{-v/2}/2  (maximum at v = 1, value 2e^{-1/2} < 2)",
           Z(sp.diff((1 + v) * sp.exp(-v / 2), v) - (1 - v) * sp.exp(-v / 2) / 2) and bool(2 * sp.exp(sp.Rational(-1, 2)) < 2))
    L0 = 2 / mu * sp.log(8 * X)          # X plays the role of M_0, mu of mu_wedge
    report("A10 exp(-mu_w L_0/2) = 1/(8 M_0)", Z(sp.exp(-mu * L0 / 2) - 1 / (8 * X)))
    # excursion comparison equations
    report("A11 (Z-1/2)' = (X-1/2) - (Z-1/2),  (Y-1/2)' = mu(Z-1/2) - mu(Y-1/2)",
           Z(dZ - ((X - sp.Rational(1, 2)) - (Zc - sp.Rational(1, 2)))) and Z(dY - (mu * (Zc - sp.Rational(1, 2)) - mu * (Y - sp.Rational(1, 2)))))
    # change of variables (eq:change): x = kX, y = kZ, p = lam(1-Y), U = log X, with k = mu lam
    report("A12 p = lam(1-Y): p' = k - mu p - y with y = kZ, k = mu lam", Z(-lam * dY - (mu * lam - mu * lam * (1 - Y) - mu * lam * Zc)))
    # equilibrium linearisation
    fmu = lambda q: q * (q + 1) * (q + 1 - mu)
    kk = g * (g - 1) * (g - mu)          # k = f_mu(gamma-1)
    z = sp.symbols('z')
    chiE = z ** 3 + A * z ** 2 + mu * z + k
    report("A13 chi_E(-gamma) = k - f_mu(gamma-1)", Z(chiE.subs(z, -g) - (k - fmu(g - 1))))
    report("A14 chi_E(z) = (z+gamma_E)(z^2 - (gamma_E-A) z + k/gamma_E) when k = f_mu(gamma_E-1)",
           Z(sp.expand((chiE - (z + g) * (z ** 2 - (g - A) * z + k / g)).subs(k, kk))))
    JE = sp.Matrix([[0, 1, 0], [0, -mu, -1], [kk, 0, -1]])
    etaE, zetaE = g - 1, (g - 1) * (g - mu)
    vE = sp.Matrix([1, -g, -g * (g - mu)])
    nE = sp.Matrix([zetaE, -etaE, -1])
    fp = sp.diff(fmu(s), s).subs(s, etaE)
    report("A15 J_E v_E = -gamma_E v_E,  J_E^T n_E = -gamma_E n_E", Z(JE * vE + g * vE) and Z(JE.T * nE + g * nE))
    report("A16 n_E . v_E = f_mu'(eta_E)", Z(nE.dot(vE) - fp))
    report("A17 f_mu'(s) = (s+1)(s+1-mu) + s(s+1-mu) + s(s+1)",
           Z(sp.diff(fmu(s), s) - ((s + 1) * (s + 1 - mu) + s * (s + 1 - mu) + s * (s + 1))))
    Pu = sp.eye(3) - vE * nE.T / fp
    report("A18 P_u^2 = P_u, P_u v_E = 0, n_E^T P_u = 0, P_u J_E = J_E P_u",
           Z(Pu * Pu - Pu) and Z(Pu * vE) and Z(nE.T * Pu) and Z(Pu * JE - JE * Pu))
    al = (g - A) / 2
    Nm = JE - al * sp.eye(3)
    om2 = kk / g - al ** 2
    report("A19 (J_E - alpha_E I)^2 P_u = -omega_E^2 P_u,  omega_E^2 = k/gamma_E - alpha_E^2", Z(Nm * Nm * Pu + om2 * Pu))
    U, p, y = sp.symbols('U p y', real=True)
    r = sp.Matrix([U, p, y - kk])
    report("A20 Lambda_E(r) = (y-k) - zeta_E U + eta_E p = -n_E . r;  r = P_u r - (Lambda_E/f') v_E",
           Z((y - kk) - zetaE * U + etaE * p + nE.dot(r)) and Z(r - (Pu * r - ((y - kk) - zetaE * U + etaE * p) / fp * vE)))
    F = sp.Matrix([p, kk - mu * p - y, kk * sp.exp(U) - y])
    report("A21 r' = J_E r + (0,0,k Rem(U))", Z(F - JE * r - sp.Matrix([0, 0, kk * (sp.exp(U) - 1 - U)])))
    sg, Up = sp.symbols('sigma U_p', positive=True)
    okrem = True
    for Us_ in (Up, -Up):                # U > 0 and U < 0 separately (at U = 0 both sides vanish)
        okrem = okrem and Z(sp.simplify(Us_ ** 2 * sp.integrate((1 - sg) * sp.exp(sg * Us_), (sg, 0, 1)) - (sp.exp(Us_) - 1 - Us_)))
    report("A22 Rem(U) = U^2 int_0^1 (1-s) e^{sU} ds  (U > 0 and U < 0)", okrem)
    tt = sp.symbols('t', real=True)
    alS, omS = sp.symbols('alpha omega', positive=True)
    # M(t) = e^{-alpha t}[cos(omega t) P - sin(omega t)/omega N P] satisfies M(0) = P and M' = -(N + alpha) M when N^2 P = -omega^2 P
    c_, s_ = sp.cos(omS * tt), sp.sin(omS * tt)
    # coefficients of P and of N P in  M' + (N + alpha) M, using N^2 P = -omega^2 P
    coefP = sp.diff(sp.exp(-alS * tt) * c_, tt) + alS * sp.exp(-alS * tt) * c_ + sp.exp(-alS * tt) * s_ / omS * omS ** 2
    coefNP = sp.diff(-sp.exp(-alS * tt) * s_ / omS, tt) + sp.exp(-alS * tt) * c_ - alS * sp.exp(-alS * tt) * s_ / omS
    report("A25 (cpt:eq:expformula): M(t) = e^{-alpha t}[cos(omega t)P_u - sin(omega t)/omega (J_E-alpha)P_u] solves M' = -J_E M, M(0) = P_u",
           Z(coefP) and Z(coefNP))
    # speed relations
    eta, zeta, x_, k_ = sp.symbols('eta zeta x k_', positive=True)
    pp = k_ - mu * p - y
    zeta_onO = (x_ - y + eta * pp) / p
    report("A23 y-k = -mu p - p';  x-k = (zeta-mu)p - (1+eta)p'  [on zeta p - eta p' = x-y]",
           Z(y - k_ + mu * p + pp) and Z((x_ - k_ - ((zeta - mu) * p - (1 + eta) * pp)).subs(zeta, zeta_onO)))
    report("A24 at p=0 on the constraint: p' = k-y = -(x-k)/(1+eta)",
           Z(((k_ - y) + (x_ - k_) / (1 + eta)).subs(y, sp.solve(sp.Eq(-eta * (k_ - y), x_ - y), y)[0])))


# ----------------------------------------------------------------------------------------------
# Part B: floating point
# ----------------------------------------------------------------------------------------------
def I_of(mu, s):
    return s * np.exp(-s) if abs(mu - 1) < 1e-12 else (np.exp(-s) - np.exp(-mu * s)) / (mu - 1)


def tail_grid():
    worst1 = worst2 = worst3 = -np.inf
    for mu in (0.02, 0.3, 0.999, 1.0, 1.001, 3.0, 40.0):
        mw = min(mu, 1.0)
        s = np.linspace(0, 80 / mw, 200001)
        I = I_of(mu, s)
        th = np.exp(-s) + I
        th_alt = np.exp(-mu * s) + mu * I
        b1 = (1 + mw * s) * np.exp(-mw * s)
        b2 = 2 * np.exp(-mw * s / 2)
        worst1 = max(worst1, np.max(I - s * np.exp(-mw * s)))
        worst2 = max(worst2, np.max(th - b1), np.max(th_alt - b1))
        worst3 = max(worst3, np.max(b1 - b2))
    report('B01 tail grid (7 values of mu): I(s) <= s e^{-mu_w s}; theta_mu <= (1+mu_w s)e^{-mu_w s} <= 2 e^{-mu_w s/2}',
           worst1 < 1e-12 and worst2 < 1e-12 and worst3 < 0, f'largest violations {worst1:.1e}, {worst2:.1e}, {worst3:.1e}')


def cycle_checks(mu, rho, n=20001):
    k = rho * whlib.kH(mu)
    lam = k / mu
    A = 1 + mu
    tag = f'[mu={mu}, rho={rho}]'
    cyc = whlib.find_cycle(mu, k)
    P = whlib.cycle_profile(mu, k, cyc, n=n, settle_periods=6)
    t, U, p, y, x, eta, zeta = (P[q] for q in ('t', 'U', 'p', 'y', 'x', 'eta', 'zeta'))
    T = cyc['T']
    pp = k - mu * p - y
    X, Yc, Zc = np.exp(U), 1 - p / lam, y / k
    avg = lambda f: np.trapezoid(f, t) / T

    # ---- Lemma lem:statebox
    M0 = (lam + mu) * (lam + 1) / mu
    mw = min(mu, 1.0)
    L0 = 2 / mw * np.log(8 * M0)
    log_beta0 = np.log(0.5) - lam * (M0 - 1) * L0
    R, S = Zc / X, Yc / X
    report(f'B02 {tag} ratio bounds: min Z/X > 1/(lam+1), min Y/X > 1/M_0',
           R.min() > 1 / (lam + 1) and S.min() > 1 / M0,
           f'min Z/X = {R.min():.4f} vs {1/(lam+1):.4f};  min Y/X = {S.min():.4f} vs {1/M0:.4f}')
    report(f'B03 {tag} upper box: max(X,Y,Z) < M_0, max Y <= max Z <= max X, min Y >= min Z >= min X, Y = 1 at max X',
           max(X.max(), Yc.max(), Zc.max()) < M0 and Yc.max() <= Zc.max() <= X.max() and Yc.min() >= Zc.min() >= X.min()
           and abs(Yc[int(np.argmax(X))] - 1) < 1e-6,
           f'max X = {X.max():.4f}, max Z = {Zc.max():.4f}, max Y = {Yc.max():.4f}, M_0 = {M0:.4f}')
    report(f'B04 {tag} means <X> = <Y> = <Z> = 1', max(abs(avg(X) - 1), abs(avg(Yc) - 1), abs(avg(Zc) - 1)) < 1e-7,
           f'{avg(X):.9f}, {avg(Yc):.9f}, {avg(Zc):.9f}')
    if X.min() < 0.5:
        im = int(np.argmin(X))
        below = X < 0.5
        i0 = im
        while below[i0 - 1]:
            i0 -= 1
        i1 = im
        while below[i1 + 1]:
            i1 += 1
        s = t[i0:i1 + 1] - t[i0]
        Zb = (Zc[i0] - 0.5) * np.exp(-s)
        Yb1 = (Yc[i0] - 0.5) * np.exp(-mu * s) + mu * (M0 - 0.5) * I_of(mu, s)
        Yb2 = (M0 - 0.5) * (np.exp(-s) + I_of(mu, s))
        okZ = np.all(Zc[i0:i1 + 1] - 0.5 <= Zb + 1e-9)
        okY = np.all(Yc[i0:i1 + 1] - 0.5 <= Yb1 + 1e-9) and np.all(Yb1 <= Yb2 + 1e-9)
        late = s >= L0
        okL = (not late.any()) or (np.all(Yc[i0:i1 + 1][late] < 0.75) and np.all((lam * X * (1 - Yc))[i0:i1 + 1][late] > 0))
        report(f'B05 {tag} excursion below 1/2: Z-1/2 <= (Z(t0)-1/2)e^(-s); Y-1/2 <= variation-of-constants bound <= (M_0-1/2) theta_mu(s); '
               f'Y < 3/4 and X\' > 0 for s >= L_0', bool(okZ and okY and okL),
               f'excursion length {s[-1]:.3f}, t_m - t_0 = {t[im]-t[i0]:.3f} < L_0 = {L0:.3f}')
        report(f'B06 {tag} t_m - t_0 < L_0 and log min X >= log beta_0', t[im] - t[i0] < L0 and np.log(X.min()) >= log_beta0,
               f'log min X = {np.log(X.min()):.3f}, log beta_0 = {log_beta0:.1f};  (log X)\' >= -lam(M_0-1): '
               f'{(lam*(1-Yc)).min():.3f} >= {-lam*(M0-1):.3f}')
    else:
        report(f'B05 {tag} min X >= 1/2 (no excursion)', True, f'min X = {X.min():.4f}')
    report(f'B07 {tag} log box (cpt:eq:logbox): U in [log beta_0, log M_0], p in [-lam(M_0-1), lam(1-beta_0)]',
           log_beta0 <= U.min() and U.max() <= np.log(M0) and -lam * (M0 - 1) <= p.min() and p.max() <= lam)

    # ---- Lemma lem:eqseparation
    ed = whlib.equilibrium_data(mu, k)
    g = ed['gammaE']
    etaE, zetaE = g - 1, (g - 1) * (g - mu)
    JE = np.array([[0, 1, 0], [0, -mu, -1], [k, 0, -1.0]])
    vE = np.array([1.0, -g, -g * (g - mu)])
    nE = np.array([zetaE, -etaE, -1.0])
    fp = (etaE + 1) * (etaE + 1 - mu) + etaE * (etaE + 1 - mu) + etaE * (etaE + 1)
    Pu = np.eye(3) - np.outer(vE, nE) / fp
    al = 0.5 * (g - A)
    om = np.sqrt(k / g - al ** 2)
    ev = np.linalg.eigvals(JE)
    report(f'B08 {tag} spectrum of J_E: -gamma_E and alpha_E +- i omega_E; gamma_E > A, gamma_E - mu > 1',
           np.min(np.abs(ev + g)) < 1e-9 and np.min(np.abs(ev - (al + 1j * om))) < 1e-9 and g > A and g - mu > 1,
           f'gamma_E = {g:.6f}, alpha_E = {al:.6f}, omega_E = {om:.6f}')
    Nm = JE - al * np.eye(3)
    # expm(-J t) grows like e^{gamma_E t} along v_E, which amplifies the roundoff of P_u v_E = 0; the comparison with the
    # matrix exponential is therefore made at small t only (the formula itself is checked exactly in A19, A25)
    errexp = 0.0
    for tt in (0.05, 0.2, 0.5, 1.0, 2.0):
        lhs = expm(-JE * tt) @ Pu
        rhs = np.exp(-al * tt) * (np.cos(om * tt) * Pu - np.sin(om * tt) / om * (Nm @ Pu))
        errexp = max(errexp, np.max(np.abs(lhs - rhs)) / max(1e-300, np.max(np.abs(rhs))))
    nrm = lambda Mx: np.linalg.norm(Mx, 2)
    Theta = nrm(Pu) + nrm(Nm @ Pu) / om
    Mform = lambda tt: np.exp(-al * tt) * (np.cos(om * tt) * Pu - np.sin(om * tt) / om * (Nm @ Pu))
    worst_decay = max(nrm(Mform(tt)) * np.exp(al * tt) / Theta for tt in np.linspace(0, 40, 801))
    report(f'B09 {tag} (cpt:eq:expformula) against expm for t <= 2, and ||e^(-J t) P_u|| <= Theta_E e^(-alpha t) (formula, t <= 40)',
           errexp < 1e-9 and worst_decay <= 1 + 1e-12, f'relative error {errexp:.1e};  max of ||e^(-Jt)P_u|| e^(alpha t)/Theta_E = {worst_decay:.4f}')
    G = solve_continuous_lyapunov(JE.T, Pu.T @ Pu)                 # J^T G + G J = P_u^T P_u
    Gq, _ = quad_vec(lambda tt: (lambda Mx: Mx.T @ Mx)(np.exp(-al * tt) * (np.cos(om * tt) * Pu - np.sin(om * tt) / om * (Nm @ Pu))),
                     0, 60.0 / al, epsabs=1e-12, epsrel=1e-11, limit=2000)
    nu = Theta ** 2 / (2 * al)
    rng = np.random.default_rng(8)
    okV = True
    worstN = -np.inf
    for _ in range(2000):
        rr = rng.normal(size=3)
        rr *= rng.uniform(0.01, 1.0) / np.linalg.norm(rr)
        w = Pu @ rr
        Vv = rr @ G @ rr
        okV = okV and (w @ w / (2 * nrm(JE)) <= Vv * (1 + 1e-9) + 1e-15) and (Vv <= nu * (w @ w) * (1 + 1e-9) + 1e-15)
        okV = okV and abs(2 * rr @ G @ (JE @ rr) - w @ w) < 1e-9 * max(1.0, w @ w)
        hh = rng.normal(size=3)
        okV = okV and abs(2 * rr @ G @ hh) <= 2 * nu * np.linalg.norm(w) * np.linalg.norm(hh) * (1 + 1e-9)
        Nv = np.array([0, 0, k * (np.exp(rr[0]) - 1 - rr[0])])
        worstN = max(worstN, abs(2 * rr @ G @ Nv) - np.e * k * nu * np.linalg.norm(w) * (rr @ rr))
    report(f'B10 {tag} V = r^T G r: Lyapunov solution equals the integral (cpt:eq:V); bounds (cpt:eq:Vbounds); DV(r)J_E r = |P_u r|^2; '
           f'|DV(r)N(r)| <= e k nu |P_u r| |r|^2 for |r| <= 1 (2000 random r)',
           np.max(np.abs(G - Gq)) < 1e-7 * np.max(np.abs(G)) and okV and worstN <= 1e-12,
           f'max|G - quadrature|/max|G| = {np.max(np.abs(G-Gq))/np.max(np.abs(G)):.1e};  nu_E = {nu:.3f}')
    rc = np.vstack([U, p, y - k]).T
    w = rc @ Pu.T
    LamE = (y - k) - zetaE * U + etaE * p
    cr = LamE / fp
    nr, nw = np.linalg.norm(rc, axis=1), np.linalg.norm(w, axis=1)
    inCp = (U >= 0) & (p <= 0) & (y - k <= 0)
    inCm = (U <= 0) & (p >= 0) & (y - k >= 0)
    C1 = 1 + np.linalg.norm(vE)
    report(f'B11 {tag} Lambda_E > 0, no cycle point in C_+ or C_-, c_r < |P_u r|, |r| <= (1+|v_E|)|P_u r|',
           LamE.min() > 0 and not inCp.any() and not inCm.any() and np.all(cr < nw) and np.all(nr <= C1 * nw),
           f'max |r|/|P_u r| = {(nr/nw).max():.4f} vs 1+|v_E| = {C1:.4f};  max c_r/|P_u r| = {(cr/nw).max():.4f}')
    Vc = np.einsum('ij,jk,ik->i', rc, G, rc)
    rdot = np.vstack([p, pp, x - y]).T
    DV = 2 * np.einsum('ij,jk,ik->i', rc, G, rdot)
    Nc = np.zeros_like(rc)
    Nc[:, 2] = k * (np.exp(U) - 1 - U)
    DV2 = nw ** 2 + 2 * np.einsum('ij,jk,ik->i', rc, G, Nc)
    fd = np.gradient(Vc, t, edge_order=2)
    eps = min(1 / C1, 1 / (2 * np.e * k * nu * C1 ** 2))
    dsep = eps / (nrm(Pu) * np.sqrt(2 * nrm(JE) * nu))
    i0 = int(np.argmin(Vc))
    small = nw <= eps
    okii = (not small.any()) or np.all(DV[small] >= Vc[small] / (2 * nu))
    report(f'B12 {tag} (cpt:eq:Vdot): DV along the flow = |P_u r|^2 + DV(r)N(r); separation with Q = one point: '
           f'|P_u r| > eps at the minimum of V, |r| > d at every point',
           np.max(np.abs(DV - DV2)) < 1e-8 * np.max(np.abs(DV)) and np.max(np.abs(fd - DV)[5:-5]) < 1e-4 * np.max(np.abs(DV))
           and nw[i0] > eps and nr.min() > dsep and okii,
           f'eps = {eps:.3e}, |P_u r| at min V = {nw[i0]:.3e};  d = {dsep:.3e}, min |r| = {nr.min():.3e};  '
           f'cycle points with |P_u r| <= eps: {int(small.sum())}')

    # ---- Lemma lem:speed
    res1 = np.max(np.abs((y - k) + mu * p + pp))
    res2 = np.max(np.abs((x - k) - ((zeta - mu) * p - (1 + eta) * pp))) / np.max(np.abs(x - k))
    xm, xM = x.min(), x.max()
    eP = 0.5 * (mu - 1 + np.sqrt((mu - 1) ** 2 + 4 * xM))
    zm = xm / (1 + eP)
    em = 0.5 * (mu - 1 + np.sqrt((mu - 1) ** 2 + 4 * zm))
    zP = xM / (1 + em)
    okpb = em <= eta.min() and eta.max() <= eP and zm <= zeta.min() and zeta.max() <= zP
    speed = np.hypot(p, pp)
    # the lemma's constant with x_- = k beta_0, x_+ = k M_0 (in logarithms: beta_0 underflows)
    eP0 = 0.5 * (mu - 1 + np.sqrt((mu - 1) ** 2 + 4 * k * M0))
    zP0 = k * M0                                    # zeta_+ = x_+/(1+eta_-) <= x_+
    log10_C2 = np.log10(max(zP0, mu) ** 2 + (1 + eP0) ** 2) - 2 * (np.log(k) + log_beta0) / np.log(10)
    # same chain with the cycle's own min X in place of beta_0 (a non-vacuous test of the three estimates)
    b_act = X.min()
    C2act = 2 + mu ** 2 + (max(zP, mu) ** 2 + (1 + eP) ** 2) / (k * min(b_act, 1.0)) ** 2
    okchain = np.all((y - k) ** 2 <= (mu ** 2 + 1) * speed ** 2 * (1 + 1e-12)) and \
        np.all((x - k) ** 2 <= (max(zP, mu) ** 2 + (1 + eP) ** 2) * speed ** 2 * (1 + 1e-9)) and \
        np.all(np.abs(U) <= np.abs(x - k) / (k * min(b_act, 1.0)) * (1 + 1e-9) + 1e-12) and np.all(nr ** 2 <= C2act * speed ** 2 * (1 + 1e-9))
    report(f'B13 {tag} speed relations (eq:speedrelations); Lemma lem:projectivebounds with the cycle\'s own x-range; '
           f'the three estimates of the proof with min X in place of beta_0',
           res1 < 1e-9 and res2 < 1e-7 and okpb and bool(okchain),
           f'residuals {res1:.1e}, {res2:.1e};  max |r|/speed = {(nr/speed).max():.4f};  min speed = {speed.min():.4e};  '
           f'log10 of the lemma\'s C_Q\'^2 (one-point Q) ~ {log10_C2:.0f}')

    # ---- Theorem thm:compact
    arc = np.trapezoid(speed, t)
    tv = 2 * (U.max() - U.min()) + 2 * (p.max() - p.min())
    tv_exact = np.sum(np.abs(np.diff(U))) + np.sum(np.abs(np.diff(p)))
    log_bound = np.log(2 * (np.log(M0) - log_beta0) + 2 * lam * M0)
    ndown = int(np.sum((p[:-1] > 0) & (p[1:] <= 0)))
    nup = int(np.sum((p[:-1] < 0) & (p[1:] >= 0)))
    ncrit_p = int(np.sum(pp[:-1] * pp[1:] < 0))
    report(f'B14 {tag} arclength <= int(|U\'|+|p\'|) = 2(max U - min U) + 2(max p - min p) <= box bound; T > 2 pi/sqrt(mu); '
           f'one downward and one upward zero of p, two critical points of p per period',
           arc <= tv * (1 + 1e-9) and abs(tv - tv_exact) < 1e-6 * tv and np.log(tv) <= log_bound and T > 2 * np.pi / np.sqrt(mu)
           and ndown + (1 if p[0] == 0 else 0) <= 1 and nup == 1 and ncrit_p == 2,
           f'arclength = {arc:.4f}, total variation = {tv:.4f}, T = {T:.4f} > {2*np.pi/np.sqrt(mu):.4f};  '
           f'T sigma <= arclength with the actual min speed: {speed.min()*T:.4f} <= {arc:.4f}')


if __name__ == '__main__':
    part_a()
    tail_grid()
    for mu_, rho_ in ((0.3, 1.5), (1.0, 2.0), (3.0, 1.5), (0.1, 1.1)):
        cycle_checks(mu_, rho_)
    print(f"\n{'ALL PASS' if FAIL == 0 else str(FAIL) + ' FAILURES'}")
    sys.exit(1 if FAIL else 0)
