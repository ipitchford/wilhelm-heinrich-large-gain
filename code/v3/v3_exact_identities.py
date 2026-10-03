#!/usr/bin/env python3
"""Exact symbolic check of every new algebraic identity used in the paper (sympy).

Only residual-equals-zero tests are counted.  Each test expands a polynomial (or, after clearing
explicit denominators, a polynomial numerator) and compares it with 0.

Conventions.  The flow is  U' = p, p' = k - mu p - y, y' = x - y, x = k e^U, so x' = x p; x is kept as
an independent symbol, which turns every identity below into a polynomial one.  The adjoint ratios
obey  eta' = zeta + (mu-1) eta - eta^2,  zeta' = x - (1+eta) zeta.  The tangent-orthogonality defect
is  O = zeta p - eta p' - (x - y); on a cycle, and on the attractor, O = 0.  An identity that holds
only on O = 0 is tested in its off-constraint form  residual = cofactor * O, with the cofactor found
by exact polynomial division and written into the receipt.

--inject-error adds a gauge identity with a wrong sign; the run must then fail (used by v3_controls).
This is finite algebra.  Positivity statements, averaging arguments and topology are in the written proofs.
"""
from __future__ import annotations
import argparse, sys
from pathlib import Path
import sympy as sp
from common import Clock, write_receipt, EVID


def run(inject_error=False, out=None, name='v3_exact_identities'):
    clock = Clock()
    U, p, y, x, eta, zeta, mu, k = sp.symbols('U p y x eta zeta mu k')
    a, b, c = sp.symbols('a b c')                 # strong-stable Floquet components, e_s = (a,-b,-c)
    U2, p2, y2, x2 = sp.symbols('U2 p2 y2 x2')    # a second orbit
    A = 1 + mu
    lam = k / mu
    pd = k - mu * p - y
    vf = {U: p, p: pd, y: x - y, x: x * p,
          eta: zeta + (mu - 1) * eta - eta**2, zeta: x - (1 + eta) * zeta,
          a: -b, b: -mu * b - c, c: -x * a - c,
          U2: p2, p2: k - mu * p2 - y2, y2: x2 - y2, x2: x2 * p2}

    def D(f, n=1):
        for _ in range(n):
            f = sum(sp.diff(f, v) * r for v, r in vf.items())
        return f

    rows = []

    def zero(group, ident, expr, note=''):
        """exact test: numerator of expr, expanded, is the zero polynomial"""
        num = sp.numer(sp.together(expr))
        r = sp.expand(num)
        rows.append(dict(group=group, id=ident, passed=bool(r == 0), residual=str(r)[:200], note=note))

    gens = (U, p, y, x, eta, zeta)
    hd = y - k
    hdd = x - y
    pdd = D(pd)
    O = zeta * p - eta * pd - hdd
    Opoly = sp.Poly(sp.expand(O), *gens)

    def modO(group, ident, expr, expected_cofactor=None, note=''):
        """exact test: expr - q*O == 0 with q the polynomial quotient; optional check of q itself"""
        num = sp.expand(expr)
        q, rem = sp.div(sp.Poly(num, *gens), Opoly)
        qe = q.as_expr()
        r = sp.expand(num - qe * O)
        ok = (r == 0)
        if expected_cofactor is not None:
            ok = ok and sp.expand(qe - expected_cofactor) == 0
        rows.append(dict(group=group, id=ident, passed=bool(ok), residual=str(r)[:200],
                         cofactor_of_O=str(sp.factor(qe)), note=note))

    K = zeta * U - eta * p - hd
    W = D(zeta) * p - D(eta) * pd
    Sg = p**2 - U * pd                    # star gauge
    Xg = pd**2 - p * pdd                  # convexity gauge
    V3 = U * hdd - p * hd                 # middle component of the radial wedge
    kpsi = x * U - x + k                  # k psi(U)
    kR = x - k - k * U                    # k R(U), R(d) = e^d - 1 - d

    # ------------------------------------------------------------------ inputs inherited from the framework
    g = 'framework identities re-derived'
    zero(g, 'O_transport', D(O) + (1 + eta) * O)
    zero(g, 'K_equation', D(K) + (1 + eta) * K - kpsi)
    modO(g, 'W_equation', D(W) + (1 + eta) * W - x * p**2)
    modO(g, 'tangent_second_order', pdd - ((eta - mu) * pd - zeta * p))
    zero(g, 'scalar_equation', pdd + A * pd + mu * p + x - k)

    # ------------------------------------------------------------------ gauges
    g = 'gauge ODEs and exponent formulas'
    modO(g, 'star_gauge_ode', D(Sg) - ((eta - mu) * Sg + p * K), expected_cofactor=-U,
         note="Sigma' = (eta-mu) Sigma + p K - U O")
    modO(g, 'convexity_gauge_ode', D(Xg) - ((eta - mu) * Xg + p * W),
         note="Xi' = (eta-mu) Xi + p W + cofactor * O")
    zero(g, 'star_gauge_ode_state_form', D(Sg) - (-mu * Sg + V3))
    zero(g, 'log_derivative_star_gauge', D(Sg) - (V3 / Sg - mu) * Sg, note="(log Sigma)' = eta_0 - mu, eta_0 = V3/Sigma")
    modO(g, 'eta_minus_eta0', eta * Sg - V3 + p * K, note='eta - eta_0 = -p K / Sigma on O = 0')
    modO(g, 'zeta_from_gauges', zeta * Sg - (Xg - pd * K), note='zeta Sigma = Xi - p\' K on O = 0')
    zero(g, 'convexity_gauge_state_form', Xg - (pd * (k - y) + p * (x - y)))
    zero(g, 'convexity_gauge_ode_state_form', D(Xg) - (x * Sg - A * Xg + kpsi * pd))
    zero(g, 'star_gauge_at_p0', Sg.subs(p, 0) - U * (y - k))
    zero(g, 'convexity_gauge_at_p0', Xg.subs(p, 0) - (k - y)**2)
    cs = sp.Symbol('c_shift')
    sc = p**2 - pd * (U + cs)
    modO(g, 'shifted_centre_gauge_ode', D(sc) - ((eta - mu) * sc + p * (K + cs * zeta)))
    zero(g, 'shifted_centre_gauge_at_p0', sc.subs(p, 0) - (y - k) * (U + cs))
    # generalised Wronskian principle
    L2 = lambda h: D(h, 2) + (mu - eta) * D(h) + zeta * h
    zero(g, 'L2_of_1', L2(sp.Integer(1)) - zeta)
    zero(g, 'L2_of_U', L2(U) - K)
    modO(g, 'L2_of_p', L2(p))
    modO(g, 'L2_of_y', L2(hd) - W)
    modO(g, 'W_state_form', W - (x * p - (1 + eta - mu) * hdd + zeta * hd))

    # ------------------------------------------------------------------ kappa
    g = 'scalar Riccati reduction'
    kap = K / Sg
    eta0 = V3 / Sg
    zeta0 = Xg / Sg
    expr = sp.together(D(kap) - (kpsi / Sg - (1 + eta + eta0 - mu) * kap))
    modO(g, 'kappa_equation', sp.expand(sp.numer(expr) * sp.cancel(Sg**2 / sp.denom(expr))),
         note="Sigma^2 {kappa' - [k psi/Sigma - (1+eta+eta_0-mu) kappa]} = cofactor * O")
    kp = sp.Symbol('kappa')
    # closed scalar Riccati: substitute eta = eta_0 - p kappa in the K-equation divided by Sigma
    Dkap_closed = (kpsi - (1 + eta0 - p * kp) * kp * Sg) / Sg - kp * (V3 - mu * Sg) / Sg
    zero(g, 'kappa_riccati_closed', Dkap_closed - (kpsi / Sg - (1 - mu + 2 * eta0) * kp + p * kp**2))
    Ks = sp.Symbol('K_s')
    zero(g, 'K_riccati_closed', (kpsi - (1 + (V3 - p * Ks) / Sg) * Ks) - (kpsi - (1 + V3 / Sg) * Ks + (p / Sg) * Ks**2))
    zero(g, 'W_affine_in_kappa', (x * p - (1 + eta - mu) * hdd + zeta * hd).subs({eta: eta0 - p * kp, zeta: zeta0 - pd * kp})
         - ((x * p - (1 + eta0 - mu) * hdd + zeta0 * hd) + kp * Xg))
    zero(g, 'kappa0_plane_contains_tangent', Sg * (zeta0 * p - eta0 * pd - hdd))
    zero(g, 'kappa0_plane_contains_E', Sg * (zeta0 * U - eta0 * p - hd))
    zero(g, 'eta0_riccati_defect', D(eta0) - (zeta0 + (mu - 1) * eta0 - eta0**2) - kpsi * p / Sg)
    zero(g, 'zeta0_riccati_defect', D(zeta0) - (x - (1 + eta0) * zeta0) - kpsi * pd / Sg)

    # ------------------------------------------------------------------ equilibrium quadratic form, Lambda_E
    g = 'equilibrium form and Lambda_E'
    gE = sp.Symbol('gamma_E')
    kE = gE * (gE - 1) * (gE - mu)
    etaE = gE - 1
    zetaE = etaE * (etaE + 1 - mu)
    al = (etaE - mu) / 2
    QE = p**2 - 2 * al * U * p + zetaE * U**2
    LamE = hd - zetaE * U + etaE * p
    zero(g, 'LambdaE_ode', (D(LamE) + gE * LamE - kR).subs(k, kE))
    zero(g, 'star_gauge_decomposition', Sg - QE - U * LamE)
    zero(g, 'QE_derivative', (D(QE) - 2 * al * QE + 2 * (p - al * U) * LamE).subs(k, kE))
    zero(g, 'sturm_form', pd - (etaE - mu) * p + zetaE * U + LamE, note="U'' - 2 sigma U' + zeta_E U = -Lambda_E")
    zz = sp.Symbol('z')
    zero(g, 'characteristic_factorisation', zz**3 + A * zz**2 + mu * zz + kE - (zz + gE) * (zz**2 + (mu - etaE) * zz + zetaE))
    zero(g, 'frequency_formula', 4 * zetaE - (gE - A)**2 - (3 * gE**2 - 2 * A * gE - (1 - mu)**2))
    zero(g, 'riccati_fixed_point_eta', (zeta + (mu - 1) * eta - eta**2).subs({zeta: zetaE, eta: etaE}))
    zero(g, 'riccati_fixed_point_zeta', (x - (1 + eta) * zeta).subs({x: kE, zeta: zetaE, eta: etaE}))

    # ------------------------------------------------------------------ polar form
    g = 'polar form and nested radii'
    th, r_, m, m1, m2 = sp.symbols('theta r m m1 m2')
    Up, pp_ = r_ * sp.cos(th), r_ * sp.sin(th)
    ppd = -mu * pp_ - m * r_
    rdot = (Up * pp_ + pp_ * ppd) / r_
    thdot = (Up * ppd - pp_ * pp_) / r_**2
    sig = sp.sin(th)**2 + mu * sp.sin(th) * sp.cos(th) + m * sp.cos(th)
    Xi_ = sp.sin(th) * (m + mu * sp.sin(th) - sp.cos(th)) / sig

    def tzero(ident, e):
        r = sp.simplify(sp.expand_trig(sp.together(e)))
        rows.append(dict(group=g, id=ident, passed=bool(r == 0), residual=str(r)[:200], note=''))
    tzero('polar_angle_rate', thdot + sig)
    tzero('polar_radial_rate', rdot / r_ / thdot - Xi_)
    tzero('polar_linear_fractional_difference',
          Xi_.subs(m, m2) - Xi_.subs(m, m1) - sp.sin(th) * (m2 - m1) / (sig.subs(m, m1) * sig.subs(m, m2)))

    # ------------------------------------------------------------------ Lotka-Volterra chain, second-order equations
    g = 'Lotka-Volterra chain and second-order equations'
    q = lam - p
    r1, r2 = x / y, y / q
    zero(g, 'LV_r1', D(r1) - r1 * (lam + 1 - r1 - q))
    zero(g, 'LV_r2', D(r2) - r2 * (mu - 1 + r1 - r2))
    zero(g, 'LV_q', D(q) - q * (r2 - mu))
    f = lambda s: s * (s + 1) * (s + 1 - mu)
    zero(g, 'eta_second_order', D(eta, 2) + (3 * eta + 2 - mu) * D(eta) + f(eta) - x)
    gg = 1 + eta
    zero(g, 'gamma_second_order', D(gg, 2) + (3 * gg - A) * D(gg) + gg * (gg - 1) * (gg - mu) - x)
    zero(g, 'chemical_ratio_second_order', D(r2, 2) + (3 * r2 + 2 - mu) * D(r2) + f(r2) - x * (A + p) / q)
    sigL = zeta + eta + eta**2 / 2
    GL = sp.Rational(3, 2) * eta**2 + (2 - mu) * eta
    zero(g, 'lienard_plane_eta', D(eta) - (sigL - GL))
    zero(g, 'lienard_plane_sigma', D(sigL) - (x - f(eta)))
    zero(g, 'log_eta_derivative', D(eta) / eta - (zeta / eta + mu - 1 - eta))
    zero(g, 'log_zeta_derivative', D(zeta) / zeta - (x / zeta - 1 - eta))
    qq = sp.Symbol('q_dev')
    zero(g, 'moment_polynomial', f(mu + qq) - (mu * A + (A**2 + mu) * qq + 2 * A * qq**2 + qq**3))
    Lm = lambda ph: -D(ph, 3) + 2 * A * D(ph, 2) - (A**2 + mu) * D(ph) + (mu * A - x) * ph
    Lc = lambda ph: D(ph, 3) + 2 * A * D(ph, 2) + (A**2 + mu) * D(ph) + (mu * A - x) * ph
    zero(g, 'backward_operator_on_x', Lm(x) - x * ((A - p) * (3 * pd - p * (A - p)) - (k - mu * A)))
    zero(g, 'compound_operator_on_q', Lc(q) - x * (2 * p + A - lam))
    zero(g, 'compound_operator_on_inverse_x', Lc(1 / x) * x - (-(A - 3 * p) * pd - p * (A - p)**2 - (k - mu * A)))
    zero(g, 'log_1_plus_eta', D(sp.log(1 + eta)) - ((zeta - mu) / (1 + eta) + mu - eta))
    zero(g, 'log_zeta_plus_mu_eta', D(sp.log(zeta + mu * eta)) - (x / (zeta + mu * eta) - (1 + eta - mu)))

    # ------------------------------------------------------------------ strong-stable vector and weak adjoint
    g = 'strong-stable ratios and weak adjoint'
    J = sp.Matrix([[0, 1, 0], [0, -mu, -1], [x, 0, -1]])
    es = sp.Matrix([a, -b, -c])
    for i in range(3):
        zero(g, f'strong_stable_vector_row{i+1}', D(es[i]) - (J * es)[i])
    r, s = b / a, c / a
    zero(g, 'ratio_r_equation', D(r) - (r**2 - mu * r - s))
    zero(g, 'ratio_s_equation', D(s) - (s * r - s - x))
    zero(g, 'ratio_r_second_order', D(r, 2) + (A - 3 * r) * D(r) + r * (r - 1) * (r - mu) - x)
    ct = c + (mu - 1) * b
    zero(g, 'conjugate_c_definition', ct + D(b) + b)
    zero(g, 'conjugate_c_equation', D(ct) - (-x * a - mu * ct))
    zero(g, 'ratio_r_conjugate_form', D(r) - (r**2 - r - ct / a))
    zero(g, 'third_order_operator_on_a', D(a, 3) + A * D(a, 2) + mu * D(a) + x * a)
    be, ka = sp.symbols('beta kappa_s')
    bed = be**2 - mu * be - ka
    kad = ka * (be - 1) - x
    zero(g, 'strong_stable_gain_primitive', kad + be * bed - bed - (be * (be - 1) * (be - mu) - x))
    zp = sp.Matrix([p, pd, hdd])
    Psi = zp.cross(es)                              # Psi e^{-At}
    for i in range(3):
        zero(g, f'weak_adjoint_equation_row{i+1}', D(Psi[i]) + A * Psi[i] + (J.T * Psi)[i],
             note="Psi = e^{At} z' x e_s solves Psi' = -J^T Psi")
    formula = a * sp.Matrix([r * hdd - s * pd, hdd + s * p, -(pd + r * p)])
    for i in range(3):
        zero(g, f'weak_adjoint_components_row{i+1}', Psi[i] - formula[i])
    zero(g, 'weak_adjoint_orthogonal_to_tangent', Psi.dot(zp))
    zero(g, 'weak_adjoint_orthogonal_to_strong_stable', Psi.dot(es))
    d = sp.Matrix([U2 - U, p2 - p, y2 - y])
    Bc = Psi.dot(d)
    zero(g, 'chord_identity', D(Bc) + A * Bc - Psi[2] * (x2 - x - x * (U2 - U)),
         note="(Psi.d)' = Psi_3 x R(dU), with x e^{dU} written as x2")
    zpp = sp.Matrix([D(v) for v in zp])
    Cz = Psi.dot(zpp)
    zero(g, 'curvature_pairing', D(Cz) + A * Cz - Psi[2] * x * p**2)
    theta = -pd / p
    zero(g, 'theta_equation', D(theta) - (theta**2 - mu * theta + hdd / p))
    zero(g, 'crossing_derivative_Psi3', D(theta - r) - ((theta - r) * (theta + r - mu) + hdd / p + s))
    psi_p = pd + r * p
    g2 = hdd / p + s
    zero(g, 'crossing_derivative_Psi2', D(g2) - s * psi_p / p + g2 * (1 + pd / p))
    g1 = hdd / pd - s / r
    zero(g, 'crossing_derivative_Psi1', D(g1) - x * psi_p / (r * pd) - g1 * (mu - 1 + hdd / pd + s / r))
    zero(g, 'supporting_line_curvature', pd**2 - p * pdd - (pd**2 + mu * p * pd + p * hdd))
    zero(g, 'supporting_line_curvature_at_contact',
         (pd**2 + mu * p * pd + p * hdd).subs(y, k - mu * p + r * p) - (p**2 * r * (r - mu) + p * (x - (k - mu * p + r * p))))
    zero(g, 'tangency_second_derivative', pdd + r * pd - ((r - mu) * pd - hdd))
    # rank-one identity: omega_phi' - J omega_phi = -L3[phi] e_3, and L3 on the standard test functions
    ph0, ph1, ph2, ph3 = sp.symbols('phi0 phi1 phi2 phi3')     # phi and its first three derivatives
    om = sp.Matrix([ph0, ph1, -ph2 - mu * ph1])
    omd = sp.Matrix([ph1, ph2, -ph3 - mu * ph2])
    L3phi = ph3 + A * ph2 + mu * ph1 + x * ph0
    res = omd - J * om + sp.Matrix([0, 0, L3phi])
    for i in range(3):
        zero(g, f'rank_one_carrier_row{i+1}', res[i])
    L3 = lambda h: D(h, 3) + A * D(h, 2) + mu * D(h) + x * h
    zero(g, 'L3_of_1', L3(sp.Integer(1)) - x)
    zero(g, 'L3_of_p', L3(p))
    zero(g, 'L3_of_pprime', L3(pd) + x * p**2)
    zero(g, 'L3_of_U', L3(U) - kpsi)
    zero(g, 'L3_of_y', L3(y) - x * (k + p**2))
    Msym = sp.Matrix(3, 3, sp.symbols('m1:10'))
    u_, v_ = sp.Matrix(sp.symbols('u1:4')), sp.Matrix(sp.symbols('v1:4'))
    cp_res = (Msym * u_).cross(v_) + u_.cross(Msym * v_) - (Msym.trace() * sp.eye(3) - Msym.T) * u_.cross(v_)
    for i in range(3):
        zero(g, f'cross_product_transport_row{i+1}', cp_res[i])

    # ------------------------------------------------------------------ compound system
    g = 'compound-system rows'
    Cm = sp.Matrix([[-mu, 1, 0], [0, -1, 1], [x, 0, -A]])
    Sd = sp.diag(1, -1, -1)
    J2 = sp.Matrix([[J[0, 0] + J[1, 1], J[1, 2], -J[0, 2]], [J[2, 1], J[0, 0] + J[2, 2], J[0, 1]],
                    [-J[2, 0], J[1, 0], J[1, 1] + J[2, 2]]])
    res = Sd * J2 * Sd - Cm
    zero(g, 'compound_matrix_conjugacy', sum(sp.expand(e)**2 for e in res))
    omR = sp.Matrix([Sg, V3, Xg])
    resR = sp.Matrix([D(v) for v in omR]) - Cm * omR - kpsi * sp.Matrix([0, p, pd])
    for i in range(3):
        zero(g, f'radial_wedge_row{i+1}', resR[i], note="(Sigma, V3, Xi)' = C (Sigma, V3, Xi) + k psi (0, p, p')")
    wedge = lambda v, u: sp.Matrix([v[1] * u[0] - v[0] * u[1], v[0] * u[2] - v[2] * u[0], v[1] * u[2] - v[2] * u[1]])
    resW = wedge(J * v_, u_) + wedge(v_, J * u_) - Cm * wedge(v_, u_)
    for i in range(3):
        zero(g, f'compound_transport_row{i+1}', resW[i])
    rvec = sp.Matrix([U, p, hd])
    wr = sp.Matrix([-(U * pd - p * p), U * hdd - hd * p, p * hdd - hd * pd])
    for i in range(3):
        zero(g, f'radial_wedge_components_row{i+1}', wr[i] - omR[i])
    P_ = lam - p
    omY = sp.Matrix([P_, y, x])
    resY = sp.Matrix([D(v) for v in omY]) - Cm * omY - sp.Matrix([0, 0, x * (2 * p + A - lam)])
    for i in range(3):
        zero(g, f'chemical_direction_defect_row{i+1}', resY[i])
    omc = sp.Matrix([1, mu, mu])
    resc = -Cm * omc - sp.Matrix([0, 0, mu * A - x])
    for i in range(3):
        zero(g, f'constant_direction_defect_row{i+1}', resc[i])
    nu = sp.Matrix([c, b, a])
    nud = sp.Matrix([A * c + D(c), A * b + D(b), A * a + D(a)])
    resn = nud + Cm.T * nu
    for i in range(3):
        zero(g, f'compound_adjoint_row{i+1}', resn[i])
    # supersolution defect rows for omega = phi (1, e, z)
    e_, z_, phi, de, dz, dphi = sp.symbols('e z phi de dz dphi')
    om = sp.Matrix([phi, phi * e_, phi * z_])
    omd = sp.Matrix([dphi, dphi * e_ + phi * de, dphi * z_ + phi * dz])
    gi = [dphi / phi + mu - e_, dphi / phi + de / e_ + 1 - z_ / e_, dphi / phi + dz / z_ + A - x / z_]
    dfc = omd - Cm * om
    for i in range(3):
        zero(g, f'supersolution_defect_row{i+1}', dfc[i] - gi[i] * om[i])
    gt = [mu - eta, D(eta) / eta + 1 - zeta / eta, D(zeta) / zeta + A - x / zeta]
    zero(g, 'supersolution_equality_at_true_pair_2', gt[1] - gt[0])
    zero(g, 'supersolution_equality_at_true_pair_3', gt[2] - gt[0])

    # ------------------------------------------------------------------ first-order rigidity at the Hopf point
    g = 'first-order rigidity'
    Dn = mu**2 + 6 * mu + 1
    al1, al2, be1, be2 = mu * A / Dn, -2 * mu / Dn, mu * (mu**2 + 4 * mu + 1) / Dn, -mu * A / Dn
    U1, p1 = sp.symbols('U1 p1')
    a1 = al1 * U1 + al2 * p1
    b1 = be1 * U1 + be2 * p1
    D1 = lambda F: sp.diff(F, U1) * p1 + sp.diff(F, p1) * (-mu * U1)
    zero(g, 'rigidity_eta_coefficients', D1(a1) - (b1 - A * a1))
    zero(g, 'rigidity_zeta_coefficients', D1(b1) - (mu * A * U1 - A * b1 - mu * a1))

    # ------------------------------------------------------------------ tau-time identities (exact primitives)
    g = 'period identities as exact derivatives'
    qv = pd                    # q = p'
    rv = pdd                   # p''
    kH = mu * A
    rho = k / kH
    zero(g, 'mean_x_primitive', D(y - p - mu * U) - (x - k), note='<x> = k, i.e. <e^U> = 1')
    Fb = p * qv + A * p**2 / 2 + x - k * U
    zero(g, 'identity_b_primitive', D(Fb) - (qv**2 - mu * p**2), note="<p'^2> = mu <p^2>")
    Fc = qv**2 / 2 + mu * p**2 / 2 + (x - k) * p
    zero(g, 'identity_c_primitive', D(Fc) - (x * p**2 - A * qv**2), note="A <p'^2> = <x p^2>")
    Fd = qv * p**2 + A * p**3 / 3 + Fc + A * Fb
    zero(g, 'identity_d_primitive', D(Fd) - ((k - kH) * p**2 - mu * p**3 + 2 * qv**2 * p),
         note="(k - k_H)<p^2> = mu <p^3> - 2 <p'^2 p>")
    Fe = -A * qv**2 / 2 - mu * p * qv - (x - k) * qv + x * p**2 / 2
    zero(g, 'identity_e_primitive', D(Fe) - (rv**2 - mu * qv**2 + x * p**3 / 2),
         note="<x p^3> = -2 (<p''^2> - mu <p'^2>)")
    F1 = qv * p**3 + A * p**4 / 4
    zero(g, 'identity_f_step1_primitive', D(F1) - ((k - x) * p**3 - mu * p**4 + 3 * qv**2 * p**2))
    F2 = -qv**2 * p / 2 - mu * p**3 / 3 - (x - k) * p**2 / 2
    zero(g, 'identity_f_step2_primitive', D(F2) - (A * qv**2 * p - qv**3 / 2 - x * p**3 / 2))
    Ff = A * Fd + F1 / rho - 2 * F2 - 2 * (1 - 1 / rho) * Fe
    If = (k - kH) * A * p**2 - (mu * p**4 - 3 * qv**2 * p**2) / rho + qv**3 - 2 * (1 - 1 / rho) * (rv**2 - mu * qv**2)
    zero(g, 'identity_f_primitive', D(Ff) - If,
         note="(k-k_H) A <p^2> = (mu<p^4> - 3<p'^2 p^2>)/rho - <p'^3> + 2(1 - 1/rho)(<p''^2> - mu <p'^2>)")
    # scaling to tau-time: P = p/sqrt(mu), Q = p'/mu, eps = sqrt(mu)/A
    sm = sp.Symbol('sqrt_mu', positive=True)
    Pm, Qm, Rm = sp.symbols('P Q Qprime')
    eps = sm / (1 + sm**2)
    sub = {p: sm * Pm, mu: sm**2}
    tform = ((rho - 1) * Pm**2 - eps**2 * ((Pm**4 - 3 * Qm**2 * Pm**2) / rho - Qm**3 + 2 * (1 - 1 / rho) * (Rm**2 - Qm**2))).subs(mu, sm**2)
    tt = sp.Symbol('pp'), sp.Symbol('ppp')
    Ifs = ((k - kH) * A * p**2 - (mu * p**4 - 3 * tt[0]**2 * p**2) / rho + tt[0]**3 - 2 * (1 - 1 / rho) * (tt[1]**2 - mu * tt[0]**2))
    Ifs = Ifs.subs({tt[0]: sm**2 * Qm, tt[1]: sm**3 * Rm}).subs(sub)
    zero(g, 'identity_f_time_rescaling', Ifs - sm**4 * (1 + sm**2)**2 * tform.subs(mu, sm**2),
         note='the t-time and tau-time forms of identity (f) differ by the factor mu^2 A^2')
    # kurtosis algebra of the wedge theorem
    t_ = sp.Symbol('t')
    m2_ = (1 + t_ + t_**2) / 3
    m4_ = (1 + t_ + t_**2 + t_**3 + t_**4) / 5
    zero(g, 'kurtosis_mixture_identity', sp.Rational(6, 5) * m2_ - sp.Rational(1, 5) - m4_ - t_ * (1 + t_) * (1 - t_**2) / 5)
    g_ = sp.Symbol('g')
    zero(g, 'kurtosis_square_identity', sp.Rational(9, 5) * g_**2 - sp.Rational(6, 5) * g_ + sp.Rational(1, 5) - (3 * g_ - 1)**2 / 5)

    # ------------------------------------------------------------------ large-gain algebra
    g = 'large-gain Lyapunov identity'
    nu_ = sp.Symbol('nu')
    aa, bb, cc = sp.symbols('a_c b_c c_c')       # cooperative system a' = -mu a + b, b' = -b + c, c' = -A c + x a
    alpha = 1 / (1 - nu_)
    beta = alpha / (A - nu_)
    Lf = aa + alpha * bb + beta * cc
    Ld = (-mu * aa + bb) + alpha * (-bb + cc) + beta * (-A * cc + x * aa)
    zero(g, 'trough_lyapunov_identity', Ld - (-nu_ * Lf + (beta * x - mu + nu_) * aa))
    zero(g, 'trough_threshold', (mu - nu_) / beta - (mu - nu_) * (1 - nu_) * (A - nu_))
    H1, H2, H3 = sp.symbols('H1 H2 H3')
    Bm = sp.Matrix([[0, 0, x], [1, mu, 0], [0, 1, 1]])
    Hd = Bm * sp.Matrix([H1, H2, H3])
    # (a,b,c) = e^{-At}(H3,H2,H1): derivative = e^{-At}(H' - A H)
    zero(g, 'cooperative_system_a', (Hd[2] - A * H3) - (-mu * H3 + H2))
    zero(g, 'cooperative_system_b', (Hd[1] - A * H2) - (-H2 + H1))
    zero(g, 'cooperative_system_c', (Hd[0] - A * H1) - (-A * H1 + x * H3))
    w_ = y - k
    zero(g, 'mass_identity_primitive', D(w_ - p - mu * U) - (x - k), note='integrated over the positive lobe')
    sg_, om_, s_t = sp.symbols('sigma omega s_time')
    Uf = sp.Function('U')(s_t)
    zf = sp.exp(-sg_ * s_t) * Uf
    lhs = sp.diff(zf, s_t, 2) + om_**2 * zf
    rhs = sp.exp(-sg_ * s_t) * (sp.diff(Uf, s_t, 2) - 2 * sg_ * sp.diff(Uf, s_t) + (sg_**2 + om_**2) * Uf)
    rows.append(dict(group=g, id='sturm_substitution', passed=bool(sp.simplify(lhs - rhs) == 0),
                     residual=str(sp.simplify(lhs - rhs)), note="z = e^{-sigma s} U:  z'' + omega^2 z = e^{-sigma s}(U'' - 2 sigma U' + zeta_E U)"))

    # ------------------------------------------------------------------ attractor / planar reduction
    g = 'planar reduction'
    Hf = sp.Function('H')
    u, v, thq = sp.symbols('u v theta_s')
    F = lambda a_, b_: sp.Matrix([b_, k - mu * b_ - Hf(a_, b_)])
    w2 = lambda a_, b_: a_[0] * b_[1] - a_[1] * b_[0]
    rows.append(dict(group=g, id='radial_scaling_wedge', note='F(theta q)^F(q) = -v[theta(H(q)-k)+k-H(theta q)]',
                     passed=bool(sp.expand(w2(F(thq * u, thq * v), F(u, v)) + v * (thq * (Hf(u, v) - k) + k - Hf(thq * u, thq * v))) == 0),
                     residual=str(sp.expand(w2(F(thq * u, thq * v), F(u, v)) + v * (thq * (Hf(u, v) - k) + k - Hf(thq * u, thq * v))))))
    F1f, F2f = [sp.Function(nm)(u, v) for nm in ('F1', 'F2')]
    a1_, a2_ = sp.symbols('a1 a2')
    DF = sp.Matrix([[sp.diff(F1f, u), sp.diff(F1f, v)], [sp.diff(F2f, u), sp.diff(F2f, v)]])
    Fv = sp.Matrix([F1f, F2f])
    Jv = sp.Matrix([a1_, a2_])
    dw = w2(DF * Fv, Jv) + w2(Fv, DF * Jv)
    rr = sp.expand(dw - (DF[0, 0] + DF[1, 1]) * w2(Fv, Jv))
    rows.append(dict(group=g, id='planar_wedge_transport', passed=bool(rr == 0), residual=str(rr), note="(F^J)' = (div F)(F^J)"))
    hh = sp.Symbol('h')
    sol = sp.solve(sp.Eq(-eta * (k - hh), x - hh), hh)[0]
    zero(g, 'crest_axis_graph_value', sol - k - (x - k) / (1 + eta), note='h(u,0) - k = k(e^u - 1)/(1+eta)')
    # two-orbit comparison (supporting-plane identity), with x e^{dU} written as x2
    Delta = zeta * (U2 - U) - eta * (p2 - p) - (y2 - y)
    zero(g, 'two_orbit_identity', D(Delta) + (1 + eta) * Delta + (x2 - x - x * (U2 - U)))
    nvec = sp.Matrix([zeta, -eta, -1])
    resn = sp.Matrix([D(v_i) for v_i in nvec]) + J.T * nvec + (1 + eta) * nvec
    for i in range(3):
        zero(g, f'normal_equation_row{i+1}', resn[i])
    zero(g, 'planar_divergence', sp.diff(k - mu * v - Hf(u, v), v) - (-mu - sp.diff(Hf(u, v), v)),
         note='div F = -mu - H_v = eta - mu when H_v = -eta')

    if inject_error:
        modO('deliberate error', 'injected_wrong_sign_star_gauge', D(Sg) - ((eta - mu) * Sg - p * K),
             note='wrong sign of p K; must fail')

    groups = {}
    for r_ in rows:
        gcount = groups.setdefault(r_['group'], [0, 0])
        gcount[0] += 1
        gcount[1] += int(r_['passed'])
    extra = {'software': {'sympy': sp.__version__},
             'groups': {gname: {'attempted': v[0], 'passed': v[1]} for gname, v in groups.items()},
             'counting_rule': 'one count per residual tested equal to zero; matrix identities are counted row by row'}
    return write_receipt(name, 'exact', rows, worst={'nonzero_residuals': sum(not r_['passed'] for r_ in rows)},
                         parameters='symbolic in mu, k and all state variables', clock=clock, extra=extra, out=out,
                         scope='Finite algebra only. Positivity, averaging and topological steps are in the written proofs.')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--inject-error', action='store_true')
    ap.add_argument('--out', type=Path, default=None)
    args = ap.parse_args()
    ok = run(args.inject_error, args.out, name='v3_exact_identities' if not args.inject_error else 'v3_exact_identities_with_injected_error')
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
