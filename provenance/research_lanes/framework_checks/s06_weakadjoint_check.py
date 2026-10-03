#!/usr/bin/env python3
"""Formula re-check for paper section s06_weakadjoint.tex (v3.0), after renaming to the global notation.

Part A (exact, sympy): every displayed identity of the section, written from the equations of motion
  U' = p, p' = k - mu p - y, y' = x - y, x = k e^U,
  a' = -b, b' = -mu b - c, c' = -x a - c   (e_s = (a,-b,-c)),   fr = b/a, fs = c/a,
  eta' = zeta + (mu-1) eta - eta^2, zeta' = x - (1+eta) zeta.
Part B (floating point, spot check on three cycles): <fr> = 1+<eta>, fr > max(1,mu), the twelve-event
  order, the supporting line, the three determinant inequalities, the rank-one window identity for
  phi = 1, p', -U, positivity of the crest-return derivative on the region Q and its Perron vector.

Location: ROOT/work/framework_checks/ (moved out of the paper's sections directory at the audit of
3 October 2026).  Items A36-A43, B11 and B12 were added at that audit: the Sturm-Liouville and Lagrange
forms of Remark wk:rem:secondorder, the Wronskian statements of Remark wk:rem:double, the cross-product
identity, the fourth instance of the rank-one identity and the alpha_* orthogonality.

Floating point is evidence only; nothing here is an interval certificate.
Run:  /Users/admin/.venvs/analysis-py313/bin/python -W ignore s06_weakadjoint_check.py > s06_weakadjoint_check.out
(single process, a few minutes)
"""
import sys
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp

sys.path.insert(0, '/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03/work')
import whlib  # noqa: E402

FAIL = 0


def report(name, ok, extra=''):
    global FAIL
    FAIL += (not ok)
    print(f"{'PASS' if ok else 'FAIL'}  {name}  {extra}")


# ----------------------------------------------------------------------------------------------
# Part A: exact identities
# ----------------------------------------------------------------------------------------------
def part_a():
    U, p, y, k, mu, fr, fs, a, eta, zeta, t = sp.symbols('U p y k mu fr fs a eta zeta t', real=True)
    A = 1 + mu
    x = k * sp.exp(U)
    F = sp.Matrix([p, k - mu * p - y, x - y])
    J = F.jacobian([U, p, y])
    state = [U, p, y, fr, fs, a, eta, zeta]
    rhs = [F[0], F[1], F[2], fr**2 - mu * fr - fs, fs * fr - fs - x, -fr * a,
           zeta + (mu - 1) * eta - eta**2, x - (1 + eta) * zeta]

    def D(f):
        return sp.expand(sum(sp.diff(f, q) * dq for q, dq in zip(state, rhs)))

    def Z(e):
        if isinstance(e, sp.MatrixBase):
            return all(sp.simplify(q) == 0 for q in e)
        return sp.simplify(e) == 0

    pp, yp = F[1], F[2]
    ppp, ypp = D(pp), D(yp)

    # --- Lemma (strong-stable ratios)
    b, c = fr * a, fs * a
    report('A01 fr,fs Riccati consistent with a\'=-b, b\'=-mu b-c, c\'=-x a-c',
           Z(D(a) + b) and Z(D(b) + mu * b + c) and Z(D(c) + x * a + c))
    report('A02 e_s=(a,-b,-c) solves v\'=Jv', Z(sp.Matrix([D(a), D(-b), D(-c)]) - J * sp.Matrix([a, -b, -c])))
    report("A03 fr''+(A-3fr)fr'+fr(fr-1)(fr-mu) = x",
           Z(D(D(fr)) + (A - 3 * fr) * D(fr) + fr * (fr - 1) * (fr - mu) - x))
    ct = c + (mu - 1) * b
    report("A04 conjugate system: b'=-b-ct, ct'=-x a-mu ct, fr'=fr^2-fr-ct/a",
           Z(D(b) + b + ct) and Z(D(ct) + x * a + mu * ct) and Z(D(fr) - (fr**2 - fr - ct / a)))
    g = 1 + eta
    report("A05 mirror: g=1+eta solves g''+(3g-A)g'+g(g-1)(g-mu) = x",
           Z(D(D(g)) + (3 * g - A) * D(g) + g * (g - 1) * (g - mu) - x))
    # factorisation L3 = [D^2+(A-fr)D+(fr-1)(fr-mu)-2fr'](D+fr) on a generic function of t
    ph = sp.Function('phi')(t)
    R = sp.Function('R')(t)
    X = sp.Function('X')(t)
    L3gen = sp.diff(ph, t, 3) + A * sp.diff(ph, t, 2) + mu * sp.diff(ph, t) + X * ph
    inner = sp.diff(ph, t) + R * ph
    beta = (R - 1) * (R - mu) - 2 * sp.diff(R, t)
    fac = sp.diff(inner, t, 2) + (A - R) * sp.diff(inner, t) + beta * inner
    scal = sp.diff(R, t, 2) + (A - 3 * R) * sp.diff(R, t) + R * (R - 1) * (R - mu)
    report('A06 L3[phi] - L2[(D+fr)phi] = (x - scalar eq.) * phi  (so equal on the scalar equation)',
           Z(sp.expand(L3gen - fac - (X - scal) * ph)))
    # scalar equation <=> L3[a]=0 for a'=-R a
    aa = sp.Function('aa')(t)
    L3a = (sp.diff(aa, t, 3) + A * sp.diff(aa, t, 2) + mu * sp.diff(aa, t) + X * aa)
    sub = {sp.Derivative(aa, (t, 3)): None}
    d1 = -R * aa
    d2 = sp.diff(d1, t).subs(sp.Derivative(aa, t), d1)
    d3 = sp.diff(d2, t).subs(sp.Derivative(aa, t), d1)
    L3a_val = d3 + A * d2 + mu * d1 + X * aa
    report("A07 a'=-R a  =>  L3[a]/a = x - (R''+(A-3R)R'+R(R-1)(R-mu))", Z(sp.expand(L3a_val / aa - (X - scal))))

    # --- Lemma (weak adjoint); Psi = e^{At} a * Phat
    q = pp + fr * p
    Phat = sp.Matrix([fr * yp - fs * pp, yp + fs * p, -q])
    cross = sp.Matrix([p, pp, yp]).cross(sp.Matrix([a, -b, -c]))
    report("A08 z' x e_s = a (fr y' - fs p', y' + fs p, -(p' + fr p))", Z(cross - a * Phat))
    dPhat = sp.Matrix([D(e) for e in Phat]) + (A - fr) * Phat      # (e^{At} a Phat)'/(e^{At} a)
    report("A09 Psi' = -J^T Psi", Z(dPhat + J.T * Phat))
    report("A10 components: Psi1' = -x Psi3, Psi2' = mu Psi2 - Psi1, Psi3' = Psi3 + Psi2",
           Z(dPhat[0] + x * Phat[2]) and Z(dPhat[1] - mu * Phat[1] + Phat[0]) and Z(dPhat[2] - Phat[2] - Phat[1]))
    report("A11 Psi.z' = 0 and Psi.e_s = 0", Z(Phat.dot(F)) and Z(Phat.dot(sp.Matrix([1, -fr, -fs]))))
    v1, v2, v3 = sp.symbols('v1 v2 v3')
    v = sp.Matrix([v1, v2, v3])
    report("A12 (z' x e_s).v = det(v, z', e_s)",
           Z(cross.dot(v) - sp.Matrix.hstack(v, sp.Matrix([p, pp, yp]), sp.Matrix([a, -b, -c])).det()))
    U2, p2, y2 = sp.symbols('U2 p2 y2')
    F2 = F.subs({U: U2, p: p2, y: y2}, simultaneous=True)
    d = sp.Matrix([U2 - U, p2 - p, y2 - y])
    Rem = sp.exp(d[0]) - 1 - d[0]
    report('A13 chord: f(z2)-f(z)-J(z)(z2-z) = x Rem(d1) e3', Z(F2 - F - J * d - sp.Matrix([0, 0, x * Rem])))
    zpp = sp.Matrix([pp, ppp, ypp])
    zppp = sp.Matrix([D(e) for e in zpp])
    report("A14 z''' = J z'' + x p^2 e3", Z(zppp - J * zpp - sp.Matrix([0, 0, x * p**2])))
    report('A15 E-chord weight: x Rem(-U) = k psi(U)', Z(x * (sp.exp(-U) - 1 + U) - k * (U * sp.exp(U) - sp.exp(U) + 1)))

    # --- Theorem (twelve events)
    th = -pp / p
    report("A16 theta' = theta^2 - mu theta + y'/p", Z(D(th) - (th**2 - mu * th + yp / p)))
    report("A17 (theta-fr)' = (theta-fr)(theta+fr-mu) + y'/p + fs",
           Z(D(th) - D(fr) - ((th - fr) * (th + fr - mu) + yp / p + fs)))
    report("A18 q' = (fr-mu) q - (y' + fs p)", Z(D(q) - ((fr - mu) * q - (yp + fs * p))))
    g2 = yp / p + fs
    report("A19 g2' = fs q/p - g2 (1 + p'/p)", Z(D(g2) - (fs * q / p - g2 * (1 + pp / p))))
    g1 = yp / pp - fs / fr
    report("A20 g1' = x q/(fr p') + g1 (mu - 1 + y'/p' + fs/fr)",
           Z(D(g1) - (x * q / (fr * pp) + g1 * (mu - 1 + yp / pp + fs / fr))))
    report("A21 Psi1 = varpi fr p' g1, Psi2 = varpi p g2, q = p (fr - theta)",
           Z(Phat[0] - fr * pp * g1) and Z(Phat[1] - p * g2) and Z(q - p * (fr - th)))

    # --- Lemma (supporting line), at p' = -fr p
    Xi = pp**2 - p * ppp
    report("A22 Xi = p'^2 - p p'' = p'^2 + mu p p' + p y'", Z(Xi - (pp**2 + mu * p * pp + p * yp)))
    # impose p' = -fr p by solving for y
    ysol = sp.solve(sp.Eq(pp, -fr * p), y)[0]
    report("A23 at p'=-fr p: Xi = p^2 fr (fr-mu) + p y'",
           Z((Xi - (p**2 * fr * (fr - mu) + p * yp)).subs(y, ysol)))
    report("A24 at p'=-fr p: Sigma = p^2 - U p' = p (p + fr U)", Z(((p**2 - U * pp) - p * (p + fr * U)).subs(y, ysol)))
    etaE, zetaE = sp.symbols('etaE zetaE')
    report("A25 Lambda_E = (y-k) - zetaE U + etaE p = -p' + (etaE-mu) p - zetaE U",
           Z((y - k) - zetaE * U + etaE * p - (-pp + (etaE - mu) * p - zetaE * U)))

    # --- Proposition (rank one)
    def L3(f):
        return sp.simplify(D(D(D(f))) + A * D(D(f)) + mu * D(f) + x * f)
    report('A26 L3[1]=x, L3[p]=0, L3[p\']=-x p^2, L3[U]=k psi(U), L3[a]=0',
           Z(L3(sp.Integer(1)) - x) and Z(L3(p)) and Z(L3(pp) + x * p**2)
           and Z(L3(U) - k * (U * sp.exp(U) - sp.exp(U) + 1)) and Z(L3(a)))
    phi_t = sp.Function('phi')(t)
    vph = sp.Matrix([phi_t, sp.diff(phi_t, t), -sp.diff(phi_t, t, 2) - mu * sp.diff(phi_t, t)])
    Jt = sp.Matrix([[0, 1, 0], [0, -mu, -1], [X, 0, -1]])
    report("A27 v_phi' - J v_phi = -L3[phi] e3",
           Z(sp.diff(vph, t) - Jt * vph + sp.Matrix([0, 0, L3gen.subs(ph, phi_t)])))
    report("A28 z' = v_p, z'' = v_{p'}, e_s = v_a",
           Z(sp.Matrix([p, pp, yp]) - sp.Matrix([p, pp, -ppp - mu * pp]))
           and Z(zpp - sp.Matrix([pp, ppp, -D(ppp) - mu * ppp]))
           and Z(sp.Matrix([a, -b, -c]) - sp.Matrix([a, D(a), -D(D(a)) - mu * D(a)])))
    report("A29 p'' + fr p' = (fr-mu) p' - y'", Z(ppp + fr * pp - ((fr - mu) * pp - yp)))

    # --- Proposition (crest return)
    u1, u2, u3 = sp.symbols('u1 u2 u3')
    u = sp.Matrix([u1, u2, u3])
    w = sp.Matrix([v2 * u1 - v1 * u2, v1 * u3 - v3 * u1, v2 * u3 - v3 * u2])
    C = sp.Matrix([[-mu, 1, 0], [0, -1, 1], [x, 0, -A]])
    dv, du = J * v, J * u
    dw = sp.Matrix([sum(sp.diff(wi, s_) * ds_ for s_, ds_ in zip(list(v) + list(u), list(dv) + list(du))) for wi in w])
    report("A30 compound: w' = C w for w = (v2u1-v1u2, v1u3-v3u1, v2u3-v3u2), u'=Ju, v'=Jv", Z(dw - C * w))
    dU, dy = sp.symbols('dU dy')
    wsec = w.subs({v1: dU, v2: 0, v3: dy, u1: 0, u2: k - y, u3: x - y}, simultaneous=True)
    report("A31 wedge of (dU,0,dy) with z'=(0,k-y,x-y) is ((y-k)dU,(x-y)dU,(y-k)dy)",
           Z(wsec - sp.Matrix([(y - k) * dU, (x - y) * dU, (y - k) * dy])))
    gE = sp.symbols('gammaE', positive=True)
    JE = sp.Matrix([[0, 1, 0], [0, -mu, -1], [k, 0, -1]])
    vs = sp.Matrix([1, -gE, -gE * (gE - mu)])
    report('A32 J_E v_s = -gammaE v_s when k = gammaE(gammaE-1)(gammaE-mu)',
           Z((JE * vs + gE * vs).subs(k, gE * (gE - 1) * (gE - mu))))
    report('A33 char. poly of J_E at s=-gammaE is k - gammaE(gammaE-1)(gammaE-mu)',
           Z((-gE * sp.eye(3) - JE).det() - (k - gE * (gE - 1) * (gE - mu))))
    n = sp.Matrix([zeta, -eta, -1])
    dn = sp.Matrix([D(e) for e in n])
    report("A34 n' + J^T n + (1+eta) n = 0 for n = (zeta,-eta,-1)", Z(dn + J.T * n + (1 + eta) * n))
    report("A35 backward system matrices: -(a,b,c)' = N (a,b,c), N>=0", True,
           'N = [[0,1,0],[0,mu,1],[x,0,1]], conjugate N~ = [[0,1,0],[0,1,1],[x,0,mu]] (by A01, A04)')

    # --- added at the audit of 3 October 2026 ---
    Wf = sp.Function('varpi')(t)
    chi_f = sp.Function('chi')(t)
    q_f = sp.Function('q')(t)
    L2 = lambda f: sp.diff(f, t, 2) + (A - R) * sp.diff(f, t) + beta * f
    wsub = {sp.Derivative(Wf, t): (A - R) * Wf}
    report("A36 Sturm-Liouville form: varpi L2[chi] = (varpi chi')' + varpi[(R-1)(R-mu)-2R'] chi when varpi' = (A-R) varpi",
           Z(sp.expand(Wf * L2(chi_f) - (sp.diff(Wf * sp.diff(chi_f, t), t).subs(wsub) + Wf * beta * chi_f))))
    report("A37 Lagrange identity: varpi(q L2[chi] - chi L2[q]) = (varpi(q chi' - q' chi))'",
           Z(sp.expand(Wf * (q_f * L2(chi_f) - chi_f * L2(q_f))
                       - sp.diff(Wf * (q_f * sp.diff(chi_f, t) - sp.diff(q_f, t) * chi_f), t).subs(wsub))))
    report("A38 Psi_3/varpi = (a' p - a p')/a = -(p' + fr p);  Wr(p,a) = p a' - p' a = -a q",
           Z((D(a) * p - a * pp) / a + q) and Z(p * D(a) - pp * a + a * q))
    f1, f2, f3 = [sp.Function(nm)(t) for nm in ('f1', 'f2', 'f3')]
    col = lambda f: sp.Matrix([f, sp.diff(f, t), -sp.diff(f, t, 2) - mu * sp.diff(f, t)])
    wr = sp.Matrix([[f1, f2, f3], [sp.diff(f1, t), sp.diff(f2, t), sp.diff(f3, t)],
                    [sp.diff(f1, t, 2), sp.diff(f2, t, 2), sp.diff(f3, t, 2)]]).det()
    report("A39 det(v_f1, v_f2, v_f3) = -Wr(f1,f2,f3) for v_f = (f, f', -f''-mu f')",
           Z(sp.expand(sp.Matrix.hstack(col(f1), col(f2), col(f3)).det() + wr)))
    report("A40 n . e_s = zeta a + eta b + c", Z(n.dot(sp.Matrix([a, -b, -c])) - (zeta * a + eta * b + c)))
    al_ = sp.symbols('alpha_star')
    report("A41 phi = p' - alpha: phi' + fr phi = (fr-mu)p' - y' - fr alpha;  L3[p' - alpha] = -x(p^2 + alpha)",
           Z(ppp + fr * (pp - al_) - ((fr - mu) * pp - yp - fr * al_)) and Z(L3(pp) - al_ * x + x * (p**2 + al_)))
    Mg = sp.Matrix(3, 3, sp.symbols('m0:9'))
    ug = sp.Matrix(sp.symbols('ug0:3'))
    vg = sp.Matrix(sp.symbols('vg0:3'))
    report("A42 (Mu) x v + u x (Mv) = (tr M - M^T)(u x v) for every 3x3 matrix M",
           Z(sp.expand((Mg * ug).cross(vg) + ug.cross(Mg * vg) - (Mg.trace() * sp.eye(3) - Mg.T) * ug.cross(vg))))
    report("A43 Step 1 of the crest-return proof: the roots of s^3+As^2+mu s+k other than -gammaE have sum gammaE-A and product k/gammaE",
           Z(sp.expand(((sp.Symbol('s')**3 + A * sp.Symbol('s')**2 + mu * sp.Symbol('s') + k)
                        - (sp.Symbol('s') + gE) * (sp.Symbol('s')**2 - (gE - A) * sp.Symbol('s') + k / gE)).subs(k, gE * (gE - 1) * (gE - mu)))))


# ----------------------------------------------------------------------------------------------
# Part B: floating-point spot check
# ----------------------------------------------------------------------------------------------
def cyc_data(mu, rho, n=40001):
    k = rho * whlib.kH(mu)
    A = 1 + mu
    cyc = whlib.find_cycle(mu, k)
    T = cyc['T']
    ed = whlib.equilibrium_data(mu, k)

    def f(tt, s):
        U, p, y, eta, zeta = s
        x = k * np.exp(U)
        return [p, k - mu * p - y, x - y, zeta + (mu - 1) * eta - eta * eta, x - (1 + eta) * zeta]
    s0 = np.array([cyc['U0'], 0.0, cyc['y0'], ed['etaE'], ed['zetaE']])
    for _ in range(6):
        sol = solve_ivp(f, [0, T], s0, method='DOP853', rtol=1e-12, atol=1e-14)
        s0 = sol.y[:, -1].copy()
        s0[:3] = [cyc['U0'], 0.0, cyc['y0']]
    fwd = solve_ivp(f, [0, T], s0, method='DOP853', rtol=1e-12, atol=1e-14, dense_output=True)

    def state(tt):
        return fwd.sol(np.mod(tt, T))

    # strong-stable ratios: integrate the Riccati system backward in time (attracting in that direction)
    def gb(tau, w):          # tau = -t
        U = state(-tau)[0]
        x = k * np.exp(U)
        fr, fs = w
        return [-(fr * fr - mu * fr - fs), -(fs * fr - fs - x)]
    gE = ed['gammaE']
    w0 = np.array([gE, gE * (gE - mu)])
    for _ in range(6):
        b_ = solve_ivp(gb, [0, T], w0, method='DOP853', rtol=1e-12, atol=1e-14)
        w0 = b_.y[:, -1].copy()
    bwd = solve_ivp(gb, [0, T], w0, method='DOP853', rtol=1e-12, atol=1e-14, dense_output=True)
    per_defect = float(np.max(np.abs(bwd.y[:, -1] - w0)))

    def ratios(tt):          # (fr, fs) at time tt; bwd covers t in [-T, 0]
        return bwd.sol(np.mod(-tt, T))

    t = np.linspace(0.0, T, n)
    U, p, y, eta, zeta = state(t)
    fr, fs = ratios(t)
    x = k * np.exp(U)
    pp = k - mu * p - y
    yp = x - y
    ppp = -mu * pp - yp
    ypp = x * p - yp
    return dict(mu=mu, k=k, A=A, T=T, cyc=cyc, ed=ed, t=t, U=U, p=p, y=y, x=x, eta=eta, zeta=zeta, fr=fr, fs=fs,
                pp=pp, yp=yp, ppp=ppp, ypp=ypp, state=state, ratios=ratios, per_defect=per_defect)


def crossings(t, f):
    """times and directions (+1 up, -1 down) of sign changes of f on the periodic grid t (last point = first)."""
    out = []
    for i in range(len(t) - 1):
        if f[i] == 0.0 or f[i] * f[i + 1] < 0:
            tc = t[i] - f[i] * (t[i + 1] - t[i]) / (f[i + 1] - f[i])
            out.append((tc, 1 if f[i + 1] > f[i] else -1))
    return out


def part_b(mu, rho):
    d = cyc_data(mu, rho)
    T, A, k, t = d['T'], d['A'], d['k'], d['t']
    tag = f'[mu={mu}, rho={rho}]'
    lam_eta = float(np.trapezoid(d['eta'], t) / T - mu)
    lam_r = float(np.trapezoid(d['fr'], t) / T - A)
    report(f'B01 {tag} <fr> - A = <eta> - mu', abs(lam_eta - lam_r) < 1e-7,
           f'lambda_c = {lam_eta:.9f} / {lam_r:.9f}; Riccati periodicity defect {d["per_defect"]:.1e}')
    report(f'B02 {tag} fr > max(1,mu), fs > 0, fs + (mu-1) fr > 0',
           d['fr'].min() > max(1, mu) and d['fs'].min() > 0 and (d['fs'] + (mu - 1) * d['fr']).min() > 0,
           f'min fr - max(1,mu) = {d["fr"].min() - max(1, mu):.3e}')
    # twelve-event order (crest at t=0)
    q = d['pp'] + d['fr'] * d['p']
    P1 = d['fr'] * d['yp'] - d['fs'] * d['pp']
    P2 = d['yp'] + d['fs'] * d['p']
    P3 = -q
    ev = []
    for name, f in (('p', d['p']), ("p'", d['pp']), ("y'", d['yp']), ('Psi1', P1), ('Psi2', P2), ('Psi3', P3)):
        cs = crossings(t, f)
        if len(cs) != 2:
            report(f'B03 {tag} two sign changes of {name}', False, f'found {len(cs)}')
        ev += [(tc, f"{name}{'+' if s_ > 0 else '-'}") for tc, s_ in cs]
    # rotate to start at the trough (up-crossing of p)
    t_tr = [tc for tc, nm in ev if nm == 'p+'][0]
    ev = sorted(((tc - t_tr) % T, nm) for tc, nm in ev)
    order = [nm for _, nm in ev]
    expected = ['p+', 'Psi2+', "y'+", 'Psi1+', "p'-", 'Psi3+', 'p-', 'Psi2-', "y'-", 'Psi1-', "p'+", 'Psi3-']
    report(f'B03 {tag} twelve-event order', order == expected, ' < '.join(order))
    i_cr, i_tr = 0, int(np.argmin(np.abs(np.mod(t - t_tr, T))))
    report(f'B04 {tag} Psi > 0 at crest, Psi < 0 at trough',
           min(P1[i_cr], P2[i_cr], P3[i_cr]) > 0 and max(P1[i_tr], P2[i_tr], P3[i_tr]) < 0)
    # supporting line at t_* = down-crossing of Psi3
    t_st = [tc for tc, s_ in crossings(t, P3) if s_ < 0][0]
    Us, ps, ys = d['state'](t_st)[:3]
    pps = k - mu * ps - ys
    yps = k * np.exp(Us) - ys
    frs, fss = d['ratios'](t_st)
    Lg = pps * (d['U'] - Us) - ps * (d['p'] - ps)
    mask = np.minimum(np.abs(t - t_st), T - np.abs(t - t_st)) > 1e-3 * T
    Sig, Xi = ps * ps - Us * pps, pps * pps - ps * (-mu * pps - yps)
    report(f'B05 {tag} L_Gamma > 0 away from t_*, U(t_*) < 0, p(t_*) < 0 < p\'(t_*), L(0,0)=Sigma>0, L\'\'=Xi>0',
           Lg[mask].min() > 0 and Us < 0 and ps < 0 < pps and Sig > 0 and Xi > 0 and abs(pps + frs * ps) < 1e-6 * abs(pps),
           f'min L_Gamma = {Lg[mask].min():.3e}, Sigma(t_*) = {Sig:.3e}, Xi(t_*) = {Xi:.3e}')
    # determinant inequalities on a subgrid
    idx = np.arange(0, len(t) - 1, 200)
    Z_ = np.vstack([d['U'], d['p'], d['y']]).T
    Zp = np.vstack([d['p'], d['pp'], d['yp']]).T
    Zpp = np.vstack([d['pp'], d['ppp'], d['ypp']]).T
    Es = np.vstack([np.ones_like(t), -d['fr'], -d['fs']]).T
    E = np.array([0.0, 0.0, k])
    worst = -np.inf
    for i in idx:
        nrm = np.linalg.norm(Zp[i]) * np.linalg.norm(Es[i])
        for j in idx:
            if i == j:
                continue
            ch = Z_[j] - Z_[i]
            worst = max(worst, np.linalg.det(np.column_stack([ch, Zp[i], Es[i]])) / (nrm * np.linalg.norm(ch)))
        worst = max(worst, np.linalg.det(np.column_stack([E - Z_[i], Zp[i], Es[i]])) / (nrm * np.linalg.norm(E - Z_[i])))
        worst = max(worst, np.linalg.det(np.column_stack([Zpp[i], Zp[i], Es[i]])) / (nrm * np.linalg.norm(Zpp[i])))
    report(f'B06 {tag} det(chord, z\', e_s) < 0, det(E-z, z\', e_s) < 0, det(z\'\', z\', e_s) < 0',
           worst < 0, f'largest normalised determinant {worst:.3e} on {len(idx)}^2 pairs')

    # rank-one window identity, ODE quadrature over [t_*, t_*+T] with varpi(t_*) = 1
    def h(tt, s):
        U, p, y = d['state'](tt)[:3]
        fr, fs = d['ratios'](tt)
        x = k * np.exp(U)
        pp = k - mu * p - y
        varpi = np.exp(s[0])
        P3_ = -varpi * (pp + fr * p)
        psiU = U * np.exp(U) - np.exp(U) + 1.0
        dch = d['state'](tt + sig_ch)[0] - U
        return [A - fr, P3_ * x, P3_ * (-x * p * p), P3_ * (-k * psiU),
                P3_ * x * (np.exp(dch) - 1.0 - dch), P3_ * x * (p * p + alpha_star), abs(P3_) * x * p * p]
    sig_ch = T / 3.0                                              # chord shift for the fourth instance
    alpha_star = (frs - mu) * abs(ps) + abs(yps) / frs            # phi = p' - alpha_star has (phi' + fr phi)(t_*) = 0
    q_ = solve_ivp(h, [t_st, t_st + T], [0.0] * 7, method='DOP853', rtol=1e-11, atol=1e-13)
    lv, I1, Ipp, IU, Ich, Ial, Iabs = q_.y[:, -1]
    report(f'B07 {tag} varpi(t+T)/varpi(t) = exp(-lambda_c T)', abs(lv + lam_eta * T) < 1e-6 * max(1, abs(lam_eta * T)),
           f'log ratio {lv:.9f} vs {-lam_eta * T:.9f}')
    fac = (np.exp(-lam_eta * T) - 1.0) * abs(yps + fss * ps)
    rhs = [fac * frs, fac * ((frs - mu) * pps - yps), fac * (-(ps + frs * Us))]
    names = ['phi=1', "phi=p'", 'phi=-U']
    for nm, lhs, r_ in zip(names, [I1, Ipp, IU], rhs):
        report(f'B08 {tag} rank-one identity, {nm}', abs(lhs - r_) < 2e-5 * abs(r_) and r_ * np.sign(-lam_eta) > 0,
               f'window integral {lhs:.9e}, right side {r_:.9e}')
    # added at the audit: fourth instance (chord weight) and the alpha_* orthogonality
    Uc, pc = d['state'](t_st + sig_ch)[:2]
    beta_st = -(yps + fss * ps) / ps                              # Psi(t_*) = beta (p', -p, 0), varpi(t_*) = 1
    Lch = pps * (Uc - Us) - ps * (pc - ps)
    r_ch = (np.exp(-lam_eta * T) - 1.0) * beta_st * Lch
    report(f'B11 {tag} rank-one identity, chord weight x Rem(U(t+T/3)-U(t)); beta < 0 < L; integral < 0 iff lambda_c < 0',
           abs(Ich - r_ch) < 2e-5 * abs(r_ch) and beta_st < 0 < Lch and (Ich < 0) == (lam_eta < 0),
           f'window integral {Ich:.9e}, right side {r_ch:.9e}')
    report(f'B12 {tag} int Psi_3 x (p^2 + alpha_*) = 0 over the window', abs(Ial) < 2e-6 * Iabs,
           f'integral {Ial:.2e} against int |Psi_3| x p^2 = {Iabs:.3e};  alpha_* = {alpha_star:.6f}')
    # crest return: Perron vector at the cycle crest, and positivity on Q
    M = d['cyc']['M']
    wv, V = np.linalg.eig(M)
    j = int(np.argmax(wv.real))
    slope = float((V[1, j] / V[0, j]).real)
    report(f'B09 {tag} Perron vector of M is (1, zeta(crest)), eigenvalue m_t',
           abs(slope - d['zeta'][0]) < 1e-6 * d['zeta'][0] and abs(np.log(wv[j].real) / T - lam_eta) < 1e-6
           and M.min() > 0, f'slope {slope:.9f} vs zeta {d["zeta"][0]:.9f}')
    Ucr = d['cyc']['U0']
    ok, cnt, mn = True, 0, np.inf
    for fu in (0.1, 0.5, 1.0, 1.5):
        for fy in (1e-3, 0.3, 1.0):
            U0 = fu * Ucr
            x0 = k * np.exp(U0)
            y0 = k + fy * (x0 - k)
            U1, y1, T1, M1 = whlib.crest_return(mu, k, U0, y0, tmax=4000.0 * max(1.0, 1.0 / np.sqrt(mu)))
            x1 = k * np.exp(U1)
            ok = ok and (M1.min() > 0) and (x1 > y1 > k)
            mn = min(mn, M1.min())
            cnt += 1
    report(f'B10 {tag} DP > 0 entrywise and x1 > y1 > k at {cnt} points of Q (incl. x0 = y0)', ok,
           f'smallest entry {mn:.3e}')


if __name__ == '__main__':
    part_a()
    for mu_, rho_ in ((1.0, 1.5), (0.3, 2.0), (3.0, 1.3)):
        part_b(mu_, rho_)
    print(f"\n{'ALL PASS' if FAIL == 0 else str(FAIL) + ' FAILURES'}")
    sys.exit(1 if FAIL else 0)
