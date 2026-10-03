#!/usr/bin/env python3
"""Re-check of the formulas of paper/sections/s11_wedge.tex after renaming to the global notation.

Part A: exact symbolic identities (sympy).  A failure here is an error in the paper.
Part B: floating-point spot checks (ordinary double precision; NOT interval-certified).

Location: ROOT/work/framework_checks/ (moved out of the paper's sections directory at the audit of
3 October 2026; ROOT = /Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03).  The script also runs from a
copy inside the bundle (BUNDLE/provenance/research_lanes/framework_checks/): the libraries whlib.py, toda.py
and ampcycle.py are looked up first in BUNDLE/code/v3/lib and then in the research-lane directories.
Blocks A13-A16 and B5-B6 were added at that audit (cross-model review of the section): the equality case of
the kurtosis bound, the polar form of the planar field, the energy drift used for the containment of the
selected cycle, the pulse identity, and the winding, crest and energy-excursion statements on computed cycles.

Single process, a few seconds.  Run:
    /Users/admin/.venvs/analysis-py313/bin/python s11_wedge.check.py
Output is printed and copied to s11_wedge.check.out next to this file.
"""
from __future__ import annotations
import sys, time, io
from pathlib import Path

HERE = Path(__file__).resolve()
ROOT_ABS = Path('/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03')


def _library_dirs():
    """Directories that may hold whlib.py, toda.py, ampcycle.py, most specific first."""
    cands = []
    for anc in HERE.parents:                                   # the bundle, if this file sits inside it
        cands.append(anc / 'code' / 'v3' / 'lib')
    for anc in HERE.parents:                                   # the bundle, seen from ROOT/work/framework_checks
        cands.append(anc / 'v3' / 'wilhelm_heinrich_v3_0' / 'code' / 'v3' / 'lib')
    par = HERE.parent.parent                                   # ROOT/work  or  BUNDLE/provenance/research_lanes
    cands += [par / 'lanes' / 'hamiltonian', par / 'hamiltonian', par,
              ROOT_ABS / 'v3' / 'wilhelm_heinrich_v3_0' / 'code' / 'v3' / 'lib',
              ROOT_ABS / 'work' / 'lanes' / 'hamiltonian', ROOT_ABS / 'work']
    out = []
    for c in cands:
        if c.is_dir() and c not in out and any((c / n).is_file() for n in ('whlib.py', 'toda.py', 'ampcycle.py')):
            out.append(c)
    return out


LIBDIRS = _library_dirs()
for cand in LIBDIRS:
    if str(cand) not in sys.path:
        sys.path.append(str(cand))

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

LOG = io.StringIO()
FAILS = []


def say(*a):
    line = ' '.join(str(x) for x in a)
    print(line, flush=True)
    LOG.write(line + '\n')


def check(name, ok, detail=''):
    say(('PASS ' if ok else 'FAIL ') + name + ('   ' + detail if detail else ''))
    if not ok:
        FAILS.append(name)


# =====================================================================================================
# Part A: exact symbolic identities
# =====================================================================================================
t0 = time.time()
say('=== Part A: exact symbolic identities (sympy) ===')
U, P, Q, eps, r, rho, mu, k = sp.symbols('U P Q varepsilon r rho mu k', real=True)
g = sp.exp(U) - 1
Rem = lambda u: sp.exp(u) - u - 1

# A1  recursion (wed:eq:recursion) reproduces the closed forms (wed:eq:wcoefficients)
w0 = -g
w1 = -P - (sp.diff(w0, U) * P + sp.diff(w0, P) * w0)
w2 = -r * g - (sp.diff(w1, U) * P + sp.diff(w1, P) * w0 + sp.diff(w0, P) * w1)
w3 = -(sp.diff(w2, U) * P + sp.diff(w2, P) * w0 + sp.diff(w1, P) * w1 + sp.diff(w0, P) * w2)
closed = [-g, g * P, -r * g - sp.exp(U) * P**2 + g**2,
          r * sp.exp(U) * P + sp.exp(U) * P**3 - 4 * g * sp.exp(U) * P - g**2 * P]
check('A1 slow-manifold coefficients w0..w3 from the recursion',
      all(sp.simplify(a - b) == 0 for a, b in zip([w0, w1, w2, w3], closed)))
wser = w0 + eps * w1 + eps**2 * w2 + eps**3 * w3
inv = sp.expand(eps * (sp.diff(wser, U) * P + sp.diff(wser, P) * wser) + wser + eps * P + (1 + r * eps**2) * g)
check('A1b invariance equation (wed:eq:invariance) holds through order eps^3',
      all(sp.simplify(inv.coeff(eps, j)) == 0 for j in range(4)))

# A2  modified energy (wed:eq:auxiliary), identity (wed:eq:drift) with phi_3 and phi_4
xi = sp.symbols('xi', positive=True)
um = sp.Function('um')                      # U_-(xi), defined by Rem(um(xi)) = xi
prim = lambda u: sp.exp(u) - u**2 / 2 - u   # a primitive of Rem
Ups = 2 * xi * (U - um(xi)) - 2 * (prim(U) - prim(um(xi)))      # Upsilon(U, xi)
dUps_dU = sp.diff(Ups, U)
dUps_dxi = sp.simplify(sp.diff(Ups, xi).subs(sp.exp(um(xi)), xi + um(xi) + 1))
check('A2a d Upsilon/dU = 2(xi - Rem(U)),  d Upsilon/dxi = 2(U - U_-(xi))',
      sp.simplify(dUps_dU - 2 * (xi - Rem(U))) == 0 and sp.simplify(dUps_dxi - 2 * (U - um(xi))) == 0)
E0 = P**2 / 2 + Rem(U)
Theta = sp.exp(2 * U) / 2 - 2 * sp.exp(U) + U + sp.Rational(3, 2)
check('A2b Theta\' = g^2 and Theta(0) = 0', sp.simplify(sp.diff(Theta, U) - g**2) == 0 and Theta.subs(U, 0) == 0)
E2_U = r * g - g**2 + (dUps_dU + dUps_dxi * g).subs(xi, E0)      # chain rule, d E0/dU = g
E2_P = (dUps_dxi * P).subs(xi, E0)                              # chain rule, d E0/dP = P
check('A2c cohomological identity (wed:eq:cohomological)',
      sp.simplify(P * E2_U - g * E2_P - (r * g * P - g**2 * P + P**3)) == 0)
wt = sp.Symbol('wtilde', real=True)
w = wser + eps**4 * wt
lhs = P * (g + eps**2 * E2_U) + w * (P + eps * P**2 + eps**2 * E2_P)
phi3 = r * P**2 - 4 * g * sp.exp(U) * P**2 + 2 * g * P**2 * (U - um(E0))
phi4 = w3 * P**2 + w2 * E2_P + eps * w3 * E2_P + wt * (P + eps * P**2 + eps**2 * E2_P)
check('A2d modified-energy identity (wed:eq:drift) with phi_3, phi_4 of (wed:eq:auxiliary), (wed:eq:phi4)',
      sp.simplify(sp.expand(lhs - eps**3 * phi3 - eps**4 * phi4)) == 0)
check('A2e d_P E_2 = 2P(U - U_-(E_0))', sp.simplify(E2_P - 2 * P * (U - um(E0))) == 0)

# A3  polynomial identities in the proof of Theorem thm:kurtosis (ii)
om, Pm, s_, gg = sp.symbols('omega P_m s g_ratio', positive=True)
m2 = sp.integrate(P**2, (P, s_, Pm)) / (Pm - s_)
m4 = sp.integrate(P**4, (P, s_, Pm)) / (Pm - s_)
m2 = sp.simplify(m2.subs(s_, om * Pm)); m4 = sp.simplify(m4.subs(s_, om * Pm))
check('A3a m_2, m_4 closed forms',
      sp.simplify(m2 - Pm**2 * (1 + om + om**2) / 3) == 0 and
      sp.simplify(m4 - Pm**4 * (1 + om + om**2 + om**3 + om**4) / 5) == 0)
check('A3b (6/5)Pm^2 m_2 - Pm^4/5 - m_4 = (Pm^4/5) omega (1+omega)(1-omega^2)  (wed:eq:polyidentity)',
      sp.simplify(sp.Rational(6, 5) * Pm**2 * m2 - Pm**4 / 5 - m4 - Pm**4 / 5 * om * (1 + om) * (1 - om**2)) == 0)
check('A3c (9/5)g^2 - (6/5)g + 1/5 = (3g-1)^2/5',
      sp.expand(sp.Rational(9, 5) * gg**2 - sp.Rational(6, 5) * gg + sp.Rational(1, 5) - (3 * gg - 1)**2 / 5) == 0)

# A4  exact linear equation for varpi = Q - w (wed:eq:varpi), w an arbitrary function
wf = sp.Function('w')(U, P)
Qdot = (-Q - eps * P - rho * g) / eps
varpi_dot = Qdot - (sp.diff(wf, U) * P + sp.diff(wf, P) * Q)
inv_res = eps * (sp.diff(wf, U) * P + sp.diff(wf, P) * wf) + wf + eps * P + rho * g
check('A4 eps*varpi\' + (1 + eps w_P)(Q - w) = -(invariance residual)',
      sp.simplify(eps * varpi_dot + (1 + eps * sp.diff(wf, P)) * (Q - wf) + inv_res) == 0)

# A5  energy balance (wed:eq:energybalance)
check('A5 d/dtau(P^2/2 + rho Rem(U)) = -eps P (Qdot + P) along (wed:eq:S1)',
      sp.simplify(P * Q + rho * g * P + eps * P * (Qdot + P)) == 0)

# A6  Lemma lem:tauidentities (f): tau-form  <=>  original-variable form (wed:eq:identitiesoriginal)
a2, a4, b22, c3, d2, e2 = sp.symbols('a2 a4 b22 c3 d2 e2', real=True)   # <p^2>,<p^4>,<p'^2p^2>,<p'^3>,<p''^2>,<p'^2>
A = 1 + mu
eps2 = mu / A**2
tau_form = (rho - 1) * a2 / mu - eps2 * ((a4 / mu**2 - 3 * b22 / mu**3) / rho - c3 / mu**3
                                         + 2 * (1 - 1 / rho) * (d2 / mu**3 - e2 / mu**2))
delta = (rho - 1) * mu * A
orig_form = delta * A * a2 - ((mu * a4 - 3 * b22) / rho - c3 + 2 * (1 - 1 / rho) * (d2 - mu * e2))
check('A6 identity (f) in original variables is (f) times mu^2 A^2', sp.simplify(tau_form * mu**2 * A**2 - orig_form) == 0)

# A7  moment identity (wed:eq:moment)
h = sp.Symbol('etahat', real=True)
fmu = lambda s: s * (s + 1) * (s + 1 - mu)
check('A7 f_mu(mu + mu*etahat) = kH + (A^2+mu) mu etahat + 2A mu^2 etahat^2 + mu^3 etahat^3',
      sp.expand(fmu(mu + mu * h) - (mu * A + (A**2 + mu) * mu * h + 2 * A * mu**2 * h**2 + mu**3 * h**3)) == 0)

# A8  identities of Subsection wed:sec:missing, with the derivative D along the flow of
#     (U,p,y,eta,zeta): U'=p, p'=k-mu p-y, y'=k e^U-y, eta'=zeta+(mu-1)eta-eta^2, zeta'=k e^U-(1+eta)zeta
p, y, eta, zeta = sp.symbols('p y eta zeta', real=True)
x = k * sp.exp(U)
dp = k - mu * p - y
deta = zeta + (mu - 1) * eta - eta**2
dzeta = x - (1 + eta) * zeta


def D(f):
    return (sp.diff(f, U) * p + sp.diff(f, p) * dp + sp.diff(f, y) * (x - y)
            + sp.diff(f, eta) * deta + sp.diff(f, zeta) * dzeta)


Oorth = zeta * p - eta * dp - (x - y)                       # tangent orthogonality (eq:orthogonality), zero on cycles
X1 = (1 + eta) * dp - deta * p - (eta - mu) * (1 + eta) * p + x - k
check('A8a (1+eta)p\' - eta\'p - (eta-mu)(1+eta)p + x - k = -(orthogonality residual)', sp.simplify(X1 + Oorth) == 0)
rho_k = k / (mu * A)
etah = (eta - mu) / mu
ss = (1 + eta) / A
Pt = p / sp.sqrt(mu)
Vs = Pt / ss
Es = Vs**2 / 2 + rho_k * Rem(U)
dE_dtau = D(Es) / sp.sqrt(mu)
rhsE = sp.sqrt(mu) * etah * Vs**2 + rho_k * g * Vs * (ss**3 - 1) / ss**2
# the difference must vanish when the orthogonality residual vanishes: substitute zeta from Oorth = 0
zeta_on_cycle = sp.solve(sp.Eq(Oorth, 0), zeta)[0]
mu_pos = sp.Symbol('mu', positive=True)
diffE = sp.simplify((dE_dtau - rhsE).subs(zeta, zeta_on_cycle))
check('A8b energy identity (wed:eq:energyidentity) on {orthogonality = 0}', sp.simplify(diffE) == 0)
AL = lambda f: D(f) + (1 + eta) * f                         # A*L, L = eps d/dtau + s, since eps d/dtau = D/A
check('A8c zeta = mu(1 + A L etahat)', sp.simplify(mu * (1 + AL(etah)) - zeta) == 0)
check('A8d spike identity (wed:eq:spike): rho e^U - s = A L^2 etahat',
      sp.simplify(AL(AL(etah)) / A - (rho_k * sp.exp(U) - ss)) == 0)

# A9  exact pulse (wed:eq:exactpulse) of x'=px, p'=-y, y'=x-y
tt = sp.Symbol('t', real=True)
xp = 6 * sp.exp(tt) / (1 + sp.exp(tt))**3
pp = -2 + 3 / (1 + sp.exp(tt))
yp = 3 * sp.exp(tt) / (1 + sp.exp(tt))**2
check('A9 explicit pulse solves the pulse system, p(-inf)=1, p(+inf)=-2',
      sp.simplify(sp.diff(xp, tt) - pp * xp) == 0 and sp.simplify(sp.diff(pp, tt) + yp) == 0
      and sp.simplify(sp.diff(yp, tt) - (xp - yp)) == 0
      and sp.limit(pp, tt, -sp.oo) == 1 and sp.limit(pp, tt, sp.oo) == -2
      and sp.limit(xp, tt, -sp.oo) == 0 and sp.limit(xp, tt, sp.oo) == 0)

# A10 quiescent phase (wed:eq:slowreturn)
lam, qm, th = sp.symbols('lambda q_minus vartheta', positive=True)
pslow = lam - (lam + qm) * sp.exp(-th)
Vslow = lam * th - (lam + qm) * (1 - sp.exp(-th))
check('A10a dV/dt~ = p, dp/dt~ = lambda - p, V(0)=0, p(0)=-q_-',
      sp.simplify(sp.diff(Vslow, th) - pslow) == 0 and sp.simplify(sp.diff(pslow, th) - (lam - pslow)) == 0
      and Vslow.subs(th, 0) == 0 and sp.simplify(pslow.subs(th, 0) + qm) == 0)
check('A10b at the return V=0: p = lambda*vartheta - q_-', sp.simplify(pslow - (lam * th - qm) + Vslow) == 0)
tmin = sp.log((lam + qm) / lam)
check('A10c min V = lambda log(1+q_-/lambda) - q_-',
      sp.simplify(Vslow.subs(th, tmin) - (lam * sp.log(1 + qm / lam) - qm)) == 0
      and sp.simplify(pslow.subs(th, tmin)) == 0)
bb = sp.Symbol('b', positive=True)
th_series = 2 * bb - sp.Rational(2, 3) * bb**2 + sp.Rational(4, 9) * bb**3
res = sp.series(th_series - (1 + bb) * (1 - sp.exp(-th_series)), bb, 0, 5).removeO()
check('A10d vartheta = 2b - (2/3)b^2 + (4/9)b^3 + O(b^4) solves vartheta = (1+b)(1-e^{-vartheta}), b = q_-/lambda',
      sp.simplify(res) == 0, 'so lambda*vartheta - q_- = q_- - (2/3) q_-^2/lambda + ...')

# A11 Hopf data of Appendix A in terms of eps = 1/sigma (Remark wed:rem:hopf)
sig = 1 / eps
transv = sig / (2 * (1 + sig**2))
l1 = -sig / ((sig**2 + 1) * (sig**2 + 4))
check('A11 -transversality*(rho-1)/l_1 = (r/2)(1+4 eps^2) with rho-1 = r eps^2',
      sp.simplify(-transv * r * eps**2 / l1 - r * (1 + 4 * eps**2) / 2) == 0)

# A12 d_P w_3 used in Remark wed:rem:scope (c)
check('A12 d_P w_3 = r e^U + 3 e^U P^2 - 4 g e^U - g^2',
      sp.simplify(sp.diff(closed[3], P) - (r * sp.exp(U) + 3 * sp.exp(U) * P**2 - 4 * g * sp.exp(U) - g**2)) == 0)
say(f'   [Part A: {time.time()-t0:.1f} s]')

# =====================================================================================================
# Part B: floating-point spot checks
# =====================================================================================================
say('')
say('=== Part B: floating-point spot checks (double precision, not interval-certified) ===')
t1 = time.time()
RemN = lambda u: np.expm1(u) - u


def Upm(xi_, rho_=1.0):
    """U_-(xi) < 0 < U_+(xi): solutions of rho*Rem(U) = xi."""
    f = lambda u: rho_ * RemN(u) - xi_
    up = brentq(f, 0.0, np.log(xi_ / rho_ + 2.0) + 2.0 + np.sqrt(2 * xi_ / rho_), xtol=1e-15, rtol=1e-15)
    um_ = brentq(f, -(xi_ / rho_ + 1.0) - 1e-12, 0.0, xtol=1e-15, rtol=1e-15)
    return um_, up


def level(xi_, rho_=1.0, r_=0.7):
    """Integrals over the closed orbit P^2/2 + rho*Rem(U) = xi of U'' + rho(e^U-1) = 0 (half orbit, doubled).
    All integrands are even in P."""
    um_, up = Upm(xi_, rho_)

    def f(t, s):
        Uv, Pv = s[0], s[1]
        e = np.exp(Uv); gv = e - 1.0
        Qv = -rho_ * gv
        phi3 = r_ * Pv**2 - 4 * gv * e * Pv**2 + 2 * gv * Pv**2 * (Uv - um_)
        dPw3 = r_ * e + 3 * e * Pv**2 - 4 * gv * e - gv**2
        return [Pv, Qv, Pv**2, Pv**4, e * Pv**2, Qv**2, gv**2, gv * Uv * Pv**2, gv * e * Pv**2, phi3,
                Qv**2 * Pv**2, Qv**3, dPw3, gv * Pv**2]

    def ev(t, s):
        return s[1]
    ev.terminal = True; ev.direction = 1
    s0 = [up, 0.0] + [0.0] * 12
    a = solve_ivp(f, [0, 1e-9], s0, method='DOP853', rtol=1e-12, atol=1e-14)
    b = solve_ivp(f, [1e-9, 60 + 30 * np.sqrt(xi_ + 1)], a.y[:, -1], method='DOP853', rtol=1e-12, atol=1e-14, events=ev)
    assert b.status == 1
    s = b.y_events[0][0]
    names = ['I0', 'I1', 'eP2', 'Q2', 'g2', 'gUP2', 'geP2', 'phi3', 'Q2P2', 'Q3', 'dPw3', 'gP2']
    out = {n: 2 * v for n, v in zip(names, s[2:])}
    out['T0'] = 2 * b.t_events[0][0]
    out['um'] = um_; out['up'] = up
    return out


worst = dict(prop=0.0, avg=0.0, phi3=0.0, kurt_gap=np.inf, i_low=np.inf, ii_first=np.inf, deriv=0.0, exact=0.0, div=0.0)
r_test = 0.7
for xi_ in [0.05, 0.7, 4.0, 25.0, 150.0]:
    d = level(xi_, 1.0, r_test)
    I0, I1, T0 = d['I0'], d['I1'], d['T0']
    worst['prop'] = max(worst['prop'], abs(d['eP2'] / I0 - 1), abs(d['Q2'] / I0 - 1), abs(d['g2'] / I0 - 1))
    worst['avg'] = max(worst['avg'], abs(d['gP2']) / I0, abs(d['gUP2'] / (I1 / 3) - 1), abs(d['geP2'] / (I1 / 3) - 1))
    worst['phi3'] = max(worst['phi3'], abs(d['phi3'] - (r_test * I0 - 2 * I1 / 3)) / (abs(r_test) * I0 + 2 * I1 / 3))
    p2, p4 = I0 / T0, I1 / T0
    Pm2 = 2 * xi_
    worst['i_low'] = min(worst['i_low'], p2 - Pm2 / 3)
    worst['ii_first'] = min(worst['ii_first'], 1.2 * Pm2 * p2 - Pm2**2 / 5 - p4)
    worst['kurt_gap'] = min(worst['kurt_gap'], 1.8 - p4 / p2**2)
    hh = 1e-4 * xi_
    dp_, dm_ = level(xi_ + hh), level(xi_ - hh)
    dI0 = (dp_['I0'] - dm_['I0']) / (2 * hh); dI1 = (dp_['I1'] - dm_['I1']) / (2 * hh)
    Rsel = lambda dd: (2.0 / 3.0) * dd['I1'] / dd['I0']
    dR = (Rsel(dp_) - Rsel(dm_)) / (2 * hh)
    worst['deriv'] = max(worst['deriv'], abs(dI0 / T0 - 1), abs(dI1 / (3 * I0) - 1), abs(dR / ((2.0 / 3.0) * (3 - p4 / p2**2)) - 1))
    # bracket of Lemma (f) at rho = 1 on the centre:  <P^4> - 3<Q^2P^2> - <Q^3> = (2/3)<P^4>
    worst['exact'] = max(worst['exact'], abs((I1 - 3 * d['Q2P2'] - d['Q3']) / (2 * I1 / 3) - 1))
    # divergence route: oint d_P w_3 = r T0 - 2 I0
    worst['div'] = max(worst['div'], abs(d['dPw3'] - (r_test * T0 - 2 * I0)) / (abs(r_test) * T0 + 2 * I0))
check('B1a Prop. wed:prop:melnikov on the centre (rho=1): oint e^U P^2 = oint Q^2 = oint g^2 = oint P^2',
      worst['prop'] < 1e-9, f"max rel. defect {worst['prop']:.1e}")
check('B1b averages (wed:eq:todaaverages)', worst['avg'] < 1e-9, f"max rel. defect {worst['avg']:.1e}")
check('B1c oint phi_3 = r I0 - (2/3) I1 (wed:eq:thirdaverage), r = 0.7', worst['phi3'] < 1e-9, f"max rel. defect {worst['phi3']:.1e}")
check('B1d Theorem thm:kurtosis (i),(ii): <P^2> >= Pm^2/3, first inequality of (ii), kurtosis < 9/5',
      worst['i_low'] > 0 and worst['ii_first'] > 0 and worst['kurt_gap'] > 0,
      f"min slacks {worst['i_low']:.3e}, {worst['ii_first']:.3e}, {worst['kurt_gap']:.3e}")
check('B1e I0\' = T0, I1\' = 3 I0, Rsel\' = (2/3)(3 - kurtosis) by central differences', worst['deriv'] < 1e-6,
      f"max rel. defect {worst['deriv']:.1e}")
check('B1f Remark wed:rem:exactselection: bracket of (f) on the centre equals (2/3)<P^4>', worst['exact'] < 1e-9,
      f"max rel. defect {worst['exact']:.1e}")
check('B1g Remark wed:rem:scope(c): oint d_P w_3 = r T0 - 2 I0', worst['div'] < 1e-9, f"max rel. defect {worst['div']:.1e}")

# general rho: Proposition wed:prop:melnikov
wr = 0.0
for rho_, xi_ in [(2.5, 0.3), (2.5, 6.0), (0.4, 2.0)]:
    d = level(xi_, rho_)
    M1 = d['I0'] - d['Q2']
    wr = max(wr, abs(d['eP2'] / d['I0'] - 1), abs(d['Q2'] / (rho_ * d['I0']) - 1), abs(M1 / ((1 - rho_) * d['I0']) - 1),
             abs(d['I0'] / d['g2'] - rho_) / rho_)
check('B1h Prop. wed:prop:melnikov for rho = 2.5 and 0.4: oint Q^2 = rho oint P^2, M_1 = (1-rho) oint P^2, ratio = rho',
      wr < 1e-9, f"max rel. defect {wr:.1e}")

# Corollary wed:cor:root and Remark wed:rem:hopf
Rs = lambda x_: (2.0 / 3.0) * level(x_)['I1'] / level(x_)['I0']
ok = True; det = []
for r_ in [0.3, 1.0, 3.0, 12.0]:
    xs = brentq(lambda x_: Rs(x_) - r_, 0.5 * r_, 1.5 * r_, xtol=1e-11, rtol=1e-11)
    d = level(xs, 1.0, r_)
    kurt = (d['I1'] / d['T0']) / (d['I0'] / d['T0'])**2
    slope = r_ * d['T0'] - 2 * d['I0']
    ok &= (0.75 * r_ <= xs <= 1.25 * r_) and abs(slope / (-d['I0'] * (2.0 / 3.0) * (3 - kurt)) - 1) < 1e-8 and slope < -0.8 * d['I0'] * (1 - 1e-12)
    det.append(f"xi*({r_})={xs:.6f}")
check('B1i root xi*(r) in [3r/4, 5r/4]; r T0 - 2 I0 = -I0 Rsel\' <= -(4/5) I0 at the root', ok, ', '.join(det))
check('B1j Rsel(xi)/xi -> 1 and T0 -> 2 pi as xi -> 0', abs(Rs(1e-4) / 1e-4 - 1) < 1e-3 and abs(level(1e-4)['T0'] - 2 * np.pi) < 1e-3,
      f"Rsel/xi = {Rs(1e-4)/1e-4:.6f}, T0 = {level(1e-4)['T0']:.6f} at xi = 1e-4")

# elementary bounds used in Subsection wed:sec:missing
ug = np.linspace(-60, 12, 200001)
slack_g = RemN(ug) + np.sqrt(2 * RemN(ug)) - np.abs(np.expm1(ug))
ok = slack_g.min() > -1e-12
mm, hhh = np.meshgrid(np.linspace(1e-3, 1.0, 200), np.concatenate((np.linspace(-0.999, 5, 400), np.geomspace(5, 1e6, 200))))
sv = 1 + mm * hhh / (1 + mm)
ok2 = np.all(np.abs(sv**3 - 1) / sv**2 <= 7 * mm * np.abs(hhh) / (1 + mm) * (1 + 1e-12))
check('B1k |g| <= Rem + sqrt(2 Rem) on [-60,12];  |s^3-1|/s^2 <= 7 mu |etahat|/A for 0<mu<=1, etahat>-1 (grid)', ok and ok2)
say(f'   [B1: {time.time()-t1:.1f} s]')

# ---- B2: true cycles of the full system
t2 = time.time()
from whlib import find_cycle, cycle_profile, kH
worst2 = {}
for mu_, rho_ in [(0.3, 1.5), (1.0, 2.0), (0.05, 2.0)]:
    k_ = rho_ * kH(mu_)
    c = find_cycle(mu_, k_, warm=80)
    pr = cycle_profile(mu_, k_, c, n=120001)
    tg, Uc, pc, yc, xc, etac, zetac = (pr[n_] for n_ in ('t', 'U', 'p', 'y', 'x', 'eta', 'zeta'))
    A_ = 1 + mu_
    e_ = np.sqrt(mu_) / A_
    dpc = k_ - mu_ * pc - yc
    ddpc = -mu_ * dpc - (xc - yc)
    Pc = pc / np.sqrt(mu_); Qc = dpc / mu_; Qd = ddpc / mu_**1.5
    av = lambda f: np.trapezoid(f, tg) / (tg[-1] - tg[0])
    E = np.exp(Uc)
    defects = {
        '(a)': abs(av(E) - 1),
        '(b)': abs(av(Qc**2) / av(Pc**2) - 1),
        '(c)': abs(av(Qc**2) / (rho_ * av(E * Pc**2)) - 1),
        '(d)': abs((rho_ - 1) * av(Pc**2) - (e_ * av(Pc**3) - 2 * e_ * av(Qc**2 * Pc))) / ((rho_ - 1) * av(Pc**2)),
        '(e)': abs(av(E * Pc**3) + (2 * e_ / rho_) * (av(Qd**2) - av(Qc**2))) / max(abs(av(E * Pc**3)), 1e-300),
        '(f)': abs((rho_ - 1) * av(Pc**2) - e_**2 * ((av(Pc**4) - 3 * av(Qc**2 * Pc**2)) / rho_ - av(Qc**3)
                                                    + 2 * (1 - 1 / rho_) * (av(Qd**2) - av(Qc**2)))) / ((rho_ - 1) * av(Pc**2)),
        'orig(b)': abs(av(dpc**2) / (mu_ * av(pc**2)) - 1),
        'orig(c)': abs(av(xc * pc**2) / (kH(mu_) * av(pc**2)) - 1),
        'orig(f)': abs((k_ - kH(mu_)) * A_ * av(pc**2) - ((mu_ * av(pc**4) - 3 * av(dpc**2 * pc**2)) / rho_ - av(dpc**3)
                                                            + 2 * (1 - 1 / rho_) * (av(ddpc**2) - mu_ * av(dpc**2)))) / ((k_ - kH(mu_)) * A_ * av(pc**2)),
    }
    # moment identity (wed:eq:moment)
    eh = (etac - mu_) / mu_
    defects['moment'] = abs((k_ - kH(mu_)) - ((A_**2 + mu_) * mu_ * av(eh) + 2 * A_ * mu_**2 * av(eh**2) + mu_**3 * av(eh**3))) / (k_ - kH(mu_))
    # energy identity (wed:eq:energyidentity): d E/d tau by the chain rule vs the formula
    sc = (1 + etac) / A_
    detac = zetac + (mu_ - 1) * etac - etac**2
    Vc = Pc / sc
    dV_dt = (dpc / np.sqrt(mu_)) / sc - Pc * (detac / A_) / sc**2
    dE_dtau_num = (Vc * dV_dt + rho_ * np.expm1(Uc) * pc) / np.sqrt(mu_)
    dE_formula = np.sqrt(mu_) * eh * Vc**2 + rho_ * np.expm1(Uc) * Vc * (sc**3 - 1) / sc**2
    defects['energy'] = np.max(np.abs(dE_dtau_num - dE_formula)) / np.max(np.abs(dE_formula))
    # log m_t = lambda_c T sign and multiplier relation m_s m_t = exp(-A T)
    defects['mprod'] = abs(np.log(c['m_s'] * c['m_t']) / (-A_ * c['T']) - 1)
    for kk, vv in defects.items():
        worst2[kk] = max(worst2.get(kk, 0.0), float(vv))
tol2 = {'energy': 1e-6, '(e)': 1e-6}
check('B2 Lemma lem:tauidentities (a)-(f), original-variable forms, moment identity, energy identity on 3 computed cycles',
      all(v < tol2.get(kk, 1e-7) for kk, v in worst2.items()),
      'worst defects: ' + ', '.join(f"{kk} {v:.1e}" for kk, v in worst2.items()))
say(f'   [B2: {time.time()-t2:.1f} s]')

# ---- B3: third-order prediction of Theorem thm:wedge at one point (Table wed:tab:selection, alpha = 1, mu = 1e-3)
t3 = time.time()
try:
    from ampcycle import amp_cycle
    mu_, al = 1e-3, 1.0
    c = amp_cycle(mu_, al)
    e_ = np.sqrt(mu_) / (1 + mu_)
    d = level(float(RemN(al)))
    kurt = (d['I1'] / d['T0']) / (d['I0'] / d['T0'])**2
    r_num = (c['rho'] - 1) / e_**2
    ratio_r = r_num / ((2.0 / 3.0) * d['I1'] / d['I0'])
    ratio_m = c['logmt'] / (-e_**3 * d['I0'] * (2.0 / 3.0) * (3 - kurt))
    check('B3 Table wed:tab:selection row alpha=1, mu=1e-3 reproduced (0.99599, 0.99261)',
          abs(ratio_r - 0.99599) < 2e-5 and abs(ratio_m - 0.99261) < 2e-5, f"ratios {ratio_r:.5f}, {ratio_m:.5f}")
except Exception as exc:                                        # pragma: no cover
    check('B3 Table wed:tab:selection row reproduced', False, f'exception {exc!r}')
say(f'   [B3: {time.time()-t3:.1f} s]')

# ---- B4: formal limit map of Subsection wed:sec:relax
t4 = time.time()
V0 = -40.0


def pulse(q_):
    """Pulse map G(q) of x'=px, p'=-y, y'=x-y in the variable v = log x."""
    def f(t, s):
        v, p_, y_ = s
        return [p_, -y_, np.exp(v) - y_]
    e0 = np.exp(V0)
    s0 = [V0, q_ - e0 / (q_ * (1 + q_)), e0 / (1 + q_)]

    def ev(t, s):
        return s[0] - V0
    ev.terminal = True; ev.direction = -1
    a = solve_ivp(f, [0, 1e-3 / max(q_, 1e-3)], s0, method='DOP853', rtol=1e-12, atol=1e-18)
    b = solve_ivp(f, [a.t[-1], 400.0 / q_ + 200.0], a.y[:, -1], method='DOP853', rtol=1e-12, atol=1e-18, events=ev)
    assert b.status == 1
    v, p_, y_ = b.y_events[0][0]
    return -(p_ - y_ - np.exp(v) / abs(p_))


def slow(lam_, qm_):
    f = lambda s: lam_ * s - (lam_ + qm_) * (-np.expm1(-s))
    hi = 2 * qm_ / lam_ + 1.0
    while f(hi) < 0:
        hi *= 2
    th_ = brentq(f, 1e-14, hi, xtol=1e-15, rtol=1e-14)
    return lam_ * th_ - qm_, th_, lam_ * np.log1p(qm_ / lam_) - qm_


G1 = pulse(1.0)
check('B4a pulse map: G(1) = 2 (exact pulse wed:eq:exactpulse)', abs(G1 - 2.0) < 1e-7, f"G(1) = {G1:.9f}")
gsm = [(pulse(q_) - q_) / q_**2 for q_ in (0.01, 0.03, 0.1)]
check('B4b (G(q)-q)/q^2 = 0.671, 0.680, 0.711 at q = 0.01, 0.03, 0.1',
      all(abs(a_ - b_) < 6e-4 for a_, b_ in zip(gsm, (0.67112, 0.68005, 0.71142))), ', '.join(f"{v:.5f}" for v in gsm))
rho_ = 1.5
Fm = lambda q_: slow(rho_, pulse(q_))[0] - q_
qst = brentq(Fm, 0.6, 1.45, xtol=1e-12, rtol=1e-12)
pp_, S_, Vmin_ = slow(rho_, pulse(qst))
dq = 1e-5 * qst
m0 = ((slow(rho_, pulse(qst + dq))[0]) - (slow(rho_, pulse(qst - dq))[0])) / (2 * dq)
check('B4c limit map at rho = 1.5: q* = 1.1502, m_0 = 0.67983, S = 2.4387, lim mu U_min = -1.0336 (Table wed:tab:limitmap)',
      abs(qst - 1.150216) < 2e-5 and abs(m0 - 0.679827) < 2e-5 and abs(S_ - 2.43869) < 2e-5 and abs(Vmin_ + 1.03364) < 2e-5,
      f"q* = {qst:.6f}, m_0 = {m0:.6f}, S = {S_:.5f}, lim mu U_min = {Vmin_:.5f}")
check('B4d matching constant (16 sqrt2/15)(5/4)^{3/2} = 2.108', abs(16 * np.sqrt(2) / 15 * 1.25**1.5 - 2.108) < 5e-4,
      f"{16*np.sqrt(2)/15*1.25**1.5:.5f}")
say(f'   [B4: {time.time()-t4:.1f} s]')


# =====================================================================================================
# Audit additions (cross-model review of the section, 3 October 2026)
# =====================================================================================================
say('')
say('=== Audit additions: exact (A13-A16) and floating point (B5-B6) ===')
t5 = time.time()

# A13 equality case of Theorem thm:kurtosis (Remark wed:rem:sharp): V-shaped potential a_+ U (U>=0), a_- |U| (U<=0).
#     On the side U>0 the velocity runs linearly from P_m to -P_m in the time 2 P_m/a_+, on the other side back in
#     2 P_m/a_-; so every time integral of an even function of P is (1/a_+ + 1/a_-) * int_{-P_m}^{P_m} (.) dP.
ap, am, PmV, Pv_ = sp.symbols('a_plus a_minus P_m P', positive=True)
wgt = 1 / ap + 1 / am
T0V = wgt * sp.integrate(1, (Pv_, -PmV, PmV))
p2V = wgt * sp.integrate(Pv_**2, (Pv_, -PmV, PmV)) / T0V
p4V = wgt * sp.integrate(Pv_**4, (Pv_, -PmV, PmV)) / T0V
xiV = PmV**2 / 2
check('A13 V-shaped potential: <P^2> = Pm^2/3, <P^4> = Pm^4/5; equality in (i) and in both inequalities of (ii); Rsel = (4/5) xi',
      sp.simplify(p2V - PmV**2 / 3) == 0 and sp.simplify(p4V - PmV**4 / 5) == 0
      and sp.simplify(sp.Rational(6, 5) * PmV**2 * p2V - PmV**4 / 5 - p4V) == 0
      and sp.simplify(sp.Rational(9, 5) * p2V**2 - p4V) == 0
      and sp.simplify(sp.Rational(2, 3) * p4V / p2V - sp.Rational(4, 5) * xiV) == 0)

# A14 polar form (wed:eq:polar) of the planar field U' = P, P' = w, with w = varsigma * what, and its consequences
vs, th_ = sp.symbols('varsigma theta', positive=True)
wh = sp.Function('what')(vs, th_)
Up_, Pp_ = vs * sp.cos(th_), vs * sp.sin(th_)
wpol = vs * wh
vs_dot = (Up_ * Pp_ + Pp_ * wpol) / vs                 # (U U' + P P')/varsigma
th_dot = (Up_ * wpol - Pp_ * Pp_) / vs**2              # (U P' - P U')/varsigma^2
ok14 = (sp.simplify(vs_dot - vs * sp.sin(th_) * (sp.cos(th_) + wh)) == 0
        and sp.simplify(th_dot - (wh * sp.cos(th_) - sp.sin(th_)**2)) == 0)
# at eps = 0, w = -g(U):  r_theta = -(U g(U) + P^2)/(U^2 + P^2)
ok14 &= sp.simplify(((U * (-g) - P**2) / (U**2 + P**2)) + (U * g + P**2) / (U**2 + P**2)) == 0
# at theta = 0:  r_theta = what = w(alpha,0)/alpha, so that P' = w(alpha,0) = alpha * r_theta   (wed:eq:crestsign)
ok14 &= sp.simplify((wh * sp.cos(th_) - sp.sin(th_)**2).subs(th_, 0) - wh.subs(th_, 0)) == 0
check('A14 polar form (wed:eq:polar); r_theta at eps=0; w(alpha,0) = alpha*r_theta(alpha,0) (wed:eq:crestsign)', ok14)

# A15 energy drift (wed:eq:energydrift): along U' = P, P' = w,  d E_0/d tau = P (w + g) = P (w - w|_{eps=0})
wgen = sp.Function('w')(U, P)
dE0 = sp.diff(E0, U) * P + sp.diff(E0, P) * wgen
ok15 = sp.simplify(dE0 - P * (wgen + g)) == 0
ok15 &= sp.simplify(sp.expand(P * (wser + g)) - eps * P * (closed[1] + eps * closed[2] + eps**2 * closed[3])) == 0
check('A15 d E_0/d tau = P (w + g) along the planar field; it is eps*P*(w_1 + eps w_2 + eps^2 w_3) for the truncated w', ok15)

# A16 pulse identity used in Subsection wed:sec:relax: d/dt (p p' + p^2/2 + x) = y^2 along x'=px, p'=-y, y'=x-y,
#     hence G(q)^2 - q^2 = 2 int y^2 dt; on the explicit pulse 2 int y^2 dt = 3 = 2^2 - 1^2.
xs, ps, ys = sp.symbols('x p y', real=True)
expr = ps * (-ys) + ps**2 / 2 + xs
dexpr = sp.diff(expr, xs) * (ps * xs) + sp.diff(expr, ps) * (-ys) + sp.diff(expr, ys) * (xs - ys)
sv_ = sp.Symbol('s', positive=True)
int_y2 = sp.integrate((3 * sv_ / (1 + sv_)**2)**2 / sv_, (sv_, 0, sp.oo))        # s = e^t, dt = ds/s
check('A16 d/dt(p p\' + p^2/2 + x) = y^2 along the pulse system; explicit pulse: 2 int y^2 dt = 3 = G(1)^2 - 1',
      sp.simplify(dexpr - ys**2) == 0 and sp.simplify(2 * int_y2 - 3) == 0)


# ---- B5: the velocity density f of the proof of Theorem thm:kurtosis for Phi = Rem: strict monotonicity, the density
#          representation (wed:eq:momentsdensity), and the Stieltjes form of the slack in the first inequality of (ii)
#          that Remark wed:rem:sharp uses for the strictness statement.
from scipy.integrate import quad


def RemAcc(u):
    """Rem(u) = e^u - 1 - u without cancellation near 0."""
    if abs(u) < 1e-3:
        return u * u / 2 * (1 + u / 3 * (1 + u / 4 * (1 + u / 5 * (1 + u / 6))))
    return float(np.expm1(u) - u)


def UpmAcc(xi_):
    f_ = lambda u: RemAcc(u) - xi_
    up = brentq(f_, 0.0, np.log(xi_ + 2.0) + 2.0 + np.sqrt(2 * xi_), xtol=1e-300, rtol=1e-15)
    um_ = brentq(f_, -(xi_ + 1.0) - 1e-12, 0.0, xtol=1e-300, rtol=1e-15)
    return um_, up


def dens(u_, xi_):
    """P, f(P) and f'(P) of (wed:eq:density) for Phi = Rem at P = P_m (1 - u^2), 0 < u <= 1.
    The level xi - P^2/2 = xi u^2 (2 - u^2) is formed without cancellation."""
    P_ = np.sqrt(2 * xi_) * (1 - u_ * u_)
    um_, up = UpmAcc(xi_ * u_ * u_ * (2 - u_ * u_))
    gp, gm = np.expm1(up), -np.expm1(um_)                       # Phi'(U_+) > 0, |Phi'(U_-)| > 0
    return P_, 1 / gm + 1 / gp, P_ * (np.exp(up) / gp**3 + np.exp(um_) / gm**3)


w5 = dict(mom=0.0, slack=0.0, mono=np.inf)
for xi_ in [0.05, 4.0, 150.0]:
    d = level(xi_)
    T0_, p2, p4 = d['T0'], d['I0'] / d['T0'], d['I1'] / d['T0']
    Pm_ = np.sqrt(2 * xi_)
    jac = lambda u_: 2 * Pm_ * u_                                # dP = -2 Pm u du; the substitution removes the endpoint singularity

    def mom(j):
        def F(u_):
            P_, f_, _ = dens(u_, xi_)
            return P_**j * f_ * jac(u_)
        return quad(F, 0, 1, epsabs=0, epsrel=1e-10, limit=200)[0]
    q0, q2, q4 = mom(0), mom(2), mom(4)
    w5['mom'] = max(w5['mom'], abs(2 * q0 / T0_ - 1), abs((2 / T0_) * q2 / p2 - 1), abs((2 / T0_) * q4 / p4 - 1))

    # slack of the first inequality of (ii) as an integral against nu_1 on (0, Pm)
    def integrand(u_):
        s__, _, fp_ = dens(u_, xi_)
        om_ = s__ / Pm_
        return (Pm_ - s__) * om_ * (1 + om_) * (1 - om_**2) * fp_ * jac(u_)
    st = (2 / T0_) * (Pm_**4 / 5) * quad(integrand, 0, 1, epsabs=0, epsrel=1e-10, limit=200)[0]
    direct = 1.2 * Pm_**2 * p2 - Pm_**4 / 5 - p4
    w5['slack'] = max(w5['slack'], abs(st / direct - 1))
    fv = np.array([dens(u_, xi_)[1] for u_ in np.linspace(1.0, 1e-3, 600)])   # increasing P in [0, Pm)
    w5['mono'] = min(w5['mono'], float(np.min(np.diff(fv) / fv[:-1])))
check('B5 density f for Phi = Rem: (wed:eq:momentsdensity) for j = 0, 2, 4; slack of (ii) = (2/T0)(Pm^4/5) int (Pm-s) omega(1+omega)(1-omega^2) f\'(s) ds > 0; f strictly increasing',
      w5['mom'] < 1e-7 and w5['slack'] < 1e-6 and w5['mono'] > 0,
      f"max rel. defects {w5['mom']:.1e} (moments), {w5['slack']:.1e} (slack); min relative increment of f on the grid {w5['mono']:.2e}")

# ---- B6: statements written out in Steps 2 and 5 of the proof of Theorem thm:wedge, on cycles of the full system
#          computed in the wedge (floating point; the theorem itself is conditional on (F)):
#          polar angle of (U,P) strictly decreasing, winding number -1, exactly two zeros of P, crest = maximum of U
#          with dP/dtau = Q < 0 there, and energy excursion max|E_0 - Rem(alpha)| = O(eps)
#          (leading term eps*max|P|^3/3 = eps*(2 Rem(alpha))^{3/2}/3, from E_1 = P^3/3).
try:
    from ampcycle import amp_cycle as _amp
    rows = []; ok6 = True
    for mu_, al in [(1e-2, 1.0), (1e-3, 1.0), (1e-4, 1.0), (1e-3, 2.0)]:
        c = _amp(mu_, al)
        pr = cycle_profile(mu_, c['k'], c, n=40001)
        Uc, pc, yc = pr['U'], pr['p'], pr['y']
        e_ = np.sqrt(mu_) / (1 + mu_)
        Pc = pc / np.sqrt(mu_); Qc = (c['k'] - mu_ * pc - yc) / mu_
        thdot = (Uc * Qc - Pc**2) / (Uc**2 + Pc**2)              # d theta/d tau
        ang = np.unwrap(np.arctan2(Pc, Uc))
        wind = (ang[-1] - ang[0]) / (2 * np.pi)
        sgn = np.sign(Pc[1:-1]); nzero = int(np.sum(sgn[1:] * sgn[:-1] < 0))   # interior sign changes of P
        E0c = 0.5 * Pc**2 + RemN(Uc)
        exc = np.max(np.abs(E0c - RemN(al))) / e_
        lead = (2 * RemN(al))**1.5 / 3
        good = (thdot.max() < 0 and abs(wind + 1) < 1e-6 and nzero == 1 and abs(Uc.max() - al) < 1e-9
                and int(np.argmax(Uc)) in (0, len(Uc) - 1) and Qc[0] < 0 and Uc[np.argmin(np.abs(Pc[10:-10])) + 10] < 0
                and abs(exc / lead - 1) < 0.25)
        ok6 &= bool(good)
        rows.append(f"(mu={mu_:g}, alpha={al:g}): max dtheta/dtau {thdot.max():.3f}, winding {wind:.6f}, Q(crest) {Qc[0]:.3f}, "
                    f"max|E0-Rem(alpha)|/eps {exc:.4f} vs (2Rem)^1.5/3 = {lead:.4f}")
    check('B6 computed wedge cycles: polar angle strictly decreasing, winding number -1, P vanishes twice, crest is the maximum of U with Q < 0, energy excursion O(eps)',
          ok6, '; '.join(rows))
except Exception as exc_:                                       # pragma: no cover
    check('B6 computed wedge cycles', False, f'exception {exc_!r}')
say(f'   [audit additions: {time.time()-t5:.1f} s]')

say('')
say(f'TOTAL TIME {time.time()-t0:.1f} s')
say('FAILURES: ' + (', '.join(FAILS) if FAILS else 'none'))
(HERE.parent / 's11_wedge.check.out').write_text(LOG.getvalue(), encoding='utf-8')
sys.exit(1 if FAILS else 0)
