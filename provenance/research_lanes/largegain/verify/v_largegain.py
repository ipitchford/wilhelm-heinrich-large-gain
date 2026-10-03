#!/usr/bin/env python3
"""Referee script C: independent large-gain cycle computation and test of Theorem 1 and the
quantitative lemma chain (L6, L8-L17) at k >= 1000.

Own hybrid integrator (not the lane's lgcycle.py, not whlib): the spike (U > -UCUT) is integrated
numerically (DOP853, state U,p,y,eta,zeta + quadratures); the trough (U < -UCUT, where
x = k e^U < k e^-60 is dropped) is solved in closed form for (U,p,w) and for the Riccati pair.
The cycle is found by iterating the return map on the section {U = -UCUT, p > 0} (it is
super-attracting).  usage: v_largegain.py [mu k]"""
import sys
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

UCUT = 60.0
PHI2 = 0.7968
FAILS = []


def chk(tag, name, cond, detail=""):
    print(f"   [{'ok  ' if cond else 'FAIL'}] {name}  {detail}")
    if not cond:
        FAILS.append((tag, name, detail))


# ---------- closed-form trough (x == 0) ----------
def trough_fns(mu, k, p0, w0, U0):
    c = 1.0 - mu

    def E(t):   # int_0^t e^{-mu(t-s)} e^{-s} ds
        return t * np.exp(-t) if abs(c) < 1e-9 else (np.exp(-mu * t) - np.exp(-t)) / c

    def IE(t):  # int_0^t E
        if abs(c) < 1e-9:
            return 1 - (1 + t) * np.exp(-t)
        return (-np.expm1(-mu * t) / mu - (-np.expm1(-t))) / c

    em = lambda t: -np.expm1(-mu * t) / mu      # (1-e^{-mu t})/mu
    w = lambda t: -k + (w0 + k) * np.exp(-t)
    p = lambda t: np.exp(-mu * t) * p0 + k * em(t) - (w0 + k) * E(t)
    U = lambda t: U0 + p0 * em(t) + k * (t - em(t)) / mu - (w0 + k) * IE(t)
    return U, p, w


def trough_riccati(mu, eta0, zeta0, t):
    """Exact solution of eta' = zeta+(mu-1)eta-eta^2, zeta' = -(1+eta)zeta (x=0); returns eta, zeta, int_0^t eta."""
    a = eta0 + zeta0 / mu
    b = zeta0 / mu
    c = mu - 1.0
    if abs(c) < 1e-9:
        G = 1 + a * t - b * (-np.expm1(-t))
        return (a - np.exp(-t) * b) / G, zeta0 * np.exp(-t) / G, np.log(G)
    if c < 0:
        g = -np.expm1(c * t) / (-c)
        G = 1 + a * g - b * (-np.expm1(-t))
        return (np.exp(c * t) * a - np.exp(-t) * b) / G, zeta0 * np.exp(-t) / G, np.log(G)
    Gs = np.exp(-c * t) * (1 - b * (-np.expm1(-t))) + a * (-np.expm1(-c * t)) / c
    return (a - np.exp(-(c + 1) * t) * b) / Gs, zeta0 * np.exp(-mu * t) / Gs, c * t + np.log(Gs)


def selftest_trough():
    rng = np.random.default_rng(1)
    worst = 0.0
    for mu in (0.02, 0.3, 1.0, 2.0):
        k = 7.0
        p0, w0, U0 = -3.0, 5.0, -60.0
        eta0, zeta0 = 2.5, 11.0
        f = lambda t, s: [s[1], -mu * s[1] - s[2], -s[2] - k, s[4] + (mu - 1) * s[3] - s[3] ** 2, -(1 + s[3]) * s[4], s[3]]
        sol = solve_ivp(f, [0, 9.0], [U0, p0, w0, eta0, zeta0, 0.0], method="DOP853", rtol=1e-12, atol=1e-14)
        U, p, w = trough_fns(mu, k, p0, w0, U0)
        e, z, I = trough_riccati(mu, eta0, zeta0, 9.0)
        ref = sol.y[:, -1]
        err = max(abs(U(9.0) - ref[0]), abs(p(9.0) - ref[1]), abs(w(9.0) - ref[2]), abs(e - ref[3]), abs(z - ref[4]), abs(I - ref[5]))
        worst = max(worst, err)
    print(f"self-test closed-form trough vs DOP853: max abs err = {worst:.2e}")
    assert worst < 1e-8


# ---------- spike (numerical) ----------
def spike(mu, k, p_in, w_in, eta_in, zeta_in):
    def f(t, s):
        U, p, y, eta, zeta = s[:5]
        x = k * np.exp(U)
        return [p, k - mu * p - y, x - y, zeta + (mu - 1) * eta - eta * eta, x - (1 + eta) * zeta,
                eta, x - k, x]

    def ev(t, s):
        return s[0] + UCUT
    ev.terminal = True
    ev.direction = -1
    s0 = [-UCUT, p_in, w_in + k, eta_in, zeta_in, 0.0, 0.0, 0.0]
    lam = k / mu
    # leave the section first
    sol = solve_ivp(f, [0, 1e-3 * UCUT / lam], s0, method="DOP853", rtol=1e-11, atol=1e-12, dense_output=True)
    t1 = sol.t[-1]
    sol2 = solve_ivp(f, [t1, t1 + 50.0], sol.y[:, -1], method="DOP853", rtol=1e-11, atol=1e-12, events=ev, dense_output=True,
                     first_step=1e-6 * UCUT / lam)
    assert sol2.status == 1, "spike did not return to the section"

    def dense(t):
        t = np.asarray(t, float)
        return np.where(t <= t1, sol.sol(np.minimum(t, t1)), sol2.sol(np.maximum(t, t1)))
    nodes = np.concatenate([sol.t, sol2.t[1:]])
    return dense, nodes, sol2.t[-1], sol2.y[:, -1]


def bracket_root(fun, nodes, direction=0):
    v = np.array([fun(t) for t in nodes])
    out = []
    for i in range(len(nodes) - 1):
        if v[i] == 0:
            out.append(nodes[i]); continue
        if v[i] * v[i + 1] < 0 and (direction == 0 or np.sign(v[i + 1] - v[i]) == direction):
            out.append(brentq(fun, nodes[i], nodes[i + 1], xtol=1e-16, rtol=1e-15))
    return out


def loop(mu, k, st):
    """one return of the section map; st = (p_in, w_in, eta_in, zeta_in).  Returns new st and a record."""
    dense, nodes, t_end, s_end = spike(mu, k, *st)
    p0, w0, eta0, zeta0 = s_end[1], s_end[2] - k, s_end[3], s_end[4]
    assert p0 < 0
    U, p, w = trough_fns(mu, k, p0, w0, -UCUT)
    lam = k / mu
    # trough: zero of p, then return to U=-UCUT
    hi = 1.0
    while p(hi) <= 0:
        hi *= 2
    s0p = brentq(p, 0.0, hi, xtol=1e-14, rtol=1e-15)
    hi = max(2 * s0p, 1.0)
    while U(hi) + UCUT <= 0:
        hi *= 2
    t_r = brentq(lambda t: U(t) + UCUT, s0p, hi, xtol=1e-12, rtol=1e-15)
    eta1, zeta1, Itr = trough_riccati(mu, eta0, zeta0, t_r)
    new = (p(t_r), w(t_r), eta1, zeta1)
    rec = dict(dense=dense, nodes=nodes, t_end=t_end, s_end=s_end, Utr=U, ptr=p, wtr=w, s0p=s0p, t_r=t_r, Itr=Itr,
               p0=p0, w0=w0, eta0=eta0, zeta0=zeta0)
    return new, rec


def Psi(mu, u):
    if abs(mu - 1) < 1e-12:
        return 1 - (1 + u) * np.exp(-u)
    return (1 - np.exp(-mu * u)) / mu - (np.exp(-u) - np.exp(-mu * u)) / (mu - 1)


ci = lambda a: (np.exp(a) - np.exp(a / 2)) * (1 - (1 + a / 2) * np.exp(-a / 2))
cii = lambda a: a * a / 2 - 1 + (1 + a) * np.exp(-a)
cb = lambda a: (1 - (1 + a / 2) * np.exp(-a / 2)) * (1 - np.exp(-a / 2))


def run(mu, k, theorem_checks=True):
    tag = f"mu={mu} k={k}"
    lam = k / mu
    A = 1 + mu
    st = (0.9 * lam, -0.9 * k, 0.0, 0.0)
    hist = []
    for it in range(8):
        st_new, rec = loop(mu, k, st)
        d = abs(st_new[0] - st[0]) / lam + abs(st_new[1] - st[1]) / k
        hist.append(d)
        st_prev = st
        st = st_new
        if d < 1e-12 and it >= 2:
            break
    # final loop evaluated from the converged section state st_prev -> rec is that loop
    print(f"== {tag}  lambda={lam:.6g} rho={k/(mu*A):.6g}  section-map increments: " + " ".join(f"{h:.1e}" for h in hist))
    dense, nodes = rec["dense"], rec["nodes"]
    Uf = lambda t: float(dense(t)[0]); pf = lambda t: float(dense(t)[1]); wf = lambda t: float(dense(t)[2]) - k
    t_up = bracket_root(Uf, nodes, +1)[0]
    t_cr = bracket_root(pf, nodes, -1)[0]
    t_dn = bracket_root(Uf, nodes, -1)[0]
    t_end = rec["t_end"]
    t_r, s0p = rec["t_r"], rec["s0p"]
    # event times in the lemma's clock (origin = down-crossing of U=0)
    fall = t_end - t_dn                       # down-crossing -> trough entry
    s0 = fall + s0p
    w0 = rec["w0"]
    if w0 > 0:
        sw = fall + np.log((w0 + k) / k)
    else:
        sw = bracket_root(wf, nodes[nodes >= t_dn], -1)[0] - t_dn
    Tm = fall + t_r + t_up
    Tp = t_dn - t_up
    T = Tm + Tp
    tau1 = t_cr - t_up
    D = -rec["Utr"](s0p)
    h = Uf(t_cr)
    pk = pf(t_up); wk = wf(t_up); Pd = -pf(t_dn); wd = wf(t_dn); Wc = wf(t_cr)
    Mplus = float(dense(t_dn)[6] - dense(t_up)[6])
    Ieta = float(rec["s_end"][5]) + rec["Itr"]
    avg_eta = Ieta / T
    print(f"   T={T:.6g}  T*mu/lambda={T*mu/lam:.5f}  T_-={Tm:.6g} T_+={Tp:.5g}  p_k/lambda={pk/lam:.6f}  D={D:.4g} h={h:.5g}")
    print(f"   s_w={sw:.5g} s_0={s0:.5g} tau_1={tau1:.5g}  M_+/p_k^2={Mplus/pk**2:.4f}  W_c={Wc:.5g}  int eta={Ieta:.5g} (spike {float(rec['s_end'][5]):.4g} + trough {rec['Itr']:.4g})")
    print(f"   <eta>={avg_eta:.6e}  <eta>/mu={avg_eta/mu:.3e}  log m_t=(<eta>-mu)T={(avg_eta-mu)*T:.6g}")
    # sanity: identities valid on an exact cycle
    chk(tag, "sanity L5a: M_+ = w_d+P_d+|w_k|+p_k (cycle closure)", abs(Mplus - (wd + Pd + abs(wk) + pk)) < 1e-6 * Mplus,
        f"rel={abs(Mplus-(wd+Pd+abs(wk)+pk))/Mplus:.1e}")
    intx = float(rec["s_end"][7])  # int x over the spike part; trough part is 0
    chk(tag, "sanity F6: <x> = k", abs(intx / T - k) < 1e-6 * k, f"<x>/k={intx/T/k:.9f}")
    # ---- L6
    r = np.roots([1, A, mu, k])
    omega = max(r.imag); gE = -min(r.real[np.abs(r.imag) < 1e-9])
    etaE = gE - 1; zetaE = etaE * (gE - mu)
    chk(tag, "L6: T_+ < pi/omega", Tp < np.pi / omega, f"T_+={Tp:.5g} pi/omega={np.pi/omega:.5g}")
    if k >= 1:
        om0 = 0.5 * np.sqrt(3 * k ** (2 / 3) - 4 * k ** (1 / 3) - 1)
        chk(tag, "L6: T_- > T_0(k)", Tm > 2 * np.pi - np.pi / om0)
    chk(tag, "F4: T > 2pi/sqrt(mu)", T > 2 * np.pi / np.sqrt(mu))
    # Lambda_E > 0 on the spike part and on the trough part
    tt = np.linspace(0, t_end, 20001)
    Ys = dense(tt)
    LamE = (Ys[2] - k) - zetaE * Ys[0] + etaE * Ys[1]
    tq = np.concatenate([np.linspace(0, min(50.0, t_r), 5001), np.geomspace(min(50.0, t_r), t_r, 2001)])
    LamT = rec["wtr"](tq) - zetaE * rec["Utr"](tq) + etaE * rec["ptr"](tq)
    chk(tag, "F2: Lambda_E > 0 (spike + trough samples)", LamE.min() > 0 and LamT.min() > 0, f"min={min(LamE.min(), LamT.min()):.3e}")
    M0 = (lam + mu) * (lam + 1) / mu
    H = np.log(M0)
    # p<lambda and y>0 hold with exponentially small margins lambda*e^{-mu T}, k*e^{-T}: below double precision, so test non-strictly
    chk(tag, "F5: h <= H = log M0, p<=lambda, y>=-1e-10 (strict margins ~e^{-mu T} are below double precision; y starts at underflowed 0)", h <= H and Ys[1].max() <= lam and Ys[2].min() >= -1e-10,
        f"h={h:.5g} H={H:.5g} lam-max p={lam-Ys[1].max():.3e} min y={Ys[2].min():.3e}")
    # ---- L7 signs
    chk(tag, "L7: w_d>0, w_k<0, W_c>0, w(trough)<0, s_w<s_0", wd > 0 and wk < 0 and Wc > 0 and rec["wtr"](s0p) < 0 and sw < s0)
    # ---- L8, L9 with the actual durations
    a = sw
    if a / 2 <= fall:
        v = -Uf(t_dn + a / 2)
    else:
        v = -rec["Utr"](a / 2 - fall)
    chk(tag, "L8: |U(s_w/2)| >= k c_i(s_w)(1-e^-v)", v >= k * ci(a) * (1 - np.exp(-v)), f"v={v:.5g} rhs={k*ci(a)*(1-np.exp(-v)):.5g}")
    if k * ci(a) >= 2:
        chk(tag, "L8: D >= phi2 k c_i(s_w)", D >= PHI2 * k * ci(a), f"D={D:.5g} rhs={PHI2*k*ci(a):.5g}")
    u = s0 - sw
    if w0 > 0:
        tw = np.log((w0 + k) / k)
        vw = -rec["Utr"](tw); Pw = -rec["ptr"](tw)
    else:
        vw = -Uf(t_dn + sw); Pw = -pf(t_dn + sw)
    phw = 1 - np.exp(-vw)
    chk(tag, "L9: P_w >= k phi(s_w)(u-1+e^-u)", Pw >= k * phw * (u - 1 + np.exp(-u)), f"P_w={Pw:.5g} rhs={k*phw*(u-1+np.exp(-u)):.5g}")
    chk(tag, "L9: D - v_w >= k phi(s_w) c_ii(u)", D - vw >= k * phw * cii(u), f"{D-vw:.5g} >= {k*phw*cii(u):.5g}")
    chk(tag, "L9: P_w < (gamma_E - mu) v_w", Pw < (gE - mu) * vw, f"ratio={Pw/((gE-mu)*vw):.4g}")
    if k * (u - 1 + np.exp(-u)) / (k ** (1 / 3) + 1) >= 2:
        chk(tag, "L9: D >= phi2 k c_ii(u)", D >= PHI2 * k * cii(u), f"D={D:.5g} rhs={PHI2*k*cii(u):.5g}")
    # ---- L10
    ell = 2.0; phl = 1 - np.exp(-ell)
    t_l = bracket_root(lambda t: Uf(t) + ell, nodes[nodes <= t_up], +1)[0]
    taub = t_up - t_l
    r_l = (t_r - s0p) + t_l
    psl = pf(t_l)
    chk(tag, "L10a: T_- - s_0 >= p_k/k ; D >= p_k^2/(2k)", Tm - s0 >= pk / k and D >= pk ** 2 / (2 * k), f"{Tm-s0:.5g}>={pk/k:.5g}; {D:.5g}>={pk**2/(2*k):.5g}")
    chk(tag, "L10b: r_l >= sqrt(2(D-l)/k)", r_l >= np.sqrt(2 * (D - ell) / k), f"{r_l:.6g} >= {np.sqrt(2*(D-ell)/k):.6g}")
    chk(tag, "L10c: p(s_l) >= k phi_l Psi_mu(r_l); p_k >= p(s_l) - mu l", psl >= k * phl * Psi(mu, r_l) and pk >= psl - mu * ell, f"p(s_l)/lam={psl/lam:.6f} rhs/lam={k*phl*Psi(mu,r_l)/lam:.6f}")
    chk(tag, "L10d: tau_b < 0.36 (and k c_b(tau_b) < l/phi_l)", taub < 0.36 and k * cb(taub) < ell / phl, f"tau_b={taub:.4g}")
    if not theorem_checks:
        return dict(mu=mu, k=k, T=T, avg_eta=avg_eta, pk=pk)
    # ---- L11
    chk(tag, "L11: p_k >= 0.374 k", pk >= 0.374 * k, f"p_k/k={pk/k:.5g}")
    print(f"   Stage-1 pigeonhole: s_w={sw:.4g} (a1=2.06)  u={u:.4g} (a2=1.99)  r_l={r_l:.6g} (a3=1.48)  tau_b={taub:.4g} (a4=0.36)")
    chk(tag, "L11: at least one of s_w>=a1, u>=a2, r_l>=a3", sw >= 2.06 or u >= 1.99 or r_l >= 1.48)
    # ---- L12
    tg = np.linspace(t_up, t_cr, 40001)
    Yg = dense(tg)
    P1 = Yg[1].max()
    eps1 = 1.1 * np.sqrt(2) * k * pk ** -1.5
    chk(tag, "L12: P_1 <= p_k(1+eps_1), P_1 < lambda", P1 <= pk * (1 + eps1) and P1 < lam, f"P1/p_k={P1/pk:.6f} 1+eps1={1+eps1:.5f}")
    # ---- L13
    th = 2 ** -0.5
    i = np.nonzero(Yg[1] <= th * pk)[0][0]
    t_th = brentq(lambda t: pf(t) - th * pk, tg[i - 1], tg[i], xtol=1e-16)
    sth = t_th - t_up
    Uth = Uf(t_th)
    Eth = th * pk * (0.5 * pk ** 2 - 2 * mu * P1 * H) / (2 * k)
    chk(tag, "L13: e^{U_theta} >= E_theta", np.exp(Uth) >= Eth, f"U_th={Uth:.5g} log E_th={np.log(max(Eth,1e-300)):.5g}")
    chk(tag, "L13: tau_1 - s_theta >= theta p_k/(mu P_1 + M_+)", tau1 - sth >= th * pk / (mu * P1 + Mplus), f"{tau1-sth:.4g} >= {th*pk/(mu*P1+Mplus):.4g}")
    Q = 0.25 * pk ** 2 * (0.5 * pk ** 2 - 2 * mu * P1 * H) - pk * k / np.sqrt(2)
    chk(tag, "L13: M_+(M_+ + mu P_1) >= Q", Mplus * (Mplus + mu * P1) >= Q, f"lhs/Q={Mplus*(Mplus+mu*P1)/Q:.4g}")
    Ht = 2 * np.log(k + 1) + 3 / np.e
    dQ = 4 * (1 + eps1) * Ht / pk + 4 * np.sqrt(2) * k / pk ** 3
    chk(tag, "L13: T_- > M_+/k >= p_k^2 sqrt(1-dQ)/(sqrt8 k) - 1/2", Tm > Mplus / k >= pk ** 2 * np.sqrt(1 - dQ) / (np.sqrt(8) * k) - 0.5,
        f"T_-={Tm:.6g} M_+/k={Mplus/k:.6g} rhs={pk**2*np.sqrt(1-dQ)/(np.sqrt(8)*k)-0.5:.6g}")
    # ---- L14
    chk(tag, "L14: T_- >= 22 and r_l >= 0.3 T_-", Tm >= 22 and r_l >= 0.3 * Tm, f"r_l/T_-={r_l/Tm:.5f}")
    z = 0.3 * Tm - 3
    chk(tag, "L14: p_k >= 0.95 phi_l k (1-e^{-mu z})/mu - l", pk >= 0.95 * phl * k * (1 - np.exp(-mu * z)) / mu - ell, f"p_k/lam={pk/lam:.5f} rhs/lam={(0.95*phl*k*(1-np.exp(-mu*z))/mu-ell)/lam:.5f}")
    # ---- L15
    d1K0 = 0.19401605
    chk(tag, "L15 (Tm): T_- > p_k^2 sqrt(1-delta1(K0))/(sqrt8 k) - 1/2", Tm > pk ** 2 * np.sqrt(1 - d1K0) / (np.sqrt(8) * k) - 0.5)
    kap0 = (0.94 * phl) ** 2 * np.sqrt(1 - d1K0) / np.sqrt(8)
    hz = (z + 3.15) * (1 + mu * z) ** 2 / z ** 2
    zs = 0.055 * k / mu ** 2
    chk(tag, "L15: h(z) > 0.3 kappa0 k and z > z#", hz > 0.3 * kap0 * k and z > zs, f"h(z)/(0.3 kappa0 k)={hz/(0.3*kap0*k):.4g}  z/z#={z/zs:.4g}")
    # ---- Theorem (i),(ii)
    chk(tag, "Thm(i): T >= 0.1833 lambda/mu", T >= 0.1833 * lam / mu, f"T mu/lambda={T*mu/lam:.5f}")
    chk(tag, "Thm(ii): p_k >= 0.819 lambda", pk >= 0.819 * lam, f"p_k/lambda={pk/lam:.6f}")
    # ---- L16 / Thm(iii)
    nu = mu / 2
    c0 = (mu - nu) * (1 - nu) * (A - nu)
    L0 = np.log(k / c0)
    K0 = 1000.0; cp = 0.819; b = 1.02
    Lh = np.log(8 * lam / 3); Hh = 3 * np.log(lam) - np.log(K0) + 2 * np.log(1 + 1 / K0); Hl = 2 * np.log(lam + 1) + 3 / np.e
    e1h = 1.1 * np.sqrt(2) * cp ** -1.5 / np.sqrt(lam)
    dp = 4 * (1 + e1h) * Hl / (cp * lam) + 4 / (th * cp ** 3 * lam ** 2)
    ac = np.sqrt(12 / (1 - dp)) + 12 / (th * cp ** 2 * lam * (1 - dp))
    A1 = np.sqrt(2) * Hh + ac
    tb1 = A1 / (cp * lam)
    What = cp ** 2 * lam ** 2 * np.exp(-tb1) / A1
    vf = b * np.sqrt(2 * (Hh + Lh) / What)
    t_L0u = bracket_root(lambda t: Uf(t) + L0, nodes[nodes <= t_up], +1)[0]
    t_L0d = bracket_root(lambda t: Uf(t) + L0, nodes[nodes >= t_dn], -1)[0]
    tauc = t_L0d - t_L0u
    chk(tag, "L16: D > L0, L0 <= Lhat, H <= Hhat", D > L0 and L0 <= Lh + 1e-12 and H <= Hh + 1e-12)
    chk(tag, "L16 rise: T_- - t_{L0} <= 2 L0/p_k <= 2Lhat/(c_p lambda)", t_up - t_L0u <= 2 * L0 / pk <= 2 * Lh / (cp * lam), f"{t_up-t_L0u:.4g} <= {2*L0/pk:.4g}")
    Xth = k * (np.exp(Uth) - 1)
    chk(tag, "L16: X_theta >= theta p_k^3 (1-delta')/4 ; s_theta <= Hhat/(theta p_k)", Xth >= 0.25 * th * pk ** 3 * (1 - dp) and sth <= Hh / (th * pk),
        f"X_th/(theta p_k^3/4)={Xth/(0.25*th*pk**3):.4g}  p_k s_th={pk*sth:.4g} (<= {Hh/th:.4g})")
    chk(tag, "L16: crest at most a_c/p_k after s_theta ; p_k tau_1 <= A_1", (tau1 - sth) <= ac / pk and pk * tau1 <= A1, f"p_k(tau1-s_th)={pk*(tau1-sth):.4g} (a_c={ac:.4g});  p_k tau_1={pk*tau1:.4g} (A_1={A1:.4g})")
    chk(tag, "L16: W_c >= p_k^2 e^{-tau1bar}/A_1 >= What", Wc >= pk ** 2 * np.exp(-tb1) / A1 >= What, f"W_c={Wc:.5g} What={What:.5g}")
    chk(tag, "L16: time crest -> U=-L0 <= v_f", t_L0d - t_cr <= vf, f"{t_L0d-t_cr:.4g} <= {vf:.4g}")
    # post-crest pointwise inequality h-U >= (W_c+k)(1-(1+v)e^-v) - k(v-1+e^-v), v <= v_+
    vv = np.linspace(0, min(t_end - t_cr, np.log(1 + Wc / k)), 4001)
    drop = h - dense(t_cr + vv)[0]
    # cancellation-free evaluation: 1-(1+v)e^-v = -expm1(-v) - v e^-v ; v-1+e^-v = v + expm1(-v)
    lowb = (Wc + k) * (-np.expm1(-vv) - vv * np.exp(-vv)) - k * (vv + np.expm1(-vv))
    gap = drop - lowb
    j = int(np.argmin(gap))
    tol = 1e-9 * np.maximum(1.0, np.abs(lowb)) + 1e-10 * abs(h)   # relative 1e-9 (integration rtol 1e-11) + abs noise in U
    chk(tag, "L16: post-crest drop >= (W_c+k)(1-(1+v)e^-v) - k(v-1+e^-v)", np.all(gap >= -tol),
        f"min(drop-lowb)={gap[j]:.3e} at v={vv[j]:.3e} (drop={drop[j]:.6e}, lowb={lowb[j]:.6e}); relative min={np.min(gap[1:]/np.maximum(lowb[1:],1e-300)):.3e}")
    Cbar = 71.79 + 12.28 * np.log(lam / K0)
    chk(tag, "Thm(iii): tau_c <= Cbar_tau(lambda)/(0.819 lambda)", tauc <= Cbar / (cp * lam), f"lambda*tau_c={lam*tauc:.5g}  bound={Cbar/cp:.5g}")
    # ---- L17 and Thm(iv)
    s_ = 2 * np.sqrt(lam)
    al = 1 / (1 - nu); be = al / (A - nu)
    rhs = -nu * (T - tauc) + (s_ - mu) * tauc + k * T / s_ ** 2 + np.log(1 + al * s_ + be * s_ ** 2)
    chk(tag, "L17: (<eta>-mu)T <= -nu(T-tau_c)+(s-mu)tau_c+kT/s^2+log(1+alpha s+beta s^2)", (avg_eta - mu) * T <= rhs, f"lhs={(avg_eta-mu)*T:.6g} rhs={rhs:.6g}")
    worst = np.inf
    for nu2 in (0.1 * mu, 0.5 * mu, 0.9 * mu):
        c02 = (mu - nu2) * (1 - nu2) * (A - nu2)
        L02 = np.log(k / c02)
        tu = bracket_root(lambda t: Uf(t) + L02, nodes[nodes <= t_up], +1)[0]
        td = bracket_root(lambda t: Uf(t) + L02, nodes[nodes >= t_dn], -1)[0]
        tc2 = td - tu
        al2 = 1 / (1 - nu2); be2 = al2 / (A - nu2)
        for s2 in (np.sqrt(2), 5.0, np.sqrt(lam), 2 * np.sqrt(lam), lam ** (2 / 3), lam):
            rhs2 = -nu2 * (T - tc2) + (s2 - mu) * tc2 + k * T / s2 ** 2 + np.log(1 + al2 * s2 + be2 * s2 ** 2)
            worst = min(worst, rhs2 - (avg_eta - mu) * T)
    chk(tag, "L17: holds for 3 nu x 6 s", worst >= 0, f"min slack={worst:.5g}")
    chk(tag, "Thm(iv): <eta> <= 0.83 mu", avg_eta <= 0.83 * mu, f"<eta>/mu={avg_eta/mu:.4e}")
    return dict(mu=mu, k=k, T=T, avg_eta=avg_eta, pk=pk)


if __name__ == "__main__":
    selftest_trough()
    if len(sys.argv) == 3:
        run(float(sys.argv[1]), float(sys.argv[2]))
    else:
        for mu in (1.0, 0.3, 0.02):
            for k in (1000.0, 3000.0):
                run(mu, k)
    print()
    print("FAILS:", FAILS if FAILS else "none")
