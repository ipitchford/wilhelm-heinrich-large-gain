#!/usr/bin/env python3
"""Referee script B: test every parameter-free lemma of theorem.tex (L5-L10, L12, L13, L17 and the
imports F2-F6) on cycles at MODERATE gain, mu in {0.02,0.3,1}, rho in {1.01,1.5,4,12}.
These cycles are outside the theorem's range (k<1000) but the lemmas are claimed for every cycle
with k>k_H, 0<mu<=1.  Uses whlib.find_cycle for the cycle and an own dense integration + event
finder; nothing from the lane's scripts.
usage: v_moderate.py [mu rho]"""
import sys
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
sys.path.insert(0, "/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03/work")
from whlib import find_cycle, equilibrium_data

PHI2 = 0.7968
FAILS = []


def chk(tag, name, cond, detail=""):
    print(f"   [{'ok  ' if cond else 'FAIL'}] {name}  {detail}")
    if not cond:
        FAILS.append((tag, name, detail))


def Psi(mu, u):
    if abs(mu - 1) < 1e-12:
        return 1 - (1 + u) * np.exp(-u)
    return (1 - np.exp(-mu * u)) / mu - (np.exp(-u) - np.exp(-mu * u)) / (mu - 1)


ci = lambda a: (np.exp(a) - np.exp(a / 2)) * (1 - (1 + a / 2) * np.exp(-a / 2))
cii = lambda a: a * a / 2 - 1 + (1 + a) * np.exp(-a)
cb = lambda a: (1 - (1 + a / 2) * np.exp(-a / 2)) * (1 - np.exp(-a / 2))


def rhs(mu, k):
    def f(t, s):
        U, p, y, eta, zeta = s[:5]
        x = k * np.exp(U)
        dp = k - mu * p - y
        return [p, dp, x - y,
                zeta + (mu - 1) * eta - eta * eta, x - (1 + eta) * zeta,
                eta, x - k, p * p, dp * dp, x * p * p]
    return f


def roots_on(fun, a, b, n=200001):
    t = np.linspace(a, b, n)
    v = fun(t)
    out = []
    idx = np.nonzero(np.sign(v[:-1]) * np.sign(v[1:]) < 0)[0]
    for i in idx:
        out.append(brentq(lambda tt: float(fun(np.array([tt]))[0]), t[i], t[i + 1], xtol=1e-14, rtol=1e-14))
    return out


def analyse(mu, k, cyc):
    tag = f"mu={mu} k={k:.6g} rho={k/(mu*(1+mu)):.4g}"
    A = 1 + mu
    lam = k / mu
    kH = mu * A
    T = cyc["T"]
    ed = equilibrium_data(mu, k)
    gE, etaE, zetaE = ed["gammaE"], ed["etaE"], ed["zetaE"]
    omega = abs(ed["lam_u"].imag)
    f = rhs(mu, k)
    s = np.array([cyc["U0"], 0.0, cyc["y0"], etaE, zetaE, 0, 0, 0, 0, 0], float)
    for _ in range(4):  # settle the Riccati pair
        so = solve_ivp(f, [0, T], s, method="DOP853", rtol=1e-11, atol=1e-13)
        s = so.y[:, -1].copy()
        s[0], s[1], s[2] = cyc["U0"], 0.0, cyc["y0"]
        s[5:] = 0
    sol = solve_ivp(f, [0, 2 * T], s, method="DOP853", rtol=1e-12, atol=1e-14, dense_output=True)
    _raw = sol.sol

    def D_(t):  # dense state with third component converted to w = y - k
        Y = np.array(_raw(t), float)
        Y[2] = Y[2] - k
        return Y
    comp = lambda i: (lambda t: D_(t)[i])
    U, P, W = comp(0), comp(1), comp(2)
    # periodicity defect
    d0, d1 = D_(0.0), D_(T)
    print(f"== {tag}: T={T:.6f} periodicity defect (U,p,w)={np.abs(d1[:3]-d0[:3]).max():.2e}  m_t={cyc['m_t']:.3e}")
    # events: crest at 0; down-crossing t_d; trough; up-crossing; crest at T
    zp = roots_on(P, 1e-9, T - 1e-9)
    chk(tag, "F3: p has exactly one zero strictly inside (crest,crest+T) [trough]", len(zp) == 1, f"zeros={len(zp)}")
    t_tr = zp[0]
    zU = roots_on(U, 0, T)
    chk(tag, "F3: U crosses 0 exactly twice per period", len(zU) == 2, f"n={len(zU)}")
    t_d, t_u = zU[0], zU[1]
    # shift: time origin at t_d
    Tm = t_u - t_d
    Tp = T - Tm
    s0 = t_tr - t_d
    tau1 = T - t_u
    at = lambda tt: D_(t_d + tt)  # state at shifted time
    Pd = -at(0)[1]; wd = at(0)[2]; pk = at(Tm)[1]; wk = at(Tm)[2]
    Dd = -at(s0)[0]; h = d0[0]
    Wc = d0[2]
    tt = np.linspace(0, T, 400001)
    Y = D_(t_d + tt)
    Uu, pp, ww = Y[0], Y[1], Y[2]
    yy = _raw(t_d + tt)[2]
    ynode = sol.y[2]            # y at accepted solver nodes (no dense-interpolation error)
    hmax = np.diff(sol.t).max()
    xx = k * np.exp(Uu)
    print(f"   T_-={Tm:.6f} T_+={Tp:.6f} s0={s0:.6f} D={Dd:.6g} h={h:.6g} p_k={pk:.6g} P_d={Pd:.6g} w_d={wd:.6g} w_k={wk:.6g}")
    # ---- imports
    LamE = ww - zetaE * Uu + etaE * pp
    chk(tag, "F2: Lambda_E > 0", LamE.min() > 0, f"min={LamE.min():.3e}")
    inCp = np.any((Uu >= 0) & (pp <= 0) & (ww <= 0)); inCm = np.any((Uu <= 0) & (pp >= 0) & (ww >= 0))
    chk(tag, "F2: no point in C_+ or C_-", not inCp and not inCm)
    chk(tag, "F4: T > 2pi/sqrt(mu)", T > 2 * np.pi / np.sqrt(mu), f"{T:.5f} vs {2*np.pi/np.sqrt(mu):.5f}")
    M0 = (lam + mu) * (lam + 1) / mu
    chk(tag, "F5: x<=k*M0, y>0 (to integrator noise 1e-9k; exact by y=int e^-(t-s)x ds), p<lambda", xx.max() <= k * M0 and ynode.min() > -1e-9 * k and pp.max() <= lam * (1 + 1e-12),
        f"max x/k={xx.max()/k:.4g} M0={M0:.4g} min y/k at nodes={ynode.min()/k:.3e} (dense interp: {yy.min()/k:.3e}, max step {hmax:.3g}) lam-max p={lam-pp.max():.3e} (strict margin ~lam e^(-mu T_-) may be below double precision)")
    Q = d1[5:] - d0[5:]
    chk(tag, "F6: <x>=k, <p'^2>=mu<p^2>", abs(Q[1]) < 1e-7 * k * T and abs(Q[3] - mu * Q[2]) < 1e-7 * abs(Q[3]),
        f"int(x-k)={Q[1]:.2e} rel={abs(Q[3]-mu*Q[2])/abs(Q[3]):.1e}")
    avg_eta = Q[0] / T
    lam_c_M = np.log(cyc["m_t"]) / T if cyc["m_t"] > 0 else np.nan
    print(f"   <eta>={avg_eta:.8f}  <eta>-mu={avg_eta-mu:.6e}  log(m_t)/T={lam_c_M:.6e}")
    # ---- L5
    qd = lambda t: D_(t_d + t)[6]
    Mplus = qd(T) - qd(Tm)
    Mneg = -(qd(Tm) - qd(0))
    chk(tag, "L5a: M_+ = int_{neg}(k-x) = w_d+P_d+|w_k|+p_k", abs(Mplus - Mneg) < 1e-7 * max(1, Mplus) and abs(Mplus - (wd + Pd + abs(wk) + pk)) < 1e-7 * max(1, Mplus),
        f"M_+={Mplus:.8g} sum={wd+Pd+abs(wk)+pk:.8g}")
    chk(tag, "L5a: M_+ < k T_-", Mplus < k * Tm, f"{Mplus:.6g} < {k*Tm:.6g}")
    chk(tag, "L5b: <(x-kH)p^2>=0 and D>log rho", abs(Q[4] - kH * Q[2]) < 1e-7 * abs(Q[4]) and Dd > np.log(k / kH),
        f"rel={abs(Q[4]-kH*Q[2])/abs(Q[4]):.1e} D={Dd:.5g} log rho={np.log(k/kH):.5g}")
    # ---- L6
    chk(tag, "L6: T_+ < pi/omega", Tp < np.pi / omega, f"T_+={Tp:.6f} pi/omega={np.pi/omega:.6f} ratio={Tp*omega/np.pi:.5f}")
    rad = 3 * k ** (2 / 3) - 4 * k ** (1 / 3) - 1
    if k >= 1 and rad > 0:
        om0 = 0.5 * np.sqrt(rad)
        chk(tag, "L6: omega>=omega0(k), T_->2pi-pi/omega0", omega >= om0 and Tm > 2 * np.pi - np.pi / om0, f"omega={omega:.4f} om0={om0:.4f} T_-={Tm:.4f}")
    chk(tag, "L6: gamma_E - mu < k^(1/3)+1", gE - mu < k ** (1 / 3) + 1)
    # ---- L7
    wtr = at(s0)[2]
    chk(tag, "L7a: w_d>0, w_k<0, w(trough)<0, w(crest)>0", wd > 0 and wk < 0 and wtr < 0 and Wc > 0)
    zw = roots_on(lambda t: D_(t_d + t)[2], 0, Tm)
    chk(tag, "L7b: unique zero s_w of w on the negative lobe, s_w<s_0", len(zw) == 1 and zw[0] < s0, f"n={len(zw)}")
    sw = zw[0]
    # (e) w(t+v)+k >= e^-v (w(t)+k): equivalent to e^t (w+k) nondecreasing
    # (e^t y)' = e^t x >= 0 is an exact identity; numerically test it at solver nodes where y is resolved
    tn = sol.t; okn = ynode > 1e-6 * k
    gn = np.log(np.where(okn, ynode, 1.0)) + tn
    dgn = np.diff(gn)[okn[1:] & okn[:-1]]
    chk(tag, "L7e: e^t(w+k) nondecreasing (solver nodes, y>1e-6k)", np.all(dgn > -1e-7), f"min increment={dgn.min():.2e}")
    pos = tt >= Tm
    g2 = np.exp(tt[pos] - Tm) * ww[pos]
    chk(tag, "L7d: e^t w nondecreasing on positive lobe", np.all(np.diff(g2) >= -1e-9 * np.abs(g2).max()))
    # ---- L8 (phase i), for a = s_w and smaller a
    for frac in (1.0, 0.7, 0.4):
        a = frac * sw
        v = -at(a / 2)[0]
        chk(tag, f"L8: v=|U(a/2)| >= k c_i(a)(1-e^-v), a={frac}*s_w", v >= k * ci(a) * (1 - np.exp(-v)) * (1 - 1e-9),
            f"a={a:.4f} v={v:.5g} rhs={k*ci(a)*(1-np.exp(-v)):.5g}")
        if k * ci(a) >= 2:
            chk(tag, f"L8: D >= phi2 k c_i(a), a={frac}*s_w", Dd >= PHI2 * k * ci(a), f"D={Dd:.5g} rhs={PHI2*k*ci(a):.5g}")
    # ---- L9 (phase ii)
    u = s0 - sw
    vw = -at(sw)[0]; Pw = -at(sw)[1]
    phw = 1 - np.exp(-vw)
    chk(tag, "L9: P_w >= k phi(s_w)(u-1+e^-u)", Pw >= k * phw * (u - 1 + np.exp(-u)) * (1 - 1e-9), f"P_w={Pw:.5g} rhs={k*phw*(u-1+np.exp(-u)):.5g}")
    chk(tag, "L9: D-v_w >= k phi(s_w) c_ii(u)", Dd - vw >= k * phw * cii(u) * (1 - 1e-9), f"D-v_w={Dd-vw:.5g} rhs={k*phw*cii(u):.5g}")
    chk(tag, "L9: P_w < (gamma_E-mu) v_w", Pw < (gE - mu) * vw, f"P_w={Pw:.5g} rhs={(gE-mu)*vw:.5g} ratio={Pw/((gE-mu)*vw):.4f}")
    kap = k * (u - 1 + np.exp(-u)) / (k ** (1 / 3) + 1)
    if kap >= 2:
        chk(tag, "L9: kappa_ii(u)>=2 => D >= phi2 k c_ii(u)", Dd >= PHI2 * k * cii(u), f"D={Dd:.5g} rhs={PHI2*k*cii(u):.5g}")
    # ---- L10 (rising phase)
    chk(tag, "L10a: T_- - s0 >= p_k/k, D >= p_k^2/2k", Tm - s0 >= pk / k and Dd >= pk * pk / (2 * k), f"{Tm-s0:.5g}>={pk/k:.5g}; {Dd:.5g}>={pk*pk/(2*k):.5g}")
    Uris = lambda t: D_(t_d + t)[0]
    for fr in (1.0, 0.5, 0.1):
        L = fr * pk * pk / (2 * k)
        if L < Dd:
            tL = brentq(lambda t: Uris(t) + L, s0, Tm, xtol=1e-14)
            chk(tag, f"L10a: T_- - t_L <= 2L/p_k (L={fr}*p_k^2/2k)", Tm - tL <= 2 * L / pk * (1 + 1e-9), f"{Tm-tL:.5g} <= {2*L/pk:.5g}")
    ells = [l for l in (2.0, Dd / 2, Dd / 10) if l < Dd]
    for l in ells:
        phl = 1 - np.exp(-l)
        sl = brentq(lambda t: Uris(t) + l, s0, Tm, xtol=1e-14)
        rl = sl - s0; taub = Tm - sl
        psl = at(sl)[1]
        chk(tag, f"L10b: r_l >= sqrt(2(D-l)/k) (l={l:.4g})", rl >= np.sqrt(2 * (Dd - l) / k), f"{rl:.5g} >= {np.sqrt(2*(Dd-l)/k):.5g}")
        chk(tag, f"L10c: p(s_l) >= k phi_l Psi_mu(r_l) (l={l:.4g})", psl >= k * phl * Psi(mu, rl) * (1 - 1e-9), f"{psl:.6g} >= {k*phl*Psi(mu,rl):.6g}")
        chk(tag, f"L10c: p_k >= p(s_l) - mu l (l={l:.4g})", pk >= psl - mu * l, f"{pk:.6g} >= {psl-mu*l:.6g}")
        chk(tag, f"L10d: k c_b(tau_b) < l/phi_l (contrapositive, l={l:.4g})", k * cb(taub) < l / phl, f"tau_b={taub:.5g} k c_b={k*cb(taub):.5g} < {l/phl:.5g}")
    # ---- L12
    loc = lambda t: D_(t_u + t)  # time after up-crossing
    tg = np.linspace(0, tau1, 200001)
    Yg = loc(tg)
    P1 = Yg[1].max()
    wneg = np.nonzero(Yg[2] > 0)[0]
    splus = tg[wneg[0]] if len(wneg) else tau1
    sg = tg[tg <= splus]
    chk(tag, "L12: p_k(s^2/2-s^3/6) < 1 on [0,s_+]", np.all(pk * (sg ** 2 / 2 - sg ** 3 / 6) < 1), f"s_+={splus:.5g} max={np.max(pk*(sg**2/2-sg**3/6)):.4f}")
    if pk >= 10:
        eps1 = 1.1 * np.sqrt(2) * k * pk ** (-1.5)
        chk(tag, "L12: P_1 <= p_k(1+eps_1)", P1 <= pk * (1 + eps1), f"P1/p_k={P1/pk:.5f} 1+eps1={1+eps1:.5f}")
        chk(tag, "L12: s_+ < 1.1 sqrt(2/p_k)", splus < 1.1 * np.sqrt(2 / pk))
    chk(tag, "L12: P_1 < lambda", P1 < lam)
    # ---- L13
    th = 2 ** -0.5
    H = np.log(M0)
    idx = np.nonzero(Yg[1] <= th * pk)[0][0]
    sth = brentq(lambda t: loc(t)[1] - th * pk, tg[idx - 1], tg[idx], xtol=1e-15)
    Uth = loc(sth)[0]
    Eth = th * pk * (0.5 * pk * pk - 2 * mu * P1 * H) / (2 * k)
    chk(tag, "L13: e^{U_theta} >= E_theta (with H=log M0)", np.exp(Uth) >= Eth, f"e^U_th={np.exp(Uth):.5g} E_th={Eth:.5g}")
    Eth2 = th * pk * (0.5 * pk * pk - 2 * mu * P1 * h) / (2 * k)
    chk(tag, "L13: e^{U_theta} >= E_theta (sharper, H:=h)", np.exp(Uth) >= Eth2, f"E_th(h)={Eth2:.5g}")
    Eth3 = th * pk * (0.5 * pk * pk - 2 * mu * P1 * Uth) / (2 * k)
    chk(tag, "L13: e^{U_theta} >= E_theta (sharpest, H:=U_theta)", np.exp(Uth) >= Eth3, f"E_th(U_th)={Eth3:.5g}")
    chk(tag, "L13: tau_1 - s_theta >= theta p_k/(mu P_1 + M_+)", tau1 - sth >= th * pk / (mu * P1 + Mplus), f"{tau1-sth:.5g} >= {th*pk/(mu*P1+Mplus):.5g}")
    Qv = 0.25 * pk * pk * (0.5 * pk * pk - 2 * mu * P1 * H) - pk * k / np.sqrt(2)
    chk(tag, "L13: M_+(M_+ + mu P_1) >= Q", Mplus * (Mplus + mu * P1) >= Qv, f"lhs={Mplus*(Mplus+mu*P1):.5g} Q={Qv:.5g}")
    Qh = 0.25 * pk * pk * (0.5 * pk * pk - 2 * mu * P1 * h) - pk * k / np.sqrt(2)
    chk(tag, "L13: M_+(M_+ + mu P_1) >= Q (sharper, H:=h)", Mplus * (Mplus + mu * P1) >= Qh, f"Q(h)={Qh:.5g}")
    if pk >= 10:
        Ht = 2 * np.log(k + 1) + 3 / np.e
        eps1 = 1.1 * np.sqrt(2) * k * pk ** (-1.5)
        dQ = 4 * (1 + eps1) * Ht / pk + 4 * np.sqrt(2) * k / pk ** 3
        if dQ < 1:
            chk(tag, "L13: M_+/k >= p_k^2 sqrt(1-dQ)/(sqrt8 k) - 1/2", Mplus / k >= pk * pk * np.sqrt(1 - dQ) / (np.sqrt(8) * k) - 0.5,
                f"M_+/k={Mplus/k:.5g} rhs={pk*pk*np.sqrt(1-dQ)/(np.sqrt(8)*k)-0.5:.5g}  M_+/p_k^2={Mplus/pk**2:.4f}")
    # ---- L17 (linear estimate) for several (nu, s)
    worst = np.inf
    for nu in (0.5 * mu, 0.25 * mu, 0.9 * mu, 0.05 * mu):
        c0 = (mu - nu) * (1 - nu) * (A - nu)
        tauc = np.count_nonzero(xx[:-1] > c0) * (tt[1] - tt[0])
        al = 1 / (1 - nu); be = al / (A - nu)
        for s_ in (np.sqrt(2), 2.0, 3.0, 2 * np.sqrt(lam), 10.0, 50.0):
            if s_ < np.sqrt(2):
                continue
            rhs_ = -nu * (T - tauc) + (s_ - mu) * tauc + k * T / s_ ** 2 + np.log(1 + al * s_ + be * s_ * s_)
            lhs_ = (avg_eta - mu) * T
            worst = min(worst, rhs_ - lhs_)
            if not lhs_ <= rhs_ + 1e-6:
                chk(tag, f"L17 nu={nu:.4g} s={s_:.4g}", False, f"lhs={lhs_:.6g} rhs={rhs_:.6g} tau_c={tauc:.5g}")
    chk(tag, "L17: (<eta>-mu)T <= bound for 4 nu x 6 s", worst >= -1e-6, f"min slack={worst:.5g}")
    # compound-system cross check of L17 mechanism: integrate (a,b,c) one period from positive data; growth<=bound
    return dict(mu=mu, k=k, T=T, Tm=Tm, Tp=Tp, pk=pk, D=Dd, avg_eta=avg_eta)


def main():
    if len(sys.argv) == 3:
        cases = [(float(sys.argv[1]), float(sys.argv[2]))]
    else:
        cases = [(mu, rho) for mu in (1.0, 0.3, 0.02) for rho in (1.01, 1.5, 4.0, 12.0)]
    guess = None
    last_mu = None
    for mu, rho in cases:
        k = rho * mu * (1 + mu)
        if mu != last_mu:
            guess = None
            last_mu = mu
        try:
            try:
                cyc = find_cycle(mu, k)
            except Exception as e0:
                print(f"   (direct find_cycle raised {e0!r})")
                cyc = dict(U0=0.0, residual=1.0)
            if cyc["U0"] < 1e-6 or cyc["residual"] > 1e-8:
                # continuation in rho from 1.5 downwards (whlib warm-up lands on E near Hopf for small mu)
                print(f"   (direct find_cycle landed on E for mu={mu} rho={rho}; continuing down from rho=1.5)")
                g = None
                for r_ in [r for r in [1.5, 1.35, 1.2, 1.12, 1.07, 1.04, 1.02, 1.01, 1.005] if r > rho] + [rho]:
                    kk = r_ * mu * (1 + mu)
                    cyc = find_cycle(mu, kk, guess=g)
                    g = (cyc["U0"], cyc["y0"] * (1.0))
                    print(f"     rho={r_}: U0={cyc['U0']:.5g} T={cyc['T']:.5f} res={cyc['residual']:.1e} m_t={cyc['m_t']:.6f}")
        except Exception as e:  # report, do not hide
            print(f"== mu={mu} rho={rho}: find_cycle FAILED: {e!r}")
            FAILS.append((f"mu={mu} rho={rho}", "find_cycle", repr(e)))
            continue
        if cyc["residual"] > 1e-8 or cyc["T"] < 1e-3:
            print(f"== mu={mu} rho={rho}: find_cycle residual {cyc['residual']:.2e} -- not a converged cycle; skipped")
            FAILS.append((f"mu={mu} rho={rho}", "find_cycle residual", str(cyc["residual"])))
            continue
        amp = cyc["U0"]
        if amp < 1e-6:
            print(f"== mu={mu} rho={rho}: find_cycle returned the equilibrium (U0={amp:.2e}); skipped")
            FAILS.append((f"mu={mu} rho={rho}", "find_cycle landed on E", ""))
            continue
        analyse(mu, k, cyc)
    print()
    print("FAILS:", FAILS if FAILS else "none")


if __name__ == "__main__":
    main()
