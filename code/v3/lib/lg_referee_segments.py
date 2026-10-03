#!/usr/bin/env python3
# --- v3.0 bundle copy of research file `lanes/largegain/verify/v_segments.py`.
# --- Changes from the original: none.
"""Referee script D: test the DERIVATION STEPS of the lemmas on orbits that are not cycles.
D1: negative lobes from random data (U=0, p<0, w>0)  -> L8, L9 (first two), L10(a)-(d).
D2: spikes from random data (U=0, p=p_k>0, -k<w<0)   -> L12, L13 (energy step, crest-time step, mass step),
                                                        L16 (w lower bound, crest time, W_c bound, post-crest drop).
D3: Lemma 17 as a statement about the linear cooperative system a'=-mu a+b, b'=-b+c, c'=-Ac+x a
    for ARBITRARY T-periodic x(t)>=0: spectral radius of the monodromy <= the claimed bound.
Independent implementation; the lane's test_segments.py / test_stage.py were not read or used."""
import sys
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

rng = np.random.default_rng(20261003)
N = int(sys.argv[1]) if len(sys.argv) > 1 else 300
PHI2 = 0.7968
nfail = 0
ntest = 0
worst = {}


def rec(name, slack, scale=1.0, info=""):
    """slack >= 0 means the inequality holds; tolerance 1e-7 relative to scale."""
    global nfail, ntest
    ntest += 1
    if not (np.isfinite(slack) and np.isfinite(scale)):   # never let a NaN pass silently
        nfail += 1
        print(f"FAIL(non-finite) {name}: slack={slack} scale={scale} {info}")
        return
    rel = slack / max(abs(scale), 1e-300)
    if name not in worst or rel < worst[name][0]:
        worst[name] = (rel, info)
    if slack < -1e-7 * max(abs(scale), 1e-12):
        nfail += 1
        print(f"FAIL {name}: slack={slack:.3e} scale={scale:.3e} {info}")


def Psi(mu, u):
    if abs(mu - 1) < 1e-9:
        return 1 - (1 + u) * np.exp(-u)
    return (1 - np.exp(-mu * u)) / mu - (np.exp(-u) - np.exp(-mu * u)) / (mu - 1)


ci = lambda a: (np.exp(a) - np.exp(a / 2)) * (1 - (1 + a / 2) * np.exp(-a / 2))
cii = lambda a: a * a / 2 - 1 + (1 + a) * np.exp(-a)
cb = lambda a: (1 - (1 + a / 2) * np.exp(-a / 2)) * (1 - np.exp(-a / 2))


def field(mu, k):
    def f(t, s):
        U, p, y, m = s
        x = k * np.exp(min(U, 700.0))
        return [p, k - mu * p - y, x - y, x - k]
    return f


def root(fun, a, b):
    return brentq(fun, a, b, xtol=1e-14, rtol=1e-14)


# ---------------- D1: negative lobes ----------------
def d1(n):
    done = 0
    tries = 0
    while done < n and tries < 20 * n:
        tries += 1
        mu = float(np.exp(rng.uniform(np.log(0.01), 0.0)))
        k = float(np.exp(rng.uniform(np.log(0.3), np.log(3000.0))))
        Pd = float(np.exp(rng.uniform(np.log(1e-2), np.log(50.0)))) * k ** rng.uniform(0, 0.7)
        wd = float(np.exp(rng.uniform(np.log(1e-3), np.log(30.0)))) * k
        f = field(mu, k)

        def up(t, s):
            return s[0]
        up.terminal = True; up.direction = 1
        sol = solve_ivp(f, [0, 80.0], [0.0, -Pd, wd + k, 0.0], method="DOP853", rtol=1e-11, atol=1e-13, events=up, dense_output=True)
        if sol.status != 1:
            continue
        Tm = sol.t[-1]
        S = sol.sol
        U = lambda t: float(S(t)[0]); P = lambda t: float(S(t)[1]); Wf = lambda t: float(S(t)[2]) - k
        tg = np.linspace(0, Tm, 4001)
        Y = S(tg)
        ip = np.nonzero(Y[1] > 0)[0]
        if len(ip) == 0 or ip[0] == 0:
            continue
        s0 = root(P, tg[ip[0] - 1], tg[ip[0]])
        iw = np.nonzero(Y[2] - k < 0)[0]
        sw = root(Wf, tg[iw[0] - 1], tg[iw[0]])
        if not (0 < sw < s0 < Tm):
            rec("D1 structure s_w<s_0<T_-", -1.0, 1.0, f"mu={mu:.3g} k={k:.4g}")
            continue
        # 0<|w|<k on the rising phase is y>0, exact for y(0)>0 (y' = x - y, x>0); numerically y decays to
        # the integrator noise floor, so only a violation beyond 1e-9 k is counted
        rec("D1 y>0 on rising phase (noise floor 1e-9k)", np.min(Y[2][tg >= s0]) + 1e-9 * k, k, f"mu={mu:.3g} k={k:.4g} min y={np.min(Y[2][tg >= s0]):.2e}")
        D = -U(s0); pk = P(Tm)
        tagi = f"mu={mu:.3g} k={k:.4g} Pd={Pd:.3g} wd={wd:.3g} D={D:.4g} T_-={Tm:.4g}"
        # L8
        for fr in (1.0, 0.5):
            a = fr * sw
            v = -U(a / 2)
            rec("L8 v>=k c_i(a)(1-e^-v)", v - k * ci(a) * (1 - np.exp(-v)), v, tagi)
            if k * ci(a) >= 2:
                rec("L8 D>=phi2 k c_i(a)", D - PHI2 * k * ci(a), D, tagi)
        # L9 (first two)
        u = s0 - sw
        vw = -U(sw); Pw = -P(sw); phw = 1 - np.exp(-vw)
        rec("L9 P_w>=k phi(s_w)(u-1+e^-u)", Pw - k * phw * (u - 1 + np.exp(-u)), Pw, tagi)
        rec("L9 D-v_w>=k phi(s_w)c_ii(u)", (D - vw) - k * phw * cii(u), D - vw, tagi)
        # L9 conclusion under its own hypothesis, with the cycle-only input P_w<(gE-mu)v_w checked as a hypothesis
        r = np.roots([1, 1 + mu, mu, k]); gE = -min(r.real[np.abs(r.imag) < 1e-9])
        if Pw < (gE - mu) * vw and k * (u - 1 + np.exp(-u)) / (k ** (1 / 3) + 1) >= 2 and k > mu * (1 + mu):
            rec("L9 D>=phi2 k c_ii(u) [given Lambda_E(s_w)>0]", D - PHI2 * k * cii(u), D, tagi)
        # L10
        rec("L10a T_- - s0 >= p_k/k", (Tm - s0) - pk / k, Tm - s0, tagi)
        rec("L10a D>=p_k^2/2k", D - pk * pk / (2 * k), D, tagi)
        for fr in (1.0, 0.3):
            L = fr * pk * pk / (2 * k)
            tL = root(lambda t: U(t) + L, s0, Tm)
            rec("L10a T_- - t_L <= 2L/p_k", 2 * L / pk - (Tm - tL), 2 * L / pk, tagi)
        for l in [x for x in (2.0, D / 2, D / 20) if x < D]:
            phl = 1 - np.exp(-l)
            sl = root(lambda t: U(t) + l, s0, Tm)
            rl = sl - s0; taub = Tm - sl; psl = P(sl)
            rec("L10b r_l>=sqrt(2(D-l)/k)", rl - np.sqrt(2 * (D - l) / k), rl, tagi)
            rec("L10c p(s_l)>=k phi_l Psi_mu(r_l)", psl - k * phl * Psi(mu, rl), psl, tagi)
            rec("L10c p_k>=p(s_l)-mu l", pk - (psl - mu * l), abs(pk) + mu * l, tagi)
            rec("L10d k c_b(tau_b) < l/phi_l", l / phl - k * cb(taub), l / phl, tagi)
        done += 1
    return done


# ---------------- D2: spikes ----------------
def d2(n):
    done = 0
    tries = 0
    th = 2 ** -0.5
    while done < n and tries < 20 * n:
        tries += 1
        mu = float(np.exp(rng.uniform(np.log(0.01), 0.0)))
        k = float(np.exp(rng.uniform(np.log(0.3), np.log(3000.0))))
        lam = k / mu
        pk = lam * float(np.exp(rng.uniform(np.log(1e-3), 0.0))) * 0.999
        wk = -k * rng.uniform(0.001, 0.999)
        f = field(mu, k)

        def crest(t, s):
            return s[1]
        crest.terminal = True; crest.direction = -1
        sol = solve_ivp(f, [0, 200.0], [0.0, pk, wk + k, 0.0], method="DOP853", rtol=1e-11, atol=1e-13, events=crest, dense_output=True)
        if sol.status != 1:
            continue
        tau1 = sol.t[-1]
        S = sol.sol
        tg = np.unique(np.concatenate([np.linspace(0, tau1, 6001), sol.t]))
        Y = S(tg)
        Ug, Pg, Wg, Mg = Y[0], Y[1], Y[2] - k, Y[3]
        h = Ug[-1]; P1 = Pg.max(); Mc = Mg[-1]; Wc = Wg[-1]
        tagi = f"mu={mu:.3g} k={k:.4g} p_k={pk:.4g} w_k={wk:.4g} tau1={tau1:.4g} h={h:.4g}"
        if P1 >= lam:      # outside the positive domain: the lemma assumes P_1<lambda only through mu*P_1<k; skip
            continue
        # (wm): e^-s (m-|w_k|) <= w <= m
        rec("wm: w<=m", np.min(Mg - Wg), max(Mc, 1.0), tagi)
        rec("wm: w>=e^-s(m-|w_k|)", np.min(Wg - np.exp(-tg) * (Mg - abs(wk))), max(Mc, 1.0), tagi)
        # L12
        ipos = np.nonzero(Wg > 0)[0]
        splus = tg[ipos[0]] if len(ipos) else tau1
        sg = tg[tg <= splus]
        rec("L12 p_k(s^2/2-s^3/6)<1 on [0,s_+]", 1 - np.max(pk * (sg ** 2 / 2 - sg ** 3 / 6)), 1.0, tagi)
        if pk >= 10:
            rec("L12 s_+ < 1.1 sqrt(2/p_k)", 1.1 * np.sqrt(2 / pk) - splus, 1.0, tagi)
            rec("L12 P_1<=p_k(1+eps_1)", pk * (1 + 1.1 * np.sqrt(2) * k * pk ** -1.5) - P1, P1, tagi)
        # L13
        i = np.nonzero(Pg <= th * pk)[0][0]
        sth = root(lambda t: float(S(t)[1]) - th * pk, tg[i - 1], tg[i])
        Uth = float(S(sth)[0])
        # energy step: p_k^2/4 <= mu P_1 U_theta + k e^{U_theta}/(theta p_k)
        rec("L13 energy: p_k^2/4 <= mu P_1 U_th + k e^{U_th}/(theta p_k)", mu * P1 * Uth + k * np.exp(Uth) / (th * pk) - pk * pk / 4, pk * pk / 4, tagi)
        rec("L13 m(s)<=k(e^U-1-U)/(theta p_k) on [0,s_th]", np.min((k * (np.exp(Ug) - 1 - Ug) / (th * pk) - Mg)[tg <= sth]), max(Mc, 1.0), tagi)
        rec("L13 tau1-s_th >= theta p_k/(mu P_1+M_c)", (tau1 - sth) - th * pk / (mu * P1 + Mc), tau1 - sth, tagi)
        Eth = th * pk * (0.5 * pk * pk - 2 * mu * P1 * h) / (2 * k)
        rec("L13 e^{U_th} >= E_theta(H:=h)", np.exp(Uth) - Eth, np.exp(Uth), tagi)
        Q = 0.25 * pk * pk * (0.5 * pk * pk - 2 * mu * P1 * h) - pk * k / np.sqrt(2)
        rec("L13 M_c(M_c+mu P_1) >= Q(H:=h)", Mc * (Mc + mu * P1) - Q, Mc * (Mc + mu * P1), tagi)
        if Q > 0:
            rec("L13 M_c >= sqrt(Q) - mu P_1/2", Mc - (np.sqrt(Q) - mu * P1 / 2), Mc, tagi)
        # L16 pieces
        Xth = k * (np.exp(Uth) - 1)
        if Xth > 0:
            sel = (tg >= sth)
            v = tg[sel] - sth
            rec("L16 w(s_th+v) >= -k + X_th(1-e^-v) up to crest", np.min(Wg[sel] - (-k + Xth * (1 - np.exp(-v)))), max(abs(Wc), k), tagi)
            vc = 1.5 / Xth * (k + np.sqrt(k * k + 4 / 3 * Xth * th * pk))
            if vc <= 1:
                rec("L16 crest <= s_th + v_c", vc - (tau1 - sth), vc, tagi)
            rec("L16 s_th <= U_th/(theta p_k)", Uth / (th * pk) - sth, sth, tagi)
        rec("L16 p_k <= W_c tau1 e^{tau1}", Wc * tau1 * np.exp(min(tau1, 700)) - pk, pk, tagi)
        # post-crest
        if Wc > 0:
            vplus = np.log(1 + Wc / k)

            def trough(t, s):
                return s[1]
            trough.terminal = True; trough.direction = 1
            sol2 = solve_ivp(f, [0, vplus], sol.y[:, -1] * np.array([1, 0, 1, 1]), method="DOP853", rtol=1e-11, atol=1e-13, events=trough, dense_output=True)
            vv = np.linspace(0, sol2.t[-1], 3001)
            Y2 = sol2.sol(vv)
            w2 = Y2[2] - k
            rec("L7e w>=e^-v(W_c+k)-k after crest", np.min(w2 - (np.exp(-vv) * (Wc + k) - k)), Wc + k, tagi)
            lowp = (Wc + k) * vv * np.exp(-vv) + k * np.expm1(-vv)
            rec("L16 |p|(v) >= (W_c+k)v e^-v - k(1-e^-v)", np.min(-Y2[1] - lowp), max(np.max(-Y2[1]), 1e-9), tagi)
            lowb = (Wc + k) * (-np.expm1(-vv) - vv * np.exp(-vv)) - k * (vv + np.expm1(-vv))
            rec("L16 h-U(v) >= (W_c+k)(1-(1+v)e^-v)-k(v-1+e^-v)", np.min((h - Y2[0]) - lowb), max(np.max(h - Y2[0]), 1e-9), tagi)
            lowc = 0.5 * Wc * vv ** 2 * (1 - 2 / 3 * (1 + k / Wc) * vv)
            rec("L16 h-U(v) >= W_c v^2/2 [1-(2/3)(1+k/W_c)v]", np.min((h - Y2[0]) - lowc), max(np.max(h - Y2[0]), 1e-9), tagi)
        done += 1
    return done


# ---------------- D3: Lemma 17 for arbitrary periodic x ----------------
def d3(n):
    done = 0
    for _ in range(n):
        mu = float(np.exp(rng.uniform(np.log(0.02), 0.0)))
        A = 1 + mu
        nu = mu * rng.uniform(0.05, 0.95)
        c0 = (mu - nu) * (1 - nu) * (A - nu)
        T = float(np.exp(rng.uniform(np.log(1.0), np.log(200.0))))
        tauc = T * float(np.exp(rng.uniform(np.log(1e-3), np.log(0.9))))
        amp = float(np.exp(rng.uniform(np.log(1.0), np.log(1e5))))
        shape = rng.integers(0, 3)
        base = c0 * rng.uniform(0.0, 1.0)

        def x(t):
            t = t % T
            if t >= tauc:
                return base
            z = t / tauc
            if shape == 0:
                return base + amp
            if shape == 1:
                return base + amp * np.sin(np.pi * z) ** 2
            return base + amp * np.exp(-((z - 0.5) / 0.08) ** 2)

        def f(t, s):
            M = s.reshape(3, 3)
            xt = x(t)
            C = np.array([[-mu, 1, 0], [0, -1, 1], [xt, 0, -A]])
            return (C @ M).ravel()
        # integrate with renormalisation at the two breakpoints to avoid overflow
        logscale = 0.0
        M = np.eye(3)
        for (t0, t1) in ((0, tauc), (tauc, T)):
            nseg = max(1, int((t1 - t0) * (np.sqrt(np.sqrt(amp)) + 1) / 5) + 1)
            edges = np.linspace(t0, t1, nseg + 1)
            for j in range(nseg):
                sol = solve_ivp(f, [edges[j], edges[j + 1]], M.ravel(), method="DOP853", rtol=1e-10, atol=1e-30, max_step=max(tauc / 20, (edges[j + 1] - edges[j]) / 4))
                M = sol.y[:, -1].reshape(3, 3)
                sc = np.abs(M).max()
                M /= sc
                logscale += np.log(sc)
        ev = np.linalg.eigvals(M)
        logr = np.log(np.max(np.abs(ev))) + logscale
        from scipy.integrate import quad
        intx = quad(x, 0, tauc, limit=200)[0]
        al = 1 / (1 - nu); be = al / (A - nu)
        for s_ in (np.sqrt(2), 2.0, 5.0, max(np.sqrt(2), amp ** (1 / 3)), max(np.sqrt(2), np.sqrt(amp))):
            bound = -nu * (T - tauc) + (s_ - mu) * tauc + intx / s_ ** 2 + np.log(1 + al * s_ + be * s_ ** 2)
            rec("L17 log(spectral radius) <= bound (arbitrary periodic x)", bound - logr, max(abs(bound), abs(logr), 1.0),
                f"mu={mu:.3g} nu={nu:.3g} T={T:.3g} tauc={tauc:.3g} amp={amp:.3g} shape={shape} s={s_:.3g} logr={logr:.5g} bound={bound:.5g}")
        done += 1
    return done


if __name__ == "__main__":
    n1 = d1(N)
    n2 = d2(N)
    n3 = d3(max(40, N // 3))
    print(f"segments: {n1} negative lobes, {n2} spikes, {n3} periodic-x linear systems; {ntest} inequality tests, {nfail} failures")
    print("tightest relative slack per inequality:")
    for name, (rel, info) in sorted(worst.items()):
        print(f"  {rel:+.3e}  {name}   [{info}]")
    print("ALL SEGMENT TESTS PASSED" if nfail == 0 else f"FAILURES: {nfail}")
