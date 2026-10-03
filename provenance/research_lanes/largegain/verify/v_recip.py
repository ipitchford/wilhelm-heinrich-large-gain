#!/usr/bin/env python3
"""Referee script E.
(1) Corollary 3 (reciprocal symmetry): cycles at (mu,k) and (1/mu, k/mu^3) must have periods
    T' = mu*T and the same weak multiplier, log m_t = (<eta>-mu)T = (<eta>'-mu')T'.
    Checked with the independent hybrid integrator of v_largegain.py at mu>1, k>=1000 mu^3,
    and at moderate gain with whlib (monodromy multipliers).
(2) Sanity for Corollary 2: the orbit leaving E along W^u(E) (innermost side) and an orbit started
    far outside land on the same cycle at k=1000 (floating-point evidence only).
"""
import sys
import numpy as np
from scipy.integrate import solve_ivp
sys.path.insert(0, "/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03/work")
from whlib import find_cycle
import v_largegain as L

ok = True


def cyc(mu, k):
    lam = k / mu
    st = (0.9 * lam, -0.9 * k, 0.0, 0.0)
    for it in range(6):
        st_new, rec = L.loop(mu, k, st)
        d = abs(st_new[0] - st[0]) / lam + abs(st_new[1] - st[1]) / k
        st = st_new
        if d < 1e-12 and it >= 2:
            break
    dense, nodes = rec["dense"], rec["nodes"]
    Uf = lambda t: float(dense(t)[0])
    t_up = L.bracket_root(Uf, nodes, +1)[0]
    t_dn = L.bracket_root(Uf, nodes, -1)[0]
    T = (rec["t_end"] - t_dn) + rec["t_r"] + t_up + (t_dn - t_up)
    Ieta = float(rec["s_end"][5]) + rec["Itr"]
    pk = float(dense(t_up)[1])
    return T, Ieta / T, pk


print("(1a) large gain, mu>1 vs reciprocal (hybrid integrator)")
for mu, kp in ((2.0, 1000.0), (2.0, 3000.0), (5.0, 1000.0), (1.0 / 0.3, 1000.0)):
    k = kp * mu ** 3
    T, e, pk = cyc(mu, k)
    Tp, ep, pkp = cyc(1 / mu, kp)
    lm = (e - mu) * T
    lmp = (ep - 1 / mu) * Tp
    good = abs(Tp - mu * T) < 1e-6 * Tp and abs(lm - lmp) < 1e-6 * abs(lm) and e <= mu - 0.17 and abs(pkp - pk / mu) < 1e-6 * pkp
    ok &= good
    print(f"  mu={mu:.4g} k={k:.6g}: T={T:.8g}  mu*T={mu*T:.8g}  T'={Tp:.8g} | log m_t={lm:.8g}  log m_t'={lmp:.8g} | <eta>={e:.6g} <= mu-0.17={mu-0.17:.6g} | p_k/mu={pk/mu:.8g} p_k'={pkp:.8g}  {'ok' if good else 'FAIL'}")

print("(1b) moderate gain, monodromy multipliers (whlib.find_cycle)")
for mu, rho in ((2.0, 1.5), (3.0, 4.0), (5.0, 2.0)):
    k = rho * mu * (1 + mu)
    c = find_cycle(mu, k)
    cp = find_cycle(1 / mu, k / mu ** 3)
    good = abs(cp["T"] - mu * c["T"]) < 1e-7 * cp["T"] and abs(c["m_t"] - cp["m_t"]) < 1e-7
    ok &= good
    print(f"  mu={mu} rho={rho}: T={c['T']:.9f} mu*T={mu*c['T']:.9f} T'={cp['T']:.9f} | m_t={c['m_t']:.10f} m_t'={cp['m_t']:.10f}  {'ok' if good else 'FAIL'}")

print("(2) W^u(E) orbit and widely spread section data all land on the same cycle (floating-point evidence only)")
from scipy.optimize import brentq
for mu, k in ((1.0, 1000.0), (0.3, 1000.0), (0.02, 1000.0), (1.0, 3000.0)):
    lam = k / mu
    # reference fixed point of the section map {U=-UCUT, p>0}
    st = (0.9 * lam, -0.9 * k, 0.0, 0.0)
    for it in range(5):
        st, _ = L.loop(mu, k, st)
    ref = st
    # (2a) orbit on W^u(E)
    ev_, vec_ = np.linalg.eig(np.array([[0, 1, 0], [0, -mu, -1], [k, 0, -1.0]]))
    i_ = int(np.argmax(ev_.imag))
    v = np.real(vec_[:, i_]); v /= np.linalg.norm(v)

    def f(t, s):
        U, p, y = s
        return [p, k - mu * p - y, k * np.exp(U) - y]

    def cut(t, s):
        return s[0] + L.UCUT
    cut.terminal = True; cut.direction = -1
    sol = solve_ivp(f, [0, 200.0], np.array([0.0, 0.0, k]) + 1e-6 * v, method="DOP853", rtol=1e-11, atol=1e-12, events=cut)
    assert sol.status == 1, "W^u(E) orbit did not reach the trough section"
    p0, w0 = sol.y[1, -1], sol.y[2, -1] - k
    U_, p_, w_ = L.trough_fns(mu, k, p0, w0, -L.UCUT)
    hi = 1.0
    while p_(hi) <= 0:
        hi *= 2
    s0p = brentq(p_, 0.0, hi)
    hi = max(2 * s0p, 1.0)
    while U_(hi) + L.UCUT <= 0:
        hi *= 2
    t_r = brentq(lambda t: U_(t) + L.UCUT, s0p, hi, xtol=1e-12, rtol=1e-15)
    stw = (p_(t_r), w_(t_r), 0.0, 0.0)
    d_first = abs(stw[0] - ref[0]) / lam + abs(stw[1] - ref[1]) / k
    for it in range(3):
        stw, _ = L.loop(mu, k, stw)
    d_w = abs(stw[0] - ref[0]) / lam + abs(stw[1] - ref[1]) / k
    # (2b) spread of section data
    dmax = 0.0
    for pf_ in (1e-3, 0.05, 0.5, 0.999999):
        for wf_ in (-0.999, -0.5, -0.001, 5.0):
            s_ = (pf_ * lam, wf_ * k, 0.0, 0.0)
            for it in range(4):
                s_, _ = L.loop(mu, k, s_)
            dmax = max(dmax, abs(s_[0] - ref[0]) / lam + abs(s_[1] - ref[1]) / k)
    good = d_w < 1e-9 and dmax < 1e-9
    ok &= good
    print(f"  mu={mu} k={k}: W^u(E) orbit reaches the trough section at t={sol.t[-1]:.3f}; distance to cycle section point: first hit {d_first:.2e}, after 3 loops {d_w:.2e}; 16 spread section data after 4 loops: max distance {dmax:.2e}  {'ok' if good else 'FAIL'}")
print("ALL OK" if ok else "SOME FAIL")
