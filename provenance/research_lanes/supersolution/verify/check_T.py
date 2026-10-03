#!/usr/bin/env python3
"""Referee check 6 (Proposition T and the imported Theorem 3.1).
 (a) Theorem 3.1 two-phase form on one cycle: Delta(s,s') = zeta(s)(U(s')-U(s)) - eta(s)(p(s')-p(s)) - (y(s')-y(s)) < 0, s != s'.
 (b) polar orbit equation dlog r/dtheta = Xi(m,theta) by finite differences of the actual trajectory.
 (c) the monotonicity mechanism with REAL nested curves: successive turns of a W^u(E) spiral and the cycle
     (conditional on brief 3(b); two distinct cycles at one parameter are not available numerically):
     log(R_out/R_in) must increase on (0,pi) and decrease on (pi,2pi); m_out > m_in at equal angle."""
import sys, time
import numpy as np
from scipy.integrate import solve_ivp
ROOT = '/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03'
sys.path.insert(0, ROOT + '/work')
from whlib import equilibrium_data
from vprofile import profile, load_cycles
C = load_cycles()

print("(a),(b)")
print(f"{'mu':>5} {'rho':>5} {'max offdiag Delta/scale':>23} {'polar eq err':>12} {'min sigma':>10}")
for key in [(0.02, 1.01), (0.02, 4.0), (0.3, 1.5), (0.3, 12.0), (1.0, 1.01), (1.0, 4.0), (1.0, 12.0)]:
    mu, rho = key; cyc = C[key]; k = cyc['k']; lam = k/mu; T = cyc['T']
    R = profile(mu, k, cyc)
    n = 700
    # phases concentrated where the solver took steps (crest region) plus uniform
    ts = np.unique(np.concatenate([np.linspace(0, T, n//2, endpoint=False), R['sol'].t[np.linspace(0, len(R['sol'].t)-2, n//2).astype(int)]]))
    Y = R['sol'].sol(ts); U, P, y, eta, zeta = Y[:5]; p = lam - P
    # p(s') - p(s) = -(P(s') - P(s)): difference P directly (p = lam - P loses the relative accuracy of P in deep troughs)
    dp = -(P[None, :] - P[:, None])
    Dl = zeta[:, None]*(U[None, :] - U[:, None]) - eta[:, None]*dp - (y[None, :] - y[:, None])
    scale = np.abs(zeta[:, None]*(U[None, :] - U[:, None])) + np.abs(eta[:, None]*dp) + np.abs(y[None, :] - y[:, None]) + 1e-300
    off = ~np.eye(len(ts), dtype=bool)
    worst = np.max((Dl/scale)[off])
    # polar equation by central differences of the trajectory
    tm = np.linspace(0.013*T, 0.987*T, 4001); h = 1e-5*T
    def pol(t):
        Uv, Pv, yv = R['sol'].sol(t)[:3]; pv = lam - Pv
        return 0.5*np.log(Uv*Uv + pv*pv), np.arctan2(pv, Uv), Uv, pv, yv
    lr1, th1, *_ = pol(tm - h); lr2, th2, *_ = pol(tm + h); lr0, th0, U0, p0, y0 = pol(tm)
    dth = np.angle(np.exp(1j*(th2 - th1)))
    lhs = (lr2 - lr1)/dth
    r0 = np.exp(lr0); m = (y0 - k)/r0; s, c = p0/r0, U0/r0
    sig = s*s + mu*s*c + m*c
    rhs = s*(m + mu*s - c)/sig
    err = np.max(np.abs(lhs - rhs)/(1 + np.abs(rhs)))
    print(f"{mu:5.2f} {rho:5.2f} {worst:23.3e} {err:12.2e} {sig.min():10.3e}", flush=True)

print("\n(c) spiral turns vs cycle: fraction of neighbouring grid angles (step pi/360) between which log(R_out/R_in) moves in the direction sign(sin theta), and min(m_out - m_in)")
print(f"{'mu':>5} {'rho':>5} {'pair':>12} {'R ratio range':>22} {'sign ok frac':>12} {'min signed increment':>22} {'min (m_out-m_in)':>16}")
for key in [(1.0, 1.5), (0.3, 1.5), (0.3, 4.0), (1.0, 4.0), (0.02, 1.5)]:
    mu, rho = key; cyc = C[key]; k = cyc['k']; lam = k/mu; T = cyc['T']
    ed = equilibrium_data(mu, k)
    v = np.real(ed['vec_u']); z0 = np.array([0.0, 0.0, k]) + 1e-5*v/np.linalg.norm(v)
    def f(t, s):
        U, p, y, ph = s
        pd = k - mu*p - y
        return [p, pd, k*np.exp(U) - y, (p*p - U*pd)/(U*U + p*p)]
    tend = (np.log(1e5) + 14.0)/ed['lam_u'].real
    sol = solve_ivp(f, [0, tend], [z0[0], z0[1], z0[2], 0.0], method='DOP853', rtol=1e-12, atol=1e-15, dense_output=True)
    th_abs0 = np.arctan2(z0[1], z0[0])
    te = np.unique(np.concatenate([sol.t] + [sol.t[:-1] + j/8*np.diff(sol.t) for j in range(1, 8)]))
    Phi = sol.sol(te)[3] - th_abs0        # = -theta (unwrapped), strictly increasing because S1 > 0
    assert np.all(np.diff(Phi) > 0)
    thg = np.linspace(0, 2*np.pi, 721)[1:-1]
    def turn(nq):                          # log R(theta), m(theta) on the turn where -theta in [2pi nq - 2pi, 2pi nq]
        target = 2*np.pi*nq - thg
        j = np.searchsorted(Phi, target)   # bracket, then vectorised bisection on the dense output (Phi is increasing)
        lo, hi = te[j-1].copy(), te[j].copy()
        for _ in range(60):
            mid_t = 0.5*(lo + hi)
            ph = sol.sol(mid_t)[3] - th_abs0
            up = ph < target
            lo = np.where(up, mid_t, lo); hi = np.where(up, hi, mid_t)
        t = 0.5*(lo + hi)
        U, p, y, ph = sol.sol(t)
        assert np.max(np.abs(ph - th_abs0 - target)) < 1e-8, np.max(np.abs(ph - th_abs0 - target))
        r = np.hypot(U, p)
        return np.log(r), (y - k)/r
    nmax = int(np.floor((Phi[-1] - 1e-9)/(2*np.pi)))
    # cycle
    Rc = profile(mu, k, cyc)
    tsc = Rc['sol'].t
    tc = np.unique(np.concatenate([tsc] + [tsc[:-1] + j/16*np.diff(tsc) for j in range(1, 16)]))
    Yc = Rc['sol'].sol(tc); Uc, Pc, yc = Yc[:3]; pc = lam - Pc
    thc = np.mod(np.arctan2(pc, Uc), 2*np.pi); o = np.argsort(thc)
    # cycle starts at the crest (theta = 0 = 2pi) and theta decreases from 2pi to 0 over one period
    thdec = np.where(tc == 0, 2*np.pi, thc); thdec[-1] = 0.0
    assert np.all(np.diff(thdec) < 0)
    j = np.searchsorted(-thdec, -thg)
    lo, hi = tc[j-1].copy(), tc[j].copy()
    for _ in range(60):
        mid_t = 0.5*(lo + hi)
        Um, Pm, ym = Rc['sol'].sol(mid_t)[:3]
        thm = np.mod(np.arctan2(lam - Pm, Um), 2*np.pi)
        up = thm > thg
        lo = np.where(up, mid_t, lo); hi = np.where(up, hi, mid_t)
    tcg = 0.5*(lo + hi)
    Uc, Pc, yc = Rc['sol'].sol(tcg)[:3]; pc = lam - Pc
    assert np.max(np.abs(np.angle(np.exp(1j*(np.arctan2(pc, Uc) - thg))))) < 1e-8
    lrc = 0.5*np.log(Uc**2 + pc**2); mc = (yc - k)/np.hypot(Uc, pc)
    sgn = np.sign(np.sin(thg))
    def report(name, lr_out, m_out, lr_in, m_in):
        q = lr_out - lr_in
        dq = np.diff(q)                    # direct monotonicity test between neighbouring grid angles on the same side of pi
        mid = 0.5*(thg[1:] + thg[:-1]); same = (thg[1:] - np.pi)*(thg[:-1] - np.pi) > 0
        sg = np.sign(np.sin(mid))
        ok = (dq*sg > 0)[same]
        wrong = np.min((dq*sg)[same])
        print(f"{mu:5.2f} {rho:5.2f} {name:>12} {np.exp(q.min()):10.5f} {np.exp(q.max()):11.5f} {ok.mean():12.4f} {wrong:22.3e} {np.min(m_out - m_in):16.3e}", flush=True)
    turns = [nq for nq in range(2, nmax+1)]   # turn 1 may start before the first sample (Phi[0] = -theta_0 > 0)
    # choose turns with amplitude between 2% and 98% of the cycle radius
    sel = []
    for nq in turns:
        lr, mm = turn(nq)
        if nq+1 <= nmax and np.exp(np.max(lr - lrc)) < 0.999 and np.exp(np.min(lr - lrc)) > 0.02:
            sel.append(nq)
    for nq in sel[::max(1, len(sel)//4)]:
        lr_in, m_in = turn(nq); lr_out, m_out = turn(nq+1)
        report(f"turn {nq+1}/{nq}", lr_out, m_out, lr_in, m_in)
        report(f"cycle/{nq}", lrc, mc, lr_in, m_in)
