#!/usr/bin/env python3
"""Referee script F: Stage-1 / Stage-2 implications (Lemmas 11 and 14) on NON-periodic negative lobes.
The proofs of L11/L14 use, besides the lobe structure, only: k>=1000, mu<=1, T_- > 5.89 (resp. >= 22),
D > 2, |w|<k on the rising phase, and Lambda_E(s_w) > 0 (i.e. P_w < (gamma_E-mu) v_w).  So on ANY orbit
segment from (U=0,p<0,w>0) to the next up-crossing satisfying those hypotheses the conclusions
  tau_b < 0.36,  [s_w>=2.06 or u>=1.99 or r_l>=1.48],  r_l >= 1.48,  p_k >= 0.374 k,
  (T_->=22:)  r_l >= 0.3 T_-,  p_k >= 0.95 phi_l k (1-e^{-mu z})/mu - 2,  z = 0.3 T_- - 3
must hold.  Independent implementation (closed-form deep trough from v_largegain)."""
import sys
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
import v_largegain as L

rng = np.random.default_rng(7)
N = int(sys.argv[1]) if len(sys.argv) > 1 else 400
UC = 60.0
phl = 1 - np.exp(-2.0)
stats = dict(lobes=0, hyp1=0, hyp2=0, noLam=0, fail=0)
mins = dict(pk_over_k=np.inf, rl=np.inf, rl_over_T=np.inf, taub=0.0, pkT_slack=np.inf)


def lobe(mu, k, Pd, wd):
    """negative lobe from (0,-Pd,wd): numerical while U>-UC, closed form below."""
    def f(t, s):
        U, p, y = s
        return [p, k - mu * p - y, k * np.exp(U) - y]

    def up(t, s):
        return s[0]
    up.terminal = True; up.direction = 1

    def cut(t, s):
        return s[0] + UC
    cut.terminal = True; cut.direction = -1
    pieces = []   # list of (t_offset, kind, object)
    sol = solve_ivp(f, [0, 400.0], [0.0, -Pd, wd + k], method="DOP853", rtol=1e-11, atol=1e-12, events=[up, cut], dense_output=True)
    if sol.status != 1:
        return None
    if len(sol.t_events[0]):          # shallow lobe: never reached -UC
        return dict(deep=False, sols=[(0.0, sol)], Tm=sol.t[-1])
    t1 = sol.t[-1]
    p0, w0 = sol.y[1, -1], sol.y[2, -1] - k
    U_, p_, w_ = L.trough_fns(mu, k, p0, w0, -UC)
    hi = 1.0
    while p_(hi) <= 0:
        hi *= 2
        if hi > 1e9:
            return None
    s0p = brentq(p_, 0.0, hi, xtol=1e-14)
    hi = max(2 * s0p, 1.0)
    while U_(hi) + UC <= 0:
        hi *= 2
    t_r = brentq(lambda t: U_(t) + UC, s0p, hi, xtol=1e-12, rtol=1e-15)
    sol2 = solve_ivp(f, [0, 400.0], [-UC, p_(t_r), w_(t_r) + k], method="DOP853", rtol=1e-11, atol=1e-12, events=[up], dense_output=True)
    if sol2.status != 1:
        return None
    return dict(deep=True, sol=sol, t1=t1, U_=U_, p_=p_, w_=w_, s0p=s0p, t_r=t_r, sol2=sol2, Tm=t1 + t_r + sol2.t[-1], p0=p0, w0=w0)


def main():
    done = 0
    tries = 0
    while done < N and tries < 40 * N:
        tries += 1
        mu = float(np.exp(rng.uniform(np.log(0.01), 0.0)))
        k = float(np.exp(rng.uniform(np.log(1000.0), np.log(6000.0))))
        wd = k * float(np.exp(rng.uniform(np.log(0.05), np.log(300.0))))
        Pd = float(np.exp(rng.uniform(np.log(1e-2), np.log(1e4))))
        lb = lobe(mu, k, Pd, wd)
        if lb is None or not lb["deep"]:
            continue
        done += 1
        stats["lobes"] += 1
        Tm = lb["Tm"]
        sol, sol2 = lb["sol"], lb["sol2"]
        t1, t_r, s0p = lb["t1"], lb["t_r"], lb["s0p"]
        # s_w, s_0, D
        if lb["w0"] > 0:
            tw = np.log((lb["w0"] + k) / k)
            if tw >= s0p:
                print("STRUCTURE VIOLATION: s_w >= s_0", mu, k); stats["fail"] += 1; continue
            sw = t1 + tw
            vw = -lb["U_"](tw); Pw = -lb["p_"](tw)
        else:
            g = lambda t: float(sol.sol(t)[2]) - k
            sw = brentq(g, 0, t1, xtol=1e-14)
            vw = -float(sol.sol(sw)[0]); Pw = -float(sol.sol(sw)[1])
        s0 = t1 + s0p
        D = -lb["U_"](s0p)
        u = s0 - sw
        # s_l on the final rising piece
        tl = brentq(lambda t: float(sol2.sol(t)[0]) + 2.0, 0, sol2.t[-1], xtol=1e-15)
        taub = sol2.t[-1] - tl
        rl = (t_r - s0p) + tl
        pk = sol2.y[1, -1]
        r = np.roots([1, 1 + mu, mu, k]); gE = -min(r.real[np.abs(r.imag) < 1e-9])
        lamE_ok = Pw < (gE - mu) * vw
        # tau_b < 0.36 needs only k>=1000 and D>2
        mins["taub"] = max(mins["taub"], taub)
        if not taub < 0.36:
            print(f"FAIL tau_b<0.36: mu={mu:.3g} k={k:.5g} tau_b={taub:.4g}"); stats["fail"] += 1
        if not lamE_ok:
            stats["noLam"] += 1
        if Tm > 5.89 and lamE_ok and D > 2:
            stats["hyp1"] += 1
            pig = (sw >= 2.06) or (u >= 1.99) or (rl >= 1.48)
            mins["pk_over_k"] = min(mins["pk_over_k"], pk / k)
            mins["rl"] = min(mins["rl"], rl)
            if not (pig and rl >= 1.48 and pk >= 0.374 * k):
                print(f"FAIL Stage1: mu={mu:.3g} k={k:.5g} T_-={Tm:.4g} s_w={sw:.4g} u={u:.4g} r_l={rl:.4g} tau_b={taub:.4g} p_k/k={pk/k:.4g}")
                stats["fail"] += 1
        if Tm >= 22 and lamE_ok and D > 2:
            stats["hyp2"] += 1
            z = 0.3 * Tm - 3
            rhs = 0.95 * phl * k * (1 - np.exp(-mu * z)) / mu - 2.0
            mins["rl_over_T"] = min(mins["rl_over_T"], rl / Tm)
            mins["pkT_slack"] = min(mins["pkT_slack"], (pk - rhs) / pk)
            if not (rl >= 0.3 * Tm and pk >= rhs):
                print(f"FAIL Stage2: mu={mu:.3g} k={k:.5g} T_-={Tm:.4g} s_w={sw:.4g} u={u:.4g} r_l={rl:.4g} p_k={pk:.5g} rhs={rhs:.5g}")
                stats["fail"] += 1
    print(stats)
    print({a: float(f"{b:.5g}") for a, b in mins.items()})
    print("ALL STAGE TESTS PASSED" if stats["fail"] == 0 else "FAILURES")


if __name__ == '__main__':
    main()
