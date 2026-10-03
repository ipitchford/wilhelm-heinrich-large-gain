#!/usr/bin/env python3
"""Referee numerics for Theorem C (NOTES 4.1).  The theorem is asymptotic (mu0 not explicit), so numerics cannot
refute it outright; what can be tested is (i) the two quantitative conclusions with their claimed error orders at
FIXED r (the lane tested at fixed amplitude), (ii) Step 1 (cycle on the slow manifold: Q - w_series = O(eps^4)),
(iii) uniqueness inside an energy ball by two independent searches.

 (P0) cross-check of the referee's cycle finder against whlib.find_cycle
 (P1) fixed r in {0.3, 1, 3}, mu -> 0:  crest energy vs h*(r),  log m_t vs -eps^3 J0 Rsel'(h*),
      sup |Q - (w0+eps w1+eps^2 w2+eps^3 w3)| / eps^4 along the cycle
 (P2) amplitude-parametrised family: for crest amplitude a on a grid, solve for (y0, rho); is rho(a) strictly
      increasing on the whole ball (one cycle per rho)?  m_t < 1 everywhere?
 (P3) planar squeeze at mu = 0.02: crest sequences started near E and at the ball boundary; monotone, same limit?
"""
import sys, time
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root
sys.path.insert(0, '/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03/work/lanes/hamiltonian/verify')
sys.path.insert(0, '/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03/work')
from whref import cycle, crest_map, rhs
from todaq import level, hstar, Phi, turning

T00 = time.time()
FAIL = []

# (P0) ---------------------------------------------------------------------------------------------------------
import whlib
print("(P0) referee finder vs whlib.find_cycle")
for mu, rho in [(1.0, 2.0), (0.1, 1.5), (0.3, 4.0)]:
    c = cycle(mu, rho)
    w = whlib.find_cycle(mu, rho * mu * (1 + mu))
    print(f"  mu={mu} rho={rho}:  T {c['T']:.9f} / {w['T']:.9f}   m_t {c['m_t']:.9e} / {w['m_t']:.9e}   U0 {c['U0']:.9f} / {w['U0']:.9f}")
    if abs(c['T'] - w['T']) > 1e-6 or abs(c['m_t'] - w['m_t']) > 1e-6 * max(1, w['m_t']):
        FAIL.append(f"P0 {mu},{rho}")


def wseries(U, P, eps, r):
    E = np.exp(U); G = E - 1
    return (-G) + eps * (G * P) + eps ** 2 * (-r * G - E * P ** 2 + G ** 2) + eps ** 3 * (r * E * P + E * P ** 3 - 4 * G * E * P - G ** 2 * P)


# (P1) ---------------------------------------------------------------------------------------------------------
print("\n(P1) fixed r, mu -> 0")
print("# r     mu         eps       h*(r)     h_crest   (h_crest-h*)/eps^2   logmt/eps^3   pred=-J0Rsel'   ratio     (ratio-1)/eps^2   sup|Q-w_ser|/eps^4   res")
for r in [0.3, 1.0, 3.0]:
    hs = hstar(r)
    d = level(hs)
    pred = -d["J0"] * d["Rselp"]
    for mu in [1e-2, 2.5e-3, 6.25e-4, 1.5625e-4]:
        eps = np.sqrt(mu) / (1 + mu)
        rho = 1 + r * eps ** 2
        k = rho * mu * (1 + mu)
        a0 = d["up"]
        c = cycle(mu, rho, guess=(a0, k + mu * rho * np.expm1(a0)), warm=6)
        hc = Phi(c["U0"])
        # Step 1 check: Q - w_series along the cycle
        sol = solve_ivp(rhs(mu, k), [0, c["T"]], [c["U0"], 0.0, c["y0"]], method="DOP853", rtol=1e-12, atol=1e-14,
                        t_eval=np.linspace(0, c["T"], 2001))
        Uc, pc, yc = sol.y
        Pc = pc / np.sqrt(mu)
        Qc = (k - mu * pc - yc) / mu
        Dser = np.max(np.abs(Qc - wseries(Uc, Pc, eps, r)))
        ratio = c["logmt"] / eps ** 3 / pred
        print(f"{r:<6g} {mu:<10.4g} {eps:<9.5f} {hs:<9.6f} {hc:<9.6f} {(hc-hs)/eps**2:<20.4f} {c['logmt']/eps**3:<13.6f} {pred:<15.6f} "
              f"{ratio:<9.6f} {(ratio-1)/eps**2:<17.3f} {Dser/eps**4:<20.3f} {c['res']:.1e}", flush=True)
        if not (c["m_t"] < 1):
            FAIL.append(f"P1 m_t>=1 r={r} mu={mu}")
print(f"   [{time.time()-T00:.0f}s]")


# (P2) ---------------------------------------------------------------------------------------------------------
def amp_cycle(mu, a, rho_guess):
    kH = mu * (1 + mu)

    def F(v):
        s, rho = v
        k = rho * kH
        y0 = k + mu * rho * s
        U1, y1, T, _, _ = crest_map(mu, k, a, y0, var=False)
        return [U1 - a, (y1 - y0) / kH]
    sol = root(F, [np.expm1(a), rho_guess], method="hybr", tol=1e-13)
    s, rho = sol.x
    k = rho * kH
    y0 = k + mu * rho * s
    U1, y1, T, M, Umin = crest_map(mu, k, a, y0)
    return dict(rho=rho, y0=y0, T=T, m_t=float(np.max(np.linalg.eigvals(M).real)), res=float(np.linalg.norm(sol.fun)), ok=sol.success)


print("\n(P2) amplitude family: r(a) = (rho(a)-1)/eps^2, compared with Rsel(Phi(a)); monotone?  m_t<1?")
for mu in [0.02, 1e-3]:
    eps = np.sqrt(mu) / (1 + mu)
    agrid = np.concatenate((np.arange(0.1, 1.0, 0.1), np.arange(1.0, 5.01, 0.25)))
    prev = -np.inf
    mono = True
    allstab = True
    rows = []
    for a in agrid:
        h = Phi(a)
        Rs = level(h)["Rsel"]
        try:
            c = amp_cycle(mu, a, 1 + eps ** 2 * Rs)
        except Exception as e:
            print(f"   mu={mu} a={a:.2f} FAILED {e}")
            continue
        rnum = (c["rho"] - 1) / eps ** 2
        if rnum <= prev:
            mono = False
        prev = rnum
        if not c["m_t"] < 1:
            allstab = False
        rows.append((a, h, rnum, Rs, c["m_t"], c["res"]))
    print(f"  mu={mu} eps={eps:.5f}:  a, h=Phi(a), r_num, Rsel(h), r_num/Rsel, m_t, residual")
    for row in rows[::2]:
        print(f"    {row[0]:<5.2f} {row[1]:<9.4f} {row[2]:<11.6f} {row[3]:<11.6f} {row[2]/row[3]:<8.5f} {row[4]:<12.9f} {row[5]:.1e}")
    print(f"    r_num strictly increasing over {len(rows)} amplitudes a in [0.1,5] (h up to {Phi(5.0):.1f}): {mono};  all m_t<1: {allstab}")
    if not mono:
        FAIL.append(f"P2 non-monotone mu={mu}")
    if not allstab:
        FAIL.append(f"P2 m_t>=1 mu={mu}")
print(f"   [{time.time()-T00:.0f}s]")

# (P3) ---------------------------------------------------------------------------------------------------------
print("\n(P3) planar squeeze, mu=0.02: crest sequences from inside (near E) and from the ball boundary")
mu = 0.02
eps = np.sqrt(mu) / (1 + mu)
for rho, a_out in [(1.01, 3.0), (1.04, 3.5)]:
    r = (rho - 1) / eps ** 2
    k = rho * mu * (1 + mu)
    c = cycle(mu, rho, guess=(turning(hstar(r))[1], k + mu * rho * np.expm1(turning(hstar(r))[1])), warm=5)
    out = {}
    for name, a0 in [("inside", 0.02), ("outside", a_out)]:
        U0, y0 = a0, k + mu * rho * np.expm1(a0)
        seq = []
        for n in range(4000):
            U0, y0, T, _, _ = crest_map(mu, k, U0, y0, var=False)
            seq.append(U0)
            if n > 5 and abs(seq[-1] - seq[-2]) < 1e-11:
                break
        seq = np.array(seq)
        dif = np.diff(seq[3:])                      # discard the fast transient onto the slow manifold
        monotone = bool(np.all(dif > 0)) if name == "inside" else bool(np.all(dif < 0))
        # ignore the last few differences, which are at round-off level
        big = np.abs(dif) > 1e-9
        monotone_big = bool(np.all(dif[big] > 0)) if name == "inside" else bool(np.all(dif[big] < 0))
        out[name] = (seq[3], seq[-1], len(seq), monotone_big)
        print(f"  rho={rho} r={r:.3f} {name:<8}: crest after transient {seq[3]:.6f} (H0={Phi(seq[3]):.4f}) -> limit {seq[-1]:.9f} "
              f"after {len(seq)} returns; monotone (|step|>1e-9): {monotone_big}")
        if not monotone_big:
            FAIL.append(f"P3 monotone rho={rho} {name}")
    print(f"     Newton cycle crest {c['U0']:.9f}, m_t={c['m_t']:.6f};  |inside-limit - outside-limit| = {abs(out['inside'][1]-out['outside'][1]):.2e}")
    if abs(out['inside'][1] - out['outside'][1]) > 1e-6 or abs(out['inside'][1] - c['U0']) > 1e-6:
        FAIL.append(f"P3 limits differ rho={rho}")
print(f"   [{time.time()-T00:.0f}s]")
print("\nFAILURES:", FAIL if FAIL else "none")
