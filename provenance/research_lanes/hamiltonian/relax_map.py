#!/usr/bin/env python3
"""Singular (mu -> 0, rho fixed) relaxation limit: 1-D return map  F_lam = G_lam o S.
Pulse (universal, v = u + log lam, zh = lam*z):  v' = p, p' = -zh, zh' = e^v - zh,
   from (v=-inf, p=q>0, zh=0) to (v=-inf, p=-S(q), zh=0).
Slow phase (s = mu t, V = mu U): p' = lam - p, V' = p, from (V=0, p=-qm) back to V=0:  p_plus = G_lam(qm).
Floating point; a diagnostic of the regime structure, not a proof."""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

V0 = -40.0


def pulse(q, with_derivative=False):
    """Return S(q) (and S'(q) by central difference if requested)."""
    def f(t, s):
        v, p, z = s
        return [p, -z, np.exp(v) - z]
    e0 = np.exp(V0)
    s0 = [V0, q - e0 / (q * (1 + q)), e0 / (1 + q)]
    def ev(t, s):
        return s[0] - V0
    ev.terminal = True; ev.direction = -1
    a = solve_ivp(f, [0, 1e-3 / max(q, 1e-3)], s0, method="DOP853", rtol=1e-12, atol=1e-18)
    b = solve_ivp(f, [a.t[-1], 400.0 / q + 200.0], a.y[:, -1], method="DOP853", rtol=1e-12, atol=1e-18, events=ev)
    assert b.status == 1
    v, p, z = b.y_events[0][0]
    pm = p - z - np.exp(v) / abs(p) * 1.0 / 1.0  # remaining int of zh: z_end + int e^v
    S = -pm
    if with_derivative:
        dq = 1e-5 * q
        return S, (pulse(q + dq) - pulse(q - dq)) / (2 * dq)
    return S


def slow(lam, qm):
    """G_lam(qm), its derivative, return time S, and V_min."""
    L = lam + qm
    g = lambda s: lam * s - L * (-np.expm1(-s))
    hi = 2 * qm / lam + 1.0
    while g(hi) < 0:
        hi *= 2
    s = brentq(g, 1e-14, hi, xtol=1e-15, rtol=1e-14)
    pp = lam * s - qm
    dG = lam * (-np.expm1(-s)) / pp - 1.0
    Vmin = lam * np.log(L / lam) - qm
    return pp, dG, s, Vmin


def fixed_point(lam):
    def F(q):
        return slow(lam, pulse(q))[0] - q
    qs = np.geomspace(1e-3, lam * (1 - 1e-6), 60)
    vals = [F(q) for q in qs]
    roots = []
    for i in range(len(qs) - 1):
        if vals[i] * vals[i + 1] < 0:
            roots.append(brentq(F, qs[i], qs[i + 1], xtol=1e-13, rtol=1e-12))
    out = []
    for q in roots:
        S, dS = pulse(q, with_derivative=True)
        pp, dG, s, Vmin = slow(lam, S)
        out.append(dict(q=q, S=S, mult=dG * dS, Tmu=s, Vmin=Vmin))
    return out


if __name__ == "__main__":
    print("# small-q check of pulse map: S(q) ~ q + (2/3) q^2")
    for q in [0.01, 0.03, 0.1, 0.3, 1.0, 3.0]:
        S = pulse(q)
        print(f"  q={q:<5g} S={S:.8f}  (S-q)/q^2={(S-q)/q**2:.5f}")
    print("# singular fixed points: rho, #roots, q*, S(q*), multiplier (-> m_t), T*mu, mu*U_min")
    for lam in [1.02, 1.1, 1.5, 3.0, 8.0, 30.0]:
        fp = fixed_point(lam)
        for d in fp:
            print(f"  rho={lam:<5g} nroots={len(fp)} q*={d['q']:.6f} S={d['S']:.6f} mult={d['mult']:.6f} T*mu={d['Tmu']:.5f} muUmin={d['Vmin']:.5f}")
