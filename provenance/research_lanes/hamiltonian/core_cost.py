#!/usr/bin/env python3
"""Floating-point sizing of the computer-assisted core (PROGRAM.md, Section 4).  NOT rigorous.

For sample (mu, rho): locate the cycle crest z* = (U*, y*), return derivative M, weak multiplier m_t, and
  kappa   = second derivative of the reduced crest map along the weak direction (central difference),
  delta*  = (1-m_t)/|kappa|        radius (in U) of the window where 'contraction' can be certified,
  dU*/drho, dU*/dmu                 (finite differences)  -> admissible parameter box widths
       d_rho = delta*/(2|dU*/drho|),  d_mu = delta*/(2|dU*/dmu|),
  N_state = 2*ln(L/((1-m_t)delta*)) / ln(1+(1-m_t)/2)   geometric cover of the crest curve outside the window
       (both sides), L = a-priori crest range log(M0),
  density = N_state/(d_rho*d_mu)    box evaluations per unit (rho, mu) area.
"""
import sys, json
import numpy as np
sys.path.insert(0, '/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03/work')
from whlib import find_cycle, crest_return, kH, equilibrium_data
from toda import level_data, h_star, turning_points


def guess_from_theory(mu, rho):
    eps = np.sqrt(mu) / (1 + mu)
    r = (rho - 1) / eps ** 2
    h = h_star(r)
    um, up = turning_points(h)
    k = rho * kH(mu)
    return (up, k + mu * rho * np.expm1(up))


def cyc(mu, rho, guess=None):
    k = rho * kH(mu)
    eps = np.sqrt(mu) / (1 + mu)
    r = (rho - 1) / eps ** 2
    if guess is None and r < 6 and mu <= 0.25:
        guess = guess_from_theory(mu, rho)
    return find_cycle(mu, k, guess=guess, warm=(0 if guess is not None else 120), newton=40)


def analyse(mu, rho):
    k = rho * kH(mu)
    c = cyc(mu, rho)
    z = np.array([c["U0"], c["y0"]])
    M = c["M"]
    w, V = np.linalg.eig(M)
    i = int(np.argmax(w.real))
    v = np.real(V[:, i]); v = v / v[0]
    wl, Vl = np.linalg.eig(M.T)
    j = int(np.argmax(wl.real))
    l = np.real(Vl[:, j]); l = l / (l @ v)
    mt = c["m_t"]
    s0 = 0.02 * min(1.0, c["U0"])
    vals = []
    for s in (+s0, -s0):
        U1, y1, T1, _ = crest_return(mu, k, z[0] + s * v[0], z[1] + s * v[1], with_derivative=False)
        vals.append(l @ (np.array([U1, y1]) - z))
    kappa = (vals[0] + vals[1]) / s0 ** 2
    dstar = (1 - mt) / abs(kappa) if kappa != 0 else np.inf
    dstar = min(dstar, 0.5 * c["U0"])
    drho = 1e-4 * (rho - 1)
    c2 = find_cycle(mu, (rho + drho) * kH(mu), guess=(c["U0"], c["y0"] * (rho + drho) / rho), warm=0, newton=40)
    dUdrho = (c2["U0"] - c["U0"]) / drho
    dmu = 1e-4 * mu
    c3 = find_cycle(mu + dmu, rho * kH(mu + dmu), guess=(c["U0"], c["y0"] * kH(mu + dmu) / kH(mu)), warm=0, newton=40)
    dUdmu = (c3["U0"] - c["U0"]) / dmu
    lam = k / mu
    M0 = (lam + mu) * (lam + 1) / mu
    L = np.log(M0)
    d_rho = dstar / (2 * abs(dUdrho))
    d_mu = dstar / (2 * abs(dUdmu))
    N_state = 2 * np.log(L / ((1 - mt) * dstar)) / np.log1p((1 - mt) / 2)
    ed = equilibrium_data(mu, k)
    mE = np.exp(2 * np.pi * ed["lam_u"].real / abs(ed["lam_u"].imag))
    return dict(mu=mu, rho=rho, U=c["U0"], T=c["T"], one_minus_mt=1 - mt, kappa=kappa, dstar=dstar, dUdrho=dUdrho,
                dUdmu=dUdmu, d_rho=d_rho, d_mu=d_mu, N_state=N_state, density=N_state / (d_rho * d_mu), L=L, mE=mE,
                res=max(c["residual"], c2["residual"], c3["residual"]))


if __name__ == "__main__":
    mus = [1.0, 0.5, 0.2, 0.1, 0.05]
    drs = [0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 4.0, 9.0, 19.0]
    rows = []
    print("# mu    rho     r=(rho-1)/eps^2  U*      T       1-m_t     kappa      delta*    dU/drho   dU/dmu    d_rho     d_mu      N_state  density   m_E     res")
    for mu in mus:
        eps = np.sqrt(mu) / (1 + mu)
        for dr in drs:
            rho = 1 + dr
            try:
                d = analyse(mu, rho)
            except Exception as e:
                print(f"{mu:<6g} {rho:<7g} FAILED {type(e).__name__}: {str(e)[:60]}", flush=True)
                continue
            rows.append(d)
            print(f"{mu:<6g} {rho:<7g} {dr/eps**2:<16.4g} {d['U']:<7.3f} {d['T']:<7.1f} {d['one_minus_mt']:<9.3e} {d['kappa']:<10.3e} "
                  f"{d['dstar']:<9.2e} {d['dUdrho']:<9.3g} {d['dUdmu']:<9.3g} {d['d_rho']:<9.2e} {d['d_mu']:<9.2e} "
                  f"{d['N_state']:<8.0f} {d['density']:<9.2e} {d['mE']:<7.3f} {d['res']:.0e}", flush=True)
    json.dump(rows, open("core_cost.json", "w"), indent=1)
