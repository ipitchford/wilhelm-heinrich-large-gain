#!/usr/bin/env python3
"""Integrate the box-evaluation density of core_cost.json over the core
   C = { mu0 <= mu <= 1,  1 + r_c*eps(mu)^2 <= rho <= rho1 }   (idealised count; see PROGRAM.md Section 4.6)."""
import json
import numpy as np

rows = json.load(open("core_cost.json"))
mus = sorted({r["mu"] for r in rows})
tab = {}
for mu in mus:
    rs = sorted([r for r in rows if r["mu"] == mu], key=lambda r: r["rho"])
    tab[mu] = rs
# mu = 1: dU*/dmu vanishes by the reciprocal symmetry, so d_mu is meaningless there; borrow d_mu from mu = 0.5
for r1, r5 in zip(tab[1.0], tab[0.5]):
    r1["d_mu"] = r5["d_mu"]
    r1["density"] = r1["N_state"] / (r1["d_rho"] * r1["d_mu"])


def per_unit_mu(mu, r_c, rho1):
    rs = tab[mu]
    x = np.log([r["rho"] - 1 for r in rs]); yv = np.log([r["density"] for r in rs])
    eps2 = mu / (1 + mu) ** 2
    lo, hi = np.log(r_c * eps2), np.log(rho1 - 1)
    g = np.linspace(lo, hi, 4000)
    # log-log interpolation, power-law extrapolation below the first grid point with the local slope
    slope0 = (yv[1] - yv[0]) / (x[1] - x[0])
    yy = np.where(g < x[0], yv[0] + slope0 * (g - x[0]), np.interp(g, x, yv))
    return float(np.trapezoid(np.exp(yy) * np.exp(g), g)), slope0


for r_c in (0.2, 0.05):
    for rho1 in (20.0,):
        print(f"# collar r_c = {r_c}, rho1 = {rho1}")
        vals = []
        for mu in mus:
            v, s = per_unit_mu(mu, r_c, rho1)
            vals.append(v)
            eps2 = mu / (1 + mu) ** 2
            print(f"   mu={mu:<5g} rho_c-1={r_c*eps2:.4f}  ideal evaluations per unit mu = {v:.3e}   (near-Hopf slope {s:.2f})")
        lm = np.log(mus); lv = np.log(vals)
        for mu0 in (0.05, 0.1, 0.2):
            g = np.linspace(np.log(mu0), 0, 2000)
            tot = np.trapezoid(np.exp(np.interp(g, lm, lv)) * np.exp(g), g)
            print(f"   mu0={mu0:<5g}: ideal total = {tot:.3e} box evaluations;  x1e2 = {tot*1e2:.1e}, x1e3 = {tot*1e3:.1e}, x1e4 = {tot*1e4:.1e}")
