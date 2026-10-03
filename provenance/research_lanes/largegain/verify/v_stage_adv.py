#!/usr/bin/env python3
"""Referee script G: adversarial search for a Stage-1 counterexample among non-periodic lobes:
minimise p_k/k over lobes with k>=1000 and T_- >= 5.89 (differential evolution over mu, k, P_d, w_d),
with and without the cycle-only hypothesis Lambda_E(s_w) > 0, i.e. P_w < (gamma_E - mu) v_w."""
import numpy as np
from scipy.optimize import differential_evolution, brentq
import v_stage as S


def quantities(x):
    mu, k, Pd, wd = np.exp(x[0]), np.exp(x[1]), np.exp(x[2]), np.exp(x[3]) * np.exp(x[1])
    lb = S.lobe(mu, k, Pd, wd)
    if lb is None:
        return None
    if not lb["deep"]:
        sol = lb["sols"][0][1]
        Tm = sol.t[-1]
        return dict(mu=mu, k=k, Pd=Pd, wd=wd, Tm=Tm, pk=sol.y[1, -1], deep=False, lam=None)
    k_ = k
    if lb["w0"] > 0:
        tw = np.log((lb["w0"] + k_) / k_)
        vw = -lb["U_"](tw); Pw = -lb["p_"](tw)
    else:
        g = lambda t: float(lb["sol"].sol(t)[2]) - k_
        sw = brentq(g, 0, lb["t1"], xtol=1e-14)
        vw = -float(lb["sol"].sol(sw)[0]); Pw = -float(lb["sol"].sol(sw)[1])
    r = np.roots([1, 1 + mu, mu, k]); gE = -min(r.real[np.abs(r.imag) < 1e-9])
    return dict(mu=mu, k=k, Pd=Pd, wd=wd, Tm=lb["Tm"], pk=lb["sol2"].y[1, -1], deep=True, lam=bool(Pw < (gE - mu) * vw))


def make_obj(require_lam):
    def obj(x):
        q = quantities(x)
        if q is None:
            return 10.0
        pen = 0.0
        if q["Tm"] < 5.89:
            pen += 5.0 + (5.89 - q["Tm"])
        if require_lam and q["lam"] is False:
            pen += 5.0
        if require_lam and q["lam"] is None:
            pen += 5.0
        return q["pk"] / q["k"] + pen
    return obj


bounds = [(np.log(0.005), 0.0), (np.log(1000.0), np.log(1e5)), (np.log(1e-4), np.log(1e7)), (np.log(1e-4), np.log(1e4))]
for req in (True, False):
    best = None
    for seed in (1, 2, 3):
        res = differential_evolution(make_obj(req), bounds, seed=seed, maxiter=60, popsize=20, tol=1e-8, polish=False)
        if best is None or res.fun < best.fun:
            best = res
    q = quantities(best.x)
    print(f"require Lambda_E(s_w)>0: {req} -> min p_k/k = {best.fun:.5f} at mu={q['mu']:.4g} k={q['k']:.5g} P_d={q['Pd']:.4g} w_d/k={q['wd']/q['k']:.4g} T_-={q['Tm']:.5f} deep={q['deep']} LambdaE_ok={q['lam']}")
    print("   claim p_k >= 0.374 k:", "consistent" if best.fun >= 0.374 else "VIOLATED")
