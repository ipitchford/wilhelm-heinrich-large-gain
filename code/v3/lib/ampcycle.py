#!/usr/bin/env python3
# --- v3.0 bundle copy of research file `lanes/hamiltonian/ampcycle.py`.
# --- Changes from the original: absolute import path removed.
"""Amplitude-fixed shooting for Wilhelm-Heinrich cycles (floating point).
Given mu and crest amplitude a = U_max, solve for (y0, k) with crest return (a,0,y0) -> (a,0,y0).
Uses whlib.crest_return (not modified)."""
import sys
import numpy as np
from whlib import crest_return, kH
from toda import level_data, Phi


def eps_of_mu(mu):
    return np.sqrt(mu) / (1.0 + mu)


def amp_cycle(mu, a, guess=None, tol=1e-10, maxit=25, verbose=False):
    eps = eps_of_mu(mu)
    if guess is None:
        h = Phi(a)
        d = level_data(h)
        r = (2.0 / 3.0) * d["J1"] / d["I"]
        rho = 1.0 + eps ** 2 * r
        k = rho * kH(mu)
        y0 = k + mu * rho * np.expm1(a)
    else:
        y0, k = guess
    for it in range(maxit):
        U1, y1, T, M = crest_return(mu, k, a, y0)
        res = np.array([U1 - a, (y1 - y0)])
        if verbose:
            print(it, res, k / kH(mu))
        if abs(res[0]) + abs(res[1]) / k < tol:
            break
        dk = 1e-6 * k
        U1b, y1b, _, _ = crest_return(mu, k + dk, a, y0, with_derivative=False)
        dFdk = np.array([(U1b - U1) / dk, (y1b - y1) / dk])
        Jm = np.array([[M[0, 1], dFdk[0]], [M[1, 1] - 1.0, dFdk[1]]])
        step = np.linalg.solve(Jm, -res)
        y0 += step[0]
        k += step[1]
    U1, y1, T, M = crest_return(mu, k, a, y0)
    ev = np.linalg.eigvals(M)
    m_t = float(np.max(ev.real))
    return dict(mu=mu, a=a, k=k, rho=k / kH(mu), y0=y0, T=T, M=M, m_t=m_t, logmt=float(np.log(m_t)),
                residual=float(abs(U1 - a) + abs(y1 - y0) / k), U0=a)


if __name__ == "__main__":
    print("# test of third-order selection  r = (rho-1)/eps^2 -> (2/3)<p^4>/<p^2>(h),  lambda_c*T/eps^3 -> T0*(2/3)(<p^4>-3<p^2>^2)/<p^2>")
    print("# mu       eps      a      h=Phi(a)  r_num      r_pred     r_num/r_pred  logmt/eps^3   pred        ratio    T*sqrt(mu)  T0       res")
    for a in [0.3, 1.0, 2.0, 3.0]:
        h = Phi(a)
        d = level_data(h)
        r_pred = (2.0 / 3.0) * d["J1"] / d["I"]
        lam_pred = d["T0"] * (2.0 / 3.0) * (d["p4"] - 3 * d["p2"] ** 2) / d["p2"]
        for mu in [3e-2, 1e-2, 3e-3, 1e-3, 3e-4, 1e-4]:
            eps = eps_of_mu(mu)
            try:
                c = amp_cycle(mu, a)
            except Exception as e:
                print(f"{mu:<8g} a={a} FAILED {e}")
                continue
            r_num = (c["rho"] - 1) / eps ** 2
            print(f"{mu:<8g} {eps:<8.5f} {a:<6g} {h:<9.5f} {r_num:<10.6f} {r_pred:<10.6f} {r_num/r_pred:<13.6f} "
                  f"{c['logmt']/eps**3:<13.6f} {lam_pred:<11.6f} {c['logmt']/eps**3/lam_pred:<8.5f} "
                  f"{c['T']*np.sqrt(mu):<11.6f} {d['T0']:<8.5f} {c['residual']:.1e}", flush=True)
