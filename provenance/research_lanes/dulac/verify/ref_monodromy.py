#!/usr/bin/env python3
"""Riccati-free confirmation of lambda_c at large rho, where double-precision monodromy fails (the lane's own
verify_theorem.out shows lam_M != <eta>-mu at rho=20).  The crest return map is differentiated by forward
differences in 400-bit arithmetic (increment 1e-18, integrator tolerance tightened to 1e-45), the weak multiplier
m_t is the larger eigenvalue of the 2x2 section derivative, and log(m_t)/T is compared with <eta>-mu, -<pK/V>,
-<pW/Q> from the case_*.json written by ref_taylor.py.

usage: ref_monodromy.py MU RHO"""
import sys, os, math, json, time
from fractions import Fraction
import numpy as np
from flint import arb, fmpq
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ref_taylor as rt
import refcyc

rt.LOGTOL = -45.0 * math.log(10.0)
M = rt.M


def ret(mu, k, U0, a0, hmax):
    y0 = M(k * a0.exp())
    st = [U0, arb(0), y0, M(k * U0.exp()), arb(1), arb(1)]      # Riccati slots are passengers here
    st1, T, n = rt.integrate_period(mu, k, st, hmax)
    return st1[0], (st1[2] / k).log(), T, n


def main(mu_fr, rho_fr, guess):
    t0 = time.time()
    mu = M(arb(fmpq(mu_fr.numerator, mu_fr.denominator)))
    k_fr = rho_fr * mu_fr * (1 + mu_fr)
    k = M(arb(fmpq(k_fr.numerator, k_fr.denominator)))
    c = refcyc.find(float(mu_fr), float(k_fr), U0=guess)
    hmax = c['T'] / 80.0
    U0 = arb(c['U0']); a0 = arb(c['a0'])
    Jm = c['J'] - np.eye(2)
    for it in range(4):                                   # polish the crest as in ref_taylor
        U1, a1, T, n = ret(mu, k, U0, a0, hmax)
        r = np.array([float((U1 - U0).mid()), float((a1 - a0).mid())])
        print(f'   polish {it}: residual {np.abs(r).max():.2e} steps={n} [{time.time()-t0:.0f}s]', flush=True)
        if np.abs(r).max() < 1e-28:
            break
        d = np.linalg.solve(Jm, -r)
        U0 = M(U0 + arb(float(d[0]))); a0 = M(a0 + arb(float(d[1])))
    dl = arb(10) ** (-18)
    Ub, ab, T, _ = ret(mu, k, U0, a0, hmax)
    Uu, au, _, _ = ret(mu, k, M(U0 + dl), a0, hmax)
    Ua, aa, _, _ = ret(mu, k, U0, M(a0 + dl), hmax)
    J11 = (Uu - Ub) / dl; J21 = (au - ab) / dl; J12 = (Ua - Ub) / dl; J22 = (aa - ab) / dl
    tr = J11 + J22; det = J11 * J22 - J12 * J21
    disc = tr * tr - 4 * det
    m_t = (tr + disc.sqrt()) / 2
    m_s = det / m_t
    lam = m_t.log() / T
    A = 1 + mu
    out = dict(mu=float(mu_fr), rho=float(rho_fr), T=float(T.mid()), entries=[rt.fl(v) for v in (J11, J12, J21, J22)],
               m_t=rt.fl(m_t), log10_m_t=rt.lg(m_t) / math.log(10), lam_monodromy=rt.fl(lam),
               log_det=rt.lg(det), minus_A_T=rt.fl(-A * T), log_m_s=rt.lg(m_s))
    fn = os.path.join(os.path.dirname(os.path.abspath(__file__)), f'case_mu{float(mu_fr):g}_rho{float(rho_fr):g}.json')
    if os.path.exists(fn):
        r = json.load(open(fn))
        out['lam_eta'] = r['lambda_c: <eta>-mu']; out['lam_V'] = r['lambda_c: -<pK/V> (GL)']; out['lam_Q'] = r['lambda_c: -<pW/Q> (GL)']
        out['rel_diff_monodromy_vs_eta'] = abs(out['lam_monodromy'] - out['lam_eta']) / abs(out['lam_eta'])
    print(json.dumps(out, indent=1))
    print(f'   entrywise positive M: {bool(J11 > 0 and J12 > 0 and J21 > 0 and J22 > 0)}   log det M + A T = {out["log_det"] - out["minus_A_T"]:+.3e}   [{time.time()-t0:.0f}s]')


if __name__ == '__main__':
    main(Fraction(sys.argv[1]), Fraction(sys.argv[2]), float(sys.argv[3]) if len(sys.argv) > 3 else 1.0)
