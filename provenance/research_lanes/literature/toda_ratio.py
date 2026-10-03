#!/usr/bin/env python3
"""Lane 'literature' numerical check (floating point, single process).

Claim tested (leading-order singular limit mu -> 0, rho fixed):
  slow time tau = sqrt(k) t, q = p/sqrt(k):   U_tautau + (mu/eps - eps e^U) U_tau + e^U - 1 = O(eps^2),  eps = sqrt(k)
  => cycles shadow level sets of H = q^2/2 + e^U - U - 1, selected by  R(h) := <e^U q^2>_h / <q^2>_h = mu A / k = 1/rho
     (the exact identity <x p^2> = mu A <p^2> of the brief).
  Zhang Zhifen's condition for the truncated Lienard equation: f/g = (1/rho' - e^U)/(e^U - 1) increasing on U<0 and U>0.
Test 1: R(h) on a grid of h (quadrature).  OUTCOME: R(h) == 1 identically (analytic reason: on H-level sets
        q_tau = -(e^U - 1), so  oint (e^U - 1) q^2 dtau = -oint q^2 dq = 0).  Hence the first-order Melnikov function is
        eps (1 - mu/k) oint q^2 dtau > 0 for every level when rho(1+mu) > 1: the truncated Lienard equation has NO limit
        cycle, and the proposed reduction is REFUTED as a route to uniqueness at fixed rho > 1.
Test 2: actual WH cycles at small mu: relative spread of H along the cycle.  OUTCOME: H is far from constant
        (spread/mean between 1 and 7; mean H ~ 1/mu), so real cycles are not near-Hamiltonian at fixed rho > 1.
"""
import sys, numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
sys.path.insert(0, '/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03/work')

def turning(h):
    F = lambda U: np.exp(U) - U - 1 - h
    Up = brentq(F, 1e-14, 50 + h)          # e^U - U - 1 = h, U>0
    Um = brentq(F, -(h + 1.0) - 1e-9, -1e-14)  # U<0 ; F(-(h+1)) = e^{-(h+1)} > 0
    return Um, Up

def R(h):
    Um, Up = turning(h)
    # substitute U = mid + half*sin(th) to kill endpoint sqrt singularities
    mid, half = 0.5*(Um+Up), 0.5*(Up-Um)
    def integ(w):
        def g(th):
            U = mid + half*np.sin(th)
            q2 = 2*(h - (np.exp(U) - U - 1))
            q = np.sqrt(max(q2, 0.0))
            return w(U)*q*half*np.cos(th)
        return quad(g, -np.pi/2, np.pi/2, limit=400, epsabs=1e-13, epsrel=1e-11)[0]
    num = integ(lambda U: np.exp(U))   # int e^U q dU  (prop. to <e^U q^2> T)
    den = integ(lambda U: 1.0)         # int q dU
    return num/den

if __name__ == '__main__':
    hs = np.concatenate(([1e-4, 1e-3, 1e-2], np.linspace(0.05, 2, 40), np.linspace(2.5, 60, 60)))
    Rs = np.array([R(h) for h in hs])
    d = np.diff(Rs)
    print('Test 1: R(h) on %d points, h in [%.0e, %.0f]' % (len(hs), hs[0], hs[-1]))
    print('  R(first) = %.8f  R(last) = %.3e   max increment = %.3e  (monotone decreasing: %s)' % (Rs[0], Rs[-1], d.max(), bool((d < 0).all())))
    for h in [1e-3, 0.1, 0.5, 1, 2, 5, 10, 30, 60]:
        print('    h=%-6g R=%.8f   1/R=%.5f' % (h, R(h), 1/R(h)))
    # small-h expansion check: R = 1 - c h + ...
    print('  (1-R(h))/h at h=1e-3,1e-2: %.6f %.6f' % ((1-R(1e-3))/1e-3, (1-R(1e-2))/1e-2))
    # Zhang ratio: d/dU [(1/rho - e^U)/(e^U - 1)] = e^U (1 - 1/rho)/(e^U - 1)^2  > 0 for rho > 1 (analytic; evaluated here)
    for rho in [1.01, 2.0, 50.0]:
        U = np.concatenate((np.linspace(-30, -1e-3, 4000), np.linspace(1e-3, 30, 4000)))
        dr = np.exp(U)*(1 - 1/rho)/np.expm1(U)**2
        print('  Zhang d(f/g)/dU > 0 on grid: %s (rho=%g)' % (bool((dr > 0).all()), rho))
    # Test 2
    from whlib import find_cycle, cycle_profile, kH
    print('Test 2: actual WH cycles (whlib.find_cycle), H along cycle in scaled variables q = p/sqrt(k)')
    print('  mu     rho    k        T        Hmean     spread(H)/Hmean   R(Hmean)   1/rho     1/(rho(1+mu))  lambda_c')
    for mu in [0.2, 0.05, 0.02]:
        for rho in [1.5, 3.0, 8.0]:
            k = rho*kH(mu)
            try:
                c = find_cycle(mu, k)
                pr = cycle_profile(mu, k, c)
                U, p = pr['U'], pr['p']
                H = 0.5*p**2/k + np.exp(U) - U - 1
                Hm = np.trapezoid(H, pr['t'])/(pr['t'][-1]-pr['t'][0])
                print('  %-6g %-6g %-8.5f %-8.3f %-9.5f %-16.4f %-10.6f %-9.6f %-13.6f %.3e' % (mu, rho, k, c['T'], Hm, (H.max()-H.min())/Hm, R(Hm), 1/rho, 1/(rho*(1+mu)), c['lambda_c']))
            except Exception as e:
                print('  %-6g %-6g FAILED: %s' % (mu, rho, repr(e)[:120]))
