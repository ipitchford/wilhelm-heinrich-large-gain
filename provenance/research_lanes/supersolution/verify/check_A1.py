#!/usr/bin/env python3
"""Referee check 5 (Lemma A1): near-Hopf cycles.
 (a) regress the TRUE (eta-mu, zeta-mu) on (U,p): coefficients must tend to
       (mu A/D, -2mu/D) and (mu(mu^2+4mu+1)/D, -mu A/D),  D = mu^2+6mu+1.
 (b) scores of affine state trials  eta~ = mu + c1 U + c2 p, zeta~ = mu + d1 U + d2 p:
       exact jet -> S = O(eps^2);  perturbed jet -> S/eps -> <m_1> < 0 with m_1 from the lemma's first-order formula
       (computed independently on the Hopf circle);  omega^Y and constants likewise negative."""
import sys, time
import numpy as np
ROOT = '/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03'
sys.path.insert(0, ROOT + '/work')
from whlib import find_cycle, kH
from vprofile import profile, load_cycles
C = load_cycles()

def score(mu, k, T, tt, U, p, y, e, z, De, Dz):
    x = k*np.exp(U); A = 1+mu
    g = np.minimum.reduce([mu - e, De/e + 1 - z/e, Dz/z + A - x/z])
    return np.trapezoid(g, tt)/T

def m1_mean(mu, c1, c2, d1, d2, n=20000):
    """<m_1> of Lemma A1 for the trial jet (per unit amplitude) on the Hopf circle U = cos(w t), p = -w sin(w t)."""
    A = 1+mu; w = np.sqrt(mu); ph = np.linspace(0, 2*np.pi, n, endpoint=False)
    u = np.cos(ph); pp = -w*np.sin(ph); up = pp; ppd = -mu*u
    a = c1*u + c2*pp; b = d1*u + d2*pp
    ad = c1*up + c2*ppd; bd = d1*up + d2*ppd
    m = np.minimum.reduce([-a, (ad + a - b)/mu, (bd + A*b - mu*A*u)/mu])
    return m.mean()

print(f"{'mu':>5} {'rho-1':>7} {'eps':>9} {'lam_c':>11} | {'aU fit':>9} {'claim':>9} {'ap fit':>9} {'claim':>9} {'bU fit':>9} {'claim':>9} {'bp fit':>9} {'claim':>9} | {'S_true/e^2':>10} {'S_jet/e^2':>10} {'S_pert/e':>9} {'<m1>pred':>9} {'S_Y/e':>9} {'<m1>_Y':>9} {'S_const/e':>9} {'<m1>_c':>9}")
for mu in (1.0, 0.3, 0.02):
    A = 1+mu; Dn = mu*mu + 6*mu + 1
    cl = (mu*A/Dn, -2*mu/Dn, mu*(mu*mu + 4*mu + 1)/Dn, -mu*A/Dn)
    c0 = C[(mu, 1.01)]; guess = (c0['U0'], c0['y0']); dprev = 0.01
    for d in (1e-2, 3e-3, 1e-3, 3e-4, 1e-4, 3e-5, 1e-5, 3e-6, 1e-6):
        k = (1+d)*kH(mu); lam = k/mu
        sc = np.sqrt(d/dprev)
        kprev = (1+dprev)*kH(mu)
        g = (guess[0]*sc, k + (guess[1] - kprev)*sc)
        try:
            cyc = find_cycle(mu, k, guess=g)
        except Exception as ex:
            print(f"{mu:5.2f} {d:7.0e} find_cycle failed: {ex}"); continue
        if cyc['U0'] < 1e-6 or cyc['residual'] > 1e-8:
            print(f"{mu:5.2f} {d:7.0e} collapsed/unconverged (U0={cyc['U0']:.2e}, res={cyc['residual']:.1e})"); continue
        guess = (cyc['U0'], cyc['y0']); dprev = d
        R = profile(mu, k, cyc)
        T = cyc['T']; tt = np.linspace(0, T, 8001)
        Y = R['sol'].sol(tt); U, P, y, eta, zeta = Y[:5]; p = lam - P; x = k*np.exp(U); pd = mu*P - y
        # first-harmonic amplitude of U
        ph = 2*np.pi*tt[:-1]/T
        eps = 2*abs(np.mean(U[:-1]*np.exp(-1j*ph)))
        # regression of the true pair
        M = np.vstack([np.ones_like(U), U, p]).T
        ce = np.linalg.lstsq(M, eta - mu, rcond=None)[0]; cz = np.linalg.lstsq(M, zeta - mu, rcond=None)[0]
        # scores
        S_true = score(mu, k, T, tt, U, p, y, eta, zeta, zeta + (mu-1)*eta - eta*eta, x - (1+eta)*zeta)
        def aff(c1, c2, d1, d2):
            e = mu + c1*U + c2*p; z = mu + d1*U + d2*p
            return score(mu, k, T, tt, U, p, y, e, z, c1*p + c2*pd, d1*p + d2*pd)
        S_jet = aff(*cl)
        pert = (1.3*cl[0], cl[1], cl[2], 0.8*cl[3])
        S_pert = aff(*pert); m1p = m1_mean(mu, *pert)
        # omega^Y direction (y/P, x/P): first-order jet on the centre plane is (mu U/A, 0) for eta~, and for zeta~: mu(1+U)(1+p/A) -> (mu, mu/A)
        eY = y/P; zY = x/P
        DeY = (x - y)/P + y*pd/P**2; DzY = x*p/P + x*pd/P**2
        S_Y = score(mu, k, T, tt, U, p, y, eY, zY, DeY, DzY)
        m1Y = m1_mean(mu, mu/A, 0.0, mu, mu/A)
        S_c = score(mu, k, T, tt, U, p, y, mu + 0*U, mu + 0*U, 0*U, 0*U); m1c = m1_mean(mu, 0, 0, 0, 0)
        print(f"{mu:5.2f} {d:7.0e} {eps:9.3e} {R['lam_eta']:+11.3e} | {ce[1]:9.5f} {cl[0]:9.5f} {ce[2]:9.5f} {cl[1]:9.5f} {cz[1]:9.5f} {cl[2]:9.5f} {cz[2]:9.5f} {cl[3]:9.5f} | "
              f"{S_true/eps**2:10.5f} {S_jet/eps**2:10.4f} {S_pert/eps:9.5f} {m1p:9.5f} {S_Y/eps:9.5f} {m1Y:9.5f} {S_c/eps:9.5f} {m1c:9.5f}", flush=True)
