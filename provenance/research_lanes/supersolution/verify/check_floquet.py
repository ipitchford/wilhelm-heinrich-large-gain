#!/usr/bin/env python3
"""Referee check 7: the Floquet objects behind Proposition B and Lemma A0, from the plain 3x3 variational monodromy
(no Riccati equations).  Multipliers {1, m_t, m_s}; strong-stable eigenvector must have sign pattern (+,-,-)
(i.e. (a,b,c) > 0), beta(0) = -v_p/v_U must agree with the backward Riccati value, and the compound monodromy
(second multiplicative compound, sign-conjugated) must be a nonnegative matrix with spectral radius m_t."""
import numpy as np
from scipy.integrate import solve_ivp
from vprofile import profile, load_cycles
C = load_cycles()
print(f"{'mu':>5} {'rho':>5} {'mult 1':>12} {'log(m_t)/T':>13} {'<eta>-mu':>13} {'-log(m_s)/T':>12} {'1+<eta>':>10} {'v_s signs':>10} {'beta0 eig':>10} {'beta0 Ricc':>10} {'kap0 eig':>10} {'kap0 Ricc':>10} {'min C-monodromy entry':>21} {'spr/m_t':>9}")
for key in [(1.0, 1.01), (1.0, 1.5), (0.3, 1.01), (0.3, 1.5), (1.0, 4.0)]:
    mu, rho = key; cyc = C[key]; k = cyc['k']; lam = k/mu; T = cyc['T']; A = 1+mu
    R = profile(mu, k, cyc)
    fwd = R['sol']
    def f(t, s):
        U = fwd.sol(t)[0]; x = k*np.exp(U)
        J = np.array([[0, 1, 0], [0, -mu, -1], [x, 0, -1.0]])
        Cm = np.array([[-mu, 1, 0], [0, -1, 1], [x, 0, -A]])
        return np.concatenate(((J @ s[:9].reshape(3, 3)).ravel(), (Cm @ s[9:].reshape(3, 3)).ravel()))
    r = solve_ivp(f, [0, T], np.concatenate((np.eye(3).ravel(), np.eye(3).ravel())), method='DOP853', rtol=1e-13, atol=1e-18)
    M = r.y[:9, -1].reshape(3, 3); MC = r.y[9:, -1].reshape(3, 3)
    w, V = np.linalg.eig(M)
    o = np.argsort(-np.abs(w)); w = w[o].real; V = V[:, o].real
    vs = V[:, 2]/V[0, 2]
    # backward Riccati value at the crest (t = 0 = T)
    gE = R['ed']['gammaE']
    def g(t, s):
        d, n = s; U = fwd.sol(t)[0]
        return [d*(1+d) - n, n*(d + 1 - mu) - k*np.exp(U)]
    s = np.array([gE - 1, gE*(gE - 1)])
    for _ in range(5):
        s = solve_ivp(g, [T, 0], s, method='DOP853', rtol=1e-12, atol=1e-250).y[:, -1]
    beta0 = 1 + s[0]; kap0 = s[1] + (1-mu)*beta0
    wc = np.linalg.eigvals(MC)
    print(f"{mu:5.2f} {rho:5.2f} {w[0]:12.9f} {np.log(w[1])/T:+13.6e} {R['lam_eta']:+13.6e} {-np.log(w[2])/T:12.7f} {R['lam_eta']+A:10.7f} {''.join('+' if v > 0 else '-' for v in vs):>10} {-vs[1]:10.6f} {beta0:10.6f} {-vs[2]:10.6f} {kap0:10.6f} {MC.min():21.3e} {np.max(np.abs(wc))/w[1]:9.6f}")
