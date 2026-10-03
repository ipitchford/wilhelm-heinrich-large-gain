#!/usr/bin/env python3
"""Referee cross-check of Lemma 1 / Lemma 2 against the FULL monodromy matrix (no Riccati in the monodromy):
at several base phases t0, Phi = fundamental matrix of v' = J v over [t0, t0+T] in (U,p,y) coordinates.
 (i)  Phi e_s = m_s e_s          with e_s=(1,-r,-s) from the backward Riccati system, m_s = exp(-(A+lambda_c)T)
 (ii) Phi^T q = m_t q            with q=(r y'-s p', y'+s p, -(p'+r p)) (the claimed weak adjoint direction), m_t = exp(lambda_c T)
 (iii) eigenvector of Phi for the smallest multiplier, normalised to first component 1, equals (1,-r,-s)
lambda_c is taken from the forward (eta,zeta) Riccati mean, i.e. from a third independent computation."""
import sys, warnings
import numpy as np
from scipy.integrate import solve_ivp
sys.path.insert(0, '/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03/work/lanes/monotone/verify')
from vlib import *
warnings.filterwarnings('ignore'); np.seterr(all='ignore')
grid = [(mu, rho) for mu in (1.0, 0.3, 0.02) for rho in (1.01, 1.5, 4.0)] + [(3.0, 1.01), (3.0, 1.5), (3.0, 4.0), (0.1, 1.2), (7.0, 2.5)]
for mu, rho in grid:
    k = rho * kH(mu); A = 1 + mu
    cyc = find_cycle(mu, k); T = cyc['T']; ed = eq_data(mu, k)

    def gb(tau, w):
        U = cyc['sol'](np.mod(T - tau, T))[0, 0]
        return [-(w[0] ** 2 - mu * w[0] - np.exp(w[1])), -((w[0] - 1.0) - k * np.exp(U - w[1]))]
    w0 = np.array([ed['gE'], np.log(ed['gE'] * (ed['gE'] - mu))])
    for _ in range(3):
        sb = solve_ivp(gb, [0, T], w0, method='DOP853', rtol=RT, atol=AT, dense_output=True); w0 = sb.y[:, -1].copy()

    def gf(t, w):
        U = cyc['sol'](np.mod(t, T))[0, 0]
        return [np.exp(w[1]) + (mu - 1) * w[0] - w[0] ** 2, k * np.exp(U - w[1]) - (1 + w[0]), w[0]]
    w0 = np.array([ed['etaE'], np.log(ed['zetaE']), 0.0])
    for _ in range(4):
        sf = solve_ivp(gf, [0, T], w0, method='DOP853', rtol=RT, atol=AT)
        lamc = sf.y[2, -1] / T - mu; w0 = np.array([sf.y[0, -1], sf.y[1, -1], 0.0])
    m_t = np.exp(lamc * T); m_s = np.exp(-(A + lamc) * T)

    def fv(t, s):
        U, p, y = s[:3]; x = k * np.exp(U)
        Phi = s[3:].reshape(3, 3)
        J = np.array([[0, 1, 0], [0, -mu, -1], [x, 0, -1.0]])
        return np.concatenate(([p, k - mu * p - y, x - y], (J @ Phi).ravel()))
    wl = wr = we = 0.0; nres_r = 0
    for frac in (0.0, 0.13, 0.31, 0.5, 0.68, 0.87):
        t0 = frac * T
        st = state(mu, k, cyc, np.array([t0])); st = {n: float(v[0]) for n, v in st.items()}
        r, ls = sb.sol(T - t0); s = np.exp(ls)
        sol = solve_ivp(fv, [0, T], np.concatenate(([st['U'], st['p'], st['y']], np.eye(3).ravel())), method='DOP853', rtol=1e-12, atol=1e-14)
        Phi = sol.y[3:, -1].reshape(3, 3)
        q = np.array([r * st['yp'] - s * st['pp'], st['yp'] + s * st['p'], -(st['pp'] + r * st['p'])])
        es = np.array([1.0, -r, -s])
        wl = max(wl, np.linalg.norm(Phi.T @ q - m_t * q) / (m_t * np.linalg.norm(q)))
        if m_s > 1e-7:
            nres_r += 1
            wr = max(wr, np.linalg.norm(Phi @ es - m_s * es) / (m_s * np.linalg.norm(es)))
            w, V = np.linalg.eig(Phi)
            i = int(np.argmin(np.abs(w))); v = V[:, i].real; v = v / v[0]
            we = max(we, np.linalg.norm(v - es) / np.linalg.norm(es))
    print(f"mu={mu} rho={rho}: m_t={m_t:.3e} m_s={m_s:.3e} | (ii) max rel residual of Phi^T q = m_t q: {wl:.1e}"
          + (f" | (i) Phi e_s = m_s e_s: {wr:.1e} | (iii) eig-vector vs (1,-r,-s): {we:.1e}" if nres_r else " | (i),(iii): m_s below roundoff of Phi, not resolvable"), flush=True)
