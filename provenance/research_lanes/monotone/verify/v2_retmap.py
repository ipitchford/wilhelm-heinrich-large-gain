#!/usr/bin/env python3
"""Referee test of claim 5 (lane R7): crest return map P on R = {p=0, x>=y>k}.
 (a) existence of the first return from every sampled point of R (the lane's proof does not address this);
 (b) P(R) in int R:  x1 > y1 > k;
 (c) DP entrywise > 0, by TWO methods: (ii) full variational equation in (U,p,y) + projection along the flow
     [independent of the lane's argument; resolvable only above a roundoff floor], and
     (iii) compound system w' = C w in log variables [the lane's/v2.0 mechanism; resolves entries of any size];
 (d) at the cycle crest: Perron eigenvector slope of DP  vs  zeta(crest) from the forward Riccati equation,
     eigenvalues vs exp(lambda_c T), det vs exp(-A T).
"""
import sys, warnings, time
import numpy as np
from scipy.integrate import solve_ivp
sys.path.insert(0, '/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03/work/lanes/monotone/verify')
from vlib import *
warnings.filterwarnings('ignore'); np.seterr(all='ignore')


def var_return(mu, k, U0, y0, tmax):
    """variational equation in standard coordinates; returns (U1,y1,tau,M2x2,|Phi|max)"""
    def f(t, s):
        U, p, y = s[:3]
        x = k * np.exp(U)
        Phi = s[3:].reshape(3, 3)
        J = np.array([[0, 1, 0], [0, -mu, -1], [x, 0, -1.0]])
        return np.concatenate(([p, k - mu * p - y, x - y], (J @ Phi).ravel()))

    def ev(t, s):
        return s[1]
    ev.direction = -1; ev.terminal = True
    s0 = np.concatenate(([U0, 0.0, y0], np.eye(3).ravel()))
    a = solve_ivp(f, [0, 1e-4], s0, method='DOP853', rtol=1e-12, atol=1e-14)
    b = solve_ivp(f, [1e-4, tmax], a.y[:, -1], method='DOP853', rtol=1e-12, atol=1e-14, events=ev)
    if b.status != 1:
        raise RuntimeError('no return (variational)')
    s = b.y_events[0][0]; tau = b.t_events[0][0]
    U1, p1, y1 = s[:3]; Phi = s[3:].reshape(3, 3)
    fz = np.array([p1, k - mu * p1 - y1, k * np.exp(U1) - y1])
    M = np.empty((2, 2))
    for j, v in enumerate((np.array([1.0, 0, 0]), np.array([0, 0, 1.0]))):
        d = Phi @ v
        d = d - fz * (d[1] / fz[1])
        M[:, j] = [d[0], d[2]]
    return U1, y1, tau, M, np.abs(Phi).max() * (1 + np.abs(fz).max() / abs(fz[1]))


def compound_return(mu, k, U0, y0):
    """log-variable compound system along the orbit; returns U1, y1, tau, ln M (2x2), all entries >0 by construction"""
    A = 1 + mu
    U1, y1, tau, sol = ret(mu, k, U0, y0, dense=True)
    x0 = k * np.exp(U0)
    C0 = np.array([[-mu, 1, 0], [0, -1, 1], [x0, 0, -A]])
    dlt = 1e-7 / max(1.0, x0, A)
    E = np.eye(3) + dlt * C0 + dlt ** 2 / 2 * C0 @ C0 + dlt ** 3 / 6 * C0 @ C0 @ C0
    lnM = np.empty((2, 2))

    def g(t, l):
        U = sol(np.array([t]))[0, 0]
        return [-mu + np.exp(l[1] - l[0]), -1.0 + np.exp(l[2] - l[1]), -A + k * np.exp(U + l[0] - l[2])]
    for j, w0 in enumerate((np.array([y0 - k, x0 - y0, 0.0]), np.array([0.0, 0.0, y0 - k]))):
        w = E @ w0
        if not (w > 0).all():
            raise RuntimeError(f'Taylor start not positive: {w}')
        sc = solve_ivp(g, [dlt, tau], np.log(w), method='DOP853', rtol=1e-11, atol=1e-12)
        l = sc.y[:, -1]
        lnM[:, j] = [l[0] - np.log(y1 - k), l[2] - np.log(y1 - k)]
        if j == 0:
            x1_minus_y1_over_dU = np.exp(l[1] - l[0]) * (y1 - k)     # w2/w1 = (x1-y1)/(y1-k)
    return U1, y1, tau, lnM, x1_minus_y1_over_dU


def perron(lnM):
    s = lnM.max()
    a, b, c, d = np.exp(lnM - s).ravel()   # [[a,b],[c,d]]
    disc = np.sqrt((a - d) ** 2 + 4 * b * c)
    lam1 = 0.5 * (a + d + disc)
    slope = (lam1 - a) / b if d >= a else 2 * c / ((a - d) + disc)
    lndet = np.log(a * d - b * c) + 2 * s if a * d > b * c else np.nan
    return np.log(lam1) + s, slope, lndet


if __name__ == '__main__':
    grid = [(mu, rho) for mu in (1.0, 0.3, 0.02) for rho in (1.01, 1.5, 4.0, 12.0)] + [(3.0, 1.5), (3.0, 4.0)]
    if len(sys.argv) > 1:
        a = sys.argv[1:]; grid = [(float(a[i]), float(a[i + 1])) for i in range(0, len(a), 2)]
    tot = dict(pts=0, noret=0, notR=0, var_neg=0, var_res=0, var_unres=0, mism=0)
    for mu, rho in grid:
        t0 = time.time()
        k = rho * kH(mu); A = 1 + mu; lam = k / mu
        cyc = find_cycle(mu, k); T = cyc['T']
        # forward Riccati for zeta(crest), lambda_c
        ed = eq_data(mu, k)

        def gf(t, w):
            U = cyc['sol'](np.mod(t, T))[0, 0]
            return [np.exp(w[1]) + (mu - 1) * w[0] - w[0] ** 2, k * np.exp(U - w[1]) - (1 + w[0]), w[0]]
        w0 = np.array([ed['etaE'], np.log(ed['zetaE']), 0.0])
        for _ in range(4):
            sf = solve_ivp(gf, [0, T], w0, method='DOP853', rtol=RT, atol=AT)
            lamc = sf.y[2, -1] / T - mu; w0 = np.array([sf.y[0, -1], sf.y[1, -1], 0.0])
        zeta_c = np.exp(w0[1])
        # (d) cycle crest
        U1, y1, tau, lnM, _ = compound_return(mu, k, cyc['U0'], cyc['y0'])
        lnl, slope, lndet = perron(lnM)
        line = (f"== mu={mu} rho={rho} T={T:.4f}: crest map by compound system: ln(m_t)/T={lnl / T:+.9f} vs lambda_c={lamc:+.9f};"
                f" Perron slope={slope:.9g} vs zeta(crest)={zeta_c:.9g} (rel diff {abs(slope / zeta_c - 1):.1e});"
                f" ln det/T={lndet / T:+.6f} vs -A={-A:+.6f}")
        print(line, flush=True)
        # (a)-(c) sample of R
        fs = (0.02, 0.2, 0.6, 1.0, 1.4, 2.0) if rho < 12 else (0.2, 0.6, 1.0, 1.1)
        ths = (1e-4, 1e-2, 0.2, 0.6, 0.95, 1.0)
        npts = 0; notR = 0; var_neg = 0; var_res = 0; var_unres = 0; mism = 0; noret = 0
        min_margin = np.inf; lg_min = np.inf; lg_max = -np.inf; worst_mis = 0.0; taumax = 0.0
        for f_ in fs:
            for th in ths:
                U0 = f_ * cyc['U0']; x0 = k * np.exp(U0); y0 = k + th * (x0 - k)
                npts += 1
                try:
                    U1, y1, tau, lnM, x1my1 = compound_return(mu, k, U0, y0)
                except Exception as e:
                    noret += 1; print('   NO RETURN / failure at', U0, y0, repr(e)); continue
                taumax = max(taumax, tau)
                x1 = k * np.exp(U1)
                if not (x1 > y1 > k):
                    notR += 1; print('   NOT IN int R:', U0, y0, '->', U1, y1, x1)
                min_margin = min(min_margin, (x1 - y1) / x1, (y1 - k) / k)
                lg_min = min(lg_min, lnM.min() / np.log(10)); lg_max = max(lg_max, lnM.max() / np.log(10))
                try:
                    Uv, yv, tv, Mv, floor = var_return(mu, k, U0, y0, tmax=tau * 1.5 + 10)
                except Exception as e:
                    print('   variational failure at', U0, y0, repr(e)); continue
                floor = 1e-9 * floor
                Mc = np.exp(lnM)
                for i in range(2):
                    for j in range(2):
                        if abs(Mv[i, j]) > 30 * floor or Mc[i, j] > 30 * floor:
                            var_res += 1
                            if Mv[i, j] <= 0:
                                var_neg += 1; print('   NEGATIVE resolved entry (variational):', U0, y0, i, j, Mv[i, j], Mc[i, j], floor)
                            rel = abs(Mv[i, j] / Mc[i, j] - 1)
                            worst_mis = max(worst_mis, rel if Mc[i, j] > 1e4 * floor else 0.0)
                            if Mc[i, j] > 1e4 * floor and rel > 1e-3:
                                mism += 1
                        else:
                            var_unres += 1
        print(f"   R sample: {npts} points (U0/U_crest in {fs}, (y-k)/(x-k) in {ths}); no-return={noret}; outside int R={notR};"
              f" min margin of (x1-y1)/x1,(y1-k)/k = {min_margin:.2e}; max return time={taumax:.1f};"
              f" log10 DP entries in [{lg_min:.1f},{lg_max:.1f}]; variational method: resolved entries={var_res} (negative={var_neg}, "
              f"mismatch vs compound>1e-3: {mism}, worst rel diff {worst_mis:.1e}), unresolved={var_unres}  [{time.time() - t0:.0f}s]", flush=True)
        for kk, vv in dict(pts=npts, noret=noret, notR=notR, var_neg=var_neg, var_res=var_res, var_unres=var_unres, mism=mism).items():
            tot[kk] += vv
    print('TOTAL', tot)
