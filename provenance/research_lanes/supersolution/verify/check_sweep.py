#!/usr/bin/env python3
"""Referee check 8: random sweep (40 parameter pairs, mu log-uniform in [0.02, 1], rho log-uniform in [1.002, 12]) plus
three very-near-Hopf cycles.  For each cycle: Theorem R (signs + 4 exponent formulas), claim 6 identities
(<S2/S1> = mu, Q_E identity, (R_K) residual), Proposition B (three exponent formulas, chi(gamma), beta > 1, kappa > (1-mu)beta)."""
import sys, time
import numpy as np
from scipy.integrate import solve_ivp
ROOT = '/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03'
sys.path.insert(0, ROOT + '/work')
from whlib import find_cycle, kH
from vprofile import profile, G1
rng = np.random.default_rng(31415)
pairs = [(float(np.exp(rng.uniform(np.log(0.02), 0.0))), float(np.exp(rng.uniform(np.log(1.002), np.log(12.0))))) for _ in range(40)]
pairs += [(0.02, 1.0001), (0.3, 1.0001), (1.0, 1.0001)]

def get_cycle(mu, rho):
    lad = [1.5]
    if rho < 1.5:
        r = 1.5
        while (r - 1) > 1.6*(rho - 1):
            r = 1 + (r - 1)/2.5; lad.append(r)
    else:
        r = 1.5
        while r*1.5 < rho:
            r *= 1.5; lad.append(r)
    lad.append(rho)
    guess = None; prev = None
    for r in lad:
        k = r*kH(mu)
        if guess is not None and r < prev:      # amplitude ~ sqrt(rho-1) when going down
            sc = np.sqrt((r - 1)/(prev - 1)); kp = prev*kH(mu)
            guess = (guess[0]*sc, k + (guess[1] - kp)*sc)
        c = find_cycle(mu, k, guess=guess)
        if c['U0'] < 1e-7 or c['residual'] > 1e-8:
            raise RuntimeError('cycle not found')
        guess = (c['U0'], c['y0']); prev = r
    return c

def backward(mu, k, cyc, R):
    lam = k/mu; A = 1+mu; T = cyc['T']; fwd = R['sol']; gE = R['ed']['gammaE']
    def f(t, s):
        d, n = s; U = fwd.sol(t)[0]
        return [d*(1+d) - n, n*(d + 1 - mu) - k*np.exp(U)]
    s = np.array([gE - 1, gE*(gE - 1)])
    for _ in range(5):
        s = solve_ivp(f, [T, 0], s, method='DOP853', rtol=1e-12, atol=1e-250).y[:, -1]
    r = solve_ivp(lambda t, s: f(t, s[:2]) + [1 + s[0]], [T, 0], [s[0], s[1], 0.0], method='DOP853', rtol=1e-12, atol=[1e-250, 1e-250, 1e-13])
    gamma = -r.y[2, -1]/T
    def f3(t, s):
        d, n, L = s[:3]
        U, P, y = fwd.sol(t)[:3]; x = k*np.exp(U); p = lam - P
        be = 1 + d; ka = n + (1-mu)*be; ah = np.exp(L + gamma*(t - T))
        return [d*(1+d) - n, n*(d + 1 - mu) - x, -be, ah*x, ah, ah*x*(lam - A - 2*p), ah*(ka*P + be*y + x), ah*(x - mu*A), ah*(ka + mu*be + mu)]
    r = solve_ivp(f3, [T, 0], [s[0], s[1], 0.0] + [0.0]*6, method='DOP853', rtol=1e-12, atol=[1e-250, 1e-250, 1e-13] + [1e-14]*6, dense_output=True)
    J = -r.y[3:, -1]/T
    tt = np.unique(np.concatenate([r.t] + [r.t[:-1] + j/4*np.diff(r.t) for j in range(1, 4)]))
    d, n = r.sol(tt)[:2]
    return dict(gamma=gamma, F1=J[2]/J[3], F2=J[4]/J[5], chi_ratio=(J[0]/J[1])/(gamma*(gamma-1)*(gamma-mu)), dmin=d.min(), nmin=n.min())

print(f"{'mu':>7} {'rho':>8} {'T':>8} {'Umin':>9} {'<eta>-mu':>13} {'devR(4)':>8} {'dev B(3)':>8} {'chi rat-1':>9} {'minS1/mx':>8} {'minS4/mx':>8} {'min sc/mx':>9} {'minK/mx':>8} {'minW/mx':>9} {'<S2/S1>-mu':>10} {'QE-id':>8} {'R_K res':>8} {'min(be-1)':>9} {'min nu':>9} {'s':>5}")
W = dict(devR=0, devB=0, chi=0, S2=0, QE=0, RK=0)
neg = []
t_all = time.time()
for mu, rho in pairs:
    t0 = time.time()
    try:
        cyc = get_cycle(mu, rho)
    except Exception as ex:
        print(f"{mu:7.4f} {rho:8.4f} cycle search failed: {ex}"); continue
    k = cyc['k']
    R = profile(mu, k, cyc)
    B = backward(mu, k, cyc, R)
    le = R['lam_eta']
    devR = max(abs(v/le - 1) for v in (R['lam_K'], R['lam_W'], *R['lam_c']))
    devB = max(abs(v/le - 1) for v in (B['gamma'] - 1 - mu, B['F1'], B['F2']))
    minsc = min(float(np.min(R[('sc', c)])/np.max(R[('sc', c)])) for c in R['cs'])
    # (R_K) residual: K' from the K-ODE vs the scalar Riccati right-hand side with explicit coefficients
    K = R['K']; lhs = k*G1(R['U']) - (1 + R['eta'])*K
    rhs = k*G1(R['U']) - (1 + R['S2']/R['S1'])*K + (R['p']/R['S1'])*K*K
    rk = float(np.max(np.abs(lhs - rhs))/np.max(np.abs(lhs)))
    W['devR'] = max(W['devR'], devR); W['devB'] = max(W['devB'], devB); W['chi'] = max(W['chi'], abs(B['chi_ratio']-1))
    W['S2'] = max(W['S2'], abs(R['mean_S2S1'] - mu)/mu); W['QE'] = max(W['QE'], abs(R['mean_QE']/R['alpha'] - 1)); W['RK'] = max(W['RK'], rk)
    mins = dict(S1=R['S1'].min(), S4=R['S4'].min(), sc=minsc, K=R['K'].min(), W=R['W'].min(), d=B['dmin'], n=B['nmin'])
    for nm, v in mins.items():
        if v <= 0: neg.append((mu, rho, nm, v))
    print(f"{mu:7.4f} {rho:8.4f} {cyc['T']:8.3f} {R['Umin']:9.3f} {le:+13.6e} {devR:8.1e} {devB:8.1e} {B['chi_ratio']-1:+9.1e} {R['S1'].min()/R['S1'].max():8.1e} {R['S4'].min()/R['S4'].max():8.1e} {minsc:9.1e} "
          f"{R['K'].min()/R['K'].max():8.1e} {R['W'].min()/R['W'].max():9.1e} {R['mean_S2S1']-mu:+10.1e} {R['mean_QE']/R['alpha']-1:+8.1e} {rk:8.1e} {B['dmin']:9.2e} {B['nmin']:9.2e} {time.time()-t0:5.1f}", flush=True)
print("\nworst over sweep:", {a: f"{b:.2e}" for a, b in W.items()})
print("non-positive minima (name, value):", neg if neg else "none")
print(f"total {time.time()-t_all:.0f}s")
