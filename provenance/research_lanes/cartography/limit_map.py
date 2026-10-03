#!/usr/bin/env python3
"""Formal singular limit mu -> 0 at fixed rho (NOT a proof; matched asymptotics checked against the full system).

Pulse (fast, O(1) time, parameter-free):   x' = p x,  p' = -y,  y' = x - y,
    heteroclinic from (0, p_in>0, 0) to (0, p_out<0, 0);  G(p_in) := -p_out.
Slow phase (time 1/mu, x = y = 0):  p' = mu (lambda - p),  with  int p dt = 0  (U returns to pulse level):
    b := -p_out/lambda,  S = (1+b)(1-e^{-S}),  p_in' = lambda * sigma(b),  sigma(b) = S - b.
Limit return map on the pre-pulse value of p:  Pi(p) = lambda * sigma(G(p)/lambda).
Fixed points for given lambda (= rho in the limit)  <=>  Lambda(p) = lambda, where Lambda(p) = G(p)/B(p) and
B(p) solves sigma(B)/B = p/G(p).   Uniqueness in the limit  <=>  Lambda strictly monotone.
Limit multiplier m_0 = sigma'(b*) G'(p*).
"""
import sys, csv
sys.path.insert(0, '/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03/work/lanes/cartography')
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
import cart as C

U_START = -46.0


def G_of(p_in, rtol=1e-12):
    x0 = np.exp(U_START)
    y0 = x0 / (1 + p_in)
    z0 = [U_START, p_in - x0 / ((1 + p_in) * p_in), y0]

    def f(t, z):
        U, p, y = z
        return [p, -y, np.exp(U) - y]

    def ev(t, z):
        return z[0] - U_START
    ev.direction = -1
    ev.terminal = True
    sol = solve_ivp(f, [0, 1e7], z0, method='DOP853', rtol=rtol, atol=1e-30, events=ev, first_step=1e-3)
    assert sol.status == 1, 'pulse did not terminate'
    U, p, y = sol.y[:, -1]
    return -(p - y), float(np.max(sol.y[0]))          # remaining int y dt = y (x negligible)


def S_of(b):
    # nonzero root of S = (1+b)(1-exp(-S))
    f = lambda S: S - (1 + b) * (1 - np.exp(-S))
    lo = 1e-12 if b > 0 else 0
    return brentq(f, max(b, 1e-14), 2 * (1 + b) + 1, xtol=1e-15, rtol=1e-15)


def sigma(b):
    if b < 1e-4:
        return b - (2.0 / 3.0) * b * b + (4.0 / 9.0) * b ** 3
    return S_of(b) - b


def dsigma(b):
    """sigma'(b) = b e^{-S} / (1 - (1+b) e^{-S})  (from implicit differentiation), returned with its log."""
    if b < 1e-4:
        return 1 - (4.0 / 3.0) * b
    S = S_of(b)
    return b * np.exp(-S) / (1 - (1 + b) * np.exp(-S))


def main():
    ps = np.concatenate((np.geomspace(0.02, 1.0, 26)[:-1], np.geomspace(1.0, 320.0, 60)))
    rows = []
    for p in ps:
        G, Umax = G_of(p)
        rows.append([p, G, Umax])
    rows = np.array(rows)
    p, G = rows[:, 0], rows[:, 1]
    print('pulse amplification G(p)/p: min', (G / p).min(), 'max', (G / p).max(), ' (small-p theory 1 + 2p/3:',
          f'G(0.02)/0.02 = {G[0]/p[0]:.6f} vs {1 + 2*0.02/3:.6f})')
    B = np.array([brentq(lambda b: sigma(b) / b - pi / Gi, 1e-9, 1e9, xtol=1e-14, rtol=1e-13) for pi, Gi in zip(p, G)])
    Lam = G / B
    dG = np.array([(G_of(pi * (1 + 1e-5))[0] - G_of(pi * (1 - 1e-5))[0]) / (2e-5 * pi) for pi in p])
    m0 = np.array([dsigma(b) for b in B]) * dG
    print('G(p)/p increasing:', bool(np.all(np.diff(G / p) > 0)), ' G\'(p) > 1 everywhere:', bool(np.all(dG > 1)),
          ' Lambda - 1 at smallest p:', Lam[0] - 1, ' ratio (Lambda-1)/p^2 at the 3 smallest p:', (Lam[:3] - 1) / p[:3] ** 2)
    print('log m_0 / (Lambda-1)^1.5 at the 5 smallest p:', np.log(m0[:5]) / (Lam[:5] - 1) ** 1.5)
    mono = bool(np.all(np.diff(Lam) > 0))
    print('Lambda(p) strictly increasing on the sampled range:', mono, '  Lambda range', Lam.min(), Lam.max())
    print('limit multiplier m_0 in (0,1) on the range:', bool(np.all((m0 > 0) & (m0 < 1))), ' max m_0 =', m0.max(), 'at Lambda =', Lam[np.argmax(m0)])
    with open('limit_map.csv', 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['p_in', 'G', 'pulse_log_x_max', 'b_star', 'Lambda=rho_limit', 'm0_limit_multiplier', 'log_m0', 'dG_dp'])
        for i in range(len(p)):
            w.writerow([repr(float(v)) for v in (p[i], G[i], rows[i, 2], B[i], Lam[i], m0[i], np.log(m0[i]), dG[i])])
    # ---- comparison with the full system at small mu
    from hunt import primary_cycle
    print('\ncomparison with the full system (prediction at mu=0 by interpolation in Lambda):')
    print('  rho     mu        log m_t(full)   log m_0(limit)   p_max/lambda(full)  p*/lambda(limit)   -p_min/lambda(full)  b*(limit)   mu*T(full)  S*(limit)')
    for rho in (1.03, 1.1, 1.5, 3.0, 10.0, 50.0, 200.0):
        lm0 = np.interp(rho, Lam, np.log(m0))
        pst = np.interp(rho, Lam, p)
        bst = np.interp(rho, Lam, B)
        for mu in (0.01, 0.002, 0.0005):
            try:
                cyc, _, _, _ = primary_cycle(mu, rho)
            except Exception as e:
                print('  ', rho, mu, 'failed', e)
                continue
            k = rho * C.kH(mu)
            pr = C.profile(mu, k, cyc['U0'], cyc['y0'], cyc['T'], nmin=4000)
            lam = k / mu
            print(f"  {rho:<6} {mu:<8} {pr['log_mt']:+.5f}        {lm0:+.5f}        {pr['p'].max()/lam:.5f}             {pst/rho:.5f}            {-pr['p'].min()/lam:.5f}             {bst:.5f}     {mu*cyc['T']:.4f}     {S_of(bst):.4f}")


if __name__ == '__main__':
    main()
