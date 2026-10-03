#!/usr/bin/env python3
"""Tolerance-tightening checks of the decisive numbers (rt = integrator relative tolerance)."""
import sys
sys.path.insert(0, '/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03/work/lanes/cartography')
import numpy as np
import cart as C

print('(i) weakest exponents on the grid: lambda_c*T at three tolerances (Newton re-converged at each)')
for mu, rho, U0, y0k in [(0.002, 1.001, 0.8646, None), (0.005, 1.001, 0.5804, None), (1.0, 1.001, 0.1232, None), (0.002, 1.0127, 2.336, None)]:
    k = rho * C.kH(mu)
    out = []
    for rt in (1e-10, 1e-12, 1e-13):
        y0g, kg = C.hopf_guess(mu, U0)
        c = C.newton_fixed_h(mu, U0, y0g, kg, rt=rt, maxit=60)
        # secant on h to hit rho
        h0, r0 = U0, c['rho']
        h1 = U0 * 1.01
        c1 = C.newton_fixed_h(mu, h1, c['y0'], c['k'], rt=rt)
        for _ in range(30):
            h2 = h1 - (c1['rho'] - rho) * (h1 - h0) / (c1['rho'] - r0)
            h0, r0 = h1, c1['rho']
            h1 = h2
            c1 = C.newton_fixed_h(mu, h1, c1['y0'], c1['k'], rt=rt)
            if abs(c1['rho'] - rho) < 1e-13 * rho:
                break
        pr = C.profile(mu, k, c1['U0'], c1['y0'], c1['T'], record=False, rt=rt)
        out.append((rt, c1['U0'], pr['log_mt'], np.log(c1['m_t'])))
    print(f"  mu={mu} rho-1={rho-1:.4g}: " + ' | '.join(f"rt={rt:g}: U*={u:.9f} logmt(Riccati)={a:+.10e} logmt(M)={b:+.10e}" for rt, u, a, b in out))

print('(ii) non-monotone r(U) = log(phi(U)/U) on W^u(E) at mu=0.0027, rho=5.02 (three tolerances)')
mu, rho = 0.0027, 5.02
k = rho * C.kH(mu)
for rt in (1e-10, 1e-12, 1e-13):
    par = np.array([mu, k])
    rtol, atol = C.tolerances(0, mu, k, rt)
    vals = []
    for th in np.linspace(0, 2 * np.pi, 9)[:-1]:
        z, ed = C.unstable_seed(mu, k, 1e-5, th)
        ev_t = np.zeros(80)
        ev_y = np.zeros((80, 3))
        s = z.copy()
        st, t, ns, nev, nrec = C.integrate(0, par, s, 0.0, 1e9, rtol, atol, 1, -2, 80, ev_t, ev_y, C._EMPTY_T, C._EMPTY_Y, np.inf, 50_000_000)
        Uc = ev_y[:nev, 0]
        vals.append(np.column_stack((Uc[2:-1], np.log(Uc[3:] / Uc[2:-1]))))
    V = np.vstack(vals)
    V = V[np.argsort(V[:, 0])]
    Ustar = V[-1, 0]
    sel = (V[:, 0] > 0.05 * Ustar) & (V[:, 0] < 0.95 * Ustar)
    i_min = np.argmin(np.where(V[:, 0] < 0.4 * Ustar, V[:, 1], np.inf))
    j = np.argmax(np.where((V[:, 0] > V[i_min, 0]), V[:, 1], -np.inf))
    print(f"  rt={rt:g}: r(0+)={2*np.pi*ed['alpha']/ed['omega']:.6f}; local min r={V[i_min,1]:.8f} at U/U*={V[i_min,0]/Ustar:.4f}; later local max r={V[j,1]:.8f} at U/U*={V[j,0]/Ustar:.4f}; rise={V[j,1]-V[i_min,1]:+.6f}")

print('(iii) S3a at the cycle crest for mu=1, rho=5: d/dU[(mu-zeta)/(k(e^U-1))] by finite differences along W^u(E) crests')
mu, rho = 1.0, 5.0
k = rho * C.kH(mu)
for rt in (1e-10, 1e-12, 1e-13):
    par = np.array([mu, k])
    rtol, atol = C.tolerances(1, mu, k, rt)
    pts = []
    for th in np.linspace(0, 2 * np.pi, 25)[:-1]:
        z, ed = C.unstable_seed(mu, k, 1e-6, th)
        s = np.zeros(12)
        s[:3] = z
        s[3], s[4] = ed['etaE'], ed['zetaE']
        s[5], s[6] = 1e-12 * k, 1e-12 * k
        ev_t = np.zeros(40)
        ev_y = np.zeros((40, 12))
        st, t, ns, nev, nrec = C.integrate(1, par, s, 0.0, 1e6, rtol, atol, 1, -2, 40, ev_t, ev_y, C._EMPTY_T, C._EMPTY_Y, np.inf, 50_000_000)
        pts.append(ev_y[4:nev, :5])
    Pp = np.vstack(pts)
    Pp = Pp[np.argsort(Pp[:, 0])]
    U, zeta = Pp[:, 0], Pp[:, 4]
    Ustar = U[-1]
    keep = np.concatenate(([True], np.diff(U) > 1e-9))
    U, zeta = U[keep], zeta[keep]
    R = (mu - zeta) / (k * np.expm1(U))
    near = (U > Ustar * (1 - 3e-6)) & (U < Ustar * (1 - 1e-8))
    sl = np.polyfit(U[near], R[near], 1)[0] if near.sum() >= 3 else np.nan
    mid = (U > 0.5 * Ustar) & (U < 0.9 * Ustar)
    slm = np.polyfit(U[mid], R[mid], 1)[0]
    print(f"  rt={rt:g}: U*={Ustar:.9f}; slope of R just inside the cycle (n={int(near.sum())}) = {sl:+.5f}; slope at mid amplitudes = {slm:+.5f}")
