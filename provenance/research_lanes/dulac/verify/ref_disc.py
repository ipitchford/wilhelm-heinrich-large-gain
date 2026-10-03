#!/usr/bin/env python3
"""Light sanity check (double precision, whlib.unstable_manifold_orbit) of the CONDITIONAL disc statements in
NOTES.md (table item 4 and O1 on D\\{E}): V>0, Q>0, tangent orthogonality O=0 and eta = eta0 - p*kappa along
spirals on W^u(E).  Also B1(ii): eta - mu -> 2 sigma at E.  This is not part of the 'proved on cycles' claims and
proves nothing about brief item 3(b); it only checks that nothing visibly contradicts it."""
import sys
import numpy as np
sys.path.insert(0, '/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03/work')
import whlib

for mu, rho in ((0.02, 1.5), (0.3, 1.01), (0.3, 4.0), (1.0, 1.5), (1.0, 4.0)):
    k = rho * mu * (1 + mu)
    worst = dict(V=np.inf, Q=np.inf, O=0.0, id=0.0)
    n_tot = 0
    for theta in (0.0, 1.0, 2.0, 3.0, 4.0, 5.0):
        o = whlib.unstable_manifold_orbit(mu, k, r0=1e-6, theta=theta, n=40001)
        U, p, y, x, eta, zeta = (o[s] for s in ('U', 'p', 'y', 'x', 'eta', 'zeta'))
        r = np.hypot(U, p)
        sel = (r > 1e-4) & np.isfinite(eta) & (U > -30)        # keep away from E (cancellation) and from deep troughs
        U, p, y, x, eta, zeta = (a[sel] for a in (U, p, y, x, eta, zeta))
        pd = k - mu * p - y; Hd = x - y; H = y - k
        V = p * p - U * pd; Q = -pd * H + p * Hd; V3 = U * Hd - p * H
        K = zeta * U - eta * p - H
        O = zeta * p - eta * pd - Hd
        worst['V'] = min(worst['V'], np.min(V / (U * U + p * p)))
        worst['Q'] = min(worst['Q'], np.min(Q / (p * p + pd * pd)))
        worst['O'] = max(worst['O'], np.max(np.abs(O) / (np.abs(zeta * p) + np.abs(eta * pd) + np.abs(Hd))))
        worst['id'] = max(worst['id'], np.max(np.abs(eta - V3 / V + p * K / V) / (np.abs(eta) + np.abs(V3 / V) + np.abs(p * K / V))))
        n_tot += len(U)
    ed = whlib.equilibrium_data(mu, k)
    o = whlib.unstable_manifold_orbit(mu, k, r0=1e-6, theta=0.0, n=40001)
    r = np.hypot(o['U'], o['p'])
    near = (r < 1e-4)
    print(f'mu={mu:<5} rho={rho:<5} samples={n_tot}  min V/r^2={worst["V"]:.3e}  min Q/(p^2+pd^2)={worst["Q"]:.3e}  max|O|rel={worst["O"]:.1e}  '
          f'max|eta-eta0+p kappa|rel={worst["id"]:.1e}  | near E: eta-mu in [{(o["eta"][near]-mu).min():.6f},{(o["eta"][near]-mu).max():.6f}]  2sigma={2*ed["lam_u"].real:.6f}  gammaE-A={ed["gammaE"]-1-mu:.6f}')
