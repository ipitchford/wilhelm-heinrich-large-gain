#!/usr/bin/env python3
"""On the limit cycle: compare kappa = K/V_rad and varpi = W/Q_tan on the upper (p>0) and lower (p<0) branch at equal U.
lambda_c*T = -oint kappa dU = -int_{Umin}^{Umax} (kappa_up - kappa_low) dU."""
import numpy as np
from lpkit import *
from explore3 import last_rev

for mu, rho in [(1.0, 1.02), (1.0, 1.5), (1.0, 5.0), (0.2, 1.05), (0.05, 1.5)]:
    u = load(mu, rho)
    ii = last_rev(u)
    U, p, pd, H, eta, zeta, K, W = (u[s][ii] for s in ('U', 'p', 'pd', 'H', 'eta', 'zeta', 'K', 'W'))
    Vrad = p * p - U * pd
    Qtan = pd * (-H) + p * (u['x'][ii] - u['y'][ii])
    kap = K / Vrad; var = W / Qtan
    up = p > 0; lo = p < 0
    Ug = np.linspace(U.min(), U.max(), 13)[1:-1]
    def interp(Uq, m, f):
        o = np.argsort(U[m]); return np.interp(Uq, U[m][o], f[m][o])
    print(f'mu={mu} rho={rho}  Umin={U.min():.3f} Umax={U.max():.3f}')
    print('    U        kappa_up  kappa_low   varpi_up  varpi_low   eta_up   eta_low')
    for Uq in Ug:
        print(f'  {Uq:+8.3f}  {interp(Uq, up, kap):9.5f} {interp(Uq, lo, kap):9.5f}   {interp(Uq, up, var):9.5f} {interp(Uq, lo, var):9.5f}  {interp(Uq, up, eta):8.4f} {interp(Uq, lo, eta):8.4f}')
