#!/usr/bin/env python3
"""Quick diagnostics on sampled attractor discs: signs of simple candidates."""
import numpy as np
from gen_data import SETS, fname, unpack

np.set_printoptions(linewidth=200, precision=4, suppress=False)


def rng(a):
    return f'[{a.min():+.4g},{a.max():+.4g}]'


for mu, rho in SETS:
    u = unpack(np.load(fname(mu, rho)))
    U, p, H, x, y, eta, zeta, K, W, L, pd = (u[s] for s in ('U', 'p', 'H', 'x', 'y', 'eta', 'zeta', 'K', 'W', 'L', 'pd'))
    A = u['A']; k = u['k']; gE = u['gE']
    Vrad = p * p - U * pd
    Qtan = pd * (-H) + p * (x - y)
    disc = (eta - mu)**2 - 4 * zeta
    print(f'mu={mu} rho={rho}: sigma={u["sig"]:.4g} omega={u["om"]:.4g} gE={gE:.4g}')
    print('   eta-mu', rng(eta - mu), ' min Vrad/r2', f'{np.min(Vrad/(U*U+p*p)):.3g}', ' min Qtan/(p^2+pd^2)',
          f'{np.min(Qtan/(p*p+pd*pd)):.3g}', ' max disc', f'{disc.max():.3g}')
    print('   xp2/W - A', rng(x * p * p / W - A), ' kpsi/K - A', rng(k * u['psi'] / K - A),
          ' x/zeta - A', rng(x / zeta - A), ' kphi/L - gE', rng(k * u['phi'] / L - gE))
    print('   K/L', rng(K / L), ' W/K', rng(W / K), ' K/Vrad', rng(K / Vrad), ' W/Qtan', rng(W / Qtan))
