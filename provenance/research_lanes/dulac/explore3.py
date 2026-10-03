#!/usr/bin/env python3
"""Sign patterns of the eta-free Dulac residuals R_i on the limit cycle (last revolution of the spiral),
by tangent quadrant:  QI p>0,p'>0 ; QII p>0,p'<0 ; QIII p<0,p'<0 ; QIV p<0,p'>0."""
import numpy as np
from lpkit import *

def last_rev(u):
    """indices of the final revolution of orbit 0 (essentially the limit cycle)."""
    m = u['orb'] == 0
    idx = np.flatnonzero(m)
    p = u['p'][idx]
    # crest crossings: p goes + to -
    cr = np.flatnonzero((p[:-1] > 0) & (p[1:] <= 0))
    return idx[cr[-2]:cr[-1] + 1]

def main():
  for mu, rho in SETS:
      u = load(mu, rho)
      ii = last_rev(u)
      U, p, pd, H, x, y, eta, zeta, K, W, L = (u[s][ii] for s in ('U', 'p', 'pd', 'H', 'x', 'y', 'eta', 'zeta', 'K', 'W', 'L'))
      t = u['t'][ii]; T = t[-1] - t[0]
      A = u['A']; k = u['k']
      Vrad = p * p - U * pd
      Qtan = pd * (-H) + p * (u['x'][ii] - u['y'][ii])
      R = {'rad': -p * K / Vrad, 'tan': -p * W / Qtan, 'zeta': x / zeta - A, 'eta': zeta / eta - 1,
           'K': k * u['psi'][ii] / K - A, 'W': x * p * p / W - A, 'q=eta-mu': eta - mu}
      quad = {'QI': (p > 0) & (pd > 0), 'QII': (p > 0) & (pd <= 0), 'QIII': (p <= 0) & (pd <= 0), 'QIV': (p <= 0) & (pd > 0)}
      lam = np.trapezoid(eta - mu, t) / T
      print(f'mu={mu} rho={rho}: T={T:.3f} lambda_c={lam:+.5f}  means: ' + ' '.join(f'{n}={np.trapezoid(v, t)/T:+.5f}' for n, v in R.items()))
      for qn, m in quad.items():
          print(f'   {qn:<5} frac={m.mean():.2f} ' + ' '.join(f'{n}:[{v[m].min():+.3g},{v[m].max():+.3g}]' for n, v in R.items()))


if __name__ == '__main__':
    main()
