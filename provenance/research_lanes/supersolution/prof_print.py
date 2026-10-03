import sys, numpy as np
from base import load
D = load()
mu = float(sys.argv[1]); rho = float(sys.argv[2]); m = int(sys.argv[3]) if len(sys.argv) > 3 else 24
P = D[(mu, rho)]; k = P['k']; lam = k/mu; A = 1+mu; n = len(P['t'])
idx = np.arange(0, n, n//m)
et = P['y']/(lam-P['p']); zt = P['x']/(lam-P['p'])
a = np.exp(P['log_ahat'])
ed = P['ed']; LamE = (P['y']-k) - ed['zetaE']*P['U'] + ed['etaE']*P['p']
print(f"mu={mu} rho={rho} k={k:.4f} lam={lam:.4f} A={A} T={P['T']:.4f} lam_c={P['lam_M']:.4e} etaE={ed['etaE']:.4f} zetaE={ed['zetaE']:.4f}")
print(f"{'t/T':>6} {'U':>8} {'p':>8} {'y/k':>7} {'x/k':>8} {'eta':>8} {'eta~':>8} {'d':>9} {'zeta':>8} {'zeta~':>8} {'e':>9} {'V~=2p+A-lam':>11} {'beta':>7} {'kappa':>8} {'ahat':>8} {'LamE':>9}")
for i in idx:
    print(f"{P['t'][i]/P['T']:6.3f} {P['U'][i]:8.3f} {P['p'][i]:8.4f} {P['y'][i]/k:7.4f} {P['x'][i]/k:8.4f} {P['eta'][i]:8.4f} {et[i]:8.4f} {P['eta'][i]-et[i]:+9.4f} {P['zeta'][i]:8.4f} {zt[i]:8.4f} {P['zeta'][i]-zt[i]:+9.4f} {2*P['p'][i]+A-lam:+11.4f} {P['beta'][i]:7.4f} {P['kappa'][i]:8.4f} {a[i]:8.4f} {LamE[i]:9.5f}")
