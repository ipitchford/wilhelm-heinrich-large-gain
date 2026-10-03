import sys, numpy as np
from base import load
D = load()
mu = float(sys.argv[1]); rho = float(sys.argv[2]); m = int(sys.argv[3]) if len(sys.argv) > 3 else 24
P = D[(mu, rho)]; k = P['k']; lam = k/mu; A = 1+mu; n = len(P['t'])
U, p, y, x, eta, zeta = P['U'], P['p'], P['y'], P['x'], P['eta'], P['zeta']
ed = P['ed']; H = y-k; pd = k-mu*p-y; yd = x-y
S1 = p*p-U*pd; K = zeta*U-eta*p-H; LamE = H-ed['zetaE']*U+ed['etaE']*p
r = np.hypot(U,p); th = np.arctan2(p,U)
idx = np.arange(0, n, n//m)
print(f"mu={mu} rho={rho} lam_c={P['lam_M']:.4e}  -<pK/S1>={-np.mean(p*K/S1):.4e}; split: p>0 part {np.sum((p*K/S1)[p>0])/n:.4e}, p<0 part {np.sum((p*K/S1)[p<0])/n:.4e}")
print(f"{'t/T':>6} {'theta/pi':>8} {'U':>8} {'p':>8} {'S1':>9} {'K':>9} {'LamE':>9} {'K/r':>8} {'K/S1':>9} {'pK/S1':>9} {'eta-mu':>8} {'DlogS1':>8}")
for i in idx:
    print(f"{P['t'][i]/P['T']:6.3f} {th[i]/np.pi:8.3f} {U[i]:8.3f} {p[i]:8.4f} {S1[i]:9.4f} {K[i]:9.5f} {LamE[i]:9.5f} {K[i]/r[i]:8.4f} {K[i]/S1[i]:9.4f} {p[i]*K[i]/S1[i]:+9.4f} {eta[i]-mu:+8.4f} {eta[i]-mu+p[i]*K[i]/S1[i]:+8.4f}")
