import numpy as np, time, warnings
warnings.filterwarnings('ignore')
from whlib import *
for mu in [0.01,0.05,0.2,1.0]:
    guess=None
    for rho in [1.02,1.1,1.3,1.6,2,3,5,8,13,20,35,60,100]:
        k=rho*kH(mu); t0=time.time()
        try:
            c=find_cycle(mu,k,guess=guess); guess=(c['U0'],c['y0'])
            pr=cycle_profile(mu,k,c,n=2001)
            Umean=np.trapezoid(pr['U'],pr['t'])/c['T']
            print(f"mu={mu:<5} rho={rho:<6} T={c['T']:9.3f} Umax={c['U0']:7.3f} Umin={pr['U'].min():9.3f} m_t={c['m_t']:.3e} lam_c={c['lambda_c']:+.4e} <U>+log rho={Umean+np.log(rho):+.3e} res={c['residual']:.0e} [{time.time()-t0:.1f}s]",flush=True)
        except Exception as e:
            print(f"mu={mu} rho={rho} FAILED {type(e).__name__}: {e} [{time.time()-t0:.1f}s]",flush=True); guess=None
