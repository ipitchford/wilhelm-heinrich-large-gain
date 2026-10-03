#!/usr/bin/env python3
"""E3: star-shapedness omega1 = p^2 - U p' > 0 on cycles; Massera identity lambda_c = -<pK/omega1>;
upper/lower half contributions; mirror comparison data for the 'next step' in NOTES section 5."""
import sys; sys.path.insert(0,'/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03/work')
import numpy as np
from whlib import *
for mu,rho in ((0.05,1.3),(0.1,1.3),(0.3,1.05),(1.0,1.2),(1.0,2.0),(0.2,2.0),(1.0,5.0),(0.2,8.0),(3.0,2.0)):
    k=rho*kH(mu); c=find_cycle(mu,k); pr=cycle_profile(mu,k,c,n=8001); T=c['T']; t=pr['t']
    U,p,y,K=pr['U'],pr['p'],pr['y'],pr['K']; pp=k-mu*p-y
    om=p*p-U*pp
    av=lambda g: np.trapezoid(g,t)/T
    integrand=p*K/om
    up=np.trapezoid(np.where(p>0,integrand,0),t)/T; lo=np.trapezoid(np.where(p<0,integrand,0),t)/T
    print(f"mu={mu} rho={rho}: min omega1/max omega1={om.min()/om.max():.3e}  lambda_c={pr['lambda_c_riccati']:.6f}  -<pK/om>={-av(integrand):.6f}"
          f"  upper(p>0)={up:.5f} lower(p<0)={lo:.5f}  minK={K.min():.2e}")
