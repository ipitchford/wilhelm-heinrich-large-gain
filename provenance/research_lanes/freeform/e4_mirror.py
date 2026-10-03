#!/usr/bin/env python3
"""E4: -lambda_c T = int_{Umin}^{Umax} [K_+/om_+ - K_-/om_-] dU  (upper arc p>0 vs lower arc p<0 at equal U).
Test: is the bracket pointwise >= 0 in U?  Also the same with other fibre choices (equal polar angle mirror)."""
import sys; sys.path.insert(0,'/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03/work')
import numpy as np
from whlib import *
for mu,rho in ((0.05,1.3),(0.3,1.05),(1.0,1.2),(1.0,2.0),(0.2,2.0),(1.0,5.0),(3.0,2.0)):
    k=rho*kH(mu); c=find_cycle(mu,k); pr=cycle_profile(mu,k,c,n=20001); T=c['T']; t=pr['t']
    U,p,y,K=pr['U'],pr['p'],pr['y'],pr['K']; pp=k-mu*p-y; om=p*p-U*pp
    g=K/om
    lo=p<0; up=p>0
    # lower arc: U decreasing in time; upper arc: U increasing
    Ul,gl=U[lo][::-1],g[lo][::-1]; Uu,gu=U[up],g[up]
    Ug=np.linspace(U.min(),U.max(),2003)[1:-1]
    gL=np.interp(Ug,Ul,gl); gU=np.interp(Ug,Uu,gu)
    d=gU-gL
    I=np.trapezoid(d,Ug)
    neg=Ug[d<0]
    print(f"mu={mu} rho={rho}: int bracket dU={I:.5f} vs -lambda_c T={-pr['lambda_c_riccati']*T:.5f};  frac(U) with bracket<0 = {np.mean(d<0):.3f}"
          + (f"  negative on U in [{neg.min():.3f},{neg.max():.3f}] (Umin={U.min():.2f},Umax={U.max():.2f}) min bracket={d.min():.3e} max={d.max():.3e}" if len(neg) else "  bracket >= 0 everywhere"))
