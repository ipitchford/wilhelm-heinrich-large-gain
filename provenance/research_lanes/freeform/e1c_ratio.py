#!/usr/bin/env python3
"""E1c: successive-crest ratio R_n=u_{n+1}/u_n and secant slopes on W^u(E); extra parameter points.
Also second differences (secant slope s_n=(u_{n+2}-u_{n+1})/(u_{n+1}-u_n) decreasing <=> consistent with concave P1)."""
import sys,time; sys.path.insert(0,'/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03/work')
import numpy as np
from scipy.integrate import solve_ivp
from whlib import *
def crests(mu,k,r0=1e-5,theta=0.0):
    ed=equilibrium_data(mu,k)
    v=np.real(np.exp(1j*theta)*ed['vec_u']); z0=np.array([0,0,k])+r0*v/np.linalg.norm(v)
    ev=lambda t,z: z[1]
    ev.direction=-1
    tend=(np.log(1/r0)+22)/ed['lam_u'].real
    sol=solve_ivp(field(mu,k),[0,tend],z0,method='DOP853',rtol=1e-12,atol=1e-16,events=ev)
    return sol.y_events[0][:,0]
t0=time.time()
for mu,rho in ((0.1,2.0),(0.1,1.3),(0.02,1.5),(0.02,3.0),(0.5,1.15),(0.5,6.0),(1.0,10.0),(0.2,8.0)):
    k=rho*kH(mu)
    for th in (0.0,1.3):
        u=crests(mu,k,theta=th); u=u[u>3e-5]
        R=u[1:]/u[:-1]; m=np.abs(R-1)>1e-9
        s=np.diff(u)[1:]/np.diff(u)[:-1]; ms=m[1:]
        okR=np.all(np.diff(R[m])<0); okS=np.all(np.diff(s[ms])<1e-9)
        print(f"mu={mu} rho={rho} th={th} n={m.sum()} R:[{R[m][0]:.5f}->{R[m][-1]:.7f}] R decreasing={okR}  secant slopes decreasing={okS}  maxdR={np.diff(R[m]).max():.2e} t={time.time()-t0:.0f}s",flush=True)
