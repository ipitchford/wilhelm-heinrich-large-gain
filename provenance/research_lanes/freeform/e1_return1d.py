#!/usr/bin/env python3
"""E1: 1-D crest return map on W^u(E) (attractor disc). Test monotonicity of P1(u)/u and concavity."""
import sys; sys.path.insert(0,'/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03/work')
import numpy as np
from scipy.integrate import solve_ivp
from whlib import *

def crests(mu,k,r0=1e-6,nmax=400):
    ed=equilibrium_data(mu,k)
    v=np.real(ed['vec_u']); z0=np.array([0,0,k])+r0*v/np.linalg.norm(v)
    ev=lambda t,z: z[1]
    ev.direction=-1
    tend=(np.log(1/r0)+25)/ed['lam_u'].real
    sol=solve_ivp(field(mu,k),[0,tend],z0,method='DOP853',rtol=1e-12,atol=1e-15,events=ev)
    Y=sol.y_events[0]
    return Y[:,0],Y[:,2],ed

for mu in (0.05,0.3,1.0):
    for rho in (1.05,1.3,2.0,4.0):
        k=rho*kH(mu)
        u,y,ed=crests(mu,k)
        u=u[u>1e-5]
        if len(u)<5: print(mu,rho,'few crests',len(u)); continue
        R=u[1:]/u[:-1]
        # P1(u_n)=u_{n+1}; ratio R_n vs u_n
        dR=np.diff(R)
        # stop once converged
        m=np.abs(R-1)>1e-9
        Rm=R[m]
        mono=np.all(np.diff(Rm)<1e-12)
        # also ratio in x-1 coordinate
        xx=np.exp(u)-1; Rx=xx[1:]/xx[:-1]; monox=np.all(np.diff(Rx[m])<1e-12)
        yy=(y[-len(u):]-k); Ry=yy[1:]/yy[:-1]; monoy=np.all(np.diff(Ry[m])<1e-12)
        print(f"mu={mu} rho={rho} ncrest={len(u)} R0={R[0]:.5f} linear={np.exp(2*np.pi*ed['lam_u'].real/ed['lam_u'].imag):.5f} "
              f"mono_U={mono} mono_x={monox} mono_y={monoy} u_last={u[-1]:.6f} n_used={m.sum()}")
