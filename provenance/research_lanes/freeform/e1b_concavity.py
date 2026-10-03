#!/usr/bin/env python3
"""E1b: is the 1-D crest return map P1 on W^u(E) concave (U-coordinate)?  Exact P1'(u_n) via
P1' = (p'(0)/p'(tau)) * exp(int (eta-mu)) with eta the pull-back adjoint ratio on W^u(E)."""
import sys; sys.path.insert(0,'/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03/work')
import numpy as np
from scipy.integrate import solve_ivp
from whlib import *

def run(mu,k,r0=1e-3):
    ed=equilibrium_data(mu,k)
    v=np.real(ed['vec_u']); z0=np.array([0,0,k])+r0*v/np.linalg.norm(v)
    def f(t,s):
        U,p,y,eta,zeta,I,Kg,Ir=s
        x=k*np.exp(U)
        K=zeta*U-eta*p-y+k
        om=p*p-U*(k-mu*p-y)
        return [p,k-mu*p-y,x-y,zeta+(mu-1)*eta-eta*eta,x-(1+eta)*zeta,eta-mu,0.0,p*K/om]
    ev=lambda t,s:s[1]
    ev.direction=-1
    tend=(np.log(1/r0)+25)/ed['lam_u'].real
    s0=np.concatenate((z0,[ed['etaE'],ed['zetaE'],0,0,0]))
    sol=solve_ivp(f,[0,tend],s0,method='DOP853',rtol=1e-10,atol=1e-13,events=ev)
    Y=sol.y_events[0]
    u=Y[:,0]; y=Y[:,2]; I=Y[:,5]; Ir=Y[:,7]
    dP=(y[:-1]-k)/(y[1:]-k)*np.exp(np.diff(I))   # P1'(u_n)
    Irad=np.diff(Ir)                               # int pK/om over the turn
    return u,dP,Irad,ed

for mu in (0.1,1.0):
    for rho in (1.2,2.0,4.0):
        k=rho*kH(mu)
        u,dP,Irad,ed=run(mu,k)
        m=(u[:-1]>1e-4)&(np.abs(u[1:]/u[:-1]-1)>1e-8)
        d=dP[m]; uu=u[:-1][m]
        conc=np.all(np.diff(d)<1e-10)
        # check consistency: u1/u0 * exp(-Irad) == P1'
        chk=np.max(np.abs(d-(u[1:]/u[:-1])[m]*np.exp(-Irad[m])))
        print(flush=True); print(f"mu={mu} rho={rho} n={m.sum():4d} P1'(first)={d[0]:.4f} P1'(last)={d[-1]:.3e} concave(P1' decreasing)={conc} "
              f"min Irad={Irad[m].min():.3e} chk={chk:.1e}")
        if not conc:
            j=np.where(np.diff(d)>=1e-10)[0]
            print("   violations at u=",uu[j][:5],"dP=",d[j][:5],d[j+1][:5])
