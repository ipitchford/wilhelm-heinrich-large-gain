#!/usr/bin/env python3
"""E2: profiles of eta vs eta~=y/q, Perron function phi*=exp(int(eta-<eta>)), dual Lienard plane data."""
import sys; sys.path.insert(0,'/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03/work')
import numpy as np
from whlib import *
np.set_printoptions(precision=4,suppress=True,linewidth=200)
for mu,rho in ((1.0,1.2),(1.0,2.0),(0.2,2.0),(1.0,5.0)):
    k=rho*kH(mu); lam=k/mu; A=1+mu
    c=find_cycle(mu,k); pr=cycle_profile(mu,k,c,n=2401)
    t,U,p,y,x,eta,zeta,K,W=[pr[s] for s in 't U p y x eta zeta K W'.split()]
    T=c['T']; q=lam-p
    avg=lambda f: np.trapezoid(f,t)/T
    et=y/q; zt=x/q
    lamc=avg(eta)-mu
    print(f"\nmu={mu} rho={rho} T={T:.4f} lam_c={lamc:.5f} m_t={c['m_t']:.4f} <eta~>-mu={avg(et)-mu:.2e}")
    print(f"  ranges: U[{U.min():.3f},{U.max():.3f}] p[{p.min():.3f},{p.max():.3f}] lam={lam:.3f} eta[{eta.min():.3f},{eta.max():.3f}] eta~[{et.min():.3f},{et.max():.3f}]")
    # extrema counts of eta, zeta, sigma
    sig=zeta+eta+eta**2/2
    def next(f): 
        d=np.diff(f); return int(np.sum(d[1:]*d[:-1]<0))
    print("  #extrema: eta",next(eta),"zeta",next(zeta),"sigma",next(sig),"eta~",next(et),"x/zeta",next(x/zeta),"zeta/eta",next(zeta/eta))
    # identities
    f=lambda s:s*(s+1)*(s+1-mu)
    print(f"  <eta x>-<eta f(eta)>={avg(eta*x)-avg(eta*f(eta)):.4f} (pred <0)   <x/zeta>-A={avg(x/zeta)-A:.5f}  <zeta/eta>-1={avg(zeta/eta)-1:.5f}")
    print(f"  pointwise: frac(zeta/eta<x/y)={np.mean(zeta/eta<x/y):.3f} frac(eta<eta~)={np.mean(eta<et):.3f} frac(x/zeta<q)={np.mean(x/zeta<q):.3f} frac(x/zeta<A)={np.mean(x/zeta<A):.3f}")
    idx=np.linspace(0,len(t)-1,13).astype(int)[:-1]
    phi=np.exp(np.concatenate(([0],np.cumsum((eta[1:]+eta[:-1])/2*np.diff(t))))-avg(eta)*t)
    print("  t/T   ",t[idx]/T)
    print("  U     ",U[idx]); print("  p     ",p[idx]); print("  y/k   ",y[idx]/k)
    print("  eta   ",eta[idx]); print("  eta~  ",et[idx]); print("  zeta/k",zeta[idx]/k); print("  zeta~/k",zt[idx]/k)
    print("  phi*  ",phi[idx]/phi.mean()); print("  phi*/q",(phi/q)[idx]/np.mean(phi/q))
    print("  x/zeta",(x/zeta)[idx])
