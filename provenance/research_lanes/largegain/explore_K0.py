#!/usr/bin/env python3
"""Exploration: which K0 does the a-priori chain support?  (design aid; the final check is constants.py)"""
import numpy as np
from scipy.optimize import brentq, minimize_scalar

ELL = 2.0; PHI = 1-np.exp(-ELL); TH = 2**-0.5
v2 = brentq(lambda v: v-2*(1-np.exp(-v)), 1, 2); PHI2 = v2/2
Psi1 = lambda u: 1-(1+u)*np.exp(-u)
ci = lambda a: (np.exp(a)-np.exp(a/2))*(1-(1+a/2)*np.exp(-a/2))
cii = lambda a: a*a/2-1+(1+a)*np.exp(-a)
cb = lambda a: Psi1(a/2)*(1-np.exp(-a/2))
om0 = lambda k: 0.5*np.sqrt(3*k**(2/3)-4*k**(1/3)-1)

def stage1(K0, verbose=False):
    T0 = 2*np.pi - np.pi/om0(K0)
    a4 = brentq(lambda a: K0*cb(a)-ELL/PHI, 1e-3, 3)*1.0001
    S = T0 - a4
    best = None
    for a1 in np.linspace(0.8, 3.0, 221):
        if K0*ci(a1) < 2: continue
        Ri2 = 2*(PHI2*ci(a1)-ELL/K0)
        if Ri2 <= 0: continue
        Ri = np.sqrt(Ri2)
        for a2 in np.linspace(0.8, 3.0, 221):
            if K0*(a2-1+np.exp(-a2))/(K0**(1/3)+1) < 2: continue
            Rii2 = 2*(PHI2*cii(a2)-ELL/K0)
            if Rii2 <= 0: continue
            a3 = S-a1-a2
            if a3 <= 0: continue
            r = min(a3, Ri, np.sqrt(Rii2))
            if best is None or r > best[0]:
                best = (r, a1, a2, a3, a4)
    r, a1, a2, a3, a4 = best
    c1 = PHI*Psi1(r) - ELL/K0
    return dict(T0=T0, a=(a1, a2, a3, a4), rmin=r, c1=c1)

def chain(K0, verbose=True):
    s1 = stage1(K0)
    c1 = s1['c1']
    Ht = lambda k: 2*np.log(k+1)+3/np.e
    eps1 = 1.1*np.sqrt(2)/(c1**1.5*np.sqrt(K0))
    d1 = 4*Ht(K0)*(1+eps1)/(c1*K0) + 4*np.sqrt(2)/(c1**3*K0**2)
    if d1 >= 1: return None
    T1 = c1**2*K0*np.sqrt(1-d1)/np.sqrt(8) - 0.5
    kap = (0.94*PHI)**2*np.sqrt(1-d1)/np.sqrt(8)
    z1 = 0.3*T1-3
    ok_T1 = T1 >= 22
    # endpoint conditions
    cond_z1 = None
    if z1 > 0:
        cond_z1 = (0.3*kap - 0.3*c1**2/np.sqrt(8)*(1+1/z1)**2)*K0 - 3.15*(1+1/z1)**2
    csharp = 0.27*kap
    cond_zs = 1/0.9 - (1+3.15/(csharp*K0))*(1+1/(csharp*K0))**2
    cT = csharp/0.3
    cp = 0.95*PHI*(1-np.exp(-csharp*K0)) - ELL/K0
    out = dict(K0=K0, **s1, eps1=eps1, d1=d1, T1=T1, kap=kap, z1=z1, ok_T1=ok_T1, cond_z1=cond_z1, cond_zs=cond_zs, cT=cT, cp=cp)
    # durations at lam = K0 (worst case), k<=lam
    lam = K0
    Hh = np.log((lam*lam/K0+1)*(lam+1)); Lh = np.log(8*lam/3)
    pk = cp*lam
    e1 = 1.1*np.sqrt(2)*lam/pk**1.5
    dE = 4*(1+e1)*Hh/pk + 4*lam/(TH*pk**3)
    A1 = np.sqrt(2)*Hh + 3.4641/np.sqrt(1-dE) + 16.97/(cp**2*lam*(1-dE))
    tau1 = A1/pk
    W = pk**2*np.exp(-tau1)/A1
    f = lambda v: (W/2)*v*v*(1-(2/3)*(1+lam/W)*v)-(Hh+Lh)
    vf = brentq(f, 1e-9, 0.5/(1+lam/W))
    tau_c = 2*Lh/pk + tau1 + vf
    s = 2*np.sqrt(lam)
    final = cT*lam/4 - (2*np.sqrt(lam)*tau_c + np.log(1+4*np.sqrt(lam)+16*lam/3))
    out.update(Hh=Hh, Lh=Lh, dE=dE, A1=A1, tau1=tau1, W=W, vf=vf, tau_c=tau_c, final=final, C_tau=tau_c*pk)
    if verbose:
        def fmt(b):
            if isinstance(b, tuple): return tuple(round(float(x), 4) for x in b)
            if isinstance(b, (bool, np.bool_)) or b is None: return b
            return round(float(b), 5)
        print({a: fmt(b) for a, b in out.items()})
    return out

if __name__ == '__main__':
    for K0 in (300, 400, 500, 600, 700, 800, 1000, 1500, 2000, 5000):
        try:
            chain(K0)
        except Exception as e:
            print(K0, 'ERR', e)
