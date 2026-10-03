#!/usr/bin/env python3
"""Exact (sympy) checks of the identities recorded in NOTES.md, plus a floating check of the Hoelder sandwich."""
import sys; sys.path.insert(0,'/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03/work')
import sympy as sp
mu,k,t=sp.symbols('mu k t',positive=True); lam=k/mu; A=1+mu
U,p,y,eta,zeta=[sp.Function(n)(t) for n in ('U','p','y','eta','zeta')]
x=k*sp.exp(U); q=lam-p
flow={sp.Derivative(U,t):p, sp.Derivative(p,t):k-mu*p-y, sp.Derivative(y,t):x-y,
      sp.Derivative(eta,t):zeta+(mu-1)*eta-eta**2, sp.Derivative(zeta,t):x-(1+eta)*zeta}
def D(e,n=1):
    for _ in range(n): e=sp.diff(e,t).subs(flow)
    return sp.simplify(e)
ok=lambda name,e: print(f"{name:55s}", "OK" if sp.simplify(e)==0 else f"FAIL {sp.simplify(e)}")
# (1) Lotka-Volterra form
r1,r2=x/y,y/q
ok("LV r1' = r1(lam+1-r1-q)", D(r1)-r1*(lam+1-r1-q))
ok("LV r2' = r2(mu-1+r1-r2)", D(r2)-r2*(mu-1+r1-r2))
ok("LV q'  = q(r2-mu)", D(q)-q*(r2-mu))
# (2) second-order Riccati (Burgers-hierarchy) equation for eta
f=lambda s:s*(s+1)*(s+1-mu)
ok("eta''+(3eta+2-mu)eta'+f(eta) = x", D(eta,2)+(3*eta+2-mu)*D(eta)+f(eta)-x)
et=y/q
ok("same operator on eta~=y/q gives x(A+p)/q", D(et,2)+(3*et+2-mu)*D(et)+f(et)-x*(A+p)/q)
# (3) dual Lienard plane
sig=zeta+eta+eta**2/2; G=sp.Rational(3,2)*eta**2+(2-mu)*eta
ok("eta' = sigma-G(eta)", D(eta)-(sig-G)); ok("sigma' = x-f(eta)", D(sig)-(x-f(eta)))
# (4) Massera identity (needs tangent orthogonality zeta p - eta p' = x-y: substitute zeta)
K=zeta*U-eta*p-y+k; om=p**2-U*(k-mu*p-y); pp=k-mu*p-y
expr=D(om)-((eta-mu)*om+p*K)
expr=expr.subs(zeta,(x-y+eta*pp)/p)
ok("(p^2-Up')' = (eta-mu)(p^2-Up') + pK  [on O=0]", expr)
# (5) backward (C-) operator on phi=x and compound operator on 1, q
Lm=lambda ph: -D(ph,3)+2*A*D(ph,2)-(A**2+mu)*D(ph)+(mu*A-x)*ph
Lc=lambda ph: D(ph,3)+2*A*D(ph,2)+(A**2+mu)*D(ph)+(mu*A-x)*ph
ok("Lc^-[x] = x[(A-p)(3p'-p(A-p)) - (k-kH)]", Lm(x)-x*((A-p)*(3*pp-p*(A-p))-(k-mu*A)))
ok("Lc[q] = x(2p+A-lam)", Lc(q)-x*(2*p+A-lam))
ok("Lc[1/x] x = -(A-3p)p' - p(A-p)^2 - (k-kH)", Lc(1/x)*x-(-(A-3*p)*pp-p*(A-p)**2-(k-mu*A)))
# floating: Hoelder sandwich  <x^(1/3)>^3 <= gamma(gamma-1)(gamma-mu) < k
import numpy as np
from whlib import *
for m_,rho in ((0.1,1.3),(1.0,1.2),(1.0,2.0),(0.2,2.0),(1.0,5.0),(0.3,1.05)):
    kk=rho*kH(m_); c=find_cycle(m_,kk); pr=cycle_profile(m_,kk,c,n=4001); T=c['T']
    av=lambda g: np.trapezoid(g,pr['t'])/T
    gam=1+av(pr['eta']); mm=gam*(gam-1)*(gam-m_)
    print(f"mu={m_} rho={rho}: exp<log x>={np.exp(av(np.log(pr['x']))):.4f} <= <x^1/3>^3={av(pr['x']**(1/3))**3:.4f} <= m={mm:.4f} < kH={kH(m_):.4f} < k={kk:.4f}   (<x p^2>/<p^2>={av(pr['x']*pr['p']**2)/av(pr['p']**2):.4f})")
