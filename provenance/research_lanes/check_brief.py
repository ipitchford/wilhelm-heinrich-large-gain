import numpy as np
from whlib import *
tz=np.trapezoid
for mu,rho in [(1.0,1.1),(0.3,1.5),(0.2,3.0),(3.0,1.4)]:
    k=rho*kH(mu); A=1+mu; lam=k/mu
    c=find_cycle(mu,k); pr=cycle_profile(mu,k,c,n=20001); t=pr['t']; T=c['T']
    U,p,y,x,eta,zeta,K,W=[pr[s] for s in 'U p y x eta zeta K W'.split()]
    av=lambda f: tz(f,t)/T
    ed=equilibrium_data(mu,k)
    LE=(y-k)-ed['zetaE']*U+ed['etaE']*p
    print(f"mu={mu} rho={rho} lam_c(M)={c['lambda_c']:+.8f}  <eta>-mu={av(eta)-mu:+.8f}  <x/zeta>-A={av(x/zeta)-A:+.8f}  <zeta/eta>-1={av(zeta/eta)-1:+.8f}  <kTheta/K>-A={av(k*(U*np.exp(U)-np.exp(U)+1)/K)-A:+.8f}  <xp^2/W>-A={av(x*p*p/W)-A:+.8f}")
    q=eta-mu
    print(f"   <y/(lam-p)>-mu={av(y/(lam-p))-mu:+.2e}  <x/y>-1={av(x/y)-1:+.2e}  minLambda_E={LE.min():.3e}  <kR/Lambda_E>-gammaE={av(k*(np.exp(U)-1-U)/LE)-ed['gammaE']:+.2e}  cubic: {(A*A+mu)*av(q)+2*A*av(q*q)+av(q**3)-(k-kH(mu)):+.2e}  min(y-zeta)={np.min(y-zeta):.3e}  min(A+p)={np.min(A+p):.3f}")
