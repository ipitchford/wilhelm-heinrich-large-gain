import sys, time
sys.path.insert(0,'/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03/work')
from whlib import *
t0=time.time()
for mu,rho in [(1.0,1.1),(1.0,2.0),(0.3,1.5),(3.0,1.3)]:
    k=rho*kH(mu)
    c=find_cycle(mu,k)
    pr=cycle_profile(mu,k,c,n=2001)
    orb=unstable_manifold_orbit(mu,k,n=4001)
    print(mu,rho,k,'T',c['T'],'m_t',c['m_t'],'lam',c['lambda_c'],pr['lambda_c_riccati'],'U0',c['U0'],'Umin',pr['U'].min(),'orb Umax',orb['U'].max(),'minK',orb['K'][10:].min(),'minW',orb['W'][10:].min(), round(time.time()-t0,1))
