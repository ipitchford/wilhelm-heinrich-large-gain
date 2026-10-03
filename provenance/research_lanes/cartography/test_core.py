import sys, time
sys.path.insert(0, '/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03/work')
sys.path.insert(0, '/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03/work/lanes/cartography')
import numpy as np
import cart as C
t0=time.time()
ref = [(0.1, 0.196052776, 38.028889, 0.43642544581),
       (1.0, 2.201497, 6.502509, 0.77861840302),
       (1.0, 5.841219, 11.121221, 0.030745762),
       (10.0, 196.052776, 3.802889, 0.43642544581)]
for mu,k,Tref,mref in ref:
    # iterate crest map from a rough start
    ed = C.equilibrium_data(mu,k)
    z0, _ = C.unstable_seed(mu,k,1e-2)
    U,Y,tt,st = C.crest_iter(mu,k,0.5,k*1.3,60)
    c = C.newton_fixed_k(mu,k,U[-1],Y[-1])
    t1=time.time()
    pr = C.profile(mu,k,c['U0'],c['y0'],c['T'])
    print(f"mu={mu} k={k} T={c['T']:.6f} (ref {Tref}) m_t={c['m_t']:.11f} (ref {mref}) res={c['residual']:.1e} it={c['iters']}"
          f" lamT_M={np.log(c['m_t']):.10f} lamT_R={pr['log_mt']:.10f} settle={pr['settle_periods']} defect={pr['ratio_defect']:.1e} closure={pr['closure']:.1e}"
          f" <x>/k={pr['mean_x']/k:.12f} Umin={pr['U_min']:.4f} minK={pr['K'].min():.2e} minW={pr['W'].min():.2e} n={len(pr['t'])} t={time.time()-t1:.2f}s")
print('total',time.time()-t0)
