#!/usr/bin/env python3
"""Joint concave programme: common EXPONENTS over all rho at fixed mu, but free per-cycle prefactors
(l1,l2) (these can be absorbed in explicit k-dependent constants such as eta_E, zeta_E).
maximise  min_rho S_rho/(-lambda_c)."""
import sys
import numpy as np
from scipy.optimize import minimize
from base import load
from cvx import make_blocks, soft_score_grad, score
D = load()
feats = tuple(sys.argv[1].split(','))
rmin = float(sys.argv[2]) if len(sys.argv) > 2 else 0.0
rmax = float(sys.argv[3]) if len(sys.argv) > 3 else 1e9
nf = 1 + len(feats); ne = nf - 1
print("features:", feats, " rho range", rmin, rmax, "(common exponents, free prefactors)")
def unpack(q, j, nc):
    ex = q[:2*ne].reshape(2, ne); cs = q[2*ne:].reshape(nc, 2)
    return np.concatenate([[cs[j, 0]], ex[0], [cs[j, 1]], ex[1]])
for mu in sorted(set(k[0] for k in D)):
    keys = [k for k in sorted(D) if k[0] == mu and rmin <= k[1] <= rmax]
    nc = len(keys)
    Bs = [make_blocks(D[k], feats, nsub=4096) for k in keys]
    norms = [abs(D[k]['lam_M']) for k in keys]
    q = np.zeros(2*ne + 2*nc); q[2*ne:] = np.log(mu)
    for s in (0.3, 0.1, 0.03, 0.01, 0.003, 0.001, 3e-4):
        def f(q):
            vals = []; grads = []
            for j, (B, nm) in enumerate(zip(Bs, norms)):
                v, g = soft_score_grad(unpack(q, j, nc), B, s*nm)
                gq = np.zeros_like(q)
                gq[:ne] = g[1:nf]; gq[ne:2*ne] = g[nf+1:]
                gq[2*ne + 2*j] = g[0]; gq[2*ne + 2*j + 1] = g[nf]
                vals.append(v/nm); grads.append(gq/nm)
            vals = np.array(vals); grads = np.array(grads)
            m = vals.min(); E = np.exp(-(vals-m)/s); Z = E.sum()
            return -(m - s*np.log(Z)), -(E/Z) @ grads
        r = minimize(f, q, jac=True, method='L-BFGS-B', options=dict(maxiter=5000, maxfun=10000, ftol=1e-15, gtol=1e-10))
        q = r.x
    vals = np.array([score(unpack(q, j, nc), B)/nm for j, (B, nm) in enumerate(zip(Bs, norms))])
    print(f"mu={mu:5.2f}  min ratio = {vals.min():+8.4f}   ratios: " + ' '.join(f"{v:+.3f}" for v in vals))
    print("      exps eta~:", np.array2string(q[:ne], precision=3), " zeta~:", np.array2string(q[ne:2*ne], precision=3))
