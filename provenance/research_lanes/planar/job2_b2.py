#!/usr/bin/env python3
"""Job 2(b) supplement: Zhang-type pointwise conditions in the other Lienard-type planes of the attractor.
 plane (yh=y-k, m=p+mu U):  yh' = X(y,m)-k-yh,  m' = -yh   (linear 'g').  div = x/kappa - 1, kappa = zeta+mu eta.
     X_yy = x/kappa^2 - x q/kappa^3,  q = (1,-mu) D^2h (1,-mu)^T.   Zhang:  Z2 = X_yy*yh - (X_y - 1) <= 0 ?
 plane (U,s): repeated from job2_b for reference (Z = P_UU g - P_U g').
Also the structural monotonicities implied by convexity (checked): d(eta)/dv <= 0, d(zeta)/du >= 0."""
import json
import numpy as np
from job2_lib import sample
from whlib import kH

def run(mu, rho):
    k = rho*kH(mu)
    S = sample(mu, k, np.linspace(0, 2*np.pi, 8, endpoint=False))
    U, y, x, eta, zeta = S['U'], S['y'], S['x'], S['eta'], S['zeta']
    huu, huv, hvv = S['huu'], S['huv'], S['hvv']
    kap = zeta+mu*eta; q = huu-2*mu*huv+mu*mu*hvv
    yh = y-k
    Xy = x/kap; Xyy = x/kap**2-x*q/kap**3
    Z2 = Xyy*yh-(Xy-1)
    m = np.abs(yh) > 1e-3*k
    res = dict(mu=mu, rho=rho, frac_Z2_violation_yh_pos=float(np.mean(Z2[m & (yh > 0)] > 0)),
               frac_Z2_violation_yh_neg=float(np.mean(Z2[m & (yh < 0)] > 0)),
               min_q=float(q.min()), min_hvv=float(hvv.min()), min_huu=float(huu.min()))
    return res

if __name__ == '__main__':
    out = []
    for mu, rho in [(1.0, 1.1), (1.0, 2.0), (0.3, 1.5), (3.0, 1.3)]:
        r = run(mu, rho); out.append(r); print(json.dumps(r), flush=True)
    json.dump(out, open('job2_b2.json', 'w'), indent=1)
