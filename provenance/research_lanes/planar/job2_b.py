#!/usr/bin/env python3
"""Job 2(b): Lienard plane U' = P(U,s), s' = -g(U), g = k(e^U-1), P_s = 1/(1+eta), P_U = (zeta-mu)/(1+eta),
P_UU = (1,P_U) D^2h (1,P_U)^T/(1+eta).  Tests on the attractor disc:
  B0  D^2h symmetric (consistency of the second-order data) and positive semidefinite (convexity)
  B1  generalised Zhang condition  d/dU|_s (P_U/g) <= 0, i.e.  Z := P_UU g - P_U g' <= 0  for U != 0
  B2  sign of P_U on {U>0}  (zeta > mu there?)
"""
import sys, json
import numpy as np
from job2_lib import sample
from whlib import kH

def run(mu, rho):
    k = rho*kH(mu)
    S = sample(mu, k, np.linspace(0, 2*np.pi, 8, endpoint=False))
    U, eta, zeta = S['U'], S['eta'], S['zeta']
    huu, huv, hvv = S['huu'], S['huv'], S['hvv']
    scale = np.abs(huu)+np.abs(hvv)+2*np.abs(huv)+1e-300
    res = dict(mu=mu, rho=rho, n=int(len(U)))
    res['B0_max_rel_asym'] = float(np.max(np.abs(S['asym'])/scale))
    tr = huu+hvv; det = huu*hvv-huv**2
    res['B0_min_trace'] = float(tr.min()); res['B0_min_det_over_tr2'] = float(np.min(det/tr**2))
    PU = (zeta-mu)/(1+eta)
    PUU = (huu+2*PU*huv+PU**2*hvv)/(1+eta)
    g = k*(np.exp(U)-1); gp = k*np.exp(U)
    Z = PUU*g-PU*gp
    m = np.abs(U) > 1e-3
    pos = m & (U > 0); neg = m & (U < 0)
    res['B1_frac_violation_Upos'] = float(np.mean(Z[pos] > 0)); res['B1_frac_violation_Uneg'] = float(np.mean(Z[neg] > 0))
    # normalised size of violations
    nz = Z/(np.abs(PUU*g)+np.abs(PU*gp)+1e-300)
    res['B1_max_normalised_Z_Upos'] = float(nz[pos].max()); res['B1_max_normalised_Z_Uneg'] = float(nz[neg].max())
    res['B2_frac_PU_negative_on_Upos'] = float(np.mean(PU[pos] < 0))
    res['B2_min_zeta_minus_mu_on_Upos'] = float((zeta-mu)[pos].min())
    # where are the violations (U range)
    if (Z[neg] > 0).any():
        res['B1_violation_U_range_neg'] = [float(U[neg][Z[neg] > 0].min()), float(U[neg][Z[neg] > 0].max())]
    if (Z[pos] > 0).any():
        res['B1_violation_U_range_pos'] = [float(U[pos][Z[pos] > 0].min()), float(U[pos][Z[pos] > 0].max())]
    res['U_range'] = [float(U.min()), float(U.max())]
    return res

if __name__ == '__main__':
    out = []
    for mu, rho in [(1.0, 1.1), (1.0, 1.5), (1.0, 2.0), (0.3, 1.5), (0.1, 1.5), (3.0, 1.3), (3.0, 2.0)]:
        r = run(mu, rho); out.append(r); print(json.dumps(r), flush=True)
    json.dump(out, open('job2_b.json', 'w'), indent=1)
