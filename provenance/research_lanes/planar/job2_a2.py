#!/usr/bin/env python3
"""Job 2(a), exact formula for the derivative of the turn integral iota(u) on the crest section:
   d iota/du = int_0^tau  w(t) * Xi / |F|^4 dt,   w(t) = g~(u) exp(int_0^t (eta-mu)),  g~(u) = h(u,0)-k = -v'(u,0),
   Xi = (h_uv v' - h_vv v)(v^2+v'^2) - delta[(1-zeta)(v^2-v'^2) + 2 delta v v'],  delta = eta - mu.
Checked against finite differences of iota over neighbouring spirals.  Also reports the sign pattern of Xi."""
import sys, json
import numpy as np
from scipy.integrate import solve_ivp
sys.path.insert(0, '/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03/work')
from whlib import equilibrium_data, kH
from job2_lib import rhs
from job2_a import turns

def ev_c(t, s): return s[1]
ev_c.direction = -1

def run(mu, rho, theta=0.3):
    k = rho*kH(mu); ed = equilibrium_data(mu, k)
    r0 = 1e-6
    tend = (np.log(1/r0)+12)/ed['lam_u'].real
    vc = np.exp(1j*theta)*ed['vec_u']; nv = np.linalg.norm(np.real(vc))
    z0 = np.array([0, 0, k])+r0*np.real(vc)/nv; w0 = r0*np.real(1j*vc)/nv
    s0 = np.concatenate((z0, [ed['etaE'], ed['zetaE'], 0, 0], w0, [0, 0]))
    sol = solve_ivp(rhs(mu, k), [0, tend], s0, method='DOP853', rtol=1e-11, atol=1e-18, events=ev_c, dense_output=True)
    tc = sol.t_events[0]; yc = sol.y_events[0]
    rows = []
    for i in range(len(tc)-1):
        u = yc[i][0]
        if u < 0.02: continue
        t = np.linspace(tc[i], tc[i+1], 6001)
        U, p, y, eta, zeta, K, W, dU, dp, dy, de, dz = sol.sol(t)
        x = k*np.exp(U); pp = k-mu*p-y
        etap = zeta+(mu-1)*eta-eta**2; zetap = x-(1+eta)*zeta
        det = p*dp-pp*dU
        G12 = (zetap*(-dU)+dz*p)/det; G21 = (-etap*dp+(-de)*(-pp))/det; G22 = (-etap*(-dU)+(-de)*p)/det
        huv = 0.5*(G12+G21); hvv = G22
        delta = eta-mu
        I = np.concatenate(([0], np.cumsum(0.5*(delta[1:]+delta[:-1])*np.diff(t))))
        gt = yc[i][2]-k
        w = gt*np.exp(I)
        Xi = (huv*pp-hvv*p)*(p*p+pp*pp)-delta*((1-zeta)*(p*p-pp*pp)+2*delta*p*pp)
        integrand = w*Xi/(p*p+pp*pp)**2
        val = np.trapezoid(integrand, t)
        rows.append((u, yc[i+1][0], I[-1], val, float(np.mean(Xi > 0))))
    return np.array(rows)

if __name__ == '__main__':
    out = []
    for mu, rho in [(1.0, 1.3), (1.0, 2.0), (0.3, 1.5), (3.0, 1.3)]:
        R = run(mu, rho)
        C, Tr, HC, HT, ed = turns(mu, rho*kH(mu), nth=16)
        o = np.argsort(C[:, 0]); uu = C[o, 0]; ii = C[o, 2]
        print(f"mu={mu} rho={rho}")
        for u, u1, iota, form, fpos in R[:: max(1, len(R)//7)]:
            # finite-difference slope of iota at u from merged spirals (local linear fit over nearest 6 points)
            j = np.argsort(np.abs(uu-u))[:6]
            if np.ptp(uu[j]) < 1e-6: continue
            fd = np.polyfit(uu[j], ii[j], 1)[0]
            print(f"   u={u:.4f}  iota={iota:+.5f}  d iota/du: formula={form:+.5f}  finite-diff={fd:+.5f}   frac(Xi>0 on turn)={fpos:.2f}")
            out.append(dict(mu=mu, rho=rho, u=float(u), iota=float(iota), formula=float(form), fd=float(fd), frac_Xi_pos=float(fpos)))
    json.dump(out, open('job2_a2.json', 'w'), indent=1)
