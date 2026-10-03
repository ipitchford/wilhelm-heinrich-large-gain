#!/usr/bin/env python3
"""Job 2 (b)/(c): Massera-type radial certificates = Dulac functions B = 1/(X ^ xi) in a chart xi centred at E.
For a chart xi=(xi1,xi2) of the attractor disc,  beta := div_xi X - D log|X ^ xi|  satisfies lambda_c T = oint beta dt
on every cycle that is star-shaped in the chart (X ^ xi != 0).  beta <= 0 on the disc minus a neighbourhood of E would
prove uniqueness (Theorem A (A8)(f)).  We evaluate beta on W^u(E) orbits for several natural charts."""
import sys, json
import numpy as np
from scipy.integrate import solve_ivp
sys.path.insert(0, '/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03/work')
from whlib import equilibrium_data, kH, find_cycle

def orbit(mu, k, theta, r0=1e-6, n=40001):
    ed = equilibrium_data(mu, k)
    tend = (np.log(1/r0)+9)/ed['lam_u'].real
    v = np.real(np.exp(1j*theta)*ed['vec_u']); z0 = np.array([0, 0, k])+r0*v/np.linalg.norm(v)
    def f(t, s):
        U, p, y, eta, zeta = s
        x = k*np.exp(U)
        return [p, k-mu*p-y, x-y, zeta+(mu-1)*eta-eta*eta, x-(1+eta)*zeta]
    t = np.linspace(0, tend, n)
    sol = solve_ivp(f, [0, tend], list(z0)+[ed['etaE'], ed['zetaE']], method='DOP853', rtol=1e-12, atol=1e-18, t_eval=t)
    return sol.t, sol.y

def charts(mu, k, U, p, y, eta, zeta):
    """each chart: (xi1, xi2, [dxi1/dU, dxi1/dp, dxi1/dy], [dxi2/dU, dxi2/dp, dxi2/dy])"""
    x = k*np.exp(U); yh = y-k; one = np.ones_like(U); zero = np.zeros_like(U)
    G = np.exp(U)-1-U
    w = np.sign(U)*np.sqrt(2*G); wU = np.where(np.abs(U) > 1e-8, (np.exp(U)-1)/np.where(w == 0, 1, w), 1.0)
    s = p-yh+mu*U; m = p+mu*U; xh = x-k
    C = {
        'C1 (U,p)': (U, p, [one, zero, zero], [zero, one, zero]),
        'C2 (U,y-k)': (U, yh, [one, zero, zero], [zero, zero, one]),
        'C3 (p,y-k)': (p, yh, [zero, one, zero], [zero, zero, one]),
        'C4 (y-k,m=p+mu U)': (yh, m, [zero, zero, one], [mu*one, one, zero]),
        'C5 (U,s)': (U, s, [one, zero, zero], [mu*one, one, -one]),
        'C6 (w,s)': (w, s, [wU, zero, zero], [mu*one, one, -one]),
        'C7 (w,p)': (w, p, [wU, zero, zero], [zero, one, zero]),
        'C8 (x-k,p)': (xh, p, [x, zero, zero], [zero, one, zero]),
        'C9 (x-k,y-k)': (xh, yh, [x, zero, zero], [zero, zero, one]),
        'C10 (w,y-k)': (w, yh, [wU, zero, zero], [zero, zero, one]),
        'C11 (1-e^-U,p)': (1-np.exp(-U), p, [np.exp(-U), zero, zero], [zero, one, zero]),
        'C12 (1-e^-U,s)': (1-np.exp(-U), s, [np.exp(-U), zero, zero], [mu*one, one, -one]),
    }
    return C

def run(mu, rho):
    k = rho*kH(mu); cyc = find_cycle(mu, k)
    res = dict(mu=mu, rho=rho, charts={})
    acc = {}
    for th in (0.0, 1.6, 3.2, 4.8):
        t, (U, p, y, eta, zeta) = orbit(mu, k, th)
        x = k*np.exp(U); zd = np.array([p, k-mu*p-y, x-y])
        r = np.sqrt(U**2+p**2)
        for name, (x1, x2, g1, g2) in charts(mu, k, U, p, y, eta, zeta).items():
            x1d = g1[0]*zd[0]+g1[1]*zd[1]+g1[2]*zd[2]; x2d = g2[0]*zd[0]+g2[1]*zd[1]+g2[2]*zd[2]
            wedge = x1d*x2-x2d*x1
            # Jacobian of chart w.r.t. (U,p) on the attractor (y=h(U,p), h_U=zeta, h_p=-eta)
            a11 = g1[0]+g1[2]*zeta; a12 = g1[1]-g1[2]*eta; a21 = g2[0]+g2[2]*zeta; a22 = g2[1]-g2[2]*eta
            jac = a11*a22-a12*a21
            with np.errstate(divide='ignore', invalid='ignore'):
                L = np.log(np.abs(jac))-np.log(np.abs(wedge))
            beta = (eta-mu)+np.gradient(L, t)
            good = (r > 0.02*cyc['U0']) & np.isfinite(beta)
            good[:5] = False; good[-5:] = False
            d = acc.setdefault(name, dict(beta=[], wedge=[], r=[]))
            d['beta'].append(beta[good]); d['wedge'].append(wedge[good]); d['r'].append(r[good])
    for name, d in acc.items():
        beta = np.concatenate(d['beta']); wedge = np.concatenate(d['wedge']); r = np.concatenate(d['r'])
        star = bool(np.all(wedge > 0) or np.all(wedge < 0))
        outer = r > 0.5*cyc['U0']
        res['charts'][name] = dict(star_shaped=star, frac_beta_pos=float(np.mean(beta > 0)),
                                   frac_beta_pos_outer=float(np.mean(beta[outer] > 0)),
                                   beta_q01=float(np.quantile(beta, 0.01)), beta_q99=float(np.quantile(beta, 0.99)))
    return res

if __name__ == '__main__':
    out = []
    for mu, rho in [(1.0, 1.1), (1.0, 2.0), (0.3, 1.5), (3.0, 1.3)]:
        r = run(mu, rho); out.append(r)
        print(f"mu={mu} rho={rho}")
        for name, v in r['charts'].items():
            print(f"   {name:20s} star={v['star_shaped']!s:5s} frac(beta>0)={v['frac_beta_pos']:.3f} outer={v['frac_beta_pos_outer']:.3f} q01={v['beta_q01']:+.3f} q99={v['beta_q99']:+.3f}")
    json.dump(out, open('job2_massera.json', 'w'), indent=1)
