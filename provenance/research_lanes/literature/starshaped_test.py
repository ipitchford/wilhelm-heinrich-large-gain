#!/usr/bin/env python3
"""Lane 'literature': test of the hypothesis of Sabatini's uniqueness theorem (arXiv:1006.2020, Theorem 3;
Massera-type 'strict star-shapedness' after the Conti-Filippov transformation) for the reduced planar flow
    u' = v,  v' = k - mu v - h(u,v)          (brief item (b); on cycles this is Prop. 5.2 of the v2.0 manuscript)
written as  u'' + u' Phi(u,u') + g(u) = 0,   g(u) = h(u,0) - k,   Phi(u,v) = mu + (h(u,v) - h(u,0))/v.
Hypothesis to test:   S(u,v) := (sgn(u) sqrt(2G)/g) [ 2G (Phi_u g - Phi g')/g^2 + Phi ] + v Phi_v   has one sign.
h(u,0), zeta(u,0) are interpolated (PCHIP) from the p=0 crossings of an orbit on W^u(E) plus the cycle's crest/trough.
Also reports the raw-plane radial-angular-monotonicity quantity alpha = -v K (Sabatini-Villari Thm 6 needs alpha >= 0).
Floating point; interpolation error not controlled -> only gross sign changes are meaningful.
"""
import sys, numpy as np
from scipy.interpolate import PchipInterpolator
sys.path.insert(0, '/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03/work')
from whlib import find_cycle, cycle_profile, unstable_manifold_orbit, equilibrium_data, kH

def run(mu, rho):
    k = rho*kH(mu)
    c = find_cycle(mu, k); pr = cycle_profile(mu, k, c)
    orb = unstable_manifold_orbit(mu, k, n=200001)
    ed = equilibrium_data(mu, k)
    U, p, y, et, ze = orb['U'], orb['p'], orb['y'], orb['eta'], orb['zeta']
    idx = np.where(np.sign(p[:-1])*np.sign(p[1:]) < 0)[0]
    us, hs, zs = [0.0], [k], [ed['zetaE']]
    for i in idx:
        s = p[i]/(p[i]-p[i+1])
        uc = U[i] + s*(U[i+1]-U[i])
        if abs(uc) < 1e-5: continue
        us.append(uc); hs.append(y[i] + s*(y[i+1]-y[i])); zs.append(ze[i] + s*(ze[i+1]-ze[i]))
    # cycle crest / trough
    j = np.where(np.sign(pr['p'][:-1])*np.sign(pr['p'][1:]) < 0)[0]
    us += [pr['U'][0]] + [pr['U'][i] for i in j]; hs += [pr['y'][0]] + [pr['y'][i] for i in j]; zs += [pr['zeta'][0]] + [pr['zeta'][i] for i in j]
    us, hs, zs = map(np.array, (us, hs, zs))
    o = np.argsort(us); us, hs, zs = us[o], hs[o], zs[o]
    keep = np.concatenate(([True], np.diff(us) > 1e-7)); us, hs, zs = us[keep], hs[keep], zs[keep]
    h0 = PchipInterpolator(us, hs); z0 = PchipInterpolator(us, zs)
    G0 = h0.antiderivative()
    def S_of(Uv, pv, yv, etav, zev):
        m = (np.abs(pv) > 0.05*np.abs(pv).max()) & (np.abs(Uv) > 0.05*np.abs(Uv).max()) & (Uv > us[0]) & (Uv < us[-1])
        u, v, yy, e, z = Uv[m], pv[m], yv[m], etav[m], zev[m]
        g = h0(u) - k; gp = z0(u); G = (G0(u) - G0(0.0)) - k*u
        Phi = mu + (yy - h0(u))/v; Phiu = (z - z0(u))/v; Phiv = (-e*v - (yy - h0(u)))/v**2
        S = (np.sign(u)*np.sqrt(2*np.abs(G))/g)*(2*G*(Phiu*g - Phi*gp)/g**2 + Phi) + v*Phiv
        return u, v, S, G, g
    u, v, S, G, g = S_of(pr['U'], pr['p'], pr['y'], pr['eta'], pr['zeta'])
    print('mu=%g rho=%g: %d p=0 crossings used for h(u,0); min G=%.2e (should be >0), u*g>0 everywhere: %s' % (mu, rho, len(us), G.min(), bool((u*g > 0).all())))
    print('   on the cycle: S>0 at %d pts, S<0 at %d pts ; S range [%.3e, %.3e]' % ((S > 0).sum(), (S < 0).sum(), S.min(), S.max()))
    for q, name in [((u > 0) & (v > 0), 'u>0,v>0'), ((u > 0) & (v < 0), 'u>0,v<0'), ((u < 0) & (v < 0), 'u<0,v<0'), ((u < 0) & (v > 0), 'u<0,v>0')]:
        if q.any(): print('      quadrant %-8s: S in [%.3e, %.3e]' % (name, S[q].min(), S[q].max()))
    alpha = -pr['p']*pr['K']
    print('   raw-plane RAM quantity alpha = -v K on the cycle: min %.3e max %.3e (sign changes with v: %s)' % (alpha.min(), alpha.max(), bool(alpha.min() < 0 < alpha.max())))

if __name__ == '__main__':
    for mu, rho in [(1.0, 1.1), (1.0, 1.3), (0.3, 1.1)]:
        try:
            run(mu, rho)
        except Exception as e:
            print('mu=%g rho=%g FAILED: %r' % (mu, rho, e))
