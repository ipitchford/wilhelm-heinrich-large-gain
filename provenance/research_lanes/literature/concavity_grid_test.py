#!/usr/bin/env python3
"""Lane 'literature': direct grid test of hypotheses (M) and (C) of H.L. Smith's 1986 monotone-concave fixed point
theorem (as quoted by Farcot-Gouze, RR-6875, Theorem 1) for the crest return map P of the Wilhelm-Heinrich system.

Chart: xi1 = (e^{aU}-1)/a, xi2 = ((y/k)^b-1)/b, origin = equilibrium E.  DP_chart = diag(phi'(P)) M diag(1/phi'(start)).
(M): all entries of DP_chart > 0.   (C): entries of DP_chart non-increasing when xi increases (tested between grid
neighbours in the +U, +y and diagonal directions).  Floating point; a finite grid; not a certificate.
"""
import sys, numpy as np
sys.path.insert(0, '/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03/work')
from whlib import find_cycle, crest_return, kH

def run(mu, rho, charts, n=9, lo=0.1, hi=1.5):
    k = rho*kH(mu)
    c = find_cycle(mu, k)
    Us, ys = c['U0'], c['y0']
    fu = np.linspace(lo, hi, n); fy = np.linspace(lo, hi, n)
    data = {}
    for i, a_ in enumerate(fu):
        for j, b_ in enumerate(fy):
            U0 = a_*Us; y0 = k + b_*(ys - k)
            try:
                U1, y1, T, M = crest_return(mu, k, U0, y0)
                data[(i, j)] = (U0, y0, U1, y1, M)
            except Exception as e:
                data[(i, j)] = None
    ok = sum(v is not None for v in data.values())
    print('mu=%g rho=%g k=%.4f crest=(%.4f, %.4f) m_t=%.5f ; returns computed %d/%d' % (mu, rho, k, Us, ys, c['m_t'], ok, n*n))
    negM = sum(1 for v in data.values() if v is not None and (v[4] <= 0).any())
    outside = sum(1 for v in data.values() if v is not None and (v[2] <= 0 or v[3] <= k))
    print('   points with a non-positive entry of M (violates (M)): %d ; returns leaving the quadrant {U>0,y>k}: %d' % (negM, outside))
    for (a, b) in charts:
        def D(v):
            U0, y0, U1, y1, M = v
            L = np.diag([np.exp(a*U1), (y1/k)**(b-1)/k]); R = np.diag([np.exp(-a*U0), k*(y0/k)**(1-b)])
            return L @ M @ R
        viol = 0; tot = 0; worst = 0.0; wloc = None
        for (i, j), v in data.items():
            if v is None: continue
            for (di, dj) in [(1, 0), (0, 1), (1, 1)]:
                w = data.get((i+di, j+dj))
                if w is None: continue
                inc = D(w) - D(v)           # should be <= 0 entrywise
                rel = inc/np.abs(D(v))
                tot += 4; viol += int((inc > 1e-9*np.abs(D(v))).sum())
                if rel.max() > worst: worst = rel.max(); wloc = (fu[i], fy[j], di, dj, int(np.argmax(rel)))
        print('   chart (a,b)=(%.4g,%.4g): (C) violations %d / %d entry-comparisons ; worst relative increase %.3e at %s' % (a, b, viol, tot, worst, wloc))

if __name__ == '__main__':
    for mu, rho in [(1.0, 2.0), (1.0, 1.1), (0.2, 1.5)]:
        astar = (mu**2 + 4*mu + 1)/(mu**2 + 6*mu + 1)   # b=1 chart with zero O(h) Hopf mismatch
        run(mu, rho, [(1.0, 1.0), (astar, 1.0), (0.5, 0.5), (0.0, 1.0), (0.0, 0.0)])


def diag_only(mu, rho, charts, n=9, lo=0.1, hi=1.5):
    """Same test restricted to the diagonal of the grid (points (s U*, k + s (y*-k)), i.e. the ray through the crest),
    comparing consecutive points along the ray only."""
    k = rho*kH(mu); c = find_cycle(mu, k); Us, ys = c['U0'], c['y0']
    ss = np.linspace(lo, hi, 2*n)
    pts = []
    for s in ss:
        U0 = s*Us; y0 = k + s*(ys-k)
        U1, y1, T, M = crest_return(mu, k, U0, y0)
        pts.append((U0, y0, U1, y1, M))
    for (a, b) in charts:
        def D(v):
            U0, y0, U1, y1, M = v
            return np.diag([np.exp(a*U1), (y1/k)**(b-1)/k]) @ M @ np.diag([np.exp(-a*U0), k*(y0/k)**(1-b)])
        viol = 0; worst = -np.inf
        for v, w in zip(pts[:-1], pts[1:]):
            rel = (D(w) - D(v))/np.abs(D(v)); viol += int((rel > 1e-9).sum()); worst = max(worst, rel.max())
        print('   ray test mu=%g rho=%g chart (%.4g,%.4g): violations %d / %d ; worst relative increase %.3e' % (mu, rho, a, b, viol, 4*(len(pts)-1), worst))

if __name__ == '__main__':
    print('\nRay-only test (ordered pairs along the ray from E through the crest):')
    for mu, rho in [(1.0, 2.0), (1.0, 1.1), (0.2, 1.5)]:
        astar = (mu**2 + 4*mu + 1)/(mu**2 + 6*mu + 1)
        diag_only(mu, rho, [(1.0, 1.0), (astar, 1.0), (0.5, 0.5), (0.0, 1.0)])
