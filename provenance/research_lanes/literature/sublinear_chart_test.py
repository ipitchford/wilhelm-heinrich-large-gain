#!/usr/bin/env python3
"""Lane 'literature': necessary-condition test for applying the Krasnosel'skii / H.L. Smith (1986) monotone-concave
fixed-point theorem to the crest return map P on {p=0, p'<0}.

If, in a chart xi = (phi1(U), phi2(y)) centred at E (phi_i(E)=0, phi_i' > 0), P were monotone and concave in Smith's
sense (DP decreasing along ordered pairs), then at every fixed point q:  DP(q) q <= q componentwise.
In original section coordinates (U,y) this reads  M v <= v  with  v_i = phi_i/phi_i'  evaluated at the crest.
Power charts:  phi1 = (e^{aU}-1)/a,  phi2 = ((y/k)^b - 1)/b   =>   v = ((1-e^{-aU})/a, (y/b)(1-(y/k)^{-b})).
(a,b)=(1/2,1/2) is the manuscript's square-root chart (Thm 8.1, fails near Hopf); a->0 gives U; b->0 gives y log(y/k).
We report  d_i = (Mv)_i/v_i - 1  (need both <= 0).  Floating point only.
"""
import sys, numpy as np
sys.path.insert(0, '/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03/work')
from whlib import find_cycle, kH

def vchart(a, b, U, y, k):
    v1 = U if abs(a) < 1e-12 else (1 - np.exp(-a*U))/a
    v2 = y*np.log(y/k) if abs(b) < 1e-12 else (y/b)*(1 - (y/k)**(-b))
    return np.array([v1, v2])

mus = [0.2, 1.0, 5.0]
rhos = [1.02, 1.05, 1.2, 1.5, 2.0, 3.0, 5.0, 8.0]
cyc = {}
for mu in mus:
    for rho in rhos:
        k = rho*kH(mu)
        try:
            c = find_cycle(mu, k)
            if not (np.isfinite(c['m_t']) and c['U0'] > 1e-6):
                raise RuntimeError('degenerate cycle')
            cyc[(mu, rho)] = (k, c)
        except Exception as e:
            print('find_cycle failed mu=%g rho=%g: %s' % (mu, rho, repr(e)[:80]))
print('cycles found:', len(cyc))
grid = [0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0]
print('\nPer-cycle defect d = max_i[(Mv)_i/v_i - 1] for named charts (negative = necessary condition satisfied):')
print('  mu    rho    m_t        A:(U,y-k)   B:(x,y)     C:sqrt      D:(U,log y)')
for (mu, rho), (k, c) in sorted(cyc.items()):
    M = np.array(c['M']); U, y = c['U0'], c['y0']
    out = []
    for (a, b) in [(0, 1), (1, 1), (0.5, 0.5), (0, 0)]:
        v = vchart(a, b, U, y, k); d = (M @ v)/v - 1
        out.append(d.max())
    print('  %-5g %-6g %-10.6f %-11.3e %-11.3e %-11.3e %-11.3e' % (mu, rho, c['m_t'], *out))
print('\nGrid search over power charts (a,b): worst defect over all cycles at each mu, and best chart')
for mu in mus:
    best = None
    for a in grid:
        for b in grid:
            worst = -np.inf
            for (m, rho), (k, c) in cyc.items():
                if m != mu: continue
                M = np.array(c['M']); v = vchart(a, b, c['U0'], c['y0'], k)
                worst = max(worst, ((M @ v)/v - 1).max())
            if best is None or worst < best[0]:
                best = (worst, a, b)
    print('  mu=%g: best (a,b)=(%g,%g), worst defect over rho grid = %.3e  -> %s' % (mu, best[1], best[2], best[0], 'PASS on grid' if best[0] <= 0 else 'no power chart on grid passes'))
# universal chart
best = None
for a in grid:
    for b in grid:
        worst = max(((np.array(c['M']) @ vchart(a, b, c['U0'], c['y0'], k))/vchart(a, b, c['U0'], c['y0'], k) - 1).max() for (m, r), (k, c) in cyc.items())
        if best is None or worst < best[0]: best = (worst, a, b)
print('  all mu: best (a,b)=(%g,%g), worst defect = %.3e' % (best[1], best[2], best[0]))
