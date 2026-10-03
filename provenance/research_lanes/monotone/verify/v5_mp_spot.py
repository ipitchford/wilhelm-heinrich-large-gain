#!/usr/bin/env python3
"""High-precision spot check of claim 5 where double precision cannot resolve the sign of DP:
Taylor-series integration (mpmath, 60 digits) of the state and the variational equation from a point of R to its first
crest return, then projection along the flow.  Independent of the compound-system argument.
usage: v5_mp_spot.py mu rho f theta [digits]    (U0 = f*U_crest(cycle), y0 = k + theta*(x0-k))"""
import sys, time
import numpy as np
import mpmath as mp
sys.path.insert(0, '/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03/work/lanes/monotone/verify')
from vlib import kH, find_cycle
from v2_retmap import compound_return
DPS = int(sys.argv[5]) if len(sys.argv) > 5 else 60
mp.mp.dps = DPS
N = max(30, int(0.72 * DPS))
TOL = mp.mpf(10) ** (-(DPS - 16))


def taylor_step(mu, k, st):
    """st = [U,x,p,y, dU1,dp1,dy1, dU2,dp2,dy2]; returns list of coefficient lists (order N)"""
    U = [st[0]]; x = [st[1]]; p = [st[2]]; y = [st[3]]
    dU = [[st[4]], [st[7]]]; dp = [[st[5]], [st[8]]]; dy = [[st[6]], [st[9]]]
    for n in range(N):
        inv = mp.mpf(1) / (n + 1)
        U.append(p[n] * inv)
        x.append(mp.fsum(x[j] * p[n - j] for j in range(n + 1)) * inv)
        p.append(((k if n == 0 else 0) - mu * p[n] - y[n]) * inv)
        y.append((x[n] - y[n]) * inv)
        for c in range(2):
            conv = mp.fsum(x[j] * dU[c][n - j] for j in range(n + 1))
            dU[c].append(dp[c][n] * inv)
            ndp = (-mu * dp[c][n] - dy[c][n]) * inv
            ndy = (conv - dy[c][n]) * inv
            dp[c].append(ndp); dy[c].append(ndy)
    return [U, x, p, y, dU[0], dp[0], dy[0], dU[1], dp[1], dy[1]]


def horner(c, h):
    s = mp.mpf(0)
    for a in reversed(c):
        s = s * h + a
    return s


def step_size(co):
    h = mp.inf
    for i, c in enumerate(co):
        # x = k e^U must keep RELATIVE accuracy (it re-grows by hundreds of orders of magnitude at the trough exit)
        sc = abs(c[0]) if i == 1 else max(abs(c[0]), mp.mpf(1))
        for m in (N, N - 1):
            if c[m] != 0:
                h = min(h, (TOL * sc / abs(c[m])) ** (mp.mpf(1) / m))
    return h * mp.mpf('0.8')


def mp_return(mu, k, U0, y0):
    mu = mp.mpf(mu); k = mp.mpf(k)
    st = [mp.mpf(U0), k * mp.exp(mp.mpf(U0)), mp.mpf(0), mp.mpf(y0), 1, 0, 0, 0, 0, 1]
    st = [mp.mpf(v) for v in st]
    t = mp.mpf(0); phase = 0; nsteps = 0      # phase 0: p<0 after crest, 1: p>0 after trough
    while True:
        co = taylor_step(mu, k, st)
        h = step_size(co)
        new = [horner(c, h) for c in co]
        nsteps += 1
        if phase == 0 and new[2] > 0:
            phase = 1
        elif phase == 1 and new[2] <= 0:
            # root of p in (0,h]
            lo, hi = mp.mpf(0), h
            for _ in range(220):
                mid = (lo + hi) / 2
                if horner(co[2], mid) > 0:
                    lo = mid
                else:
                    hi = mid
            s = (lo + hi) / 2
            fin = [horner(c, s) for c in co]
            return fin, t + s, nsteps
        st = new; t += h


if __name__ == '__main__':
    mu, rho, f_, th = [float(a) for a in sys.argv[1:5]]
    k = rho * kH(mu)
    cyc = find_cycle(mu, k)
    U0 = f_ * cyc['U0']; x0 = k * np.exp(U0); y0 = k + th * (x0 - k)
    U1, y1, tau, lnM, _ = compound_return(mu, k, U0, y0)
    t0 = time.time()
    fin, tmp, nsteps = mp_return(mu, k, U0, y0)
    Um, xm, pm, ym = fin[:4]
    fz = [pm, mp.mpf(k) - mp.mpf(mu) * pm - ym, xm - ym]
    M = mp.matrix(2, 2)
    for c in range(2):
        d = fin[4 + 3 * c: 7 + 3 * c]
        corr = d[1] / fz[1]
        M[0, c] = d[0] - fz[0] * corr
        M[1, c] = d[2] - fz[2] * corr
    print(f"mu={mu} rho={rho} point of R: U0={U0:.6f} (={f_} x crest), y0={y0:.6f}, x0={x0:.6f}")
    print(f"  double/log integrator : U1={U1:.12f} y1={y1:.12f} tau={tau:.10f}")
    print(f"  mp Taylor ({mp.mp.dps} digits, order {N}, {nsteps} steps, {time.time() - t0:.0f}s): U1={mp.nstr(Um, 14)} y1={mp.nstr(ym, 14)} tau={mp.nstr(tmp, 13)}  p at event={mp.nstr(pm, 3)}")
    print(f"  consistency x/(k e^U) - 1 at the event = {mp.nstr(xm / (mp.mpf(k) * mp.exp(Um)) - 1, 3)}")
    print(f"  in int R:  x1-y1={mp.nstr(xm - ym, 6)} > 0,  y1-k={mp.nstr(ym - mp.mpf(k), 6)} > 0")
    for i in range(2):
        for j in range(2):
            print(f"  DP[{i},{j}]  mp variational+projection = {mp.nstr(M[i, j], 12)}   compound(log) = {np.exp(lnM[i, j]):.10e}   (log10 {lnM[i, j] / np.log(10):.2f})")
    print(f"  all four entries positive in {DPS}-digit arithmetic:", all(M[i, j] > 0 for i in range(2) for j in range(2)))
