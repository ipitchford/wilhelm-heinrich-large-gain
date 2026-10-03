#!/usr/bin/env python3
"""Sampling of the attractor disc (W^u(E) inside the cycle, post-transient orbits outside) with pulled-back
eta, zeta, K, W and the full Hessian of h (grad h = (zeta, -eta)) obtained from one transversal variation.

Turn = crest (p=0, p'<0) to next crest.  Per-turn and per-point arrays are returned."""
import sys
sys.path.insert(0, '/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03/work/lanes/cartography')
import numpy as np
import cart as C


def _flow5(mu, k, s):
    U, p, y, eta, zeta = s[0], s[1], s[2], s[3], s[4]
    x = k * np.exp(U)
    return np.array([p, k - mu * p - y, x - y, zeta + (mu - 1) * eta - eta * eta, x - (1 + eta) * zeta])


def _hess(mu, k, R):
    """R: array (n,13) of mode-3 states.  Returns hUU, hUp, hpp, asym, detR."""
    U, p, y, eta, zeta = R[:, 0], R[:, 1], R[:, 2], R[:, 3], R[:, 4]
    dU, dp, de, dz = R[:, 7], R[:, 8], R[:, 10], R[:, 11]
    x = k * np.exp(U)
    pp = k - mu * p - y
    ed = zeta + (mu - 1) * eta - eta * eta
    zd = x - (1 + eta) * zeta
    det = p * dp - pp * dU
    with np.errstate(divide='ignore', invalid='ignore'):
        H11 = (zd * dp - dz * pp) / det
        H12 = (-zd * dU + dz * p) / det
        H21 = (-ed * dp + de * pp) / det
        H22 = (ed * dU - de * p) / det
        asym = np.abs(H12 - H21) / (np.abs(H11) + np.abs(H22) + np.abs(H12) + 1e-300)
    return H11, 0.5 * (H12 + H21), H22, asym, det


def sample_orbit(mu, k, s0, n_skip, n_rec, T_est, Ustar=None, pts_per_turn=400, max_turns=100000,
                 stop_tol=1e-9, rt=C.RT, maxsteps=20_000_000):
    """Integrate a mode-3 orbit from s0 (13-vector).  Skips n_skip crest-to-crest turns, then records up to n_rec
    turns.  Returns (turns, points): list of per-turn dicts and array of recorded points (with turn index)."""
    par = np.array([mu, k])
    rtol, atol = C.tolerances(3, mu, k, rt)
    s = s0.copy()
    ev_t = np.zeros(1)
    ev_y = np.zeros((1, 13))
    nrec_max = 60 * pts_per_turn + 20000
    rec_t = np.zeros(nrec_max)
    rec_y = np.zeros((nrec_max, 13))
    turns = []
    pts = []
    # go to the first crest
    ncrest = 0
    cur = None          # current turn record
    seg_pts = []
    t_abs = 0.0
    hmax = np.inf
    T_last = T_est
    for it in range(4 * max_turns):
        recording = ncrest > n_skip and len(turns) < n_rec
        hmax = T_last / pts_per_turn if recording else np.inf
        s[12] = s[12]
        st, t, ns, nev, nrec = C.integrate(3, par, s, 0.0, 1e3 * T_est + 1e7, rtol, atol, 1, 3, 1, ev_t, ev_y,
                                           rec_t if recording else C._EMPTY_T, rec_y if recording else C._EMPTY_Y,
                                           hmax, maxsteps)
        if st != 1:
            break
        t_abs += t
        pprime = k - mu * s[1] - s[2]
        is_crest = pprime < 0
        if recording and cur is not None:
            seg = rec_y[:nrec].copy()
            tt = rec_t[:nrec] + cur['_t0']
            seg_pts.append((tt, seg))
        if cur is not None:
            cur['_t0'] = cur['_t0'] + t
        # section-projected variation at this event
        f5 = _flow5(mu, k, s)
        d5 = np.array([s[7], s[8], s[9], s[10], s[11]])
        c = d5[1] / f5[1]
        dsec = d5 - c * f5          # dp = 0
        scale = dsec[0]
        if is_crest:
            # close the current turn
            if cur is not None:
                cur.update(U1=s[0], y1=s[2], eta1=s[3], zeta1=s[4], K1=s[5], W1=s[6], T=cur['_t0'],
                           I_eta=s[12] - cur['_I0'], zU1=dsec[4] / scale, eU1=dsec[3] / scale, yU1=dsec[2] / scale)
                if recording:
                    turns.append(cur)
                    if seg_pts:
                        tt = np.concatenate([a for a, b in seg_pts])
                        seg = np.concatenate([b for a, b in seg_pts])
                        pts.append((len(turns) - 1, tt, seg))
                T_last = cur['_t0']
            ncrest += 1
            # renormalise the variation inside the section: dU = 1, dp = 0
            s[7:12] = dsec / scale
            s[12] = 0.0
            cur = dict(U0=s[0], y0=s[2], eta0=s[3], zeta0=s[4], K0=s[5], W0=s[6], zU0=s[11], eU0=s[10], yU0=s[9],
                       _t0=0.0, _I0=0.0)
            seg_pts = []
            if len(turns) >= n_rec:
                break
            if Ustar is not None and len(turns) >= 1 and abs(s[0] - Ustar) < stop_tol * max(1.0, Ustar):
                break
        else:
            if cur is not None:
                cur.update(Ut=s[0], yt=s[2], etat=s[3], zetat=s[4], Kt=s[5], Wt=s[6], zUt=dsec[4] / scale,
                           eUt=dsec[3] / scale, t_trough=cur['_t0'])
        if ncrest > max_turns:
            break
    return turns, pts


def sample_attractor(mu, rho, n_theta=12, r0_list=None, n_out=24, out_span=None, n_rec_in=None, pts_per_turn=400,
                     verbose=False):
    """Returns dict with the primary cycle, per-turn table and per-point table for the attractor disc
    (inside: W^u(E)) and for post-transient outside orbits."""
    from hunt import primary_cycle
    k = rho * C.kH(mu)
    A = 1 + mu
    cyc, mono, hb, rb = primary_cycle(mu, rho)
    Us, ys, T = cyc['U0'], cyc['y0'], cyc['T']
    pr = C.profile(mu, k, Us, ys, T, record=False)
    ed = C.equilibrium_data(mu, k)
    turn_E = 2 * np.pi / ed['omega']
    growth = ed['alpha'] * turn_E                     # log expansion per turn near E
    contr = (ed['gammaE'] + ed['alpha']) * turn_E      # pull-back contraction per turn near E
    n_skip = int(np.ceil(27.0 / contr)) + 1
    if r0_list is None:
        if growth > 0.35:
            r0_list = [1e-6]
            n_theta_use = max(n_theta, 16)
        else:
            # slow spiral: seed the linear eigenplane at many amplitudes (O(r0^2) offset decays in the skipped turns)
            r0_list = list(Us * np.geomspace(2e-3, 0.6, 22)) + list(Us * np.linspace(0.65, 1.0, 12))
            n_theta_use = 1
    else:
        n_theta_use = n_theta
    if n_rec_in is None:
        n_rec_in = 400 if growth > 0.35 else 3
    turns, pts, src = [], [], []
    for r0 in r0_list:
        for ith in range(n_theta_use):
            th = 2 * np.pi * ith / n_theta_use
            z, _ = C.unstable_seed(mu, k, r0, th)
            zdot, _ = C.unstable_seed(mu, k, r0, th + np.pi / 2)
            s0 = np.zeros(13)
            s0[:3] = z
            s0[3], s0[4] = ed['etaE'], ed['zetaE']
            s0[5] = 1e-12 * k
            s0[6] = 1e-12 * k
            s0[7:10] = (zdot - np.array([0.0, 0.0, k])) / r0
            tr, pt = sample_orbit(mu, k, s0, n_skip, n_rec_in, turn_E, Ustar=Us, pts_per_turn=pts_per_turn)
            base = len(turns)
            for d in tr:
                d['region'] = 0
                d['orbit'] = len(src)
            turns += tr
            pts += [(base + i, tt, seg) for i, tt, seg in pt]
            src.append(('Wu', r0, th))
    # outside seeds on the crest section
    if out_span is None:
        out_span = min(max(1.0 * Us, 1.0), 3.0)
    T_hopf = 2 * np.pi / np.sqrt(mu)
    n_skip_out = int(np.ceil(32.0 / (A * T_hopf))) + 1
    for U0 in Us + np.concatenate((Us * np.geomspace(2e-3, 0.08, n_out // 3), np.linspace(0.1 * Us, out_span, n_out - n_out // 3))):
        s0 = np.zeros(13)
        s0[0] = U0
        s0[2] = k + k * np.expm1(U0) / (1 + pr['eta0'])
        s0[3], s0[4] = pr['eta0'], pr['zeta0']
        s0[5], s0[6] = pr['K0'], pr['W0']
        s0[7], s0[9] = 1.0, pr['zeta0']
        tr, pt = sample_orbit(mu, k, s0, n_skip_out, 2, T, Ustar=None, pts_per_turn=pts_per_turn)
        base = len(turns)
        for d in tr:
            d['region'] = 1
            d['orbit'] = len(src)
        turns += tr
        pts += [(base + i, tt, seg) for i, tt, seg in pt]
        src.append(('out', U0, 0.0))
    keys = ['U0', 'y0', 'eta0', 'zeta0', 'K0', 'W0', 'zU0', 'eU0', 'yU0', 'U1', 'y1', 'eta1', 'zeta1', 'K1', 'W1',
            'zU1', 'eU1', 'yU1', 'T', 'I_eta', 'Ut', 'yt', 'etat', 'zetat', 'Kt', 'Wt', 'zUt', 'eUt', 't_trough',
            'region', 'orbit']
    turns = [t for t in turns if all(kk in t for kk in keys)]
    TT = {kk: np.array([t[kk] for t in turns], dtype=float) for kk in keys}
    if pts:
        tid = np.concatenate([np.full(len(tt), i) for i, tt, seg in pts])
        tt = np.concatenate([tt for i, tt, seg in pts])
        seg = np.concatenate([seg for i, tt, seg in pts])
    else:
        tid = np.zeros(0)
        tt = np.zeros(0)
        seg = np.zeros((0, 13))
    hUU, hUp, hpp, asym, detR = _hess(mu, k, seg)
    P = dict(turn=tid, t=tt, U=seg[:, 0], p=seg[:, 1], y=seg[:, 2], eta=seg[:, 3], zeta=seg[:, 4], K=seg[:, 5],
             W=seg[:, 6], hUU=hUU, hUp=hUp, hpp=hpp, asym=asym)
    return dict(mu=mu, rho=rho, k=k, cycle=cyc, profile=pr, eq=ed, turns=TT, points=P, growth=growth, n_skip=n_skip,
                branch_monotone=mono)


if __name__ == '__main__':
    import time
    for mu, rho in [(1.0, 2.0), (0.1, 1.05), (0.0027, 5.02)]:
        t0 = time.time()
        S = sample_attractor(mu, rho)
        TT, P = S['turns'], S['points']
        print(f"mu={mu} rho={rho}: U*={S['cycle']['U0']:.4f} growth/turn={S['growth']:.3f} n_skip={S['n_skip']} turns={len(TT['U0'])} "
              f"(in {int((TT['region']==0).sum())}, out {int((TT['region']==1).sum())}) points={len(P['U'])} time={time.time()-t0:.1f}s")
        ins = TT['region'] == 0
        print('   U0/U* inside range', TT['U0'][ins].min() / S['cycle']['U0'], TT['U0'][ins].max() / S['cycle']['U0'])
        print('   consistency |yU0 - zeta0|/zeta0 max', np.max(np.abs(TT['yU0'] - TT['zeta0']) / TT['zeta0']),
              ' Hessian asym: median', np.nanmedian(P['asym']), 'max', np.nanmax(P['asym']))
        # divergence formula vs direct finite-difference of the return map
        phi_p = (TT['y0'] - S['k']) / (TT['y1'] - S['k']) * np.exp(TT['I_eta'] - mu * TT['T'])
        o = np.argsort(TT['U0'][ins])
        U = TT['U0'][ins][o]; U1 = TT['U1'][ins][o]; ph = phi_p[ins][o]
        fd = np.gradient(U1, U)
        sel = slice(2, -2)
        print('   phi\' (divergence formula) vs finite difference: max rel diff', np.max(np.abs(fd[sel] / ph[sel] - 1)))
