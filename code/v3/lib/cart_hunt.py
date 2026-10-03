#!/usr/bin/env python3
# --- v3.0 bundle copy of research file `lanes/cartography/hunt.py`.
# --- Changes from the original: file renamed; absolute import path removed; imports cart_grid (the bundle copy of grid.py).
"""Part 2: hunt for non-uniqueness at random parameter points.

For each (mu, rho):
  0. primary cycle by amplitude continuation from Hopf (records whether rho(h) is increasing along the branch);
  A. displacement scan of the crest return map on the slow manifold: seeds on the crest section, transient
     returns discarded (normal contraction e^{-AT} per return), then (U_n, U_{n+1}) pairs recorded, inside and
     outside the primary crest, each seed launched from two different off-manifold y0 values;
     pass  <=>  sign(U_{n+1}-U_n) = sign(U* - U_n) everywhere and r(U) = log(phi(U)/U) strictly decreasing;
  B. long runs from W^u(E), far-outside crests and random 3-D initial conditions: all must converge to the
     primary crest (Newton polish at the end when contraction is slow);
  C. Newton on the crest map from perturbed guesses: any fixed point distinct from the primary one is reported.
"""
import sys, time, json
import numpy as np
import cart as C
from cart_grid import solve_h


def primary_cycle(mu, rho_t):
    """Amplitude continuation to rho_t.  Returns (cycle dict at k_t, branch monotone flag, branch arrays)."""
    hs = list(np.geomspace(0.01, 1.5, 22)) + list(np.arange(1.75, 60.01, 0.25))
    branch = []
    prev = prev2 = None
    for h in hs:
        c = solve_h(mu, h, prev, prev2)
        branch.append(c)
        prev2, prev = prev, c
        if c['rho'] > rho_t:
            break
    rb = np.array([b['rho'] for b in branch])
    hb = np.array([b['U0'] for b in branch])
    if rb[-1] < rho_t:
        raise RuntimeError('branch did not reach target rho')
    if rb[0] > rho_t:
        # need smaller amplitudes
        hs2 = list(np.geomspace(1e-4, 0.01, 12))[:-1]
        br2 = []
        prev = prev2 = None
        for h in hs2:
            c = solve_h(mu, h, prev, prev2)
            br2.append(c)
            prev2, prev = prev, c
        branch = br2 + branch
        rb = np.array([b['rho'] for b in branch])
        hb = np.array([b['U0'] for b in branch])
    mono = bool(np.all(np.diff(rb) > 0))
    i = int(np.searchsorted(rb, rho_t))
    i = max(i, 1)
    ca, cb = branch[i - 1], branch[i]
    ha, hb_ = ca['U0'], cb['U0']
    ga, gb = np.log(ca['rho'] - 1) - np.log(rho_t - 1), np.log(cb['rho'] - 1) - np.log(rho_t - 1)
    pa, pb = ca, cb
    c = cb
    for it in range(40):
        if gb == ga:
            break
        hn = hb_ - gb * (hb_ - ha) / (gb - ga)
        c = solve_h(mu, hn, pb, pa)
        gn = np.log(c['rho'] - 1) - np.log(rho_t - 1)
        ha, ga, pa = hb_, gb, pb
        hb_, gb, pb = hn, gn, c
        if abs(c['rho'] - rho_t) < 2e-13 * rho_t or abs(hb_ - ha) < 1e-14 * hb_:
            break
    k_t = rho_t * C.kH(mu)
    cyc = dict(c)
    cyc['k'] = k_t
    if c['m_t'] < 0.9:
        try:
            cf = C.newton_fixed_k(mu, k_t, c['U0'], c['y0'], maxit=6)
            if cf['residual'] < 1e-10:
                cyc = cf
        except Exception:
            pass
    return cyc, mono, hb, rb


def hunt_point(mu, rho, rng, budget_s=25.0, verbose=False):
    t_start = time.time()
    out = dict(mu=mu, rho=rho)
    k = rho * C.kH(mu)
    A = 1 + mu
    cyc, mono, hb, rb = primary_cycle(mu, rho)
    Us, ys, T = cyc['U0'], cyc['y0'], cyc['T']
    pr = C.profile(mu, k, Us, ys, T, record=False)
    out.update(k=k, U_star=Us, y_star_over_k=ys / k, T=T, log_mt=pr['log_mt'], lambda_c=pr['lambda_c'],
               m_t=float(np.exp(pr['log_mt'])), branch_monotone=mono, primary_residual=cyc['residual'])
    eta_c = pr['eta0']
    # time one return
    t0 = time.time()
    C.crest_return(mu, k, Us, ys)
    t_ret = max(time.time() - t0, 2e-5)
    T_hopf = 2 * np.pi / np.sqrt(mu)
    n_trans = int(np.ceil(32.0 / (A * T_hopf))) + 1
    n_rec = 2

    def y_slow(U0, f=1.0):
        return k + f * k * np.expm1(U0) / (1.0 + eta_c)

    # ---------------- A. displacement scan (time-boxed; seeds visited coarse-to-fine)
    n_in, n_out = 70, 70
    span_out = min(max(1.5 * Us, 2.0), 4.0)
    seeds_in = Us * np.concatenate((np.geomspace(1e-4, 0.3, n_in // 3), np.linspace(0.3, 0.999, n_in - n_in // 3 + 1)[1:]))
    seeds_out = Us + np.concatenate((Us * np.geomspace(1e-3, 0.1, n_out // 3), np.linspace(0.1 * Us, span_out, n_out - n_out // 3 + 1)[1:]))
    seeds = np.concatenate((seeds_in, seeds_out))
    order = rng.permutation(len(seeds))
    pairs = []      # (U_n, U_{n+1}, y_n, source)   source 0/1: seed variants, 2: W^u(E)
    nfail = 0
    nused = 0
    t_box = time.time()
    for idx in order:
        if time.time() - t_box > 0.45 * budget_s and nused >= 12:
            break
        U0 = seeds[idx]
        nused += 1
        for var, f in enumerate((1.0, 0.6 + 0.8 * rng.random())):
            try:
                U, Y, tt, st = C.crest_iter(mu, k, U0, y_slow(U0, f), n_trans + n_rec, maxsteps=3_000_000)
            except Exception:
                nfail += 1
                continue
            if st != 1:
                nfail += 1
                continue
            for j in range(n_trans, n_trans + n_rec):
                pairs.append((U[j], U[j + 1], Y[j], var))
    # W^u(E) orbits: every crest lies on the attractor disc, no transient needed
    n_wu_pairs = 0
    t_box = time.time()
    for ith in range(8):
        if time.time() - t_box > 0.15 * budget_s and ith >= 2:
            break
        z, ed = C.unstable_seed(mu, k, 1e-4, theta=2 * np.pi * (ith + rng.random()) / 8)
        n_wu = 60
        par = np.array([mu, k])
        rtol, atol = C.tolerances(0, mu, k)
        ev_t = np.zeros(n_wu)
        ev_y = np.zeros((n_wu, 3))
        s = z.copy()
        st, t, ns, nev, nrec = C.integrate(0, par, s, 0.0, 60.0 * T + 200.0 / max(ed['alpha'], 1e-12) if ed['alpha'] * T > 0.02 else 60.0 * T,
                                           rtol, atol, 1, -2, n_wu, ev_t, ev_y, C._EMPTY_T, C._EMPTY_Y, np.inf, 3_000_000)
        for j in range(n_trans, nev - 1):
            if ev_y[j, 0] > 0:
                pairs.append((ev_y[j, 0], ev_y[j + 1, 0], ev_y[j, 2], 2))
                n_wu_pairs += 1
    P = np.array(pairs)
    out['scan_seeds'] = int(nused)
    out['scan_fail'] = int(nfail)
    out['scan_pairs'] = int(len(P))
    out['scan_wu_pairs'] = int(n_wu_pairs)
    out['_pairs'] = P
    if len(P):
        Un, Un1 = P[:, 0], P[:, 1]
        d = Un1 - Un
        tolU = 1e-8 * max(1.0, Us) + 50 * abs(pr['closure'])
        far = np.abs(Un - Us) > 1e-6 * Us + tolU
        wrong = far & (np.sign(d) != np.sign(Us - Un)) & (np.abs(d) > tolU)
        out['scan_wrong_sign'] = int(wrong.sum())
        out['scan_U_over_Ustar_range'] = [float(Un.min() / Us), float(Un.max() / Us)]
        # monotonicity of r(U) = log(phi(U)/U)
        o = np.argsort(Un)
        Uo, r = Un[o], np.log(Un1[o] / Un[o])
        dU = np.diff(Uo)
        dr = np.diff(r)
        ok = dU > 1e-7 * max(1.0, Us)
        noise = 1e-9 + 20 * abs(pr['closure']) + 2e-10 / Uo[:-1]
        viol = ok & (dr > noise)
        out['scan_r_max_increase_over_noise'] = float((dr[ok] / noise[ok]).max()) if ok.any() else float('nan')
        out['scan_r_monotone_violations'] = int(viol.sum())
        out['scan_r_max_increase'] = float(dr[ok].max()) if ok.any() else float('nan')
        out['scan_r_slope_max'] = float((dr[ok] / dU[ok]).max()) if ok.any() else float('nan')
        out['scan_r_at_E'] = float(2 * np.pi * ed['alpha'] / ed['omega'])
        if wrong.any():
            out['scan_wrong_examples'] = [[float(a_), float(b_)] for a_, b_ in zip(Un[wrong][:5], Un1[wrong][:5])]
    # ---------------- B. long runs (time-boxed)
    ics = []
    z, ed = C.unstable_seed(mu, k, 1e-3, theta=2 * np.pi * rng.random())
    ics.append(('Wu', z))
    lam = k / mu
    for j in range(4):
        U0 = rng.uniform(-6.0, Us + 2.0)
        p0 = lam * (1 - np.exp(rng.normal(0.0, 1.5)))            # any p < lambda
        y0 = k * np.exp(rng.normal(0.0, 2.0))
        ics.append((f'rand{j}', np.array([U0, p0, y0])))
    ics.append(('out+3', np.array([Us + 3.0, 0.0, y_slow(Us + 3.0, 0.5 + rng.random())])))
    ics.append(('out+6', np.array([Us + 6.0, 0.0, y_slow(Us + 6.0, 0.5 + rng.random())])))
    long_res = []
    t_box = time.time()
    per_ic = 0.3 * budget_s / len(ics)
    for name, z in ics:
        try:
            par = np.array([mu, k])
            s = z.copy()
            rtol, atol = C.tolerances(0, mu, k)
            chunk = 8
            ev_t = np.zeros(chunk)
            ev_y = np.zeros((chunk, 3))
            t_ic = time.time()
            tot = 0
            Uhist = []
            yf = None
            while True:
                st, t, ns, nev, nrec = C.integrate(0, par, s, 0.0, C.tmax_default(mu, k) * chunk, rtol, atol, 1, -2,
                                                   chunk, ev_t, ev_y, C._EMPTY_T, C._EMPTY_Y, np.inf, 4_000_000)
                tot += nev
                Uhist += list(ev_y[:nev, 0])
                if nev:
                    yf = ev_y[nev - 1, 2]
                if st != 1:
                    break
                if len(Uhist) >= 2 and abs(Uhist[-1] - Us) < 1e-9 * max(1.0, Us) and abs(Uhist[-2] - Us) < 1e-9 * max(1.0, Us):
                    break
                if time.time() - t_ic > per_ic or tot >= 4000:
                    break
            if tot < 2:
                long_res.append(dict(ic=name, status='too_expensive' if st == -1 else 'too_few_returns', nev=int(tot)))
                continue
            Uf = Uhist[-1]
            gap = abs(Uf - Us) / max(1.0, Us)
            rec = dict(ic=name, nev=int(tot), gap_after_iter=float(gap), min_U_seen=float(min(Uhist)), max_U_seen=float(max(Uhist)))
            if gap > 1e-7:
                # still moving: is it at least heading towards U*?
                rec['heading_ok'] = bool(np.sign(Uhist[-1] - Uhist[-2]) == np.sign(Us - Uhist[-1]))
                try:
                    cf = C.newton_fixed_k(mu, k, Uf, yf, maxit=60, maxsteps=3_000_000)
                    rec['newton_gap'] = float(abs(cf['U0'] - Us) / max(1.0, Us))
                    rec['newton_res'] = cf['residual']
                    if abs(cf['U0']) < 1e-5 * max(1.0, Us) or cf['U0'] < 0:
                        rec['status'] = 'slow_newton_to_E'
                    else:
                        rec['status'] = 'same' if (rec['newton_gap'] < 1e-6 and cf['residual'] < 1e-8) else 'DIFFERENT'
                    if rec['status'] == 'DIFFERENT':
                        rec['newton_U0'] = cf['U0']
                        rec['newton_y0'] = cf['y0']
                        rec['newton_mt'] = cf['m_t']
                except Exception as e:
                    rec['status'] = 'newton_failed'
                    rec['err'] = repr(e)
            else:
                rec['status'] = 'same'
            long_res.append(rec)
        except Exception as e:
            long_res.append(dict(ic=name, status='error', err=repr(e)))
    out['long'] = long_res
    out['long_same'] = int(sum(r['status'] == 'same' for r in long_res))
    out['long_different'] = int(sum(r['status'] == 'DIFFERENT' for r in long_res))
    out['long_other'] = int(sum(r['status'] not in ('same', 'DIFFERENT') for r in long_res))
    out['long_heading_bad'] = int(sum(r.get('heading_ok') is False for r in long_res))
    # ---------------- C. Newton from perturbed guesses
    newt = dict(same=0, equilibrium=0, failed=0, new=[])
    for j in range(12):
        U0 = min(Us * np.exp(rng.normal(0.0, 0.8)), Us + 3.0)
        y0 = k + (ys - k) * np.expm1(U0) / np.expm1(Us) * np.exp(rng.normal(0.0, 0.5))
        try:
            cf = C.newton_fixed_k(mu, k, U0, y0, maxit=60, maxsteps=3_000_000)
        except Exception:
            newt['failed'] += 1
            continue
        if cf['residual'] > 1e-8:
            newt['failed'] += 1
        elif abs(cf['U0'] - Us) < 1e-6 * max(1.0, Us):
            newt['same'] += 1
        elif abs(cf['U0']) < 1e-5 * max(1.0, Us) or cf['U0'] < 0:
            newt['equilibrium'] += 1
        else:
            newt['new'].append(dict(U0=cf['U0'], y0=cf['y0'], T=cf['T'], m_t=cf['m_t'], res=cf['residual']))
    out['newton'] = newt
    out['seconds'] = time.time() - t_start
    return out


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    seed = int(sys.argv[2]) if len(sys.argv) > 2 else 20261003
    outname = sys.argv[3] if len(sys.argv) > 3 else 'hunt_results.jsonl'
    budget = float(sys.argv[4]) if len(sys.argv) > 4 else 25.0
    rng = np.random.default_rng(seed)
    allpairs = {}
    with open(outname, 'a') as f:
        for i in range(n):
            mu = float(10 ** rng.uniform(-3, 3))
            rho = float(1 + 10 ** rng.uniform(-3, np.log10(299.0)))
            try:
                r = hunt_point(mu, rho, rng, budget_s=budget)
            except Exception as e:
                r = dict(mu=mu, rho=rho, error=repr(e))
            r['index'] = i
            r['seed'] = seed
            P = r.pop('_pairs', None)
            if P is not None:
                allpairs[f'p{i}'] = P
            f.write(json.dumps(r) + '\n')
            f.flush()
            if 'error' in r:
                print(f"[{i}] mu={mu:.4g} rho={rho:.5g} ERROR {r['error']}", flush=True)
            else:
                print(f"[{i}] mu={mu:.4g} rho-1={rho-1:.4g} U*={r['U_star']:.3f} log_mt={r['log_mt']:+.4e} mono={r['branch_monotone']}"
                      f" scan: pairs={r['scan_pairs']} wrong={r.get('scan_wrong_sign')} rviol={r.get('scan_r_monotone_violations')} wu={r.get('scan_wu_pairs')} rslope_max={r.get('scan_r_slope_max'):.3e}"
                      f" long: same={r['long_same']} diff={r['long_different']} other={r['long_other']}"
                      f" newton: same={r['newton']['same']} E={r['newton']['equilibrium']} fail={r['newton']['failed']} NEW={len(r['newton']['new'])} t={r['seconds']:.1f}s", flush=True)


    np.savez_compressed(outname.replace('.jsonl', '_pairs.npz'), **allpairs)


if __name__ == '__main__':
    main()
