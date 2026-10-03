#!/usr/bin/env python3
"""Part 1: weak-exponent map.  Amplitude continuation (crest U_max = h fixed, unknowns y0,k) from Hopf,
then secant on h to hit each target rho.  Writes branch_table.csv, grid_table.csv, grid_cycles.npz."""
import sys, time, csv, json
sys.path.insert(0, '/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03/work/lanes/cartography')
import numpy as np
import cart as C

MUS = [0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.85, 1.0]
RHOS = 1.0 + np.logspace(-3, np.log10(199.0), 25)
RHO_MAX = 200.0


def predict(prev, prev2, mu, h):
    if prev is None:
        return C.hopf_guess(mu, h)
    if prev2 is None:
        s = (h / prev['U0']) ** 2
        k = C.kH(mu) * (1 + (prev['rho'] - 1) * s)
        y0 = k * (1 + (prev['y0'] / prev['k'] - 1) * (np.expm1(h) / np.expm1(prev['U0'])))
        return y0, k
    if h < 1.5:
        a = (np.log(h) - np.log(prev2['U0'])) / (np.log(prev['U0']) - np.log(prev2['U0']))
    else:
        a = (h - prev2['U0']) / (prev['U0'] - prev2['U0'])
    lr = np.log(prev2['rho'] - 1) + a * (np.log(prev['rho'] - 1) - np.log(prev2['rho'] - 1))
    ly = np.log(prev2['y0'] / prev2['k'] - 1) + a * (np.log(prev['y0'] / prev['k'] - 1) - np.log(prev2['y0'] / prev2['k'] - 1))
    k = C.kH(mu) * (1 + np.exp(lr))
    return k * (1 + np.exp(ly)), k


def solve_h(mu, h, prev, prev2):
    y0, k = predict(prev, prev2, mu, h)
    return C.newton_fixed_h(mu, h, y0, k, maxit=80)


def describe(mu, c):
    """Profile + derived scalars for a cycle dict (from newton_fixed_h / newton_fixed_k)."""
    k = c['k']
    pr = C.profile(mu, k, c['U0'], c['y0'], c['T'], record=False)
    A = 1 + mu
    kh = C.kH(mu)
    q1 = pr['lambda_c']
    ident = ((A * A + mu) * q1 + 2 * A * pr['q2'] + pr['q3']) - (k - kh)
    lm_M = np.log(c['m_t']) if c['m_t'] > 1e-7 else np.nan
    row = dict(mu=mu, rho=k / kh, k=k, U_max=c['U0'], U_min=pr['U_min'], y0_over_k=c['y0'] / k, T=c['T'],
               T_scaled=c['T'] * np.sqrt(mu) / (2 * np.pi), log_mt=pr['log_mt'], log_mt_M=lm_M,
               lambda_c=pr['lambda_c'], lambda_c_T=pr['log_mt'], mean_eta=pr['mean_eta'],
               meanU_plus_log_rho=pr['mean_U'] + np.log(k / kh), mean_U=pr['mean_U'],
               q2=pr['q2'], q3=pr['q3'], moment_identity_residual=ident / (k - kh),
               mean_x_rel_err=pr['mean_x'] / k - 1, newton_residual=c['residual'], closure=pr['closure'],
               ratio_defect=pr['ratio_defect'], settle_periods=pr['settle_periods'],
               eta_crest=pr['eta0'], zeta_crest=pr['zeta0'])
    return row


def run_mu(mu, log):
    # --- coarse branch in amplitude
    hs = list(np.geomspace(0.01, 1.5, 22)) + list(np.arange(1.75, 40.01, 0.25))
    branch = []
    prev = prev2 = None
    for h in hs:
        c = solve_h(mu, h, prev, prev2)
        if c['residual'] > 1e-9:
            log(f"  WARNING branch residual {c['residual']:.1e} at mu={mu} h={h}")
        branch.append(c)
        prev2, prev = prev, c
        if c['rho'] > 1.25 * RHO_MAX:
            break
    rb = np.array([b['rho'] for b in branch])
    hb = np.array([b['U0'] for b in branch])
    mono = bool(np.all(np.diff(rb) > 0))
    log(f"mu={mu}: branch {len(branch)} pts, h in [{hb[0]:.3g},{hb[-1]:.3g}], rho-1 in [{rb[0]-1:.3e},{rb[-1]-1:.3e}], rho(h) increasing: {mono}")
    brows = []
    for b in branch:
        r = describe(mu, b)
        brows.append(r)
    # --- targets
    rows = []
    for rho_t in RHOS:
        if rho_t > rb[-1] or rho_t < rb[0]:
            log(f"  target rho={rho_t} outside branch range; skipped")
            continue
        i = int(np.searchsorted(rb, rho_t))
        ca, cb = branch[i - 1], branch[i]
        # secant in h on g(h) = log(rho(h)-1) - log(rho_t-1), seeded from the bracket
        ha, hb_ = ca['U0'], cb['U0']
        ga, gb = np.log(ca['rho'] - 1) - np.log(rho_t - 1), np.log(cb['rho'] - 1) - np.log(rho_t - 1)
        pa, pb = ca, cb
        c = None
        for it in range(40):
            hn = hb_ - gb * (hb_ - ha) / (gb - ga)
            c = solve_h(mu, hn, pb, pa)
            gn = np.log(c['rho'] - 1) - np.log(rho_t - 1)
            ha, ga, pa = hb_, gb, pb
            hb_, gb, pb = hn, gn, c
            if abs(c['rho'] - rho_t) < 2e-13 * rho_t or abs(hb_ - ha) < 1e-14 * hb_:
                break
        k_t = rho_t * C.kH(mu)
        # final object: cycle at exactly k_t; polish with one fixed-k Newton if well conditioned
        cyc = dict(c)
        cyc['k'] = k_t
        if c['m_t'] < 0.9:
            try:
                cf = C.newton_fixed_k(mu, k_t, c['U0'], c['y0'], maxit=6)
                if cf['residual'] < 1e-10:
                    cyc = cf
            except Exception as e:
                log(f"  fixed-k polish failed at mu={mu} rho={rho_t}: {e}")
        else:
            U1, y1, T, M, _ = C.crest_return_var(mu, k_t, c['U0'], c['y0'])
            cyc.update(T=T, M=M, m_t=C.multipliers(M, mu, T)[0], residual=float(np.hypot(U1 - c['U0'], (y1 - c['y0']) / k_t)))
        r = describe(mu, cyc)
        r['secant_iters'] = it + 1
        rows.append(r)
        log(f"  rho-1={rho_t-1:10.4e} h={r['U_max']:8.4f} Umin={r['U_min']:12.4f} T={r['T']:12.4f} log_mt={r['log_mt']:+.6e} (M {r['log_mt_M']:+.6e})"
            f" lam_c={r['lambda_c']:+.4e} <U>+log rho={r['meanU_plus_log_rho']:+.4e} res={r['newton_residual']:.1e} clos={r['closure']:.1e} set={r['settle_periods']}")
    return brows, rows


def main():
    mus = MUS if len(sys.argv) < 2 else [float(a) for a in sys.argv[1:]]
    tag = '' if len(sys.argv) < 2 else '_partial'
    logf = open(f'grid_log{tag}.txt', 'w')

    def log(s):
        print(s, flush=True)
        logf.write(s + '\n')
        logf.flush()
    allb, allr = [], []
    for mu in mus:
        t0 = time.time()
        try:
            b, r = run_mu(mu, log)
        except Exception as e:
            log(f"mu={mu}: FAILED {e!r}")
            continue
        allb += b
        allr += r
        log(f"mu={mu}: done in {time.time()-t0:.1f}s")
    for name, rows in (('branch_table', allb), ('grid_table', allr)):
        if not rows:
            continue
        keys = list(rows[-1].keys())
        with open(f'{name}{tag}.csv', 'w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=keys, extrasaction='ignore', restval='')
            w.writeheader()
            for r in rows:
                w.writerow({kk: (repr(float(v)) if isinstance(v, (float, np.floating)) else v) for kk, v in r.items()})
    log('written')


if __name__ == '__main__':
    main()
