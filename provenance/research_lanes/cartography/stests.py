#!/usr/bin/env python3
"""Part 3: structural statements S1-S6 on the attractor disc, with pass/fail and worst cases.
Usage: stests.py [out.json]   (parameter list is fixed below)"""
import sys, json, time
sys.path.insert(0, '/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03/work/lanes/cartography')
import numpy as np
import cart as C
from attractor import sample_attractor

COORDS = ['U', 'y-k', 'e^U-1', 'sqrtK', 'sqrtW', 'Lambda_E']


def chi(U):
    U = np.asarray(U, dtype=float)
    out = np.empty_like(U)
    small = np.abs(U) < 1e-3
    out[small] = U[small] ** 2 / 6.0
    Ub = U[~small]
    out[~small] = Ub - 2.0 * (np.expm1(Ub) - Ub) / np.expm1(Ub)
    return out


def g2(U):
    U = np.asarray(U, dtype=float)
    out = np.where(np.abs(U) < 1e-3, U * U / 2 + U ** 3 / 3, U * np.exp(U) - np.expm1(U))
    return out


def coords(mu, k, ed, U, y, eta, zeta, K, W, zU, eU):
    """value and derivative along the crest section for each coordinate"""
    etad = zeta + (mu - 1) * eta - eta * eta
    etadU = zU + (mu - 1) * eU - 2 * eta * eU
    Wc = etad * (y - k)                         # W on the section p = 0
    WU = etadU * (y - k) + etad * zeta
    Kc = zeta * U - y + k
    KU = U * zU
    with np.errstate(invalid='ignore', divide='ignore'):
        val = {'U': U, 'y-k': y - k, 'e^U-1': np.expm1(U), 'sqrtK': np.sqrt(Kc), 'sqrtW': np.sqrt(Wc),
               'Lambda_E': (y - k) - ed['zetaE'] * U}
        der = {'U': np.ones_like(U), 'y-k': zeta, 'e^U-1': np.exp(U), 'sqrtK': KU / (2 * np.sqrt(Kc)),
               'sqrtW': WU / (2 * np.sqrt(Wc)), 'Lambda_E': zeta - ed['zetaE']}
    return val, der, Kc, Wc


def worst(arr, where, n=1):
    if len(arr) == 0:
        return None
    i = int(np.nanargmax(arr))
    return dict(value=float(arr[i]), at=float(where[i]))


def analyse(S):
    mu, k, rho = S['mu'], S['k'], S['rho']
    A = 1 + mu
    ed = S['eq']
    Us = S['cycle']['U0']
    TT, P = S['turns'], S['points']
    res = dict(mu=mu, rho=rho, k=k, U_star=Us, log_mt=S['profile']['log_mt'], lambda_c=S['profile']['lambda_c'],
               n_turns=int(len(TT['U0'])), n_points=int(len(P['U'])), growth_per_turn_at_E=S['growth'])
    U0, U1 = TT['U0'], TT['U1']
    tolU = 1e-7 * max(1.0, Us)
    inside = U0 < Us - tolU
    outside = U0 > Us + tolU
    res['n_inside'] = int(inside.sum())
    res['n_outside'] = int(outside.sum())
    res['U0_over_Ustar_range'] = [float(U0.min() / Us), float(U0.max() / Us)]
    # ---- consistency of the pulled-back geometry (brief (b))
    res['b_check_max_rel_dy/dU_minus_zeta'] = float(np.max(np.abs(TT['yU0'] - TT['zeta0']) / TT['zeta0']))
    zmax = np.max(P['zeta'])
    cond = (P['zeta'] > 1e-9 * zmax) & (P['eta'] - max(mu - 1.0, 0.0) > 1e-7 * max(mu, 1.0))
    res['hessian_well_conditioned_fraction'] = float(cond.mean())
    res['b_check_hessian_max_asym_on_conditioned'] = float(np.nanmax(P['asym'][cond]))
    with np.errstate(invalid='ignore', divide='ignore'):
        conv_turn = (np.abs(TT['yU0'] / TT['zeta0'] - 1) < 1e-6) & (np.abs(TT['yU1'] / TT['zeta1'] - 1) < 1e-6)
    res['variation_converged_turn_fraction'] = float(conv_turn.mean())
    T_c = S['cycle']['T']
    res['smoothness_exponent_r=log(m_s)/log(m_t)'] = float(A * T_c / abs(S['profile']['log_mt']) - 1)
    good = np.isfinite(P['asym']) & (P['asym'] < 1e-6) & cond & conv_turn[P['turn'].astype(int)]
    res['b_check_hessian_symmetric_fraction'] = float(good.mean()) if len(good) else None
    res['b_check_min_K_over_k'] = float(np.min(P['K']) / k)
    res['b_check_min_W'] = float(np.min(P['W']))
    # ---- S1: one-turn average of eta - mu
    a = TT['I_eta'] / TT['T'] - mu
    b = TT['I_eta'] - mu * TT['T']
    s1 = {}
    # literal: consecutive turns of the same orbit
    same = (TT['orbit'][1:] == TT['orbit'][:-1]) & (np.abs(TT['U0'][1:] - TT['U1'][:-1]) < 1e-12 * max(1.0, Us))
    cons_in = same & inside[:-1] & inside[1:]
    da = a[1:] - a[:-1]
    s1['consecutive_pairs_inside'] = int(cons_in.sum())
    if cons_in.any():
        viol = cons_in & (da > 1e-11 * (1 + np.abs(a[:-1])))
        s1['literal_violations'] = int(viol.sum())
        j = np.where(cons_in)[0][np.argmax(da[cons_in])]
        s1['literal_worst_increase'] = dict(da=float(da[j]), U_over_Ustar=float(U0[j] / Us), a_n=float(a[j]), a_np1=float(a[j + 1]))
    # as a function of crest amplitude
    for name, sel in (('inside', inside), ('outside', outside)):
        if sel.sum() < 3:
            continue
        o = np.argsort(U0[sel])
        Uo, ao, bo = U0[sel][o], a[sel][o], b[sel][o]
        keep = np.concatenate(([True], np.diff(Uo) > 1e-6 * max(1.0, Us)))
        Uo, ao, bo = Uo[keep], ao[keep], bo[keep]
        if len(Uo) < 3:
            continue
        d = np.diff(ao)
        s1[f'avg_monotone_decreasing_in_U_{name}'] = bool(np.all(d < 1e-11 * (1 + np.abs(ao[:-1]))))
        s1[f'avg_n_increases_{name}'] = int((d > 1e-11 * (1 + np.abs(ao[:-1]))).sum())
        i = int(np.argmax(d))
        s1[f'avg_worst_increase_{name}'] = dict(da=float(d[i]), U_over_Ustar=float(Uo[i] / Us), a=float(ao[i]))
        d = np.diff(bo)
        s1[f'integral_monotone_decreasing_in_U_{name}'] = bool(np.all(d < 1e-11 * (1 + np.abs(bo[:-1]))))
    s1['a_at_E'] = float(ed['etaE'] - mu)
    s1['a_at_cycle'] = float(S['profile']['lambda_c'])
    s1['pass_literal'] = bool(s1.get('literal_violations', 0) == 0) if cons_in.any() else None
    s1['pass_function_of_U_inside'] = s1.get('avg_monotone_decreasing_in_U_inside')
    res['S1'] = s1
    # ---- S2: return map in several coordinates
    phi_p = (TT['y0'] - k) / (TT['y1'] - k) * np.exp(TT['I_eta'] - mu * TT['T'])
    v0, d0, K0c, W0c = coords(mu, k, ed, TT['U0'], TT['y0'], TT['eta0'], TT['zeta0'], TT['K0'], TT['W0'], TT['zU0'], TT['eU0'])
    v1, d1, K1c, W1c = coords(mu, k, ed, TT['U1'], TT['y1'], TT['eta1'], TT['zeta1'], TT['K1'], TT['W1'], TT['zU1'], TT['eU1'])
    res['b_check_max_rel_K_crest_identity'] = float(np.max(np.abs(K0c - TT['K0']) / np.abs(TT['K0'])))
    res['b_check_max_rel_W_crest_identity'] = float(np.max(np.abs(W0c - TT['W0']) / np.abs(TT['W0'])))
    s2 = {}
    for cname in COORDS:
        c0, c1, dc0, dc1 = v0[cname], v1[cname], d0[cname], d1[cname]
        with np.errstate(invalid='ignore', divide='ignore'):
            phic = dc1 / dc0 * phi_p
            elast = c0 / c1 * phic
        e = dict(coordinate_monotone=bool(np.all(dc0[np.isfinite(dc0) & conv_turn] > 0)))
        for name, sel in (('inside', inside), ('outside', outside)):
            sel = sel & np.isfinite(phic) & np.isfinite(elast)
            if cname in ('sqrtK', 'sqrtW'):
                sel = sel & conv_turn
            if sel.sum() < 4:
                continue
            o = np.argsort(U0[sel])
            Uo, ph, el = U0[sel][o], phic[sel][o], elast[sel][o]
            keep = np.concatenate(([True], np.diff(Uo) > 1e-5 * max(1.0, Us)))
            Uo, ph, el = Uo[keep], ph[keep], el[keep]
            if len(Uo) < 3:
                continue
            dph = np.diff(ph)
            thr = 1e-9 * (1 + np.abs(ph[:-1]))
            inc = dph > thr
            dec = dph < -thr
            e[f'concave_{name}'] = bool(not inc.any())
            e[f'convex_{name}'] = bool(not dec.any())
            e[f'frac_convex_intervals_{name}'] = float(inc.mean())
            i = int(np.argmax(dph / np.diff(Uo)))
            e[f'max_phi2_{name}'] = dict(dphi1_dU=float(dph[i] / (Uo[i + 1] - Uo[i])), U_over_Ustar=float(Uo[i] / Us))
            e[f'elasticity_lt_1_{name}'] = bool(np.all(el < 1 + 1e-9))
            i = int(np.argmax(el))
            e[f'max_elasticity_{name}'] = dict(value=float(el[i]), U_over_Ustar=float(Uo[i] / Us))
            away = Uo > 0.01 * Us
            if away.any():
                e[f'max_elasticity_{name}_U>0.01U*'] = float(el[away].max())
                e[f'elasticity_lt_1_{name}_U>0.01U*'] = bool(el[away].max() < 1 - 1e-12)
            # where is phi'' > 0 : report U/U* ranges
            if inc.any():
                e[f'convex_U_over_Ustar_range_{name}'] = [float(Uo[:-1][inc].min() / Us), float(Uo[1:][inc].max() / Us)]
        s2[cname] = e
    res['S2'] = s2
    # ---- S3: Zhang-type monotonicity on the section p = 0
    s3 = {}
    for name, U, zeta, eta, zU, eU in (('crest_U>0', TT['U0'], TT['zeta0'], TT['eta0'], TT['zU0'], TT['eU0']),
                                       ('trough_U<0', TT['Ut'], TT['zetat'], TT['etat'], TT['zUt'], TT['eUt'])):
        ok = np.isfinite(zU) & np.isfinite(eU) & (zeta > 1e-150) & conv_turn
        U, zeta, eta, zU, eU = U[ok], zeta[ok], eta[ok], zU[ok], eU[ok]
        Us_here = TT['U0'][ok]
        eUx = np.exp(U)
        em1 = np.expm1(U)
        # d/dU [(mu - zeta)/(k(e^U-1))] >= 0   <=>   Q := (zeta - mu) e^U - zU (e^U - 1) >= 0
        Q = (zeta - mu) * eUx - zU * em1
        Qn = Q / (np.abs(zeta - mu) * eUx + np.abs(zU * em1) + 1e-300)
        s3[f'S3a_{name}_pass'] = bool(np.all(Qn > -1e-9))
        s3[f'S3a_{name}_frac_fail'] = float((Qn < -1e-9).mean())
        i = int(np.argmin(Qn))
        s3[f'S3a_{name}_worst'] = dict(normalised=float(Qn[i]), U=float(U[i]), crest_U_over_Ustar=float(Us_here[i] / Us))
        # Lienard f = (mu - zeta)/(1+eta) along the section
        f = (mu - zeta) / (1 + eta)
        fU = (-zU * (1 + eta) - (mu - zeta) * eU) / (1 + eta) ** 2
        QL = fU * em1 - f * eUx
        QLn = QL / (np.abs(fU * em1) + np.abs(f * eUx) + 1e-300)
        s3[f'S3L_{name}_pass'] = bool(np.all(QLn > -1e-9))
        s3[f'S3L_{name}_frac_fail'] = float((QLn < -1e-9).mean())
        i = int(np.argmin(QLn))
        s3[f'S3L_{name}_worst'] = dict(normalised=float(QLn[i]), U=float(U[i]), crest_U_over_Ustar=float(Us_here[i] / Us))
    # over the disc, derivative at fixed Lienard variable s
    okp = good & (np.abs(P['U']) > 1e-9) & (P['zeta'] > 1e-150)
    Up, pp_, eta, zeta = P['U'][okp], P['p'][okp], P['eta'][okp], P['zeta'][okp]
    hUU, hUp, hpp = P['hUU'][okp], P['hUp'][okp], P['hpp'][okp]
    PU = (zeta - mu) / (1 + eta)
    dz = hUU + PU * hUp
    de = -hUp - PU * hpp
    f = (mu - zeta) / (1 + eta)
    df = (-dz * (1 + eta) - (mu - zeta) * de) / (1 + eta) ** 2
    em1 = np.expm1(Up)
    QD = df * em1 - f * np.exp(Up)
    QDn = QD / (np.abs(df * em1) + np.abs(f * np.exp(Up)) + 1e-300)
    turn_U0 = TT['U0'][P['turn'][okp].astype(int)]
    for name, sel in (('inside', turn_U0 < Us - tolU), ('outside', turn_U0 > Us + tolU)):
        if sel.sum() == 0:
            continue
        s3[f'S3b_disc_{name}_pass'] = bool(np.all(QDn[sel] > -1e-7))
        s3[f'S3b_disc_{name}_frac_fail'] = float((QDn[sel] < -1e-7).mean())
        i = int(np.argmin(QDn[sel]))
        s3[f'S3b_disc_{name}_worst'] = dict(normalised=float(QDn[sel][i]), U=float(Up[sel][i]), p_over_sqrt_mu=float(pp_[sel][i] / np.sqrt(mu)),
                                            crest_U_over_Ustar=float(turn_U0[sel][i] / Us))
    # necessary for S3b near U = 0: zeta >= mu on the line U = 0
    nearU0 = (np.abs(P['U']) < 0.02 * max(1.0, Us)) & (P['zeta'] > 1e-150)
    if nearU0.any():
        zz = P['zeta'][nearU0] - mu
        pz = P['p'][nearU0]
        s3['zeta_minus_mu_on_U=0: min_over_mu'] = float(zz.min() / mu)
        s3['zeta_minus_mu_on_U=0: negative_fraction_p>0'] = float((zz[pz > 0] < 0).mean()) if (pz > 0).any() else None
        s3['zeta_minus_mu_on_U=0: negative_fraction_p<0'] = float((zz[pz < 0] < 0).mean()) if (pz < 0).any() else None
    res['S3'] = s3
    # ---- S4: Hessian rank
    okh = good
    hUU, hUp, hpp = P['hUU'][okh], P['hUp'][okh], P['hpp'][okh]
    # scale-free: use p/sqrt(mu)-scaled Hessian so both directions are comparable
    sm = np.sqrt(mu)
    a11, a12, a22 = hUU, hUp * sm, hpp * mu
    tr = a11 + a22
    det = a11 * a22 - a12 * a12
    disc = np.sqrt(np.maximum(tr * tr - 4 * det, 0))
    lmin = 0.5 * (tr - disc)
    lmax = 0.5 * (tr + disc)
    with np.errstate(invalid='ignore', divide='ignore'):
        ratio = np.where(lmax > 0, lmin / lmax, -1.0)
    turn_U0 = TT['U0'][P['turn'][okh].astype(int)]
    s4 = dict(min_eig_ratio=float(ratio.min()), frac_negative_det=float((det < 0).mean()),
              median_eig_ratio=float(np.median(ratio)), min_trace=float(tr.min()))
    i = int(np.argmin(ratio))
    s4['worst'] = dict(eig_ratio=float(ratio[i]), U=float(P['U'][okh][i]), p_over_sqrt_mu=float(P['p'][okh][i] / sm),
                       crest_U_over_Ustar=float(turn_U0[i] / Us))
    s4['strictly_convex_pass'] = bool(ratio.min() > 0)
    near = np.abs(turn_U0 / Us - 1) < 2e-3
    if near.any():
        s4['min_eig_ratio_near_cycle'] = float(ratio[near].min())
        s4['median_eig_ratio_near_cycle'] = float(np.median(ratio[near]))
    res['S4'] = s4
    # ---- S5: sign of (1+eta)E = K - (zeta - mu) chi(U)
    Eh = P['K'] - (P['zeta'] - mu) * chi(P['U'])
    turn_U0 = TT['U0'][P['turn'].astype(int)]
    s5 = {}
    for name, sel in (('inside', turn_U0 < Us - tolU), ('outside', turn_U0 > Us + tolU), ('near_cycle', np.abs(turn_U0 / Us - 1) < 2e-3)):
        if sel.sum() == 0:
            continue
        e = Eh[sel]
        s5[f'frac_negative_{name}'] = float((e < 0).mean())
        s5[f'min_over_k_{name}'] = float(e.min() / k)
        s5[f'max_over_k_{name}'] = float(e.max() / k)
        neg = e < 0
        if neg.any():
            Un, pn = P['U'][sel][neg], P['p'][sel][neg]
            s5[f'negative_region_{name}'] = dict(U_range=[float(Un.min()), float(Un.max())],
                                                 frac_p_positive=float((pn > 0).mean()), frac_U_positive=float((Un > 0).mean()))
        pos = e > 0
        if pos.any():
            Un, pn = P['U'][sel][pos], P['p'][sel][pos]
            s5[f'positive_region_{name}'] = dict(U_range=[float(Un.min()), float(Un.max())],
                                                 frac_p_positive=float((pn > 0).mean()), frac_U_positive=float((Un > 0).mean()))
    s5['one_signed'] = bool((Eh >= 0).all() or (Eh <= 0).all())
    res['S5'] = s5
    # ---- S6: pointwise ratios
    with np.errstate(divide='ignore', invalid='ignore'):
        ze = P['zeta'] / P['eta']
        R1 = k * np.exp(P['U']) / P['zeta']
        R2 = ze + mu
        R3 = k * g2(P['U']) / P['K']
        R4 = k * np.exp(P['U']) * P['p'] ** 2 / P['W']
    okr = (P['zeta'] > 1e-150) & (P['eta'] > 1e-150)
    s6 = dict(min_zeta_over_eta=float(np.nanmin(ze[okr])), max_zeta_over_eta=float(np.nanmax(ze[okr])),
              frac_zeta_over_eta_le_1=float((ze[okr] <= 1).mean()), zetaE_over_etaE=float(ed['zetaE'] / ed['etaE']))
    for nm, R in (('x/zeta', R1), ('zeta/eta+mu', R2), ('k g2(U)/K', R3), ('x p^2/W', R4)):
        r = R[okr] - A
        s6[f'{nm}: min-A'] = float(np.nanmin(r))
        s6[f'{nm}: max-A'] = float(np.nanmax(r))
        s6[f'{nm}: sign_changing'] = bool(np.nanmin(r) < 0 < np.nanmax(r))
    res['S6'] = s6
    return res


POINTS = [(mu, 1 + d) for mu in (0.002, 0.01, 0.05, 0.2, 1.0, 5.0, 50.0) for d in (0.003, 0.03, 0.3, 1.5, 2.5, 4.0, 8.0, 20.0)] + [(0.0027, 5.02), (0.002, 13.92), (500.0, 34.0), (0.5, 5.0), (2.0, 5.0), (0.5, 12.0)]


def main():
    outname = sys.argv[1] if len(sys.argv) > 1 else 'stests_results.json'
    pts = POINTS
    if len(sys.argv) > 2:
        pts = [tuple(float(v) for v in a.split(',')) for a in sys.argv[2:]]
    allres = []
    for mu, rho in pts:
        t0 = time.time()
        try:
            S = sample_attractor(mu, rho)
            r = analyse(S)
            r['seconds'] = time.time() - t0
        except Exception as e:
            import traceback
            traceback.print_exc()
            r = dict(mu=mu, rho=rho, error=repr(e))
        allres.append(r)
        if 'error' in r:
            print(f"mu={mu} rho={rho}: ERROR {r['error']}", flush=True)
            continue
        s2 = r['S2']
        print(f"mu={mu:<7} rho-1={rho-1:<7.4g} U*={r['U_star']:.3f} log_mt={r['log_mt']:+.3e} turns in/out={r['n_inside']}/{r['n_outside']} symH={r['b_check_hessian_symmetric_fraction']:.3f}"
              f" | S1 lit={r['S1'].get('pass_literal')} fn={r['S1'].get('pass_function_of_U_inside')}"
              f" | S2 concave(in) " + ' '.join(f"{c}:{'Y' if s2[c].get('concave_inside') else 'n'}" for c in COORDS)
              + " | elast<1(in,U>.01U*/out) " + ' '.join(f"{c}:{'Y' if s2[c].get('elasticity_lt_1_inside_U>0.01U*') else 'n'}{'Y' if s2[c].get('elasticity_lt_1_outside') else ('-' if s2[c].get('elasticity_lt_1_outside') is None else 'n')}" for c in COORDS)
              + f" | S3a c/t={r['S3']['S3a_crest_U>0_pass']}/{r['S3']['S3a_trough_U<0_pass']} S3L={r['S3']['S3L_crest_U>0_pass']}/{r['S3']['S3L_trough_U<0_pass']} S3b={r['S3'].get('S3b_disc_inside_pass')}"
              f" | S4 minratio={r['S4']['min_eig_ratio']:.2e} | S5 neg in={r['S5'].get('frac_negative_inside', float('nan')):.2f} cyc={r['S5'].get('frac_negative_near_cycle', float('nan')):.2f}"
              f" | S6 min z/e={r['S6']['min_zeta_over_eta']:.3f} t={r['seconds']:.1f}s", flush=True)
    json.dump(allres, open(outname, 'w'), indent=1)


if __name__ == '__main__':
    main()
