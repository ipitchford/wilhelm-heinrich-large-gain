#!/usr/bin/env python3
"""Referee test of lane 'monotone' claims 1-4 (DS planes, twelve-event order, rank-one window identity,
strong-stable ratios).  Independent of the lane's scripts.  All sign decisions are made in SIGNED-LOG
arithmetic (monomials with log-magnitudes), so they survive troughs with U ~ -8000, y ~ e^{-900}.
usage: v1_floquet.py [mu rho] ...   (default: the 12-point extreme grid + mu>1 cases)
"""
import sys, time, warnings, itertools
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
sys.path.insert(0, '/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03/work/lanes/monotone/verify')
from vlib import *
warnings.filterwarnings('ignore')
np.seterr(all='ignore')


# ---------- signed-log polynomials -------------------------------------------------------------
class SL:
    def __init__(self, terms):
        self.terms = terms          # list of (sign, log|.|) broadcastable arrays

    @staticmethod
    def pos(lg):
        return SL([(1.0, lg)])

    @staticmethod
    def real(v):
        v = np.asarray(v, float)
        return SL([(np.sign(v), np.log(np.abs(v)))])

    def __neg__(self):
        return SL([(-s, l) for s, l in self.terms])

    def __add__(self, o):
        return SL(self.terms + o.terms)

    def __sub__(self, o):
        return self + (-o)

    def __mul__(self, o):
        return SL([(s1 * s2, l1 + l2) for (s1, l1), (s2, l2) in itertools.product(self.terms, o.terms)])

    def ev(self):
        """ratio = sum / sum|.| in [-1,1]  and  log(sum|.|)"""
        shp = np.broadcast_shapes(*[np.broadcast(np.asarray(s), np.asarray(l)).shape for s, l in self.terms])
        S = np.stack([np.broadcast_to(np.asarray(s, float), shp) for s, _ in self.terms])
        Lg = np.stack([np.broadcast_to(np.asarray(l, float), shp) for _, l in self.terms])
        Lg = np.where(S == 0, -np.inf, Lg)
        Lm = np.max(Lg, axis=0)
        Lm = np.where(np.isfinite(Lm), Lm, 0.0)
        w = np.exp(Lg - Lm)
        tot = np.sum(S * w, axis=0); ab = np.sum(w, axis=0)
        return np.where(ab > 0, tot / np.where(ab > 0, ab, 1), 0.0), Lm + np.log(ab)


def atoms(mu, k, U, V, Q, r, ls):
    lam = k / mu
    x = SL.pos(np.log(k) + U); y = SL.pos(V); q = SL.pos(Q)
    R = SL.pos(np.log(r)); S = SL.pos(ls)
    LAM = SL.pos(np.log(lam) + 0 * U); MU = SL.pos(np.log(mu) + 0 * U)
    p = LAM - q
    pp = MU * q - y
    yp = x - y
    q1 = R * yp - S * pp
    q2 = yp + S * p
    q3 = -(pp + R * p)
    ppp = -(MU * pp) - yp
    ypp = x * p - yp
    return dict(x=x, y=y, q=q, p=p, pp=pp, yp=yp, q1=q1, q2=q2, q3=q3, ppp=ppp, ypp=ypp)


def adapted_grid(cyc, n, nfine=200001):
    T = cyc['T']
    tf = np.linspace(0, T, nfine)
    w = cyc['sol'](tf)
    sc = np.ptp(w, axis=1)[:, None]
    m = np.sqrt(np.sum((np.diff(w, axis=1) / sc) ** 2, axis=0)) + 1e-3 / nfine
    c = np.concatenate(([0], np.cumsum(m))); c /= c[-1]
    return np.interp((np.arange(n) + 0.5) / n, c, tf)


def sign_changes(tg, ratio):
    """list of (t_lo, t_hi, direction) for sign changes of a sampled signed quantity"""
    sg = np.sign(ratio)
    out = []
    idx = np.nonzero(sg[:-1] * sg[1:] < 0)[0]
    for i in idx:
        out.append((tg[i], tg[i + 1], +1 if sg[i + 1] > 0 else -1))
    return out, int(np.sum(sg == 0))


def run(mu, rho, quick=False):
    t00 = time.time()
    k = rho * kH(mu); lam = k / mu; A = 1 + mu
    cyc = find_cycle(mu, k)
    T = cyc['T']
    out = dict(mu=mu, rho=rho, T=T)
    # ---- strong-stable ratios (backward) and adjoint ratios (forward), with dense interpolants
    ed = eq_data(mu, k)

    def gb(tau, w):
        U = cyc['sol'](np.mod(T - tau, T))[0, 0]
        r, ls, Ir = w
        return [-(r * r - mu * r - np.exp(ls)), -((r - 1.0) - k * np.exp(U - ls)), r]
    w0 = np.array([ed['gE'] * 1.7, np.log(0.3 * ed['gE'] * (ed['gE'] - mu)), 0.0])   # deliberately NOT the equilibrium values
    for _ in range(3):
        sb = solve_ivp(gb, [0, T], w0, method='DOP853', rtol=RT, atol=AT, dense_output=True)
        w1 = sb.y[:, -1].copy(); rs_defect = np.abs(w1[:2] - w0[:2]).max(); mean_r = w1[2] / T
        w0 = np.array([w1[0], w1[1], 0.0])
    rs = sb.sol

    def r_ls(t):
        z = rs(T - np.mod(t, T))
        return z[0], z[1]

    def gf(t, w):
        U = cyc['sol'](np.mod(t, T))[0, 0]
        eta, lz, Ie = w
        return [np.exp(lz) + (mu - 1) * eta - eta * eta, k * np.exp(U - lz) - (1 + eta), eta]
    w0 = np.array([ed['etaE'], np.log(ed['zetaE']), 0.0])
    for _ in range(4):
        sf = solve_ivp(gf, [0, T], w0, method='DOP853', rtol=RT, atol=AT)
        w1 = sf.y[:, -1].copy(); ez_defect = np.abs(w1[:2] - w0[:2]).max(); mean_eta = w1[2] / T
        w0 = np.array([w1[0], w1[1], 0.0])
    zeta_crest = float(np.exp(w0[1])); eta_crest = float(w0[0])
    lam_eta = mean_eta - mu
    lam_r = mean_r - A
    # return-map multiplier where finite differences can resolve it
    try:
        M = jac_fd(mu, k, cyc['U0'], cyc['y0'])
        ev = np.linalg.eigvals(M)
        m_t_fd = float(np.max(ev.real)); lam_fd = np.log(m_t_fd) / T if m_t_fd > 1e-5 else np.nan
        i = int(np.argmax(ev.real)); w_, V_ = np.linalg.eig(M); vv = V_[:, int(np.argmax(w_.real))].real
        perron_slope = vv[1] / vv[0] if m_t_fd > 1e-5 else np.nan
        Mpos = bool((M > 0).all()) if m_t_fd > 1e-5 else None
    except Exception:
        lam_fd = np.nan; perron_slope = np.nan; Mpos = None
    out.update(lam_eta=lam_eta, lam_r=lam_r, lam_fd=lam_fd, mean_r=mean_r, mean_eta=mean_eta,
               rs_defect=rs_defect, ez_defect=ez_defect, zeta_crest=zeta_crest, perron_slope=perron_slope, Mpos=Mpos)

    # ---- claim 4: r > max(1,mu), s > 0 ; margins in log variables
    tg = np.unique(np.concatenate((np.linspace(0, T, 40001)[:-1] + 0.5 * T / 40000, adapted_grid(cyc, 20000))))
    tg = tg[(tg > 0) & (tg < T)]
    U, V, Q = cyc['sol'](tg)
    r, ls = r_ls(tg)
    out['min_r_minus_max'] = float(np.min(r - max(1.0, mu)))
    out['min_ln_s'] = float(np.min(ls))
    # exact margin: delta = r - max(1,mu), sigma = s~ (mu<=1) or s (mu>=1), backward in log variables
    m0 = max(1.0, mu)

    def gm(tau, w):
        U_ = cyc['sol'](np.mod(T - tau, T))[0, 0]
        ld, lsg = w
        d = np.exp(ld)
        if mu <= 1:   # delta' = delta(1+delta) - sig ; sig' = sig(1-mu+delta) - x
            return [-(1 + d - np.exp(lsg - ld)), -((1 - mu + d) - k * np.exp(U_ - lsg))]
        else:         # delta' = delta(mu+delta) - s ; s' = s(mu-1+delta) - x
            return [-(mu + d - np.exp(lsg - ld)), -((mu - 1 + d) - k * np.exp(U_ - lsg))]
    w0 = np.array([0.0, 0.0])
    for _ in range(3):
        sm = solve_ivp(gm, [0, T], w0, method='DOP853', rtol=1e-11, atol=1e-12, dense_output=True)
        w1 = sm.y[:, -1].copy(); mg_defect = np.abs(w1 - w0).max(); w0 = w1
    ld = sm.sol(T - tg)[0]
    out['margin_min_log10'] = float(np.min(ld) / np.log(10)); out['mg_defect'] = mg_defect
    # consistency of the margin system with the (r, ln s) system where double precision resolves it
    dd = r - m0
    msk = dd > 1e-6
    out['margin_consistency'] = float(np.max(np.abs(np.exp(ld[msk]) / dd[msk] - 1))) if msk.any() else np.nan

    # ---- claim 2: twelve-event order
    at = atoms(mu, k, U, V, Q, r, ls)
    sig = {}
    for nm in ['p', 'pp', 'yp', 'q1', 'q2', 'q3']:
        sig[nm] = at[nm].ev()[0]
    ev_list = []
    counts = {}
    for nm in ['p', 'pp', 'yp', 'q1', 'q2', 'q3']:
        sc, nz = sign_changes(np.concatenate((tg, [tg[0] + T])), np.concatenate((sig[nm], [sig[nm][0]])))
        counts[nm] = len(sc)
        for lo, hi, d in sc:
            ev_list.append((0.5 * (lo + hi), nm, d))
    ev_list.sort()
    seq = [f"{nm}{'+' if d > 0 else '-'}" for _, nm, d in ev_list]
    claimed = ['p+', 'q2+', 'yp+', 'q1+', 'pp-', 'q3+', 'p-', 'q2-', 'yp-', 'q1-', 'pp+', 'q3-']
    ok12 = False
    if len(seq) == 12 and 'p+' in seq:
        i0 = seq.index('p+')
        ok12 = (seq[i0:] + seq[:i0]) == claimed
    out['counts'] = counts; out['order_ok'] = ok12; out['seq'] = seq
    ev_t = {f"{nm}{'+' if d > 0 else '-'}": tt for tt, nm, d in ev_list}
    # minimal gap between consecutive events relative to grid spacing there (resolution check)
    tt_sorted = np.array([e[0] for e in ev_list])
    gaps = np.diff(np.concatenate((tt_sorted, [tt_sorted[0] + T])))
    out['min_event_gap'] = float(gaps.min())
    out['events'] = [(f'{tt / T:.9f}', nm, d) for tt, nm, d in ev_list]
    # signs at crest (t=0) and at trough
    U0_, V0_, Q0_ = cyc['sol'](np.array([0.0])); r0_, ls0_ = r_ls(np.array([0.0]))
    a0_ = atoms(mu, k, U0_, V0_, Q0_, r0_, ls0_)
    out['Psi_crest_pos'] = bool(a0_['q1'].ev()[0][0] > 0 and a0_['q2'].ev()[0][0] > 0 and a0_['q3'].ev()[0][0] > 0)
    # exact trough: root of Q - ln(lam) (decreasing), then impose p = 0 exactly
    lo_, hi_ = [(a_, b_) for a_, b_, d_ in sign_changes(np.concatenate((tg, [tg[0] + T])), np.concatenate((sig['p'], [sig['p'][0]])))[0] if d_ > 0][0]
    ttr = brentq(lambda t_: float(cyc['sol'](np.mod(np.array([t_]), T))[2, 0]) - np.log(lam), lo_, hi_, xtol=1e-14 * T)
    Ut_, Vt_, Qt_ = cyc['sol'](np.mod(np.array([ttr]), T)); rt_, lst_ = r_ls(np.array([ttr]))
    # at p = 0 exactly:  q1 = r y' - s p',  q2 = y',  q3 = -p'   with  y' = x - y,  p' = k - y   (signed-log, no absorbed terms)
    x_ = SL.pos(np.log(k) + Ut_); y_ = SL.pos(Vt_); K_ = SL.pos(np.log(k) + 0 * Ut_)
    yp_ = x_ - y_; pp_ = K_ - y_
    q1_ = SL.pos(np.log(rt_)) * yp_ - SL.pos(lst_) * pp_
    out['Psi_trough_neg'] = bool(q1_.ev()[0][0] < 0 and yp_.ev()[0][0] < 0 and (-pp_).ev()[0][0] < 0)
    out['yp_trough_neg'] = bool(yp_.ev()[0][0] < 0)
    out['yp_crest_pos'] = bool(a0_['yp'].ev()[0][0] > 0)
    # ties: events sharing a grid cell (unresolvable in double-precision TIME); only (p,q2) ties are decidable analytically:
    # at a zero of p, Psi_2 = varpi*y', so 'p+ before q2+' <=> y'<0 at the trough, 'p- before q2-' <=> y'>0 at the crest.
    ties = []
    cyc_ev = ev_list + [(ev_list[0][0] + T,) + ev_list[0][1:]]
    for (ta_, na_, da_), (tb_, nb_, db_) in zip(cyc_ev[:-1], cyc_ev[1:]):
        if tb_ - ta_ < 1e-9 * T:
            ties.append((na_, nb_))
    out['ties'] = ties
    out['ties_ok'] = all(set(tp) == {'p', 'q2'} for tp in ties) and out['yp_trough_neg'] and out['yp_crest_pos']
    # ---- refine t1 (zero of q3 in F) and t2
    def q3f(t):
        U_, V_, Q_ = cyc['sol'](np.mod(np.array([t]), T))
        r_, ls_ = r_ls(np.array([t]))
        return float(atoms(mu, k, U_, V_, Q_, r_, ls_)['q3'].ev()[0][0])
    lo, hi = [(a, b) for a, b, d in sign_changes(np.concatenate((tg, [tg[0] + T])), np.concatenate((sig['q3'], [sig['q3'][0]])))[0] if d < 0][0]
    t1 = brentq(q3f, lo, hi, xtol=1e-13 * T)
    out['t1_over_T'] = t1 / T
    s1 = state(mu, k, cyc, np.array([t1])); r1, ls1 = r_ls(np.array([t1])); r1 = float(r1[0]); sv1 = float(np.exp(ls1[0]))
    U1, p1, pp1, yp1, y1_, x1_ = [float(s1[n][0]) for n in ['U', 'p', 'pp', 'yp', 'y', 'x']]
    out['t1_in_F'] = bool(p1 < 0 and pp1 > 0 and yp1 < 0)
    out['U_t1'] = U1
    # Lemma 5: supporting line  l(t) = p'(t1)(U - U(t1)) - p(t1)(p - p(t1))  (plain doubles; p = lam - q)
    st = state(mu, k, cyc, tg)
    ell = pp1 * (st['U'] - U1) - p1 * (st['p'] - p1)
    dt1 = np.abs(((tg - t1 + T / 2) % T) - T / 2)
    out['ell_neg_count'] = int(np.sum(ell[dt1 > 1e-7 * T] <= 0))
    out['min_ell_over_scale'] = float(np.min(ell[dt1 > 0.01 * T]) / np.max(np.abs(ell)))
    out['L_at_E'] = float(p1 * p1 - pp1 * U1)
    out['ell_dd_t1'] = float(p1 * p1 * r1 * (r1 - mu) + p1 * yp1)

    # ---- claim 1: DS determinants in signed-log arithmetic
    nb, nc = (120, 500) if quick else (240, 1200)
    tb = np.unique(np.concatenate((np.linspace(0, T, nb // 2 + 1)[:-1], adapted_grid(cyc, nb // 2), [t1, 0.0])))
    tc = np.unique(np.concatenate((np.linspace(0, T, nc // 2 + 1)[:-1], adapted_grid(cyc, nc // 2))))
    Ub, Vb, Qb = cyc['sol'](tb); rb, lsb = r_ls(tb)
    Uc, Vc, Qc = cyc['sol'](tc)
    B = atoms(mu, k, Ub[:, None], Vb[:, None], Qb[:, None], rb[:, None], lsb[:, None])

    def det_with(dU, Qother, Vother):
        """det(d, z', e_s) with e_s=(1,-r,-s):  d_U q1 + d_p q2 + d_y q3 ;  d_p = q_base - q_other, d_y = y_other - y_base"""
        d_U = SL.real(dU)
        d_p = SL.pos(Qb[:, None] + 0 * Qother) - SL.pos(Qother + 0 * Qb[:, None])
        d_y = SL.pos(Vother + 0 * Vb[:, None]) - SL.pos(Vb[:, None] + 0 * Vother)
        return (d_U * B['q1'] + d_p * B['q2'] + d_y * B['q3']).ev()[0]
    ratio = det_with(Uc[None, :] - Ub[:, None], Qc[None, :], Vc[None, :])
    dist = np.abs(((tc[None, :] - tb[:, None] + T / 2) % T) - T / 2)
    # log-coordinate distance between the two points, to classify 'near-diagonal'
    dlog = np.sqrt((Uc[None, :] - Ub[:, None]) ** 2 + (Vc[None, :] - Vb[:, None]) ** 2 + (Qc[None, :] - Qb[:, None]) ** 2)
    off = dlog > 1e-2
    out['chord_pairs'] = int(off.sum())
    out['chord_viol'] = int(np.sum(ratio[off] > 0))
    out['chord_max_ratio'] = float(np.max(ratio[off]))
    out['chord_unres'] = int(np.sum(np.abs(ratio[off]) < 1e-8))
    ii, jj = np.nonzero(off & (np.abs(ratio) < 1e-7))
    out['chord_small_maxdt'] = float(np.max(dist[ii, jj]) / T) if len(ii) else 0.0
    out['chord_small_n'] = int(len(ii))
    out['chord_small'] = [(f'{tb[i] / T:.6f}', f'{tc[j] / T:.6f}', f'{ratio[i, j]:+.2e}', f'dlog={dlog[i, j]:.3f}') for i, j in zip(ii[:12], jj[:12])]
    # near diagonal: sign must still be negative wherever resolved
    nd = (~off) & (dist > 0)
    out['neardiag_pos_resolved'] = int(np.sum(ratio[nd] > 1e-7))
    # E
    rE = det_with(-Ub[:, None], np.log(lam) + 0 * Qb[:, None], np.log(k) + 0 * Vb[:, None])
    out['E_max_ratio'] = float(np.max(rE))
    # z''
    rC = (B['pp'] * B['q1'] + B['ppp'] * B['q2'] + B['ypp'] * B['q3']).ev()[0]
    out['curv_max_ratio'] = float(np.max(rC))
    # cross-check of the signed-log determinant against a plain numpy determinant (only where no underflow)
    stb = state(mu, k, cyc, tb); stc = state(mu, k, cyc, tc)
    chk = 0.0; nchk = 0
    rng = np.random.default_rng(1)
    for _ in range(400):
        i = rng.integers(len(tb)); j = rng.integers(len(tc))
        if not off[i, j]:
            continue
        d = np.array([stc['U'][j] - stb['U'][i], stc['p'][j] - stb['p'][i], stc['y'][j] - stb['y'][i]])
        zp = np.array([stb['p'][i], stb['pp'][i], stb['yp'][i]])
        es = np.array([1.0, -rb[i], -np.exp(lsb[i])])
        dt = np.linalg.det(np.array([d, zp, es]).T)
        if abs(ratio[i, j]) > 1e-3 and np.isfinite(dt) and dt != 0:
            nchk += 1
            chk = max(chk, float(np.sign(dt) != np.sign(ratio[i, j])))
    out['det_crosscheck_mismatch'] = chk; out['det_crosscheck_n'] = nchk

    # ---- claim 3: rank-one window identity   (varpi(t1)=1)
    om = 2 * np.pi / T
    tests = {
        '1': lambda t, s_: (1.0, 0.0, s_['x']),
        "p'": lambda t, s_: (s_['pp'], s_['ppp'], -s_['x'] * s_['p'] ** 2),
        '-U': lambda t, s_: (-s_['U'], -s_['p'], -k * (1 - np.exp(s_['U']) + s_['U'] * np.exp(s_['U']))),
        'y': lambda t, s_: (s_['y'], s_['yp'], s_['x'] * (k + s_['p'] ** 2)),
        'cos': lambda t, s_: (np.cos(om * t), -om * np.sin(om * t),
                              om ** 3 * np.sin(om * t) - A * om ** 2 * np.cos(om * t) - mu * om * np.sin(om * t) + s_['x'] * np.cos(om * t)),
        'sin3': lambda t, s_: (np.sin(3 * om * t), 3 * om * np.cos(3 * om * t),
                               -27 * om ** 3 * np.cos(3 * om * t) - 9 * A * om ** 2 * np.sin(3 * om * t) + 3 * mu * om * np.cos(3 * om * t) + s_['x'] * np.sin(3 * om * t)),
    }
    names = list(tests)

    def gw(t, w):
        s_ = state(mu, k, cyc, np.array([t])); s_ = {n: float(v[0]) for n, v in s_.items()}
        r_, ls_ = r_ls(np.array([t])); r_ = float(r_[0])
        q3_ = -(s_['pp'] + r_ * s_['p'])
        vp = np.exp(w[0])
        return [A - r_] + [vp * q3_ * tests[n](t, s_)[2] for n in names]
    sw = solve_ivp(gw, [t1, t1 + T], np.zeros(1 + len(names)), method='DOP853', rtol=1e-11, atol=1e-13 )
    lnvp_T = sw.y[0, -1]                      # = A T - int r = -lambda_c T by the lemma
    out['lam_window'] = -lnvp_T / T
    s1d = {n: float(v[0]) for n, v in s1.items()}
    Psi2_t1 = yp1 + sv1 * p1                  # varpi(t1)=1
    out['Psi2_t1_neg'] = bool(Psi2_t1 < 0)
    worst = 0.0; signs = {}; errs = {}
    for j, n in enumerate(names):
        ph, dph, _ = tests[n](t1, s1d)
        lhs = sw.y[1 + j, -1]
        for nm_l, lamc in (('eta', lam_eta),):
            rhs_ = (np.exp(-lamc * T) - 1) * abs(Psi2_t1) * (dph + r1 * ph)
            worst = max(worst, abs(lhs - rhs_) / max(abs(lhs), abs(rhs_)))
            errs[n] = (abs(lhs - rhs_) / max(abs(lhs), abs(rhs_)), lhs, rhs_)
        signs[n] = (float(np.sign(lhs)), float(np.sign(dph + r1 * ph)))
    out['window_rel_err'] = worst; out['window_signs'] = signs; out['window_errs'] = errs
    out['secs'] = time.time() - t00
    return out


if __name__ == '__main__':
    args = sys.argv[1:]
    quick = '--quick' in args
    args = [a for a in args if a != '--quick']
    if args:
        grid = [(float(args[i]), float(args[i + 1])) for i in range(0, len(args), 2)]
    else:
        grid = [(mu, rho) for mu in (1.0, 0.3, 0.02) for rho in (1.01, 1.5, 4.0, 12.0)] + [(3.0, 1.01), (3.0, 1.5), (3.0, 4.0), (3.0, 12.0), (50.0, 1.5), (50.0, 4.0), (0.1, 1.2), (7.0, 2.5)]
    for mu, rho in grid:
        o = run(mu, rho, quick=quick)
        print(f"== mu={mu} rho={rho} T={o['T']:.4f}  [{o['secs']:.0f}s]")
        print(f"   lambda_c: <eta>-mu={o['lam_eta']:+.9f}  <r>-A={o['lam_r']:+.9f}  window={o['lam_window']:+.9f}  FD return map={o['lam_fd']:+.9f}"
              f"   | periodicity defects rs={o['rs_defect']:.1e} etazeta={o['ez_defect']:.1e}")
        print(f"   C4: min(r-max(1,mu))={o['min_r_minus_max']:+.3e} (log-margin system: min log10={o['margin_min_log10']:.2f}, consistency {o['margin_consistency']:.1e}, defect {o['mg_defect']:.1e}); min ln s={o['min_ln_s']:.2f}")
        print(f"   C2: counts={o['counts']} order_ok={o['order_ok']} min event gap={o['min_event_gap']:.3e} Psi>0@crest={o['Psi_crest_pos']} Psi<0@trough={o['Psi_trough_neg']} ties={o['ties']} ties_ok={o['ties_ok']}")
        if not o['order_ok']:
            print('      events =', o['events'])
        print(f"   L5: t1/T={o['t1_over_T']:.5f} in F={o['t1_in_F']} U(t1)={o['U_t1']:.4g} #(l<=0)={o['ell_neg_count']} min l/max|l| (|t-t1|>1%T)={o['min_ell_over_scale']:+.2e} L(E)={o['L_at_E']:+.3e} l''(t1)={o['ell_dd_t1']:+.3e}")
        print(f"   C1: chords: pairs={o['chord_pairs']} violations={o['chord_viol']} max ratio={o['chord_max_ratio']:+.3e} unresolved(|ratio|<1e-8)={o['chord_unres']} near-diag resolved positives={o['neardiag_pos_resolved']};"
              f" E max ratio={o['E_max_ratio']:+.3e}; z'' max ratio={o['curv_max_ratio']:+.3e}; numpy-det cross-check mismatches={o['det_crosscheck_mismatch']:.0f}/{o['det_crosscheck_n']}")
        print(f"      chord pairs with |ratio|<1e-7: {o['chord_small_n']}, all within |dt|/T <= {o['chord_small_maxdt']:.2e} of the diagonal")
        print(f"   C3: window identity worst rel err={o['window_rel_err']:.2e}; Psi2(t1)<0={o['Psi2_t1_neg']}; signs (integral, (phi'+r phi)(t1)) = {o['window_signs']}")
        print('      per-test (relerr, lhs, rhs):', {n: (f'{e[0]:.1e}', f'{e[1]:+.4e}', f'{e[2]:+.4e}') for n, e in o['window_errs'].items()})
        print(f"   R7: Perron slope of FD crest map={o['perron_slope']:.8f} vs zeta(crest)={o['zeta_crest']:.8f}; M>0 (FD)={o['Mpos']}", flush=True)
