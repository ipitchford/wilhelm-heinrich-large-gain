#!/usr/bin/env python3
"""Referee test battery for Theorem A (A2)-(A8) and Prop 8.1 at one parameter pair.
usage: v1_attractor.py mu rho [nphase]
Writes verify/out_v1_<mu>_<rho>.json and prints a summary.  Attempts to FALSIFY; a pass is not a proof."""
import sys, json, math, time
import numpy as np
from scipy.linalg import expm
sys.path.insert(0, '/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03/work/lanes/planar/verify')
from vlib import *

mu = float(sys.argv[1]); rho = float(sys.argv[2])
nphase = int(sys.argv[3]) if len(sys.argv) > 3 else 6
k = rho * mu * (1 + mu); A = 1 + mu; lam = k / mu
t_start = time.time()
ed = eq_data(mu, k)
out = dict(mu=mu, rho=rho, k=k, gamE=ed['gam'], etaE=ed['eta'], zetaE=ed['zeta'],
           lam_u=[ed['lam'].real, ed['lam'].imag])
print(f"== mu={mu} rho={rho} k={k:.6g} gamE={ed['gam']:.6f} lam_u={ed['lam']:.6f}")

# ---------- spirals on W^u(E) ----------
t_skip = 35.0 / ed['gam']
r0 = min(1e-3, 1e-3 * math.exp(-ed['lam'].real * t_skip))
per0 = 2 * math.pi / ed['lam'].imag
NS = 2600                       # samples per run in total
pts = []                        # rows (U,p,y,eta,zeta)
crest_tabs = []
phases = [2 * math.pi * j / nphase for j in range(nphase)]
tconv0 = None
nphase_eff = nphase
for j, th in enumerate(phases):
    if j >= nphase_eff:
        break
    if tconv0 is None:
        tconv, cr = spiral_until_converged(mu, k, r0, th)
        tconv0 = tconv
        out['spiral_turns'] = int(len(cr)); out['t_conv'] = tconv
        print(f"   spiral: {len(cr)} crests, t_conv={tconv:.1f}, r0={r0:.2e}, t_skip={t_skip:.2f}  [{time.time()-t_start:.0f}s]")
        if len(cr) > 400:        # slow spiral: one phase is already dense
            nphase_eff = 2
        else:
            nphase_eff = nphase
    if j == 0:
        crest_tabs.append(cr[cr[:, 0] > t_skip])
    n_j = NS // nphase_eff
    te = np.linspace(t_skip, tconv, 4 * n_j)
    sol = run_aug(mu, k, spiral_start(mu, k, r0, th), 0.0, tconv, t_eval=te)
    S = sol.y[:5].T
    # half uniform in time, half uniform in planar arclength (scaled)
    sc = np.array([np.ptp(S[:, 0]) + 1e-300, np.ptp(S[:, 1]) + 1e-300])
    d = np.hypot(*(np.diff(S[:, :2], axis=0) / sc).T)
    arc = np.concatenate(([0], np.cumsum(d)))
    idx_t = np.arange(0, len(te), 8)
    idx_a = np.unique(np.searchsorted(arc, np.linspace(0, arc[-1], n_j // 2)))
    idx = np.unique(np.concatenate((idx_t, np.clip(idx_a, 0, len(te) - 1))))
    pts.append(S[idx])
spiral_pts = np.vstack(pts)
print(f"   spiral samples: {len(spiral_pts)}  [{time.time()-t_start:.0f}s]")

# ---------- the cycle (own Newton) ----------
cr = crest_tabs[0]
cyc = newton_cycle(mu, k, cr[-1, 1], cr[-1, 2])
ev = np.linalg.eigvals(cyc['M'])
prof = cycle_riccati(mu, k, cyc)
T = cyc['T']
eta_avg = float(np.trapezoid(prof['eta'], prof['t']) / T)
lamc_ric = eta_avg - mu
lamc_I = float(prof['I'][-1] / T)
qr = qr_exponents(mu, k, cyc, nsub=300, periods=3)
lamc_qr = float(np.sum(qr[-1]) / T)
out.update(T=T, U_crest=cyc['U0'], y_crest=cyc['y0'], newton_res=cyc['res'], M_eigs=[complex(e).real for e in ev],
           lamc_riccati=lamc_ric, lamc_I=lamc_I, lamc_qr=lamc_qr, lamc_qr_prev=float(np.sum(qr[-2]) / T),
           Umin=float(prof['U'].min()), Umax=float(prof['U'].max()), ric_periods=prof['periods_used'],
           ric_defect=prof['defect'])
mt_M = max(e.real for e in ev)
out['logmt_M'] = float(math.log(mt_M)) if mt_M > 0 else None
out['logmt_ric'] = lamc_ric * T
print(f"   cycle: T={T:.6f} crest U={cyc['U0']:.6f} Umin={prof['U'].min():.3f} res={cyc['res']:.1e}")
print(f"   A8(a): lambda_c  Riccati <eta>-mu = {lamc_ric:.10f} | integral I/T = {lamc_I:.10f} | QR (no Riccati) = {lamc_qr:.10f}"
      f" | log(max eig M)/T = {out['logmt_M']/T if out['logmt_M'] is not None else float('nan'):.10f}")
# left eigenvector test  M3^T n = m_s n  and  n.f = 0 on the cycle
n0 = np.array([prof['zeta'][0], -prof['eta'][0], -1.0])
ms_pred = math.exp(-(1 + eta_avg) * T)
lhs = cyc['Phi'].T @ n0
out['MTn_relerr'] = float(np.linalg.norm(lhs - ms_pred * n0) / (np.linalg.norm(cyc['Phi']) * np.linalg.norm(n0)))
xc = k * np.exp(prof['U'])
O_cyc = prof['zeta'] * prof['p'] - prof['eta'] * (k - mu * prof['p'] - prof['y']) - (xc - prof['y'])
O_sc = np.abs(prof['zeta'] * prof['p']) + np.abs(prof['eta'] * (k - mu * prof['p'] - prof['y'])) + np.abs(xc - prof['y']) + 1e-300
out['O_cycle_rel'] = float(np.max(np.abs(O_cyc) / O_sc))
print(f"   M^T n = m_s n rel.err {out['MTn_relerr']:.2e};  n.f on cycle rel {out['O_cycle_rel']:.2e}")

cyc_pts = np.column_stack([prof[q] for q in ('U', 'p', 'y', 'eta', 'zeta')])[::4][:-1]
P = np.vstack((spiral_pts, cyc_pts))
U, p, y, eta, zeta = P.T
x = k * np.exp(U)
N = len(P)
out['N_points'] = N

# ---------- (A5): O = n.f = 0, K>0, W>0 ; (A6) cones, Lambda_E ----------
pp = k - mu * p - y
O = zeta * p - eta * pp - (x - y)
Osc = np.abs(zeta * p) + np.abs(eta * pp) + np.abs(x - y) + 1e-300
K = zeta * U - eta * p - y + k
etap = zeta + (mu - 1) * eta - eta ** 2
zetap = x - (1 + eta) * zeta
W = zetap * p - etap * pp
LamE = (y - k) - ed['zeta'] * U + ed['eta'] * p
rad = np.hypot(U / (np.ptp(U)), p / np.ptp(p))
Ksc = np.abs(zeta * U) + np.abs(eta * p) + np.abs(y - k) + 1e-300
Wsc = np.abs(zetap * p) + np.abs(etap * pp) + 1e-300
Lsc = np.abs(y - k) + np.abs(ed['zeta'] * U) + np.abs(ed['eta'] * p) + 1e-300
inCp = (U >= 0) & (p <= 0) & (y - k <= 0)
inCm = (U <= 0) & (p >= 0) & (y - k >= 0)
out.update(O_rel_max=float(np.max(np.abs(O) / Osc)), K_rel_min=float(np.min(K / Ksc)), W_rel_min=float(np.min(W / Wsc)),
           LamE_rel_min=float(np.min(LamE / Lsc)), n_in_Cplus=int(inCp.sum()), n_in_Cminus=int(inCm.sum()),
           K_min=float(K.min()), W_min=float(W.min()), LamE_min=float(LamE.min()))
print(f"   A5: max|O|rel={out['O_rel_max']:.2e}  min K/sc={out['K_rel_min']:.3e}  min W/sc={out['W_rel_min']:.3e}"
      f"  A6: min LamE/sc={out['LamE_rel_min']:.3e}  #C+={inCp.sum()} #C-={inCm.sum()}")

# ---------- (A2) bounds with x_-/x_+ over the attractor (= over the cycle, by (A7)) ----------
xm, xM = k * math.exp(prof['U'].min()), k * math.exp(prof['U'].max())
eta_p = (mu - 1 + math.sqrt((mu - 1) ** 2 + 4 * xM)) / 2
zeta_m = xm / (1 + eta_p)
eta_m = (mu - 1 + math.sqrt((mu - 1) ** 2 + 4 * zeta_m)) / 2
zeta_p = xM / (1 + eta_m)
out.update(bounds=dict(eta_m=eta_m, eta_p=eta_p, zeta_m=zeta_m, zeta_p=zeta_p,
                       eta_min=float(eta.min()), eta_max=float(eta.max()), zeta_min=float(zeta.min()), zeta_max=float(zeta.max())),
           bounds_ok=bool(eta.min() >= eta_m * (1 - 1e-9) and eta.max() <= eta_p * (1 + 1e-9)
                          and zeta.min() >= zeta_m * (1 - 1e-9) and zeta.max() <= zeta_p * (1 + 1e-9)),
           eta_gt_mu_minus_1=bool(eta.min() > max(0.0, mu - 1)))
print(f"   A2 bounds: eta in [{eta.min():.4g},{eta.max():.4g}] vs [{eta_m:.4g},{eta_p:.4g}];"
      f" zeta in [{zeta.min():.4g},{zeta.max():.4g}] vs [{zeta_m:.4g},{zeta_p:.4g}]  ok={out['bounds_ok']}")

# ---------- (A3)/(A4): supporting planes for ALL ordered pairs ----------
viol = 0; worst = -np.inf; worst_pair = None; nstrict = 0; worst_abs = -np.inf
close_worst = -np.inf
dq_sc = np.array([np.ptp(U), np.ptp(p)])
for i0 in range(0, N, 400):
    sl = slice(i0, min(N, i0 + 400))
    dU = U[None, :] - U[sl, None]; dp = p[None, :] - p[sl, None]; dy = y[None, :] - y[sl, None]
    t1 = zeta[sl, None] * dU; t2 = -eta[sl, None] * dp
    D = t1 + t2 - dy                      # must be < 0 for a != b
    sc = np.abs(t1) + np.abs(t2) + np.abs(dy)
    sc[sc == 0] = 1.0
    R = D / sc
    ii = np.arange(sl.start, sl.stop)
    R[ii - sl.start, ii] = -1.0
    D[ii - sl.start, ii] = -np.inf
    j = np.unravel_index(np.argmax(R), R.shape)
    if R[j] > worst:
        worst = float(R[j]); worst_pair = (int(j[0] + sl.start), int(j[1]))
    worst_abs = max(worst_abs, float(D.max()))
    viol += int((R > 1e-9).sum())
    nstrict += int((R >= 0).sum())
out.update(support_pairs=int(N * (N - 1)), support_worst_rel=worst, support_worst_abs=worst_abs,
           support_viol_gt_1e9=viol, support_nonneg=nstrict)
a, b = worst_pair
out['support_worst_pair'] = dict(a=P[a].tolist(), b=P[b].tolist())
print(f"   A3: {N*(N-1)} ordered pairs; worst Delta/scale = {worst:.3e} (abs max {worst_abs:.3e});"
      f" #(rel>1e-9) = {viol}; #(rel>=0) = {nstrict}")

# ---------- crest sequences: Prop 8.1 and M (1,zeta) || (1,zeta') ----------
tab = crest_tabs[0]
gt = tab[:, 2] - k
g_formula = k * (np.exp(tab[:, 1]) - 1) / (1 + tab[:, 3])
out['gtilde_formula_relerr'] = float(np.max(np.abs(gt - g_formula) / np.abs(gt)))
dcr = np.diff(tab[:, 1])
mono = bool(np.all(dcr > -1e-11 * abs(tab[-1, 1])))
out['crests_increasing'] = mono
out['crest_max_decrease'] = float(max(0.0, -dcr.min()))
out['crest_first_last'] = [float(tab[0, 1]), float(tab[-1, 1])]
# choose up to 8 crests spread over the spiral
sel = np.unique(np.linspace(0, len(tab) - 2, 8).astype(int))
rows = []
for j in sel:
    u0, y0, e0, z0, I0 = tab[j, 1], tab[j, 2], tab[j, 3], tab[j, 4], tab[j, 5]
    u1, y1, e1, z1, I1 = tab[j + 1, 1], tab[j + 1, 2], tab[j + 1, 3], tab[j + 1, 4], tab[j + 1, 5]
    iota = I1 - I0
    pred = (y0 - k) / (y1 - k) * math.exp(iota)
    try:
        U1, Y1, tau, M, Phi = crest_map(mu, k, u0, y0)
        w = M @ np.array([1.0, z0])
        var = w[0]
        par = (w[1] - z1 * w[0]) / (abs(w[1]) + abs(z1 * w[0]) + 1e-300)
        Mnorm = float(np.abs(M).max())
        ret_err = abs(U1 - u1)
    except Exception as e:
        var = float('nan'); par = float('nan'); Mnorm = float('nan'); ret_err = float('nan')
    rows.append(dict(u=u0, R1=u1, iota=iota, pred=pred, var=float(var), rel=float(abs(var - pred) / abs(pred)),
                     parallel_defect=float(par), Mmax=Mnorm, ret_err=float(ret_err)))
out['prop81'] = rows
print("   Prop 8.1:  u -> R1(u);  R1' predicted g(u)e^iota/g(R1) vs variational (M(1,zeta))_U ; parallel defect")
for r in rows:
    print(f"      u={r['u']:.6g} R1={r['R1']:.6g} iota={r['iota']:+.5f} pred={r['pred']:.6e} var={r['var']:.6e}"
          f" rel={r['rel']:.1e} par={r['parallel_defect']:.1e}")
print(f"   gtilde formula rel.err {out['gtilde_formula_relerr']:.2e}; crests increasing: {mono}")

# finite-difference version of R1' between two neighbouring orbits of W^u(E) (no variational equation, no zeta)
if len(tab) >= 3:
    dth = 2e-3
    tc_a, cr_a = spiral_until_converged(mu, k, r0, phases[0] + dth)
    tc_b, cr_b = spiral_until_converged(mu, k, r0, phases[0] - dth)
    cr_a = cr_a[cr_a[:, 0] > t_skip]; cr_b = cr_b[cr_b[:, 0] > t_skip]
    m = min(len(cr_a), len(cr_b), len(tab))
    fd_rows = []
    for j in sel:
        if j + 1 >= m:
            continue
        du0 = cr_a[j, 1] - cr_b[j, 1]; du1 = cr_a[j + 1, 1] - cr_b[j + 1, 1]
        if abs(du0) < 1e-7 * abs(tab[j, 1]) or abs(du1) < 1e-10:
            continue
        fd = du1 / du0
        iota = tab[j + 1, 5] - tab[j, 5]
        pred = (tab[j, 2] - k) / (tab[j + 1, 2] - k) * math.exp(iota)
        fd_rows.append(dict(u=float(tab[j, 1]), fd=float(fd), pred=float(pred), rel=float(abs(fd - pred) / abs(pred)), du0=float(du0)))
    out['prop81_fd'] = fd_rows
    for r in fd_rows:
        print(f"      FD: u={r['u']:.6g} du0={r['du0']:.2e} fd={r['fd']:.6e} pred={r['pred']:.6e} rel={r['rel']:.1e}")

# ---------- (A2): pull-back from extreme data on the cycle, and the stated rate ----------
Bm = np.array([[0, 0, xm], [1, mu, 0], [0, 1, 1.0]]); Bp = np.array([[0, 0, xM], [1, mu, 0], [0, 1, 1.0]])
m0 = expm(Bm).min(); M0 = expm(Bp).max()
Delta0 = 2 * math.log(M0 / m0) if m0 > 0 else float('inf')
tau0 = math.tanh(Delta0 / 4)
out.update(m0=float(m0), M0=float(M0), Delta0=Delta0, one_minus_tau0=float(1 - tau0))
gaps = []
for data in ((0.0, 0.0), (1e6, 1e6), (0.0, 1e6), (1e6, 0.0)):
    s = np.array([cyc['U0'], 0.0, cyc['y0'], data[0], data[1], 0.0])
    tt = np.array([1.0, 2.0, 5.0, 10.0, 20.0, 40.0])
    # reference periodic solution evaluated at the same times
    sr = np.array([cyc['U0'], 0.0, cyc['y0'], prof['eta'][0], prof['zeta'][0], 0.0])
    try:
        a_ = run_aug(mu, k, s, 0.0, tt[-1], t_eval=tt)
        b_ = run_aug(mu, k, sr, 0.0, tt[-1], t_eval=tt)
        if a_.y.shape[1] < len(tt):
            gaps.append(dict(data=data, failed=True)); continue
        g = np.maximum(np.abs(np.log(a_.y[3] / b_.y[3])), np.abs(np.log(a_.y[4] / b_.y[4])))
        gaps.append(dict(data=data, t=tt.tolist(), gap=g.tolist(), minvals=[float(a_.y[3].min()), float(a_.y[4].min())]))
    except Exception as e:
        gaps.append(dict(data=data, failed=str(e)))
out['pullback'] = gaps
for g in gaps:
    if 'gap' in g:
        print(f"   A2 pull-back from {g['data']}: gap at t=1,2,5,10,20,40: " + ' '.join(f"{v:.1e}" for v in g['gap']))
    else:
        print(f"   A2 pull-back from {g['data']}: FAILED {g}")
print(f"   A2 rate constants: m0={m0:.3e} M0={M0:.3e} Delta0={Delta0:.3f} 1-tau0={1-tau0:.3e}")

out['seconds'] = time.time() - t_start
fn = f"/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03/work/lanes/planar/verify/out_v1_{mu:g}_{rho:g}.json"
json.dump(out, open(fn, 'w'), indent=1, default=float)
print(f"   done in {out['seconds']:.0f}s -> {fn}")
