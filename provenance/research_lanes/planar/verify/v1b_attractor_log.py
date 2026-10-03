#!/usr/bin/env python3
"""Referee battery for Theorem A (A2)-(A8) and Prop 8.1 at one parameter pair, LOG-VARIABLE version
(valid through deep troughs).  usage: v1b_attractor_log.py mu rho [nphase]
Every sign test is evaluated as (sum of signed terms)/(largest |term|), computed without underflow."""
import sys, json, math, time
import numpy as np
import mpmath as mp
sys.path.insert(0, '/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03/work/lanes/planar/verify')
from vlog import *

np.seterr(divide='ignore')
mu = float(sys.argv[1]); rho = float(sys.argv[2])
nphase = int(sys.argv[3]) if len(sys.argv) > 3 else 6
k = rho * mu * (1 + mu); A = 1 + mu; lam = k / mu
lk, ll = math.log(k), math.log(lam)
t_start = time.time()
ed = eq_data(mu, k)
out = dict(mu=mu, rho=rho, k=k, gamE=ed['gam'], etaE=ed['eta'], zetaE=ed['zeta'], lam_u=[ed['lam'].real, ed['lam'].imag])
print(f"== mu={mu} rho={rho} k={k:.6g} gamE={ed['gam']:.6f} lam_u={ed['lam']:.6f}")

# ---------- spirals on W^u(E) ----------
t_skip = 35.0 / ed['gam']
r0 = min(1e-3, 1e-3 * math.exp(-ed['lam'].real * t_skip))
NS = 2700
phases = [2 * math.pi * j / nphase for j in range(nphase)]
tconv, cr = spiral_until_converged(mu, k, r0, phases[0])
tab = cr[cr[:, 0] > t_skip]
nphase_eff = 2 if len(cr) > 400 else nphase
out['spiral_turns'] = int(len(cr)); out['t_conv'] = tconv
print(f"   spiral: {len(cr)} crests, t_conv={tconv:.1f}, r0={r0:.2e}, t_skip={t_skip:.2f}  [{time.time()-t_start:.0f}s]")
pts = []
for j in range(nphase_eff):
    n_j = NS // nphase_eff
    te = np.linspace(t_skip, tconv, 6 * n_j)
    sol = run_log(mu, k, spiral_start(mu, k, r0, phases[j]), 0.0, tconv, t_eval=te)
    assert sol.status == 0 and np.all(np.isfinite(sol.y)), 'spiral integration failed / log-Riccati not finite'
    S = sol.y[:5].T
    pS = lam - np.exp(S[:, 1])

    def arc_idx(Ucoord, n):
        sc = np.array([np.ptp(Ucoord) + 1e-300, np.ptp(pS) + 1e-300])
        d = np.hypot(np.diff(Ucoord) / sc[0], np.diff(pS) / sc[1])
        arc = np.concatenate(([0], np.cumsum(d)))
        return np.clip(np.searchsorted(arc, np.linspace(0, arc[-1], n)), 0, len(arc) - 1)
    idx = np.unique(np.concatenate((np.arange(0, len(te), 18), arc_idx(S[:, 0], n_j // 3),
                                    arc_idx(np.maximum(S[:, 0], -8.0), n_j // 3))))
    pts.append(S[idx])
spiral_pts = np.vstack(pts)
print(f"   spiral samples: {len(spiral_pts)}  [{time.time()-t_start:.0f}s]")

# ---------- the cycle ----------
cyc = newton_cycle(mu, k, tab[-1, 1], math.exp(tab[-1, 2]))
T = cyc['T']
ev = np.linalg.eigvals(cyc['M'])
prof = cycle_riccati(mu, k, cyc)
Sc = prof['S']
lamc_ric = float(Sc[5, -1] / T)                         # <eta> - mu  from  I' = eta - mu
eta_avg = lamc_ric + mu
qr = qr_exponents(mu, k, cyc)
lamc_qr = float(np.sum(qr[-1]) / T)
mt_M = max(e.real for e in ev)
out.update(T=T, U_crest=cyc['U0'], y_crest=cyc['y0'], newton_res=cyc['res'], lamc_riccati=lamc_ric, lamc_qr=lamc_qr,
           lamc_qr_prev=float(np.sum(qr[-2]) / T), lamc_M=(float(math.log(mt_M) / T) if mt_M > 0 else None),
           Umin=float(Sc[0].min()), Umax=float(Sc[0].max()), ric_periods=prof['periods_used'], ric_defect=prof['defect'],
           ric_hist=prof['hist'][:12])
print(f"   cycle: T={T:.6f} crest U={cyc['U0']:.6f} Umin={Sc[0].min():.3f} res={cyc['res']:.1e}; Riccati periodic after"
      f" {prof['periods_used']} periods (defect {prof['defect']:.1e})")
print(f"   A8(a): lambda_c: Riccati {lamc_ric:.10f} | QR, no Riccati {lamc_qr:.10f} (prev period {out['lamc_qr_prev']:.10f})"
      f" | log max eig(M)/T {out['lamc_M'] if out['lamc_M'] is not None else float('nan'):.10f}")
n0 = np.array([math.exp(Sc[4, 0]), -math.exp(Sc[3, 0]), -1.0])
ms_pred = math.exp(-(1 + eta_avg) * T)
out['MTn_relerr'] = float(np.linalg.norm(cyc['Phi'].T @ n0 - ms_pred * n0) / (np.linalg.norm(cyc['Phi']) * np.linalg.norm(n0)))
print(f"   M3^T n(a0) = m_s n(a0): rel.err {out['MTn_relerr']:.2e}   (m_s = {ms_pred:.3e})")

cyc_pts = Sc[:5].T[::3][:-1]
P = np.vstack((spiral_pts, cyc_pts))
U, w, c, a, b = P.T
N = len(P)
p = lam - np.exp(w)
y = np.exp(c)
out['N_points'] = N; out['N_cycle_points'] = int(len(cyc_pts))
one = np.ones(N)
lgabs = lambda z: np.log(np.abs(z))

# ---------- (A5) O = n.f = 0 ; K > 0 ; W > 0 ----------
# O = zeta p - eta (mu q - y) - (x - y)
O, _ = nsum([(np.sign(p), b + lgabs(p)), (-one, a + math.log(mu) + w), (one, a + c), (-one, lk + U), (one, c)])
Kn, _ = nsum([(np.sign(U), b + lgabs(U)), (-one, a + ll), (one, a + w), (-one, c), (one, lk + 0 * U)])
terms = [(one, lk + U + ll), (-one, lk + U + w), (-one, b + ll), (one, b + w), (-one, a + b + ll), (one, a + b + w),
         (-one, math.log(mu) + b + w), (one, b + c), (one, math.log(mu) + 2 * a + w), (-one, 2 * a + c)]
if mu != 1.0:
    sg = math.copysign(1.0, mu - 1)
    terms += [(-sg * one, math.log(abs(mu - 1) * mu) + a + w), (sg * one, math.log(abs(mu - 1)) + a + c)]
Wn, _ = nsum(terms)
LamE = (y - k) - ed['zeta'] * U + ed['eta'] * p
Lsc = np.abs(y - k) + np.abs(ed['zeta'] * U) + np.abs(ed['eta'] * p) + 1e-300
inCp = (U >= 0) & (p <= 0) & (y - k <= 0)
inCm = (U <= 0) & (p >= 0) & (y - k >= 0)
out.update(O_rel_max=float(np.max(np.abs(O))), K_rel_min=float(Kn.min()), W_rel_min=float(Wn.min()),
           LamE_rel_min=float(np.min(LamE / Lsc)), n_in_Cplus=int(inCp.sum()), n_in_Cminus=int(inCm.sum()))
print(f"   A5: max |O|/maxterm = {out['O_rel_max']:.2e};  min K/maxterm = {out['K_rel_min']:.3e};  min W/maxterm = {out['W_rel_min']:.3e}")
print(f"   A6: min Lambda_E/scale = {out['LamE_rel_min']:.3e};  points in E+C_+: {inCp.sum()}, in E+C_-: {inCm.sum()}")

# ---------- (A2) bounds (mpmath, since x_- underflows) ----------
mp.mp.dps = 40
xm = mp.mpf(k) * mp.e ** mp.mpf(float(Sc[0].min())); xM = mp.mpf(k) * mp.e ** mp.mpf(float(Sc[0].max()))
m_ = mp.mpf(mu)
eta_p = (m_ - 1 + mp.sqrt((m_ - 1) ** 2 + 4 * xM)) / 2
zeta_m = xm / (1 + eta_p)
eta_m = (m_ - 1 + mp.sqrt((m_ - 1) ** 2 + 4 * zeta_m)) / 2
if eta_m <= max(0, mu - 1):                      # cancellation guard: series
    eta_m = (zeta_m / (1 - m_)) if mu < 1 else mp.sqrt(zeta_m) if mu == 1 else (m_ - 1) + zeta_m / (m_ - 1)
zeta_p = xM / (1 + eta_m)
lb = [float(mp.log(eta_m)), float(mp.log(eta_p)), float(mp.log(zeta_m)), float(mp.log(zeta_p))]
ok = bool(a.min() >= lb[0] - 1e-9 and a.max() <= lb[1] + 1e-9 and b.min() >= lb[2] - 1e-9 and b.max() <= lb[3] + 1e-9)
out.update(bounds_log=dict(eta_m=lb[0], eta_p=lb[1], zeta_m=lb[2], zeta_p=lb[3], a_min=float(a.min()), a_max=float(a.max()),
                           b_min=float(b.min()), b_max=float(b.max())), bounds_ok=ok,
           eta_minus_max0mu1_min=float(np.min(np.exp(a) - max(0.0, mu - 1))))
print(f"   A2 bounds (logs): log eta in [{a.min():.4g},{a.max():.4g}] vs [{lb[0]:.4g},{lb[1]:.4g}];"
      f" log zeta in [{b.min():.4g},{b.max():.4g}] vs [{lb[2]:.4g},{lb[3]:.4g}]  ok={ok};"
      f" min(eta - max(0,mu-1)) = {out['eta_minus_max0mu1_min']:.3e}")

# ---------- (A3)/(A4): supporting planes, ALL ordered pairs, underflow-free ----------
worst = -np.inf; worst_pair = None; n_pos = 0; n_gt = 0; n_sig = 0; margin = -np.inf; n_close = 0; far_worst = -np.inf
hist_edges = np.array([-3.1, -1e-3, -1e-6, -1e-9, -1e-12, 0.0, 1e-12, 1e-9, 1e-6, 3.1])
hist = np.zeros(len(hist_edges) - 1, dtype=np.int64)
for i0 in range(0, N, 300):
    sl = slice(i0, min(N, i0 + 300))
    dU = U[None, :] - U[sl, None]
    dw = w[None, :] - w[sl, None]
    dc = c[None, :] - c[sl, None]
    with np.errstate(divide='ignore', invalid='ignore'):
        L1 = b[sl, None] + np.log(np.abs(dU))                                         # zeta_a (U_b-U_a)
        L2 = a[sl, None] + np.maximum(w[None, :], w[sl, None]) + np.log(-np.expm1(-np.abs(dw)))   # eta_a (q_b-q_a)
        L3 = np.maximum(c[None, :], c[sl, None]) + np.log(-np.expm1(-np.abs(dc)))    # -(y_b-y_a)
    M = np.maximum(np.maximum(L1, L2), L3)
    with np.errstate(invalid='ignore'):
        R = np.sign(dU) * np.exp(L1 - M) + np.sign(dw) * np.exp(L2 - M) - np.sign(dc) * np.exp(L3 - M)
    ii = np.arange(sl.start, sl.stop)
    R[ii - sl.start, ii] = -1.0
    R[~np.isfinite(R)] = 9.0             # would be a coincident pair of distinct samples: flag
    # rounding amplification: coordinates carry abs. error ~1e-12(1+|coord|); differences of near-coincident
    # points lose digits.  delta = relative separation;  a violation is 'significant' if R > 10*(1e-11/delta+1e-11)
    delta = np.maximum(np.maximum(np.abs(dU) / (1 + np.abs(U[sl, None])), np.abs(dw) / (1 + np.abs(w[sl, None]))),
                       np.abs(dc) / (1 + np.abs(c[sl, None])))
    delta[ii - sl.start, ii] = 1.0
    with np.errstate(divide='ignore'):
        thr = 10 * (1e-11 / delta + 1e-11)
    n_sig += int((R > thr).sum())
    margin = max(margin, float(np.max(R - thr)))
    close = delta < 1e-6
    n_close += int(close.sum()) - 0
    far_worst = max(far_worst, float(np.max(np.where(close, -9.0, R))))
    j = np.unravel_index(np.argmax(R), R.shape)
    if R[j] > worst:
        worst = float(R[j]); worst_pair = (int(j[0] + sl.start), int(j[1]))
    n_pos += int((R >= 0).sum()); n_gt += int((R > 1e-9).sum())
    hist += np.histogram(R, bins=hist_edges)[0]
out.update(support_pairs=int(N * (N - 1)), support_worst_rel=worst, support_nonneg=n_pos, support_viol_gt_1e9=n_gt,
           support_hist=dict(edges=hist_edges.tolist(), counts=hist.tolist()))
ia, ib = worst_pair
out['support_worst_pair'] = dict(a=P[ia].tolist(), b=P[ib].tolist(), ia=ia, ib=ib, n_spiral=int(len(spiral_pts)))
out.update(support_significant=n_sig, support_margin=margin, support_n_close=n_close, support_worst_rel_far=far_worst)
print(f"   A3: {N*(N-1)} ordered pairs; worst Delta/maxterm = {worst:.3e}; #(>=0) = {n_pos}; #(>1e-9) = {n_gt};"
      f" SIGNIFICANT (beyond rounding amplification) = {n_sig}; worst among pairs with rel. separation >1e-6 = {far_worst:.3e}"
      f" ({n_close} closer pairs)")
print(f"       histogram of Delta/maxterm over {hist_edges.tolist()}: {hist.tolist()}")

# ---------- crest sequences: Remark 5.1 formula, Prop 8.1(ii), invariance of the section tangent (1,zeta) ----------
uC = tab[:, 1]; yC = np.exp(tab[:, 2]); eC = np.exp(tab[:, 3]); zC = np.exp(tab[:, 4]); IC = tab[:, 5]
gt = yC - k
out['gtilde_formula_relerr'] = float(np.max(np.abs(gt - k * np.expm1(uC) / (1 + eC)) / np.abs(gt)))
dcr = np.diff(uC)
out['crests_increasing'] = bool(np.all(dcr > -1e-11 * abs(uC[-1])))
out['crest_first_last'] = [float(uC[0]), float(uC[-1])]
sel = np.unique(np.linspace(0, len(tab) - 2, 8).astype(int))
rows = []
for j in sel:
    iota = IC[j + 1] - IC[j]
    pred = gt[j] / gt[j + 1] * math.exp(iota)
    try:
        U1, Y1, tau, M, Phi = crest_map(mu, k, uC[j], yC[j])
        wv = M @ np.array([1.0, zC[j]])
        var = float(wv[0]); par = float((wv[1] - zC[j + 1] * wv[0]) / (abs(wv[1]) + abs(zC[j + 1] * wv[0]) + 1e-300))
        noise = float(1e-12 * np.abs(Phi).max() * (1 + zC[j]))
    except Exception as e:
        var = par = noise = float('nan')
    rows.append(dict(u=float(uC[j]), R1=float(uC[j + 1]), iota=float(iota), pred=float(pred), var=var,
                     rel=float(abs(var - pred) / abs(pred)), parallel_defect=par, noise_floor=noise))
out['prop81'] = rows
print("   Prop 8.1(ii): R1'(u) predicted g(u)e^iota/g(R1 u)  vs  U-component of M(1,zeta(u)) [variational, 3D];"
      " parallel defect of M(1,zeta) vs (1,zeta(R1 u))")
for r in rows:
    print(f"      u={r['u']:.6g} R1={r['R1']:.6g} iota={r['iota']:+.5f} pred={r['pred']:.6e} var={r['var']:.6e}"
          f" rel={r['rel']:.1e} par={r['parallel_defect']:.1e} (abs noise ~{r['noise_floor']:.0e})")
print(f"   Remark 5.1: h(u,0)-k = k(e^u-1)/(1+eta): rel.err {out['gtilde_formula_relerr']:.2e};"
      f" crests increasing: {out['crests_increasing']} ({uC[0]:.4g} -> {uC[-1]:.6g})")

# finite differences between neighbouring orbits of W^u(E) (no variational equation, no zeta)
dth = 2e-3
_, cr_a = spiral_until_converged(mu, k, r0, phases[0] + dth)
_, cr_b = spiral_until_converged(mu, k, r0, phases[0] - dth)
cr_a = cr_a[cr_a[:, 0] > t_skip]; cr_b = cr_b[cr_b[:, 0] > t_skip]
m = min(len(cr_a), len(cr_b), len(tab))
fd_rows = []
for j in sel:
    if j + 1 >= m:
        continue
    du0 = cr_a[j, 1] - cr_b[j, 1]; du1 = cr_a[j + 1, 1] - cr_b[j + 1, 1]
    if abs(du0) < 1e-8 * abs(uC[j]) or abs(du1) < 3e-11:
        continue
    iota = IC[j + 1] - IC[j]
    pred = gt[j] / gt[j + 1] * math.exp(iota)
    fd_rows.append(dict(u=float(uC[j]), fd=float(du1 / du0), pred=float(pred), rel=float(abs(du1 / du0 - pred) / abs(pred)),
                        du0=float(du0), du1=float(du1)))
out['prop81_fd'] = fd_rows
for r in fd_rows:
    print(f"      FD: u={r['u']:.6g} du0={r['du0']:.2e} du1={r['du1']:.2e} fd={r['fd']:.6e} pred={r['pred']:.6e} rel={r['rel']:.1e}")

# ---------- (A2) pull-back on the cycle from extreme positive data ----------
gaps = []
for data in ((1e-8, 1e-8), (1e8, 1e8), (1e-8, 1e8), (1e8, 1e-8)):
    res = cycle_riccati(mu, k, cyc, a0=math.log(data[0]), b0=math.log(data[1]), max_periods=400)
    gap = max(abs(res['S'][3, 0] - Sc[3, 0]), abs(res['S'][4, 0] - Sc[4, 0]))
    gaps.append(dict(data=data, periods=res['periods_used'], final_gap=float(gap), hist=res['hist'][:6]))
    print(f"   A2 pull-back on the cycle from (eta,zeta)={data}: periodic after {res['periods_used']} periods,"
          f" log-gap to reference {gap:.1e}; per-period log-changes {['%.1e' % h for h in res['hist'][:5]]}")
out['pullback'] = gaps

out['seconds'] = time.time() - t_start
fn = f"/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03/work/lanes/planar/verify/out_v1b_{mu:g}_{rho:g}.json"
json.dump(out, open(fn, 'w'), indent=1, default=float)
print(f"   done in {out['seconds']:.0f}s -> {fn}")
