#!/usr/bin/env python3
"""Referee extras (attempted falsifications of side statements).  Reads the cycle data written by v1b.
 (i)   Lemma 2.1 upper bounds: limsup x <= lam/c_1, y <= M_x-type bounds, checked against the cycle extremes.
 (ii)  (A6): the stable manifold of E is not in the attractor: backward orbits along +-v_s leave D_chem / blow up.
 (iii) (A2) uniqueness: backward Riccati solutions from perturbed positive data leave the positive quadrant.
 (iv)  (A7): all sampled points of W^u(E) project strictly inside the projected cycle (polygon test).
 (v)   (A1): orbits from random initial data in D_chem (including far away) end on the same cycle.
usage: v6_extras.py mu rho"""
import sys, math, json
import numpy as np
from scipy.integrate import solve_ivp, dblquad
from matplotlib.path import Path
sys.path.insert(0, '/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03/work/lanes/planar/verify')
from vlog import *

mu = float(sys.argv[1]); rho = float(sys.argv[2])
k = rho * mu * (1 + mu); lam = k / mu
ed = eq_data(mu, k)
o = json.load(open(f"/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03/work/lanes/planar/verify/out_v1b_{mu:g}_{rho:g}.json"))
cyc = newton_cycle(mu, k, o['U_crest'], o['y_crest'])
T = cyc['T']
print(f"== mu={mu} rho={rho} k={k:.5g}  T={T:.5f}")
out = dict(mu=mu, rho=rho)

# cycle in original variables via the log system
prof = cycle_riccati(mu, k, cyc, n=4001)
S = prof['S']
Uc = S[0]; qc = np.exp(S[1]); pc = lam - qc; lyc = S[2]
# (i) Lemma 2.1 upper bounds
c1 = dblquad(lambda b, a: math.exp(-mu * a - b - lam * (a + b)), 0, 1, lambda a: 0, lambda a: 1 - a)[0]
Mx_limsup = lam / c1
lxmax = math.log(k) + Uc.max()
out.update(c1=c1, log_xmax=lxmax, log_limsup_bound=math.log(Mx_limsup), log_ymax=float(lyc.max()), qmax=float(qc.max()))
print(f" (i) Lemma 2.1: x_max on cycle = {math.exp(lxmax):.5g} <= lam/c_1 = {Mx_limsup:.5g}: {math.exp(lxmax) <= Mx_limsup};"
      f"  y_max = {math.exp(lyc.max()):.5g} <= lam/c_1: {math.exp(lyc.max()) <= Mx_limsup};"
      f"  q_max = {qc.max():.5g} <= lam/(c_1 mu) = {Mx_limsup/mu:.5g}: {qc.max() <= Mx_limsup/mu}")
out['lemma21_ok'] = bool(math.exp(lxmax) <= Mx_limsup and math.exp(lyc.max()) <= Mx_limsup and qc.max() <= Mx_limsup / mu)


def f3(t, z):
    U, p, y = z
    return [p, k - mu * p - y, k * ex(U) - y]

# (ii) stable manifold backward
vs = np.array([1.0, -ed['gam'], -ed['gam'] * (ed['gam'] - mu)])
res = []
for sgn in (+1, -1):
    z0 = np.array([0.0, 0.0, k]) + sgn * 1e-6 * vs / np.linalg.norm(vs)

    def leave(t, z):
        return min(z[2], lam - z[1], 1e3 - abs(z[0]))     # leaves D_chem (y<=0 or p>=lam) or |U|>1e3
    leave.terminal = True
    sol = solve_ivp(f3, [0, -400.0], z0, method='DOP853', rtol=1e-11, atol=1e-13, events=leave)
    zf = sol.y[:, -1]
    res.append(dict(sign=sgn, t_exit=float(sol.t[-1]), exited=bool(sol.status == 1), U=float(zf[0]), p=float(zf[1]), y=float(zf[2])))
    print(f" (ii) W^s(E), branch {sgn:+d}: backward orbit {'LEAVES D_chem / |U|>1e3' if sol.status == 1 else 'stays (?)'} at t={sol.t[-1]:.2f}:"
          f" U={zf[0]:.3g} p={zf[1]:.3g} (lam={lam:.3g}) y={zf[2]:.3g}")
out['stable_manifold'] = res

# (iii) backward Riccati from perturbed positive data, along the cycle (x periodic)
from scipy.interpolate import CubicSpline
Ucp = Uc.copy(); Ucp[-1] = Ucp[0]
cs = CubicSpline(prof['t'], Ucp, bc_type='periodic')
eta0, zeta0 = math.exp(S[3, 0]), math.exp(S[4, 0])


def ric(t, s):
    e, z = s
    x = k * math.exp(float(cs(t % T)))
    return [z + (mu - 1) * e - e * e, x - (1 + e) * z]


def neg(t, s):
    return min(s[0], s[1])
neg.terminal = True


def big(t, s):
    return 1e8 - max(abs(s[0]), abs(s[1]))
big.terminal = True
rng = np.random.default_rng(3)
left = 0; rows = []
for trial in range(12):
    fac = np.exp(rng.normal(0, 0.5 if trial < 8 else 0.01, 2))
    s0 = [eta0 * fac[0], zeta0 * fac[1]]
    sol = solve_ivp(ric, [0, -40 * T], s0, method='DOP853', rtol=1e-10, atol=1e-13, events=[neg, big], max_step=T / 20)
    ex_ = (len(sol.t_events[0]) > 0) or (len(sol.t_events[1]) > 0)
    left += int(ex_)
    rows.append(dict(fac=fac.tolist(), exited=bool(ex_), t=float(sol.t[-1]), end=sol.y[:, -1].tolist()))
ref = solve_ivp(ric, [0, -3 * T], [eta0, zeta0], method='DOP853', rtol=1e-10, atol=1e-13, events=[neg, big], max_step=T / 20)
print(f" (iii) backward Riccati from 12 perturbed positive data: {left}/12 leave the closed positive quadrant or blow up"
      f" (times/T: {[round(r['t']/T, 2) for r in rows]}); unperturbed solution followed backward for"
      f" {abs(ref.t[-1])/T:.2f} periods before rounding ejects it")
out['riccati_backward'] = dict(left=left, rows=rows)

# (iv) polygon containment
poly = Path(np.column_stack((Uc, pc)))
t_skip = 35.0 / ed['gam']
r0 = min(1e-3, 1e-3 * math.exp(-ed['lam'].real * t_skip))
tconv, cr = spiral_until_converged(mu, k, r0, 0.3)
# stop sampling when the spiral is within 1e-6 (relative) of the cycle crest
crU = cr[:, 1]
jstop = int(np.argmax(crU > crU[-1] * (1 - 1e-6))) if np.any(crU > crU[-1] * (1 - 1e-6)) else len(crU) - 1
tstop = cr[max(jstop - 1, 0), 0]
sol = run_log(mu, k, spiral_start(mu, k, r0, 0.3), 0.0, tstop, t_eval=np.linspace(t_skip, tstop, 6000))
Pts = np.column_stack((sol.y[0], lam - np.exp(sol.y[1])))
inside = poly.contains_points(Pts)
out['inside_fraction'] = float(inside.mean()); out['n_spiral_pts'] = int(len(Pts))
print(f" (iv) (A7): {inside.sum()}/{len(Pts)} spiral points (up to 1e-6 from the cycle crest) project inside the projected cycle")

# (v) random initial data
rng = np.random.default_rng(11)
ends = []
for trial in range(8):
    X0 = k * math.exp(rng.uniform(-6, 6)); y0 = k * math.exp(rng.uniform(-6, 6)); q0 = lam * math.exp(rng.uniform(-6, 3))
    s0 = np.array([math.log(X0 / k), math.log(q0), math.log(y0), 0.0, 0.0, 0.0])
    tend = 60 * T + 40.0 / max(1e-3, -o['lamc_riccati'])
    crest = make_crest(mu, k)
    sol = run_log(mu, k, s0, 0.0, tend, events=[crest], rtol=1e-11)
    if len(sol.t_events[0]) < 2:
        ends.append(None); continue
    ends.append(float(sol.y_events[0][-1][0]))
out['random_ic_crests'] = ends
dev = max(abs(e - cyc['U0']) for e in ends if e is not None) if any(e is not None for e in ends) else None
print(f" (v) (A1): final crest U of 8 random orbits in D_chem: max deviation from the cycle crest {cyc['U0']:.6f} is {dev}"
      f"  (None entries: {sum(e is None for e in ends)})")
json.dump(out, open(f"/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03/work/lanes/planar/verify/out_v6_{mu:g}_{rho:g}.json", 'w'), indent=1, default=float)
