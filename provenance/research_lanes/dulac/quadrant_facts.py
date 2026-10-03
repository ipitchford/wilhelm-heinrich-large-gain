#!/usr/bin/env python3
"""Numerical status of the quadrant sign facts on cycles (20 fresh sets) and on the whole disc (10 LP sets):
  F1: eta < mu on QI (p>0,p'>0);  F2: eta > mu on QIII (p<0,p'<0);  F3: eta(crest) > mu > eta(trough);
  F4: eta0 = V_UH/V > 0 (third component of r^z' positive);  F5: sup K/V, sup W/Q."""
import numpy as np
from verify_theorem import cycle_log
from gen_data import eqdata
from lpkit import load, SETS

def facts(mu, U, p, pd, H, x, y, eta, K, W):
    V = p * p - U * pd; Q = pd * (-H) + p * (x - y); VUH = U * (x - y) - p * H
    q = eta - mu
    QI = (p > 0) & (pd > 0); QIII = (p < 0) & (pd < 0)
    return dict(F1=q[QI].max(), F2=q[QIII].min(), F4=(VUH / V).min(), kap=(K / V).max(), var=(W / Q).max())

print('--- cycles ---')
for mu in (0.1, 0.3, 0.7, 1.0):
    guess = None
    for rho in (1.001, 1.1, 2.5, 8.0, 20.0):
        k = rho * mu * (1 + mu)
        if rho < 1.01:
            ed = eqdata(mu, k); A = 1 + mu; D = mu * mu + 6 * mu + 1
            eps = np.sqrt(4 * ed['sigma'] * (A * A + mu) * D / (A * mu * mu)); guess = (eps, k + mu * eps)
        c = cycle_log(mu, k, guess=guess); guess = c['crest']
        f = facts(mu, c['U'], c['p'], c['pd'], c['H'], c['x'], c['y'], c['eta'], c['K'], c['W'])
        i_c = np.argmax(c['U']); i_t = np.argmin(c['U'])
        print(f'mu={mu:<4} rho={rho:<6} max q on QI={f["F1"]:+.3e}  min q on QIII={f["F2"]:+.3e}  q(crest)={c["eta"][i_c]-mu:+.3e} q(trough)={c["eta"][i_t]-mu:+.3e} '
              f'min eta0={f["F4"]:+.3e}  max K/V={f["kap"]:.3f} max W/Q={f["var"]:.3f}')
print('--- whole disc (W^u(E) samples) ---')
for s in SETS:
    u = load(*s)
    f = facts(u['mu'], u['U'], u['p'], u['pd'], u['H'], u['x'], u['y'], u['eta'], u['K'], u['W'])
    print(f'mu={s[0]:<5} rho={s[1]:<5} max q on QI={f["F1"]:+.3e}  min q on QIII={f["F2"]:+.3e}  min eta0={f["F4"]:+.3e}  max K/V={f["kap"]:.3f} max W/Q={f["var"]:.3f}  (2 sigma={2*u["sig"]:.3e})')
