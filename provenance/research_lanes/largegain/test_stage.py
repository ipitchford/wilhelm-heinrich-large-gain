#!/usr/bin/env python3
"""Direct test of the Stage-1 / Stage-2 lemmas (Lemmas 11 and 14 of theorem.tex) on random
NON-periodic negative-lobe segments.  The proofs of those lemmas use only:
   sign structure of the lobe, |w| < k, D > ell, and Lambda_E(s_w) > 0  (i.e. P_w < (gamma_E - mu) v_w).
So for ANY orbit segment with these properties and k >= 1000:
   T_- >= 5.89  ==>  p_k >= 0.374 k                (Stage 1)
   T_- >= 22    ==>  r_ell >= 0.3 T_-              (Stage 2)
Segments violating a hypothesis are skipped and counted."""
import sys, numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

ELL = 2.0
def rhs(mu, k):
    def f(t, s):
        U, p, w = s
        return [p, -mu*p-w, -w+k*np.expm1(min(U, 50.0))]
    return f
def upx(t, s): return s[0]
upx.direction = 1; upx.terminal = True

def one(mu, k, Pd, wd):
    A = 1+mu
    g = brentq(lambda g: g*(g-1)*(g-mu)-k, A, k**(1/3)+A+1)
    sol = solve_ivp(rhs(mu, k), [0, 3000.0], [-1e-12, -Pd, wd], method='DOP853', rtol=1e-10, atol=1e-12,
                    events=upx, dense_output=True)
    if sol.status != 1: return 'no upcross', None
    Tm = sol.t_events[0][0]; S = sol.sol
    tt = np.unique(np.concatenate([sol.t, np.linspace(0, Tm, 20001)])); Y = S(tt)
    Ug, pg, wg = Y
    if np.any(wg <= -k): return 'y<0', None
    iw = np.where(wg < 0)[0]; ip = np.where(pg > 0)[0]
    if len(iw) == 0 or iw[0] == 0 or len(ip) == 0: return 'degenerate', None
    s_w = brentq(lambda t: S(t)[2], tt[iw[0]-1], tt[iw[0]])
    s_0 = brentq(lambda t: S(t)[1], tt[ip[0]-1], tt[ip[0]])
    if np.any(pg[ip[0]:] < 0): return 'p not unimodal', None
    D = -S(s_0)[0]; pk = S(Tm)[1]
    if D <= ELL: return 'D<=ell', None
    Pw = -S(s_w)[1]; vw = -S(s_w)[0]
    if not (Pw < (g-mu)*vw): return 'LambdaE(s_w)<=0', None
    ir = np.where((tt > s_0) & (Ug > -ELL))[0][0]
    s_l = brentq(lambda t: S(t)[0]+ELL, tt[ir-1], tt[ir])
    return 'ok', dict(Tm=Tm, pk=pk, r_l=s_l-s_0, s_w=s_w, u=s_0-s_w, tau_b=Tm-s_l, D=D)

if __name__ == '__main__':
    rng = np.random.default_rng(7)
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 400
    out = {}; n1 = n2 = f1 = f2 = 0; worst1 = np.inf; worst2 = np.inf; maxtb = 0
    for i in range(n):
        mu = float(rng.choice([1.0, 0.6, 0.25]))
        k = float(np.exp(rng.uniform(np.log(1000), np.log(6000))))
        wd = k*np.exp(rng.uniform(np.log(0.02), np.log(3e4)))
        Pd = wd/k**(1/3)*np.exp(rng.uniform(np.log(1e-3), np.log(0.9)))
        r, d = one(mu, k, Pd, wd)
        out[r] = out.get(r, 0)+1
        if d is None: continue
        maxtb = max(maxtb, d['tau_b'])
        if d['Tm'] >= 5.89:
            n1 += 1; ratio = d['pk']/k; worst1 = min(worst1, ratio)
            if ratio < 0.374: f1 += 1; print('STAGE-1 FAIL', mu, k, d)
        if d['Tm'] >= 22:
            n2 += 1; ratio = d['r_l']/d['Tm']; worst2 = min(worst2, ratio)
            if ratio < 0.3: f2 += 1; print('STAGE-2 FAIL', mu, k, d)
    print('outcomes:', out)
    print(f'Stage 1 (T_- >= 5.89 => p_k >= 0.374 k): tested {n1}, failures {f1}, min p_k/k = {worst1:.4f}')
    print(f'Stage 2 (T_- >= 22 => r_l >= 0.3 T_-):   tested {n2}, failures {f2}, min r_l/T_- = {worst2:.4f}')
    print(f'max strip time tau_b over all segments = {maxtb:.4f}  (lemma: < 0.36)')
    print('ALL STAGE CHECKS PASSED' if f1 == 0 and f2 == 0 and maxtb < 0.36 else 'FAILURES')
