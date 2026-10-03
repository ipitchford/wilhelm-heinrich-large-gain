#!/usr/bin/env python3
"""Test the *derivation steps* of the a-priori lemmas on random NON-periodic orbit segments.
Every inequality tested here was derived using only the sign structure of the segment
(not periodicity), so it must hold on arbitrary orbit segments with that sign structure.
This exercises the 'medium-size' regimes that no actual cycle visits."""
import sys, numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

ELL = 2.0; PHI = 1-np.exp(-ELL); TH = 2**-0.5
def Psi(mu, u):
    if abs(mu-1) < 1e-12: return 1-(1+u)*np.exp(-u)
    return (1-np.exp(-mu*u))/mu - (np.exp(-mu*u)-np.exp(-u))/(1-mu)
ci = lambda a: (np.exp(a)-np.exp(a/2))*(1-(1+a/2)*np.exp(-a/2))
cii = lambda a: a*a/2-1+(1+a)*np.exp(-a)

def rhs(mu, k):
    def f(t, s):
        U, p, w, m = s
        ex = np.expm1(min(U, 50.0))
        return [p, -mu*p-w, -w+k*ex, k*ex]
    return f

def upx(t, s): return s[0]
upx.direction = 1; upx.terminal = True
def dnx(t, s): return s[0]
dnx.direction = -1; dnx.terminal = True

TOL = 1e-7
stats = {}
def chk(name, lhs, op, rhs):
    ok = (lhs <= rhs + TOL*(1+abs(rhs))) if op == '<=' else (lhs >= rhs - TOL*(1+abs(rhs)))
    a = stats.setdefault(name, [0, 0, np.inf])
    a[0] += 1; a[1] += (not ok)
    slack = (rhs-lhs) if op == '<=' else (lhs-rhs)
    a[2] = min(a[2], slack/(1+abs(rhs)))
    if not ok:
        print('  FAIL', name, lhs, op, rhs)
    return ok

def one(mu, k, Pd, wd, rng):
    f = rhs(mu, k)
    # ---- negative lobe from the down-crossing
    s0 = [-1e-12, -Pd, wd, 0.0]
    sol = solve_ivp(f, [0, 400.0], s0, method='DOP853', rtol=1e-10, atol=1e-12, events=upx, dense_output=True)
    if sol.status != 1: return 'no upcross'
    Tm = sol.t_events[0][0]; S = sol.sol
    U = lambda t: S(t)[0]; p = lambda t: S(t)[1]; w = lambda t: S(t)[2]
    tt = np.unique(np.concatenate([sol.t, np.linspace(0, Tm, 4001)])); Y = S(tt)
    Ug, pg, wg = Y[0], Y[1], Y[2]
    if np.any(wg <= -k): return 'y<0'
    iw = np.where(wg < 0)[0]
    if len(iw) == 0 or iw[0] == 0: return 'no s_w'
    s_w = brentq(w, tt[iw[0]-1], tt[iw[0]])
    ip = np.where(pg > 0)[0]
    if len(ip) == 0: return 'no s_0'
    s_0 = brentq(p, tt[ip[0]-1], tt[ip[0]])
    if np.any(pg[ip[0]:] < 0): return 'p not unimodal'
    D = -U(s_0); pk = p(Tm); wk = w(Tm)
    # (i)
    for a in (s_w, 0.6*s_w):
        v = -U(a/2)
        chk('(i) |U(a/2)| >= k phi c_i(a)', v, '>=', k*(1-np.exp(-v))*ci(a))
    # (ii)
    u = s_0-s_w; Pw = -p(s_w); vw = -U(s_w); phw = 1-np.exp(-vw)
    chk('(ii) P_w >= k phi_w (u-1+e^-u)', Pw, '>=', k*phw*(u-1+np.exp(-u)))
    chk('(ii) D - v_w >= k phi_w c_ii(u)', D-vw, '>=', k*phw*cii(u))
    # rising
    chk('G1 s_r >= pk/k', Tm-s_0, '>=', pk/k)
    chk('G1 D >= pk^2/2k', D, '>=', pk*pk/(2*k))
    for uu in np.linspace(0, min(pk/k, Tm-s_0), 7)[1:]:
        chk('G1 p(T- - u) >= pk - k u', p(Tm-uu), '>=', pk-k*uu)
        chk('G1 U(T- - u) <= -(pk u - k u^2/2)', U(Tm-uu), '<=', -(pk*uu-k*uu*uu/2))
    if D > ELL:
        ir = np.where((tt > s_0) & (Ug > -ELL))[0][0]
        s_l = brentq(lambda t: U(t)+ELL, tt[ir-1], tt[ir])
        r_l = s_l-s_0
        chk('(iii) r_l >= sqrt(2(D-l)/k)', r_l, '>=', np.sqrt(2*(D-ELL)/k))
        chk('(iii) p(s_l) >= k phi_l Psi(r_l)', p(s_l), '>=', k*PHI*Psi(mu, r_l))
        chk('(iii) pk >= p(s_l) - mu l', pk, '>=', p(s_l)-mu*ELL)
        for fr in (0.25, 0.5, 0.75):
            s = s_l+fr*(Tm-s_l); v = -U(s)
            chk('(iv) p(s) >= k phi(s) Psi(s-s_l)', p(s), '>=', k*(1-np.exp(-v))*Psi(mu, s-s_l))
            chk('(iv) |U(s)| >= p(s)(1-e^{-mu(T- - s)})/mu', v, '>=', p(s)*(1-np.exp(-mu*(Tm-s)))/mu)
    # ---- positive lobe
    if not (wk < 0 and pk > 0): return 'bad upcross'
    s1 = [1e-12, pk, wk, 0.0]
    sol2 = solve_ivp(f, [0, 50.0], s1, method='DOP853', rtol=1e-11, atol=1e-13, events=dnx, dense_output=True)
    if sol2.status != 1: return 'no downcross'
    Tp = sol2.t_events[0][0]; S2 = sol2.sol
    U2 = lambda t: S2(t)[0]; p2 = lambda t: S2(t)[1]; w2 = lambda t: S2(t)[2]; m2 = lambda t: S2(t)[3]
    t2 = np.unique(np.concatenate([sol2.t, np.linspace(0, Tp, 4001)])); Y2 = S2(t2)
    ic = np.where(Y2[1] < 0)[0][0]
    tau1 = brentq(p2, t2[ic-1], t2[ic])
    h = U2(tau1); Mp = m2(Tp); P1 = Y2[1][:ic+1].max()
    chk('C: M+ = w_d - w_k + p_k - p_d (segment)', abs(Mp - (w2(Tp)-wk+pk-p2(Tp)))/(1+Mp), '<=', 1e-6)
    for s in np.linspace(0, Tp, 9)[1:]:
        chk('F: w(s) <= m(s)', w2(s), '<=', m2(s))
    if pk >= 10:
        chk('F: P1 <= pk + 1.1 k sqrt(2/pk)', P1, '<=', pk+1.1*k*np.sqrt(2/pk))
    Q = TH**2*pk**2*((1-TH**2)*pk**2-2*mu*P1*h)/2 - TH*pk*k
    chk('F: M+(M+ + mu P1) >= Q (H:=h)', Mp*(Mp+mu*P1), '>=', Q)
    j = np.where(Y2[1][:ic+1] <= TH*pk)[0][0]
    s_th = brentq(lambda t: p2(t)-TH*pk, t2[j-1], t2[j]); Uth = U2(s_th)
    chk('F: e^{U_th} >= th pk((1-th^2)pk^2 - 2 mu P1 U_th)/(2k)', np.exp(Uth), '>=', TH*pk*((1-TH**2)*pk**2-2*mu*P1*Uth)/(2*k))
    chk('F: tau1 - s_th >= th pk/(mu P1 + M+)', tau1-s_th, '>=', TH*pk/(mu*P1+Mp))
    chk('G2: s_th <= U_th/(th pk)', s_th, '<=', Uth/(TH*pk))
    Xth = k*np.expm1(Uth)
    if Xth > 0:
        vc = (k+np.sqrt(k*k+(4/3)*Xth*TH*pk))/((2/3)*Xth)
        if vc <= 1:
            chk('G2: tau1 - s_th <= v_c', tau1-s_th, '<=', vc)
    Wc = w2(tau1)
    chk('G3: W_c >= pk e^{-tau1}/tau1', Wc, '>=', pk*np.exp(-tau1)/tau1)
    # fall after the crest: continue through the next negative lobe until p=0
    sc = [h, 0.0, Wc, 0.0]
    def pz(t, s): return s[1]
    pz.direction = 1; pz.terminal = True
    sol3 = solve_ivp(f, [0, 400.0], sc, method='DOP853', rtol=1e-10, atol=1e-12, events=pz, dense_output=True, first_step=1e-6)
    tb = sol3.t_events[0][0] if sol3.status == 1 else 400.0
    vplus = np.log1p(Wc/k)
    for v in np.linspace(0, min(vplus, tb), 7)[1:]:
        drop = h - sol3.sol(v)[0]
        chk('G3: drop(v) >= (Wc+k)(1-(1+v)e^-v) - k(v-1+e^-v)', drop, '>=', (Wc+k)*(1-(1+v)*np.exp(-v))-k*(v-1+np.exp(-v)))
        chk('G3: w(crest+v)+k >= e^-v (Wc+k)', sol3.sol(v)[2]+k, '>=', np.exp(-v)*(Wc+k))
    return 'ok'

if __name__ == '__main__':
    rng = np.random.default_rng(20261003)
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 300
    outcomes = {}
    for i in range(n):
        mu = float(rng.choice([1.0, 0.7, 0.3, 0.1, 0.03]))
        k = float(np.exp(rng.uniform(np.log(max(3*mu*(1+mu), 0.5)), np.log(3000))))
        scale = np.exp(rng.uniform(np.log(0.05), np.log(30)))
        Pd = scale*k**(1/3)*np.exp(rng.normal(0, 1.0))
        wd = scale*k**(2/3)*np.exp(rng.normal(0, 1.5))
        try:
            r = one(mu, k, Pd, wd, rng)
        except Exception as e:
            r = 'exc:'+type(e).__name__
        outcomes[r] = outcomes.get(r, 0)+1
    print('segment outcomes:', outcomes)
    bad = 0
    for name, (cnt, nf, slack) in stats.items():
        print(f'{name:62s} tests={cnt:5d} fails={nf:3d} min rel slack={slack:.3e}')
        bad += nf
    print('ALL SEGMENT CHECKS PASSED' if bad == 0 else f'{bad} FAILURES')
