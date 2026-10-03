#!/usr/bin/env python3
"""Numerical sanity test of every a-priori lemma of theorem.tex on computed cycles (floating point).
Each line prints  name: lhs  (<= or >=)  rhs  -> OK/FAIL.   A FAIL means the lemma as stated is wrong."""
import sys, numpy as np
from scipy.optimize import brentq
from lgcycle import find_cycle_lg, cycle_dense

ELL = 2.0
PHI = 1 - np.exp(-ELL)
TH = 1/np.sqrt(2.0)

def Psi(mu, u):
    if abs(mu-1) < 1e-12:
        return 1 - (1+u)*np.exp(-u)
    return (1-np.exp(-mu*u))/mu - (np.exp(-mu*u)-np.exp(-u))/(1-mu)

def eq_data(mu, k):
    A = 1+mu
    g = brentq(lambda g: g*(g-1)*(g-mu)-k, A, k**(1/3)+A+1)
    sig = (g-A)/2
    zE = (g-1)*(g-mu)
    om = np.sqrt(zE - sig*sig)
    return dict(g=g, sig=sig, om=om, etaE=g-1, zetaE=zE)

def root(f, a, b):
    return brentq(f, a, b, xtol=1e-13, rtol=1e-13)

def analyse(mu, k, verbose=True):
    lam = k/mu; A = 1+mu
    cyc = find_cycle_lg(mu, k)
    T = cyc['T']
    sol = cycle_dense(mu, k, cyc)
    S = sol.sol
    U = lambda t: S(t)[0]; p = lambda t: S(t)[1]; w = lambda t: S(t)[2]
    # coarse scan for events
    tt = np.unique(np.concatenate([sol.t, np.linspace(0, T, 200001)]))
    Y = S(tt)
    Ug, pg, wg = Y[0], Y[1], Y[2]
    def first_cross(arr, sgn_from, start=0):
        idx = np.where((np.sign(arr[start:-1]) == sgn_from) & (np.sign(arr[start+1:]) == -sgn_from))[0]
        return idx + start
    i_w = first_cross(wg, 1)[0];  s_w = root(w, tt[i_w], tt[i_w+1])
    i_0 = first_cross(pg, -1)[0]; s_0 = root(p, tt[i_0], tt[i_0+1])
    i_k = first_cross(Ug, -1, 1)[0]; Tm = root(U, tt[i_k], tt[i_k+1])
    i_c = first_cross(pg, 1, i_k)[0]; t_c = root(p, tt[i_c], tt[i_c+1])
    Tp = T - Tm; tau1 = t_c - Tm
    Pd = -p(0.0); wd = w(0.0); pk = p(Tm); wk = w(Tm); D = -U(s_0); h = U(t_c)
    ed = eq_data(mu, k)
    res = []
    def chk(name, lhs, op, rhs):
        ok = (lhs <= rhs) if op == '<=' else (lhs >= rhs)
        res.append(ok)
        if verbose:
            print(f'  {name:42s} {lhs:14.6g} {op} {rhs:14.6g}  {"OK" if ok else "FAIL"}')
    if verbose:
        print(f'mu={mu} k={k} lam={lam:.4g} T={T:.6g} T-={Tm:.6g} T+={Tp:.4g} s_w={s_w:.4g} s_0={s_0:.4g} '
              f'D={D:.5g} h={h:.4g} pk={pk:.5g} (pk/lam={pk/lam:.4f}) Pd={Pd:.5g} wd={wd:.5g} wk={wk:.5g}')
    # Lemma A
    chk('A: T+ < pi/omega', Tp, '<=', np.pi/ed['om'])
    chk('B: T > 2pi/sqrt(mu)', T, '>=', 2*np.pi/np.sqrt(mu))
    # sign pattern
    chk('sign: w_d>0', wd, '>=', 0); chk('sign: w_k<0', -wk, '>=', 0); chk('order s_w<s_0', s_w, '<=', s_0)
    chk('Lambda_E(0)>0: wd - etaE*Pd', wd - ed['etaE']*Pd, '>=', 0)
    chk('Lambda_E(T-)>0: etaE*pk - |wk|', ed['etaE']*pk + wk, '>=', 0)
    chk('x_min<k_H: D > log rho', D, '>=', np.log(k/(mu*A)))
    # mass
    tl = np.linspace(Tm, T, 400001)   # fine grid on the positive lobe only (the spike is ~1/lam wide)
    Mp = float(np.trapezoid(k*np.expm1(S(tl)[0]), tl))
    chk('C: mass identity (rel err, quadrature)', abs(Mp - (wd+Pd-wk+pk))/Mp, '<=', 1e-6)
    chk('C: M+ < k T-', Mp, '<=', k*Tm)
    # rising phase
    i_l = np.where((Ug[i_0:i_k+1] < -ELL))[0]
    s_l = root(lambda t: U(t)+ELL, tt[i_0+i_l[-1]], tt[i_0+i_l[-1]+1])
    r_l = s_l - s_0; tau_b = Tm - s_l
    chk('E: p(s_l) >= k phi Psi(r_l)', p(s_l), '>=', k*PHI*Psi(mu, r_l))
    chk('E: p_k >= p(s_l) - mu*ell', pk, '>=', p(s_l) - mu*ELL)
    chk('E: r_l >= sqrt(2(D-ell)/k)', r_l, '>=', np.sqrt(2*(D-ELL)/k))
    if Tm >= 22:
        chk('Stage2: r_l >= 0.3 T-', r_l, '>=', 0.3*Tm)
    # descending-phase lemmas (i),(ii) evaluated at the actual durations
    ci = lambda a: (np.exp(a)-np.exp(a/2))*(1-(1+a/2)*np.exp(-a/2))
    cii = lambda a: a*a/2 - 1 + (1+a)*np.exp(-a)
    if k*ci(s_w) >= 2:
        chk('D(i): D >= 0.7968 k c_i(s_w)', D, '>=', 0.7968*k*ci(s_w))
    u = s_0 - s_w
    kap = k*(u-1+np.exp(-u))/(ed['g']-mu)
    Pw = -p(s_w); vw = -U(s_w)
    chk('D(ii): P_w >= k phi(v_w)(u-1+e^-u)', Pw, '>=', k*(1-np.exp(-vw))*(u-1+np.exp(-u)))
    chk('D(ii): P_w < (gE-mu) v_w', Pw, '<=', (ed['g']-mu)*vw)
    if kap >= 2:
        chk('D(ii): D >= 0.7968 k c_ii(u)', D, '>=', 0.7968*k*cii(u))
    # Lemma F (spike mass)
    mrise = (tt >= Tm) & (tt <= t_c)
    P1 = pg[mrise].max()
    if pk >= 10:
        chk('F: P1 <= pk(1+eps1)', P1, '<=', pk*(1+1.1*np.sqrt(2)*k/pk**1.5))
    M0 = (lam+mu)*(lam+1)/mu; H = np.log(M0)
    chk('box: h <= H', h, '<=', H)
    Q = TH**2*pk**2*((1-TH**2)*pk**2 - 2*mu*P1*H)/2 - TH*pk*k
    chk('F: M+(M+ + mu P1) >= Q', Mp*(Mp+mu*P1), '>=', Q)
    # s_theta
    ir = np.where(mrise)[0]
    j = ir[np.where(pg[ir] <= TH*pk)[0][0]]
    s_th = root(lambda t: p(t)-TH*pk, tt[j-1], tt[j]); U_th = U(s_th)
    Eth = TH*pk*((1-TH**2)*pk**2 - 2*mu*P1*H)/(2*k)
    chk('F: exp(U_theta) >= E_theta', np.exp(U_th), '>=', Eth)
    Xth = k*(Eth-1)
    chk('F: tau1 - s_theta >= theta pk/(mu P1+M+)', t_c - s_th, '>=', TH*pk/(mu*P1+Mp))
    # Lemma G
    nu = mu/2; c0 = (mu-nu)*(1-nu)*(A-nu); L0 = np.log(k/c0)
    i_r = np.where(Ug[i_0:i_k+1] < -L0)[0]
    t_r = root(lambda t: U(t)+L0, tt[i_0+i_r[-1]], tt[i_0+i_r[-1]+1])
    chk('G1: tau_rise <= 2 L0/pk', Tm - t_r, '<=', 2*L0/pk)
    chk('G1: D >= pk^2/(2k)', D, '>=', pk**2/(2*k))
    vc = (k + np.sqrt(k*k + (4/3)*Xth*TH*pk))/((2/3)*Xth) if Xth > 0 else np.inf
    chk('G2: s_theta <= H/(theta pk)', s_th - Tm, '<=', H/(TH*pk))
    chk('G2: tau1 - s_theta <= v_c', t_c - s_th, '<=', vc)
    Wc = w(t_c)
    chk('G3: W_c >= pk e^-tau1/tau1', Wc, '>=', pk*np.exp(-tau1)/tau1)
    # fall time: crest -> U=-L0 (wraps around the period end)
    i_f = np.where(Ug[:i_0+1] > -L0)[0]
    t_f = root(lambda t: U(t)+L0, tt[i_f[-1]], tt[i_f[-1]+1])
    fall = (T - t_c) + t_f
    Wlow = pk*np.exp(-tau1)/tau1
    vf0 = np.sqrt(2*(H+L0)/Wlow)
    # solve (W/2) v^2 [1-(2/3)(1+k/W)v] = H+L0
    g = lambda v: (Wlow/2)*v*v*(1-(2/3)*(1+k/Wlow)*v) - (H+L0)
    try:
        vf = root(g, vf0*0.5, min(1.0/((1+k/Wlow)), 2*vf0+0.5))
    except Exception:
        vf = np.inf
    chk('G3: fall time <= v_f', fall, '<=', vf)
    tau_c = (Tm - t_r) + Tp + t_f
    # final Lyapunov bound
    alpha = 1/(1-nu); beta = alpha/(A-nu); s = 2*np.sqrt(lam)
    mcrest = (Ug > -L0)
    xg = k*np.exp(np.minimum(Ug, 60))
    Mc = float(np.trapezoid(np.where(mcrest, xg, 0.0), tt))
    logmt = S(T)[3]
    bound = -nu*(T-tau_c) + (s-mu)*tau_c + Mc/s**2 + np.log(1+alpha*s+beta*s*s)
    chk('H: log m_t <= Lyapunov bound (actual tau_c, M_c)', logmt, '<=', bound)
    bound2 = -nu*(T-tau_c) + (s-mu)*tau_c + k*T/s**2 + np.log(1+alpha*s+beta*s*s)
    chk('H: log m_t <= bound with M_c<=kT', logmt, '<=', bound2)
    if k >= 1000 and mu <= 1:
        chk('THM (i): T >= 0.1833 lam/mu', T, '>=', 0.1833*lam/mu)
        chk('THM (ii): p_k >= 0.819 lam', pk, '>=', 0.819*lam)
        chk('THM (iii): tau_c <= (71.79+12.28 log(lam/1000))/(0.819 lam)', tau_c, '<=', (71.79+12.28*np.log(lam/1000))/(0.819*lam))
        chk('THM (iv): <eta> <= 0.83 mu', mu + logmt/T, '<=', 0.83*mu)
        chk('Lemma boot: T- >= 43.89', Tm, '>=', 43.89)
        chk('Lemma stage1: p_k >= 0.374 k', pk, '>=', 0.374*k)
    if verbose:
        print(f'  [info] tau_c={tau_c:.4g} tau_rise={Tm-t_r:.4g} tau1={tau1:.4g} fall={fall:.4g} v_f={vf:.4g} '
              f'log m_t={logmt:.5g} bound={bound:.5g} bound2={bound2:.5g}  M+/pk^2={Mp/pk**2:.4f} T mu/lam={T*mu/lam:.4f} r_l/T-={r_l/Tm:.3f}')
    return all(res)

if __name__ == '__main__':
    cases = [(1.0, 16.0), (1.0, 60.0), (1.0, 300.0), (1.0, 1000.0), (1.0, 3000.0), (0.5, 400.0), (0.5, 1000.0), (0.3, 15.6), (0.3, 1000.0), (0.1, 100.0), (0.02, 20.0)]
    if len(sys.argv) > 2: cases = [(float(sys.argv[1]), float(sys.argv[2]))]
    allok = True
    for mu, k in cases:
        ok = analyse(mu, k)
        allok &= ok
        print('  => all OK' if ok else '  => SOME FAIL')
    print('ALL LEMMA CHECKS PASSED' if allok else 'THERE WERE FAILURES')
