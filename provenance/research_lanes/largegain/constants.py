#!/usr/bin/env python3
"""Check every numerical constant and scalar inequality used in theorem.tex.

All quantities are evaluated with mpmath at 40 significant digits.  Every check prints the two
sides and the margin; all margins are many orders of magnitude above the evaluation error, so no
interval arithmetic is needed.  Exit status 1 if any check fails.

Usage:  python constants.py            (K0 = 1000, the value used in the theorem)
        python constants.py 800        (try another K0 with the same fixed constants)
"""
import sys
from mpmath import mp, mpf, exp, log, sqrt, pi, findroot, e as E

mp.dps = 40
K0 = mpf(sys.argv[1]) if len(sys.argv) > 1 else mpf(1000)

FAILS = []
def chk(label, lhs, op, rhs):
    ok = (lhs >= rhs) if op == '>=' else (lhs <= rhs)
    print(f"  [{'ok' if ok else 'FAIL'}] {label}:  {mp.nstr(lhs, 8)} {op} {mp.nstr(rhs, 8)}")
    if not ok:
        FAILS.append(label)

# ---------------------------------------------------------------- fixed constants of the proof
ELL = mpf(2)                       # strip level
PHI = 1 - exp(-ELL)                # phi_ell
TH = 1/sqrt(2)                     # theta
A1W, A2W, A3W, A4W = mpf('2.06'), mpf('1.99'), mpf('1.48'), mpf('0.36')   # Stage-1 windows
C1 = mpf('0.374')                  # p_k >= C1 k           (Lemma 8)
CSHARP = mpf('0.055')              # z^sharp = CSHARP k/mu^2 (Lemma 11)
CT = CSHARP/mpf('0.3')             # T >= CT lam/mu
CP = mpf('0.819')                  # p_k >= CP lam
BF = mpf('1.02')                   # b in v_f
NU_FRAC = mpf('0.5')               # nu = mu/2

Psi1 = lambda u: 1 - (1+u)*exp(-u)
ci = lambda a: (exp(a)-exp(a/2))*(1-(1+a/2)*exp(-a/2))
cii = lambda a: a*a/2 - 1 + (1+a)*exp(-a)
cb = lambda a: Psi1(a/2)*(1-exp(-a/2))

print(f"K0 = {K0}")
print("== elementary constants")
v2 = findroot(lambda v: v - 2*(1-exp(-v)), 1.6)
PHI2 = 1 - exp(-v2)
print(f"  v2 = {mp.nstr(v2, 12)},  phi2 = 1-exp(-v2) = v2/2 = {mp.nstr(PHI2, 12)}")
chk("phi2 >= 0.7968", PHI2, '>=', mpf('0.7968'))
PHI2 = mpf('0.7968')              # the rounded-down value used in the text
chk("1-exp(-3) >= 0.95", 1-exp(-3), '>=', mpf('0.95'))
chk("c_i(3) e^{-3} >= 0.34", ci(3)*exp(-3), '>=', mpf('0.34'))
chk("max_{0<mu<=1} 3 mu log(1/mu) = 3/e <= 1.1037", 3/E, '<=', mpf('1.1037'))
chk("ell/phi_ell (strip ratio) <= 2.3131", ELL/PHI, '<=', mpf('2.3131'))
chk("s_+ lemma: 1.21 (1 - 1.1 sqrt(2/10)/3) >= 1", mpf('1.21')*(1-mpf('1.1')*sqrt(mpf(2)/10)/3), '>=', 1)

print("== Lemma 2/3: frequency and first foothold")
om0 = sqrt(3*K0**(mpf(2)/3) - 4*K0**(mpf(1)/3) - 1)/2
T0 = 2*pi - pi/om0
print(f"  omega_0(K0) = {mp.nstr(om0, 10)},  T_0(K0) = 2 pi - pi/omega_0 = {mp.nstr(T0, 10)}")
chk("windows fit: a1+a2+a3+a4 <= T_0(K0)", A1W+A2W+A3W+A4W, '<=', T0)

print("== Lemma 7/8: Stage 1")
chk("K0 c_i(a1) >= 2", K0*ci(A1W), '>=', 2)
chk("kappa_ii(a2) = K0 (a2-1+e^{-a2})/(K0^{1/3}+1) >= 2", K0*(A2W-1+exp(-A2W))/(K0**(mpf(1)/3)+1), '>=', 2)
chk("K0 c_b(a4) >= ell/phi_ell", K0*cb(A4W), '>=', ELL/PHI)
Ri = sqrt(2*(PHI2*ci(A1W) - ELL/K0)); Rii = sqrt(2*(PHI2*cii(A2W) - ELL/K0))
print(f"  R_i = {mp.nstr(Ri, 8)}, R_ii = {mp.nstr(Rii, 8)}, a3 = {A3W}")
chk("R_i >= a3", Ri, '>=', A3W); chk("R_ii >= a3", Rii, '>=', A3W)
chk("c1: phi Psi_1(a3) - ell/K0 >= C1", PHI*Psi1(A3W) - ELL/K0, '>=', C1)
chk("C1 K0 >= 10 (needed for the s_+ bound)", C1*K0, '>=', 10)

print("== Lemma 9/10: spike mass with p_k >= C1 k")
Ht = 2*log(K0+1) + 3/E
eps1 = mpf('1.1')*sqrt(2)/(C1**mpf('1.5')*sqrt(K0))
d1 = 4*Ht*(1+eps1)/(C1*K0) + 4*sqrt(2)/(C1**3*K0**2)
print(f"  Htilde(K0) = {mp.nstr(Ht, 8)}, eps1bar = {mp.nstr(eps1, 6)}, delta_1 = {mp.nstr(d1, 6)}")
chk("delta_1 < 1", d1, '<=', mpf('0.999'))
T1 = C1**2*K0*sqrt(1-d1)/sqrt(8) - mpf('0.5')
print(f"  T_1(K0) = {mp.nstr(T1, 8)}")
chk("T_1(K0) >= 22", T1, '>=', 22)

print("== Lemma 11: Stage 2 (continuum bootstrap)")
chk("Stage-2 r_l: 0.89/1.89 >= 0.4708", mpf('0.89')/mpf('1.89'), '>=', mpf('0.4708'))
chk("Stage-2 r_l: 1/1.89 <= 0.5292", 1/mpf('1.89'), '<=', mpf('0.5292'))
chk("Stage-2 r_l: 0.4708*a4 + 0.5292 <= 0.6987", mpf('0.4708')*A4W + mpf('0.5292'), '<=', mpf('0.6987'))
chk("Stage-2 r_l: 0.1708*22 - 0.9416 log 22 - 0.6987 >= 0", mpf('0.1708')*22 - mpf('0.9416')*log(22) - mpf('0.6987'), '>=', 0)
chk("Stage-2 r_l: derivative 0.1708 - 0.9416/22 >= 0", mpf('0.1708') - mpf('0.9416')/22, '>=', 0)
chk("Stage-2 r_l: T' = 22 - 2 log 22 - a4 >= 2.79 (case u<2)", 22 - 2*log(22) - A4W, '>=', mpf('2.79'))
chk("Stage-2 r_l: 0.5292*2.79 >= 2 - 0.5292", mpf('0.5292')*mpf('2.79'), '>=', 2 - mpf('0.5292'))
chk("Stage-2: 2 log 22 >= 3 (so c_i(a) >= 0.34 e^a)", 2*log(22), '>=', 3)
chk("Stage-2: sqrt(2*phi2*0.34 - 2 ell/(K0*22^2)) >= 0.73", sqrt(2*PHI2*mpf('0.34') - 2*ELL/(K0*22**2)), '>=', mpf('0.73'))
chk("Stage-2: (0.89u-1)^2 <= 2 phi2 (u^2/2-1) - 2 ell/K0 at u=2", (mpf('0.89')*2-1)**2, '<=', 2*PHI2*(2-1) - 2*ELL/K0)
chk("Stage-2: coefficient phi2 - 0.89^2 >= 0 (monotone in u)", PHI2 - mpf('0.89')**2, '>=', 0)
chk("Stage-2: kappa_ii(2) >= 2", K0*(2-1+exp(-2))/(K0**(mpf(1)/3)+1), '>=', 2)
z10 = mpf('0.3')*T1 - 3
print(f"  z_1(K0) = 0.3 T_1 - 3 = {mp.nstr(z10, 8)}")
chk("z_1(K0) >= 3.6 (i.e. T_1 >= 22)", z10, '>=', mpf('3.6'))
chk("0.01 phi K0 z10/(1+z10) >= ell (0.95 -> 0.94)", mpf('0.01')*PHI*K0*z10/(1+z10), '>=', ELL)
kap0 = (mpf('0.94')*PHI)**2*sqrt(1-d1)/sqrt(8)
print(f"  kappa_0 = {mp.nstr(kap0, 8)},  0.3 kappa_0 = {mp.nstr(0.3*kap0, 8)}")
slope = mpf('0.3')*kap0 - mpf('0.3')*C1**2*(1+1/z10)**2/sqrt(8)
chk("endpoint z_1: 0.3 kappa_0 - 0.3 C1^2 (1+1/z10)^2/sqrt8 >= 0", slope, '>=', 0)
hs = CSHARP*(1+mpf('3.15')/(CSHARP*K0))*(1+1/(CSHARP*K0))**2
chk("endpoint z^sharp: c#(1+3.15/(c# K0))(1+1/(c# K0))^2 <= 0.3 kappa_0", hs, '<=', mpf('0.3')*kap0)
chk("c_p: 0.95 phi (1-exp(-c# K0)) - ell/K0 >= CP", mpf('0.95')*PHI*(1-exp(-CSHARP*K0)) - ELL/K0, '>=', CP)
chk("interval placement: 0.3 C1^2/sqrt8 <= 0.015 <= c#  (z_1(k) < z^sharp)", mpf('0.3')*C1**2/sqrt(8), '<=', mpf('0.015'))
chk("interval placement: 0.015 <= c#", mpf('0.015'), '<=', CSHARP)
print(f"  c_T = c#/0.3 = {mp.nstr(CT, 8)}")

print("== Lemma 12: crest duration at lambda = K0 (all majorants are decreasing in lambda)")
lam = K0
Hhat = 3*log(lam) - log(K0) + 2*log(1+1/K0)
Lhat = log(8*lam/3)
Htl = 2*log(lam+1) + 3/E
e1 = mpf('1.1')*sqrt(2)/(CP**mpf('1.5')*sqrt(lam))
dp = 4*(1+e1)*Htl/(CP*lam) + 4/(TH*CP**3*lam**2)
chk("delta' < 1", dp, '<=', mpf('0.999'))
ac = sqrt(12)/sqrt(1-dp) + 12/(TH*CP**2*lam*(1-dp))
A1 = sqrt(2)*Hhat + ac
tau1 = A1/(CP*lam)
W = CP**2*lam**2*exp(-tau1)/A1
vf = BF*sqrt(2*(Hhat+Lhat)/W)
print(f"  Hhat={mp.nstr(Hhat,8)} Lhat={mp.nstr(Lhat,8)} eps1={mp.nstr(e1,6)} delta'={mp.nstr(dp,6)} a_c={mp.nstr(ac,8)}")
print(f"  A_1={mp.nstr(A1,8)} tau1bar={mp.nstr(tau1,6)} What={mp.nstr(W,8)} v_f={mp.nstr(vf,6)}")
chk("(a) c_p^2 lam >= 2 Lhat  (p_k^2 >= 2 k L_0)", CP**2*lam, '>=', 2*Lhat)
chk("(a') C_P lam >= 10", CP*lam, '>=', 10)
chk("(c) v_c <= a_c/(c_p lam) <= 1", ac/(CP*lam), '<=', 1)
chk("(d) v_f <= log(1+What/lam)", vf, '<=', log(1+W/lam))
chk("(e) b^2 [1-(2/3)(1+lam/What) v_f] >= 1", BF**2*(1-mpf(2)/3*(1+lam/W)*vf), '>=', 1)
chk("(e') v_f <= 1", vf, '<=', 1)
chk("monotonicity: Hhat(K0) >= 3", Hhat, '>=', 3)
chk("monotonicity: Hhat + Lhat >= 4", Hhat + Lhat, '>=', 4)
chk("monotonicity: 2 < 2 log(2) + 3/e (Htilde_lam/lam decreasing for lam>=1)", mpf(2), '<=', 2*log(2) + 3/E)
tauc = (2*Lhat + A1)/(CP*lam) + vf
print(f"  tau_c_hat(K0) = {mp.nstr(tauc, 8)}")
# majorant C_tau(lam) = c_a + c_b log(lam) with tau_c <= C_tau/(c_p lam):
kf = BF*exp(tau1/2)/sqrt(2)
Ctau = 2*Lhat + A1 + kf*(Hhat + Lhat + A1)
chk("AM-GM majorant: c_p lam v_f <= kf (Hhat+Lhat+A1)", CP*lam*vf, '<=', kf*(Hhat+Lhat+A1))
cb_log = 2 + 3*sqrt(2) + kf*(3 + 1 + 3*sqrt(2))       # coefficient of log(lam) in C_tau
print(f"  C_tau(K0) = {mp.nstr(Ctau, 8)}, log-coefficient c_b = {mp.nstr(cb_log, 8)}")
chk("C_tau(K0) <= 71.79 and c_b <= 12.28 (majorant Cbar_tau)", max(Ctau - mpf('71.79'), cb_log - mpf('12.28')), '<=', 0)
chk("C_tau/sqrt(lam) decreasing for lam>=K0:  c_b <= C_tau(K0)/2", cb_log, '<=', Ctau/2)
chk("C_tau/lam decreasing (same criterion, weaker): c_b <= C_tau(K0)", cb_log, '<=', Ctau)

print("== Lemma 13 + final inequality (nu = mu/2, s = 2 sqrt(lam))")
chk("s = 2 sqrt(K0) >= sqrt(2)", 2*sqrt(K0), '>=', sqrt(2))
logterm = log(1 + 4*sqrt(lam) + mpf(16)/3*lam)
Cbar = mpf('71.79')
Phi_ = CT*lam/4 - 2*Cbar/(CP*sqrt(lam)) - logterm
print(f"  (c_T/4) lam = {mp.nstr(CT*lam/4, 8)},  2 C_tau/(c_p sqrt(lam)) = {mp.nstr(2*Ctau/(CP*sqrt(lam)), 8)},  log term = {mp.nstr(logterm, 8)}")
chk("FINAL: (c_T/4) K0 - 2 C_tau(K0)/(c_p sqrt(K0)) - log(1+4 sqrt(K0)+16 K0/3) > 0", Phi_, '>=', mpf('1e-6'))
chk("monotonicity of the final expression: c_T/4 >= 1/K0", CT/4, '>=', 1/K0)
# quantitative by-product: <eta> <= (1 - 1/4 + q) mu
q = (2*Cbar/(CP*sqrt(lam)) + logterm)/(CT*lam)
chk("q(K0) < 0.0772", q, '<=', mpf('0.0772'))
chk("majorant monotone: 12.28 <= Cbar/2", mpf('12.28'), '<=', Cbar/2)
print(f"  by-product: <eta> - mu <= -(1/4 - q) mu with q = {mp.nstr(q, 6)}  ->  <eta> <= {mp.nstr(mpf(3)/4+q, 6)} mu")
chk("by-product: 3/4 + q <= 0.86", mpf(3)/4 + q, '<=', mpf('0.86'))

print()
if FAILS:
    print("FAILED CHECKS:", FAILS); sys.exit(1)
print(f"ALL CONSTANT CHECKS PASSED for K0 = {K0}")
