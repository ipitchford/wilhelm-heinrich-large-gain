#!/usr/bin/env python3
# --- v3.0 bundle copy of research file `lanes/largegain/verify/v_consts.py`.
# --- Changes from the original: chk() also appends each check to RECORDS.
"""Referee script A: independent recomputation of every numerical constant and scalar
inequality used in theorem.tex (mpmath, 50 digits), plus grid checks of the monotonicity
claims that reduce side conditions to lambda = K0.  Written from the .tex only; does not
import or read the lane's constants.py."""
import sys
import mpmath as mp
mp.mp.dps = 50
F = mp.mpf
fails = []
RECORDS = []


def chk(name, cond, detail=""):
    tag = "ok  " if cond else "FAIL"
    print(f"[{tag}] {name} {detail}")
    if not cond:
        fails.append(name)
    RECORDS.append(dict(name=name, passed=bool(cond), detail=str(detail)))


K0 = F(1000)
ell = F(2)
phl = 1 - mp.e ** (-ell)
theta = 1 / mp.sqrt(2)
v2 = mp.findroot(lambda v: v - 2 * (1 - mp.e ** (-v)), 1.6)
ph2 = F("0.7968")
chk("v2 root, phi2<=1-exp(-v2)", ph2 <= 1 - mp.e ** (-v2), f"v2={mp.nstr(v2,8)} 1-e^-v2={mp.nstr(1-mp.e**(-v2),8)}")
# monotonicity of v/(1-e^-v)
g = lambda v: v / (1 - mp.e ** (-v))
chk("v/(1-e^-v) increasing (grid)", all(g(F(i) / 100) < g(F(i + 1) / 100) for i in range(1, 3000)))

Psi1 = lambda u: 1 - (1 + u) * mp.e ** (-u)


def Psi(mu, u):
    mu = F(mu)
    if mu == 1:
        return Psi1(u)
    # int_0^u e^{-mu(u-r)}(1-e^{-r}) dr
    return (1 - mp.e ** (-mu * u)) / mu - (mp.e ** (-u) - mp.e ** (-mu * u)) / (mu - 1)


ci = lambda a: (mp.e ** a - mp.e ** (a / 2)) * (1 - (1 + a / 2) * mp.e ** (-a / 2))
cii = lambda a: a * a / 2 - 1 + (1 + a) * mp.e ** (-a)
cb = lambda a: Psi1(a / 2) * (1 - mp.e ** (-a / 2))

# cross-check closed forms against quadrature
for a in [F("0.36"), F("1.48"), F("1.99"), F("2.06"), F(3)]:
    q = mp.quad(lambda s: s * mp.e ** (-s), [0, a / 2]) * (mp.e ** a - mp.e ** (a / 2))
    chk(f"c_i closed form a={a}", abs(q - ci(a)) < F(10) ** -30)
    q = mp.quad(lambda r: r * (1 - mp.e ** (-r)), [0, a])
    chk(f"c_ii closed form a={a}", abs(q - cii(a)) < F(10) ** -30)
    for mu in [F("0.02"), F("0.3"), F(1)]:
        q = mp.quad(lambda r: mp.e ** (-mu * (a - r)) * (1 - mp.e ** (-r)), [0, a])
        chk(f"Psi_mu closed form mu={mu} a={a}", abs(q - Psi(mu, a)) < F(10) ** -30)
        chk(f"Psi_mu>=Psi_1 mu={mu} a={a}", Psi(mu, a) >= Psi1(a))

# ---- Lemma 6
om0 = lambda k: mp.sqrt(3 * k ** (F(2) / 3) - 4 * k ** (F(1) / 3) - 1) / 2
T0 = lambda k: 2 * mp.pi - mp.pi / om0(k)
print("omega0(K0)=", mp.nstr(om0(K0), 10), " T0(K0)=", mp.nstr(T0(K0), 10))
a1, a2, a3, a4 = F("2.06"), F("1.99"), F("1.48"), F("0.36")
chk("T0(K0) >= a1+a2+a3+a4 = 5.89", T0(K0) >= a1 + a2 + a3 + a4 and a1 + a2 + a3 + a4 == F("5.89"),
    f"margin={mp.nstr(T0(K0)-F('5.89'),6)}")

# omega >= omega0(k), gamma_E - mu < k^{1/3}+1 on a grid (true roots)
import numpy as np
bad = 0
minratio = 1e9
for mu in [1e-4, 1e-3, 0.02, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0]:
    for k in np.concatenate([np.geomspace(max(1.0, 1.0001 * mu * (1 + mu)), 1e7, 200)]):
        if k <= mu * (1 + mu):
            continue
        r = np.roots([1, 1 + mu, mu, k])
        gam = -min(r.real[np.abs(r.imag) < 1e-9])
        om = max(r.imag)
        rad = 3 * k ** (2 / 3) - 4 * k ** (1 / 3) - 1
        if rad > 0:
            minratio = min(minratio, om / (0.5 * np.sqrt(rad)))
            if om < 0.5 * np.sqrt(rad) * (1 - 1e-12):
                bad += 1
        if not (gam - mu < k ** (1 / 3) + 1) or not (gam ** 3 > k):
            bad += 1
chk("omega>=omega0(k), gamma_E-mu<k^(1/3)+1, gamma_E^3>k on grid", bad == 0, f"min omega/omega0={minratio:.5f}")

# ---- Lemma 11 (Stage 1)
chk("K0 c_i(a1) >= 2", K0 * ci(a1) >= 2, mp.nstr(K0 * ci(a1), 8))
kap2 = lambda k, a: k * (a - 1 + mp.e ** (-a)) / (k ** (F(1) / 3) + 1)
chk("kappa_ii(a2) >= 2 at K0", kap2(K0, a2) >= 2, mp.nstr(kap2(K0, a2), 8))
chk("K0 c_b(a4) >= ell/phi_ell", K0 * cb(a4) >= ell / phl, f"{mp.nstr(K0*cb(a4),8)} vs {mp.nstr(ell/phl,8)}")
r1 = mp.sqrt(2 * (ph2 * ci(a1) - ell / K0))
r2 = mp.sqrt(2 * (ph2 * cii(a2) - ell / K0))
chk("case1 r_l >= 1.48", r1 >= a3, mp.nstr(r1, 8))
chk("case2 r_l >= 1.48", r2 >= a3, mp.nstr(r2, 8))
c1 = F("0.374")
chk("phi_l Psi_1(1.48) - l/K0 >= c1", phl * Psi1(a3) - ell / K0 >= c1, mp.nstr(phl * Psi1(a3) - ell / K0, 8))
chk("rho >= K0/2 > e^l", K0 / 2 > mp.e ** ell)

# ---- Lemma 12
sbar = F("1.1") * mp.sqrt(F(2) / 10)
chk("sbar<=0.492 at p_k=10", sbar <= F("0.492"), mp.nstr(sbar, 8))
chk("1.21(1-sbar/3) >= 1.0115 > 1", F("1.21") * (1 - sbar / 3) >= F("1.0115"), mp.nstr(F("1.21") * (1 - sbar / 3), 8))

# ---- Lemma 14 (Stage 2)
f3 = lambda a: (1 - mp.e ** (-a / 2)) * (1 - (1 + a / 2) * mp.e ** (-a / 2))
chk("c_i(a)e^-a = f3(a)", abs(ci(F(3)) * mp.e ** -3 - f3(F(3))) < F(10) ** -40)
chk("f3(3) >= 0.34", f3(F(3)) >= F("0.34"), mp.nstr(f3(F(3)), 8))
chk("sqrt(0.68 phi2 T^2 - 2l/K0) >= 0.73T at T=22", mp.sqrt(F("0.68") * ph2 * 22 ** 2 - 2 * ell / K0) >= F("0.73") * 22)
chk("0.68 phi2 > 0.73^2", F("0.68") * ph2 - F("0.73") ** 2 > 0 and (F("0.68") * ph2 - F("0.73") ** 2) * 484 >= 2 * ell / K0)
Tp22 = 22 - 2 * mp.log(22) - F("0.36")
chk("T' >= 15.4 at T=22", Tp22 >= F("15.4"), mp.nstr(Tp22, 8))
chk("kappa_ii(2) >= 2 at K0", kap2(K0, F(2)) >= 2)
du = lambda u: (ph2 - F("0.89") ** 2) * u * u + F("1.78") * u - 1 - 2 * ph2 - 2 * ell / K0
chk("quadratic difference positive at u=2 and phi2>0.89^2", du(F(2)) > 0 and ph2 > F("0.89") ** 2, mp.nstr(du(F(2)), 6))
chk("(0.89T'-1)/1.89 >= 0.4708T'-0.5292", F("0.89") / F("1.89") >= F("0.4708") and F(1) / F("1.89") <= F("0.5292"))
chk("T'-2 >= 0.4708T'-0.5292 for T'>=2.79", (1 - F("0.4708")) * F("2.79") >= 2 - F("0.5292"))
chk("0.4708*0.36+0.5292 <= 0.6987 and 2*0.4708=0.9416", F("0.4708") * F("0.36") + F("0.5292") <= F("0.6987") and 2 * F("0.4708") == F("0.9416"))
gT = lambda T: F("0.1708") * T - F("0.9416") * mp.log(T) - F("0.6987")
chk("0.1708T-0.9416logT-0.6987 >0 at 22, increasing", gT(F(22)) > 0 and F("0.1708") - F("0.9416") / 22 > 0, mp.nstr(gT(F(22)), 6))
chk("1-e^-3 >= 0.95", 1 - mp.e ** -3 >= F("0.95"))

# ---- Lemma 15 (bootstrap)
Ht = lambda k: 2 * mp.log(k + 1) + 3 / mp.e
eb1 = lambda k: F("1.1") * mp.sqrt(2) * c1 ** (-F(3) / 2) / mp.sqrt(k)
d1 = lambda k: 4 * Ht(k) * (1 + eb1(k)) / (c1 * k) + 4 * mp.sqrt(2) / (c1 ** 3 * k ** 2)
print("Htilde(K0)=", mp.nstr(Ht(K0), 8), "eps1bar=", mp.nstr(eb1(K0), 8), "delta1(K0)=", mp.nstr(d1(K0), 8))
chk("delta1(K0) < 1 (=0.19402)", d1(K0) < 1 and abs(d1(K0) - F("0.19402")) < F("1e-5"))
chk("delta1 decreasing on grid", all(d1(K0 * F(10) ** (F(i) / 20)) > d1(K0 * F(10) ** (F(i + 1) / 20)) for i in range(0, 200)))
# mu*H <= Htilde check on grid
bad = 0
for mu in [F(10) ** (-F(i) / 10) for i in range(0, 80)]:
    for k in [K0, K0 * 10, K0 * 1000]:
        M0 = (k / mu + mu) * (k / mu + 1) / mu
        if mu * mp.log(M0) > Ht(k):
            bad += 1
chk("mu*H <= Htilde(k) on grid", bad == 0)
sq = mp.sqrt(1 - d1(K0))
T1 = lambda k: c1 ** 2 * k * sq / mp.sqrt(8) - F(1) / 2
print("T1(K0)=", mp.nstr(T1(K0), 8))
chk("T1(K0) > 22", T1(K0) > 22)
z10 = F("10.169")
chk("0.3 T1(K0) - 3 >= z10", F("0.3") * T1(K0) - 3 >= z10, mp.nstr(F("0.3") * T1(K0) - 3, 8))
chk("0.01 phi_l K0 z10/(1+z10) >= l", F("0.01") * phl * K0 * z10 / (1 + z10) >= ell)
kap0 = (F("0.94") * phl) ** 2 * sq / mp.sqrt(8)
print("kappa0=", mp.nstr(kap0, 8), "0.3kappa0=", mp.nstr(F("0.3") * kap0, 8))
e1 = F("0.3") * c1 ** 2 * (1 + 1 / z10) ** 2 / mp.sqrt(8)
chk("h(z1) endpoint coefficient <= 0.3 kappa0", e1 <= F("0.3") * kap0, f"{mp.nstr(e1,6)} <= {mp.nstr(F('0.3')*kap0,6)}")
cs = F("0.055")
e2 = lambda k: cs * (1 + F("3.15") / (cs * k)) * (1 + 1 / (cs * k)) ** 2
chk("h(z#) endpoint coefficient <= 0.3 kappa0 at K0", e2(K0) <= F("0.3") * kap0, f"{mp.nstr(e2(K0),6)} <= {mp.nstr(F('0.3')*kap0,6)}")
chk("0.3 c1^2/sqrt8 < 0.015 < c#", F("0.3") * c1 ** 2 / mp.sqrt(8) < F("0.015") < cs)
chk("c_p: 0.95 phi_l (1-e^-55) - l/K0 >= 0.819", F("0.95") * phl * (1 - mp.e ** -55) - ell / K0 >= F("0.819"),
    mp.nstr(F("0.95") * phl * (1 - mp.e ** -55) - ell / K0, 8))

# brute force: the implication  [z >= z1(k), h(z) > 0.3 kappa0 k]  =>  z > z#   on a (mu,k,z) grid
viol = 0
for mu in [F("1e-4"), F("0.02"), F("0.1"), F("0.3"), F("0.6"), F(1)]:
    for k in [K0, K0 * F("1.5"), K0 * 3, K0 * 30, K0 * 1000]:
        z1 = F("0.3") * T1(k) - 3
        zs = cs * k / mu ** 2
        hh = lambda z: (z + F("3.15")) * (1 + mu * z) ** 2 / z ** 2
        for i in range(0, 2001):
            z = z1 * (zs / z1) ** (F(i) / 2000)
            if hh(z) > F("0.3") * kap0 * k:
                viol += 1
chk("no z in [z1,z#] with h(z) > 0.3 kappa0 k (grid 6x5x2001)", viol == 0)

# ---- Lemma 16 (crest) at K0 and for all lambda>=K0 (grid)
cp = F("0.819")
b = F("1.02")


def crest(lam):
    Lh = mp.log(8 * lam / 3)
    Hh = 3 * mp.log(lam) - mp.log(K0) + 2 * mp.log(1 + 1 / K0)
    Hl = 2 * mp.log(lam + 1) + 3 / mp.e
    e1h = F("1.1") * mp.sqrt(2) * cp ** (-F(3) / 2) / mp.sqrt(lam)
    dp = 4 * (1 + e1h) * Hl / (cp * lam) + 4 / (theta * cp ** 3 * lam ** 2)
    ac = mp.sqrt(12 / (1 - dp)) + 12 / (theta * cp ** 2 * lam * (1 - dp))
    A1 = mp.sqrt(2) * Hh + ac
    tb1 = A1 / (cp * lam)
    What = cp ** 2 * lam ** 2 * mp.e ** (-tb1) / A1
    vf = b * mp.sqrt(2 * (Hh + Lh) / What)
    return dict(Lh=Lh, Hh=Hh, Hl=Hl, e1h=e1h, dp=dp, ac=ac, A1=A1, tb1=tb1, What=What, vf=vf)


c0 = crest(K0)
print({k: mp.nstr(v, 8) for k, v in c0.items()})
tb1K0 = c0["tb1"]
Ctau = lambda lam, c: 2 * c["Lh"] + c["A1"] + b / mp.sqrt(2) * mp.e ** (tb1K0 / 2) * (c["Hh"] + c["Lh"] + c["A1"])
print("C_tau(K0)=", mp.nstr(Ctau(K0, c0), 8))
chk("C_tau(K0) <= 71.79", Ctau(K0, c0) <= F("71.79"))
cbd = 2 + 3 * mp.sqrt(2) + b / mp.sqrt(2) * mp.e ** (tb1K0 / 2) * (4 + 3 * mp.sqrt(2))
chk("slope c_b <= 12.28", cbd <= F("12.28"), mp.nstr(cbd, 8))
bad = []
prev = None
for i in range(0, 601):
    lam = K0 * F(10) ** (F(i) / 50)  # up to 1e15
    c = crest(lam)
    conds = dict(
        cp_lam_ge10=cp * lam >= 10,
        rise=cp ** 2 * lam >= 2 * c["Lh"],
        dp_lt1=c["dp"] < 1,
        vc_le1=c["ac"] / (cp * lam) <= 1,
        vf_le1=c["vf"] <= 1,
        vf_le_vplus=c["vf"] <= mp.log(1 + c["What"] / lam),
        drop=b * b * (1 - F(2) / 3 * (1 + lam / c["What"]) * c["vf"]) >= 1,
        Cbar=2 * c["Lh"] + c["A1"] + cp * lam * c["vf"] <= F("71.79") + F("12.28") * mp.log(lam / K0),
        Ctau_form=cp * lam * c["vf"] <= b * mp.e ** (c["tb1"] / 2) / mp.sqrt(2) * (c["Hh"] + c["Lh"] + c["A1"]) * (1 + F(10) ** -40),
        tb1_dec=c["tb1"] <= tb1K0,
    )
    for n, v in conds.items():
        if not v:
            bad.append((mp.nstr(lam, 6), n))
chk("crest side conditions for lambda in [1e3,1e15] (601 pts)", not bad, str(bad[:5]))
print("at K0: drop factor=", mp.nstr(b * b * (1 - F(2) / 3 * (1 + K0 / c0['What']) * c0['vf']), 8),
      " v_plus=", mp.nstr(mp.log(1 + c0['What'] / K0), 8))
# elementary inequalities used in Lemma 16
chk("v-1+e^-v >= v^2/3 on [0,1]", all((F(i) / 1000) - 1 + mp.e ** (-F(i) / 1000) >= (F(i) / 1000) ** 2 / 3 for i in range(0, 1001)))
chk("1-(1+v)e^-v >= v^2/2 - v^3/3 on [0,3]", all(1 - (1 + F(i) / 1000) * mp.e ** (-F(i) / 1000) >= (F(i) / 1000) ** 2 / 2 - (F(i) / 1000) ** 3 / 3 for i in range(0, 3001)))
# L0 <= Lhat, H <= Hhat, on a (mu,k) grid
badg = 0
for mu in [F(10) ** (-F(i) / 8) for i in range(0, 40)]:
    for k in [K0, K0 * 2, K0 * 50, K0 * 10 ** 4]:
        lam = k / mu
        cc0 = mu / 2 * (1 - mu * mu / 4)
        M0 = (lam + mu) * (lam + 1) / mu
        c = crest(lam)
        tol_ = F(10) ** -40   # equality is attained at mu=1 (L0=Lhat) and at mu=1,k=K0 (H=Hhat)
        if mp.log(k / cc0) > c["Lh"] + tol_ or mp.log(M0) > c["Hh"] + tol_ or mu * mp.log(M0) > c["Hl"] + tol_:
            badg += 1
chk("L0<=Lhat, H<=Hhat, mu H<=Htilde_lambda on grid", badg == 0)

# ---- Theorem (iv)
cT = cs / F("0.3")
Cbar = lambda lam: F("71.79") + F("12.28") * mp.log(lam / K0)
q = lambda lam: (2 * Cbar(lam) / (cp * mp.sqrt(lam)) + mp.log(1 + 4 * mp.sqrt(lam) + F(16) / 3 * lam)) / (cT * lam)
print("q(K0)=", mp.nstr(q(K0), 8), " parts:", mp.nstr(2 * Cbar(K0) / (cp * mp.sqrt(K0)), 6), mp.nstr(mp.log(1 + 4 * mp.sqrt(K0) + F(16) / 3 * K0), 6))
chk("q(K0) < 0.0772", q(K0) < F("0.0772"))
chk("q decreasing on [1e3,1e15] grid", all(q(K0 * F(10) ** (F(i) / 50)) > q(K0 * F(10) ** (F(i + 1) / 50)) for i in range(0, 600)))
chk("0.75 + 0.0772 < 0.83", F("0.75") + F("0.0772") < F("0.83"))
marg = cT / 4 * K0 - 2 * Ctau(K0, c0) / (cp * mp.sqrt(K0)) - mp.log(1 + 4 * mp.sqrt(K0) + F(16) / 3 * K0)
print("final margin at K0:", mp.nstr(marg, 8))
# alpha<=2, beta<=4/3 for nu=mu/2
chk("alpha<=2, beta<=4/3, 1/beta<2 on mu grid", all((1 / (1 - m / 2) <= 2) and (1 / ((1 - m / 2) * (1 + m / 2)) <= F(4) / 3) for m in [F(i) / 100 for i in range(1, 101)]))

print()
print("FAILS:", fails if fails else "none")
sys.exit(1 if fails else 0)
