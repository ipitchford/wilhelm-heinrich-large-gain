#!/usr/bin/env python3
"""Referee check of Theorem A (NOTES 3.2).  Independent of the lane's toda.py:
level integrals are computed by tanh-sinh quadrature in U (mpmath, 30 digits), not by ODE integration.

  T0 = 2 int dU/P,  J0 = 2 int P dU,  J1 = 2 int P^3 dU,   P = sqrt(2(h - Phi(U))),  U in [U_-, U_+].
  <P^2> = J0/T0, <P^4> = J1/T0, kurt = <P^4>/<P^2>^2, Rsel = (2/3) J1/J0.

Tests:
 (T1) polynomial identity used in (ii)
 (T2) Toda potential, h in [1e-6, 1e6]: kurt <= 9/5, (ii) first inequality, (i), 4h/5 <= Rsel <= 4h/3,
      J0' = T0, J1' = 3 J0 and Rsel' = (2/3)(3-kurt) by central differences, Rsel' >= 4/5
 (T3) other strictly convex potentials (try to break 9/5)
 (T4) power check: a NON-convex potential must violate the bound (shows the test can detect failure)
 (T5) monotonicity of the velocity time-density f(P) on a grid (the single structural step of the proof)
"""
import mpmath as mp
import sympy as sp

mp.mp.dps = 30
FAIL = []


def turning(Phi, h, side):
    # bracket by doubling, then bisection + secant polish
    a, b = mp.mpf(0), mp.mpf(side)
    while Phi(b) < h:
        b *= 2
    for _ in range(200):
        m = (a + b) / 2
        if Phi(m) < h:
            a = m
        else:
            b = m
    return (a + b) / 2


def level(Phi, h):
    h = mp.mpf(h)
    um, up = turning(Phi, h, -1), turning(Phi, h, +1)
    P = lambda u: mp.sqrt(max(2 * (h - Phi(u)), mp.mpf(0)))
    def integ(g):
        return 2 * (mp.quad(g, mp.linspace(um, 0, 9)) + mp.quad(g, mp.linspace(0, up, 9)))
    T0 = integ(lambda u: 1 / P(u) if P(u) > 0 else mp.mpf(0))
    J0 = integ(P)
    J1 = integ(lambda u: P(u) ** 3)
    return dict(h=h, um=um, up=up, T0=T0, J0=J0, J1=J1, p2=J0 / T0, p4=J1 / T0,
                kurt=(J1 / T0) / (J0 / T0) ** 2, Rsel=mp.mpf(2) / 3 * J1 / J0)


# (T1) ---------------------------------------------------------------------------------------------------------
t = sp.symbols('t')
m2 = (1 + t + t**2) / 3
m4 = (1 + t + t**2 + t**3 + t**4) / 5
ident = sp.expand(sp.Rational(6, 5) * m2 - sp.Rational(1, 5) - m4 - t * (1 + t) * (1 - t**2) / 5)
g = sp.symbols('g')
ident2 = sp.expand(sp.Rational(9, 5) * g**2 - sp.Rational(6, 5) * g + sp.Rational(1, 5) - (3 * g - 1)**2 / 5)
print("(T1) (6/5)m2 - 1/5 - m4 - t(1+t)(1-t^2)/5 =", ident, ";  (9/5)g^2-(6/5)g+1/5-(3g-1)^2/5 =", ident2)
if ident != 0 or ident2 != 0:
    FAIL.append("T1")

# (T2) ---------------------------------------------------------------------------------------------------------
toda = lambda u: mp.expm1(u) - u
print("\n(T2) Toda potential e^U-U-1")
print("# h        kurt        9/5-kurt    (ii)slack      <P^2>/(Pm^2/3)  Rsel/h    Rsel'(fd)   (2/3)(3-kurt)  J0'/T0-1   J1'/(3J0)-1")
worst = mp.mpf(0)
minRp = mp.inf
for h in ['1e-6', '1e-4', '1e-2', '0.1', '0.5', '1', '3', '10', '30', '100', '1e3', '1e4', '1e5', '1e6']:
    d = level(toda, h)
    hh = d['h']
    dh = hh * mp.mpf('1e-5')
    dp, dm = level(toda, hh + dh), level(toda, hh - dh)
    Rp_fd = (dp['Rsel'] - dm['Rsel']) / (2 * dh)
    Rp_formula = mp.mpf(2) / 3 * (3 - d['kurt'])
    J0p = (dp['J0'] - dm['J0']) / (2 * dh)
    J1p = (dp['J1'] - dm['J1']) / (2 * dh)
    Pm2 = 2 * hh
    slack2 = mp.mpf(6) / 5 * Pm2 * d['p2'] - Pm2 ** 2 / 5 - d['p4']        # (ii) first inequality, must be >= 0
    print(f"{mp.nstr(hh,6):<9} {mp.nstr(d['kurt'],10):<11} {mp.nstr(mp.mpf(9)/5-d['kurt'],4):<11} "
          f"{mp.nstr(slack2/Pm2**2,4):<14} {mp.nstr(d['p2']/(Pm2/3),8):<15} {mp.nstr(d['Rsel']/hh,7):<9} "
          f"{mp.nstr(Rp_fd,9):<11} {mp.nstr(Rp_formula,9):<14} {mp.nstr(J0p/d['T0']-1,2):<10} {mp.nstr(J1p/(3*d['J0'])-1,2)}")
    minRp = min(minRp, Rp_fd, Rp_formula)
    if d['kurt'] > mp.mpf(9) / 5 or slack2 < 0 or d['p2'] < Pm2 / 3 or not (mp.mpf(4) / 5 <= d['Rsel'] / hh <= mp.mpf(4) / 3):
        FAIL.append(f"T2 h={h}")
    if abs(Rp_fd - Rp_formula) > 1e-6 or abs(J0p / d['T0'] - 1) > 1e-6 or abs(J1p / (3 * d['J0']) - 1) > 1e-6:
        FAIL.append(f"T2-derivative h={h}")
print("min Rsel' over grid =", mp.nstr(minRp, 10), " (claim: >= 0.8)")
if minRp < mp.mpf(4) / 5:
    FAIL.append("T2 Rsel' < 4/5")

# fine monotonicity scan of Rsel (bijection claim) on 400 log-spaced h
prev = None
mono = True
for i in range(0, 121):
    h = mp.mpf(10) ** (mp.mpf(-6) + mp.mpf(i) / 10)
    R = level(toda, h)['Rsel']
    if prev is not None and R <= prev:
        mono = False
    prev = R
print("Rsel strictly increasing on 121 log-spaced h in [1e-6,1e6]:", mono)
if not mono:
    FAIL.append("T2 monotone")

# (T3) ---------------------------------------------------------------------------------------------------------
pots = {
    "U^2/2": lambda u: u * u / 2,
    "U^4/4": lambda u: u ** 4 / 4,
    "cosh U - 1": lambda u: mp.cosh(u) - 1,
    "sqrt(1+U^2)-1": lambda u: mp.sqrt(1 + u * u) - 1,
    "log cosh U": lambda u: mp.log(mp.cosh(u)) if abs(u) < 30 else abs(u) - mp.log(2),
    "(e^{3U}-3U-1)/9": lambda u: (mp.expm1(3 * u) - 3 * u) / 9,
    "U^2/2 for U<0, 50U^2 for U>0 (C^1 only)": lambda u: u * u / 2 if u < 0 else 50 * u * u,
    "e^U-U-1 + U^6": lambda u: mp.expm1(u) - u + u ** 6,
    "|U|^{2.5}": lambda u: abs(u) ** mp.mpf('2.5'),
    "|U|^{1.2} (convex, C^1 only)": lambda u: abs(u) ** mp.mpf('1.2'),
}
print("\n(T3) other convex potentials: max kurtosis over h in {1e-3,...,1e4}")
for name, Phi in pots.items():
    ks = []
    for h in ['1e-3', '1e-1', '1', '10', '1e2', '1e4']:
        if name == "sqrt(1+U^2)-1" or name == "log cosh U":
            pass
        d = level(Phi, h)
        ks.append(d['kurt'])
    mk = max(ks)
    print(f"  {name:<42} kurt: " + " ".join(mp.nstr(x, 7) for x in ks) + f"   max={mp.nstr(mk,8)}  ok={mk <= mp.mpf(9)/5 + mp.mpf('1e-12')}")
    if mk > mp.mpf(9) / 5 + mp.mpf('1e-12'):
        FAIL.append("T3 " + name)

# (T4) ---------------------------------------------------------------------------------------------------------
print("\n(T4) power check with NON-convex potentials (the bound should fail somewhere)")
nonconv = {
    "1-cos U (h near separatrix 2)": (lambda u: 1 - mp.cos(u), ['1.9', '1.99', '1.999']),
    "U^2/(1+U^2)*... log(1+U^2) (concave tails)": (lambda u: mp.log(1 + u * u), ['1', '5', '10']),
}
for name, (Phi, hs) in nonconv.items():
    ks = [level(Phi, h)['kurt'] for h in hs]
    print(f"  {name:<45} kurt: " + " ".join(mp.nstr(x, 7) for x in ks), "  -> exceeds 9/5:", any(k > mp.mpf(9) / 5 for k in ks))

# (T5) ---------------------------------------------------------------------------------------------------------
print("\n(T5) density f(P) = 1/|Phi'(U_-)| + 1/Phi'(U_+) at phi = h - P^2/2: increasing in P?  (Toda)")
for h in ['0.01', '1', '100']:
    hh = mp.mpf(h)
    Pm = mp.sqrt(2 * hh)
    prev = None
    ok = True
    for i in range(1, 200):
        P = Pm * i / 200
        phi = hh - P * P / 2
        um, up = turning(toda, phi, -1), turning(toda, phi, +1)
        f = 1 / abs(mp.expm1(um)) + 1 / mp.expm1(up)
        if prev is not None and f <= prev:
            ok = False
        prev = f
    print(f"  h={h}: f strictly increasing on 199 grid points: {ok}")
    if not ok:
        FAIL.append("T5")

print("\nFAILURES:", FAIL if FAIL else "none")
