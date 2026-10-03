# Lane `dulac` — adversarial referee report (2026-10-03)

Scope: the five statements marked "proved" in `NOTES.md` (Theorem D(a)–(e), Prop. CD, Prop. B1, O1–O2).
The referee tried to refute each one numerically with independent code and then checked every proof line.
Nothing in this report is an interval-arithmetic certificate; "holds" means: every proof step was re-derived by
hand and by an independent exact-algebra script, and no numerical counterexample was found on 21 cycles
(400-bit arithmetic) plus the double-precision side checks listed below.

## 0. Verdicts

| # | statement | verdict | one-line reason |
|---|---|---|---|
| 1 | V>0, Q>0 on every cycle; V' = (η−μ)V + pK, Q' = (η−μ)Q + pW; projection strictly convex and star-shaped about E | **holds** | proof complete; inputs are only v2.0 Thm 3.1, Cor. 3.2 (K>0), Lemma 3.3 (W>0), Thm 4.1; no use of stability or uniqueness; 21/21 cycles pass, including ρ=20 where min Q/(p²+p'²) = 1.2e-20 |
| 2 | λ_c = −⟨pK/V⟩ = −⟨pW/Q⟩; ⟨η₀⟩ = ⟨η₁⟩ = μ; η − η₀ = −pK/V; (η,ζ) affine in κ; compound system for (V,V₃,Q) | **holds** | exact algebra re-derived; three formulas agree to 16 digits on 21 cycles and with a Riccati-free 400-bit monodromy at ρ=12 and ρ=20 |
| 3 | ⟨(αp+βp')(αη+βζ)⟩ < 0 for all (α,β)≠0; branch orderings; Prop. CD | **holds** | proof complete given claim 1 and Thm 3.1 for two phases of one cycle; quadratic form negative definite on 21 cycles; pointwise step verified for 200 random directions |
| 4 | Prop. B1: only R ≤ 0 possible; no Dulac function regular at E; B must blow up at least like r⁻² | **gap** | (i) the cycle-average argument needs a *hyperbolic* attracting cycle, Smith gives only orbital asymptotic stability (λ_c ≤ 0); (ii) rests on unproved brief item 3(b) (disclosed) and the "at least r⁻²" clause is argued only for the pure power B = r^{-2s}. Both are repairable (Sec. 3.4); no counterexample |
| 5 | O1–O2: R = E_Λ(z) − pκ; certificate ⟺ upper bound on κ where p<0 and lower bound where p>0; supporting planes within one cycle give only lower bounds | **holds** (on cycles) | exact identity; F₁∧(q₂−q₁) < 0 and Δ<0 on 1.9 million ordered pairs; the disc version is conditional on brief 3(b), as the lane's own table item 4 says |

Overall: the core by-product (Theorem D, Prop. CD, O1–O2) survived every attack. The only defects are in the
auxiliary Prop. B1. No silent weakening of the lane's target was found: the headline says plainly that no
Dulac function was proved and the conjecture is not proved.

## 1. Independent numerics

All scripts are in `verify/`; none imports a lane script. Python: `/Users/admin/.venvs/analysis-py313/bin/python`.

| script | what it does | output |
|---|---|---|
| `refcyc.py` | own double-precision Newton shooting for the cycle in log variables (U, log Z, log Y), central-difference section derivative | seeds the next script |
| `ref_taylor.py` | 400-bit Taylor-series integrator (order 36, python-flint `arb_series` as a multiprecision float, tolerance 1e-32) for state + Riccati pair; crest polished to residual ≤ 1e-28; every quantity (V, Q, V₃, K, W, O, η₀, η₁, κ) computed from its *definition*; derivatives for the identities taken from the flow's own Taylor coefficients; cycle means by composite Gauss–Legendre (24 nodes per step) and, independently, by termwise series integration | `case_mu*_rho*.json` |
| `run_all.py` | task grid μ ∈ {0.02, 0.3, 1} × ρ ∈ {1.01, 1.5, 4, 12}, plus μ = 5 (Theorem D is claimed for all μ>0) | `ref_taylor.out` |
| `run_extra.py` | μ=1 and 0.7 at ρ=20 (the two cycles on which the lane first saw Q<0), μ=0.005 (ρ=1.01, 3), μ=50 (ρ=2), and a re-run of μ=0.02, ρ=12 | `ref_extra.out`, `ref_extra2.out` |
| `ref_monodromy.py` | Riccati-free λ_c: forward-difference section derivative of the crest map in 400-bit arithmetic (increment 1e-18, tolerance 1e-45) | `ref_monodromy.out` |
| `ref_directions.py` | pointwise step of D(d) for 40 random directions × 49 levels on 5 cycles (double precision, whlib) | `ref_directions.out` |
| `ref_disc.py` | sanity check of the conditional disc statements on W^u(E) spirals (double precision, whlib) | `ref_disc.out` |
| `ref_sym.py` | referee's own exact (sympy) derivation of every identity, with the explicit cofactor of the constraint O | `ref_sym.out` |
| `summarize.py` | table below | `summary.out` |

Commands: `python refcyc.py`; `research-run run --cpu 1 --memory-gb 1 --label wh-dulac-referee -- python run_all.py`
(≈ 17 min, of which 11.5 min is μ=0.02, ρ=12); `python run_extra.py [mu,rho ...]`; `python ref_monodromy.py 1 12 5.5`;
`python ref_monodromy.py 1 20 6.5`; `python ref_directions.py`; `python ref_disc.py`; `python ref_sym.py`; `python summarize.py`.

### 1.1 Results on 21 cycles (400-bit)

    mu   rho          T      Umin   lambda_c=<eta>-mu   <eta0>-mu <eta1>-mu   minV/r2   minQ/s2    <eta p>  <zeta pd>       det   idmax
 0.005  1.01   104.7413     -3.11  -4.31850813105e-05     +5e-37    +1e-36  1.55e-03  2.03e-03 -9.408e-05 -4.910e-05  2.34e-09   3e-33
 0.005     3  1059.2925  -1587.99  -1.59913117761e-03     +2e-35    +4e-28  5.99e-06  7.38e-07 -1.772e-02 -3.782e-02  4.06e-04   2e-33
  0.02  1.01    46.4096     -1.31  -1.81932169554e-04     +1e-39    +2e-39  1.14e-02  1.41e-02 -3.910e-04 -2.155e-04  4.38e-08   2e-36
  0.02   1.5   130.0667    -53.44  -3.25990023299e-03     -2e-33    +3e-30  5.09e-04  2.00e-04 -1.606e-02 -1.615e-02  1.46e-04   2e-34
  0.02     4   343.0578   -793.11  -7.72365255768e-03     +3e-34    -2e-28  5.20e-05  3.32e-06 -1.127e-01 -3.306e-01  2.36e-02   9e-34
  0.02    12   877.7237  -8356.54  -1.28112922004e-02     -4e-34    -9e-30  4.08e-06  1.80e-10 -4.831e-01 -4.405e+00  1.54e+00   8e-33
   0.3  1.01    11.5256     -0.48  -1.94228053127e-03     +3e-37    +3e-37  2.50e-01  2.62e-01 -5.996e-03 -6.895e-03  2.61e-05   2e-36
   0.3   1.5    14.7034     -7.46  -5.63597952746e-02     +1e-33    +2e-31  6.30e-02  4.40e-02 -2.830e-01 -4.409e-01  8.12e-02   1e-35
   0.3     4    29.7593    -88.62  -1.26570339699e-01     +1e-33    +9e-28  7.21e-03  3.67e-04 -1.905e+00 -7.936e+00  1.04e+01   2e-35
   0.3    12    74.6786   -926.10  -2.03960370317e-01     -4e-33    -8e-27  5.00e-04  9.76e-10 -8.180e+00 -1.067e+02  6.64e+02   9e-34
   0.7    20    67.1585  -1986.25  -5.57643127907e-01     +5e-34    -7e-30  4.12e-04  8.83e-19 -3.955e+01 -1.314e+03  4.29e+04   2e-33
     1  1.01     6.3042     -0.44  -3.98399812506e-03     +7e-38    +1e-37  8.48e-01  8.68e-01 -2.003e-02 -6.019e-02  8.05e-04   2e-36
     1   1.5     7.4907     -6.03  -1.46389340857e-01     +2e-43    +3e-44  2.69e-01  2.01e-01 -1.033e+00 -3.698e+00  2.68e+00   5e-38
     1     4    13.8416    -64.18  -3.95940343043e-01     -2e-33    +4e-27  3.31e-02  8.32e-04 -7.390e+00 -6.237e+01  3.42e+02   8e-35
     1    12    34.4883   -668.28  -6.74126202215e-01     +2e-33    +4e-28  2.07e-03  5.14e-12 -3.307e+01 -8.412e+02  2.24e+04   7e-34
     1    20    55.3084  -1927.51  -7.70636848475e-01     +5e-33    -3e-28  5.97e-04  1.16e-20 -6.029e+01 -2.599e+03  1.31e+05   2e-33
     5  1.01     2.8269     -0.52  -7.22360679253e-03     -1e-39    -2e-39  9.99e-01  9.95e-01 -9.975e-02 -2.309e+00  1.40e-01   4e-38
     5   1.5     3.8374     -8.89  -1.85190037894e-01     +4e-32    +4e-31  7.31e-01  5.14e-01 -4.515e+00 -1.520e+02  4.30e+02   4e-35
     5     4     8.2029   -112.14  -4.11342699859e-01     +7e-33    +2e-26  9.05e-02  5.48e-03 -3.051e+01 -2.827e+03  5.80e+04   1e-35
     5    12    20.6734  -1175.44  -6.67883678728e-01     -5e-33    +1e-25  6.75e-03  4.16e-08 -1.309e+02 -3.794e+04  3.72e+06   1e-33
    50     2     3.7012   -141.69  -2.32751224508e-01     +1e-32    +1e-26  5.05e-01  1.43e-01 -8.316e+01 -2.935e+05  1.42e+07   6e-34

(minV/r2 = min V/(U²+p²), minQ/s2 = min Q/(p²+p'²) over ≥ 1 900 samples per cycle, up to 265 000;
det = ⟨ηp⟩⟨ζp'⟩ − ⟨ζp⟩²; idmax = largest relative residual of all pointwise identities.)

On every one of the 21 cycles:

- **V>0, Q>0, K>0, W>0** at every sample (sign decided in ball arithmetic, not from floats).
  At ρ=20 the minimum of Q/(p²+p'²) is 1.16e-20 (μ=1) and 8.83e-19 (μ=0.7): the lane's values are reproduced
  to three digits, so its "numerical trap" fix was correct and Q is genuinely positive there.
- **λ_c**: −⟨pK/V⟩ and −⟨pW/Q⟩ equal ⟨η⟩−μ to all 16 digits stored; Gauss–Legendre and termwise-series
  quadrature of the two rational integrands agree to ≤ 1e-15 relative. Where the double-precision monodromy is
  usable (m_t ≥ 1e-5, 16 cycles) it agrees within its finite-difference noise.
- **Riccati-free check at large ρ** (`ref_monodromy.out`): μ=1, ρ=12: log(m_t)/T = −0.6741262022145094;
  μ=1, ρ=20: −0.7706368484748485 (m_t = 3.08e-19); both equal ⟨η⟩−μ, −⟨pK/V⟩, −⟨pW/Q⟩ to 16 digits; the section
  derivative is entrywise positive and log det M = −AT to 1e-13.
- **⟨η₀⟩−μ, ⟨η₁⟩−μ** vanish to ≤ 1e-25. η₀ is negative on part of the cycle for μ = 0.005 and 0.02 (every ρ tested)
  and for μ = 0.3 at ρ ≥ 4, as the lane notes; so ω = (V,V₃,Q) is indeed not componentwise positive in general.
- **Identities** V' = (η−μ)V + pK, Q' = (η−μ)Q + pW, V' = −μV + V₃, V₃' = −V₃ + Q + kψp, Q' = xV − AQ + kψp',
  κ' = kψ/V − (1+η+η₀−μ)κ, η = η₀ − pκ, η = η₁ − pW/Q, ζ = ζ₀ − p'κ, V₃ = pK + ηV, W = W₀ + κQ: relative
  residual ≤ 1e-32 (the size of the tangent-orthogonality defect O of the relaxed Riccati pair).
- **Sign pattern** (η−η₀)p < 0, (η−η₁)p < 0: zero violations. In deep troughs the direct difference is below
  400-bit resolution (|pW/Q| ~ e^{-8000}·η); those samples are counted as unresolved (μ=0.02, ρ=12:
  422 048 resolved, 108 832 unresolved, 0 violations) and are covered by K/V>0, W/Q>0 plus the identity.
- **Monotone integrals**: ⟨ηp⟩ < 0, ⟨ζp'⟩ < 0, det > 0, ⟨ηp'⟩ = ⟨ζp⟩ to 1e-12 relative.
- **Branch orderings** at 167 exactly matched U-levels and 167 exactly matched p-levels per cycle (levels
  solved to 1e-48 on the Taylor polynomials, clustered at crest, trough and the extrema of p): η_low > η_up,
  y_low > y_up, η_up Δp < y_low − y_up < η_low Δp; U_r > U_l, y_r > y_l, ζ_r > ζ_l, ζ_l ΔU < y_r − y_l < ζ_r ΔU.
  All strictly positive; smallest relative margin 1.4e-34 (μ=50), resolved at 400 bits.
- **Pairs** (300 points per cycle, 89 700 ordered pairs, 1.9 million in total): Δ(base A; B) < 0 (Thm 3.1 for two
  phases of one cycle), F_A∧(q_B−q_A) < 0 (strict convexity, clockwise), gradient monotonicity
  (ζ_A−ζ_B)(U_A−U_B) − (η_A−η_B)(p_A−p_B) > 0, and the O2 identity Δ = [ζ₀ΔU − η₀Δp − Δy] + κ_A F_A∧Δq to 1e-32.
- **Unimodality inputs**: exactly one interior sign change of p and two of p' per period; phase order
  crest < p_min < trough < p_max < crest.
- Moment identities of v2.0 (⟨x⟩ = k, ⟨p⟩ = 0, ⟨p'²⟩ = μ⟨p²⟩, ⟨xp²⟩ = A⟨p'²⟩) hold to 5e-34 and the crest closes to
  ≤ 2e-34, so the cycles are cycles.
- O2 side remarks: y > ζ on all 21 cycles; η < λ − p holds on all 10 cycles with μ ≤ 1/2 and *fails* on
  (μ,ρ) = (0.7,20), (1,4), (1,12), (1,20), (5,1.5), (5,4), (5,12), (50,2), so the lane's restriction μ ≤ 1/2 is not
  cosmetic and is correctly stated.
- **Refutation of an unproved side remark** (NOTES table item 10, "K/V < 0.6 on D, numerically supported, not
  used"): max κ on the cycle is 0.627 at (0.005, 3), 0.607 at (0.005, 1.01), 0.606 at (0.02, 12). The bound 0.6 is
  false for small μ. Max W/Q over the 21 cycles is 1.118 (< 1.15). Item 10 is not among the proved claims.

`ref_directions.out`: for 40 random directions × 49 levels on each of 5 cycles, ℓ = αU + βp has exactly two
critical points, each level has exactly two points, s < 0 and (αη+βζ)_A < (αη+βζ)_B in every instance.

`ref_disc.out` (conditional statements only): on W^u(E) spirals for 5 parameter sets, V>0, Q>0,
η = η₀ − pκ to 1e-8, and η − μ → 2σ at E (e.g. μ=0.02, ρ=1.5: η−μ ∈ [0.009441, 0.009452], 2σ = 0.009447).
The entry max|O| = 1 for (0.02, 1.5) is a whlib artifact (absolute tolerance 1e-16 against ζ ~ 1e-20 in the trough).

No numerical counterexample to any of the five statements was found.

## 2. Exact algebra (`ref_sym.out`)

Written from the statements, with x an independent symbol (D x = xp) so that every residual is polynomial.
All checks pass (28 identity lines plus the wedge and equilibrium lines). Cofactors of the constraint O = ζp − ηp' − (x−y):
V' − [(η−μ)V + pK] = −U·O; V₃ − [pK + ηV] = −U·O; Q' − [(η−μ)Q + pW] = (ηp + k − μp + p − y)·O;
ηV − (V₃ − pK) = U·O; ζV − (Q − p'K) = p·O; L₂[p] = O; L₂[y−k] − W = (η−μ+1)·O.
Identically zero (no constraint needed): V' = −μV + V₃, Q = p'(k−y) + p(x−y), Q' = xV − AQ + kψp',
V₃' = −V₃ + Q + kψp, the κ-equation, the Riccati defects kψ(p,p')/V of (η₀,ζ₀), the wedge identification
(V, V₃, Q) = (−(12), (13), (23)) components of r∧z', W(η₀−pκ, ζ₀−p'κ) = W₀ + κQ, the two O2 plane facts
(the κ=0 plane contains z' and E), (y−ζ)' + (y−ζ) = ηζ, (λ−p−η)' + μ(λ−p−η) = (y−ζ) + η(η+1−2μ),
γ_E − A = 2σ, and 𝔻Q₀ = 2σQ₀ for the planar linearisation. These agree with the lane's `sym_check.out`.

## 3. Proof audit, line by line

### 3.1 Theorem D(a),(b) (claim 1) — holds
- (b): V' = pp' − Up''; inserting p'' = (η−μ)p' − ζp gives V' − (η−μ)V = p[p' + μp + ζU − ηp] = pK. Correct.
  Q' = p'p'' − pp'''; with p''' = (η−μ)p'' + η'p' − ζ'p − ζp' the remainder is p(ζ'p − η'p') = pW. Correct.
- (a): at a zero of p, O = 0 gives x − k = (1+η)(y − k), so sign(y−k) = sign(U) and U ≠ 0 (else the point is E);
  hence V = U(y−k) > 0, Q = (y−k)² > 0 there. Between the trough t_a and crest t_b, p ≥ 0 and K>0 make V/E
  nondecreasing; from t_b to t_a+T, p ≤ 0 makes it nonincreasing, and V(t_a+T) = V(t_a) > 0 closes the argument.
  The only inputs are K>0 (v2.0 Cor. 3.2), W>0 (Lemma 3.3) and "p has exactly two simple zeros, trough U<0,
  crest U>0" (Thm 4.1). No compactness, no uniformity, no stability, no uniqueness is used. Correct.
- Curvature (pp'' − p'²)/(p²+p'²)^{3/2} = −Q/(…) < 0 on a regular simple closed C² curve gives strict convexity
  (standard); d arg(U+ip)/dt = −V/(U²+p²) < 0 gives strict star-shape. Correct.
- Hypotheses of imported results: v2.0 Thm 3.1 is stated for a cycle base and any distinct log-bounded entire
  solution, "distinct phases of the same cycle are allowed"; its proof (re-derived: n' + Jᵀn + (1+η)n = 0) does
  not use k > k_H or stability. Satisfied.

### 3.2 Theorem D(c),(e) (claim 2) — holds
- 0 = ⟨(log V)'⟩ = ⟨η−μ⟩ + ⟨pK/V⟩ needs only V>0 and periodicity; λ_c = ⟨η⟩−μ is v2.0 Prop. 2.1. Correct.
- (log V)' = η₀ − μ and (log Q)' = η₁ − μ are state identities, so ⟨η₀⟩ = ⟨η₁⟩ = μ exactly. Correct.
- The 2×2 system {ζU − ηp = K + H, ζp − ηp' = x − y} has determinant ±V ≠ 0; solving gives η = (V₃ − pK)/V,
  ζ = (Q − p'K)/V. Correct. κ' follows from K' + (1+η)K = kψ and V' = (η₀−μ)V; it is an exact identity.
- (e): rows checked by hand and by `ref_sym.py`; the wedge sign convention matches v2.0 (4.6). Correct.

### 3.3 Theorem D(d), Prop. CD (claim 3) — holds
- For ℓ = αU + βp, ℓ' = (α,β)·(p,p') vanishes exactly twice per period because the tangent angle is strictly
  monotone with total turning 2π (claim 1). Two points A, B per intermediate level, q_A − q_B = s(β,−α), s ≠ 0;
  the set of (direction, level) pairs is connected and s is continuous, so sign s is constant; s < 0 at (1,0).
  Adding Thm 3.1 with base A and with base B gives s[(αη+βζ)_A − (αη+βζ)_B] > 0. Then
  ∮(αη+βζ)dℓ = ∫[g_A − g_B]dℓ < 0. Correct; ⟨ηp'⟩ = ⟨ζp⟩ from ηp' = ζp − y'. Correct.
- Prop. CD: λ_cT = ∮F_U dU + ∮F_p dp + ∮(3)dt is the chain rule plus periodicity of F along the cycle; the first
  two integrals are ≤ 0 by (1),(2) with the orientation checked (upper branch dU>0; trough-side branch dp>0).
  Correct. Minor wording point: the parenthetical "true whenever F_U is nondecreasing in η and y, …" needs F
  defined and monotone on a set containing the coordinate-wise path between the two compared points, not merely
  "near the cycle". It is not used in any conclusion.
- No circularity: neither statement uses λ_c < 0, uniqueness, or the existence of a Dulac function.

### 3.4 Prop. B1 (claim 4) — gap
1. **B1(i), undisclosed gap.** The argument is ⟨R⟩ = λ_c on any cycle, "an attracting cycle exists (Smith)",
   hence R ≥ 0, R ≢ 0 is impossible. The imported theorem (v2.0 Sec. 1.3: "at least one of which is orbitally
   asymptotically stable") gives λ_c ≤ 0, not λ_c < 0. If that cycle is neutral (λ_c = 0, which is exactly the
   case the conjecture has not excluded), R ≥ 0 forces only R ≡ 0 on that cycle and does not contradict R ≢ 0 on
   D∖{E}. What the written argument proves: no Λ has R ≥ 0 with R ≢ 0 *on the Smith cycle*; in particular R > 0 is
   impossible, and the full statement holds at every parameter where the stable cycle is hyperbolic.
   Repair (needs brief 3(b), which the lane says (i) does not need): on the planar flow take a transversal Σ at
   the asymptotically stable cycle γ, a point q ∈ Σ close to γ, and the region Ω between γ and the closed curve
   {orbit arc q → P(q)} ∪ [P(q), q]. Then ∬_Ω B·R = flux of Bf through [P(q), q] < 0, so R ≥ 0 on a one-sided
   neighbourhood of γ is impossible whether or not γ is hyperbolic.
2. **B1(ii), disclosed dependence.** R(E) = η_E − μ = γ_E − A = 2σ is correct algebra (checked), and η − μ → 2σ
   on W^u(E) is confirmed numerically. The conclusion "no Λ that is C¹ near E works" needs η to be defined and
   continuous on D at E, i.e. brief item 3(b), which is unproved. The lane states this dependence.
3. **B1(ii), under-argued clause.** "B must blow up at E at least like r⁻²" is justified only by the example
   B ~ r^{-2s}. A general two-line argument exists: along an orbit of W^u(E), R ≤ 0 gives
   d/dt log B ≤ −(η−μ), and η − μ − 2σ = O(r) is integrable at t → −∞, so
   log B(z(t)) ≥ log B(z(t₀)) + 2σ(t₀−t) − C, i.e. B ≥ c·r⁻². This also covers "any B vanishing at E".
4. Reporting detail: "best margin exactly −2.000σ for exp(poly₃(U,p))" — `byproduct.out` shows −2.073σ for
   (μ,ρ) = (0.05, 5); consistent with "≤ −2σ", not with "exactly".

No numerical counterexample to B1; all 21 cycles computed here have λ_c < 0.

### 3.5 O1–O2 (claim 5) — holds on cycles
- R = 𝔻Λ + η − μ = (𝔻Λ + η₀ − μ) − pκ is claim 2. The equivalence with one-sided bounds on κ is immediate.
  W = W₀(z) + κQ and K = κV are exact, so (z,κ) is closed under the flow. Correct.
- O2: Δ = [ζ₀ΔU − η₀Δp − Δy] + κ₁F₁∧(q₂−q₁) (checked exactly and on 1.9 million pairs), F₁∧(q₂−q₁) < 0 by
  strict convexity, so Δ < 0 is a lower bound on κ₁. Each supporting-plane inequality involves κ only at its
  base point, so no combination of them bounds κ from above. Correct.
- On D∖{E} the same statements need V>0 on the disc and O = 0 on the attractor, i.e. brief 3(b). The referee
  checked the lane's sketch (at p = 0, ζ_E U < H < ζU gives UH > 0; the monotonicity argument for V/E does not
  need periodicity): it is correct given 3(b). The claim as worded in the lane summary does not repeat this
  condition; NOTES table item 4 does.

## 4. Other observations

- `verify_theorem.out` (lane): the monodromy column disagrees with ⟨η⟩−μ at ρ=20 (μ=0.3: −0.246 vs −0.232;
  μ=0.7: −0.308 vs −0.558; μ=1: −0.564 vs −0.771; μ=0.1: NaN) and in the sixth digit at ρ=8. NOTES Sec. 2 reports
  only the three Riccati-based values. The cause is double-precision loss in whlib's variational return map when
  m_t ≲ 1e-12; the referee's 400-bit monodromy removes the doubt (Sec. 1.1). The statement is true; the
  disclosure was incomplete.
- Near-Hopf asymptotics (NOTES Sec. 2, "numerically supported"): at μ=1, ρ=1.01 the referee's λ_c = −3.984e-3 and
  ⟨ηp⟩ = −2.003e-2 are consistent with −ε²/40 and −ε²/8 for ε ≈ 0.40 (mean of crest 0.369 and trough 0.435).
- The LP results of NOTES Sec. 3 (existence of per-parameter Dulac functions in a 26-feature family) are
  numerical claims and were not re-run; they are not among the statements marked proved.
- Referee-side corrections during this review (recorded so that nobody mistakes them for failures):
  the first grid pass judged signs from floats that underflow in deep troughs and flagged μ=0.02, ρ=12 and
  μ=0.005, ρ=3; both were re-run with exact sign tests and pass (`ref_extra2.out`). The first-pass double
  monodromy tolerance was also too tight for m_t ~ 1e-5 (finite-difference noise).

## 5. Limits of this review

- Floating/multiprecision evidence only; no interval enclosure of the cycles.
- All computed cycles are attracting; Theorem D for neutral or repelling cycles rests on the proof alone
  (which does not use stability).
- v2.0 results (Thm 3.1, Cor. 3.2, Lemma 3.3, Thm 4.1, Prop. 2.1) were taken as given; their identities were
  re-derived and their conclusions (Δ<0, K>0, W>0, two zeros of p) were confirmed numerically on every cycle.
- Brief item 3(b) was not verified.
