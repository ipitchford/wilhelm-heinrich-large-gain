# Cartography lane — NOTES (2026-10-03)

Pure numerics. No statement below is a proof except the one-line identity in §4.2, which is marked.
All numbers are ordinary floating point (numba DOP853, rtol 1e-12, events located by Newton on the step length).
File and column descriptions: `README.md`.

**Bottom line.** No second cycle and no cycle with lambda_c ≥ 0 was found anywhere (300 grid cycles, 827 branch
points, 240 + 240 random parameter points). Status `partial`. Of the candidate structural statements,
**S1 and S4 survive every test**; **S2 (all six coordinates), S3b and one-signedness in S5 are refuted**;
**S3a survives except in a band around mu ≈ 1 with rho ≳ 4.5, where it is refuted on the crest side**.

---

## 1. Tools (what was built)

- `cart.py`: compiled DOP853 with per-component tolerances (U is controlled absolutely because troughs reach
  U ≈ −2·10⁶; eta, zeta, K, W relatively). Crest return map with 2×2 derivative M and d/dk.
- Cycle finder that works through Hopf at small mu: **amplitude-fixed shooting** (crest U_max = h prescribed,
  unknowns y0 and k), continued in h from the Hopf point, then a secant on h to hit each target rho.
  This is well conditioned where fixed-k Newton is not (m_t → 1), and it yields rho(h) along the whole primary
  branch, so a fold of the branch would be seen directly as non-monotone rho(h).
- log m_t is taken from the pulled-back Riccati quadrature T(⟨eta⟩ − mu), which stays accurate when
  m_t ~ e^{−500}; the eigenvalue of M is used as a cross-check where m_t > 1e-7.
- Self-test: `python test_core.py` reproduces the four bundle receipts (T to 7 digits, m_t to 6–9 digits; the receipts carry a rounded k).
- Tolerance tightening (`python tol_check.py`, rt = 1e-10, 1e-12, 1e-13): the weakest exponents, the
  non-monotone return-map ratio and the S3a counterexample are unchanged to ≥ 8 significant digits.

Two numerical traps found and fixed (both would have produced false alarms):
1. `K' = −(1+eta)K + k(Ue^U − e^U + 1)`: the source cancels catastrophically near E; a series is used for |U| < 0.5.
2. For mu ≫ 1 in direct coordinates the explicit integrator sits at its stability limit in troughs and p wiggles
   around its quasi-static value; near the trough bottom this produced one spurious "crest" (a pair
   U_n = 12.6 → U_{n+1} = −670 in the first hunt run). Crest events are now accepted only at U > 0 and trough
   events only at U < 0 (on the attractor crest ⇔ U > 0, by tangent orthogonality). The first hunt run is kept
   in `superseded/`.

---

## 2. Part 1 — map of the weak exponent

Command: `python grid.py` (≈ 100 s). Output `grid_table.csv` (12 mu × 25 rho = 300 cycles,
mu ∈ {0.002 … 1}, rho − 1 ∈ [1e-3, 199] log-spaced) and `branch_table.csv` (827 amplitude-branch points).

| question | answer (numerical) |
|---|---|
| all cycles located? | 300/300; Newton residual ≤ 9.6e-11; ⟨x⟩ = k to 4.3e-6; moment identity k − k_H = (A²+mu)⟨q⟩ + 2A⟨q²⟩ + ⟨q³⟩ to 4.3e-6 relative (worst at rho = 200) |
| is lambda_c ever ≥ 0? | **No.** max lambda_c = −1.93e-6 (mu = 0.002, rho − 1 = 1e-3, log m_t = −2.829e-4, stable under tolerance change to 10 digits). On branches the maximum is −2.0e-10 at h = 0.01 (the Hopf end) |
| is log m_t monotone in rho at fixed mu? | **Yes, strictly decreasing** at all 12 mu on the grid, and strictly decreasing in amplitude along all 12 branches; lambda_c is also strictly decreasing in rho |
| fold of the primary branch? | **None**: rho(h) strictly increasing on all branches (and on the branches of all 240 hunt points); min d log(rho−1)/dh = 0.503 (→ ½ at large h) |
| how close to neutral? | \|log m_t\| ≥ 0.548 × the Hopf-linear value 2πσ(rho−1)/(1+σ²), σ = (1+mu)/√mu, over the whole grid (minimum at mu = 0.2, rho − 1 = 3.4); ≥ 0.556 over the hunt points. m_t > 0.999 occurs only for rho − 1 < 0.0035 and m_t > 0.99 only for rho − 1 < 0.022 (grid and hunt points together). **No near-neutral behaviour away from Hopf.** |
| ⟨U⟩ + log rho | < 0 on every cycle (max −3.0e-3): brief item (g) holds |
| large rho | lambda_c → −0.955 mu … −0.969 mu at rho = 200 (⟨eta⟩ ≈ 0.04 mu); U_max ≈ 2 log rho + c(mu); U_min down to −1.7·10⁶; log m_t ≈ −1.26 rho (mu → 0) to −2.54 rho (mu = 1) |
| near Hopf | log m_t / [−2πσ(rho−1)/(1+σ²)] → 1 (0.9999 at mu = 1, rho − 1 = 1e-3); the Hopf-linear regime requires rho − 1 ≪ mu |

log m_t = lambda_c·T (selected rows; full table in `grid_table.csv`):

| rho − 1 | mu=0.002 | 0.01 | 0.05 | 0.2 | 0.5 | 1 |
|---|---|---|---|---|---|---|
| 1.0e-3 | −2.83e-4 | −6.16e-4 | −1.28e-3 | −2.06e-3 | −2.42e-3 | −2.51e-3 |
| 1.27e-2 | −4.17e-3 | −7.83e-3 | −1.60e-2 | −2.59e-2 | −3.07e-2 | −3.19e-2 |
| 0.161 | −0.103 | −0.116 | −0.171 | −0.292 | −0.373 | −0.396 |
| 1.23 | −1.02 | −1.04 | −1.14 | −1.47 | −2.00 | −2.31 |
| 9.42 | −9.20 | −9.30 | −9.79 | −11.6 | −15.2 | −19.5 |
| 199 | −251 | −253 | −264 | −302 | −380 | −505 |

### 2.1 Limit mu → 0 (the "near-Hamiltonian" question) — three regimes

The answer to "what is m_t as mu → 0 at fixed rho" is **not** m_t → 1.

- **(R1) Hopf-linear, rho − 1 ≪ mu:** log m_t ≈ −2π√mu·(rho − 1).
- **(R2) near-Hamiltonian window, rho − 1 = O(mu):** in slow time the equation is
  U_ττ + rho(e^U − 1) = −ε(U_τττ + U_τ), ε = √mu/(1+mu). The first-order Melnikov function has **no amplitude
  selection** (identity in §4.2), so the amplitude is fixed at higher order and the natural parameter is
  (rho − 1)/mu. Numerically, at fixed crest amplitude h, (rho − 1)/mu → F(h) and log m_t/mu^{3/2} → c(h) as mu → 0:
  F(0.5) ≈ 0.145, F(1) ≈ 0.70, F(2) ≈ 4.1, F(3) ≈ 13.8; c(0.5) ≈ −0.94, c(1) ≈ −4.6, c(2) ≈ −28, c(3) ≈ −120.
  These are the values at mu = 0.002 and are still drifting (F by 4–10 %, c by 6–30 % between mu = 0.01 and
  0.002), so they are indicative, not converged limits. F is increasing in h on the computed range.
  So in this window lambda_c·T = O(mu^{3/2}), negative.
- **(R3) relaxation limit, fixed rho > 1:** the amplitude grows like U_max ≈ log(1/mu) + const, T ≈ S*(rho)/mu,
  and **m_t → m_0(rho) ∈ (0,1)**, i.e. lambda_c·T → log m_0(rho) < 0 and lambda_c ≈ mu·log m_0/S*.
  Formal matched asymptotics (a sketch, not a proof; `limit_map.py`):
  - pulse (fast, parameter-free): x' = px, p' = −y, y' = x − y, connecting (x,y,p) = (0,0,p_in>0) to (0,0,−G(p_in));
  - slow phase (x = y = 0): p' = mu(lambda − p) with ∫p dt = 0, giving p_in' = lambda·σ(b), b = G/lambda,
    σ(b) = S − b, S = (1+b)(1 − e^{−S});
  - limit return map Π(p) = lambda·σ(G(p)/lambda); its fixed points are the solutions of Λ(p) = rho with
    Λ(p) = p/σ(B(p)), σ(B)/B = p/G(p); limit multiplier m_0 = σ'(b*)·G'(p*), σ'(b) = b e^{−S}/(1 − (1+b)e^{−S}).
  - Computed (`limit_map.csv`, p ∈ [0.02, 320]): **Λ is strictly increasing** (one limit cycle for each
    rho ∈ (1.0002, 320)), G(p)/p increasing from 1.0135, G' > 1, 0 < m_0 < 1; small amplitude
    Λ − 1 ≈ 0.405 p², log m_0 ≈ −2.106 (rho − 1)^{3/2}.
  - Agreement with the full system (O(mu) convergence):

    | rho | log m_t at mu = 0.01 | 0.002 | 0.0005 | limit log m_0 |
    |---|---|---|---|---|
    | 1.03 | −0.01863 | −0.01198 | −0.01061 | −0.01022 |
    | 1.1 | −0.06746 | −0.05656 | −0.05448 | −0.05391 |
    | 1.5 | −0.40512 | −0.38978 | −0.38688 | −0.38593 |
    | 3 | −1.70932 | −1.68471 | −1.68008 | −1.67858 |
    | 10 | −8.8261 | −8.7336 | −8.7163 | −8.7127 |
    | 50 | −58.321 | −57.811 | −57.715 | −57.683 |
    | 200 | −253.09 | −251.02 | −250.63 | −250.50 |

    p_max/lambda, −p_min/lambda and mu·T converge to p*/rho, b*, S* in the same way (`limit_map_stdout.txt`).
  - Consequence for the proof lanes: small mu is **not** a near-neutral regime at fixed rho; uniformly on the
    data |log m_t| ≳ max(2π√mu(rho−1), 2.1(rho−1)^{3/2})·O(1). The mu → 0 uniqueness question reduces (formally)
    to monotonicity of Λ, a statement about one parameter-free scattering function G. By reciprocal symmetry
    the same holds for mu → ∞.

---

## 3. Part 2 — hunt for non-uniqueness

Command: `python hunt.py 240 20261003 hunt_results.jsonl 8` (run through `research-run`, 290 s of compute).
240 random points, mu log-uniform in [1e-3, 1e3] (130 with mu < 1, 109 with mu > 1, simulated directly in the
original coordinates), rho − 1 log-uniform in [1e-3, 299]. Summary in `hunt_summary.txt`.

| test | result |
|---|---|
| primary cycle by amplitude continuation | 239/240 found; all lambda_c < 0 (max log m_t = −4.17e-4); rho(h) increasing on every branch. The one failure (mu = 141.65, rho = 195.66: stiff trough, step limit in direct coordinates) was rerun through the reciprocal rate mu = 1/141.65 with the same result (`hunt_results_supplement.jsonl`) |
| A. displacement scan of the crest return map on the attractor (W^u(E) inside; post-transient orbits outside, each seed launched from two different y0) | 239,103 pairs (106,877 on W^u(E)); seeds from 1e-4·U* inside to U* + min(max(1.5U*, 2), 4) outside: **0 pairs with the wrong sign** of U_{n+1} − U_n, i.e. the 1-D return map has exactly one crossing of the diagonal at every point |
| B. long runs from W^u(E), far-outside crests (U* + 3, U* + 6) and random 3-D initial conditions | 1673 runs: 1648 reach the primary crest; 22 too slow to arrive (near Hopf; all heading towards U*; Newton from the last point falls onto E or fails at E); 3 too expensive; **0 converge elsewhere** |
| C. Newton on the crest map from perturbed guesses | 2868 starts: 2028 → primary cycle, 681 → equilibrium, 159 no convergence, **0 new fixed points** |

The superseded first run (`superseded/`, a different set of 240 random points because the random stream is
consumed differently) gave the same conclusions: 0 new fixed points, 0 different limits, and one wrong-sign
pair that was the event artefact described in §1.

Near-tangency: the ratio r(U) = log(phi(U)/U) is **not** monotone for min(mu, 1/mu) ≲ 0.004 (24 of the 240
points; see S2), so its interior dip is the place where a second cycle would be born. `python dipscan.py`
(mu ∈ {0.03 … 0.0003} × rho − 1 ∈ {0.01 … 30}, `dipscan_stdout.txt`): where the interior dip exists it sits
near U/U* ≈ 0.2 and r there is 0.77·r(0+) at mu = 0.0027, rho = 5.02 and 0.81·r(0+) at mu = 0.0003,
rho − 1 ∈ {0.3, 1, 3}; it does not deepen as mu decreases, and r stays positive up to the cycle in all 40 cases.
**No near-tangency of the return map was found.**

---

## 4. Part 3 — structural statements on the attractor disc

Commands: `python stests.py` (62 parameter points, incl. mu > 1), `python cycle_tests.py` (300 grid cycles),
`python s3scan.py`, `python portraits.py`. Inside the cycle the disc is sampled on W^u(E) (orbits from the
linear eigenplane, transient turns discarded); outside by post-transient orbits. eta, zeta, K, W by pull-back;
the Hessian of h from one transversal variation of (U, p, y, eta, zeta).

Consistency of the geometry used (brief item (b)), at the same 62 points:
dy/dU along the section equals the pulled-back zeta to 1e-9 … 1e-12; the Hessian computed from the two
directions is symmetric to ≤ 1e-9 on well-conditioned points; K ≥ 0 and W ≥ 0 everywhere on the disc
(min K/k = 1.1e-15 next to E; W bottoms out at underflow level in deep troughs).
The section derivative zU = dζ/dU obtained from the variation agrees with finite differences across
independent W^u(E) orbits to 4–6 digits (mu = 1, rho = 3.5 and 5; mu = 0.2, rho = 5).

| statement | verdict | detail / worst case |
|---|---|---|
| **S1** one-turn mean of eta − mu decreases turn by turn along the W^u(E) spiral | **PASS 62/62** | 0 violations among all consecutive-turn pairs; the mean is also a strictly decreasing function of crest amplitude inside (62/62) and outside (31/31 points that have outside samples). It runs from eta_E − mu > 0 at E to lambda_c < 0 at the cycle. Tightest literal case: mu = 0.002, rho = 1.003, where consecutive turns differ by −2e-14 (round-off level, the spiral barely moves); there the amplitude-monotone form is the meaningful test and passes. The one-turn *integral* ∫(eta − mu)dt is **not** monotone (fails at small mu, large rho) |
| **S2** return map concave on the crest section | **REFUTED for all six coordinates** | y − k, e^U − 1, √K, √W: convex near E at **all 62** points (window starting at E and ending at U/U* between 0.003 and 0.65, median 0.07–0.23). U: concave at 49/62, convex window U/U* ≈ 0.12–0.5 at mu = 0.002 (rho − 1 ≥ 0.3), mu = 0.01 (rho − 1 ≥ 2.5), mu = 0.0027, mu = 500. Lambda_E: concave at 43/62, fails at mu ≤ 0.05 and mu ≥ 50 with large rho. The weaker property "phi_c(c)/c decreasing" (elasticity < 1) also fails: U at 9/62 (max elasticity 1.079 at mu = 0.002, rho = 9), Lambda_E at 17/62 (max 1.44), the other four at ≥ 59/62 (max 1.41–1.69). Outside the cycle the map is concave in all six coordinates (30/30) |
| **S3a** (mu − zeta(U,0))/(k(e^U − 1)) increasing along the attractor section | **REFUTED on the crest side in a band**; survives elsewhere | trough side U < 0: pass 62/62. Crest side U > 0: pass 60/62; fails at mu = 1, rho = 5 and 9, in a thin layer next to the cycle and at the cycle itself (slope of the ratio = −0.135 just inside the cycle at mu = 1, rho = 5, confirmed by finite differences at three tolerances; +0.053 at mid amplitudes). `s3scan.py`: at the cycle crest it fails for mu ∈ [0.7, 1.5] with rho ≥ 4.5–8 (mu = 0.85, 1, 1.2: from rho = 4.5; mu = 0.7, 1.5: from rho = 8) and passes for mu ≤ 0.5 and mu ≥ 2 up to rho = 10 |
| **S3L** same with f = (mu − zeta)/(1 + eta) | same pattern | fails at the same two points |
| **S3b** f/g monotone in U at fixed Liénard variable s, over the disc | **REFUTED everywhere (0/62)** | necessary condition near U = 0 is zeta ≥ mu on the line U = 0; it fails on the upstroke (p > 0): min (zeta − mu)/mu = −0.40 (mu = 0.002, rho = 1.03) |
| **S4** D²h strictly positive definite on W^u(E) | **PASS 62/62** on reliably computed points | smallest eigenvalue ratio 2.5e-12 … 0.10 (the surface is nearly ruled in troughs but never indefinite). Caveat below about the cycle itself |
| **S4b** projected cycle convex in (U, p) | **PASS 300/300** | p·p'' − p'² has no sign change on any grid cycle |
| **S5** sign of (1+eta)E = K − (zeta − mu)χ(U) | **not one-signed in general** | positive on the whole disc at 22/62 points and on 174/300 grid cycles. On the grid a negative part first appears at rho − 1 = 0.74 (mu = 0.002), 1.23 (mu = 0.005–0.05), 2.05 (mu ≥ 0.1). It is a patch (1–7 % of sampled points) on the downstroke, U < 0 and p < 0 (for example U from −5.6 to −163 at mu = 0.0027, rho = 5.02), where zeta is still above mu |
| **S6** zeta/eta ≤ 1 pointwise? | somewhere yes, everywhere no | zeta/eta ≤ 1 on part of every cycle and > 1 on another part (300/300; at E it equals gamma_E − mu > 1). x/zeta − A, zeta/eta + mu − A, k(Ue^U − e^U + 1)/K − A and x p²/W − A are **all sign-changing on every cycle**; their cycle means agree with lambda_c to 6 digits (three portrait sets, `portrait_*.npz/.png`) |

### 4.1 A caveat for anyone using D²h at the cycle
The attractor section approaches the cycle crest like a generic curve into a stable node of the return map
(multipliers m_s < m_t), so it is C^{r} there with r = log m_s / log m_t = A·T/|log m_t| − 1, and h is C² at
the cycle only if r ≥ 2, i.e. |lambda_c| ≤ A/3. This is the standard node argument, not something proved here;
the numerical symptom is clear: at mu = 1, rho = 21 (r = 1.57) the pulled-back transversal variation no longer
converges (|dy/dU − zeta|/zeta = 0.16 at the cycle) while it converges to 1e-9 at rho = 5 (r = 3.4).
On the grid r < 2 happens only for mu between about ½ and 2 at large rho (mu = 1: rho ≳ 12). The S3a failure
band coincides with small r (r ≲ 3.8): dζ/dU on the section rises steeply in a boundary layer next to the cycle.
So "h is a convex C¹ function" is consistent with everything computed; "h is C²" is not safe at the cycle.

### 4.2 The only proved item (elementary)
On every periodic orbit of U_ss = −(e^U − 1): ∮U_ss² ds = ∮U_s² ds.
Proof: ∮(e^U − 1)² ds = −∮(e^U − 1)U_ss ds = ∮e^U U_s² ds = ∮e^U U_s dU, and
∮(e^U − 1)U_s dU = ∮√(2(E − G)) dG = 0 with G = e^U − U − 1, because G takes the same values on the way up
and down. Numerical check `melnikov_check.py`: relative difference ≤ 8e-14 for h = 0.1 … 6.
Hence for U_ττ + rho(e^U − 1) = −ε(U_τττ + U_τ) the first-order energy drift per turn is
ε(rho − 1)∮U_τ² dτ: it selects no amplitude, which is why regime (R2) lives at rho − 1 = O(mu).

### 4.3 Why S2 fails at small mu (structural, not a numerical accident)
At small mu the map has two regimes inside the cycle: quasi-harmonic near E, where the expansion per turn
r(0+) = 2πα/ω ≈ π√mu(rho−1)/√rho is small and decreases with amplitude, and pulse-relaxation at larger
amplitude, where the limit map gives r ≈ (2/3)(rho − 1)b + … increasing in amplitude before it turns down to
0 at the cycle. Example (mu = 0.0027, rho = 5.02): r = 0.2871 at E, 0.2202 at U/U* = 0.25, 0.2452 at
U/U* = 0.54, 0 at U*; identical from W^u(E) orbits and from section seeds, and at rt = 1e-10, 1e-12, 1e-13.
Any coordinate that is a fixed function of crest amplitude inherits this, so concavity or sub-homogeneity of
the return map cannot be the mechanism of uniqueness for extreme rate ratios.

---

## 5. Checks of brief items

- (b) convex graph: supported (C¹, convex, gradient (zeta, −eta), K, W ≥ 0); C² fails at the cycle when r < 2 (§4.1).
- (c) four mean formulas for lambda_c: agree to 6 digits; moment identity to 4e-6.
- (f) ⟨y/(lambda − p)⟩ = mu and ⟨x/y⟩ = 1: confirmed to 8e-5 and 1.6e-8 for rho − 1 < 0.5 (trapezoid on the mesh;
  not checked to useful accuracy at larger rho, where the integrands are spikes).
- (g) ⟨U⟩ < −log rho: holds on all 300 cycles.
- (l) m_t decreasing in rho: confirmed on the full grid; near Hopf at small mu the cycles are found without difficulty
  by amplitude-fixed shooting.

## 6. Sharpest obstruction met

None against the conjecture. Against the proof strategies: (i) return-map concavity/star-shapedness is false in
every tested coordinate once min(mu, 1/mu) ≲ 0.01; (ii) Zhang-type f/g monotonicity is false over the disc for
all parameters and false even on the section for mu ≈ 1, rho ≳ 4.5; (iii) every pointwise candidate in S5/S6 is
sign-changing. What is left standing is integral or order-type: the one-turn mean of eta − mu is monotone in
amplitude (S1), h is strictly convex (S4), and rho(h), log m_t are monotone along the branch.

## 7. What I would try next

1. S1 in its amplitude form — "the one-turn mean of the planar divergence eta − mu is strictly decreasing in
   crest amplitude on the attractor disc" — implies uniqueness given (b) (every cycle would be hyperbolic
   attracting, and two nested attracting cycles need a non-attracting one between them). It passed 62/62 inside
   and 31/31 outside, including the parameter sets that break S2. It deserves a dedicated attack and a wider
   stress test (finer amplitude sampling near E, mu < 0.002, rho > 21).
2. The mu → 0 relaxation limit as a model problem: prove Λ increasing (equivalently, a monotonicity property of
   the scattering function G of x' = px, p' = −y, y' = x − y). This is parameter-free.
3. The window rho − 1 = O(mu): compute the second-order averaged drift to get F(h) and c(h) in closed form and
   check F' > 0; the first order is identically degenerate (§4.2).
4. Interval validation of lambda_c < 0 on a compact (mu, rho) box is feasible away from the three asymptotic
   corners; the corners are exactly R1–R3 above.

## 8. Limits of this evidence

Finite sampling; floating point; mu ∈ [1e-3, 1e3] and rho ≤ 300 only. Outside the cycle "the return map" is
defined only asymptotically (post-transient orbits), so outside results describe the approach to the cycle,
not an invariant object. For mu > 1 the Hessian is unreliable in deep troughs (cancellation in eta'); those
points are masked (`hessian_reliable_fraction` in `stests_summary.csv`), so S3b/S4 statements for mu > 1 cover
the crest half of the disc. Isolas of cycles far outside U* + 6, or thinner than the scan spacing
(≈ 140 seeds × 2 per point), would not be seen.
