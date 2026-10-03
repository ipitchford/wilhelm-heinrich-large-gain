# Lane `planar` — adversarial referee report (2026-10-03)

Target: the five statements the lane marks as proved (Theorem A (A1)–(A8) and Proposition 8.1 in
`attractor_theorem.tex`; the identity in `NOTES.md` §2(b)). Brief: refute them if possible.

All referee scripts are in `verify/` and were written from scratch. They share no code with the lane's
`job*.py`; `whlib` is imported only in two cross-checks (`v8`, `v9`). Python:
`/Users/admin/.venvs/analysis-py313/bin/python`. Everything is ordinary floating point (DOP853); nothing
below is an interval certificate. Batches went through `research-run` (1 CPU, 1 GB).

## 0. Outcome

| # | Statement | Verdict | Basis |
|---|---|---|---|
| 1 | Theorem A (A1)–(A5), all mu, k > 0 | **holds** | every step re-derived (§1); no numerical violation on 15 parameter pairs, 1.2e8 ordered pairs of attractor points (§2) |
| 2 | Theorem A (A6)–(A7), k > k_H | **holds** | every step re-derived; cone, stable-manifold, containment and limit-set tests pass |
| 3 | Theorem A (A8) | **holds** | every step re-derived; lambda_c from the Riccati mean, from a Riccati-free QR computation and from the monodromy agree to 1e-10 |
| 4 | Proposition 8.1 (conditional criterion) | **holds** | every step re-derived; R_1' formula confirmed to 1e-9..1e-12 by a 3D variational computation and by finite differences |
| 5 | Identity of NOTES §2(b) | **gap** | both identities are exact (sympy). The sign claim "bracket >= 0 by convexity" is stated without a range for theta and is **false for theta > 1** (bracket <= −(theta−1)K(q) < 0; negative at 100 % of test points). It holds for 0 <= theta <= 1. Corrected statement in §3 |

Theorem A and Proposition 8.1 survived every attack attempted here. The only defect found in a
"proved" statement is the missing hypothesis in claim 5; it has no effect on Theorem A, on
Proposition 8.1, or on the lane's qualitative conclusion in §3.6 (the comparison is odd under the half-turn
for every fixed theta).

Structured verdict returned to the workflow: `refuted = true`, on account of claim 5 alone (a claimed-proved
sign statement with a numerical counterexample for theta > 1). Claims 1–4 are not refuted.

Four numerical alarms appeared during the work and were each traced to the numerics, not the theorem (§2.3).

## 1. Proof audit, line by line

Notation as in the tex. "OK" means the step was re-derived independently; algebraic identities were also
checked symbolically (`verify/v3_identities.py`, 19 identities, all residuals 0).

### (A1) Lemma 2.1 and the attractor
- Coordinates: q = lambda − p gives q' = y − mu q (sympy). Forward invariance and global existence: OK.
- Upper bound: x(r) >= x(t)e^{−lambda(t−r)} needs only q > 0; the double integral c(t) is the stated one
  (a = t−s, b = s−r); c nondecreasing; for t >= 1, x' <= x(lambda − c_1 x). OK.
  Numerical test of the constant: x_max on the cycle <= lambda/c_1 in all 12 pairs, margin only a factor
  2.3–100 (e.g. (1,12): 6349 vs 15000), so the bound is not vacuous and is respected.
- Lower bound: the estimate of the convolution by s e^{−ms}, the choice of L, x(t_0) = a on a maximal
  interval inside [t_4, ∞), x' / x >= −M_q on the first L time units, monotone afterwards. OK. The times
  t_1..t_4 depend on the solution; the constants do not. OK.
- (a) absorbing set B_1, (b) invariance (injectivity of phi_s used correctly), (c) characterisation (a bounded
  entire solution in R^3 has y >= x_−, q >= x_−/mu, hence lies in a compact subset of D_chem), (d) connectedness. OK.

### (A2) Riccati
- (1.2) ⇔ (1.3) and Riccati ⇔ ratios of H' = B(t)H: OK (sympy).
- Lemma 3.1: Peano–Baker monotonicity, irreducibility of B_− (cycle 1→2→3→1), entries of Phi(t+1,t) in
  [m_0, M_0]; Hilbert-diameter bound Delta_0 = 2 log(M_0/m_0); the factorisation into n−1 contracting factors,
  one diameter-bounding factor and one nonnegative invertible factor. OK. Birkhoff–Hopf is used with its
  correct hypothesis (entrywise positive matrix, finite projective diameter).
- Lemma 3.2 (penalised maximum): OK.
- Proposition 3.3: nested compacta C(t,s), d-diameter → 0, equivariance for t' > t and for t < 0; uniqueness
  (any nonnegative entire solution has its normalisation in every C(t,s)); bounds via Lemma 3.2 in the order
  zeta-sup, eta-sup, zeta-inf, eta-inf, zeta-sup again (each uses only bounds already obtained); attraction,
  including data on the boundary of the quadrant (H_3' >= H_3 keeps H_3 > 0). OK.
- Continuity on A: R_s is continuous by finite-time continuous dependence (data (1,1) stay positive since
  Phi >= I entrywise), convergence is uniform because Delta_0, tau_0 depend only on x_∓ over A. OK.
- Remark (not an error): the rate in (3.4) is true but quantitatively empty at large rho:
  1 − tau_0 ≈ 2e^{−Delta_0/2} with Delta_0 of the order of 2 × (U-range), e.g. 1 − tau_0 ≈ 1e-31 at (1,4). The
  actual contraction happens once per turn, at the peak; through a trough (x ≈ 0) the Hilbert distance of two
  solutions stays constant (`verify/v10_stall.py`, log variables: at (0.3,12) the log-gap is 2.92e-2 from
  t = 10 to t = 67 and falls to 1e-11 across the peak at t ≈ 74). Anyone who wants to use (3.4)
  quantitatively must not expect more.

### (A3) Supporting planes
d' = J_1 d + (0,0,r), Delta' + (1+eta)Delta = −r (sympy), Delta bounded, 1+eta >= 1, limit s → −∞,
strictness by uniqueness. OK. No circularity: only (A1), (A2) are used.

### Lemma 4.2 (gauges)
O' = −(1+eta)O; K' + (1+eta)K = k(Ue^U − e^U + 1) (holds identically, O not needed);
W' + (1+eta)W − xp² = −eta'·O, so the stated equation holds on A where O = 0 (sympy). W > 0 off E. OK.

### (A4), (A5)
H is a max of affine functions with compact slope set; H(pi b) = y_b by (A3); two-sided inequality (1.6) from
(A3) applied both ways; differentiability at interior points and the Whitney C^1 condition from uniform
continuity of g. F globally Lipschitz; conjugacy; (1.7) is O = 0; div F = eta − mu on int D. OK.
"Convex C^1 graph" is to be read as stated in the theorem: h = H|_D with H convex on R², h ∈ C^1(int D),
(h, g) a Whitney C^1 jet on D. This is said explicitly (NOTES §1.3); nothing is silently weakened.

### (A6)
Characteristic polynomial, gamma_E > A, complex pair with positive real part, v_s ∈ int C_+, J_E^T n_E =
−gamma_E n_E, E^u = n_E^⊥ (sympy). Cone exclusion from Lemma 4.2(ii),(iv). No forward convergence to E:
stable-manifold theorem + cone exclusion. W^u(E) ⊂ A: backward orbit converges to E ∈ D_chem (open), forward
orbit bounded by Lemma 2.1. N_E ⊂ D by the inverse function theorem. Limit sets: generalised
Poincaré–Bendixson for the Lipschitz field F on the compact invariant set D with the single equilibrium 0; the
alternatives "{0}" and "contains a homoclinic orbit" both give a bounded entire orbit ≠ E converging to E
forwards, excluded. alpha ≠ omega for non-periodic points: monotonicity on a transversal. OK.

### Lemmas 6.1, 6.2 and (A7)
One crest/one trough, clockwise, nesting: OK (uses only v' = −k(e^u−1)/(1+eta) on {v = 0} and Jordan).
Compactness of the crest set: OK (the argument that a limit of crests is a crest uses that a non-periodic orbit
lies in the open annulus between its alpha- and omega-limit and cannot cross a separating periodic orbit).
D ⊂ closure(G): OK. closure(G) ⊂ D: the degree argument is correct — S_0 is a compact subset of D_chem
(y >= y_min > 0; p < lambda because the linear function v attains its maximum over closure(G) on J_out ⊂ D),
the boundary loop is carried into J_out by Theta_s for every s, so its winding number about q_0 stays ±1,
whereas it would be 0 if q_0 were missed; attraction of the compact set S_0 then puts q_0 in D. OK.

### (A8)
- M^T n(a_0) = m_s n(a_0); third multiplier from Liouville; m_s < m_t since eta > max(0, mu−1); (7.1). OK
  (numerically: relative error of M^T n = m_s n is 1e-12..1e-15 in 15 cases).
- (a) planar multiplier = exp∮(eta − mu); Whitney extension at Gamma_out (the multiplier depends only on
  div on J_out, so it is extension-independent). OK.
- (b) m_t > 1: stable manifold of the Poincaré map is a curve tangent to e_s-bar, but points of A approach a_0
  tangentially to n(a_0)^⊥ by (1.6), and n(a_0)·e_s ≠ 0. m_t < 1: Lyapunov stability. OK.
- (c) Gamma_in via an orbit of N_E; Gamma_out via W^u_loc(Gamma_out) ⊂ A projecting onto a neighbourhood of
  J_out, contradicting (A7). OK.
- (d) isolation of crests: one detail is left to the reader and was filled here: a nearby cycle's crest point
  is a fixed point of the *first-return* map near a_c because its orbit meets {p = 0, U > 0} only at its own
  crest (Lemma 6.1 + injectivity of pi on A), so the return time is near T. Then hyperbolic ⇒ isolated,
  compact + isolated ⇒ finite, adjacent pair, contradiction with (b). OK.
- (e), (f): OK.

### Proposition 8.1
(i) return exists, C^1; "increasing because orbits do not cross" is terse, but R_1' > 0 follows from (ii)
since g̃ > 0 on u > 0. (ii) wedge transport w' = (div F)w and the two boundary evaluations: OK (sympy +
numerics below). (iii) both criteria: OK (concave with R̃(0+) >= 0 ⇒ R̃(e)/e nonincreasing, strictly if strictly
concave; value 1 at the right end).

### Imported theorems
Birkhoff–Hopf, generalised Poincaré–Bendixson and the transversal lemma for Lipschitz planar fields,
stable/unstable manifolds for hyperbolic equilibria and hyperbolic fixed points of C^1 maps, orbital
asymptotic stability, Jordan/Schoenflies/winding number, Whitney C^1 extension: each is used within its
standard hypotheses. Theorem numbers in Hartman and Smith were not re-opened; the statements used are
the standard ones. Smith's competitive Poincaré–Bendixson theorem is indeed not used. The only external
model-specific input is the remark A = {E} for k <= k_H, which is outside the five claims.

### Circularity
None found. The order is (A1) → (A2) → (A3) → gauges → (A4) → (A5) → (A6) → (A7) → (A8); stability or
uniqueness of cycles is never assumed. (A8)(d) derives uniqueness from lambda_c < 0 on all cycles using only
(A6), (A7), (A8)(b).

## 2. Numerical attempts at falsification

### 2.1 Main battery — `verify/v1b_attractor_log.py mu rho` (library `verify/vlog.py`)

Log-variable formulation (state U, log(lambda−p), log y, log eta, log zeta), which keeps relative accuracy
through troughs as deep as U = −8356. Attractor sample = orbits of W^u(E) at several phases (or one dense
slow spiral near Hopf) plus the cycle (own Newton on the crest section). Every sign test is evaluated as
(signed sum)/(largest |term|) without underflow.

Run: `verify/run_batch2.sh` (log `verify/batch2.log`, JSON `verify/out_v1b_*.json`, table `verify/v7_summary.py`).

| mu | rho | U_min | lambda_c Riccati | lambda_c QR (no Riccati) | ordered pairs | Delta >= 0 | beyond rounding | worst Delta/maxterm (rel. sep. > 1e-6) | max abs O | min K | min W | pts in E+C_± | P8.1 variational | P8.1 FD |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.02 | 1.01 | −1.3 | −0.000181932 | −0.000181932 | 6.5e6 | 0 | 0 | −1.4e-06 | 3e-11 | >0 | >0 | 0 | 2e-10 | 7e-06 |
| 0.02 | 1.5 | −53.4 | −0.003259900 | −0.003259900 | 1.0e7 | 0 | 0 | −1.4e-06 | 3e-11 | >0 | >0 | 0 | 7e-11 | 3e-07 |
| 0.02 | 4 | −793.1 | −0.007723653 | −0.007723653 | 8.4e6 | 0 | 0 | −2.0e-07 | 3e-09 | >0 | >0 | 0 | 8e-10 | 2e-08 |
| 0.02 | 12 | −8356.5 | −0.012811292 | −0.012811292 | 7.1e6 | 0 | 0 | −4.6e-06 | 3e-10 | >0 | >0 | 0 | 4e-11 | 1e-08 |
| 0.3 | 1.01 | −0.5 | −0.001942281 | −0.001942281 | 6.5e6 | 1 | 0 | −6.7e-07 | 2e-12 | >0 | >0 | 0 | 3e-12 | 2e-07 |
| 0.3 | 1.5 | −7.5 | −0.056359795 | −0.056359795 | 6.5e6 | 0 | 0 | −1.2e-05 | 4e-12 | >0 | >0 | 0 | 4e-11 | 2e-07 |
| 0.3 | 4 | −88.6 | −0.126570340 | −0.126570340 | 1.0e7 | 1 | 0 | −1.2e-05 | 1e-11 | >0 | >0 | 0 | 2e-10 | 1e-06 |
| 0.3 | 12 | −926.1 | −0.203960370 | −0.203960370 | 8.6e6 | 0 | 0 | −3.3e-06 | 2e-11 | >0 | >0 | 0 | 6e-11 | 1e-04 |
| 1 | 1.01 | −0.4 | −0.003983998 | −0.003983998 | 6.4e6 | 0 | 0 | −6.8e-07 | 4e-13 | >0 | >0 | 0 | 7e-12 | 3e-07 |
| 1 | 1.1 | −1.7 | −0.038211636 | −0.038211636 | 6.4e6 | 0 | 0 | −1.3e-06 | 6e-13 | >0 | >0 | 0 | 3e-11 | 6e-08 |
| 1 | 1.5 | −6.0 | −0.146389341 | −0.146389341 | 6.5e6 | 0 | 0 | −3.3e-06 | 4e-12 | >0 | >0 | 0 | 8e-11 | 2e-07 |
| 1 | 4 | −64.2 | −0.395940343 | −0.395940343 | 1.1e7 | 0 | 0 | −6.0e-07 | 4e-12 | >0 | >0 | 0 | 1e-09 | 2e-07 |
| 1 | 12 | −668.3 | −0.674126202 | −0.674126202 | 9.1e6 | 0 | 0 | −4.7e-07 | 9e-12 | >0 | >0 | 0 | 2e-11 | 2e-07 |
| 3 | 1.5 | −7.2 | −0.187751996 | −0.187751996 | 6.5e6 | 0 | 0 | −3.2e-07 | 2e-12 | >0 | >0 | 0 | 1e-10 | 2e-06 |
| 3 | 12 | −878.5 | −0.683406959 | −0.683406959 | 8.9e6 | see §2.3(c) | see §2.3(c) | see §2.3(c) | 1e-11 | >0 | see §2.3(c) | 0 | 8e-11 | 1e-03 |

Columns:
- **(A3)/(A4)**: Delta_ab = zeta_a(U_b−U_a) − eta_a(p_b−p_a) − (y_b−y_a) < 0 for all ordered pairs of sampled
  attractor points (both orders, hence also the two-sided inequality (1.6) and grad h = (zeta, −eta)).
  "Delta >= 0" counts pairs with a nonnegative value; the two such pairs are near-coincident points
  (relative separation < 1e-6) whose value is below the rounding amplification of the difference (§2.3(b)).
- **(A5)**: O = n·f = 0 to 1e-9..1e-13 relative; K > 0, W > 0 on every sample.
- **(A6)**: no sample in E + C_+ or E + C_−; Lambda_E > 0 on every sample.
- **(A2)**: log eta, log zeta finite on every sample (so eta, zeta > 0) and inside the bounds (1.4) in all 15
  cases; periodic Riccati pair reached from data (1e-8,1e-8), (1e8,1e8), (1e-8,1e8), (1e8,1e-8) in 2–5 periods
  with log-gap <= 5e-12. Zero data and the rate (3.4) were tested in original variables (`verify/v1_attractor.py`)
  where that formulation is valid: from (0,0), (1e6,1e6), (0,1e6), (1e6,0) the gaps at t = 1, 2, 5, 10, 20, 40
  respect Delta_0 tau_0^{floor(t)−1} at (1,1.01), (0.3,1.01), (0.02,1.01), (1,1.5), (0.3,1.5), (3,1.5)
  (e.g. (1,1.01): Delta_0 = 3.435, tau_0 = 0.696; gaps 0.38, 0.055, 1.9e-4, 1e-8).
- **(A8)(a)**: lambda_c = <eta> − mu (Riccati) against a computation that never touches the Riccati equation
  (discrete QR of the 3D variational flow, sum of the two leading exponents): agreement to <= 2e-11 in all
  15 cases; log of the leading eigenvalue of the crest-section derivative agrees wherever m_t is above rounding.
  Left-eigenvector relation M^T n(a_0) = e^{−(1+<eta>)T} n(a_0): relative error <= 2e-12.
- **Prop 8.1(ii)** ("P8.1 variational"): with M the derivative of the 3D crest-return map at a crest point of
  W^u(E), the vector M(1, zeta(u)) is parallel to (1, zeta(R_1 u)) (defect <= 4e-10: the section curve of A is
  invariant with tangent (1, h_u) = (1, zeta)) and its U-component equals g̃(u)e^{iota(u)}/g̃(R_1 u) with the
  stated relative error (8 crests per case, from the first turns to the cycle). "P8.1 FD": the same derivative
  from finite differences between neighbouring orbits of W^u(E) (no variational equation, no zeta).
  Remark 5.1 (g̃ = k(e^u−1)/(1+eta) on the crest axis): relative error <= 2e-10. Crests of successive turns
  increase monotonically to the cycle crest in all cases.

### 2.2 Other tests

| script | attack | result |
|---|---|---|
| `v3_identities.py` | sympy check of 19 identities: (1.3), O, K, W equations, Thm 4.1, H' = BH ⇔ Riccati, eigen-relations at E, Remark 5.1, Remark 8.2, both identities of claim 5 (and the general formula D omega − C omega = conj(f ∧ (Dnu f − J nu))), wedge transport of Prop 8.1 | all residuals 0 (`out_v3_identities.txt`) |
| `v2_closepairs.py` | consecutive crests of one orbit, section hit exactly (Hénon), three tolerances; curvature estimates c_fwd, c_bwd = −2 Delta/d² both ways | (1,1.01), (1,1.1), (0.3,1.5), (3,1.3): c_fwd ≈ c_bwd > 0 and tolerance-independent wherever 1e-12/d² is small (e.g. 1.7113 at (1,1.1)); sign anomalies only at separations d below 1e-6 (rtol 1e-10), 1.1e-7 (rtol 1e-12), 4e-8 (rtol 3e-14), with abs(Delta) at most 7e-12, 8e-14, 1.2e-14 respectively |
| `v5_mu_gt1.py` | mu > 1 troughs in cancellation-free variables | (3,12), (3,1.5), (10,4): no pair with relative separation > 1e-6 has Delta >= 0; worst −8e-7; W > 0 |
| `v4_bracket.py` | sign of the bracket of claim 5, theta from 0.05 to 3 | §3 |
| `v6_extras.py` (12 pairs) | (i) Lemma 2.1 upper bounds; (ii) W^s(E) backward; (iii) backward Riccati from perturbed positive data; (iv) W^u(E) inside the projected cycle; (v) random initial data in D_chem | (i) respected 12/12; (ii) both branches leave D_chem within 15 time units, 24/24 (so no bounded entire orbit on W^s(E)); (iii) 144/144 perturbed solutions leave the closed quadrant within half a period; (iv) all inside after the chord artefact is removed (§2.3(d)); (v) all orbits end on the same cycle (crest deviation <= 6e-11; at (0.02,1.01) the gap is 1.2e-10 after 2.5e5 time units and decays at the rate lambda_c, `v6d_slow_convergence.py`) |
| `v8_crosscheck.py` | referee integrators vs the v2.0 bundle receipts and whlib | T to 6 decimals, m_t to 7–9 digits in 4 cases |

### 2.3 Alarms raised and resolved (none is a counterexample)

(a) *Original-variable integration at rho >= 4* (`v1_attractor.py`, `batch1.log`): thousands of apparent
violations of (A3), W < 0, eta < 0 at (0.3,12), (0.02,1.5), (0.02,4), (1,12). Cause: in the trough eta and zeta
fall below any absolute tolerance (true minima e^{−17}..e^{−800}); the integrated values go slightly
negative and everything built on them is noise. In log variables all of these disappear (table above).
`v1` is kept only for the zero-data pull-back tests at moderate parameters.

(b) *Near-coincident pairs*: +5e-14 at (1,1.01) in `v1`, one pair each at (0.3,1.01) and (0.3,4) in `v1b`
(Delta/maxterm 2e-7 and 8e-7, absolute 1e-19 in the second), and the lane's own +3e-10 at (1,1.1).
For the (0.3,4) pair the two points are 1.7e-7 apart at U = −62 (a cycle point and a spiral point that has
converged onto the cycle); the symmetric part of the two directions is −1.35e-21 each, as predicted by
−W s²/2, and the antisymmetric part is the 1e-5 relative error of differences of numbers agreeing to 9
digits. `v2_closepairs.py` shows the pattern directly: the anomaly threshold moves down with the tolerance and
the resolved curvature is strictly positive. The lane's reading of its +3e-10 pair as rounding is supported.

(c) *mu > 1, deep trough* ((3,12) in `v1b`: 96849 pairs with Delta/maxterm in [0, 5e-10], one flagged,
min W = −6e-15, eta − (mu−1) = −4e-14). Cause: for mu > 1 the trough has eta → mu−1 and y → (mu−1)(lambda−p),
so −eta Δp − Δy cancels to leading order and the sign is below double precision in (y, eta). With
r = y − (mu−1)(lambda−p) (r' = x − mu r, r > 0) and eps = eta − (mu−1) (eps' = zeta − eta eps, eps > 0) one has
exactly Delta_ab = zeta_a ΔU + eps_a Δq − Δr; in these variables (`v5_mu_gt1.py`) all violations vanish.
mu > 1 is outside the requested parameter set; it was added as an extra and is clean.

(d) *Polygon containment* (`v6` (iv)): up to 38 of 6000 spiral points flagged outside the 4001-vertex polygon
of the projected cycle. With 160001 vertices all are inside except one point at (0.02,12); an exact
comparison (cycle integrated backward from its crest to the same U, `v6c_flagged_point.py`) puts that point
3.0e-4 *inside* the cycle. Chord artefact (`out_v6b.txt`, `out_v6c.txt`).

## 3. Claim 5 — what is true

Identity (exact, any function H): F(theta q) ∧ F(q) = −v·B(theta, q), B = theta(H(q) − k) + k − H(theta q).

Sign of B. With H convex and H(0) = k:
- 0 <= theta <= 1: B >= 0. If moreover theta ∈ (0,1) and theta q ∈ D \ {0}, the supporting plane at the point of A
  above theta q gives B > (1 − theta) K(theta q) > 0.
- every theta >= 0: the supporting plane at q gives H(theta q) >= theta h(q) + (1−theta)k − (1−theta)K(q), i.e.
  B <= (1 − theta) K(q). Hence **for theta > 1 and q ≠ 0: B <= −(theta − 1) K(q) < 0.**

So "the bracket is >= 0 by convexity" is correct only for theta <= 1, a restriction the notes do not state.
`v4_bracket.py` at (1,1.5), (0.3,1.5), (1,4), (0.02,1.5), (1,1.01), 3000 points each: B > 0 at 100 % of points for
theta ∈ {0.05, 0.3, 0.6, 0.9, 0.99} (lower bound of B from Delaunay interpolation >= −3e-9, the mesh error;
−1e-4 on the coarser mesh at mu = 0.02), B < 0 at 100 % of points for theta ∈ {1.01, 1.1, 1.5, 3}, and
B <= (1−theta)K(q) to rounding throughout.

Correct geometric statement: where p > 0 the flow crosses the homothetic copies theta·Gamma towards theta = 1
from both sides (outward through the shrunk copies, inward through the enlarged ones); where p < 0 it
crosses them away from theta = 1. For the Massera range theta < 1 this is the lane's sentence. If h is used
instead of H (as in the notes), theta q ∈ D is needed as well, and star-shapedness of D about 0 is only
numerical. The 3D identity D omega − C omega = k(Ue^U − e^U + 1)(0, p, p') for nu = (U, p, y−k) is exact.

## 4. Side findings (outside the five claims)

- `whlib.cycle_profile` loses positivity of the Riccati pair in deep troughs (`v9_whlib_trough.py`):
  min eta = −3.6e-9, min W = −3.4e-9 at (0.02,4); min eta = −1.6e-21 at (0.3,12). lambda_c from whlib stays
  correct to 9 digits (its monodromy value is off by 4e-5 at (0.3,12), where m_t is near rounding). Pointwise sign tests of eta, zeta,
  K, W built on whlib are unreliable for rho >= 4, mu < 1; a log-variable integrator is in `verify/vlog.py`.
- Crest-section curvature h_uu is strictly positive in every resolved case (`v2`), consistent with the lane's
  statement that h is strictly convex inside (numerical only).
- Not checked: the lane's Job 2 numerics (job2_*.py), which are reported as numerical and not as proved.

## 5. Files in `verify/`

`vlog.py` (log-variable library), `vlib.py` + `v1_attractor.py` (original variables; valid near Hopf only),
`v1b_attractor_log.py`, `v2_closepairs.py`, `v3_identities.py`, `v4_bracket.py`, `v5_mu_gt1.py`, `v6_extras.py`,
`v6b_polygon.py`, `v6c_flagged_point.py`, `v6d_slow_convergence.py`, `v7_summary.py`, `v8_crosscheck.py`,
`v9_whlib_trough.py`, `v10_stall.py`, `run_batch{1,2,3}.sh`, logs `batch{1,2,3}.log`, outputs `out_*.json`, `out_*.txt`.
