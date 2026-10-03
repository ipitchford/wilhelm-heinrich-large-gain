# VERIFY — lane `monotone` (adversarial referee pass, 2026-10-03)

Brief: refute the five statements the lane marks "proved". Materials read in full: `work/BRIEF.md`,
`NOTES.md`, `ds_theorem.tex`, the lane's output files, and the imported v2.0 statements with their
proofs (`manuscript.tex` Prop. 2.1, Prop. 2.2, Thm 3.1, Cor. 3.2, Thm 4.1, Thm 4.2, Prop. 4.4).
All referee scripts are in `verify/`; none imports the lane's scripts or `mono_lib.py`.

## 0. Verdicts

| # | Statement | Verdict | One-line reason |
|---|-----------|---------|-----------------|
| 1 | Dual supporting planes (Thm 6) | **holds** | Every step checked; 5.77 M chord determinants on 20 cycles, 0 violations |
| 2 | Twelve-event order (Thm 4) | **holds** | Three crossing derivatives re-derived exactly; order found on 20/20 cycles |
| 3 | Rank-one window identity (Prop 8) | **holds** | Identity re-derived; reproduced to ≤ 6e-6 on 20 cycles for well-conditioned test functions |
| 4 | Strong-stable ratios (Lemma 1) | **holds**, with a wording defect | Lemma correct. "Equivalent to the conjecture" is an overstatement (see 3.2) |
| 5 | Crest return map on R (NOTES R7) | **gap, repaired here** | Written proof never shows P is defined and C^1 on all of R; closed in section 4 |

No statement is false. No numerical counterexample was found. One gap (statement 5) was found and
repaired; one sentence (statement 4) claims more than is proved.

## 1. Line-by-line audit of `ds_theorem.tex`

Imported facts. (F1) = v2.0 Prop. 2.1 (stated exactly as used; it allows m_t = 1). (F2) = Thm 4.1.
(F3) = Thm 4.2 together with the zero counts inside the proof of Thm 4.1 (p has exactly two simple
zeros; p' and y' each have exactly one sign change per half-cycle), which is what the sign table in
(F3) needs. (F4) = Cor. 3.2 eq. (3.6). None of these uses stability or uniqueness, so there is no
circularity. The nestedness of projections of distinct cycles used in Thm 6 is the last sentence of
Thm 4.1; the tex cites it as "(F2)" although its own (F2) item does not list it (citation
imprecision only).

**Lemma 1.** Sign change (a,b,c) = (v1,−v2,−v3) and τ = −t give a_τ = b, b_τ = μb + c,
c_τ = xa + c: cooperative, irreducible (cycle a←b←c←a, x > 0), so the backward period map is a
positive matrix. Its spectrum is {1, 1/m_t, 1/m_s}; by (F1) m_s < min(1, m_t), so the Perron root
is 1/m_s, simple, positive eigenvector. Riccati equations and "r > μ at the minimum of r" checked.
For r > 1: c̃ = c + (μ−1)b gives a' = −b, b' = −b − c̃, c̃' = −xa − μc̃ (checked by hand and by
sympy). This system is conjugate to the first by a constant matrix, so it has the same
multipliers; its positive Perron solution maps back to a Floquet solution with the simple
multiplier m_s, hence to a positive multiple of e_s. So c̃ > 0 and r' = r² − r − c̃/a gives r > 1.
No gap.

**Lemma 2.** (Mu)×v + u×(Mv) = (tr M − Mᵀ)(u×v) verified symbolically for a general 3×3 M;
tr J = −A; quasi-periodicity from m_s e^{AT} = 1/m_t; component formula, chord identity and
(Ψ·z'')' = Ψ₃xp² verified exactly. No gap.

**Theorem 4.** Checked interval by interval against (F3):
- ψ_p = p' + rp > 0 on [t_tr, t_Ymin] and < 0 on [t_cr, t_Ymax] (both terms of one sign, not both zero).
- On C and F: (θ − r)' = (θ − r)(θ + r − μ) + y'/p + s, so every zero is a strict up-crossing
  (y'/p > 0 on C ∪ F, s > 0); θ runs 0 → +∞; exactly one zero.
- g₂ = y'/p + s: exact identity g₂' = sψ_p/p − g₂(1 + p'/p). On the two intervals ψ_p/p > 0 by the
  first bullet (no circularity: the Ψ₃ analysis does not use Ψ₂).
- g₁ = y'/p' − s/r: exact identity g₁' = xψ_p/(rp') + g₁(μ − 1 + y'/p' + s/r); same argument.
- The four closed intervals plus two open intervals cover the period for each component, so
  "exactly two zeros" is justified, and the stated interlacing follows. No gap.

**Lemma 5.** ℓ' = 0 iff (p,p')(t) ∥ (p,p')(t₁), which forces pp' < 0, i.e. t ∈ C ∪ F (C and F are
exactly the set pp' < 0), and θ(t) = r(t₁) > μ; each level crossing has θ' > 0. ℓ''(t₁) =
p²r(r−μ) + py' > 0. Two critical points, one a strict local minimum, so it is the unique global
minimum. Convex-hull argument for interior points is correct. No gap.

**Theorem 6.** B' = Ψ₃w, w ≥ 0, and the sign pattern of Ψ₃ make t₁ + jT the maximiser of B on each
window; Ψ(t₁) = κ(p', −p, 0) with κ = −Ψ₂/p < 0; B(t₁ + jT) = κ_j L(·) < 0 by Lemma 5. Applies
verbatim to z(· + s), to E, to z'' (with ℓ''(t₁)), and to any entire solution with L > 0 at the
times t₁ + jT. For an inner cycle, every (point, t) pair is reached by choosing the phase of the
inner solution. No gap.

**Proposition 8.** ω_φ' − Jω_φ = −L₃[φ]e₃ verified; boundary term evaluated with Ψ₃(t₁) = 0,
p'(t₁) = −rp(t₁), Ψ₂(t₁) < 0. Three special cases: (φ' + rφ)(t₁) equals r, (r−μ)p' − y', and
−(p + rU) respectively, all > 0 at t₁ (the last needs U(t₁) < 0, which follows from (F4) as
written). No gap.

## 2. Numerical falsification attempts

Own integrator: log coordinates (U, ln y, ln(λ − p)), DOP853, rtol 1e-12; own cycle finder
(slow-manifold bisection then Newton with finite-difference Jacobian). All sign decisions use
signed-log arithmetic (each quantity is a sum of signed monomials with stored log-magnitudes), so
they remain valid when y ≈ e^{−900} and U ≈ −8356.

Parameter set (20 cycles): μ ∈ {1, 0.3, 0.02} × ρ ∈ {1.01, 1.5, 4, 12} (the requested grid), plus
(3, 1.01), (3, 1.5), (3, 4), (3, 12), (50, 1.5), (50, 4), (0.1, 1.2), (7, 2.5). The lane tested
ρ ≤ 5 and μ ≥ 0.05 only. Extremes reached: T = 877.7, U_min = −8356 (μ = 0.02, ρ = 12).
λ_c ∈ [−0.683, −0.00018]: all 20 cycles attracting.

### Statement 1 (`v1_floquet.py`, `v1b_inner.py`)
- Chords: 5,770,615 (base, chord) pairs, non-uniform grids resolving both crest and trough.
  Violations: 0. Every pair with |normalised det| < 1e-7 (303 pairs) lies within |Δt|/T ≤ 1.2e-3
  of the diagonal, where det → 0 quadratically; all of those are still negative.
- E: worst normalised det −4.5e-3. z'': worst −1.7e-5 (μ = 0.02, ρ = 12). All negative.
- Signed-log determinant cross-checked against `numpy.linalg.det` on ~7,600 random pairs: 0 sign mismatches.
- Inner entire solutions (no second cycle exists at any tested parameter, so the test bodies are
  the W^u(E) spiral and interior orbits after two periods): 16 parameter pairs, ~13,000 inner
  points × 400 base phases each. 0 violations; the worst ratio scales with the crest-amplitude
  deficit δ (−9e-8 … −3e-6 for δ ∈ [1e-7, 1e-5)) and is ≈ 0⁻ only for points numerically on the cycle.
- Lemma 5 inputs: t₁ ∈ F, U(t₁) < 0, L(E) > 0, ℓ''(t₁) > 0, ℓ > 0 away from t₁ on 20/20.
- Monodromy cross-check (`v3_monodromy.py`, no Riccati inside Φ): Φᵀq = m_t q with relative
  residual 2.5e-12 … 8e-6 on 14 cycles; the eigenvector of Φ for the smallest multiplier equals
  (1, −r, −s) to ≤ 6e-12 on the 4 cycles where m_s is above the roundoff of Φ.

### Statement 2 (`v1_floquet.py`)
- Exactly two sign changes per period of each of p, p', y', Ψ₁, Ψ₂, Ψ₃ on 20/20 cycles, in the
  claimed cyclic order on 20/20. Ψ > 0 at the crest and Ψ < 0 at the trough on 20/20.
- At four deep-trough cycles the trough and the up-crossing of Ψ₂ fall in the same double-precision
  time cell (the gap is ~y(trough)/p', below 1e-300). The order there is decided exactly: at p = 0,
  Ψ₂ = ϖy', and y' < 0 at the trough in signed-log arithmetic. No other ties.

### Statement 3 (`v1_floquet.py`)
- Window integral by adaptive ODE quadrature from t₁ to t₁ + T against the right-hand side, with
  λ_c taken from the independent forward (η, ζ) system. Test functions 1, p', −U, cos(ωt),
  sin(3ωt): worst relative error 8e-9, 6e-6, 7e-9, 2e-8, 6e-6 (all at μ = 0.02, ρ = 12; ≤ 2e-7
  on the other 19 cycles). Signs of the three special integrals as claimed on 20/20.
- φ = y agrees to ≤ 4e-7 on the 15 cycles where it is well conditioned. On five deep-trough
  cycles with μ < 1 (μ = 0.3: ρ = 4, 12; μ = 0.02: ρ = 1.5, 4, 12) the right-hand side carries the
  factor (x + (r−1)y)(t₁) and is 1e-8 … 1e-70 times the size of the two integrals whose difference
  it is, so the comparison is pure roundoff. This is conditioning of the test, not a counterexample.
- λ_c from three sources (⟨η⟩ − μ, ⟨r⟩ − A, ln ϖ over the window) agree to 9 digits on 20/20;
  the finite-difference return map agrees wherever it resolves m_t.

### Statement 4 (`v1_floquet.py`, `v4_sym.py`)
- ⟨r⟩ = 1 + ⟨η⟩ to 9 digits on 20/20; s > 0.
- r − max(1, μ): plain double precision resolves it only down to ~1e-10 (it then returns values
  like −2e-11, pure roundoff). The margin was therefore integrated in log variables in the
  conjugate (a, b, c̃) system (μ ≤ 1) or the (a, b, c) system (μ ≥ 1); minimum margins range from
  10^{−0.06} to 10^{−369}. Where both are resolvable the two systems give the same r (relative
  difference of the margin 1e-11 … 3e-4, consistent with roundoff). This confirms the one
  non-trivial step of the lemma (the conjugate system has the same dominant solution). Positivity
  itself is structural in log variables and so cannot be refuted by this computation; the proof is
  what carries it.
- Scalar equation r'' + (A − 3r)r' + r(r−1)(r−μ) = x and its mirror for g = 1 + η: exact (sympy).

### Statement 5 (`v2_retmap.py`, `v5_mp_spot.py`)
- 468 points of R on 14 parameter pairs: U₀/U_crest ∈ {0.02 … 2.0}, (y−k)/(x−k) ∈ {1e-4 … 1},
  including the boundary x = y. A first return exists in every case (return times up to 1.9e5);
  every image has x₁ > y₁ > k.
- DP by the full variational equation plus projection (independent of the compound-system
  argument): 1,373 entries above the roundoff floor, all positive, agreeing with the
  compound-system values to ≤ 2.2e-8. The other 499 entries are below the floor in double precision
  (down to 1e-1670).
- High precision (own Taylor integrator in mpmath, 80–160 digits, variational equation plus
  projection): 4 points, including one with entries 1e-64 … 1e-69 (μ = 1, ρ = 12) and one on
  x = y. All 16 entries positive and equal to the compound-system values to 11 digits.
- Cycle crest: Perron slope of DP equals ζ(crest) to ≤ 1e-12 and ln(m_t)/T equals λ_c to 9 digits
  on 14/14.

### Exact identities (`v4_sym.py`)
31/31 pass, written from the equations of motion only.

## 3. Defects

3.1 **Statement 5, gap.** NOTES R7 proves: *if* the orbit from a point of R returns to
{p = 0, y > k}, then the image is in int R and DP > 0. It does not show that the return exists for
every point of R, nor that the first return is transversal (needed for P to be C^1 there). v2.0
Prop. 4.4 only treats the local return at a cycle crest, so it cannot be cited. Repaired in
section 4.

3.2 **Statement 4, wording.** What is proved is ⟨r⟩ = 1 + ⟨η⟩, hence ⟨r⟩ < A ⇔ λ_c < 0. "The
conjecture is equivalent to ⟨r⟩ < A" is correct only in the brief's shorthand (obligation
"⟨η⟩ < μ on every cycle"). As a statement about uniqueness it is one-directional: ⟨r⟩ < A on every
cycle ⇒ uniqueness (brief 3(a)); uniqueness ⇒ ⟨r⟩ ≤ A only, since a unique cycle with m_t = 1 is
not excluded. Lemma 1 itself is unaffected.

3.3 **Citation.** Nestedness is cited as (F2) but is not in the tex's (F2) item; it is in v2.0 Thm 4.1.

3.4 **Coverage of the lane's numerics.** The lane's 23 cycles stop at ρ = 5 and μ = 0.05, and its
`ds_independent.py` reports up to 8/25 return-map points as unresolved. The referee runs close
both (ρ = 12, μ = 0.02 and 50; high-precision DP).

## 4. Repair of statement 5

Claim: for every z₀ = (U₀, 0, y₀) with x₀ ≥ y₀ > k, the forward orbit has a first later zero of p
with p' > 0 and a second with p' < 0, both transversal; hence P is defined and C^1 near z₀.

(i) *No zero of p on W^s(E) \ {E}.* E is hyperbolic with stable eigenvector
v_s = (1, −γ_E, −γ_E(γ_E − μ)), i.e. (1, γ_E, γ_E(γ_E − μ)) > 0 in the coordinates (a, b, c). On
W^s(E), z'(t) = −γ_E c₀ e^{−γ_E t}(v_s + o(1)) with c₀ ≠ 0, so z'(t) lies in the interior of ±K
(K the positive cone of (a, b, c)) for large t. z' solves the variational equation, which is
cooperative in (a, b, c) in backward time along any orbit; hence z'(t) ∈ ± int K for all t, and in
particular p ≠ 0.

(ii) *No tangency.* For the variation (1, 0, 0) the wedge w(0) = (y₀ − k, x₀ − y₀, 0) is
nonnegative and nonzero, so w(t) > 0 for t > 0 (cooperative irreducible compound system). At a
point with p = p' = 0 the first wedge component v₂p − v₁p' vanishes. So p and p' never vanish
together for t > 0.

(iii) *Existence.* p'(0) = k − y₀ < 0, so p < 0 just after 0. If p < 0 for all t > 0, U decreases.
If U → −∞ then x → 0, y → 0 and p' ≥ k − y > k/2 eventually, contradicting p < 0. If U stays
bounded the orbit is bounded with U monotone, its ω-limit set is {E}, the orbit lies on W^s(E),
contradicting (i) because p(0) = 0. So there is a first zero t_a, with p'(t_a) > 0 by (ii). The
same dichotomy on (t_a, ∞) (if U → +∞ then y → ∞ and p' ≤ k − y → −∞) gives a next zero t_b with
p'(t_b) < 0, i.e. y(t_b) > k.

(iv) Both zeros are transversal and p(t) = p'(0)t + O(t²) near 0, so the zero pattern persists for
nearby section points and P is C^1 at z₀. The lane's wedge computation then gives DP > 0 and
x₁ > y₁ > k. The Perron-eigenvector statement follows from n = (ζ, −η, −1) annihilating z' and e_t
(v2.0 Prop. 2.2) and m_t > m_s.

With (i)–(iv) statement 5 is fully proved; Smith's convergence theorem is not needed.

## 5. Reproduction

    PY=/Users/admin/.venvs/analysis-py313/bin/python ; cd work/lanes/monotone/verify
    $PY v1_floquet.py      > v1_grid.out       # statements 1-4, 20 cycles, ~50 s
    $PY v1b_inner.py       > v1b_inner.out     # statement 1, inner entire solutions, ~60 s
    $PY v2_retmap.py       > v2_retmap.out     # statement 5, 468 points of R, ~3.5 min
    $PY v3_monodromy.py    > v3_monodromy.out  # e_s and Psi against the full monodromy, ~20 s
    $PY v4_sym.py          > v4_sym.out        # 31 exact identities
    $PY v5_mp_spot.py mu rho f theta digits    # high-precision DP at one point of R (v5_mp_spot.out)

Single process, each under 5 minutes and well under 1 GB.

## 6. Limits of this pass

- Floating point throughout, apart from the sympy identities and the mpmath spot checks; nothing
  is an interval certificate.
- Statements about "every cycle" and "any cycle nested inside" were tested on the one cycle that
  exists at each parameter and on interior entire solutions; a second cycle was never available.
- Three referee-side artefacts were found and fixed during the pass and are not findings about the
  lane: a first inner-solution filter that mistook on-cycle points for interior points; the
  ill-conditioned φ = y window test; and an absolute-error step control in the Taylor integrator
  that lost x in deep troughs (first run produced a false return; fixed by relative control on x,
  after which x/(ke^U) − 1 ≤ 1e-71 at the event).
- The imported v2.0 facts were read and their use checked, but v2.0 itself was not re-refereed
  beyond the proofs of Prop. 2.1, Thm 3.1, Cor. 3.2, Thm 4.1, Thm 4.2 and Prop. 4.4.
