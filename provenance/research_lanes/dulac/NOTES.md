# Lane `dulac` — notes (2026-10-03)

Target (brief items 3(b),(c),(d)): find Λ on the attractor disc D∖{E} with `R := 𝔻Λ + η − μ` of one sign, then prove it.

**Outcome in one paragraph.** No provable Dulac function was found and the conjecture is NOT proved in this lane.
Numerically a Dulac function exists for each of the 10 sampled parameter sets (26-feature family, sign R ≤ 0), but
no small family works for all sets, no closed-form coefficients emerged, and every candidate found reduces to an
*upper* bound on the tangent-plane tilt `K/V` on the descending half `p<0` of a cycle, which none of the available
convexity inequalities supplies (they give lower bounds there). On the way, a rigorous intermediate theorem was
obtained (Theorem D below): every cycle has a strictly convex, strictly star-shaped (U,p)-projection, two new
exact formulas for the weak exponent, an explicit state function η₀ with cycle mean exactly μ that brackets η
with a sign fixed by the sign of p, and three "monotone-integral" inequalities. All statements of Theorem D were
tested on 20 fresh cycles before the proofs were written.

Notation as in the brief: `H = y−k`, `x = k e^U`, `p' = k − μp − y`, `A = 1+μ`, `ψ(U) = Ue^U − e^U + 1 ≥ 0`,
`𝔻` = derivative along the flow, `⟨·⟩` = cycle mean. Sign convention used throughout: **R ≤ 0 is the only possible
sign** (see Prop. B1).

---------------------------------------------------------------------------------------------------------------

## 1. Rigorous results (status: PROVED, every step written; inputs are only results proved in v2.0)

Inputs from v2.0 used: tangent equation (2.8) `p'' = (η−μ)p' − ζp`; orthogonality (2.7); `K>0` (3.5); `W>0` (3.8);
Theorem 3.1 for two phases of one cycle; Theorem 4.1 (regular injective projection, p has exactly two simple zeros —
trough with U<0, p'>0 and crest with U>0, p'<0 — and exactly one max and one min) and (4.3).

Define on a cycle

    V  := p² − U p'            = p² + U(μp + y − k)                 (state function)
    Q  := p'² − p p''          = p'(k − y) + p(x − y)               (state function)
    V₃ := U(x − y) − p(y − k)                                        (state function)
    η₀ := V₃ / V ,   η₁ := (xV + kψ(U)p')/Q − 1 ,   κ := K/V ,   ϖ := W/Q .

### Theorem D
For every positive nonconstant cycle of (1.5), at every (μ,k) ∈ 𝒫 (indeed for every (μ,k) at which a cycle exists):

**(a) Positivity.** `V > 0` and `Q > 0` at every phase.
Consequently the projection t ↦ (U,p) is a *strictly convex* closed curve (signed curvature `−Q/(p²+p'²)^{3/2} < 0`),
and it is *strictly star-shaped about the equilibrium* (polar angle has derivative `−V/(U²+p²) < 0`).

**(b) Wronskian identities.**
`V' = (η−μ)V + pK = −μV + V₃`, `Q' = (η−μ)Q + pW = xV − AQ + kψ(U)p'`, and `V₃ = pK + ηV`.

**(c) Exponent formulas and an explicit neutral comparison function.**

    λ_c = −⟨pK/V⟩ = −⟨pW/Q⟩ ,        λ_c T = −∮ (K/V) dU = −∮ (W/Q) dU ,
    ⟨η₀⟩ = ⟨η₁⟩ = μ ,     η − η₀ = −pK/V ,     η − η₁ = −pW/Q .

So `η < min(η₀,η₁)` on the ascending half `p>0`, `η > max(η₀,η₁)` on the descending half `p<0`, and
`λ_c = ⟨η − η₀⟩`. Moreover the whole gradient is an affine function of the single unknown κ>0:
`η = η₀ − pκ`, `ζ = ζ₀ − p'κ` with `ζ₀ := Q/V`, and κ solves the scalar equation
`κ' = kψ(U)/V − (1 + η + η₀ − μ)κ`.

**(d) Monotone integrals.** For every (α,β) ≠ (0,0):  `⟨(αp + βp')(αη + βζ)⟩ < 0`. In particular
`⟨ηp⟩ < 0`, `⟨ζp'⟩ < 0`, and `⟨ηp⟩⟨ζp'⟩ > ⟨ζp⟩²` (note `⟨ηp'⟩ = ⟨ζp⟩`).
At equal U the lower-branch point has larger η and larger y; at equal p the crest-side point has larger U, y and ζ.

**(e) Compound form.** `ω := (V, V₃, Q)` is the sign-conjugated wedge `r ∧ z'` (r = (U,p,H)) and satisfies
`ω' = Cω + kψ(U)·(0, p, p')` with C the compound matrix of v2.0 (4.6). Equivalently `(η₀,ζ₀)` solves the Riccati
system (2.5) with defects `kψ(U)(p, p')/V` in the two equations. (For the supersolution lane: the defect is ≥0
exactly on the quadrant p,p' ≥ 0. Note V₃, hence η₀, is negative on part of the cycle for μ ≤ 0.3 (`quadrant_facts.out`), so ω is not
componentwise positive in general.)

#### Proof
*(b)* `V' = 2pp' − pp' − Up'' = pp' − Up''`. Insert (2.8):
`V' − (η−μ)V = pp' + ζUp − (η−μ)p² = p[p' + μp + ζU − ηp] = p[k − y + ζU − ηp] = pK`.
Using instead `p'' = −μp' − (x−y)`: `V' = pp' + μUp' + U(x−y) = −μV + p(p'+μp) + U(x−y) = −μV + V₃`.
Subtracting, `V₃ = pK + ηV`.
`Q' = p'p'' − pp'''`; differentiate (2.8): `p''' = (η−μ)p'' + η'p' − ζ'p − ζp'`, so
`Q' − (η−μ)Q = p'p'' − (η−μ)p'² − η'pp' + ζ'p² + ζpp' = p'[(η−μ)p' − ζp] − (η−μ)p'² + ζpp' + p(ζ'p − η'p') = pW`.
The state forms `Q = p'(k−y) + p(x−y)` and `Q' = xV − AQ + kψ p'` follow by substituting
`p'' = −μp' − (x−y)`, `p''' = −μp'' − xp + (x−y)`.
(All of these are verified symbolically, including the explicit multiples of the constraint O off the cycle:
`sym_check.py`.)

*(a)* Let t_a be the trough and t_b ∈ (t_a, t_a+T) the crest, the only zeros of p. At a zero of p, (4.3) gives
`y − k = (x−k)/(1+η)`, which has the sign of U ≠ 0; hence `V = −Up' = U(y−k) > 0` and `Q = p'² = (y−k)² > 0`
at both zeros. Put `E(t) = exp ∫_{t_a}^t (η−μ)`. By (b), `(V/E)' = pK/E`.
On [t_a,t_b], p ≥ 0 and K > 0, so V/E is nondecreasing and `V(t) ≥ E(t)V(t_a) > 0`.
On [t_b, t_a+T], p ≤ 0, so V/E is nonincreasing and `V(t)/E(t) ≥ V(t_a+T)/E(t_a+T) > 0`.
The same argument with W>0 gives Q>0.
Curvature: velocity (p,p'), acceleration (p',p''), so `U'p'' − p'U'' = pp'' − p'² = −Q < 0`. A regular simple closed
C² plane curve whose curvature never vanishes is strictly convex (the tangent angle is strictly monotone and turns
by exactly 2π by the Umlaufsatz). Star-shape: (U,p) ≠ (0,0) on the cycle because V>0, and
`d/dt arg(U+ip) = (Up' − p²)/(U²+p²) = −V/(U²+p²) < 0`.

*(c)* log V and log Q are smooth and T-periodic, so `0 = ⟨(log V)'⟩ = ⟨η−μ⟩ + ⟨pK/V⟩`, likewise for Q; λ_c = ⟨η⟩−μ
is Prop. 2.1 of v2.0. `dU = p dt` gives the line-integral forms. `(log V)' = η₀ − μ` and `(log Q)' = η₁ − μ` are
state identities by (b), so their means vanish; subtracting the two expressions for V'/V, Q'/Q gives
`η − η₀ = −pK/V`, `η − η₁ = −pW/Q`. Solving the linear system {`ζU − ηp = K + H`, `ζp − ηp' = x − y`}
(determinant −V ≠ 0) gives `η = (V₃ − pK)/V`, `ζ = (Q − p'K)/V`. The κ-equation is `κ'/κ = K'/K − V'/V` with
(3.5) and `1 + 2η₀ − μ − pκ = 1 + η + η₀ − μ`.

*(d)* Fix (α,β) ≠ 0 and ℓ := αU + βp. By strict convexity ℓ has exactly one maximum and one minimum on the cycle
and every intermediate level is attained at exactly two points A (where ℓ' > 0) and B (where ℓ' < 0).
`q_A − q_B = s(β, −α)` with s ≠ 0; s depends continuously on the level and on the direction and never vanishes, and
for (α,β) = (1,0) the point A is the upper one, i.e. s < 0; hence s < 0 always. Theorem 3.1 with base A and with
base B, added: `(ζ_A − ζ_B)(U_A − U_B) − (η_A − η_B)(p_A − p_B) > 0`, i.e. `s[(αη+βζ)_A − (αη+βζ)_B] > 0`, so
`(αη+βζ)_A < (αη+βζ)_B`. Therefore `T⟨(αη+βζ)(αp+βp')⟩ = ∮(αη+βζ) dℓ = ∫ [(αη+βζ)_A − (αη+βζ)_B] dℓ < 0`.
The quadratic form in (α,β) is negative definite; `⟨ηp'⟩ = ⟨ζp⟩` because `ηp' = ζp − y'`. The individual
orderings (y, and U at equal p) follow from the two one-sided inequalities: at equal U,
`η_+(p_+ − p_−) < y_− − y_+ < η_−(p_+ − p_−)`; at equal p, `ζ_l(U_r − U_l) < y_r − y_l < ζ_r(U_r − U_l)` and
`U_r > U_l` because the two points never coincide and at p=0 they are the crest and the trough.

*(e)* Rows 1–3: `V' = −μV + V₃` (b); `V₃' = −V₃ + Q + kψ(U)p` and `Q' = xV − AQ + kψ(U)p'` are state identities
(`sym_check.py`). Dividing by V gives the Riccati defects. ∎

**General principle behind (b),(c)** (also proved in `sym_check.py`): with `L₂ = D² + (μ−η)D + ζ` (the second-order
factor of (2.9)) and `V_{f,g} := f g' − f' g`,
`V_{f,g}' = (η−μ)V_{f,g} + f L₂[g] − g L₂[f]`, and `L₂[1] = ζ`, `L₂[U] = K`, `L₂[p] = 0`, `L₂[y−k] = W`.
V is the pair (p,U), Q the pair (p,H). Whenever `V_{f,g} > 0` on the cycle, `λ_c = ⟨(g L₂[f] − f L₂[g])/V_{f,g}⟩`.

### Proposition CD ("convex-Dulac" sufficient criterion; PROVED, but no certificate found)
Let F(U,p,y,η,ζ) be C¹ near a cycle and suppose on the cycle
(1) at equal U, `F_U(upper point) ≤ F_U(lower point)` (true whenever F_U is nondecreasing in η and y, nonincreasing
    in p and independent of ζ);
(2) at equal p, `F_p(trough-side point) ≤ F_p(crest-side point)` (true whenever F_p is nondecreasing in U, y, ζ
    and independent of η);
(3) `η − μ + F_y·(x−y) + F_η·η' + F_ζ·ζ' ≤ 0`, not identically 0.
Then λ_c < 0. *Proof.* `λ_cT = ∮(η−μ+𝔻F)dt = ∮F_U dU + ∮F_p dp + ∮(3)dt`; the first two are ≤ 0 by (1),(2)
exactly as in D(d). ∎
This lets the favourable part of 𝔻F be discarded; in Green form the discarded terms are `−∬ tr(M·D²h)` with M
positive semidefinite. Tested ansatz `F = N(U)η + M(p)ζ + (n₁(U)+m₁(p))y + S₀(y)` with N,M,n₁,m₁ nondecreasing:
feasible on 2 of 10 cycles only (`cdkit.out`); so this ansatz is **refuted as a universal certificate**.

### Proposition B1 (by-product asked for in the lane prompt; PROVED modulo brief item (b) for the statement at E)
(i) For any Λ that is C¹ on a neighbourhood of a cycle, `⟨𝔻Λ + η − μ⟩ = λ_c`. Since an attracting cycle exists
(Smith), no Λ on D∖{E} can have `R ≥ 0` everywhere with R ≢ 0: **the sign must be R ≤ 0**.
(ii) If Λ is differentiable at E then `R(E) = η_E − μ = γ_E − A = 2σ > 0` (σ = Re λ_u). So no Dulac function that is
regular at E exists: in particular **`B = exp(β(U))` and `B = B(U,p)` polynomial (or any C¹ positive B, or any B
vanishing at E) cannot work**; B must blow up at E at least like r^{-2} (B ~ r^{-2s} gives R → 2σ(1−s) at E).
For `B = exp(β(U))` there is the additional failure `R = η − μ` at every point of the axis p=0; on all 30 computed
cycles `η(crest) > μ` (`quadrant_facts.out`), so it fails on the cycle itself as well.
LP confirmation on sampled h: `byproduct.out` (best margin exactly −2.000σ for exp(poly₃(U,p)); both signs infeasible
for exp(β(U))). Inverse polynomials `B = 1/Q(U,p)`, deg Q ≤ 4, are feasible on the annulus r ≥ 0.01 for
(μ,ρ) = (0.5,1.2), (1,1.02) only; not found for the other 5 sets tried (`byproduct2.out`; the LP is poorly
conditioned, so this is weak evidence either way).

---------------------------------------------------------------------------------------------------------------

## 2. Numerical tests of §1 (all run BEFORE the proofs were written)

| test | command | result |
|---|---|---|
| all algebraic identities, exact | `python sym_check.py` | 19/19 residuals 0 (`sym_check.out`) |
| Theorem D on 20 fresh cycles, μ∈{0.1,0.3,0.7,1}, ρ∈{1.001,1.1,2.5,8,20} | `python verify_theorem.py` | ALL OK (`verify_theorem.out`): minV/r² ∈ [2e-5, 0.95] > 0, minQ/(p²+p'²) > 0 (as small as 1e-20 at ρ=20 but positive), identity residuals ≤ 2e-8, the three values of λ_c agree to ≥ 6 digits, ⟨η₀⟩−μ, ⟨η₁⟩−μ ≤ 1e-9, sign pattern exact, ⟨ηp⟩<0, ⟨ζp'⟩<0, det>0, branch orderings >0 |
| V>0, Q>0 on the whole disc (W^u(E) spirals, 10 sets, 42k–156k samples each) | `python explore1.py` | min V/r² ≥ 1.9e-4, min Q/(p²+p'²) ≥ 5.6e-6 (`explore1.out`) |
| quadrant sign facts | `python quadrant_facts.py` | see §5 |

Numerical trap found and fixed: `p' = −μp − H` loses all relative accuracy when p ≈ λ (deep troughs, ρ ≳ 8);
use `p' = k(Y − Z)` in the log variables, and `Q = p'(k−y) + p(x−y)`. Before the fix two ρ=20 cycles falsely
showed Q<0 at the 1e-18 level.

Near-Hopf closed forms (derived from the centre-manifold coefficients of v2.0 Thm 8.1 and the appendix l₁;
checked at ρ=1.001 for μ = 0.3 and 1 to 3–4 digits; status: numerically supported asymptotics, not written as a proof):
with ε = crest amplitude and D = μ²+6μ+1,
`λ_c ≈ −Aμ²ε²/(2(A²+μ)D)`, `⟨ηp⟩ ≈ −μ²ε²/D`, `⟨ζp'⟩ ≈ −μ²(μ²+4μ+1)ε²/(2D)`,
so `λ_c/⟨ηp⟩ → A/(2(A²+μ))` at Hopf. This ratio is NOT a lower bound away from Hopf (μ=1: 0.200 → 0.189 at ρ=1.1
→ 0.013 at ρ=20), so "λ_c ≤ c(μ)⟨ηp⟩" with the Hopf constant is **refuted**.

---------------------------------------------------------------------------------------------------------------

## 3. Step 1 — LP discovery on the attractor disc

Data (`gen_data.py`): 10 sets (μ,ρ) ∈ {(0.05,1.02),(0.05,1.5),(0.05,5),(0.2,1.05),(0.2,2),(0.5,1.2),(0.5,3),(1,1.02),
(1,1.5),(1,5)}; W^u(E) spirals from r=1e-7 to the cycle (1–16 orbits per set so the spirals fill the disc),
integrated in log variables (U, log Z, log Y, log η, log ζ, log K, log W, log Λ_E) so every gauge keeps relative
accuracy next to E and in troughs (U down to −550). Consistency of pulled-back K, W, Λ_E with their state formulas:
≤ 5e-8 relative.
LP (`cplp.py`): maximise t subject to `Σ c_i 𝔻Λ_i + η − μ ≤ −t` on *all* samples (cutting-plane loop, so the
reported margin is exact on the full sample; margins in units of σ = Re λ_u; note t ≤ |λ_c| ≈ 2σ near Hopf).
An earlier attempt that solved on a subsample only over-fitted badly (full-sample margins −3σ…−40σ); discarded.

Feasibility per set (t*/σ, sign R ≤ 0; `lp_run2.out`):

| family (features = Λ_i) | # | (.05,1.02) | (.05,1.5) | (.05,5) | (.2,1.05) | (.2,2) | (.5,1.2) | (.5,3) | (1,1.02) | (1,1.5) | (1,5) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| G: log ζ, η, K, W, Λ_E, x, y, λ−p, 1+η | 9 | −2.00 | −2.00 | −2.00 | −2.00 | −2.00 | −2.00 | −2.00 | −2.00 | −2.00 | −2.00 |
| a: log V, U, p, H | 4 | −2.69 | −2.15 | −1.51 | −1.92 | −1.32 | −1.44 | −0.78 | −1.90 | −0.96 | −0.72 |
| b: a + log Q | 5 | −1.41 | −2.05 | −1.46 | −0.46 | −1.01 | −0.44 | −0.77 | −0.25 | −0.56 | −0.72 |
| c: b + log ζ, log η | 7 | −1.19 | −2.00 | −1.46 | **+0.32** | −0.99 | **+0.76** | −0.16 | **+1.12** | **+0.49** | −0.35 |
| d: c + log y, log(λ−p), log(1+η), log Q₀ | 11 | −0.30 | −0.86 | −0.31 | **+1.17** | **+0.22** | **+1.43** | **+0.41** | **+1.62** | **+1.27** | **+0.34** |
| e: d + log K, log W, log Λ_E | 14 | −0.28 | −0.80 | −0.27 | +1.18 | +0.27 | +1.56 | +0.41 | +1.64 | +1.29 | +0.35 |
| f: e + {p/(λ−p), H/Λ_E, K/Λ_E, W/K, K/V, W/Q, W/V, Q/V, Λ_E/V, η, ζ, x} | 26 | **+1.69** | **+0.26** | **+0.05** | +1.88 | +0.64 | +1.83 | +0.86 | +1.99 | +1.66 | +0.78 |

(Q₀ = p² − 2σUp + ζ_E U², the quadratic form with 𝔻Q₀ = 2σQ₀ for the linearisation at E.)
- **Opposite sign (R ≥ 0): infeasible in every family and set**, as it must be by Prop. B1(i).
- **Pure gauge family G is infeasible on all 10 sets with margin exactly −2.0000σ** (full samples, `lp_G.out`): for
  every coefficient choice max R ≥ 2σ, attained in the limit at E. Reason (Farkas dual, `farkas.out`): there is a
  probability measure supported on ≤ 6 directions in the linear regime (r ~ 1e-8…1e-5) under which every
  `𝔻 log F_i` averages to 0 while η − μ = 2σ. The gauge logs K, W, Λ_E have O(1) angular oscillation `S/G − γ_E` at E
  (S/G vanishes on a ray and exceeds γ_E elsewhere) and the remaining logs have 𝔻 → 0 at E, so no combination
  lowers R below 2σ in every direction. The Wronskian logs (log V, log Q, log Q₀) have *no* O(1) oscillation at E
  (each tends to 2σ), which is why they carry the solution: in every feasible solution
  `−(c_V + c_Q + c_{Q₀}) = 1.3…3.0 > 1`, i.e. B ~ r^{-2s} with s>1, and `|c_K|,|c_W|,|c_{Λ_E}| ≤ 0.23`, usually with
  c_K and c_{Λ_E} of opposite sign (the combination log(K/Λ_E) → 0 at E). The remark in brief 3(d) that B "must behave
  like r^{-2}" holds only inside the η-cancelling ansatz; for fixed ρ>1 a stronger singularity is allowed and is
  what the LP picks.
- Greedy forward selection from the 26-feature pool starting from family a (`lp_greedy.out`): after 9 additions
  (13 features) the three μ=0.05 sets are still infeasible (−0.16, −0.13, −0.05). **Smallest feasible family found
  for all 10 sets: 26 features.** For μ ≥ 0.2 five features {log V, log Q₀, U, p, H} already work on 6 of 7 sets
  (greedy run used every second sample).
- **Single coefficient vector** (family f, `lp_shared.out`): all 10 sets −0.95σ (infeasible); μ=0.05 (3 sets) −0.43;
  μ=1 (ρ = 1.02, 1.5, 5) **+0.59** (feasible); the four near-Hopf sets +0.01 (marginal). No simple dependence of
  the per-set optimal coefficients on (μ,k) is visible (e.g. c_V ranges −1.2…−2.7, c_{log ζ} from −0.3 to +23).
- Polynomial corrections in (U,p,H) up to degree 2–3 do not rescue the small families (`lp_run1.out`).
- Adding the two-point monotone-cone terms of Prop. CD to the disc LP changes nothing (`lp_run3.out`): the binding
  constraints sit next to E where those terms vanish.

Caveat recorded so nobody repeats it: **cycle-only LPs are meaningless as evidence.** On a single known cycle with
λ_c<0 the cohomological equation `𝔻Λ = λ_c − (η−μ)` is solvable, and already 10 quadratic state monomials reach
R ≤ −0.99|λ_c| on each of the 10 cycles (`cd_run2.out`). Only universality in (μ,k) or validity on a region
carries information.

**Step 2 was therefore not reached**: no structured family with guessable closed-form coefficients appeared, so
there was no candidate to re-test on 20 new sets (the 20 new sets were used for Theorem D instead).

---------------------------------------------------------------------------------------------------------------

## 4. Step 3 — why the approach stalls (sharpest obstruction)

**O1. Every state-function Dulac candidate needs a two-sided bound on κ = K/V, and the missing side is an upper
bound on the descending half.** For any C¹ function Λ of the state (U,p,y) (this includes log V, log Q, log Q₀,
log Λ_E, polynomials, …), Theorem D(c) gives exactly

    R = 𝔻Λ + η − μ = E_Λ(z) − p·κ ,        E_Λ(z) := 𝔻Λ + η₀ − μ   (explicit) .

So R ≤ 0 ⟺ `κ ≥ E_Λ/p` where p>0 and `κ ≤ −E_Λ/|p|` where p<0. Functions of (z,η,ζ,K,W) do not escape:
(η,ζ) are affine in κ and W = W₀(z) + κQ, so the pair (z,κ) is closed under the flow and any such Λ gives
R = R(z,κ), which must be ≤ 0 for the true κ(z).

**O2. Convexity only bounds κ from below on a cycle.** For two points z₁,z₂ of one cycle, Theorem 3.1 with base z₁
reads `κ₁ > (Π₁⁰(q₂) − y₂)/|F₁∧(q₂−q₁)|`, where Π₁⁰ is the plane through z₁, z₁' and E (the κ=0 plane) and
F₁∧(q₂−q₁) < 0 for every q₂ because the projection is strictly convex and clockwise (Theorem D(a)). K>0 is the case
z₂ = E, and W>0 is `κ > −W₀/Q`. An upper bound on κ₁ would need an attractor point on the *outer* side of the
tangent line at q₁; on the outermost cycle there is none. The only upper bounds available are dynamic and crude:
`ζ < y` (new, one line: `(y−ζ)' + (y−ζ) = ηζ > 0`), i.e. `κ < (y−ζ₀)/|p'|` where p'<0; ζ>0, i.e. `κ < ζ₀/p'` where
p'>0; and `η < λ − p` for μ ≤ 1/2 (`(λ−p−η)' + μ(λ−p−η) = (y−ζ) + η(η+1−2μ)`); both bounds hold on all disc
samples (`lp_G.out`). None is asymptotically sharp at E:
they give κ = O(1/r) while the truth is κ → q(U,p)/Q₀(U,p) = O(1), so they cannot reproduce the first Lyapunov
coefficient and fail in the whole near-Hopf regime.

**O3. Quadrant picture on cycles** (`explore3.out`, `quadrant_facts.out`; 30 cycles). With quadrants QI (p>0,p'>0),
QII (p>0,p'<0), QIII (p<0,p'<0), QIV (p<0,p'>0):

| residual R for Λ = | QI | QII | QIII | QIV |
|---|---|---|---|---|
| 0 (R = η−μ) | **<0** | mixed | **>0** | mixed |
| −log V, −log Q (R = −pκ, −pϖ) | ≤0 (proved) | ≤0 (proved) | **>0** (proved) | >0 (proved) |
| −log\|p'\| (R = ζp/p', exact) | ≥0 | ≤0 | ≥0 | ≤0 |
| log ζ (R = x/ζ − A) | mixed | >0 | mixed | <0 |
| log η (R = ζ/η − 1) | mixed | >0 (mostly) | mixed | <0 |
| log K, log W | mixed | mixed | mixed | mixed |

Every Wronskian `V_f` with f ∈ span{1,U,H,p} yields a *lower* bound for η on QIII (source `p·L₂[f]` with
L₂[f] ∈ cone{ζ,K,W}); there η exceeds μ on all 30 cycles (η−μ up to +21 at ρ=20) and nothing in the toolbox
bounds it from above. QIII (just after the X crest) is where any certificate has to spend the credit accumulated
on QI; "η<μ on QI", "η>μ on QIII", "η(crest) > μ > η(trough)" hold on all 30 cycles but are **unproved**
(and false on the disc near E, where η−μ ≈ 2σ > 0 in every quadrant).

**O4. Tests of specific simple ideas — all refuted.**
- `B = e^{G(U)}/V`: needs `sup_{p<0} κ(U,·) ≤ inf_{p>0} κ(U,·)`; violated in 26–60 of 60 U-bins on every set, inside
  the linear regime as well (`explore2.out`). On cycles κ_up > κ_low for U ≲ 0 but κ_up < κ_low near the crest
  (`explore4.out`).
- `B = e^{G(U)} V^{-a} Q^{-(1−a)}`: on cycles a = 2/3 gives the right ordering at every tabulated U for 4 of 5
  cycles examined, but for (0.05,1.5) no constant a works (needs a>0.845 and a<0.837); on the disc it cannot work
  near E (the angle function aκ+(1−a)ϖ would have to be constant, impossible for σ≠0 because Q₀ has a cross term).
- "⟨η(λ−p)⟩ ≤ k" (would imply the conjecture via ⟨ηp⟩<0): **false near Hopf** (requires 2(A²+μ) ≤ A²).
- Half-plane ordering of K/V, pure-gauge families, regular B, exp(β(U)): refuted as above.

**What a richer class would have to contain.** (i) Near E/Hopf: the exact 2-jet of the attractor (the tight class
`B = e^{ℓ(U,p)}/(V − C₃(U,p) + …)` with C₃ determined by D²h(E)), since R must vanish to first order and be
negative at second order with the Lyapunov constant; (ii) globally: a genuinely new upper bound for η (or K) after
the crest. The LP shows such B exists numerically; the obstruction is provability, not existence.

---------------------------------------------------------------------------------------------------------------

## 5. Statements attempted — status table

| # | statement | status |
|---|---|---|
| 1 | Theorem D (a)–(e) | **proved** (§1), tested on 20+10 cycles and symbolically |
| 2 | Prop. CD (convex-Dulac criterion) | **proved**; no certificate satisfying it found |
| 3 | Prop. B1: sign must be R ≤ 0; no Dulac function regular at E; exp(β(U)), polynomial B fail | **proved** (part (ii) uses brief 3(b) for the value of η at E) |
| 4 | V>0, Q>0 on the whole attractor disc D∖{E} | proof is the same given brief 3(b) (K, W > 0 on bounded entire orbits; at p=0, `ζ_E U < H < ζU`); numerically true on 10 discs; **conditional on 3(b)** |
| 5 | A Dulac function exists in a 26-feature explicit family for each of 10 parameter sets | **numerical** (exact on the samples) |
| 6 | A Dulac function in ≤ 14 log-type features for μ = 0.05 | **refuted numerically** (families a–e, greedy) |
| 7 | One coefficient vector for all (μ,k) in family f | **refuted numerically**; feasible for μ=1, ρ∈{1.02,1.5,5} |
| 8 | Pure gauge ansatz of brief 3(d) (logs of W, ζ, η, K, Λ_E, x, y, λ−p) | **refuted** (margin −2σ, dual measure at E) |
| 9 | η<μ on QI, η>μ on QIII, η(crest)>μ>η(trough) on cycles | numerically supported (30 cycles), **unproved** |
| 10 | K/V < 0.6 and W/Q < 1.15 on D | numerically supported, unproved, not used |
| 11 | near-Hopf asymptotics of λ_c, ⟨ηp⟩, ⟨ζp'⟩ | numerically supported (3–4 digits), sketch only |
| 12 | half-plane monotonicity of K/V; e^{G(U)}/V; ⟨η(λ−p)⟩ ≤ k; λ_c ≤ c_H(μ)⟨ηp⟩; ansatz A' | **refuted** |
| 13 | the uniqueness conjecture | **open**; nothing here proves ⟨η⟩<μ |

---------------------------------------------------------------------------------------------------------------

## 6. What to try next

1. Upper bound for η on QIII. Theorem D reduces the conjecture to `∫(κ_up − κ_low)dU > 0`; κ obeys the *linear*
   filter `κ' = kψ(U)/V − ν κ`, ν = 1+η+η₀−μ, with an explicit nonnegative source. A comparison argument on this
   filter (source symmetric in U, weight 1/V and rate ν asymmetric between the branches) is the most concrete
   route to the sign.
2. Use Theorem D(e) in the supersolution lane: `ω = r∧z'` has V,Q>0 and source `kψ(U)(0,p,p')`; modify the vector
   field n = r by a state-dependent correction so that the source `[F,n]∧F` becomes nonnegative (it already is on
   p,p' ≥ 0).
3. In a hypothetical two-cycle configuration the outer cycle *does* give upper bounds on κ for the inner one
   (points on the outer side of inner tangent lines). A contradiction argument that plays an inner against an
   outer cycle could use both sides; not attempted here.
4. Prove items 9 (quadrant signs); with them `Λ=0` on QI and exact `−log|p'|` on QII∪QIV leave only QIII.
5. If a computer-assisted proof is acceptable, the LP certificates of family f (per-parameter, exact on samples)
   are the natural objects to validate with interval arithmetic on a parameter box; near Hopf the tight class
   of §4 must be used instead.

---------------------------------------------------------------------------------------------------------------

## 7. Files (all in this directory) and commands

Python: `/Users/admin/.venvs/analysis-py313/bin/python <script>`; every run < 2 min except `lp_greedy.py`
(≈ 12 min, run through `research-run --cpu 1 --memory-gb 1 --label wh-dulac-greedy`).

- `gen_data.py` → `data_mu*_rho*.npz` (disc samples; `unpack()` gives derived arrays)
- `lpkit.py` (feature library with exact 𝔻Λ_i), `cplp.py` (cutting-plane LP), `cdkit.py` (cycle profiles, ansatz A')
- `lp_run1.py/.out` (subsample LP + full recheck; shows over-fitting), `lp_G.py/.out` (family G, both signs,
  a-priori bounds), `lp_run2.py/.out` (families a–f),
  `lp_run3.py/.out` (with monotone cone), `lp_shared.py/.out` (shared coefficients; the file also contains a
  re-print of lp_run2), `lp_greedy.py/.out`, `farkas.py/.out` (dual obstruction measures)
- `byproduct.py/.out`, `byproduct2.py/.out` (exp(β(U)), exp(poly), 1/poly)
- `sym_check.py/.out` (exact identities), `verify_theorem.py/.out` (Theorem D on 20 fresh cycles),
  `quadrant_facts.py/.out`
- `explore1–4.py/.out` (diagnostics: positivity, half-plane test, quadrant signs, branch tables),
  `cd_run2.py/.out` (cycle-only LPs — see caveat in §3)
