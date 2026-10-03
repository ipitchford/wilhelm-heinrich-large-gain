# VERIFY — adversarial referee report on lane `hamiltonian` (2026-10-03)

Scope: the three statements the lane marks as proved (Theorem A, Lemma 1, Theorem C) plus Prop. 2, which carries the
headline. Method: (1) independent numerics (own integrator `verify/whref.py`, own quadrature `verify/todaq.py`,
30-digit tanh-sinh in `verify/v1_theoremA.py`; none of the lane's scripts are imported; `whlib` is used only once,
as a cross-check target); (2) line-by-line reading of NOTES §1, §2, §3.1, §3.2, §4.1; (3) search for circularity
and for silently weakened statements. All numerics are floating point; nothing here is an interval certificate.

Python: `/Users/admin/.venvs/analysis-py313/bin/python verify/<script>`; each run is single-process and under 90 s.

## Verdicts

| # | statement | verdict | residual caveat |
|---|---|---|---|
| 1 | Theorem A (kurtosis ≤ 9/5, `Rsel' ≥ 4/5`, bijection, unique simple root `h*(r) ≤ 5r/4`) | **holds** | none |
| 2 | Lemma 1 identities (b), (c), (f) on every cycle, all parameters | **holds** | none |
| 3 | Theorem C (local uniqueness in an energy ball for small `mu`, `rho = 1 + r eps^2`) | **holds**, conditional on Fenichel's theorem exactly as the lane states | one regularity slip (C^1 written where C^2 is needed and available), repaired below; `mu0` not explicit |
| – | Prop. 2 (`M(h) = (1−rho) oint U'^2`, no zero) | **holds** | none |

No numerical counterexample was found; no unrepaired gap was found. No circularity. The scope of Theorem C is
stated honestly by the lane (local in phase space, not global uniqueness).

## 1. Theorem A

### Proof audit (NOTES §3.2)
- *Velocity density.* On the half orbit `P ≥ 0` the two arcs are graphs over `P`; `dtau = dP/|Phi'(U)|`,
  `U = U_∓(h − P^2/2)`. As `P` grows, `phi = h − P^2/2` falls, `U_∓ → 0` monotonically, and `|Phi'(U_∓)|` falls
  because `Phi'` is increasing with `Phi'(0) = 0`. So `f` is increasing. Integrable because the period is finite
  (turning points are simple: `Phi'(U_±) ≠ 0`). Correct. Only convexity is used.
- *(i)* Chebyshev's integral inequality for the two increasing functions `P^2`, `f` on `(0,Pm)`. Correct.
- *(ii)* `f(P) = nu([0,P])` with `nu` = atom `f(0+) > 0` at 0 plus `df`; `nu` is sigma-finite, Tonelli applies;
  `int (Pm − a) dnu = int_0^{Pm} f = T0/2`, so `pi` is a probability measure; `m_2, m_4` as stated. The polynomial
  identity `(6/5) m_2 − 1/5 − m_4 = t(1+t)(1−t^2)/5` (units `Pm = 1`) and `(9/5)g^2 − (6/5)g + 1/5 = (3g−1)^2/5`
  were re-derived by hand and by sympy (`v1`, T1: both residuals 0). Even moments over the full period equal those
  over the half orbit by reversibility. Correct. Equality case = V-shaped potential (`nu` = atom at 0), so 9/5 is
  the supremum.
- *(iii)* `J1' = 3 J0` (integrand vanishes at the endpoints, derivative bounded). `J0' = T0` differentiates across an
  integrable `1/sqrt` endpoint singularity; this is the area formula `d Area{H0 ≤ h}/dh = T0`, standard.
  `(J1/J0)' = 3 − <P^4>/<P^2>^2`. Correct.
- *(iv)* `Rsel ≤ (2/3) Pm^2 = (4/3)h` gives `Rsel(0+) = 0`; integrating (iii) gives `Rsel ≥ (4/5)h`; hence a
  bijection and `h*(r) ≤ 5r/4`. `M3'(h*) = −J0 Rsel'` since `r − Rsel(h*) = 0`. Correct.
- The reduction `oint F3 = r J0 − (2/3) J1` uses `oint G P^2 = 0`, `oint G U P^2 = (1/3) oint P^4`,
  `oint G E P^2 = (1/3) oint E P^4 = (1/3) oint P^4`; each re-derived by hand (substitute `G = −P'`, integrate by
  parts) and checked numerically (`v3`, block D/E).

### Numerical attempts to break it (`verify/v1_theoremA.py` → `v1_theoremA.out`, `verify/v5_extras.py`)
- Toda potential, `h ∈ [1e-6, 1e6]` (14 values, 30-digit quadrature in `U`, a different method from the lane's ODE
  integration): kurtosis rises from 1.500000 to 1.7999988, always below 9/5; first inequality of (ii) has positive
  slack throughout; `<P^2> ≥ Pm^2/3`; `Rsel/h ∈ [0.80001, 1.0]`.
- `Rsel'` by central differences agrees with `(2/3)(3 − kurt)` to 9 digits; minimum over the grid 0.8000008 (at
  `h = 1e6`). `J0'/T0 − 1` and `J1'/(3J0) − 1` are below `1.2e-11`.
- `Rsel` strictly increasing on 121 log-spaced `h ∈ [1e-6, 1e6]`.
- Ten other convex potentials (including `U^4/4` with `Phi''(0) = 0`, the asymptotically linear `sqrt(1+U^2) − 1`
  and `log cosh U`, a strongly asymmetric one, and two that are only C^1): maximum kurtosis `9/5 − 1.5e-12`
  (`log cosh`, `h = 1e6`), never above. Margins printed in `v5_extras.out` (X1).
- Power check: two non-convex potentials (pendulum near the separatrix; `log(1+U^2)`) give kurtosis 2.3–3.5, so the
  test does detect violations when convexity fails.
- `f(P)` strictly increasing on 199 grid points for `h = 0.01, 1, 100`.

Verdict: **holds**.

## 2. Lemma 1 (and Prop. 2)

### Proof audit (NOTES §1, §2)
Reduction to (S) checked: eliminating `y` gives `U''' + A U'' + mu U' + k(e^U − 1) = 0`; with `tau = sqrt(mu) t`
and division by `A mu` this is (S). Every identity was re-derived by hand:
- (a) average of (S). (b) multiplier `P`: `<PQ'> = −<Q^2>`, `<PQ> = 0`, `<GP> = 0`. (c) multiplier `Q`:
  `<GQ> = −<e^U P^2>`. (d) `(rho−1)<P^2> = −rho<G P^2>`, substitute `rho G = −(eps Q' + Q + eps P)`,
  `<Q'P^2> = −2<Q^2 P>`, `<QP^2> = 0`. (e) multiplier `Q'`: `<GQ'> = <e^U P^3>/2`, `<PQ'> = −<Q^2>`.
- (f) `−rho<G P^3> = eps(<P^4> − 3<Q^2P^2>)`; `<Q^2 P> = (eps/2)<Q^3> + (rho/2)<e^U P^3>`; inserting in (d) leaves
  `eps(1−rho)<e^U P^3>`, which (e) turns into `2 eps^2 (1 − 1/rho)(<Q'^2> − <Q^2>)`. Matches the stated (f).
- Conversion to the brief's variables (`P = p/sqrt(mu)`, `Q = p'/mu`, `Q' = p''/mu^{3/2}`, `eps^2 = mu/A^2`) checked.
- Prop. 2: on `P^2/2 + rho Phi = h`, `G = −Q/rho` gives `oint G P^2 = 0`; `oint Q^2 = rho oint e^U P^2 = rho oint P^2`.
  Correct; `M(h) = (1−rho) oint P^2` has no zero for `rho ≠ 1`.
Only periodicity and smoothness of the solution are used. No hidden hypothesis.

### Numerical attempts to break it (`verify/v2_identities.py` → `v2_identities.out`)
Own cycle finder (two-leg return map with variational equations, Newton), period averages by running integrals
inside the ODE (rtol 1e-13), so spiky cycles are handled. Grid `mu ∈ {0.02, 0.3, 1} × rho ∈ {1.01, 1.5, 4, 12}` plus
`(3, 2)` and `(0.02, 1.003)`; troughs down to `U_min = −8357`, `m_t` from 0.997 down to `4.8e-11`.
Defects normalised by the largest term of each identity:

| identity | worst defect over 14 cycles |
|---|---|
| (a) `<e^U> = 1` | 2.4e-11 |
| (b) `<Q^2> = <P^2>` | 3.9e-12 |
| (c) `<Q^2> = rho<e^U P^2>` | 3.9e-12 |
| (d) | 2.6e-13 |
| (e) | 3.2e-11 |
| (f) | 4.8e-13 |
| (f) in the brief's variables | 4.8e-13 |

The worst cases are the deepest cycle (`mu = 0.02, rho = 12`, periodic closure `1.9e-8`). The lane's own figure
("1e-11 on 8 cycles") is consistent with its `exact_identity.out` (max `1.04e-11`). The lane tested only (b), (c),
(f); (a), (d), (e) are now tested too.

Extra (X2 in `v5_extras.out`): evaluated on Toda orbits, the right side of (f) equals `(2/3) eps^2 <P^4>` to
`1e-13`; this derives the selection rule `r = Rsel(h)` from the exact identity alone, independently of the
modified-energy computation of §3.1.

Verdict: **holds**.

## 3. Theorem C

### Proof audit (NOTES §4.1)
*Algebra (re-derived, not re-used; `verify/v3_normal_form.py` → `v3_normal_form.out`).*
- `w_0..w_3` obtained by recursion from the invariance equation agree with the lane's formulas (A).
- `D' = −(1/eps + w_P) D` for `D = Q − w` is an exact consequence of the invariance equation (B); it holds for any
  `Q` because the fast equation is linear in `Q`.
- `X(H~)` has zero coefficients at `eps^0, eps^1, eps^2` and `F3` at `eps^3` (C); `sigma_h = 2(U − U_−(h))` follows
  from the definition because the boundary term carries the factor `h − Phi(U_−(h)) = 0`.
- `oint F3 dtau` by ODE integration equals `r J0 − (2/3) J1` from quadrature, relative difference `≤ 8e-13`, for
  `h ∈ {0.05, 0.7, 4, 25, 150}` and three values of `r` (D).

*Step 0.* Fenichel's hypotheses hold: analytic vector field; critical manifold `Q = −G(U)` is a graph over the
compact disc `K` with smooth boundary; the single fast eigenvalue is `−1` everywhere; `r` enters as a parameter
(trivial slow variable on a slightly enlarged interval). Joint `C^k` dependence up to `eps = 0` is part of the
theorem as stated by Fenichel/Jones. Taylor coefficients are forced by the invariance equation (needs mixed
derivatives of order ≤ 5, fine for `k_s ≥ 6`). I did not re-derive Fenichel's theorem; the lane says the same.

*Step 1.* Correct. A periodic `D` that decays exponentially forward in time vanishes. Consequently every cycle with
projection in `K` is a cycle of the planar field, and conversely by local invariance. The linearisation is block
triangular in `(U, P, D)`, so the three Floquet multipliers are `1`, the planar multiplier, and
`exp(−oint(1/eps + w_P))`; the planar one is the weak multiplier `m_t`.

*Step 2.* Correct, with one repair. The text says `R4` and `Dsp` are `C^1`. That suffices for existence,
uniqueness of the zero on the annulus, `c* = h* + O(eps)` and `∂_c Dsp(c*) < 0`, i.e. for conclusions (1), (2) and
"hyperbolic attracting". The remainder `O(eps^4)` in (3) needs `∂_c Dsp` Lipschitz in `(c, eps)`, i.e. `Dsp ∈ C^2`.
This is available under the hypothesis already in force: `R4` is the fourth-order Taylor remainder of a `C^{k_s}`
function, hence `C^{k_s−4} ⊇ C^2` for `k_s ≥ 6`, and the flow and return time are `C^{k_s}`. **Repair: replace
"C^1" by "C^2" in Step 2.** A second cosmetic point: the coordinate `c = H~|_S` depends on `eps`, so the argument
should be run on a fixed interval slightly larger than the image of `[h_a, h_b]`; nothing changes.

*Step 3.* Correct. In polar coordinates the field is `rho' = rho R`, `theta' = Theta` with `R, Theta ∈ C^{k_s−1}`
up to `rho = 0` and `eps = 0`, `Theta(0, ·) = −1` at `eps = 0`, so `g ∈ C^{k_s−1}` on the closed set. Vanishing of
`∂_eps^j g(a, 0)`, `j ≤ 2`, for each `a > 0` extends to `a = 0` by continuity; Taylor's formula with the uniform
bound on `∂_eps^4 g` (needs `k_s ≥ 5`) gives `g = eps^3(g3 + eps g4)`; `g3 > 0` on the compact set because
`Phi(a) ≤ h_a < h*(r1) ≤ h*(r)` and `g3(0) = pi r`. Checks: `g3(a)/pi → 1` as `a → 0` for `r = 1` (0.905, 0.978,
0.995, 0.998 at `a = 0.3, 0.1, 0.03, 0.01`; `v3` block F); and the exact linearisation at `E` gives
`Re z = eps^3 r / (2(1+eps^2))`, i.e. `Pi'(0) − 1 = pi r eps^3 + O(eps^5)`, the same number.

*Step 4.* Correct. Two facts are used without being stated, both standard: a closed orbit of a planar flow meets a
transversal arc exactly once; and `S` is a transversal all the way down to the origin (`w(U,0,eps) < 0` for
`0 < U ≤ U_+(h_b+1)`, from `w(0,0,eps) = 0`, `∂_U w(0,0,0) = −1` and `C^1` closeness). The only equilibrium of the
planar field in `K` is the origin (a planar equilibrium lifts to an equilibrium of the full system).

*Independent derivation of the multiplier formula.* The planar multiplier is `exp(oint w_P dtau)`. On a true cycle
`eps oint G = 0` (Lemma 1(a)), `eps^2 oint ∂_P w_2 = −2 eps^2 oint e^U U' = 0`, and to leading order
`oint ∂_P w_3 = r T0 − 2 J0 = M3'(h)`, which equals `−J0 Rsel'` at `r = Rsel(h)`. Verified numerically (`v3`, block
G, to 8 digits). This confirms (3) by a route that does not use the modified energy.

*Circularity.* None. The proof uses neither the v2.0 uniqueness near Hopf nor any stability statement; Theorem A
and Lemma 1 are self-contained.

*Silent weakening.* None found. The statement is local (cycles inside `{H0 ≤ h_b}`), and the lane says so in the
theorem, in §4.2 and in the headline. `mu0` depends on `r1, r2, h_b` and is not explicit.

### Numerical tests (`verify/v4_theoremC.py` → `v4_theoremC.out`; `verify/v5_extras.py` → `v5_extras.out`)
The theorem is asymptotic, so numerics cannot refute it outright; they can contradict its two quantitative
conclusions or exhibit a second cycle in the ball.

- (P0) The referee's cycle finder agrees with `whlib.find_cycle` on three parameter pairs (`T` to 9 digits, `m_t`
  to 8).
- (P1) Fixed `r ∈ {0.3, 1, 3}` (the lane tested at fixed amplitude), `mu = 1e-2 … 1.56e-4`:

  | r | (h_crest − h*)/eps^2 | log m_t / (−eps^3 J0 Rsel') − 1, divided by eps^2 | sup\|Q − w_series\|/eps^4 on the cycle |
  |---|---|---|---|
  | 0.3 | 1.145 → 1.161 | −1.87 → −2.00 | 8.7 → 8.8 (the script prints 11.6 at the smallest `mu`; rerun with `atol = 1e-17` gives 8.82, so that entry is an integration-tolerance artefact amplified by `1/(mu eps^4)`) |
  | 1 | 4.19 → 4.37 | −3.85 → −4.32 | 65 → 68 |
  | 3 | 17.7 → 19.4 | −8.9 → −10.8 | 582 → 640 |

  All three columns converge to constants: the crest energy is `h* + O(eps^2)` and the multiplier formula has
  relative error `O(eps^2)`, both better than claimed (`O(eps)`, `O(eps^4)` absolute), and the cycle sits on the
  slow manifold to the predicted order (Step 1). The constants grow with `r`, as the lane notes.
- (P2) Amplitude-parametrised family (own implementation): for 26 crest amplitudes `a ∈ [0.1, 5]` (crest energy up
  to 142), `rho(a)` is strictly increasing and `m_t < 1`, at `mu = 0.02` and at `mu = 1e-3`. So each `rho` in the
  range has one cycle with crest in the ball.
- (P3) Planar squeeze at `mu = 0.02`: crest sequences started near `E` (increasing) and near the ball boundary
  (`H0 ≈ 8.6` for `r = 0.52`, `H0 ≈ 12.0` for `r = 2.08`; decreasing) are monotone and converge to the same crest
  (difference `2e-9`), which is the Newton cycle. On a two-dimensional invariant manifold this excludes any other
  cycle between the two starting points.
- (X3) On the brief's grid `mu ∈ {0.02, 0.3, 1} × rho ∈ {1.01, 1.5, 4, 12}` every computed cycle has `m_t < 1`.
  The leading-order formula is accurate only where it should be (`mu = 0.02, rho = 1.01`: ratio 0.955) and is off
  by a factor 1.5–7 at fixed `rho ≥ 1.5`, where `r` is 2–570 and the regime is not Hamiltonian. This is outside the
  theorem's hypotheses and agrees with the lane's own caveat; it is not a contradiction.

Verdict: **holds**, as a theorem conditional on Fenichel's persistence theorem, with the one-word repair above.

## 4. Things a reader should not over-read
- Theorem C is not a uniqueness theorem for small `mu`. A large cycle (crest energy between `h_b` and `~2/mu`) is
  not excluded by anything in the lane; the lane states this.
- Theorem A is a statement about the unperturbed Toda centre; by itself it says nothing about the full system.
- Lemma 1 holds on every cycle but is a selection rule, not an amplitude bound.
- Sections 4.2a, 4.2b, 5 and 6 of NOTES are outlines or numerics and were not refereed here.

## 5. Files
| file | content | run time |
|---|---|---|
| `verify/whref.py` | referee's integrator, return map with variational equations, cycle finder, running-integral averages | – |
| `verify/todaq.py` | Toda level integrals by Gauss–Kronrod in `U`, `h*(r)` | – |
| `verify/v1_theoremA.py`, `.out` | Theorem A: identities, Toda scan, other potentials, power check | 52 s |
| `verify/v2_identities.py`, `.out` | Lemma 1 (a)–(f) on 14 cycles | 9 s |
| `verify/v3_normal_form.py`, `.out` | slow-manifold recursion, `D` equation, `F3`, averages, `g3(0)`, divergence route | 1 s |
| `verify/v4_theoremC.py`, `.out` | fixed-`r` asymptotics, amplitude family, planar squeeze | 81 s |
| `verify/v5_extras.py`, `.out` | sharpness margins, (f) on Toda orbits, brief's grid | 11 s |
