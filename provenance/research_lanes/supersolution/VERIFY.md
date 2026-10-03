# Lane `supersolution` — adversarial verification (2026-10-03)

Referee brief: refute the six statements the lane marks as proved. Method: (1) independent numerics (own integrator,
own variables, none of the lane's scripts or its `cycles.pkl`); (2) line-by-line proof check, including every imported
result; (3) search for circularity and silently weakened statements.

**Outcome: not refuted.** All six statements hold as stated. No numerical counterexample on 65 independently computed
cycles; every algebraic identity verified symbolically; every proof step checked. One omitted justification
(positivity of the strong-stable Floquet vector in Proposition B) is repaired below in one line from the paper's own
Perron argument. Four reporting defects are listed in section 4; none affects a proved statement.

All numerics are floating point. Nothing here is an interval certificate.

---------------------------------------------------------------------------------------------------------------------
## 1. Verdicts

| # | statement | verdict | basis |
|---|---|---|---|
| 1 | Theorem R (S1 > 0, S4 > 0, both ODEs, star-shape, convexity, four exponent formulas) | holds | proof steps re-derived; `sym_identities.py`; `check_R.py`, `check_R2.py`, `check_sweep.py` (65 cycles) |
| 2 | Proposition T (radial-ratio monotonicity for nested cycles) | holds | proof steps re-derived, including the shifted centre; `sym_identities.py`; `check_T.py` |
| 3 | Lemma A0 (direction-field bound, phi drops, equality at truth, concavity) | holds | proof steps re-derived; `check_A0.py` (14 000 random trials + adversarial optimisation); `check_floquet.py` |
| 4 | Lemma A1 (first-order rigidity at Hopf, under the stated expansion hypothesis) | holds | proof steps re-derived; coefficients solved by hand and by sympy; `check_A1.py` |
| 5 | Proposition B (Target B identities; beta > 1, kappa > (1-mu) beta for mu <= 1) | holds (one justification supplied) | proof re-derived; `check_B.py`, `check_floquet.py`, `check_sweep.py` |
| 6 | Scalar reduction ((R_K), constraint line, S1 = Q_E + U Lam_E, alpha identity) | holds | `sym_identities.py`; `check_R.py`, `check_sweep.py` |

---------------------------------------------------------------------------------------------------------------------
## 2. Proof audit

Notation as in the lane's NOTES. `f = (p, p', y')`, `H = y - k`, `p' = -mu p - H`, `O = zeta p - eta p' - y'`.

### Claim 1 — Theorem R
- `S1' = pp' - U p''`, `p'' = -mu p' - y'` gives `S1' = pp' + mu U p' + U y'`. Expanding `(eta-mu)S1 + pK` and using
  `O = 0` and `-H = p' + mu p` gives the same expression. Correct. (Also `S1' = S2 - mu S1` identically, with no eta.)
- Zeros of p: `S1 = U(y-k)`. Paper (4.1) `x - k = (1+eta)(y-k)` with `1+eta > 0` gives `sign(y-k) = sign(U)`; the paper's
  Theorem 4.1 proof gives `U != 0` there. So `S1 > 0` at crest and trough.
- Two-arc argument: variation of constants forward from the trough on `{p > 0}` and backward from the next trough on
  `{p < 0}`. Both integrands have the sign needed because `K > 0` (paper Cor. 3.2). The two arcs exhaust the period
  because p has exactly two simple zeros (paper Thm 4.1). No gap.
- `S4' - (eta-mu)S4 = p'[p'' - (eta-mu)p' + zeta p] + p(zeta' p - eta' p') = pW` by the tangent equation. At `p = 0`,
  `S4 = (k-y)^2 > 0`. Same two-arc argument with `W > 0` (paper Lemma 3.3). Curvature `= -S4/|(p,p')|^3`, and the
  projection is a regular Jordan curve (paper Thm 4.1), so nonvanishing curvature gives strict convexity.
- `d theta/dt = -S1/r^2 < 0` and winding number -1 give strict star-shape about E.
- Shifted centre: `s_c = S1 - c p'`; `s_c' = (eta-mu)s_c + p(K + c zeta)` by the tangent equation; at `p = 0`,
  `s_c = (y-k)(U+c)`, positive at the crest, and at the trough because `U_min + c < 0`. `K + c zeta > 0`. Correct.
- Averages: `S1`, `S4`, `s_c` are positive and periodic, so the log-derivative averages vanish. Correct.
- Imported hypotheses: K > 0 and W > 0 come from paper Thm 3.1 (no stability input). No circularity.

### Claim 2 — Proposition T
- Polar system: `r'/r = sin(cos - mu sin - m)`, `theta' = -sigma`, `sigma = S1/r^2`. Verified by hand and by sympy.
- `dXi/dm = sin/sigma^2`, and `Xi(m_2) - Xi(m_1) = sin (m_2 - m_1)/(sigma_1 sigma_2)` exactly (numerator and
  denominator affine in m; cross term `s0 - n0 cos = 1`). Verified by sympy.
- `m_2 > m_1`: with base `z_1` on the inner cycle, paper Thm 3.1 against E gives `d_1 > m_1`; against `z_2` on the same
  ray gives `(R_2 - R_1) d_1 < y_2 - y_1`. Then `y_2 - k > (R_2 - R_1) d_1 + R_1 m_1 > R_2 m_1` uses `R_2 > R_1`, which
  follows from nesting plus star-shape (Theorem R). Arbitrary phase pairs are allowed by Thm 3.1. Correct.
- Shifted centre: for `c > 0`, `(h - h(-c,0))/r` is nondecreasing along rays (convexity of the envelope h, paper
  Thm 5.1 with S = E + two cycles) and `-(k - h(-c,0))/r` is strictly increasing since `h(-c,0) < h(0,0) = k`. A direct
  proof also works: the same two inequalities with `K_1 + c zeta_1 > 0` in place of `K_1 > 0`. Star-shape about
  `(-c,0)` for both cycles needs `c < |U_min,1| < |U_min,2|`, as stated. Correct.
- The statement is conditional on two cycles existing; it derives no contradiction, and the lane says so.

### Claim 3 — Lemma A0
- Defect of `omega = phi(1, eta~, zeta~)` in `w' = Cw` is `diag(g) omega` with the stated `g_i`. Verified by sympy.
- `C = S J^[2] S` (verified), so the compound monodromy has eigenvalues `{m_s, m_t, m_s m_t}`; since `m_s < min(1, m_t)`
  its spectral radius is `m_t`.
- C is Metzler, so its evolution is nonnegative; `w = e^{-int m} omega` satisfies `w' >= Cw`, hence
  `Phi_C(T,0) omega(0) <= e^{-T<m>} omega(0)` with `omega(0) > 0`; Collatz–Wielandt (upper bound form) gives
  `m_t <= e^{-T<m>}`. Correct direction.
- At the true pair all three `g_i` equal `mu - eta` (verified), so `S = -lambda_c`.
- Each `g_i` is linear minus a positive multiple of an exponential of a linear form in the exponents; concave.
  Minimum and average preserve concavity. Correct.

### Claim 4 — Lemma A1
- Order 0: `S -> min(mu - e0, 1 - z0/e0, A - k_H/z0)`; nonnegativity forces `mu <= z0 <= e0 <= mu`. Correct.
- Order eps: brackets `(-a, (a'+a-b)/mu, (b' + Ab - mu A u)/mu)`; uses `k/mu = A + O(eps^2)`. Averaging the three
  inequalities `bracket_i >= m_1` with `<m_1> >= 0`, `<u> = 0` gives `0 <= <b> <= <a> <= 0`, hence all brackets equal
  `m_1` identically, hence the linearised Riccati system `a' = b - Aa`, `b' = mu A u - Ab - mu a`. Homogeneous matrix has
  eigenvalues `-A ± i sqrt(mu)`, so the periodic solution is unique. Correct.
- Coefficients on the centre plane (`U'' = -mu U`): the four linear equations for `(alpha_1, alpha_2, beta_1, beta_2)`
  were solved by hand and by sympy; the displayed values are correct.
- The hypothesis (constant order-0 term, C^1-uniform expansion) is stated in the lemma. A C^1 state function satisfies it.
  The radial direction also does: `S1 = mu eps^2 + O(eps^3)` uniformly, `S2/S1 = mu + O(eps)`.
- Scope remark (not a gap in the lemma): section 5 of NOTES says "any state-function certificate". The written proof
  covers directions with a constant order-0 limit. Directions singular at E with a phase-dependent order-0 term
  `(e0(t), z0(t))` are covered by a short extension: Lemma A0 for the constant matrix `C_H` (Perron exponent 0 at Hopf)
  gives `S_0 <= 0`; `S_0 = 0` is the equality case of Collatz–Wielandt for a positive matrix, which forces all `g_i`
  equal and `(e0, z0)` the Perron direction `(mu, mu)`.
- The four "dead" directions are also dead by direct arguments independent of A1: `omega^Y` has
  `S = -<(lambda - A - 2p)_+> < 0` because `<p> = 0`; constants need `e < mu`, `z < e`, `z > k/A > mu`; radial has
  `<mu - eta~> = 0` exactly.

### Claim 5 — Proposition B
- `(a, -b, -c)` solves the variational equation; `nu = e^{At}(c,b,a)` solves `nu' = -C^T nu` (sympy);
  `A - gamma = -lambda_c`; hence `lambda_c <nu^.omega> = -<nu^.(omega' - C omega)>` for periodic C^1 omega.
- Defects of `omega^Y` and `omega^c` verified (sympy). Denominators are positive because `P = lambda - p > 0` on
  positive cycles.
- `chi(gamma) = <x a^>/<a^>` from `a''' + A a'' + mu a' + x a = 0` and `(D-gamma)(D-gamma+1)(D-gamma+mu)`.
- `(kappa + beta^2/2 - beta)' = chi(beta) - x` (sympy), so `k = <chi(beta)>`.
- `v = c - (1-mu) b`: `v' = -mu v - x a`. `e^{mu t} v` strictly decreasing and `v = e^{-gamma t} v^` with `gamma > mu`
  force `v^ > 0`. Then `(b - a)' = -v < 0` and the same argument gives `b > a`. Correct.
- **Omitted justification, supplied here.** The proof starts from "(a,b,c) > 0" without proof (the brief flags it as
  unverified). It follows from the paper: the backward variational flow in the coordinates `(a,b,c) = S delta`
  is `d/d(-t)(a,b,c) = B^T (a,b,c)`, whose period matrix is the transpose of the paper's strictly positive adjoint
  period matrix; Perron–Frobenius gives a positive eigenvector with eigenvalue `e^{gamma T}`, i.e. multiplier `m_s`.
  Confirmed numerically (`check_floquet.py`: strong-stable eigenvector sign pattern `(+,-,-)` from the plain 3x3 monodromy).

### Claim 6 — scalar reduction
- `S2 = eta S1 + pK`, `S3 = zeta S1 + p'K` on `O = 0`; `K' = kG1 - (1+eta)K`; substitution gives (R_K); `N = K/S1`
  equation; `S2/S1 = mu + (log S1)'`; `S1 = Q_E + U Lam_E`; `DQ_E = 2 alpha Q_E - 2(p - alpha U) Lam_E`; `Q_E` positive
  definite above Hopf because the unstable pair is complex (a real cubic with positive coefficients has no
  nonnegative real root). All verified by sympy and numerically.
- The scalar comparison for (R_K) is labelled "only SKETCHED" by the lane and is not counted as proved.

### Circularity and weakening
- No step uses stability, uniqueness or the sign of the Lyapunov coefficient. Imported results: paper Prop 2.1, 2.2,
  Thm 3.1, Cor 3.2, Lemma 3.3, Thm 4.1, Thm 5.1, local Hopf expansion. Their hypotheses hold in each use.
- No statement was found weaker than its label, apart from the scope remark under Claim 4 and the items in section 4.

---------------------------------------------------------------------------------------------------------------------
## 3. Numerical falsification attempts

Python: `/Users/admin/.venvs/analysis-py313/bin/python -W ignore <script>` inside `verify/`. Single process, each run
under two minutes. Cycle location uses `whlib.find_cycle` (shared library) with continuation; every profile quantity
is recomputed by `vprofile.py`: state variables `(U, P = lambda - p, y)` with relative error control, forward Riccati,
`K` and `W` from their own ODEs and again from their definitions, and all period averages integrated as ODE components.

| script | receipt | what it tests | result |
|---|---|---|---|
| `sym_identities.py` | `sym_identities_out.txt` | 33 identities used in the six proofs | 0 failures |
| `vcycles.py` | (stdout) | cycles at mu in {0.02, 0.3, 1} x rho in {1.01, 1.5, 4, 12} | found; residual <= 1e-11 |
| `check_R.py` | `check_R_out.txt` | Theorem R at the 12 extremes | signs positive; four formulas agree with `<eta>-mu` to 8e-11 relative |
| `check_R2.py` | `check_R2_out.txt` | centres `c = 0.99, 0.999 |U_min|`; mu = 5, 20 | positive; agreement 4e-10 |
| `check_B.py` | `check_B_out.txt` | Proposition B at the 12 extremes | three exponent formulas ratio 1.0000000; `chi(gamma)` to 8 digits; `<chi(beta)>/k = 1.000000000` |
| `check_floquet.py` | `check_floquet_out.txt` | 3x3 variational monodromy, no Riccati | multipliers, `(+,-,-)` eigenvector, `beta(0)`, `kappa(0)` agree; compound monodromy positive with spectral radius `m_t` |
| `check_A0.py` | `check_A0_out.txt` | random and optimised trial directions on 7 cycles | no trial exceeds `-lambda_c`; optimum returns to the true pair; no concavity violation |
| `check_A1.py` | `check_A1_out.txt` | near-Hopf regression and scores, rho - 1 from 1e-2 down to 1e-6 | coefficients converge to the claimed values; `S/eps` of inexact jets converges to the lemma's `<m_1> < 0` |
| `check_T.py` | `check_T_out.txt` | Thm 3.1 two-phase form; polar equation; nested spiral turns | Delta < 0 on about 4.9e5 phase pairs per cycle (7 cycles); polar equation to 7e-7; monotone ratio in 38 of 38 nested pairs |
| `check_sweep.py` | `check_sweep_out.txt` | 40 random (mu, rho) + 3 cycles at rho = 1.0001 | worst deviation 4e-9 (near Hopf, absolute 1e-13); no nonpositive minimum |
| `probe_delta.py`, `probe_underflow.py` | `*_out.txt` | diagnosis of the only sign anomalies seen | both are artefacts (below) |

Details.

- **Extremes (claim 1).** At rho = 12 the trough reaches `U_min = -8357` (mu = 0.02), `-926` (0.3), `-668` (1).
  `min S1/max S1 >= 4.1e-3`, `min S4/max S4 >= 6e-15`, `min s_c/max s_c >= 1.4e-6`. `whlib`'s return-map exponent is
  noise-limited there (`m_t` down to 1e-13, reported as NaN at mu = 1, rho = 12); the Riccati exponent and the
  backward exponent `<beta> - A` agree to 7 digits in every case, and the 3x3 monodromy agrees wherever `m_t > 1e-3`.
- **Anomaly 1 (not a counterexample).** At mu = 0.02, rho = 12 the integrated `W` and `beta - 1` reach about `-1e-200`
  and `-1e-245`. Their true size in the trough is of order `e^{-800}`, below double range. Moving the absolute
  tolerance floor moves the negative value with it (`probe_underflow.py`), so it is tolerance noise.
- **Anomaly 2 (not a counterexample).** A first version of the Thm 3.1 two-phase test showed 10 positive Delta values
  of size 1e-17 at mu = 1, rho = 12. Cause: forming `p = lambda - P` before differencing loses the relative accuracy
  of `P ~ 2e-10`. Differencing `P` directly removes all ten (`probe_delta.py`, corrected `check_T.py`).
- **Lemma A0.** Trials are `true pair x exp(th . F)` with bounded state features and pure phase functions, 2000 per
  cycle at five scales, and an absolute monomial family. `S(true) + lambda_c` is the grid quadrature error
  (<= 4e-6 absolute). Largest random score is below `-lambda_c` by >= 4.7e-5. Powell from three starts never improves
  on the true pair.
- **Lemma A1.** At mu = 1, eps = 0.04 the fitted coefficients are 0.24998, -0.24999, 0.74994, -0.24999 against
  1/4, -1/4, 3/4, -1/4. For a jet perturbed by the factors (1.3, 1, 1, 0.8), at the smallest amplitude
  (rho - 1 = 1e-6, eps = 0.004 to 0.011) the measured `S/eps` is -0.06360 (mu = 1), -0.06507 (mu = 0.3), -0.09437
  (mu = 0.02) against the lemma's first-order prediction `<m_1>` = -0.06381, -0.06514, -0.09267, computed
  independently on the Hopf circle; the gap halves when eps halves. `omega^Y` and the constant direction converge to
  their predicted negative `<m_1>` in the same way. The exact affine jet gives `S/eps^2 -> -0.375` (mu = 1), i.e. the
  slack is second order, as the lemma says.
- **Proposition T.** Two cycles at one parameter are not available, so the proposition cannot be tested directly.
  The same mechanism was tested on real nested curves: successive turns of a `W^u(E)` spiral and the cycle (this uses
  brief 3(b), which is not proved). In all 38 pairs `log(R_out/R_in)` increases on `(0, pi)` and decreases on
  `(pi, 2 pi)` across every one of the 717 neighbouring grid-angle pairs (step pi/360), and `m_out > m_in`. A first version using linear interpolation showed
  spurious wrong-sign slopes at mu = 0.02; exact evaluation at the target angles by bisection removes them.

---------------------------------------------------------------------------------------------------------------------
## 4. Reporting defects in the lane's notes (none changes a verdict)

1. **Proposition B ratio line.** NOTES says "ratio check 1.000000 on 54/55 cycles, 0.9989 at mu = 1, rho = 10".
   `targetB_out.txt` shows 50 rows at 1.000000, four at 0.999997–1.000014 and one at 0.998865, all at rho = 10. The
   deviation comes from dividing by `whlib`'s noise-limited `lam_M`; with the Riccati exponent the ratio is 1.0000000
   at all twelve extremes, including mu = 1, rho = 12.
2. **beta > 1 receipt.** `radial_check_out.txt` shows `min(beta - 1) < 0` in 10 of 55 rows (down to -6e-10) and
   `min(kappa - (1-mu)beta) < 0` in 10 rows. This is round-off in the `(beta, kappa)` variables where `beta -> 1+`.
   In the variables `delta = beta - 1`, `nu = kappa - (1-mu)beta` (relative error control) both minima are positive
   wherever they are representable (down to 1e-159). NOTES does not mention the negative rows next to Theorem R.
3. **Equation numbers.** NOTES cites the paper's orthogonality and tangent equations as (2.9) and (2.10); they are
   (2.7) and (2.8). (4.1), Thm 3.1, Cor. 3.2, Lemma 3.3, Thm 4.1, Thm 5.1 are cited correctly.
4. **Wording.** "Stability therefore reduces to one explicit upper bound for eta on the falling arc" describes a
   sufficient route, not an equivalence. Obstruction 1 ("any state-function certificate") is wider than Lemma A1 as
   written; see the scope remark under Claim 4 for the extension that closes it.

---------------------------------------------------------------------------------------------------------------------
## 5. Limits of this verification

- Floating point only; agreement to 1e-9 is evidence, not proof. The verdicts rest on the proof audit.
- Proposition T has no direct numerical test (no second cycle exists numerically at any tested parameter).
- The paper's own results (Thm 3.1, 4.1, 5.1, Lemma 3.3) were taken as given, as the brief directs; only the
  two-phase inequality of Thm 3.1 and positivity of K, W were spot-checked numerically.
- `find_cycle` supplies the crest point; closure of the referee's own integration over one period is <= 1e-10.
