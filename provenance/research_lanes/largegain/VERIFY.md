# Lane `largegain` — adversarial referee report (2026-10-03)

Brief: refute the five statements claimed as proved in `NOTES.md` / `theorem.tex`.
Outcome: **not refuted.** No numerical counterexample, no unjustified step found. Verdict per statement at the end.

All scripts are in `verify/` and were written from `theorem.tex` only (the lane's `constants.py`,
`check_lemmas.py`, `test_segments.py`, `test_stage.py`, `lgcycle.py` were not read or reused; `constants.py`
was executed once, unread, to confirm the lane's reproducibility claim — it prints "ALL CONSTANT CHECKS PASSED").
Python: `/Users/admin/.venvs/analysis-py313/bin/python`, single process, every run under 3 minutes.

## 1. What was run

| script | what it tests | result |
|---|---|---|
| `v_consts.py` | every numerical constant and scalar inequality of the proof, mpmath 50 digits; closed forms of `c_i, c_ii, Psi_mu` against quadrature; `omega >= omega_0(k)`, `gamma_E - mu < k^{1/3}+1` against true roots on a (mu,k) grid; the crest side conditions and `C_tau <= 71.79 + 12.28 log(lambda/1000)` on 601 values of lambda in [1e3, 1e15]; `q` decreasing; brute-force check of the bootstrap implication `h(z) > 0.3 kappa_0 k, z >= z_1 => z > z#` on a 6x5x2001 grid | 86 checks, 0 failures |
| `v_moderate.py` | every parameter-free lemma (L5–L10, L12, L13, L17) and the imports F2–F6 on cycles **outside** the theorem's range: mu in {1, 0.3, 0.02} x rho in {1.01, 1.5, 4, 12} (cycle from `whlib.find_cycle`, own dense integration and event finder) | 12 cycles, 579 checks, 0 failures |
| `v_moderate.py mu rho` scan | same, mu in {1, 0.6, 0.3, 0.1, 0.05} x rho in {1.0005, 1.002, 1.005, 1.02, 1.05, 1.1, 1.2, 1.35, 2, 3, 6, 8, 20, 40} | 70 cycles, 3308 checks ok, 2 flagged (integration noise, see 3) |
| `v_largegain.py` | Theorem 1 (i)–(iv) and the whole quantitative chain (L6, L8–L17, all intermediate bounds of L16) on cycles with k in {1000, 3000}, mu in {1, 0.3, 0.02}; **own hybrid integrator** (numerical spike, closed-form trough for the state and for the Riccati pair) | 6 cycles, 252 checks, 0 failures |
| `v_segments.py 400` | the *derivation steps* on non-periodic data: 400 negative lobes (L8, L9, L10), 400 spikes (L12, L13, L16 pieces), and Lemma 17 as a statement about the linear cooperative system for 133 **arbitrary** periodic x(t) (spectral radius of the monodromy vs the bound) | 15 870 inequality tests, 0 failures |
| `v_stage.py 1500` | the Stage-1/Stage-2 implications (L11, L14) on 1500 non-periodic lobes with k >= 1000 | 937 lobes satisfy the L11 hypotheses (min p_k/k = 1.013 >= 0.374), 613 the L14 hypotheses (min r_l/T_- = 0.488 >= 0.3); max tau_b = 0.0069 < 0.36; 0 failures |
| `v_stage_adv.py` | adversarial search: differential evolution minimising p_k/k over lobes with k >= 1000, T_- >= 5.89 | minimum 0.916 (at T_- = 5.890, mu = 1), with or without the Lambda_E hypothesis; claim is 0.374 |
| `v_recip.py` | Cor. 3: `T' = mu T` and `log m_t` equal at (mu,k) and (1/mu, k/mu^3) — 4 large-gain pairs (hybrid integrator) and 3 moderate pairs (monodromy); plus: the orbit on W^u(E) and 16 widely spread section data land on the same cycle at 4 parameter sets | all agree to >= 8 digits; ok |
| `v_rho40_hybrid.out` | re-computation of the two flagged rho = 40 cycles with the hybrid integrator | ok |

Outputs: `verify/*.out`.

Key numbers (independent integrator; they reproduce the lane's table):

| mu | k | lambda | T mu/lambda | p_k/lambda | lambda tau_c (bound) | <eta>/mu (bound 0.83) | log m_t |
|---|---|---|---|---|---|---|---|
| 1 | 1000 | 1000 | 1.30699 | 0.999999 | 27.9 (87.7) | 1.70e-2 | -1284.77 |
| 1 | 3000 | 3000 | 1.30492 | 1.000000 | 31.7 (104.1) | 6.5e-3 | -3889.26 |
| 0.3 | 1000 | 3333 | 1.30450 | 1.000000 | 33.1 (105.7) | 4.1e-3 | -4330.51 |
| 0.3 | 3000 | 10000 | 1.30409 | 1.000000 | 36.8 (122.2) | 1.5e-3 | -13020.9 |
| 0.02 | 1000 | 50000 | 1.30392 | 1.000000 | 45.1 (146.3) | 3.5e-4 | -65173.3 |
| 0.02 | 3000 | 150000 | 1.30390 | 1.000000 | 48.7 (162.8) | 1.3e-4 | -195560 |

Tightest places found (all on the correct side):

- Lemma 6, `T_+ omega / pi`: 0.984 at (mu, rho) = (1, 1.0005), 0.982 at (0.3, 1.0005), 0.933 at rho = 1.01, decreasing in rho. The inequality is sharp at Hopf, as it must be.
- Lemma 9, `P_w / ((gamma_E - mu) v_w)` (the single place where `Lambda_E > 0` enters the negative lobe): up to 0.986 near Hopf.
- Lemma 17 on arbitrary periodic x: smallest relative slack 0.25; on cycles the smallest absolute slack is 8.0.
- Lemma 10(c) on segments: relative slack down to 4e-8 (the bound is attained in the limit).
- On all 82 moderate cycles `<eta> - mu < 0` (largest value -2.3e-5, next to Hopf).

## 2. Line-by-line check of the proof

Every lemma was re-derived by hand before any numerics.

**Imports.** Re-derived rather than trusted:

- (F1) Used only in Lemma 17. `(a,b,c) = e^{-At}(H_3,H_2,H_1)` solves the second additive compound of the variational system (checked: `J^{[2]}` with one sign flip is `C`). Its monodromy is entrywise positive, the Perron root is the spectral radius, and the compound multipliers are the pairwise products of `{1, m_s, m_t}`. So the bound of Lemma 17 controls both transverse multipliers even without the identification `m_t = exp((<eta>-mu)T)`; with it (v2.0 Prop. 2.1, checked) one gets the statement on `<eta>`.
- (F2) `Delta = zeta_E U - eta_E p - w` satisfies `Delta' + gamma_E Delta = -k(e^U - 1 - U)` (checked by direct differentiation), hence `Lambda_E = -Delta > 0` on any periodic solution. `C_+` is excluded by `Lambda_E > 0`; `C_-` is excluded independently because it is backward invariant (the backward system in `(-U, p, w)` is cooperative with the origin fixed) and U is monotone there.
- (F3) v2.0 Thm 4.1 read in full: injectivity from Thm 3.1 in both directions, regularity from the orthogonality identity, winding-number count. No gap found. Numerically: exactly one trough and two zero crossings of U per period on all 82 moderate-gain cycles.
- (F4) Wirtinger plus `<p'^2> = mu <p^2>`; strictness argument correct.
- (F5) Ratio argument for `Z/X`, `Y/X` re-derived: `X <= M_0`.
- (F6) Direct averaging.
- (F7) **Checked against the source**, arXiv:1109.5204v1, pages 2–3. Theorem 1.1(c) reads: "If k > k3 + k5, then every solution not starting on S or on the one-dimensional stable manifold of E, converges to a nontrivial periodic orbit." `S = {x = 0}` is outside the open domain D, and the statement is for all positive `k3, k5`, i.e. every mu > 0. The quotation in v2.0 and in `theorem.tex` is accurate. (Prop. 1.2 there also says `W^s(E)` in the octant is E plus two monotone orbits.)

**Lemmas.** For each: the differential identities, the sign of every dropped term, every kernel comparison (`e^{-mu s} >= e^{-s}` is used only against non-negative integrands), and every constant.

- L5: both identities correct; `x - k_H` changes sign because p has only two zeros.
- L6: `U'' - 2 sigma U' + zeta_E U = -Lambda_E` correct; the Sturm integral has the right boundary term `omega z(pi/omega) >= 0`; `4 omega^2 = 3 gamma_E^2 - 2A gamma_E - (1-mu)^2` correct; `T_0(1000) = 5.892767 > 5.89`.
- L7: correct; (d), (e) are exact integrating-factor identities.
- L8, L9, L10: correct. The backward formula for `|p|` in L9 and the forward one in L10 are the right ones; L10(d) contradiction is valid because `s_l + a/2 < T_-`.
- L11: `T_- = s_w + u + r_l + tau_b`, thresholds sum to exactly 5.89, margins 0.0028 (window), 2.3690 vs 2.3130 (strip), 1.4863 vs 1.48 (both cases), 0.37453 vs 0.374. All recomputed.
- L12: the contradiction is evaluated at `s_bar` only when `s_bar <= s_+`; `-w <= |w_k|` on `[0, s_+]` follows from L7(d).
- L13: energy identity, `w <= m`, change of variables, and the quadratic step `M_+ >= sqrt(Q) - mu P_1/2` correct; the case `E_theta < 1` is harmless.
- L14: both cases correct; the minimax `(0.89 T' - 1)/1.89` and the final `0.1708 T - 0.9416 log T - 0.6987 > 0` (0.148 at T = 22, increasing) recomputed.
- L15: no continuity or continuation argument is actually used (the notes call it a "continuum argument"; the proof is a direct exclusion): Step 1 gives `z >= z_1(k)` unconditionally, `h` is quasi-convex (`sign h' = sign(mu z^2 - z - 6.3)`, re-derived), both endpoint values are below `0.3 kappa_0 k` for all k >= 1000 and mu <= 1, so `z > z#`. Not circular.
- L16: every sub-step re-derived; the monotonicity paragraph (the repair after the cross-model review) is correct: each side condition is monotone in the right direction, and the 601-point grid agrees.
- L17: `L' = -nu L + (beta x - mu + nu) a`; Dini bound for the max-norm; matching constants `1/alpha < 1`, `1/beta < 2 <= s^2`. Correct for all `0 < nu < mu <= 1`, `s >= sqrt 2`.
- Theorem: arithmetic re-done, `q(1000) = 0.077177 < 0.0772`.
- Cor. 2: basins open and disjoint, union is `D \ W^s(E)` by (F7), and the measure-zero cone argument for connectedness is valid.
- Cor. 3: `sigma`, `rho` invariance and `t' = mu t` re-derived; every nonconstant periodic solution of the scalar equation lies in the positive domain (so "cycle" means the same on both sides); only cycles and multipliers are transported, and uniqueness is then argued at `(mu,k)` itself. The earlier false sentence about conjugate positive domains is gone.

**Circularity / silent weakening.** None found. Stability is never assumed: the chain is F2–F6 → L6 → L11 → L13 → L15 → L16 → L17. The statement is openly restricted to `k >= 1000` (not uniform in mu); the notes and the remark say so.

## 3. Things that looked like failures and were not

Reported so that nobody mistakes them later. All were defects of my checkers.

- `v_consts.py`, first run: `H <= H_hat` "failed" by 2e-50 at mu = 1, k = 1000 — it is an equality case; tolerance added.
- `v_moderate.py`: `y > 0` and `p < lambda` "failed" by 1e-9 k and 2e-10..4e-9 in deep troughs. Both are exact for periodic solutions (`y = int e^{-(t-s)} x ds`, `(lambda-p)' + mu(lambda-p) = y`), with true margins `~ e^{-T}`, far below double precision; the explicit integrator and its dense interpolant produce noise of that size. The two remaining flags in the scan (rho = 40, mu = 1 and 0.6) were recomputed with the closed-form trough: `lambda - p >= 0` there (`v_rho40_hybrid.out`), and the periods agree with whlib to 6 digits.
- `v_largegain.py`: the post-crest inequality "failed" by 1e-6 because `1-(1+v)e^{-v}` was evaluated with cancellation against `W_c ~ 1e11`; with `expm1` the relative slack is >= 5e-4 everywhere.
- `v_segments.py`: a NaN guard was added so that a non-finite value can never count as a pass.

## 4. Remarks that do not affect the verdicts

1. Theorem (iii): with the stated time origin the set `{t in [0,T): x > c_0}` wraps around t = 0; it is one interval modulo T.
2. "Psi_mu is increasing in u" is used in L14 and asserted without proof; one line suffices: `Psi' = (1-e^{-u}) - mu Psi >= (1-e^{-u}) e^{-mu u}`.
3. The numerical tests cannot probe the quantifier "every cycle" (only one cycle exists numerically); the nearest surrogates are the non-periodic lobe tests and the adversarial search of `v_stage_adv.py`. The verdicts rest on the hand check.
4. Nothing here is an interval-arithmetic certificate; none is needed, because the proof's only numerical content is the list of scalar inequalities recomputed in `v_consts.py`.
5. Imports from v2.0 were checked by me for this report, but v2.0 itself is unrefereed; the theorem inherits that status for (F1)–(F6).

## 5. Verdicts

| # | statement | verdict |
|---|---|---|
| 1 | Theorem 1 (i)–(iv), 0 < mu <= 1, k >= 1000 | holds |
| 2 | `T_+ < pi/omega`; `T_- > 2pi - pi/omega_0(k)`; `p_k >= 0.374 k` for k >= 1000 | holds |
| 3 | spike mass lemma and bootstrap to `T_- > 0.1833 lambda/mu`, `p_k >= 0.819 lambda` | holds |
| 4 | linear estimate (Lemma 17) | holds |
| 5 | Corollaries 2 and 3 (uniqueness, global attraction off `W^s(E)`, reciprocal range) | holds, conditional on Smith 2012 Thm 1.1(c), which was read in the source and matches |
