# Lane `largegain` — notes (2026-10-03)

## 0. Outcome in one paragraph

**Theorem (proved, full proof in `theorem.tex` / `theorem.pdf`, 8 pages).**
For `0 < mu <= 1` and `k >= K0 = 1000` (i.e. `rho >= 1000/(mu(1+mu))`), *every* positive cycle satisfies

- (i) `T >= 0.1833 * lambda/mu`,
- (ii) `p >= 0.819 * lambda` at the up-crossing of `U = 0`,
- (iii) `|{x > c0}| <= (71.79 + 12.28 log(lambda/1000)) / (0.819 lambda)`, `c0 = (mu/2)(1 - mu^2/4)`,
- (iv) `<eta> <= 0.83 mu < mu`, so `0 < m_s < m_t < 1`.

Corollary (uses Smith's convergence theorem as quoted in v2.0 + the basin-partition argument of brief 3(a)):
the cycle is unique and attracts every positive orbit off `W^s(E)`.
Reciprocal symmetry `(mu,k) -> (1/mu, k/mu^3)` (same scalar equation in frequency-one time, `rho` invariant,
`t' = mu t`): for `mu >= 1` the same holds when `k >= 1000 mu^3`, with `<eta> <= mu - 0.17`.
In terms of `rho`: `rho >= rho_0(mu) = 1000/(m(1+m))`, `m = min(mu, 1/mu)`; `rho_0(1) = 500`.

**Uniformity as mu -> 0 is NOT achieved** (`rho_0 ~ 1000/mu`). The threshold is poor (numerically every
computed cycle is attracting for all `rho > 1`), but it is explicit and needs no continuation.

Status of the proof: every step written out; imports (F1)–(F6) from the v2.0 manuscript and (F7) = Smith 2012
for the uniqueness corollary only. Not refereed. A cross-model adversarial review was requested (section 6).

## 1. Files (all in this directory)

| file | purpose |
|---|---|
| `theorem.tex`, `theorem.pdf` | the theorem and its complete proof (17 numbered statements) |
| `constants.py` | recomputes **every** numerical constant / scalar inequality used in the proof (mpmath, 40 digits); exit 1 on failure |
| `check_lemmas.py` (+ `check_lemmas.out`) | tests every lemma inequality, and the theorem's conclusions, on computed cycles |
| `test_segments.py` | tests the *derivation steps* of the lemmas on random **non-periodic** orbit segments (covers the "medium" regimes no real cycle visits) |
| `lgcycle.py` | large-gain cycle finder (forward integration, the cycle is super-attracting) |
| `scale1.py` | first scaling survey (uses `whlib`) |
| `explore_K0.py` | design aid: which `K0` the chain supports |
| `test_stage.py` | direct test of the Stage-1/Stage-2 lemmas (L11, L14) on random non-periodic negative lobes with `k >= 1000` |
| `sol_review.sh`, `sol_review.md` | cross-model adversarial review request / reply (reply refers to the pre-repair text) |

Commands (python = `/Users/admin/.venvs/analysis-py313/bin/python`, run inside this directory):

    python constants.py            # ALL CONSTANT CHECKS PASSED for K0 = 1000   (<1 s)
    python test_segments.py 400    # ALL SEGMENT CHECKS PASSED                  (~5 s)
    python check_lemmas.py         # ALL LEMMA CHECKS PASSED                    (~4 min, 11 cycles)
    python check_lemmas.py 1.0 1000.0   # single case
    python test_stage.py 500       # ALL STAGE CHECKS PASSED (69 lobes with T_- >= 5.89: min p_k/k = 0.917 >= 0.374;
                                   #   14 lobes with T_- >= 22: min r_l/T_- = 0.678 >= 0.3; max strip time 0.098 < 0.36)
    pdflatex theorem.tex           # 9 pages

## 2. What the numerics say (floating point)

Computed cycle, `lambda = k/mu` (from `check_lemmas.out`, `scale1.py`):

| mu | k | lambda | T | T mu/lambda | T_+ (U>0) | p_k/lambda | D = -U_min | h = U_max | log m_t |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 16 | 16 | 24.11 | 1.507 | 0.568 | 0.9965 | 285 | 4.81 | -13.98 |
| 1 | 60 | 60 | 81.37 | 1.356 | 0.203 | 0.9997 | 4.43e3 | 7.36 | -67.5 |
| 1 | 300 | 300 | 394.3 | 1.314 | 0.053 | 1.0000 | 1.16e5 | 10.55 | -375.7 |
| 1 | 1000 | 1000 | 1307.0 | 1.307 | 0.0189 | 1.0000 | 1.30e6 | 12.95 | -1284.8 |
| 1 | 3000 | 3000 | 3914.8 | 1.305 | 0.0071 | 1.0000 | 1.17e7 | 15.15 | -3889.3 |
| 0.5 | 1000 | 2000 | 5220.2 | 1.305 | 0.0106 | 1.0000 | 1.04e7 | 15.03 | — |
| 0.3 | 1000 | 3333 | 14494 | 1.304 | 0.0069 | 1.0000 | 4.82e7 | 16.56 | — |
| 0.1 | 100 | 1000 | 13056 | 1.306 | 0.0215 | 1.0000 | 1.30e7 | 15.25 | -1290.4 |
| 0.02 | 20 | 1000 | 65274 | 1.306 | 0.0234 | 1.0000 | 6.49e7 | 16.86 | -1290.4 |

Empirical laws: `T ≈ 1.305 lambda/mu`, spike mass `M_+ ≈ 1.305 p_k^2`, `p_k ≈ lambda`, `∫_crest x^{1/3} dt ≈ 4.6`
(universal spike), `∫ eta dt` over a period grows only like `log lambda` (≈ 3.5 on the crest + ≈ `log`-growing
trough transient), while `mu T ≈ 1.3 lambda`. So `log m_t ≈ -(mu) T (1 - O(log(lambda)/lambda))`: the true
margin is enormous; the difficulty is purely the *a-priori* control of **every** cycle.

## 3. Structure of the proof (statement → status)

Notation: `w = y - k`; time origin at the down-crossing of `U=0`; negative lobe `(0,T_-)`, positive lobe
`(T_-,T)`; `s_w` = zero of `w`, `s_0` = trough, `s_l` = rising passage of `U=-l` (`l=2`), `p_k = p(T_-)`,
`M_+ = ∫_{U>0}(x-k) = ∫_{U<0}(k-x)`.

| # | statement | status | numerical test |
|---|---|---|---|
| L5 | `M_+ = w_d + P_d + |w_k| + p_k < k T_-`; `<(x-k_H)p^2> = 0` so `x_min < k_H` | proved | cycles ✓ |
| L6 | **`T_+ < pi/omega`** (Sturm, uses `Lambda_E > 0`); `omega >= ½ sqrt(3k^{2/3}-4k^{1/3}-1)`; hence `T_- > 2pi - pi/omega_0(k)` | proved | cycles ✓ |
| L7 | sign structure of `w, p, phi` on both lobes (cone exclusion) | proved | cycles ✓ |
| L8 | phase (i): `s_w >= a` ⇒ `|U(a/2)| >= k c_i(a)(1-e^{-|U(a/2)|})` ⇒ `D >= 0.7968 k c_i(a)` | proved | cycles ✓, segments ✓ |
| L9 | phase (ii): `P_w >= k phi(s_w)(u-1+e^{-u})`, `D - v_w >= k phi(s_w) c_ii(u)`, `P_w < (gamma_E - mu) v_w` ⇒ `D >= 0.7968 k c_ii(a)` | proved | cycles ✓, segments ✓ |
| L10 | rising phase: `p' < k` (⇒ `D >= p_k^2/2k`, flank time `<= 2L/p_k`); `p(s_l) >= k phi_l Psi_mu(r_l)`; strip time `tau_b < a` when `k c_b(a) >= l/phi_l` | proved | cycles ✓, segments ✓ |
| L11 | **Stage 1**: `k >= 1000` ⇒ `p_k >= 0.374 k` (pigeonhole on `T_- = s_w + u + r_l + tau_b > 5.89`) | proved | cycles ✓ |
| L12 | `P_1 <= p_k (1 + 1.1 sqrt2 k p_k^{-3/2})` | proved | cycles ✓, segments ✓ |
| L13 | **spike mass**: `M_+(M_+ + mu P_1) >= ¼ p_k^2[½ p_k^2 - 2 mu P_1 H] - p_k k/sqrt2`, so `T_- > M_+/k >= p_k^2 sqrt(1-delta)/(sqrt8 k) - ½` | proved | cycles ✓, segments ✓ |
| L14 | **Stage 2**: `T_- >= 22` ⇒ `r_l >= 0.3 T_-` ⇒ `p_k >= 0.95 phi_l k (1-e^{-mu z})/mu - 2`, `z = 0.3T_- - 3` | proved | cycles ✓ |
| L15 | **bootstrap** (continuum argument, quasi-convex `h(z)`): `T_- > 0.1833 lambda/mu`, `p_k >= 0.819 lambda` | proved | cycles ✓ |
| L16 | crest duration `tau_c <= (71.79 + 12.28 log(lambda/1000))/(0.819 lambda)` | proved | cycles ✓, segments ✓ |
| L17 | linear estimate: `(<eta>-mu)T <= -nu(T-tau_c) + (s-mu)tau_c + kT/s^2 + log(1+alpha s+beta s^2)` via `L = a+alpha b+beta c` (trough) and `N = max(s^2 a, s b, c)` (crest) on the cooperative system `a'=-mu a+b, b'=-b+c, c'=-Ac+xa` | proved | cycles ✓ |
| Thm 1 | (i)–(iv) above | proved | cycles ✓ (k = 1000, 3000; mu = 1, .5, .3) |
| Cor 2 | uniqueness + global attraction off `W^s(E)` | proved **modulo (F7)** (Smith) | — |
| Cor 3 | `mu >= 1` by reciprocal symmetry | proved | — |

"cycles ✓": inequality checked on the 11 computed cycles of `check_lemmas.py` (all pass). "segments ✓":
derivation step checked on 365 random non-periodic segments (all pass, 0 failures in ~18 000 inequality tests).
These tests can only catch gross errors; the lemmas concern hypothetical cycles that do not exist.

### The mechanism in words

1. `Lambda_E > 0` makes the positive lobe shorter than half a linear period, `T_+ < pi/omega = O(k^{-1/3})`,
   while Wirtinger gives `T > 2pi/sqrt(mu)`. So **every** cycle has a long negative lobe (`T_- > 5.89` at `k>=1000`).
2. A long negative lobe is impossible for a "small/medium" cycle: split `T_-` into four durations; each long
   duration forces either a deep trough (`D >~ 0.8 k`) or a steep rising flank; the only escape (hovering in the
   strip `-2 < U < 0`) is excluded by the cooperative structure of the rising phase and by `Lambda_E > 0`.
   Result: `p_k >= 0.374 k`.
3. The spike converts slope into mass quadratically (`M_+ >~ p_k^2/sqrt8`), and mass into time (`T_- > M_+/k`).
   Feeding this back (`T_-` long ⇒ rising phase long ⇒ `p_k` close to `lambda`) and closing with a
   quasi-convexity argument gives `T >= 0.18 lambda/mu`, `p_k >= 0.82 lambda`: every cycle is a relaxation cycle.
4. Then the crest lasts `O(log(lambda)/lambda)`; the adjoint system can grow there by at most
   `exp(2 sqrt(lambda) tau_c + mu T/4)` times a polynomial factor, and decays like `exp(-(mu/2)·)` on the trough.

## 4. Sharpest obstructions / limitations

- **Non-uniformity in mu.** The chain needs the fast scale `k^{-1/3}` to be short compared with the window
  `T_0 ≈ 2pi` supplied by `T > 2pi/sqrt(mu)`. Two steps need `k` (not `rho`) large: the strip estimate
  `k c_b(a_4) >= l/phi_l` and, more seriously, phase (ii), where the only information that keeps the cycle away
  from `E` is `Lambda_E > 0`, i.e. `P_w < (gamma_E - mu) v_w`; for small `k` the factor `gamma_E - mu ≈ 1 + k`
  no longer beats `k(a-1+e^{-a})` unless `a ~ mu^{-1/2}` and `rho sqrt(mu) >~ 1`. Re-doing the estimates with
  windows `∝ mu^{-1/2}` would plausibly give `rho_0(mu) = O(mu^{-1/2})`; a `mu`-uniform `rho_0` needs a new idea
  for phase (ii) (e.g. using `K > 0`/the cycle's own adjoint instead of `Lambda_E`). **Not done, not claimed.**
- **The threshold `K0 = 1000` is limited by Stage 1** (window sum `5.89 <= T_0(1000) = 5.8928`) and by the
  endpoint check `h(z^sharp) <= 0.3 kappa_0 k`. With re-tuned windows the same proof closes for `K0 ≈ 700`
  (`explore_K0.py`), not lower. The final linear estimate has huge slack (margin 31.7 at `K0`).
- Losses: `r_l >= 0.3 T_-` (truth ≈ 0.99), `M_+ >= 0.354 p_k^2` (truth 1.305), `phi_l = 0.865`. A sharper
  spike lemma (the exact scale-invariant spike `V''' = -k e^V` has `M = c_* p_k^2`, `c_* ≈ 1.305`) and a
  computer-assisted treatment of the O(1) problem `U''' + e^U - 1 = 0` in scaled variables would be needed to
  bring `rho_0` down to O(10).
- The moderate-`rho` regime (`1 < rho < rho_0`) is untouched by this lane; near Hopf the mechanism is absent.

## 5. Approaches considered and discarded (for the record)

- `<eta> <= sqrt<zeta>`, `f(<eta>) < k` (Jensen), `L' <= (beta x - nu)L`: all linear-in-`∫x` bounds are useless
  (`∫x = kT`); the spike must be handled by a norm adapted to `x^{1/3}`.
- Area identity `∬_Omega (rho e^U - 1) dU dp = 0` (equivalent to `<x p^2> = k_H <p^2>`): gives
  `∫_{x<k_H} p^2 >= (rho-1) ∫_{U>0} p^2` but no amplitude lower bound by itself.
- Virial-type identities found on the way (true on every cycle, not needed in the end):
  `<w^2> = k_H <p^2>`, `<(x-k)^2> = <w^2> + <w'^2>`, `A<p^2> = <xU> - k<U>`.
- Exact formula `U(t) = -∫_0^∞ e^{-sigma r} sin(omega r)/omega * Lambda_E(t+r) dr` (anti-causal Green function of
  the unstable pair): only controls the fast scale.

## 6. Cross-model adversarial review

`sol_review.sh` sent `theorem.tex` to `gpt-5.6-sol` with a refute-it brief (reply in `sol_review.md`; treated
as hypotheses and each checked by hand). Outcome:

- Lemmas 5–14 (identities, short lobe, signs, phases (i)/(ii), rising phase, Stage 1, `P_1`, spike mass,
  Stage 2) and L17 (linear estimate): **"no error found"**, each re-derived by the reviewer.
- GAP (L15, Step 3): non-emptiness of `[z_1(k), z^sharp]` not stated — **fixed** (`z_1(k) < 0.015k < 0.055k <= z^sharp`;
  checked in `constants.py`).
- GAP (L16): the blanket sentence "side conditions need only be checked at `lambda = K0`" was not a proof —
  **fixed**: explicit monotonicity paragraph (each majorant is a product of positive decreasing functions of
  `lambda`), extra checks in `constants.py`.
- GAP (L16): kernel comparison `e^{-mu(v-r)} >= e^{-(v-r)}` after the crest only valid where the majorant of `w`
  is nonnegative — **fixed** (made explicit; it was already restricted to `v <= v_+`).
- GAP (Cor 2): connectedness of `D \ W^s(E)` used unproved dimension theory — **fixed**: elementary argument
  (cones over countably many compact `C^1` arcs have measure zero, so two-segment paths avoid `W^s(E)`).
- **FATAL for one sentence of Cor 3 (confirmed by hand):** I had written that the time change conjugates the two
  flows "on their positive domains". False: in scalar coordinates `D = {p<lambda, p' < mu(lambda-p)}` while
  `D' = {p<lambda, p' < lambda-p}`. Only periodic orbits and their multipliers are transported. **Fixed**:
  transport cycles + multipliers, then run the Cor 2 argument at `(mu,k)` itself (Smith's theorem holds for all
  `mu > 0`). The main theorem is unaffected.

No web searches were used by this lane.

## 7. What I would try next

1. `mu`-uniform version: redo Stage 1 with windows scaled by `mu^{-1/2}`; replace the `Lambda_E` step in phase (ii).
2. Lower `K0`: sharpen `r_l >= 0.3T_-` (true ≈ 1) by bounding `s_w` and `u` logarithmically
   (`s_w <= log(1 + w_d/(k phi))`, `u <= (1/mu) log(1 + 2 mu T_-)`), and prove `M_+ >= c p_k^2` with `c` near 1.3
   from the scale-invariant spike.
3. For moderate `rho`: the same four-duration decomposition gives `p_k >~ k^{2/3}`-type footholds for any `k`
   with `T_0(k) > 0`; the missing piece is a non-asymptotic version of step 3 (computer-assisted in the scaled
   variables `(U, p k^{-1/3}, w k^{-2/3})`).
