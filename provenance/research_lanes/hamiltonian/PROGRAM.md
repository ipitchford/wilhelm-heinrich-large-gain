# PROGRAM — a complete rigorous route to the Wilhelm–Heinrich uniqueness conjecture (design, not executed)

Lane `hamiltonian`, 2026-10-03. Nothing in this file is a proof. Numbers are floating point unless marked;
the only interval computations are the Taylor-model probes of Section 4.7. Proved inputs are cited by location.
Notation: `mu ∈ (0,1]`, `rho = k/k_H > 1`, `A = 1+mu`, `lambda = k/mu = rho A`, `eps = sqrt(mu)/(1+mu) ∈ (0,1/2]`,
`r = (rho−1)/eps^2`, `tau = sqrt(mu) t`.

## 0. Verdict

- The architecture of the proof of Jones' conjecture (asymptotic region + near-bifurcation region + computer-assisted
  compact core, glued by "every cycle is attracting ⇒ unique") transfers. The ODE core is much easier than the
  delay-equation core was.
- **Compact core: feasible.** Estimated `10^5–10^8` rigorous return-map evaluations for
  `mu ∈ [0.05,1]`, `rho ≤ 20`, outside a Hopf collar `r ≤ 0.2` (ideal count `7·10^3`, times an overestimation factor
  `10^2–10^4` that I cannot pin down without a prototype). At an assumed 10 ms per evaluation (CAPD-class C++; not
  measured here) that is 2 CPU-hours to 8 CPU-days. At the **measured** 1–10 s per evaluation of TaylorModels.jl it is
  feasible only at the optimistic end.
- **The new difficulty relative to Jones is the second parameter.** The parameter set is not compact at `mu → 0`, and
  there the problem has *two* singular regimes with a matching zone (NOTES §3–5): a Toda/Hamiltonian wedge
  `rho − 1 = O(mu)` selected at third order, and a relaxation regime for `rho − 1 = O(1)`. No amount of computation at
  `mu ≥ mu0` removes the need for an analytic theorem on `0 < mu < mu0`, all `rho > 1`.
- **What exists today:** the local Toda-wedge theorem (NOTES Thm A, Thm C). **What does not exist:** an exclusion
  of large non-attracting cycles for small `mu` (NOTES §4.2b reduces it to two a-priori estimates); any theorem in the relaxation regime; an explicit large-gain threshold;
  an explicit Hopf collar. Each is a well-posed project; the small-`mu` wedge is the largest and the riskiest.
- Without the wedge, the program would still deliver "uniqueness for `mu0 ≤ mu ≤ 1/mu0`, all `rho > 1`", provided the
  large-gain and collar pieces are done. That is the realistic first theorem.

## 1. Logical skeleton

- (L1) *Symmetry.* In `tau`-time the scalar equation depends on `(sigma, rho)` only and `sigma(mu) = sigma(1/mu)`;
  cycles correspond bijectively (`p = U'`, `y = k − mu p − p'`). So `mu ∈ (0,1]` suffices for the cycle count.
- (L2) *Pointwise criterion* (brief 3(a); rests on Smith 2012: every positive orbit off `W^s(E)` converges to a
  cycle). If at `(mu,rho)` every cycle has `m_t < 1`, the cycle is unique. No continuation and no existence proof
  is needed, and local uniqueness inside a box never has to be shown.
- (L3) *Where the cycles are* (v2.0 Thm 4.1, 4.2, eq. (4.1), (6.1), Cor. 3.2, Lemma 2.3). Every cycle crosses the
  crest section `Sigma = {p = 0, y > k}` exactly once per least period, so cycles ↔ fixed points of the **first**
  return map `Pm` on `Sigma`. In section coordinates `(U,y)` every crest lies in the explicit region

      R(mu,rho) = { 0 < U ≤ L,   max( zeta_E U,  k(e^U − 1)/(1+eta_+) )  <  y − k  <  k(e^U − 1) },
      L = log M0,  M0 = (lambda+mu)(lambda+1)/mu,   eta_+ = (mu−1+sqrt((mu−1)^2+4 k M0))/2.

  (`y−k > zeta_E U` is `Lambda_E > 0`; the other bounds are `x − k = (1+eta)(y−k)` with `0 < eta ≤ eta_+`.)
  The line `y − k = zeta_E U` is the trace of the unstable eigenplane of `E`.
- (L4) *The multiplier test.* For every `z ∈ R` whose return lies in `Sigma`, `DPm(z)` is entrywise positive (the
  wedge argument of v2.0 Prop. 4.4 uses only `x > y > k` at the start, `y > k` at the return, and the cooperative
  compound system `w' = C w`). `det DPm = e^{−A T} < e^{−4 pi}` at cycles. Hence at a crest `z`:
  `m_t < 1 ⟺ exists v > 0 with DPm(z) v < v` (Collatz–Wielandt).
- (L5) *Return times are bounded on `R` away from `E`:* `W^s(E)` never meets `{p = 0}` except at `E` (both branches
  stay in the backward-invariant octants `{U>0,p<0,y<k}` and `{U<0,p>0,y>k}`), so `R minus B_delta(E)` is at positive
  distance from it. Sign check: on the faces `U = 0`, `p = 0`, `y = k` of each octant the backward field points
  inward (`−p`, `−(k−y)`, `−k(e^U−1)` have the required signs), and the stable eigenvector
  `(1, −gamma_E, −k/(gamma_E−1))` points into the first octant. Numerical test on 5 parameter pairs:
  `ws_check.py` → `ws_check.out`, all 10 branches stay in their octant. Since every orbit off `W^s(E)` converges to
  a cycle and cycles cross `Sigma` transversally, first-return times are locally bounded on the compact set.

## 2. The parameter plane: scales, established in this lane

| regime | scaling | governing object | weak multiplier | source |
|---|---|---|---|---|
| Hopf end | `r → 0`, any `mu` | normal form, `a^2 ≈ 2r(1+4 eps^2)` | `log m_t ≈ −2 pi eps^3 r` | v2.0 App. A |
| Toda wedge | `mu → 0`, `r = O(1)` | `r = Rsel(h) = (2/3)<p^4>/<p^2>` on `H0 = h` | `log m_t = −eps^3 J0 Rsel'(h*)`, `Rsel' ≥ 4/5` | NOTES Thm A, C (proved, local) |
| matching | `1 ≪ r ≪ 1/eps^2` | `h* ≈ (5/4) r` | `log m_t ≈ −2.108 (rho−1)^{3/2}` | formal + numerics |
| relaxation | `mu → 0`, `rho` fixed | 1-D map `F_lambda = G_lambda∘S`, branch `rho = Lambda(q)` | `m_t → F'(q*) ∈ (0,1)`: 0.947 (`rho`=1.1), 0.680 (1.5), 0.187 (3), 0.0015 (8) | NOTES §5 (numerical) |
| large gain | `rho → ∞` | pulse + deep trough | `log m_t ≈ −2.5 rho` at `mu = 1` (−505 at `rho = 200`) | cartography grid |

Consequences for the design. (1) The natural near-Hopf variable is `r`, and in it the collar is **uniform** in
`mu ∈ (0,1]`; in `rho` its width is `≈ r_c·mu` and cannot be uniform. (2) At fixed `rho > 1` the weak contraction does
**not** degenerate as `mu → 0` (it tends to the relaxation multiplier), so a core computation at small `mu` is limited
by the period `T ~ Sigma(rho)/mu`, not by loss of hyperbolicity. (3) The only place where `1 − m_t → 0` is the
Hopf curve `r → 0`.

## 3. The partition

With placeholders `mu0 = 0.05`, `r_c = 0.2`, `rho1 = 20` (to be fixed by the collar and large-gain work):

    H  (Hopf collar)   mu ∈ [mu0,1],  0 < r ≤ r_c                 [rho − 1 ≤ r_c mu/(1+mu)^2; 0.0091 … 0.05]
    C  (core)          mu ∈ [mu0,1],  r ≥ r_c,  rho ≤ rho1
    L  (large gain)    mu ∈ [mu0,1],  rho ≥ rho1(mu)
    W  (small-mu wedge) 0 < mu < mu0, all rho > 1, split as
       W_H: r ≤ r_c     W_A: r_c ≤ r ≤ r2 (Toda)     W_M: r ≥ r2, rho − 1 ≤ delta (matching)
       W_B: delta ≤ rho − 1, rho ≤ rho1 (relaxation)   W_L: rho ≥ rho1

Claim to prove in every piece: *every cycle is hyperbolic attracting.* By (L2) the pieces need not be glued.

## 4. The core: a branch-and-bound over all crests

### 4.1 Certificates
For a parameter box `Q = [mu]×[rho]` and a state box `b ⊂ Sigma`, either of the following disposes of `b`:

- **(E) exclusion.** A rigorous enclosure of `Pm(b)` (for all parameters in `Q`) is disjoint from `b`. Then `b`
  contains no crest. Needs only a `C^0` tube enclosure and a crossing-time enclosure.
- **(C) contraction.** There is `v > 0` with `Mbar v < v`, where `Mbar ≥ DPm(z)` entrywise for all `z ∈ b`, all
  parameters in `Q`. Then every crest in `b` has `m_t < 1`. `Mbar` can be obtained **without** a `C^1` integrator:
  `DPm` comes from the compound system `w' = C(t) w`, `C = [[−mu,1,0],[0,−1,1],[x(t),0,−A]]`, which is cooperative and
  depends on the orbit only through `x(t)`; an upper envelope `xbar(t) ≥ x(t; z)` over the tube gives an upper
  solution `wbar ≥ w`, and `Mbar` follows by dividing by the lower bound of `y − k` at the return (take the sup of
  `wbar` over the crossing-time window).

If `R(mu,rho) ⊂ ⋃ b_i` for every `(mu,rho) ∈ Q` and each `b_i` is disposed of, then every cycle with parameters in
`Q` is hyperbolic attracting, hence unique. Existence is not checked (Smith 2012 gives it), and two crests in one
(C)-box are not excluded directly — (L2) excludes them.

Near `E`: use the blow-up `(U, yh)` with `yh = (y−k)/U ∈ [zeta_E, k(e^U−1)/U]`. The blown-up return map
`(U,yh) ↦ (U·m(U,yh), yh')` extends continuously to `U = 0` (linear return map of `J_E`). Its only invariant
direction is `yh = zeta_E` (an invariant ray is an eigenvector of `e^{J_E T}`; the stable eigenvector is not in the
section), with factor `m_E = exp(2 pi Re/Im) > 1` (1.006 at `mu = 0.05, rho = 1.01` up to 9.7 at `mu = 1, rho = 20`;
column `m_E` of `core_cost.out`). So boxes `[0,d0]×[yh]` are excluded by (E) in the `yh`- or the `U`-coordinate; this replaces the
non-explicit constant `d_Q` of v2.0 Thm 6.1 Step 2 by a finite computation.

### 4.2 Why the proved structure makes this small
- One crest per cycle ⇒ fixed points of the *first* return only; no period-`n` search.
- `det DPm < e^{−4 pi}`: after one return every box is a sliver along the crest curve `Gamma` of the attractor; the
  search is effectively one-dimensional. Coarse boxes away from `Gamma` are excluded with margins of order their
  distance to `Gamma`.
- The contraction window is wide. With `pi` the reduced map along the weak direction, `delta* = (1−m_t)/|pi''|`:
  near Hopf `delta* ≈ U*/3`, and it grows to `> 1` for `rho ≥ 2` (column `delta*`). (C) only has to hold on
  `|U − U*| ≲ delta*`; outside, (E) has margin `≳ (1−m_t)|U − U*|`.
- Nested crests are ordered in both coordinates (v2.0 Cor. 4.3): once the attracting crest is enclosed, further
  crests could only lie in the NE/SW cones `zeta_i (U−U*) < y−y* < zeta_o (U−U*)`; this halves the region again.
- The quiescent phase can be integrated in closed form. For `U ≤ U_s` the term `x = k e^U ≤ k e^{U_s}` is a bounded
  perturbation of the linear system `p' = k − mu p − y`, `y' = −y`, `U' = p`; the exit time solves an explicit scalar
  equation. Troughs of depth `10^2–10^5` (`rho` up to 20, `mu` down to 0.05) therefore cost nothing and produce no
  wrapping.

### 4.3 Algorithm
1. Parameter boxes `Q` of size `(d_mu, d_rho)` (Section 4.5), `mu`, `k` carried as interval parameters.
2. Floating-point: locate the crest `z*`, the weak eigenvector, `delta*`.
3. (C) on a box around `z*` of half-length `≈ delta*/2` along the weak direction (a few sub-boxes if needed).
4. (E) on a geometric cover of the two remaining arcs of a tubular neighbourhood of `Gamma` (box length
   `≈ (1−m_t)·dist/2`), and on the blown-up boxes at `E`.
5. (E) on a coarse cover of the rest of `R` (adaptive bisection on failure).

### 4.4 Sizing data (floating point; `core_cost.py` → `core_cost.out`)
`kappa = pi''(U*)` by central differences along the weak eigenvector; `d_rho = delta*/(2|dU*/drho|)`,
`d_mu = delta*/(2|dU*/dmu|)`; `N_state = 2 ln(L/((1−m_t) delta*)) / ln(1+(1−m_t)/2)`.

| mu | rho | r | U* | T | 1−m_t | delta* | d_rho | d_mu | N_state |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 1.05 | 0.2 | 0.752 | 6.4 | 0.118 | 0.25 | 0.020 | (≈1, symmetric point) | 153 |
| 1 | 1.2 | 0.8 | 1.300 | 6.7 | 0.385 | 0.51 | 0.11 | – | 29 |
| 1 | 2 | 4 | 2.347 | 8.8 | 0.855 | 1.17 | 0.67 | – | 7 |
| 1 | 20 | 76 | 6.565 | 55 | 1.000 | 3.3 | 17 | – | 4 |
| 0.2 | 1.05 | 0.36 | 0.880 | 14.5 | 0.094 | 0.29 | 0.020 | 0.17 | 200 |
| 0.2 | 2 | 7.2 | 2.812 | 24 | 0.707 | 1.41 | 0.68 | 0.25 | 9 |
| 0.05 | 1.01 | 0.22 | 0.642 | 28.6 | 0.0125 | 0.21 | 0.0036 | 0.023 | 2335 |
| 0.05 | 1.1 | 2.2 | 1.684 | 34 | 0.106 | 0.57 | 0.043 | 0.027 | 162 |
| 0.05 | 2 | 22 | 3.821 | 78 | 0.604 | 1.34 | 0.57 | 0.039 | 14 |
| 0.05 | 20 | 419 | 8.271 | 580 | 1.000 | 4.1 | 21 | 0.11 | 4 |

(55 parameter points in the file; all Newton residuals `≤ 6e-11`.)

### 4.5 Idealised count (`core_total.py` → `core_total.out`)
Integrating `N_state/(d_rho d_mu)` over the core (log–log interpolation; the `mu = 1` row borrows `d_mu` from
`mu = 0.5` because `dU*/dmu = 0` there by symmetry):

| collar | `mu0` | ideal number of box evaluations | ×10^2 | ×10^4 |
|---|---|---|---|---|
| `r_c = 0.2` | 0.2 | 5.7e2 | 5.7e4 | 5.7e6 |
| `r_c = 0.2` | 0.1 | 1.9e3 | 1.9e5 | 1.9e7 |
| `r_c = 0.2` | 0.05 | 6.7e3 | 6.7e5 | 6.7e7 |
| `r_c = 0.05` | 0.05 | 4.0e4 | 4.0e6 | 4.0e8 |

The density behaves like `(rho−1)^{−2.3}` near Hopf and like `mu^{−2.8}` per unit `mu`, so the count is dominated by
the corner (`mu0`, collar edge). "Ideal" means: perfect enclosures, boxes as large as the margins allow, a
one-dimensional state search. The multiplier for a real implementation (two section dimensions, wrapping, the
compound upper bound in (C), safety) is unknown; `10^2–10^4` is my bracket, not a measurement.

### 4.6 What I did not size
The coarse cover of `R` far from `Gamma` (expected cheap, large margins), the blown-up neighbourhood of `E`
(a few hundred boxes per parameter box, `≈ ln(U*/d0)/ln m_E`), and the growth of the cover when `[mu]` is an
interval (the period changes by `dT/dmu·d_mu`, which the crossing-time enclosure must absorb).

### 4.7 Rigorous-integration probe (`tm_probe.jl`, `tm_probe2.jl` → `.out`)
TaylorModels.jl 0.12 `validated_integ` (Taylor models, order 2 in the initial box, order 18 in time,
`abstol = 1e-14`, `parse_eqs=false`), time-`T` flow of `(U,p,y)` from a crest box of half-width `delta` in `U`,
`delta·k` in `y`. Run under `research-run --cpu 1 --memory-gb 2`. This is an enclosure of the time-`T` flow, not yet
of the return map (the crossing-time step is missing), and the library is not audited.

| mu | rho | T | U_min | steps | time per evaluation | widths at T for `delta = 1e-3` (U, p, y) | widths for `delta = 1e-2` |
|---|---|---|---|---|---|---|---|
| 1 | 1.2 | 6.7 | −2.8 | 29 | 1.5 s | 6.3e-3, 3.0e-2, 8.5e-2 | 5.7e-2, 0.28, 0.83 |
| 1 | 2 | 8.8 | −13 | 40 | 1.0 s | 1.2e-2, 0.15, 0.71 | 0.12, 1.5, 7.2 |
| 1 | 8 | 24.1 | −285 | 75 | 1.9 s | 5.7e-2, 5.3, 89 | 1.5, 63, 1.1e3 |
| 0.2 | 1.5 | 19.2 | −8.9 | 49 | 2.7 s | 9.3e-3, 5.1e-2, 0.11 | 9.1e-2, 0.50, 1.1 |
| 0.05 | 1.5 | 56.7 | −23 | 120 | 7.3 s | 1.9e-2, 0.10, 0.18 | 0.18, 0.96, 1.8 |
| 0.05 | 3 | 112 | −174 | 182 | 9.3 s | 9.7e-2, 1.1, 3.5 | 1.3, 11, 39 |

Findings. (a) The enclosure survives a full period, including a trough at `U = −285`, with no blow-up. Widths scale
linearly in `delta` from `1e-6` to `1e-3` in all six cases; at `delta = 1e-2` the two deep-trough cases are
superlinear by a factor 1.3–2.6 (so usable boxes there are `≲ 3e-3`). The large `p`, `y` widths are the phase shear
along the flow, which the section projection removes. (b) 30–190 steps per period at order 18. (c) 25–60 ms per
step in this unoptimised Julia set-up.
A doubles-based `C^0`/`C^1` Lohner integrator (CAPD) typically costs tens of microseconds per step for a 3-D field of
this size; I assume `≈ 10 ms` per evaluation there. **That figure is an assumption.**

### 4.8 Cost summary and failure modes
- With CAPD-class speed: `6.7e5–6.7e7` evaluations × 10 ms = 2 CPU-hours to 8 CPU-days for `mu0 = 0.05`, `r_c = 0.2`.
- With TaylorModels.jl as measured (≈3 s): `2e6–2e8` s; only the optimistic end (3 weeks) is usable.
- Failure modes: (i) the (C) upper bound from the tube envelope may be too crude near Hopf, forcing a true `C^1`
  integrator (×5 cost); (ii) spikes at `rho ≳ 10` (`x_max ≈ k e^{U*}` up to `3·10^4` at `mu = 1`, `rho = 20`) make
  `p' = k − mu p − y` fast and may force `10^3` steps per pulse; (iii) parameter intervals in `mu` shift the period;
  (iv) if the collar cannot reach `r_c = 0.2`, the count grows like `r_c^{−1.3}`.

## 5. The Hopf collar `H` (`r ≤ r_c`, `mu ∈ [mu0,1]`)

Target: every cycle with `r ∈ (0,r_c]` is hyperbolic attracting. Two layers.

1. *Confinement.* For `rho ∈ [1, 1+r_c eps^2]` no crest has `U ≥ a_c`. By (E) on `R ∩ {U ≥ a_c}` exactly as in the
   core; the margin there is `≈ eps^3 U (Rsel-like(U) − r)` and stays positive if `a_c^2 ≳ 4 r_c(1+4 eps^2)`
   (`a_c ≈ 1.0–1.3` for `r_c = 0.2`).
2. *Small amplitude.* In the blown-up variables `(a, yh, r)` on `[0,a_c]×[yh]×[0,r_c]` the scalar factor is
   `m − 1 = eps^3 (c1 r − c2 a^2 + …)`. Away from the corner `(a,r) = (0,0)` ordinary (E)/(C) boxes with relative
   accuracy suffice. At the corner one needs a **validated third-order expansion in `a`** of the return map
   (Taylor models in `(a,r)` with interval remainder, or CAPD `C^3` jets). The leading coefficients are known in
   closed form (`l_1 = −sigma/((sigma^2+1)(sigma^2+4))`, `d Re z/d rho = sigma/(2(1+sigma^2))`, v2.0 App. A), are
   nonzero uniformly on `mu ∈ [mu0,1]`, and only the remainder has to be bounded.

Alternative with a track record: Fourier-space Newton–Kantorovich with desingularised amplitude
(van den Berg–Jaquette for Wright near `pi/2`; van den Berg–Lessard–Queirolo for Hopf in ODEs) plus a-priori bounds
on Fourier coefficients of all small cycles. Either way this is a contained project; expected cost is small compared
with the core (three variables `(mu, r, a)`), but the corner analysis is where the care goes. Not sized by a probe.

## 6. Large gain `L` — interface only

Needed from the large-gain lane: an explicit continuous `rho1(mu)` on `[mu0,1]` with "every cycle is hyperbolic
attracting for `rho ≥ rho1(mu)`". Remarks for whoever proves it:
- (C) needs only an *upper* envelope of `x(t)` along the cycle (cooperative compound system). So a-priori shape
  bounds on an arbitrary cycle (pulse height, pulse width, quiescent length), not asymptotics of a particular one,
  are what is required — the same logic as Xie's large-`alpha` argument for Wright.
- The core can be extended towards larger `rho` cheaply in box count (`N_state ≈ 4`, `d_rho ≈ rho`), the limit being
  pulse stiffness. `rho1 ≈ 20–50` looks reachable; `rho1 ≈ 200` (`U_min ≈ −2·10^5`, `T ≈ 525` at `mu = 1`) only with the
  closed-form quiescent phase.

## 7. The small-`mu` wedge `W` — the part that decides whether the full conjecture is reachable

| piece | statement needed | available | missing |
|---|---|---|---|
| `W_H ∪ W_A` (`r ≤ r2`) | every cycle attracting | NOTES Thm A (proved), Thm C (proved, local, `mu0` not explicit), Claim C′ (outline) | (a) Claim C′ written out; (b) **(APB-N)**: no *non-attracting* cycle leaves `{H0 ≤ h_b}` — reduced in NOTES §4.2b to two a-priori estimates (period `≤ C/mu`, `p ≥ −C`) via the v2.0 variance window; (c) explicit constants |
| `W_M` (`r ≥ r2`, `rho−1 ≤ delta`) | every cycle attracting | formal law `h* ≈ 5r/4`, `log m_t ≈ −2.108(rho−1)^{3/2}`; one numerical match | uniform third-order averaging on an unbounded energy range; overlap with `W_B` |
| `W_B` (`delta ≤ rho−1`, `rho ≤ rho1`) | every cycle attracting | singular map `F_lambda`, universal pulse map `S` with `S(1) = 2` exactly; `Lambda` increasing on a grid (numerical); 3-digit agreement with true cycles | (B1) `Lambda' > 0` (1-D interval proof + two asymptotic ends); (B2) the true return map is a `C^1`-small perturbation of `F_lambda`, error `O(mu log 1/mu)`, explicit; (B3) every cycle is of pulse type |
| `W_L` (`rho ≥ rho1`) | every cycle attracting | – | large-gain argument uniform as `mu → 0` |

Comments.
- Explicit constants in `W_A` are attainable without Fenichel: the fast equation is linear in `Q`, so for a bounded
  entire solution `Q = −g + eps g' − eps^2 g'' + eps^3 g''' − eps^3 ∫ e^{−(tau−s)/eps} g''''(s) ds` with
  `g = eps P + rho G(U)`, and every `sup|g^{(n)}|` is bounded by a polynomial in `sup|U|, sup|P|` by induction.
  The Toda integrals `J0, J1` and the remainder integrals can then be enclosed by interval quadrature.
- The true gap in `W_A` is large cycles (NOTES §4.2). Excluding them outright (APB) would need `W_M`-type uniform
  estimates at the lower end of the energy window `h_b < h ≲ 1/eps^2` and relaxation estimates at its upper end.
  The cheaper target is (APB-N): a non-attracting cycle has `<(eta−mu)^2> ≤ r mu^2/2` (variance window), which
  through the exact relation `rho e^U − s = A L^2 qh` forbids spikes; what is missing is a period bound and a lower
  bound on `p` for arbitrary cycles. This route works only for bounded `r`, i.e. exactly in `W_H ∪ W_A`.
- (B2) is a regular perturbation once the section is placed before the pulse and `q ≥ delta'`: the pulse lasts
  `O(log(1/mu)/q)` and is perturbed by `O(mu)` terms; the quiescent phase is explicit up to `∫ x dt = O(1/q)`.
  It degenerates as `q → 0`, i.e. exactly at the matching zone.
- Pushing `mu0` down by computation is possible but does not close the wedge: per unit `mu` the ideal count grows
  like `mu^{−2.8}` and the period like `1/mu`.

Effort. This block is several papers' worth of singular-perturbation analysis with explicit constants. I would
not start it before the core + collar + large-gain theorem exists for `mu ≥ mu0`.

## 8. Alternatives considered

- *Amplitude formulation.* Parametrise cycles by crest amplitude `a` and solve for `rho` (as `ampcycle.py` does).
  Both singular limits are graphs over the amplitude (`r = Rsel(h)` proved monotone; `rho = Lambda(q)` numerically
  monotone), and the Hopf point becomes regular. The claim "(AMP) one `rho` per amplitude, increasing in `a`" implies
  the conjecture but is stronger; a computer-assisted proof of (AMP) would have to search `rho` for each `a`, which
  is the same dimension count as Section 4. I see no saving, except in the collar.
- *One contraction metric on a trapping annulus.* The reduced crest map expands near `E` (`m_E > 1`) and contracts
  near the cycle, so no single metric contracts on an annulus reaching `E`; the usable form is exactly
  "(C) on the window, (E) outside it" of Section 4.
- *Pointwise certificates* (Dulac function on the attractor disc, or a compound supersolution `omega`, brief 3(d),(e)).
  They would replace integration by interval evaluation of one differential inequality. The dulac and planar lanes
  report numerical existence at sampled parameters but no family that is robust in parameters, and a function of
  `(eta, zeta)` lives on an attractor that is not known in closed form. Not ready as a computational route. In the
  Toda wedge the leading-order Dulac function is `1/J0(Ht)` (NOTES §6).
- *Concavity of the reduced crest map* (planar lane, numerical, 14/14): if proved it gives the conjecture outright
  and makes this program unnecessary; it is the cheapest thing to try to refute numerically at more parameters.

## 9. Order of work and what each step would establish

1. Prototype at `mu = 1`: certify `rho ∈ [1.05, 20]` (ideal count ≈ 300; `3·10^4–3·10^6` evaluations). First rigorous
   uniqueness statement on a segment; calibrates the overhead factor that Section 4.5 cannot.
2. Two-parameter core for `mu ∈ [0.2,1]`, then `[0.05,1]`.
3. Collar in blown-up coordinates; validated third-order corner.
4. Large-gain theorem with explicit `rho1(mu)`; extend the core to meet it.
   → theorem: uniqueness for `mu ∈ [mu0, 1/mu0]`, all `rho > 1`.
5. Wedge: `W_A` with explicit constants and the Hopf end; then (B1); then (B2)–(B3); `W_M` and `W_L` last.

Step 1 is days of work for someone with a CAPD set-up. Steps 2–3 are weeks. Step 4 depends on the large-gain lane.
Step 5 is open-ended; I do not have an argument that its corners (`mu → 0` with `rho → 1` beyond the Toda scale, and
`mu → 0` with `rho → ∞`) are free of surprises.
