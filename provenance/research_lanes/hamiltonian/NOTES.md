# Lane `hamiltonian` — NOTES (2026-10-03)

Everything for this lane is in this directory. Python: `/Users/admin/.venvs/analysis-py313/bin/python <script>`.
All numerics are ordinary floating point unless stated; the only interval computation is the Taylor-model probe
(`tm_probe*.jl`), which is a cost probe, not a certificate. The program for the full conjecture is in `PROGRAM.md`.

## 0. Bottom line

1. **The first-order Melnikov statement in the lane prompt is false, and this is proved.** On every level curve of
   the Toda centre, `oint e^U U'^2 = oint U'^2`, so
   `M(h) = oint (U'^2 − U''^2) dtau = (1 − rho) oint U'^2 dtau < 0` for every `h > 0` and every `rho > 1`.
   `M` has **no** positive zero and the ratio `oint U'^2 / oint (e^U−1)^2` is the constant `rho`, not a monotone
   function (Prop. 2). (The literature lane found the same identity numerically, independently.) Consequently, at
   fixed `rho > 1` there is no cycle of bounded energy as `mu → 0`; numerically the crest energy grows like
   `c(rho)/mu` (Section 2.2).
2. **The correct Hamiltonian regime is the wedge `rho − 1 = r·eps^2`, `eps = sqrt(mu)/(1+mu)`** (so `rho − 1 ≈ r·mu`),
   and the selection happens at **third** order (first order is sign-definite, second order vanishes identically).
   The cycle energy is the root of `r = Rsel(h) := (2/3) <p^4>/<p^2>` on the level `H0 = p^2/2 + e^U − U − 1 = h`
   (Section 3). Checked against true cycles down to `mu = 1e-4`: relative error `4e-4`, decaying like `eps^2`.
3. **Theorem A (proved, elementary, complete):** `Rsel` is a real-analytic bijection `(0,∞) → (0,∞)` with
   `Rsel'(h) ≥ 4/5` for every `h > 0`. So for every `r > 0` the third-order selection equation has exactly one
   root `h*(r)`, it is simple, and the sign is the attracting one. The proof is a two-line moment inequality
   (kurtosis of the velocity is at most 9/5 because the potential is convex); no Abelian-integral criterion is needed.
4. **Theorem C (proved given Fenichel's theorem; local in phase space):** for `r` in a compact subset of `(0,∞)`
   and `mu` small, the original system has exactly one cycle inside any fixed energy ball `{H0 ≤ h_b}`, it is
   hyperbolic attracting, and `log m_t = −eps^3 J0(h*) Rsel'(h*) + O(eps^4) < 0`.
5. **Sharpest obstruction:** Theorem C is not global. Nothing proved here excludes, for small `mu` and
   `rho − 1 = O(mu)`, an additional *large non-attracting* cycle (crest energy between `h_b` and the a-priori ceiling
   of order `1/mu`). Excluding all large cycles would need third-order averaging uniform on an unbounded energy range
   plus the relaxation regime. A cheaper reduction is given in §4.2b: a non-attracting cycle obeys the v2.0 variance
   window `<(eta−mu)^2> ≤ r mu^2/2`, and through three exact formulas (proved, checked numerically) this forbids
   spiky cycles — but closing the argument still needs two a-priori estimates for arbitrary cycles
   (period `≤ C/mu`, `p ≥ −C`) that I did not prove.
6. **The other half of the small-`mu` limit is a relaxation oscillation, not a Hamiltonian one** (Section 5).
   At fixed `rho > 1`, `mu → 0`: the cycle is a pulse (universal scattering map `S` of `v''' + v'' + e^v = 0`)
   followed by an explicit slow return. The singular 1-D return map reproduces `m_t`, `mu·T` and `mu·U_min` of true
   cycles to 3 digits. Its fixed points are parametrised by the pulse amplitude `q`, with `rho = Lambda(q)`;
   `Lambda` is increasing on a 33-point grid `q ∈ [0.02, 12]` (numerical, not proved). The two regimes match:
   `log m_t ≈ −2.108 (rho−1)^{3/2}` in the overlap.
7. **By-product (proved): an exact identity for `k − k_H` on every cycle at every parameter** (Lemma 1(f)).

Status of the lane: `theorem` (Theorem A, Lemma 1, Prop. 2 complete; Theorem C complete given a cited standard
theorem). The full conjecture is **not** proved, not even for small `mu`.

---------------------------------------------------------------------------------------------------------------

## 1. Setting and exact identities

Frequency-one time `tau = sqrt(mu) t`, `sigma = (1+mu)/sqrt(mu) ≥ 2`, `eps = 1/sigma ∈ (0, 1/2]`,
`rho = k/k_H`. The scalar equation (paper eq. (1.7) in `tau`-time; Appendix A) is

    (S)   eps U''' + U'' + eps U' + rho (e^U − 1) = 0          (' = d/dtau).

`mu ∈ (0,1] ↔ eps ∈ (0,1/2]` is a bijection. Write `P = U'`, `Q = U''`, `G(U) = e^U − 1`, `Phi(U) = e^U − U − 1`,
`<·>` = mean over one period. In `t`-variables `P = p/sqrt(mu)`, `Q = p'/mu`, `Q' = p''/mu^{3/2}`.

**Lemma 1 (exact identities on every cycle, all parameters). PROVED.**

    (a) <e^U> = 1
    (b) <Q^2> = <P^2>
    (c) <Q^2> = rho <e^U P^2>,   hence   <(rho e^U − 1) P^2> = 0
    (d) (rho − 1) <P^2> = eps <P^3> − 2 eps <Q^2 P>
    (e) <e^U P^3> = −(2 eps/rho) ( <Q'^2> − <Q^2> )
    (f) (rho − 1) <P^2> = eps^2 [ (1/rho)(<P^4> − 3<Q^2 P^2>) − <Q^3> + 2(1 − 1/rho)(<Q'^2> − <Q^2>) ]

In the brief's variables (f) reads
`(k − k_H) A <p^2> = (1/rho)(mu<p^4> − 3<p'^2 p^2>) − <p'^3> + 2(1 − 1/rho)(<p''^2> − mu<p'^2>)`.

*Proof.* All by multiplying (S) by a multiplier, averaging, and integrating by parts on the period.
(a): multiplier 1. (b): multiplier `P`; `<PQ'> = −<Q^2>`, `<PQ> = 0`, `<G P> = <(e^U − U)'> = 0`, leaving
`eps(<P^2> − <Q^2>) = 0`. (c): multiplier `Q`; `<QQ'> = <PQ> = 0` and `<G Q> = <e^U P'> = −<e^U P^2>`.
(d): by (b),(c), `(rho−1)<P^2> = −rho <G P^2>`; substitute `rho G = −(eps Q' + Q + eps P)` and use
`<Q P^2> = 0`, `<Q' P^2> = −2<Q^2 P>`. (e): multiplier `Q'`; `<e^U Q'> = −<e^U P Q> = <e^U P^3>/2`.
(f): `<P^3> = <e^U P^3> − <G P^3>` and `−rho<G P^3> = <(eps Q' + Q + eps P)P^3> = eps(<P^4> − 3<Q^2P^2>)`;
`<Q^2 P> = −<(eps Q' + eps P + rho G) Q P> = (eps/2)<Q^3> + (rho/2)<e^U P^3>`; insert both and (e) into (d). ∎

*Numerical test* (`exact_identity.py` → `exact_identity.out`): 8 cycles, `mu ∈ {0.05,…,3}`, `rho ∈ {1.1,…,5}`;
(f) holds with relative error `≤ 1.1e-11`, (b) and (c) to `≤ 4e-12`. The algebraic core of (c) is also checked
symbolically (`normal_form_check.py`, item (3)).

Remark. (c) is the exact form of the "first-order Melnikov condition"; (f) is the exact form of the third-order
selection rule of Section 3 (on Toda orbits its right side reduces to `(2/3) eps^2 <P^4>`). Neither bounds the
amplitude by itself.

## 2. The first-order statement in the lane prompt is false

**Prop. 2. PROVED.** Let `rho > 0` and let `gamma_h` be the closed orbit `P^2/2 + rho Phi(U) = h > 0` of
`U'' + rho(e^U − 1) = 0`. Then `oint e^U P^2 dtau = oint P^2 dtau` and `oint Q^2 dtau = rho oint P^2 dtau`. Hence

    M(h) := oint (U'^2 − U''^2) dtau = oint [U'^2 − rho^2 (e^U−1)^2] dtau = (1 − rho) oint U'^2 dtau,

which for `rho > 1` is strictly negative for every `h > 0`. The ratio `oint U'^2 / oint (e^U−1)^2` equals `rho`
identically. In particular `M` has no positive zero, and the ratio is not monotone (it is constant).

*Proof.* On `gamma_h`, `G = −Q/rho`, so `oint G P^2 dtau = −(1/rho) oint P^2 P' dtau = 0`. And
`oint Q^2 = −oint P Q' = rho oint e^U P^2 = rho oint P^2`. ∎

Meaning: the energy balance `dH/dtau = −eps P (U''' + U')` gives a first-order drift `−eps M(h) = eps(rho−1) oint P^2 > 0`
at *every* energy. With `rho − 1 = O(1)` the Hamiltonian approximation predicts unbounded growth and no cycle; the
actual cycle sits where the perturbation is no longer small.

### 2.1 Numerical confirmation on the unperturbed centre
`toda.py` → `toda.out`: for `h ∈ [1e-3, 1e3]`, `oint e^U P^2 / oint P^2 = 1.000000000` and
`oint (e^U−1)^2 / oint P^2 = 1.000000000` (rho = 1 normalisation), 12 values of `h`.

### 2.2 Numerical confirmation on true cycles at fixed `rho` (refutes "energy tends to a root of M")
`fixed_rho.py` → `fixed_rho.out` (`find_cycle` from whlib; crest energy `h = Phi(U_max)`):

| rho | mu | U_max | U_min | h | h·eps^2 | mu·U_min | mu·T | m_t |
|---|---|---|---|---|---|---|---|---|
| 1.5 | 0.1 | 2.553 | −13.5 | 9.29 | 0.768 | −1.352 | 3.193 | 0.5692 |
| 1.5 | 0.01 | 4.480 | −104.5 | 82.7 | 0.811 | −1.045 | 2.521 | 0.6669 |
| 1.5 | 0.001 | 6.737 | −1032.5 | 835.1 | 0.833 | −1.032 | 2.447 | 0.6785 |
| 3 | 0.1 | 4.093 | −96.4 | 54.8 | 4.531 | −9.639 | 5.943 | 0.1379 |
| 3 | 0.01 | 6.195 | −801.7 | 483.1 | 4.736 | −8.017 | 5.331 | 0.1810 |
| 3 | 0.001 | 8.476 | −7882 | 4789 | 4.779 | −7.882 | 5.269 | 0.1861 |

The energy diverges like `1/eps^2`; `m_t` tends to a limit in `(0,1)` (not to 1). The exact identity (c) holds on
all of these cycles (`|<(rho e^U−1)P^2>|/<P^2> ≤ 6e-7`, typically `1e-13`).

## 3. Third-order selection in the wedge `rho = 1 + eps^2 r`

### 3.1 Derivation (exact algebra; machine-checked in `normal_form_check.py`)
First-order system `U' = P, P' = Q, eps Q' = −Q − eps P − rho G(U)`, `rho = 1 + eps^2 r`, `E = e^U`.
Slow-manifold series `Q = w0 + eps w1 + eps^2 w2 + eps^3 w3 + O(eps^4)`:

    w0 = −G,   w1 = G P,   w2 = −r G − E P^2 + G^2,   w3 = r E P + E P^3 − 4 G E P − G^2 P.

Reduced field `X = P ∂_U + w ∂_P`. With `H0 = P^2/2 + Phi(U)`, `h1 = P^3/3`,
`h2 = r Phi(U) − Psi(U) + s(U,P)`, `Psi' = G^2`, `s = sigma(U, H0)`,
`sigma(U,h) = int_{U_−(h)}^{U} 2(h − Phi(u)) du` (smooth on `H0 > 0`; `X0 s = P^3`), the modified energy
`Ht = H0 + eps h1 + eps^2 h2` satisfies

    X(Ht) = eps^3 F3 + O(eps^4),    F3 = r P^2 − 4 G E P^2 + 2 G P^2 (U − U_−(H0)).

Averaging over the Toda level `H0 = h` (using `oint G P^2 = 0`, `oint G U P^2 = (1/3) oint P^4`,
`oint G E P^2 = (1/3) oint E P^4 = (1/3) oint P^4`):

    M3(h) := oint F3 dtau = r J0(h) − (2/3) J1(h),    J0 = oint P dU = oint P^2 dtau,   J1 = oint P^3 dU = oint P^4 dtau.

So the first-order Melnikov function is `eps(rho−1)J0 = eps^3 r J0` (sign-definite, Prop. 2), the second-order
function vanishes identically, and the selection equation is

    r = Rsel(h) := (2/3) J1(h)/J0(h) = (2/3) <P^4>/<P^2>      on H0 = h.

Sanity: as `h → 0`, `Rsel ≈ h = a^2/2`, i.e. `a^2 = 2r`, which is the Hopf normal form of Appendix A of v2.0
(`l_1 ≈ −eps^3`, `d Re z/d rho ≈ eps/2`) at leading order in `eps`.

### 3.2 Theorem A (PROVED — complete, elementary)

**Theorem A.** Let `Phi` be any `C^2` strictly convex function on `R` with `Phi(0) = Phi'(0) = 0`, `Phi → ∞` at
`±∞` (in particular `Phi = e^U − U − 1`). For `h > 0` let `<·>_h` be the time average over the periodic orbit
`P^2/2 + Phi(U) = h` of `U'' + Phi'(U) = 0`, and `Pm = sqrt(2h) = max P`. Then

    (i)   Pm^2/3 ≤ <P^2>_h < Pm^2;
    (ii)  <P^4>_h ≤ (6/5) Pm^2 <P^2>_h − Pm^4/5 ≤ (9/5) <P^2>_h^2;
    (iii) with J0(h) = oint P dU, J1(h) = oint P^3 dU and Rsel = (2/3) J1/J0:
          J0' = T0(h) (the period), J1' = 3 J0, and
          Rsel'(h) = (2/3) (3 − <P^4>_h/<P^2>_h^2) ≥ 4/5;
    (iv)  (4/5) h ≤ Rsel(h) ≤ (4/3) h; Rsel is a strictly increasing bijection of (0,∞) onto itself.

Hence for every `r > 0` the equation `M3(h) = J0(h)(r − Rsel(h)) = 0` has exactly one root `h*(r) > 0`,
`h*(r) ≤ (5/4) r`, it is simple, and `M3'(h*) = −J0(h*) Rsel'(h*) ≤ −(4/5) J0(h*) < 0` (energy is restored from
both sides: the attracting sign).

*Proof.* Let `U_−(phi) < 0 < U_+(phi)` be the two solutions of `Phi(U) = phi`, `phi ∈ (0,h]`.
*Time density of the velocity.* By the reversibility `(U,P,tau) → (U,−P,−tau)` it suffices to look at the half
orbit `P ≥ 0`, of duration `T0/2`. It consists of the arc from `(U_−(h),0)` to `(0,Pm)`, on which `P' = −Phi'(U) > 0`,
and the arc from `(0,Pm)` to `(U_+(h),0)`, on which `P' < 0`. On each arc `P` is a strictly monotone coordinate and
`dtau = dP/|Phi'(U)|` with `U = U_∓(h − P^2/2)`. As `P` increases, `phi = h − P^2/2` decreases, `U_∓(phi) → 0`
monotonically, and `|Phi'(U_∓)|` decreases because `Phi'` is increasing and vanishes at 0. So for every continuous `F`,

    int_0^{T0/2} F(P(tau)) dtau = int_0^{Pm} F(P) f(P) dP,   f(P) := 1/|Phi'(U_−(h − P^2/2))| + 1/Phi'(U_+(h − P^2/2)),

with `f > 0` **strictly increasing** on `(0,Pm)`, integrable, `int_0^{Pm} f = T0/2`.
*(i).* Chebyshev's integral inequality for the two increasing functions `P^2` and `f` on `(0,Pm)`:
`(1/Pm) int P^2 f ≥ (1/Pm) int P^2 · (1/Pm) int f`, i.e. `<P^2> ≥ Pm^2/3`. (Proof of Chebyshev:
`int int (P^2 − S^2)(f(P) − f(S)) dP dS ≥ 0`; all integrals converge since `f` is integrable.) The upper bound is
trivial.
*(ii).* An increasing integrable `f` on `(0,Pm)` is a mixture of indicators: `f(P) = nu([0,P])` for a nonnegative
measure `nu` on `[0,Pm)` (take `f` right-continuous). By Tonelli, for `j = 2,4`,
`int_0^{Pm} P^j f dP = int (int_a^{Pm} P^j dP) dnu(a)`. With the probability measure
`dpi(a) = (2/T0)(Pm − a) dnu(a)` and `m_j(a) = (Pm − a)^{-1} int_a^{Pm} P^j dP`, this is `<P^j> = int m_j dpi`.
Put `t = a/Pm ∈ [0,1)`: `m_2 = Pm^2 (1+t+t^2)/3`, `m_4 = Pm^4 (1+t+t^2+t^3+t^4)/5`, and

    (6/5) Pm^2 m_2 − Pm^4/5 − m_4 = (Pm^4/5) t (1+t)(1−t^2) ≥ 0.

Integrating against `pi` gives the first inequality of (ii). The second is
`(9/5) g^2 − (6/5) g + 1/5 = (1/5)(3g − 1)^2 ≥ 0` with `g = <P^2>/Pm^2`.
*(iii).* `J0(h) = 2 int_{U_−}^{U_+} sqrt(2(h − Phi))dU`, `J1(h) = 2 int (2(h − Phi))^{3/2} dU`; the integrands
vanish at the endpoints, so `J1' = 3 J0` and `J0' = 2 int dU/sqrt(2(h − Phi)) = T0`. With `J0 = T0<P^2>`,
`J1 = T0<P^4>`: `(J1/J0)' = (3 J0^2 − J1 T0)/J0^2 = 3 − <P^4>/<P^2>^2 ≥ 3 − 9/5` by (ii).
*(iv).* `Rsel(0+) = 0` because `Rsel ≤ (2/3) Pm^2 = (4/3) h`; integrate (iii). Analyticity of `J0, J1` in `h > 0`
is standard for analytic `Phi` (not needed for the conclusions). ∎

Remarks. (1) The constant `4/5` is sharp: `Rsel'(h) → 4/5` as `h → ∞` (velocity density → uniform) and
`Rsel'(0+) = 1` (arcsine density). (2) The theorem uses only convexity of the potential; neither the Li–Zhang nor
the Grau–Mañosas–Villadelprat criterion is needed. Applied to the *first-order* pair `(oint P dU, oint e^U P dU)`
the GMV criterion degenerates, as it must: `l_1 ≡ l_0` because the two integrals coincide (Prop. 2).
(3) `Rsel' > 0 ⟺ J1 = oint P^3 dU` is strictly log-concave in `h`.

*Numerical test, run before the proof was written* (`toda.py` → `toda.out`): `<P^4>/<P^2>^2` increases from
1.5000 (`h = 1e-3`) to 1.7988 (`h = 1e3`), always `< 9/5`; `Rsel/h` decreases from 1.0000 to 0.8045;
`<P^2>/(2h/3)` decreases to 1.006 at `h = 1e3` (so (i) is asymptotically sharp too).

### 3.3 Numerical verification against true cycles (job 1(a), in corrected form)
`ampcycle.py` → `ampcycle.out`. Cycles are computed by amplitude-fixed shooting (crest `U_max = a` prescribed,
unknowns `y0, k`; uses `whlib.crest_return`), which is well conditioned at Hopf. Columns: `r_num = (rho−1)/eps^2`
of the true cycle, prediction `Rsel(Phi(a))`; `log m_t / eps^3` of the true cycle, prediction
`T0 (2/3)(<P^4> − 3<P^2>^2)/<P^2> = −J0 Rsel'`.

| a | mu | r_num / r_pred | (log m_t/eps^3) / pred |
|---|---|---|---|
| 0.3 | 1e-2 | 0.96373 | 0.95298 |
| 0.3 | 1e-3 | 0.99619 | 0.99501 |
| 0.3 | 1e-4 | 0.99962 | 0.99950 |
| 1 | 1e-3 | 0.99599 | 0.99261 |
| 1 | 1e-4 | 0.99960 | 0.99925 |
| 2 | 1e-3 | 0.99302 | 0.97829 |
| 2 | 1e-4 | 0.99929 | 0.99775 |
| 3 | 1e-3 | 0.98318 | 0.93877 |
| 3 | 1e-4 | 0.99821 | 0.99329 |

Errors scale like `eps^2` (factor ≈ 10 per decade of `mu`), consistent with the `(tau,eps) → (−tau,−eps)` symmetry
of (S). Also `sqrt(mu)·T → T0(h)` (last columns of the output). The error grows with `a`: the expansion parameter
is `eps·P ~ eps·sqrt(2h)`, not `eps`.

## 4. What is rigorous for the original system at small `mu > 0` (job 1(c))

### 4.1 Theorem C (PROVED, given Fenichel's theorem; local in phase space)

**Theorem C.** Fix `0 < r1 ≤ r2 < ∞` and `h_b > h*(r2)`. There is `mu0 > 0` such that for every `mu ∈ (0, mu0)` and
every `rho = 1 + r eps^2` with `r ∈ [r1,r2]` (`eps = sqrt(mu)/(1+mu)`):
(1) the system has **exactly one** cycle whose projection `(U, U_tau) = (U, p/sqrt(mu))` stays in
`D_b = {H0 ≤ h_b}`;
(2) its crest energy is `h*(r) + O(eps)`;
(3) it is hyperbolic attracting, with `log m_t = −eps^3 J0(h*(r)) Rsel'(h*(r)) + O(eps^4) ≤ −(4/5) eps^3 J0(h*) + O(eps^4)`.

*Proof.* Throughout, `k_s ≥ 6` is a fixed smoothness index and `K = {H0 ≤ h_b + 1}`.

*Step 0 (slow manifold).* (S) is the fast–slow system `U' = P, P' = Q, eps Q' = −Q − eps P − (1+eps^2 r)G(U)`
with critical manifold `Q = −G(U)`, a graph over the `(U,P)` plane with fast eigenvalue `−1` (normally hyperbolic,
attracting); `r` is a parameter. By Fenichel's theorem (Fenichel 1979, Thm 9.1; Jones 1995, Thms 1–2) there are
`eps0 > 0` and `w ∈ C^{k_s}(K × [0,eps0] × [r1,r2])`, `w|_{eps=0} = −G`, such that `M_eps = {Q = w(U,P,eps,r)}` is
locally invariant over `K` for `eps ∈ (0,eps0]`. The `eps`-Taylor coefficients of `w` at `eps = 0` are forced by the
invariance equation `eps(w_U P + w_P w) = −w − eps P − rho G`, so they are the `w_j` of Section 3.1 and
`w = w0 + eps w1 + eps^2 w2 + eps^3 w3 + eps^4 w~`, `w~ ∈ C^{k_s−4}`.

*Step 1 (every cycle over `K` lies on `M_eps`).* For any solution with `(U,P) ∈ K`, `D := Q − w(U,P,eps,r)`
satisfies, by subtracting the invariance equation, the exact linear equation `D' = −(1/eps + w_P) D`. With
`C = sup_K |w_P|` and `eps < 1/(2C)`, `|D(tau)| ≤ |D(tau0)| e^{−(tau−tau0)/(2 eps)}`. For a periodic solution with
projection in `K`, letting `tau0 → −∞` gives `D ≡ 0`. The same argument applied to the equilibrium gives
`w(0,0,eps,r) = 0`. So cycles with projection in `K` are exactly the cycles in `K` of the planar `C^{k_s}` field
`X_eps: U' = P, P' = w`, which at `eps = 0` is the Toda centre. Their second transverse multiplier is
`exp(−oint (1/eps + w_P) dtau) = e^{−T/eps + O(eps^2)}` (since `w_P = eps G + O(eps^2)` and `oint G dtau = 0` by
Lemma 1(a)); together with Step 2 this is consistent with `m_s m_t = e^{−sigma T}`.

*Step 2 (annulus).* Fix `0 < h_a < h*(r1)`. On `K minus {H0 < h_a/2}` the functions of Section 3.1 are smooth and
`X_eps Ht = eps^3 F3 + eps^4 R4` with `R4` of class `C^1` up to `eps = 0` (Hadamard's lemma; `k_s ≥ 6`). The section
`S = {P = 0, U > 0}` is transverse (`P' = w = −G + O(eps) < 0`), `Ht` is a regular coordinate on it, and the return
map in this coordinate is `c ↦ c + eps^3 Dsp(c,eps,r)`,
`Dsp = int_0^{T_eps(c)} (F3 + eps R4) dtau`, a `C^1` function on `[c_a,c_b] × [0,eps0] × [r1,r2]` with
`Dsp(c,0,r) = M3(c;r) = J0(c)(r − Rsel(c))`. By Theorem A, `M3(·;r)` has exactly one zero in `[h_a,h_b]`, with
derivative `≤ −(4/5)J0(h_a) < 0` there. By uniform convergence and the implicit function theorem, for `eps` small
`Dsp(·,eps,r)` has exactly one zero `c*(eps,r)`, `∂_c Dsp(c*) < 0`, `c* → h*(r)`. This gives (2), a cycle with
`Pi'(c*) = 1 + eps^3 ∂_c Dsp ∈ (0,1)`, and
`log m_t = eps^3 M3'(h*) + O(eps^4) = −eps^3 J0(h*) Rsel'(h*) + O(eps^4)`.

*Step 3 (no cycle with crest energy `< h_a`).* In polar coordinates about the origin the angular speed of `X_0` is
`−(U G(U) + P^2)/(U^2+P^2) < 0` on `{0 < H0 ≤ h_a}` and tends to `−1` at the origin, so for `eps` small the return
map `Pi_eps` on the ray `{P = 0, U = a > 0}` is defined on `(0,a_a]`, extends by `Pi_eps(0) = 0`, and
`Pi_eps(a) − a = a·g(a,eps,r)` with `g ∈ C^{k_s−1}([0,a_a] × [0,eps0] × [r1,r2])` (integrate `da/dtheta`).
For every fixed `a > 0`, Step 2 gives `g(a,eps) = eps^3 g3(a) + O(eps^4)`, `g3(a) = M3(Phi(a))/(a G(a)) > 0`
(because `Phi(a) ≤ h_a < h*(r)`). Since the `eps`-derivatives of `g` at `eps = 0` are continuous in `a` on `[0,a_a]`,
`∂_eps^j g(a,0) = 0` for `j ≤ 2` and `∂_eps^3 g(a,0)/6 = g3(a)` hold on the closed interval, with
`g3(0) = lim_{a→0} M3(Phi(a))/(a G(a)) = pi r > 0`. Taylor's formula in `eps` gives
`g = eps^3 (g3 + eps g4)` with `g4` bounded and `g3 ≥ c > 0` on the compact set, so `g > 0`: no fixed point in `(0,a_a]`.

*Step 4 (conclusion).* A cycle of the planar field in `K` encircles the only equilibrium, so it crosses `S`; its
crest has energy in `(0,h_a)` (excluded by Step 3) or in `[h_a,h_b]` (unique by Step 2). The cycle of Step 2 lies in
`D_b` because it is `O(eps)`-close to `gamma_{h*(r)}` and `h*(r) < h_b`. ∎

Status and caveats. Every model-specific step is written out and the algebra is machine-checked. The proof cites
Fenichel's persistence theorem (with `C^k` dependence on `eps` up to `eps = 0`) and uses standard smoothness of
planar return maps; I have not re-derived those. `mu0` is not explicit. The numerical tests of Section 3.3 were run
first and agree to `4e-4`.

### 4.2 What Theorem C does NOT give — the sharpest obstruction of this lane

Theorem C is a statement about cycles inside a fixed energy ball. **Global uniqueness in the wedge
`rho − 1 ≤ r2 eps^2` would follow from Theorem C plus the following missing statement:**

> (APB) For every `r2` there are `h_b` and `mu0` such that for `mu < mu0`, `rho − 1 ≤ r2 eps^2`, no cycle leaves
> `{H0 ≤ h_b}`.

What is known a priori (v2.0 eq. (6.1)): `U_max ≤ log M0`, `M0 = (lambda+mu)(lambda+1)/mu ≈ 2/mu` in the wedge, so
the crest energy is at most `rho Phi(log M0) ≈ 2/mu ≈ 2/eps^2`. (No comparable bound on `|p|` at negative `p` is
available: `p < lambda` is one-sided.) So the gap to close is the crest-energy window `h_b < h ≲ 2/eps^2`. It splits
in two:

- `h_b ≤ h ≤ delta^2/eps^2` (small `delta`): the slow-manifold expansion is still valid with effective small
  parameter `eps·sqrt(2h) ≤ delta·sqrt 2`, and formally the third-order drift is
  `eps^3 J0(h)(r − Rsel(h)) < 0` there (since `Rsel(h) ≥ (4/5)h > r2`). Making this rigorous needs error bounds for
  the normal form that are uniform on a **non-compact** energy range (spike of height `log h`, trough of depth `h`).
  Not done. Fenichel on a compact set does not cover it.
- `h ≥ delta^2/eps^2`: pulses with `p = O(1)` in the original time. The Hamiltonian approximation is wrong at
  leading order (the pulse map `S` of Section 5 has `S(q) − q = O(q^2)`, not small). One needs the relaxation
  analysis of Section 5 at `rho → 1`, where the singular fixed-point equation has no solution with `q ≥ delta` —
  true numerically, unproved.

I tried to get (APB) cheaply from the exact identities of Lemma 1 and from the v2.0 obstruction theorem
(neutral cycle on the fixed-`mu` path). Neither works: (c) and (f) are selection rules, satisfied to leading order by
*every* Toda orbit or pulse, and the obstruction theorem permits a large semi-stable cycle to be born away from the
Hopf branch ("out of the blue"), which is exactly what (APB) must exclude.

### 4.2a The Hopf end `r → 0` (detailed outline; same standard inputs as Theorem C; not counted as proved)

Claim C′: the conclusion of Theorem C holds for all `r ∈ (0, r2]`, i.e. for `1 < rho ≤ 1 + r2 eps^2`, `mu < mu0`.
Outline. Let `eta_h` be the half-turn map of `X_eps` on the `U`-axis (smooth through `a = 0`, `eta_h(a) = −e0 a + e1 a^2 + …`);
`Pi = eta_h∘eta_h` gives `Pi(a) − a = a·g`, `g = c0·u(a,eps,r) + a^2 v(a,eps,r)` with `c0 = e0^2 − 1 = Pi'(0) − 1`,
`u(0) = 1` (the `a^1` coefficient of `g` is a multiple of `c0`). The linearisation at `E` is exactly neutral at
`rho = 1` for every `eps` (`(eps z+1)(z^2+1)`), so `c0 = eps^3 r c^0(eps,r)`, `c^0(0,0) = pi`; and by Step 3,
`v = eps^3 v^`, `v^(0,0,0) = lim −kappa(a) Rsel(Phi(a))/a^2 = −pi/2`. Zeros of `g` with `a > 0` solve
`r = a^2 w(a,eps,r)`, `w(0,0,0) = 1/2`; the implicit function theorem gives one branch `r = r^(a,eps)`, increasing in
`a` on `(0,a1]`, so one simple attracting zero for each small `r`. On `[a1,a_b]` and `r ≤ Rsel(Phi(a1))/2` the leading
term `kappa(a)(r − Rsel(Phi(a)))` is negative and bounded away from zero. I have not written the smoothness
bookkeeping (`k_s ≥ 8`) line by line.

### 4.2b A reduction of the obstruction: only *non-attracting* large cycles must be excluded

By the pointwise criterion (brief 3(a)) the wedge needs "every cycle is attracting", not (APB). With Theorem C it is
enough to prove

> (APB-N) for `mu < mu0`, `rho − 1 ≤ r2 eps^2`: every cycle with `lambda_c ≥ 0` stays in `{H0 ≤ h_b}`.

A non-attracting cycle is strongly constrained. Put `q = eta − mu > −mu`, `qh = q/mu`, `s = (1+eta)/A`, `c = zeta/mu`.

**Facts for a cycle with `lambda_c ≥ 0` (PROVED; consequences of the v2.0 moment identity and `eta > 0`).**
`k − k_H = (A^2+mu)<q> + 2A<q^2> + <q^3>` and `<q^3> ≥ −mu<q^2>` give

    <qh^2> ≤ r/(A(2+mu)) ≤ r/2,     0 ≤ <qh> ≤ r mu/(A(A^2+mu)),     <|qh|^3> ≤ r/(A mu) + 2<qh^2>,     <c> = 1 + A<qh> + mu<qh^2> ≥ 1.

**Exact formulas on every cycle (PROVED; `variance_route.py` checks (X1),(X2) to `1e-8` on 6 cycles).**

    (X1)  U'' − (log(1+eta))' U' − sqrt(mu) qh U' + rho G(U)/s = 0            (tau-time; from orthogonality + Riccati)
    (X2)  with V = U'/s and E1 = V^2/2 + rho Phi(U):   E1' = sqrt(mu) qh V^2 + rho G(U) V (s^3 − 1)/s^2
    (X3)  rho e^U − s = A L^2 qh,   L = eps d/dtau + s      (equivalently c = L^{-1}(rho e^U), A qh = L^{-1}(c − 1))

(X2) contains no derivative of `qh`, and `|(s^3−1)/s^2| ≤ 7 mu |qh|/A` (as `s > 1/A`), `|G| ≤ Phi + sqrt(2 Phi)`, so

    |(log E1)'| ≤ |qh| [ 2 sqrt(mu) + (7 mu/A)(sqrt(2 E1) + 2 sqrt(rho)) ].

**What follows formally (SKETCH).** (i) Testing (X3) against a bump `w` around the crest interval
`J = {rho e^U > 1}` (one interval per period, v2.0 Thm 4.1) and using Cauchy–Schwarz with `<qh^2> ≤ r/2`:
`<(rho e^U − 1)_+> ≤ A sqrt(r/2) sqrt(theta + O(eps/T)) + o(1)`, `theta = |J|/T`: a non-attracting cycle cannot be
spiky. (ii) If `sqrt(mu) T` and `mu T sqrt(max E1)` are bounded, (X2) keeps `E1` within a bounded factor along the
cycle; a cycle with large `min E1` has a long deep trough and a short crest (`theta → 0`, `<(rho e^U−1)_+> ≥ 0.6`),
contradicting (i). So (APB-N) would follow from two a-priori estimates for arbitrary cycles in the wedge:

    (P-T)  T ≤ C1/mu in the original time,         (P-E)  p ≥ −C2 in the original time (so E1 ≤ C/mu).

Both hold for the relaxation cycles of Section 5 (`mu T → Sigma`, `p ≥ −S(q)`), neither is proved for arbitrary
cycles. Data: on true (attracting) cycles in the wedge `<qh^2>` is 1.7–2.0 times the ceiling
(`variance_route.out`, `mu ≤ 0.01`) and `<qh^2>/<G^2> = 0.97–0.99` at `mu = 1e-3`, i.e. `qh ≈ e^U − 1`; for a Toda-like
cycle the ceiling forces `<P^2>_h ≤ r/2`, hence `h ≤ 3r/4` by Theorem A(i). The route does not help at fixed
`rho` (`r → ∞`): at `mu = 0.05, rho = 3` the true cycle has `<qh^2> = 14.6` below the ceiling 20.5.

### 4.3 Uniformity (the three questions of job 1(c))

| regime | what holds | status |
|---|---|---|
| `rho` in a compact subset of `(1,∞)`, `mu → 0` | **Not** a Hamiltonian regime (`r = (rho−1)/eps^2 → ∞`, `h* ≈ (5/4) r → ∞`). Relaxation oscillation, Section 5. Theorem C says nothing. | numerical only |
| near `rho = 1`: `r ∈ [r1,r2] ⊂ (0,∞)` | Theorem C, uniform in `r` on compacts. Width of this wedge in `rho` is `≈ r2·mu`. | proved (local in phase space) |
| `r → 0` (Hopf end) | Claim C′ (§4.2a): uniqueness in `{H0 ≤ h_b}` for all `r ∈ (0,r2]`, `mu < mu0`. | detailed outline |
| `r → ∞` with `r eps^2 → 0` (matching) | formal: `h* ≈ (5/4) r`, `log m_t ≈ −(16√2/15) eps^3 h*^{3/2} = −2.108 (rho−1)^{3/2}`; agrees with the relaxation map at `rho = 1.02` (−0.00566 vs −0.00596). Not uniform in Theorem C. | numerical/formal |
| `rho → ∞` | outside both regimes (large gain); the relaxation map has `q* → lambda^−`, multiplier → 0 (0.00145 at `rho = 8`). | numerical only |

So a correct statement of "the small-`mu` theorem that would follow" is: *Theorem C + Claim C′ + (APB-N) give
uniqueness for `mu < mu0(r2)` in the wedge `1 < rho ≤ 1 + r2 eps^2`* (and (APB-N) reduces to the two a-priori
estimates (P-T), (P-E) plus the estimates sketched in §4.2b); and separately, *for `rho` in a compact
subset of `(1,∞)`, uniqueness for small `mu` would follow from monotonicity of `Lambda` (Section 5) plus a
singular-perturbation theorem of relaxation (entry–exit) type.* Neither is uniform as `rho → ∞`.

## 5. The relaxation regime (`rho` fixed, `mu → 0`) — numerical, with one exact solution

Original time, `lambda = k/mu = rho(1+mu)`, `y = kZ`: `U' = p`, `p' = mu(lambda(1−Z) − p)`, `Z' = e^U − Z`.

- *Quiescent phase* (`U < 0`, `Z ≈ 0`): slow time `s = mu t`, `V = mu U`: `V' = p`, `p' = lambda − p`. From
  `p = −q_− < 0` at `V = 0` the return to `V = 0` happens at `s = Sigma` with `lambda Sigma = (lambda+q_−)(1 − e^{−Sigma})`,
  with `p_+ = lambda Sigma − q_− =: G_lambda(q_−)`; depth `mu U_min = lambda log(1 + q_−/lambda) − q_−`.
- *Pulse* (`U = log(1/mu) + u`, `z = mu Z`; put `v = u + log lambda`, `zh = lambda z`): `v' = p`, `p' = −zh`,
  `zh' = e^v − zh`, i.e. `v''' + v'' + e^v = 0`, from `(v,p,zh) = (−∞, q, 0)` to `(−∞, −S(q), 0)`. `S` is universal.
  Exact: multiplying `p'' + p' = −e^v` by `p` and integrating gives `S(q)^2 − q^2 = 2 int zh^2 dt > 0` (the pulse
  amplifies `|p|`); `S(q) = q + (2/3)q^2 + O(q^3)` (Liouville soliton `e^v = (q^2/2) sech^2(qt/2)` plus first correction);
  **exact pulse** `p = −2 + 3/(1+e^t)`, `e^v = 6 e^t/(1+e^t)^3`, i.e. `S(1) = 2` (found numerically, then verified by
  substitution).
- *Singular return map* `F_lambda = G_lambda ∘ S` on `(0,lambda)`; `F_lambda(q) − q = (2/3)(1 − 1/lambda) q^2 + O(q^3)`,
  the relaxation counterpart of Prop. 2.

`relax_map.py` → `relax_map.out` (~100 s): on a 60-point sign scan of `F_lambda(q) − q` over `(1e-3, lambda)` there
is exactly one sign change for each `rho ∈ {1.02, 1.1, 1.5, 3, 8}` (at `rho = 30` the scan resolves none: the root
is within `1e-6·lambda` of `lambda`):

| rho | q* | multiplier F'(q*) | Sigma (= lim mu·T) | lim mu·U_min | true cycle at mu = 1e-3 (m_t, mu·T, mu·U_min) |
|---|---|---|---|---|---|
| 1.02 | 0.2127 | 0.99436 | 0.4508 | −0.0258 | – |
| 1.1 | 0.4717 | 0.94764 | 1.0322 | −0.1444 | – |
| 1.5 | 1.1502 | 0.67983 | 2.4387 | −1.0336 | 0.67851, 2.4470, −1.0325 |
| 3 | 2.9177 | 0.18665 | 5.2615 | −7.8699 | 0.18607, 5.2685, −7.8823 |
| 8 | 7.9994 | 0.00145 | 11.998 | −68.105 | – |

*Amplitude parametrisation of the singular limit.* The fixed-point equation is `S(q) = lambda g(q/lambda)`,
`g(a) = s − a` with `(e^s − 1)/s = 1/(1−a)`. Since `g(a)/a` increases from 1 to ∞ (checked on a 500-point grid, not
proved) and `S(q) > q`, for each `q > 0` there is exactly one `lambda = Lambda(q) > q`. **Uniqueness in the singular limit for every `rho > 1` ⟺ `Lambda` is strictly
increasing.** `lambda_of_q.py` → `lambda_of_q.out`: `Lambda` is increasing on all 33 grid points `q ∈ [0.02, 12]`
(`Lambda(0.02) = 1.00016`, `Lambda(1) = 1.39592`, `Lambda(12) = 12.000007`); small-`q` behaviour `Lambda − 1 ≈ 0.4 q^2`
matches the Toda large-`h` law `h* ≈ (5/4) r`. Status: numerically supported; not proved. It is a one-variable
statement about a universal function and is a natural target for a 1-D interval proof.

## 6. Other observations (not proved; for the other lanes)

- *Dulac hint from the Toda regime.* In the wedge the averaged divergence of `B·X_eps` with `B = b(Ht)` is
  `eps^3 (b M3)'`; the choice `b = 1/J0` (reciprocal of the area enclosed by the level through the point) gives
  `(r − Rsel)' = −Rsel' ≤ −4/5`. So to leading order a Dulac function is `1/J0(Ht)` (times a near-identity
  correction), and `1/J0 ~ 1/(pi(U^2+P^2))` at `E`, consistent with brief 3(d).
- *Amplitude formulation.* In both limits the branch is a graph over the amplitude: `r = Rsel(h)` (proved monotone)
  and `rho = Lambda(q)` (numerically monotone). The statement "(AMP): for each `mu` and each crest amplitude `a`
  there is exactly one `rho` with a cycle of that amplitude, and it increases with `a`" implies the conjecture and is
  the form in which both limits are clean. It is stronger than the conjecture and untested beyond the primary branch.

## 7. Files and commands

| file | purpose | run time |
|---|---|---|
| `toda.py` / `toda.out` | Toda level integrals, identities of Prop. 2, kurtosis, `Rsel` | 2 s |
| `normal_form_check.py` / `.out` | sympy: slow-manifold series, modified energy, `F3`, identity (c) | 5 s |
| `ampcycle.py` / `ampcycle.out` | amplitude-fixed cycles at `mu = 3e-2 … 1e-4` vs third-order prediction | 7 s |
| `fixed_rho.py` / `fixed_rho.out` | fixed `rho`, `mu → 0`: energy diverges; identity (c) on cycles | 14 s |
| `exact_identity.py` / `.out` | Lemma 1 (b),(c),(f) on 8 true cycles | 20 s |
| `relax_map.py` / `relax_map.out` | singular relaxation map, fixed points, multipliers | 100 s |
| `variance_route.py` / `.out` | exact formulas (X1),(X2) and the variance ceiling on true cycles | 15 s |
| `ws_check.py` / `.out` | `W^s(E)` stays in two octants (PROGRAM (L5)) | 3 s |
| `lambda_of_q.py` / `.out` | `Lambda(q)` monotonicity scan | 60 s |
| `crest_data.py`, `core_cost.py`, `core_total.py`, `tm_probe*.jl` and outputs | sizing probes for `PROGRAM.md` | see there |

## 8. What I would try next

1. (APB-N) by the variance-window route of §4.2b: prove (P-T) and (P-E) for arbitrary cycles in the wedge, then make
   the spike inequality (i) and the trough/crest timing rigorous. This avoids uniform third-order averaging on an
   unbounded energy range, and it is the shortest path I see to a *global* small-`mu` theorem near the Hopf curve.
   (Fallback: uniform averaging for `h_b ≤ h ≤ delta^2/eps^2`, where only the sign of the drift is needed, with
   margin `Rsel(h) − r ≥ (4/5)h − r2`.)
2. `Lambda' > 0` by interval arithmetic on the pulse equation for `q` in a compact range, plus the two asymptotic
   ends (`q → 0`: Liouville soliton with `S = q + (2/3)q^2 + c3 q^3`; `q → ∞`: `u''' ≈ −e^u`).
3. A singular-perturbation theorem for the relaxation cycle with explicit `O(mu log 1/mu)` errors (the map is a
   regular perturbation of `F_lambda` once `q ≥ delta`); this is also what the small-`mu` wedge of `PROGRAM.md` needs.
4. Write out the Hopf-end blow-up (`r → 0`) so that Theorem C covers `r ∈ (0, r2]`.
