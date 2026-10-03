# Lane `freeform` — notes (2026-10-03)

Status: **partial**. No proof of the conjecture, no counterevidence. A handful of small rigorous
lemmas (full proofs below), one numerically supported structural lemma (the 1-D return map on the
attractor is concave and sublinear in the crest coordinate `U`), and a log of routes that died, each
with the reason. Python: `/Users/admin/.venvs/analysis-py313/bin/python`. All scripts are in this directory.

Notation as in BRIEF: `x = k e^U`, `q = lambda - p`, `A = 1+mu`, `gamma = 1 + <eta>`, `delta = k - k_H`,
`f(s) = s(s+1)(s+1-mu)`, `F(s) = s(s-1)(s-mu)` (so `f(s) = F(s+1)`), `m := gamma(gamma-1)(gamma-mu)`.
Stability (`lambda_c <= 0`) is equivalent to `gamma <= A`, i.e. `m <= k_H`.

---

## 1. Rigorous results (each step written out)

### L1. Hölder sandwich (improves brief 3(g))

**Statement.** On every nonconstant cycle,

    < x^{1/3} >^3  <  gamma (gamma-1)(gamma-mu)  <  <x> = k .

Consequently every attracting or neutral cycle satisfies `< e^{U/3} > < rho^{-1/3}` (strict).
Since `exp<U> <= <e^{U/3}>^3` (Jensen), this implies and strengthens brief 3(g).

**Proof.** The three functions `a = eta`, `b = zeta/eta`, `c = x/zeta` are positive and periodic and `abc = x`.
From the Riccati system, `(log eta)' = zeta/eta + mu - 1 - eta` and `(log zeta)' = x/zeta - 1 - eta`; averaging,
`<b> = 1 + <eta> - mu = gamma - mu`, `<c> = gamma`, and `<a> = gamma - 1`. Generalised Hölder with exponents
(3,3,3): `<(abc)^{1/3}> <= (<a><b><c>)^{1/3}`. Equality would force `zeta = c1 eta^2`, `x = c2 eta^3`, and then the
first Riccati equation becomes the autonomous scalar ODE `eta' = (c1-1) eta^2 + (mu-1) eta`, whose periodic
solutions are constant; but `eta` is nonconstant on a nonconstant cycle (paper, after eq. (2.9)). So the left
inequality is strict. The right inequality is the paper's strict Jensen bound `f(<eta>) < k` (`f(<eta>) = m`).
`m` is increasing in `gamma` on `gamma > max(1,mu)`, and `m = k_H` iff `gamma = A`, which gives the corollary. ∎

Reading: the "effective mean" `m = <x â>/<â>` (brief 3(f)) lies strictly between the power mean of order 1/3 and
the arithmetic mean of `x`, and the conjecture says it lies below the `p^2`-weighted mean `<x p^2>/<p^2> = k_H`.
Numerical check (`verify_identities.py`, last block), e.g. `mu=1, rho=2`: `0.3195 < 1.0841 < 2 (=k_H) < 4 (=k)`.

### L2. Exact Lotka–Volterra form

With `r1 = x/y`, `r2 = y/q`, `q = lambda - p` (all positive on the chemical domain):

    r1' = r1 (lambda + 1 - r1 - q),   r2' = r2 (mu - 1 + r1 - r2),   q' = q (r2 - mu).

Proof: direct differentiation (sympy: `verify_identities.py`, first three lines "OK"). The interior equilibrium
is `(1, mu, lambda)`; the classical LV time-average property reproduces `<x/y> = 1`, `<y/q> = mu`, `<q> = lambda`.
The interaction matrix `[[-1,0,-1],[1,-1,0],[0,1,0]]` is a cyclic chain, not a competitive LV system, so the
carrying-simplex literature (Zeeman classes, Hofbauer–So) does not apply. WH95 describe the model as a 3-D
S-system, so this form is very likely folklore; it is recorded because it gives L5 below in one line.

### L3. Second-order Riccati (Cole–Hopf / Burgers-hierarchy) form of the stability problem

Let `𝓛[w] := w'' + (3w + 2 - mu) w' + f(w)`. For any positive `C^3` function `phi` and constant `c`,

    𝓛[ phi'/phi + c - 1 ] = (D+c)(D+c-1)(D+c-mu) phi / phi .

Proof: put `H3 = e^{ct} phi`, `H2 = H3' - H3`, `H1 = H2' - mu H2`; then `w = H2/H3`, `z = H1/H3` satisfy
`w' = z + (mu-1)w - w^2` identically and `z' + (1+w) z = H1'/H3`; substituting `z = w' + w^2 + (1-mu) w` gives
`𝓛[w] = H1'/H3`. ∎  Special cases (sympy-checked):

* `𝓛[eta] = x`  (the adjoint ratio solves a *forced Liénard equation* with positive damping `3 eta + 2 - mu`);
* `𝓛[y/q] = x (A+p)/(lambda-p)`  (the chemical ratio solves the same equation with a different forcing);
* in the plane `(eta, sigma)`, `sigma = zeta + eta + eta^2/2`:  `eta' = sigma - G(eta)`, `sigma' = x - f(eta)`,
  `G(s) = 3s^2/2 + (2-mu)s` — a Liénard plane whose first equation is explicit;
* time-reversed twin: `beta := -a'/a` for the strong-stable Floquet solution `a` satisfies
  `beta'' - (3 beta - A) beta' + F(beta) = x`, `<beta> = gamma`, while `1+eta` satisfies the same with `+(3·-A)`.

So: `gamma[x] <= A` iff there is a periodic mean-zero `w` with
`x <= w'' + (3w + 2A) w' + (w+A)(w+1)(w+mu)` pointwise. This is brief 3(e) in logarithmic form (`phi = exp ∫w`);
positivity of `phi` alone suffices (the Green kernel `G3` of `(D+A)(D+1)(D+mu)` on periodic functions is positive, so
`phi >= G3*(x phi)` and Collatz–Wielandt gives `spr(G3∘M_x) <= 1`; the equivalence of `spr <= 1` with a nonpositive
Perron–Floquet exponent is the standard periodic next-generation threshold theorem, imported, not reproved here), the
side conditions `(D+mu)phi > 0`, `(D+1)(D+mu)phi > 0` in 3(e) are not needed. Not a new criterion.

Backward form (C−), also equivalent: a positive periodic `phi` with `(A-D)(1-D)(mu-D) phi >= x phi`.
Explicit values: `Lc^-[1] = k_H - x`,  `Lc^-[x] = x[(A-p)(3p' - p(A-p)) - delta]` (no `x` inside the bracket),
`Lc[1/x]·x = -(A-3p)p' - p(A-p)^2 - delta`. None is one-signed on a cycle.

### L4. The Gauss map of the cycle family is injective

**Statement.** If `z1 != z2` are points of cycles (or `E`) at the same parameters, then `(eta1,zeta1) != (eta2,zeta2)`.
Hence each cycle is a simple closed curve in the `(eta, zeta)` plane, distinct cycles have disjoint images,
and in the `(eta, sigma)` plane each is traversed clockwise.

**Proof.** Theorem 3.1 with base `z1`: `n1·(z2 - z1) < 0`; with base `z2`: `n2·(z1 - z2) < 0`, `n_i = (zeta_i, -eta_i, -1)`.
If `n1 = n2` the two inequalities contradict each other. Orientation: `eta' = sigma - G(eta)` and `G` is increasing on
`eta > max(0,mu-1)`; if the topmost point `(eta_t, sigma_max)` had `sigma_max < G(eta_t)` then, at the point where `eta`
is maximal, `sigma = G(eta_max) >= G(eta_t) > sigma_max`, impossible. So `eta' >= 0` at the top. ∎

The area inequalities this yields (`∮ b(eta) dsigma < 0` for increasing `b`) are *trivial* identities
(`= -∮ b'(eta) eta'^2 dt`), see §3.6.

### L5. No certificate can keep the chemical ratios `x : y : q`

**Statement.** There is no positive periodic scalar `M(t)` and constants `c1,c2,c3 > 0` such that
`G = M·(c1 x, c2 y, c3 q)` satisfies the vector supersolution inequalities
`G3' + mu G3 >= G2`, `G2' + G2 >= G1`, `G1' + A G1 >= x G3` along a cycle with `k > k_H`.

**Proof.** Divide the i-th inequality by `G_i` and average; `<(log G_i)'> = 0` and L2 give
`mu >= (c2/c3) mu`, `1 >= c1/c2`, `A >= (c3/c1) lambda`, i.e. `c3/c1 <= 1/rho < 1` and `c1 <= c2 <= c3`. ∎

In particular all monomial Dulac/supersolution ansätze `x^a y^b q^c (c1 x, c2 y, c3 q)` are dead for every cycle,
by the LV averaging property — any working certificate must distort the ratios `x:y:q` non-trivially.
(The set of vector supersolutions is closed under sums and under weighted geometric means — Hölder — so partial
certificates can be combined; the two obvious ones, `(mu,mu,1)` for input `k_H` and `(x,y,q)` for input
`x(A+p)/q`, do not combine to cover the crest: one would need `x/k_H <= ((A+p)/q)^{theta/(1-theta)}` at `p=0`.)

### L6. Massera identity (exact; reduces the planar question to one sign)

Let `omega1 := p^2 - U p'`. On every cycle (sympy-checked using tangent orthogonality)

    omega1' = (eta - mu) omega1 + p K .

So wherever `omega1 > 0` on the whole cycle (projected cycle star-shaped about the origin; true in every computed
case, not proved), `lambda_c = -< p K / omega1 >`, and the conjecture is `∮ p K/omega1 dt > 0`: `K > 0` and
`omega1 > 0`, but `p` changes sign, so the radial scaling contracts on the half `p > 0` and expands on `p < 0`.
In polar form `d log r / d(-theta) = sin(theta)(cos(theta) - mu sin(theta) - m)/Omega`, `m = (h-k)/r` increasing in `r`
(convexity of `h`), `Omega = -theta'`, and `∂/∂m` of the right side is `-sin(theta)/Omega^2`: exact sign `-sgn(p)`.
Any other transversal field is `g·(U,p)` modulo the flow direction and changes the integrand by `(log g)'`, i.e.
this is the Dulac problem of brief 3(d) with `B = 1/omega1`; near Hopf the integrand is cubic to leading order and
its mean is `O(h^2)`, so no pointwise sign is available without normal-form information at `E`.

---

## 2. Numerically supported (NOT proved)

### N1. The crest return map on the attractor is concave and sublinear in `U`

Let `P1` be the first-return map of the planar flow on the attractor to `{p = 0, U > 0}`, in the coordinate
`u = U_crest`. Along `W^u(E)` (which fills the disc inside the innermost cycle):

* `R(u) = P1(u)/u` is strictly decreasing from `exp(2π Re λ_u/Im λ_u)` to 1 — every case tested;
* `P1'(u) = (p'(0)/p'(τ)) exp ∫(eta - mu)` is decreasing (concave `P1`) — every case tested except two
  flagged samples that are noise (see below);
* equivalently (L6), `∫_{crest}^{next crest} p K/omega1 dt > 0` on every turn (min over turns `3e-5 … 5e-4`, it is
  `O(u^2)` near `E`).

Tests: `e1_return1d.py` (mu ∈ {0.05,0.3,1}, rho ∈ {1.05,1.3,2,4}); `e1c_ratio.py` (mu ∈ {0.02,0.1,0.2,0.5,1},
rho up to 10, two starting angles); `e1b_concavity.py` (exact `P1'` from the pull-back adjoint, mu ∈ {0.1,1},
rho ∈ {1.2,2,4}). Caveats: (i) `P1` is sampled only on the crest sequence of one or two orbits per parameter
(7–1000 points); (ii) one "violation" in `e1b` is at the first sample `u ≈ 1e-3` where the pull-back `eta` had not
converged (consistency check `5.6e-4`); the two "secant slope" failures in `e1c` occur where successive crests differ
by `~1e-10` (cancellation). In the coordinates `e^u - 1` or `y - k` the ratio is *not* monotone, so `U` is the
right coordinate. If N1 held on the whole attractor it would give uniqueness at once (`P1(0)=0`, `P1(u)/u`
strictly decreasing ⇒ at most one positive fixed point). I found no proof: the backward flow is monotone and
convex (brief 3(j)), but convex combinations of attractor points lie in the forward-invariant epigraph, on the
wrong side for the convexity inequality to bound the attractor's own return map (§3.3).

---

## 3. Routes tried and why each died

1. **Scaling/commutator (Massera) in 3-D.** For `S` with `[F,S]` explicit: `S = z - E` gives `[F,S] = (0,0,-k g(U))`;
   chemical scaling `S = (1, p - lambda, y)` gives `[F,S] = (q,0,0)`; `S = z''` gives `(0,0,x p^2)`. In each case the
   compound source `conj(F ∧ [F,S])` has components proportional to `p`, `p'` or `y'`, which change sign. Any
   wedge `F ∧ u` is a compound vector with source `F ∧ (u' - Ju)`; one-signed sources would need `u' - Ju` to rotate with
   `F`. Reduces to brief 3(e).
2. **Lyapunov–Schmidt / fold type.** `N(k)` := number of cycles. If cycles can never be *born* as `k` increases
   (local count near every neutral cycle non-increasing in `k`), then compactness + Hopf collar + existence give
   `N ≡ 1`; if they can never *die*, uniqueness for large `k` propagates down to all `k`. The fold coefficients
   are `b = ∫Ψ·∂_k F` and `a ∝ ∫Ψ_3 x v_1^2` with `Ψ = z' × v_s` the weak adjoint, whose components change sign
   (`Ψ_3 ∝ -(b p + a p')`), so neither sign is available. Moreover N1 (concave `P1`, return map increasing in `k`)
   says the generic sign is the *birth* type, so "death-only" is the wrong target and the route cannot close
   without large-`k` uniqueness plus a proof of the sign of `a`.
3. **Concavity of the 2-D crest map (Krasnoselskii).** The crest map is monotone on `{y>k, x>=y}`, and a concave
   monotone map fixing `E` has a unique positive fixed point. But in the linear coordinates `(U, y-k)` concavity
   forces `M ξ <= ξ` at a fixed point `ξ`, and by the paper's Theorem 7.1 formulas the second component is
   `(m_t - 1) r_e - (r_v - r_e)(1 - M_22)` with `r_v - r_e = -(mu^2+4mu+1) h/(2 A D_mu) < 0`: positive at `O(h)`.
   So 2-D concavity fails near Hopf in these coordinates; it can hold only in coordinates that straighten the
   (unknown) attractor curve to second order.
4. **Cross-parameter comparison 3(h) as a rotated family.** For a projection along a direction `nu`, the planar
   fields at `k1 < k2` have cross product `(k2-k1)(nu_3 p - nu_1 y') + s·det[nu, F, J nu]`; since `<p> = <y'> = 0` the first
   term changes sign for every constant `nu`. No linear projection makes the family rotated in Duff's sense; (h) is
   one-directional and is consistent with both birth and death folds (second-order expansion gives no constraint
   for `k2 < k1`).
5. **Spectral-radius bounds without a test function.** `spr(G3 ∘ M_x) < 1` is the conjecture (`G3` the kernel of
   `(D+A)(D+1)(D+mu)`). Trace bounds fail: at Hopf `tr N ≈ 1 - 2mu/(A^2+mu) < 1 = spr`, the other eigenvalues are
   complex with negative real parts. Kingman log-convexity gives only lower bounds on `gamma` (L1). Upper bounds
   on a Perron root need a test vector; everything pointwise collapses to 3(e).
6. **Moment/area identities.** Simple-closed-curve area inequalities in the `(eta,sigma)` plane reduce to
   `∮ b'(eta) eta'^2 > 0`. All Riccati moment identities hold for arbitrary periodic input `x(t)` with `<x> = k`,
   for which `<eta> < mu` is false (e.g. `x ≡ k`), so they cannot decide the sign.
7. **Darboux/Bäcklund commutation** of the factorisation `(D+1+eta)(D^2+(mu-eta)D+zeta)`: the commuted operator
   is `D^3 + A D^2 + (mu + 3 eta') D + …`, not of the same form. No self-map.
8. **Dual (Legendre) Hamilton–Jacobi form.** `K = k + h^*(zeta,-eta)` satisfies
   `(1+eta)(K - zeta K_zeta) + (zeta - g(eta)) K_eta + k(e^{K_zeta} - 1) = 0`, with Lagrangian the relative entropy
   `x log(x/k) - x + k`; `K` is the discounted relative entropy of the past of `x`. Pretty, but the divergence in the
   dual plane is `x U_zeta - (3 eta + 2 - mu)` and `<x U_zeta> = 2 gamma + 2 lambda_c` is again a tautology.

9. **Mirror comparison in `U` (tested, refuted as a pointwise statement).** Parametrising both half-cycles by `U`,
   `-lambda_c T = ∫_{Umin}^{Umax} [K_+/omega1_+ - K_-/omega1_-] dU` (upper arc `p>0` minus lower arc `p<0` at equal `U`).
   `e4_mirror.py`: the identity reproduces `-lambda_c T` to 5 digits (mu=0.3,1,0.2,3; rho=1.05…2), and `omega1 > 0`
   with `lambda_c = -<pK/omega1>` to 6 digits on nine cycles (`e3_starshape.py`). But the bracket is **negative** on
   essentially the whole range `0 < U < U_max` (12–40 % of the `U`-range) and positive for `U < 0`: just after the
   crest `K/omega1` is larger than just before it. So there is no pointwise-in-`U` Massera/Zhang argument; the
   positive balance comes from the long excursion `U < 0`. (Two rows of that script, mu=0.05/rho=1.3 and
   mu=1/rho=5, print a wrong integral: an interpolation artefact from a non-monotone `U` sample at `p≈0`; ignore them.)

## 4. Sharpest obstruction

Every route ends at the same place: the weak transverse direction is the *second* direction of the cooperative
(backward) structure, so its adjoint has no sign, and all sign-definite objects (supporting planes, `K`, `W`,
Hölder/Kingman) control only the Perron (strong) direction or give *lower* bounds on `gamma`. An upper bound
`gamma < A` needs a positive test function singular or degenerate at `E` that matches the Hopf normal form to
second order; near Hopf the required inequality is an `O(h^2)` effect hidden under `O(h)` sign-changing terms
(L6), and for large amplitude it is a different mechanism (`x ≈ 0` most of the period). No listed or new
structural identity bridges the two.

## 5. What I would try next

* Prove N1 in the form `∫_turn pK/omega1 dt > 0`. The equal-`U` mirror comparison is refuted (§3.9); a comparison
  at mirror polar angles `±theta` (using `h(U,p) < h(U,-p)` for `p > 0`, so `m_up < m_low`, and
  `Omega_up - Omega_low = cos(theta)(2 mu sin(theta) + m_up - m_low)`) is untested and would have to be non-pointwise too.
* Use L1 quantitatively: a neutral cycle needs `<e^{U/3}> < rho^{-1/3}`, which with `<e^U> = 1` forces an explicit
  lower bound on the `U`-range; combine with the paper's variance window to shrink the admissible neutral set.
* Hand L5 and the (C−) operator values to the supersolution lane; the crest is the binding point for every
  explicit candidate tried (`1`, `q`, `x`, `1/x`, `x^{-a}`), the trough for `1/x`.

Web searches used: 2 of 5 (3-D LV limit-cycle literature; WH as S-system). Nothing found that settles the
cyclic-chain LV case.
