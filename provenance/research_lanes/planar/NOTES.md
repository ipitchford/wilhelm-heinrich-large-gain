# Lane `planar` — NOTES (2026-10-03)

Everything for this lane is in this directory. Python: `/Users/admin/.venvs/analysis-py313/bin/python <script>`.
All numerics are ordinary floating point (DOP853, rtol 1e-11..1e-12); none is an interval certificate.

## 0. Outcome in one paragraph

**Job 1 is done:** the attractor-graph theorem of brief item 3(b) is proved in full in
`attractor_theorem.tex` (compiled: `attractor_theorem.pdf`, 13 pages). No gap was found in the sketch; three points
of the sketch needed sharpening (listed in §1.3). The proof is self-contained apart from classical theorems
(Birkhoff–Hopf contraction, planar Poincaré–Bendixson for Lipschitz fields, stable-manifold theorems, Jordan/
Schoenflies, Whitney C^1 extension). It does **not** use Smith's competitive Poincaré–Bendixson theorem, Smith's
finiteness of cycles, or anything from v2.0 (the two v2.0 facts needed are re-derived in two short paragraphs);
permanence is re-proved with explicit constants.
**Job 2 is not solved.** No uniqueness proof. For each of (a)–(d) there is either a proved conditional statement
with a numerically supported missing inequality, or a numerical refutation of the natural pointwise candidate.
Two numerical leads came out of it (§3.5): the crest return map is concave in the coordinate U (14/14 parameter
pairs), and a low-degree *state-explicit* compound supersolution exists on the sampled attractor at some
parameters (LP), but not robustly across parameters in the basis tried.

## 1. Job 1 — the theorem (status: PROVED, full proof in `attractor_theorem.tex`)

### 1.1 Statement (verbatim from the tex, Theorem 1.1)

Let mu>0, k>0, system U'=p, p'=k−mu p−y, y'=k e^U−y, D_chem={y>0, p<lambda}, pi(U,p,y)=(U,p).

- (A1) D_chem is forward invariant; the semiflow has a compact connected invariant global attractor 𝒜 ⊂ D_chem
  attracting every compact subset; 𝒜 = union of all log-bounded entire orbits of the system in R^3.
- (A2) Along every log-bounded entire orbit the Riccati system eta'=zeta+(mu−1)eta−eta², zeta'=x−(1+eta)zeta has
  exactly one solution on R with eta≥0, zeta≥0. It satisfies eta_−≤eta≤eta_+, zeta_−≤zeta≤zeta_+ with the explicit
  constants of v2.0 Lemma 2.3 evaluated at x_∓ = min/max of k e^U over 𝒜 (so eta > max(0,mu−1) uniformly), attracts
  all nonnegative data with |log eta~ − log eta|, |log zeta~ − log zeta| ≤ Delta_0 tau_0^{floor(t−s)−1}
  (Delta_0 = 2 log(M_0/m_0), tau_0 = tanh(Delta_0/4), m_0/M_0 = min/max entry of exp(B_∓)), and (eta,zeta) is a
  continuous flow-equivariant function on 𝒜. On cycles it is the Perron–Floquet pair of v2.0; at E it is (eta_E,zeta_E).
- (A3) For all a≠b in 𝒜:  y_b − y_a > zeta(a)(U_b−U_a) − eta(a)(p_b−p_a).
- (A4) pi|𝒜 is a homeomorphism onto D=pi(𝒜); H(u,v)=max_{a∈𝒜}{y_a+zeta(a)(u−U_a)−eta(a)(v−p_a)} is convex and globally
  Lipschitz on R²; 𝒜 = graph of h=H|_D; for q,q' in D: 0 ≤ h(q')−h(q)−g(q)·(q'−q) ≤ (g(q')−g(q))·(q'−q) with
  g=(zeta,−eta) continuous on D; hence h is C^1 on int D with grad h=(zeta,−eta), and (h,g) is a Whitney C^1 jet on D.
- (A5) F(u,v)=(v, k−mu v−H(u,v)) is globally Lipschitz, its flow leaves D invariant and is conjugate by pi to the flow
  on 𝒜; zeta v − eta(k−mu v−h) = k e^u − h on all of D (generalised tangent orthogonality n·f=0 on 𝒜); F is C^1 on
  int D with div F = eta − mu; K>0 and W>0 on 𝒜\{E} with K'+(1+eta)K = k(Ue^U−e^U+1), W'+(1+eta)W = x p².
- For k>k_H:
- (A6) 𝒜 ∩ (E+C_±) = {E}; no log-bounded entire orbit other than E converges to E forwards; W^u(E) ⊂ 𝒜;
  (0,0) ∈ int D; every orbit in 𝒜\{E} has a cycle as omega-limit set, and {E} or a cycle as alpha-limit set, and no
  non-periodic orbit has the same cycle as alpha- and omega-limit; cycles exist.
- (A7) The set of crests is compact in (0,∞); D is exactly the closed Jordan domain bounded by the projection of
  the outermost cycle; 𝒜 is a closed topological disc with boundary Gamma_out.
- (A8) (a) for cycles in int D the planar Floquet multiplier of F equals m_t = exp∮(eta−mu); for Gamma_out the same
  for a Whitney extension; (b) m_t>1 ⇒ no orbit of 𝒜\Gamma has omega-limit Gamma; m_t<1 ⇒ none has alpha-limit Gamma;
  (c) lambda_c(Gamma_in) ≤ 0 and lambda_c(Gamma_out) ≤ 0; (d) all cycles lambda_c<0 ⇒ exactly one cycle
  (⇔ non-uniqueness ⇒ some cycle has lambda_c ≥ 0); (e) all cycles hyperbolic ⇒ finitely many, odd number, alternating
  −,+,…,−; (f) Dulac form: Lambda on the union of cycles with D Lambda + eta − mu ≤ 0, not a.e. zero on any cycle ⇒
  unique hyperbolic attracting cycle; the reverse sign is impossible.

Rigorous by-products (same file, §8): Proposition 8.1 (crest return map: R_1' = g~(u) e^{iota(u)} / g~(R_1 u);
in the energy coordinate e = ∫_0^u (h(s,0)−k) ds the return map has derivative exactly e^{iota}; uniqueness if
R_1(u)/u is strictly decreasing or iota strictly decreasing — CONDITIONAL), and Remark 8.2 (divergence of the reduced
flow in the planes (U,y), (p,y), (y,m), (U,s): zeta/eta−1, x/zeta−A, x/(zeta+mu eta)−1, (zeta−mu)/(1+eta); the
component functions are restrictions of convex/concave functions).

### 1.2 Imported inputs (all classical; listed in the tex §1.1)

Birkhoff–Hopf theorem; Poincaré–Bendixson and transversal lemmas for planar locally Lipschitz fields (Hartman VII.4);
stable/unstable manifold theorems and Lyapunov stability of attracting cycles (Hartman IX); Jordan, Schoenflies,
winding number; Whitney C^1 extension (only for (A8)(a) at the outermost cycle).
Published model results (Smith 2012 Thms 1.1, 1.3, 1.4 — checked against the arXiv PDF 1109.5204) are consistent with
but NOT used by the proof: permanence is re-proved (Lemma 2.1, explicit box), convergence to cycles on 𝒜 follows from
classical planar PB, finiteness of cycles is not needed (compactness of the crest set is proved directly).
For k ≤ k_H the remark "𝒜={E}" does cite Smith 1.1(b) / Boros–Hofbauer Thm 13.

### 1.3 Where the brief's sketch needed sharpening (no errors, three precisions)

1. "convex C^1 function over the closed disc": D need not be convex; the correct statement is h = H|_D with H convex
   on R², h C^1 on int D and a Whitney C^1 jet on D (one-sided at the boundary cycle). The key inequality is the
   two-sided bound 0 ≤ h(q')−h(q)−g(q)(q'−q) ≤ (g(q')−g(q))(q'−q), which turns continuity of (eta,zeta) into C^1.
2. "D is the disc bounded by the outermost cycle" has two halves: D ⊂ closed interior (planar PB + monotone
   transversal lemma) and closed interior ⊂ D. The second half is not a planar statement; it is proved by a degree
   argument (flow a spanning disc of Gamma_out forward, use attraction of compact sets). No shape theory needed.
3. Brief (d) says a Dulac weight "must behave like r^{-2}" at E. Not so: B = Q^{-c} with Q the invariant quadratic form
   and any c>1 gives div(BF)<0 near E. Only weights built from log K, log W with unit total weight are forced to r^{-2}.
   There is no local obstruction at E above Hopf; the difficulty is global.

### 1.4 Numerical checks of Job 1 (all pass; commands)

| script | what | result |
|---|---|---|
| `job1_checks.py` | pull-back uniqueness; O=0 on W^u(E); supporting planes for 9.0e6 pairs; bounds; cone exclusion; K,W,Lambda_E>0 | gap 360 → ≤5e-9 after t=25; rel O ≤ 7e-11; 0 violations in 3 cases, 1 pair at +3e-10 (near-coincident points, rounding) in (mu,rho)=(1,1.1); bounds hold; no point in E+C_± |
| `job1_section.py` | sharp section test: crest points of 7 different W^u(E) spirals merged and sorted; convexity bracket zeta_i ≤ slope ≤ zeta_{i+1}; same on {U=0,p>0} with −eta; 6 parameter pairs | y strictly increasing in U in all cases; largest bracket violation 1.2e-11 in y-units; strict for all pairs with separation >1e-4 |
| `job1_offattractor.py` | held-out spiral vs Hermite section curves; off-attractor orbits (same (U,p), y shifted) | held-out gap ≤ 5e-13 k … 4e-7 k; off-attractor gap/dy matches exp(−∫(1+eta)) to ~1% for dy=±0.01k over four section crossings |
| `job1_checks2.py` | plane-fit gradient error vs stencil radius; crests of successive turns | error decreases with stencil; crests increase monotonically to the cycle crest. (Its polyline-crossing test X1 is INVALID — chord artefacts — and is superseded by `job1_section.py`.) |

Parameter pairs (mu,rho): (1,1.1), (1,2), (0.3,1.5), (3,1.3) [+ (0.1,1.3), (1,3) in the section test].

## 2. Job 2 — uniqueness for the reduced planar system: NOT PROVED

Reduced system on D: u'=v, v'=k−mu v−h(u,v); h convex, h_u=zeta>0, h_v=−eta<−max(0,mu−1); div = eta−mu;
D = {u_min ≤ u ≤ u_max, v_−(u) ≤ v ≤ v_+(u)} (vertically convex). Useful exact form:
v' = −g~(u) − Psi(u,v), g~(u)=h(u,0)−k=k(e^u−1)/(1+eta(u,0)) (increasing, convex), Psi(u,v)=∫_0^v(mu−eta(u,s))ds
(convex in v), energy E=v²/2+G~(u) with E' = −v Psi(u,v).

### 2(a) Zhang-type comparison of divergence integrals — CONDITIONAL PROOF + numerically supported inequality

- Proved (Prop. 8.1 of the tex): on the crest half-axis {v=0,u>0}, R_1'(u) = g~(u) e^{iota(u)}/g~(R_1(u)),
  iota(u) = ∫ over the turn of (eta−mu) dt; in the energy coordinate the return map has derivative e^{iota}.
  If iota is strictly decreasing (energy return map strictly concave), or if R_1(u)/u is strictly decreasing,
  the cycle is unique.
- Numerics `job2_a.py` (14 pairs: mu ∈ {0.03,0.1,0.3,1,3,10}, rho from 1.02 to 3; 12 spirals each, event-refined):
  iota on the CREST section is strictly decreasing from 2 pi (eta_E−mu)/omega to log m_t in all 14 cases.
  **The same statement on the TROUGH section is false** (iota_trough is not monotone; increments up to +7e-3), and the
  half-turn integrals behave oppositely: crest→trough is increasing, trough→crest decreasing.
- `job2_a3.py`: R_1 is concave and sublinear in the coordinate u (14/14) and in the energy coordinate (14/14), in
  the Filippov coordinate w=sqrt(2(e^u−1−u)) in 13/14 (fails at mu=1, rho=3), and is NOT concave/sublinear in x−k or
  y−k (0/14). So the property is coordinate- and section-specific.
- Exact missing inequality (derived; verified against finite differences in `job2_a2.py` to ~1e-3 relative; the
  derivation is a routine Jacobi-field computation and is NOT written at referee standard):
  d iota/du = ∫_0^{tau(u)} g~(u) e^{I(t)} Xi / |F|^4 dt < 0, with I(t)=∫_0^t(eta−mu), F=(v,v'),
  Xi = (h_uv v' − h_vv v)(v²+v'²) − (eta−mu)[(1−zeta)(v²−v'²) + 2(eta−mu) v v'].
  Xi changes sign on every turn (positive on 19–80% of the turn), so this is a genuinely integral inequality.
- Status: supported, unproved. A Zhang-type comparison restricted to closed cycles (outer cycle more stable) would
  not suffice by itself anyway: it leaves a semistable innermost cycle; the statement for non-closed turns is needed.

### 2(b) Liénard plane U'=P(U,s), s'=−k(e^U−1) — pointwise candidates REFUTED numerically

- Proved: P is (the restriction of) a jointly convex function, P_s=1/(1+eta), P_U=(zeta−mu)/(1+eta),
  P_UU=(1,P_U)D²h(1,P_U)^T/(1+eta) ≥ 0; the U-nullcline s = mu U − g~(U) is concave. (Remark 8.2 + elementary.)
- `job2_lib.py` computes D²h on W^u(E) from a transversal Jacobi field and the linearised Riccati equations
  (symmetric to ≤1e-3 except close to the cycle at (1,1.1); positive definite on the sample: h is strictly convex inside).
- `job2_b.py`: the generalised Zhang condition d/dU|_s (P_U/g) ≤ 0, i.e. Z=P_UU g − P_U g' ≤ 0, is VIOLATED on 4–18% of
  the disc for U>0 and 2–18% for U<0 (7 parameter pairs), with full-size violations. Also zeta>mu on {U>0} is false
  (1–5% of points), so P_U changes sign on the right half.
- `job2_b2.py`: the analogous condition in the plane (y−k, m=p+mu U), where the "g" is linear
  (y'=k e^{𝒰(y,m)}−y, m'=−(y−k)), is violated on 34–60% of the disc.
- `job2_massera.py`: Massera-type radial certificates = Dulac weights 1/(X∧xi) in 12 charts
  ((U,p),(U,y),(p,y),(y,m),(U,s),(w,s),(w,p),(x,p),(x,y),(w,y),(1−e^{−U},p),(1−e^{−U},s)): the cycle and all
  spirals are star-shaped in essentially every chart, but beta = div − D log|X∧xi| is positive on 21–53% of the
  disc in every chart. Best are the Filippov charts (w,s),(w,p) (21–42%).
- Exact reason (proved, one line): for radial scaling in (U,p), F(theta q) ∧ F(q) = −v [theta(h(q)−k) + k − h(theta q)],
  and the bracket is ≥0 by convexity. So the flow crosses the scaled orbit arcs outward where p>0 and inward where
  p<0. In 3D: for nu=(U,p,y−k), D omega − C omega = k(Ue^U−e^U+1)·(0,p,p') for omega=conj(f∧nu). Convexity gives a
  comparison that is odd under the half-turn, whereas uniqueness needs a one-signed one.

### 2(c) Curves without contact / Dulac functions — explicit candidates REFUTED; LP leads

- `job2_c.py`: D Lambda + eta − mu is sign-changing on the disc (positive on 20–65% of points) for
  Lambda ∈ {0, −log W+2 log zeta, −log K+2 log zeta, −log Lambda_E+log zeta, −log W+log zeta+log eta,
  −log K+log zeta+log eta, log zeta−log y−log(lambda−p), log zeta, log eta, −log Q+2 log zeta, −2 log Q+3 log zeta}
  (Q = W − 2 sigma K' + 2 zeta_E K, the gauge whose leading part at E is the invariant quadratic form).
  An LP over the span of 14 such functions (logs of W,K,zeta,eta,y,lambda−p,Lambda_E,Q,1+eta and U,p,y,eta,zeta)
  has optimum t* = +0.07…+0.15 > 0 in all 8 runs: no pointwise Dulac function in that span.
- `job2_c2.py`: LP over transversal fields nu with components of degree ≤2 in (u,v,y−k,e^u−1) (28 unknowns), Dulac
  weight 1/(X∧nu): FEASIBLE (tau*<0) in all 8 runs on r>0.05·crest and r>0.3·crest; validation on unsampled
  points: beta ≤ 1e-4. The inequality involved is affine in (eta,zeta) with state-polynomial coefficients.
- `job2_c3.py`, `job2_c4.py`, `job2_c5.py`, `job2_c6.py`: state-explicit version (no eta, zeta): nu a vector field in
  (U,p,y−k) with components of degree ≤2 in (U,p,y−k,e^U−1) (42 unknowns), omega=conj(f∧nu),
  D omega − C omega = conj(f∧(D nu f − J nu)). Results: see §3.5 (mixed; degree 1 infeasible).
- The energy family: E' = v²(etabar−mu), etabar = vertical mean of eta; {eta>mu} is the region BELOW a graph
  v=v_0(u) (because eta is nonincreasing in v). Level curves of E localise cycles but cannot give uniqueness.

### 2(d) Rotated vector fields in k — natural candidate REFUTED; and insufficient even if true

- `job2_d.py` D1: cycles at different k (mu fixed, rho from 1.05 to 3) are nested (0 crossings) in the planes
  (U,p), (U+log k,p), (U,s/k), (U,p/lambda).
- D2: in the plane (U, s/k) (where s' does not depend on k after scaling) the family is NOT rotated:
  (e^U−1) ∂_k P~_k is positive on only 69–74% of the disc, with negative values of comparable magnitude.
- Brief (h) (attractor at k2 lies above the convex envelope at k1 in V=log x coordinates) is one-sided: it gives
  ∂_k h ≥ 0 but no upper bound, so 1 − ∂_k h has no determined sign.
- Even a perfect rotated family only gives: each point of a section is the crest of a cycle for at most one k
  (amplitude curve k = kappa(s) single-valued in s). Uniqueness is the statement that kappa is monotone; an S-shaped
  curve (fold pair) is compatible with rotation. So (d) cannot close the conjecture without an extra monotonicity.

## 3. Leads, obstruction, next steps

### 3.5 LP leads (numerical only — NOT a proof of anything)

Setting: nu a vector field on R^3 with components polynomial of degree ≤2 in (U, p, y−k, e^U−1) (42 unknowns);
omega = conj(f∧nu) (conj(W)=(W12,−W13,−W23)); then D omega − C omega = conj(f∧(Dnu f − J nu)) exactly, C the compound
matrix of v2.0 (4.6). omega>0 and D omega − C omega ≥ 0 (≢0) along a cycle imply lambda_c<0 (brief (e)); nothing unknown
(eta, zeta, h) appears. The LP asks for this on sampled attractor points with planar radius > 0.05·crest.

- Radial field nu=(U,p,y−k): omega>0 on the whole sampled attractor (all three planar projections are star-shaped),
  but D omega − C omega = k(Ue^U−e^U+1)(0,p,p') changes sign. Degree-1 span (+e^U−1): LP infeasible (t*<0) at all 4
  parameter pairs tried (`job2_c3.py`).
- Degree-2 span, `job2_c5.py` (full sample, raw normalisation, coefficient bounds ±2000; validation on the cycle itself
  and on held-out spirals):
  feasible AND validated out of sample at (mu,rho) = (1,1.02), (1,1.1), (1,1.5), (0.3,1.1), (0.3,1.5);
  feasible on the fitted sample but NOT validated (negative values on the cycle / held-out spirals) at (1,2), (1,3),
  (1,5), (3,1.3) [marginal], (3,2.5); LP infeasible within the bounds at (0.3,3), (0.1,1.3), (0.1,2); t*≈0 at (10,1.5).
- `job2_c6.py` (scale-invariant normalisation omega ≥ omega_rad, unbounded coefficients, cycle points included,
  5001-point LP): feasible with positive margin at (1,1.02), (1,1.1), (1,1.5), (1,2); held-out spirals violate the
  second inequality on 0.06–0.07% of points in three of the four cases (0 at (1,1.5)); coefficients are huge
  (1e5–1e8), i.e. the certificate is badly conditioned. The run was stopped at (1,3) after ~7 CPU-minutes without an
  LP solution; the remaining grid was not run. `job2_c6.json` holds the four completed cases.
- Sparse (L1-minimal) solutions (`job2_c4.py`) have the common shape nu ≈ c·[(U,p,y−k) + quadratic correction], with
  the U²-coefficient of the first component ≈ −0.21…−0.25 in the four mu ≥ 1 cases.
- Planar relaxation (`job2_c2.py`, inequality affine in eta, zeta): feasible in all 8 runs.

Reading: near Hopf and for mu of order one there is numerically a low-degree, state-explicit compound supersolution
along the attractor; it is not robust across parameters in this basis and these coordinates, and nothing here is
validated on a neighbourhood of the attractor, let alone proved. It is a lead for the supersolution lane, not a result.

### 3.6 Sharpest obstruction

The reduction keeps exactly one structural fact about h beyond the first-order PDE: convexity. Every comparison
that convexity yields (scaled orbits, supporting planes, convex backward flow) is one-sided in y, and after
projection it reverses with the sign of p: outward in the upper half-plane, inward in the lower (exact identity in
2(b)). Zhang/Massera/Dulac criteria need a one-signed quantity on a whole annulus; all natural pointwise candidates
tested change sign on 20–60% of the disc. What is true numerically is an integrated, section-specific statement
(concavity of the crest return map in U; false on the trough section and in the coordinates x, y), whose
derivative has a sign-changing integrand. So a proof must use something asymmetric between crest and trough —
the exponential nonlinearity itself — not convexity of h alone.

### 3.7 What I would try next

1. State-explicit compound supersolutions (brief (e)) with nu = radial + quadratic correction, in scaled variables
   (U, p/sqrt(mu), (y−k)/k); fit the correction as a function of (mu,rho); then try to verify the three polynomial-
   exponential inequalities on an explicit region containing 𝒜\N_E (cone exclusions, Lambda_E>0, state box).
   This is the only route found here in which nothing unknown (eta, zeta, h) appears in the inequality.
2. Prove sublinearity R_1(theta u) ≥ theta R_1(u) of the crest return map directly in 3D: the scaled orbit theta·z is an
   exact subsolution (theta k(e^U−1) ≥ k(e^{theta U}−1)); what is missing is a comparison principle on the half-turn
   crest→trough→crest that survives the sign change of p.
3. Large-rho asymptotics (heuristically: in the deep trough eta→max(0,mu−1), suggesting lambda_c→−min(mu,1); not
   checked) to close the non-compact end, as Xie's analysis does for Wright's equation; then the problem is a compact parameter set plus the Hopf collar.
