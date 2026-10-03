# Lane `supersolution` — NOTES (2026-10-03)

Scope: brief items 3(e),(f). Explicit functions of the state only (no eta, zeta inside a certificate).
All numerics are ordinary floating point on 55 cycles
(mu in {0.05,0.1,0.2,0.5,1}, rho in {1.001,1.003,1.01,1.03,1.1,1.3,1.6,2,3,5,10}); nothing is an interval certificate.

**Bottom line.** No proof of the conjecture and no counter-evidence (all 55 cycles attracting; three independent
exponents `log m_t/T`, `<eta>-mu`, `<beta>-A` agree). New rigorous intermediate results: Theorem R (every cycle's
(U,p)-projection is strictly convex and strictly star-shaped about E, with exact formulas
`lambda_c = -<p K/S1> = -<p W/S4>`) and Proposition T (for two nested cycles the radial ratio R_2/R_1 is strictly
monotone on each half-plane, maximal on the trough axis and minimal on the crest axis). Plus exact reductions, a
convex-programming engine that decides whether a log-linear certificate family is alive, and sharp obstructions for
Targets A, B, C.

Notation beyond the brief: `H = y-k`, `P = lambda - p`, `p' = k - mu p - y`, `y' = x - y`,
`G(U) = e^U-1-U`, `G1(U) = U e^U - e^U + 1` (both >= 0), `2 alpha = eta_E - mu = gamma_E - A > 0` (alpha = real part of
the unstable pair of E), (r, theta) = polar coordinates in the (U,p)-plane,

    S1 = p^2 - U p'            S2 = U y' - p H            S3 = p y' - p' H        ( (S1,S2,S3) = f ^ r_E, r_E = z - E )
    S4 = p'^2 - p p''          Q_E = p^2 - 2 alpha U p + zeta_E U^2  (positive definite above Hopf)

---------------------------------------------------------------------------------------------------------------------
## 1. Rigorous results (complete proofs; they use only v2.0 results: K>0, W>0, orthogonality, crest identities, Thm 3.1)

### Theorem R (star-shape, strict convexity, radial exponent formulas) — PROVED; verified numerically on 55 cycles

On every nonconstant positive cycle:

1. `S1 = p^2 - U p' > 0` at every phase (the (U,p)-projection is strictly star-shaped about (0,0); the polar angle is
   strictly decreasing), and `S1' = (eta - mu) S1 + p K`.
2. `S4 = p'^2 - p p'' > 0` at every phase (the (U,p)-projection has nowhere-vanishing curvature: it is a strictly convex
   Jordan curve), and `S4' = (eta - mu) S4 + p W`.
3. Consequently
       eta - mu = (log S1)' - p K / S1 = (log S4)' - p W / S4,
       **lambda_c = - < p K / S1 > = - < p W / S4 >**,
       -lambda_c T = \oint (K/S1) dU = \int_0^{2pi} sin(theta) K r^3 / S1^2  d theta.
   Stability <=> `<pK/S1> > 0`, with K > 0 and S1 > 0: only the sign of p is indefinite.
4. Shifted centres: for every `0 <= c < |U_min|`, `s_c = p^2 - p'(U+c) > 0` and `lambda_c = - < p (K + c zeta) / s_c >`.

*Proof.* (1) `S1' = 2pp' - pp' - U p'' = pp' - U p''` and `p'' = -mu p' - y'`, so `S1' = pp' + mu U p' + U y'`.
On the other hand, using `K = zeta U - eta p - H`, the orthogonality `zeta p - eta p' = y'` (paper (2.9)) and
`-H = p' + mu p`:
`(eta-mu)S1 + pK = U(zeta p - eta p') - mu p^2 + mu U p' - pH = U y' + mu U p' + p p'`. Hence `S1' = (eta-mu)S1 + pK`.
At the two zeros of p, `S1 = -U p' = U (y-k)`, and the crest identity `x - k = (1+eta)(y-k)` (paper (4.1)) gives
`sign(y-k) = sign(U) != 0`, so `S1 > 0` at the crest and at the trough. Let `E(t,s) = exp int_s^t (eta - mu)`.
On the rising arc (trough t_tr to crest, p > 0): `S1(t) = E(t,t_tr) S1(t_tr) + int_{t_tr}^t E(t,s) p K ds > 0` because K > 0
(paper Cor. 3.2). On the falling arc (crest to the next trough t_tr', p < 0) solve backwards from the trough:
`S1(t) = E(t_tr',t)^{-1} [ S1(t_tr') - int_t^{t_tr'} E(t_tr',s) p K ds ] > 0` because `-pK > 0`.
(2) `S4' = p'p'' - p p'''`; with the tangent equation `p'' = (eta-mu)p' - zeta p` (paper (2.10)) and its derivative,
`S4' - (eta-mu)S4 = p'[p'' - (eta-mu)p' + zeta p] + p(zeta' p - eta' p') = p W`. At the zeros of p, `S4 = p'^2 = (k-y)^2 > 0`.
W > 0 (paper Lemma 3.3); the same two-arc argument gives S4 > 0. The curvature of t -> (U,p) is `-S4/|(p,p')|^3`.
(3) Divide by S1 (resp. S4) and average over the period; `d theta/dt = -S1/r^2` and `p = r sin(theta)` give the angular form.
(4) `s_c = S1 - c p'`, so `s_c' = (eta-mu)S1 + pK - c[(eta-mu)p' - zeta p] = (eta-mu)s_c + p(K + c zeta)`; at p = 0,
`s_c = (y-k)(U+c) > 0` at the crest and, because `U_min + c < 0`, at the trough; `K + c zeta > 0`. Same argument. ∎

Numerical check (`radial_check.py`, receipt `radial_check_out.txt`): on all 55 cycles `-<pK/S1>`, `-<pW/S4>` and the
shifted-centre values (c = 0.5|U_min|, 0.9|U_min|) agree with `<eta> - mu` to 6 significant digits; ODE residuals <= 3e-9;
min S1/max S1 in [6e-3, 0.95]; min S4 > 0 everywhere.

*Corollaries (immediate, proved).*
- `(p/U)' = -S1/U^2 < 0` and `(p'/p)' = -S4/p^2 < 0`: `|U| = |log(x/k)|` is strictly log-concave on each of the two
  arcs where U != 0, and `|p|` (i.e. `|1-Y|`) is strictly log-concave on each half-cycle.
- General rule: for any vector v along the cycle, `f ^ v = s_v (1,eta,zeta) + N_v (0,p,p')` with `s_v = (f^v)_1`,
  `N_v = n.v`, `n = (zeta,-eta,-1)`; and `s_v' = (eta-mu)s_v + p N_v + (f^[f,v])_1`. In particular `S2 = eta S1 + p K`,
  `S3 = zeta S1 + p' K`, so `S2 > 0` on the rising arc and `S3 > 0` where `p' >= 0`.
- One-sided explicit bounds: on the rising arc `eta - mu < (log S1)'` and `eta - mu < (log S4)'`; on the falling arc both
  inequalities are reversed. So the only thing missing for a proof along this line is an explicit *upper* bound for eta
  on the falling arc; every pointwise constraint known there (K>0, W>0, zeta>0, the constraint line) bounds eta from below.
- `zeta < y` pointwise: `(zeta - y)' + (zeta - y) = -eta zeta < 0`, so the periodic function zeta - y is negative at its maximum.
- Explicit identities: `S1 = Q_E + U Lam_E`; `D Q_E = 2 alpha Q_E - 2 (p - alpha U) Lam_E`, hence
  `<(p - alpha U) Lam_E / Q_E> = alpha` on every cycle (verified to 9 digits, `attractor_check.py`). The instability
  of E and the stability of the cycle thus have parallel exact forms,
  `alpha = <(p - alpha U) D_h(z,E)/Q_E>` and `-lambda_c = <p D_h(E,z)/S1>`, where `Lam_E = D_h(z,E)` and `K = D_h(E,z)`
  are the two Bregman divergences of the convex graph h between z and E. I could not turn the parallel into an inequality.

### Proposition T (two nested cycles: radial-ratio monotonicity) — PROVED

Let Gamma_1, Gamma_2 be distinct cycles at the same parameters, Gamma_1 inside Gamma_2 (paper Thm 4.1), and let R_i(theta) be
their polar radii in the (U,p)-plane (single-valued by Theorem R). Put `m_i(theta) = (y_i - k)/R_i` and
`sigma_i = S1/R_i^2 = sin^2 + mu sin cos + m_i cos > 0`. Then `m_2 > m_1` for every theta and

    d/d theta  log( R_2/R_1 ) = sin(theta) (m_2 - m_1) / (sigma_1 sigma_2).

Hence R_2/R_1 is strictly increasing on (0,pi), strictly decreasing on (pi,2pi):

    U_min,2/U_min,1  >  R_2(theta)/R_1(theta)  >  U_max,2/U_max,1  > 1        (theta not in {0, pi}),

and `\int_0^{2pi} sin(theta) (m_2 - m_1)/(sigma_1 sigma_2) d theta = 0`. The same holds for polar coordinates centred at any
(-c, 0) with 0 <= c < |U_min,1|.

*Proof.* On any orbit with S1 != 0, from `U' = p`, `p' = -mu p - (y-k) = -r(mu sin + m)`:
`r'/r = sin (cos - mu sin - m)`, `theta' = -sigma`, so `d log r/d theta = Xi(m,theta) := sin (m + mu sin - cos)/sigma(m,theta)`
(checked numerically, `polar_check.py`, 1e-14). Xi is linear-fractional in m with `dXi/dm = sin/sigma^2` (the numerator is
`sin^2 + cos^2 = 1`), and since sigma is affine in m, `Xi(m_2) - Xi(m_1) = sin (m_2-m_1)/(sigma_1 sigma_2)`.
`m_2 > m_1`: fix theta, let `d_1 = zeta_1 cos - eta_1 sin`. Theorem 3.1 of the paper with base z_1 and second solution E gives
K_1 > 0, i.e. `d_1 > m_1`; with base z_1 and second solution z_2 (any phases are allowed) it gives
`(R_2 - R_1) d_1 < y_2 - y_1`. Therefore `y_2 - k > (R_2-R_1) d_1 + R_1 m_1 > R_2 m_1`. For a shifted centre use the convex
envelope h of paper Thm 5.1 for {E, Gamma_1, Gamma_2}: `(h - k)/r_c = (h - h(-c,0))/r_c - (k - h(-c,0))/r_c` is increasing
along rays because h is convex and strictly increasing in U; the polar equation has the same form. ∎

Remarks. (i) Scaled copies of a cycle are crossed inwards on the upper half-plane and outwards on the lower one — this is the
precise reason Massera's radial argument does not close. (ii) The infinitesimal version of the integral identity is the
neutrality condition `\int sin(theta) K r^3/S1^2 d theta = 0`. (iii) I did not find a second relation contradicting the
monotonicity; `<x> = k`, `<y> = k` on both cycles read `\int (e^{R_i cos} - 1)/sigma_i = 0`, `\int R_i m_i/sigma_i = 0`.

### Lemma A0 (the certificate is a direction field; concavity) — PROVED

For C^1 positive `phi, eta~, zeta~` along a cycle put `omega = phi (1, eta~, zeta~)`. Then `omega' = C omega + diag(g) omega` with
`g1 = Dlog phi + mu - eta~`, `g2 = Dlog phi + Dlog eta~ + 1 - zeta~/eta~`, `g3 = Dlog phi + Dlog zeta~ + A - x/zeta~`, and
`lambda_c <= -<min_i g_i>`. Moreover `<min_i g_i>` does not depend on phi:

    S[eta~,zeta~] := < min( mu - eta~ ,  Dlog eta~ + 1 - zeta~/eta~ ,  Dlog zeta~ + A - x/zeta~ ) > ,   lambda_c <= -S,

with equality (S = -lambda_c, the maximum over all positive directions) at the true (eta, zeta).
Equivalent: `< eta~ + max(0, -R1/eta~, -R2/zeta~) > < mu` with R1, R2 the two Riccati residuals of the trial.
If `log eta~ = th1.F(z)`, `log zeta~ = th2.F(z)` for any fixed feature functions F_j(z) of the state, each g_i is a
*concave* function of theta = (th1, th2); hence S(theta) is concave, and so are the joint max–min over several cycles and
the pointwise version `min_{t,i} g_i` with `log phi = th0.F`.

*Proof.* With `m = min_i g_i` and `E(t) = exp(-int_0^t m)`, `w = E omega` satisfies `w' = C w + diag(g-m) w >= C w` componentwise.
C is Metzler, so its evolution Phi_C(T,0) is a nonnegative matrix and `Phi_C(T,0) omega(0) <= w(T) = e^{-T<m>} omega(0)`,
omega(0) > 0. Collatz–Wielandt gives `spr Phi_C(T,0) <= e^{-T<m>}`; the spectrum of the compound monodromy is
{m_s, m_t, m_s m_t}, so spr = m_t = e^{lambda_c T}. phi drops because `<Dlog phi> = 0`. Each g_i is "linear minus a positive
multiple of exp(linear)" in theta. ∎

So the search inside a log-linear family is a convex programme: a negative optimum is a genuine refutation of the family
on that cycle (up to solver tolerance), not a failed search.

### Lemma A1 (first-order rigidity at Hopf) — PROVED (under the stated expansion hypothesis)

Let cycles of amplitude eps -> 0 (k = k_H + O(eps^2), x = k(1 + eps u + o(eps)), <u> = 0) carry trial directions with
C^1-uniform expansions `eta~ = e0 + eps a(t) + o(eps)`, `zeta~ = z0 + eps b(t) + o(eps)`. If `S >= 0` for all small eps then
`e0 = z0 = mu` and `(a,b) = (eta_1, zeta_1)`, the first-order part of the true pair. In state variables any C^1 certificate must
therefore satisfy, on the centre plane at Hopf (D = mu^2 + 6 mu + 1),

    eta~ = mu + (mu A/D) U - (2 mu/D) p + ...,      zeta~ = mu + (mu(mu^2+4mu+1)/D) U - (mu A/D) p + ...

i.e. it must encode the 2-jet of W^u(E); only the second-order terms carry slack (of size |lambda_c| ~ eps^2).

*Proof.* Order 0: `min(mu - e0, 1 - z0/e0, A - k_H/z0) >= 0` forces `e0 <= mu`, `z0 <= e0`, `z0 >= k_H/A = mu`.
Order eps: `g = eps (-a, (a'+a-b)/mu, (b' + A b - mu A u)/mu) + o(eps)`. Let m_1 be the minimum of the three brackets;
`<m_1> >= 0` is necessary. Averaging the three inequalities `bracket_i >= m_1`: `<a> <= 0`, `<a> - <b> >= 0`, `A<b> >= 0`,
hence `<a> = <b> = <m_1> = 0` and each bracket equals m_1 identically. The resulting equations `a' = b - A a`,
`b' = mu A u - A b - mu a` are the linearised Riccati system, whose periodic solution is unique (matrix eigenvalues -A ± i sqrt(mu)).
The displayed coefficients solve it (checked by hand against the paper's c20, c11, c02; numerically the optimal exponents
of the monomial family at mu=1, rho=1.001 give U-coefficient 0.2496 and p-coefficient -0.2497 against 1/4, -1/4). ∎

Consequences: `omega^Y`, `omega^c`, every constant direction, the constraint-line bounds and the radial direction
(S2/S1, S3/S1) are first-order inexact and therefore dead near Hopf, each for this reason.

### Proposition B (Target B identity, quantitative; strong-stable ratios) — PROVED; verified numerically

Let (a,b,c) > 0 be the strong-stable Floquet solution (`a' = -b, b' = -mu b - c, c' = -c - x a`, decay e^{-gamma t}), hats =
periodic parts. Then `nu = e^{At}(c,b,a)` solves `nu' = -C^T nu`, so for every C^1 periodic vector omega(t):
`lambda_c <nu^.omega> = -<nu^.(omega' - C omega)>`. With omega^Y = (P, y, x) and omega^c = (1,mu,mu):

    lambda_c = < a^ x (lambda - A - 2p) > / < c^ P + b^ y + a^ x >      (ratio check 1.000000 on 54/55 cycles, 0.9989 at mu=1, rho=10)
    lambda_c = < a^ (x - k_H) > / < c^ + mu b^ + mu a^ >,    chi(gamma) = <x a^>/<a^>,  chi(s) = s(s-1)(s-mu),
    lambda_c = <beta> - A,   beta = b/a,   k = <chi(beta)>.

For mu <= 1: `beta > 1` and `kappa = c/a > (1-mu) beta` pointwise.
*Proof of the last line.* `v = c - (1-mu) b` has `v' = -mu v - x a`, so `e^{mu t} v` is strictly decreasing; `v = e^{-gamma t} v^`
with gamma > 1 >= mu; if `v^(t0) <= 0` then `e^{(mu-gamma)(t0+T)} v^(t0) < e^{(mu-gamma) t0} v^(t0)` forces `v^(t0) > 0`,
contradiction. Then `(b-a)' = -v < 0` and the same periodicity argument gives b > a. ∎
(Analogue of the paper's eta > max(0, mu-1). Jensen on k = <chi(beta)> only gives gamma < gamma_E, i.e. lambda_c < 2 alpha.)

### Reductions (proved identities; no sign information by themselves)

- *Constraint line.* The true (eta,zeta) lies on the explicit line `zeta p - eta p' = y'`. Parametrised by K:
  `eta = (S2 - pK)/S1`, `zeta = (S3 - p'K)/S1`, and K is a positive periodic solution of the **scalar** Riccati equation

      (R_K)   K' = k G1(U) - (1 + S2/S1) K + (p/S1) K^2 ,     S2/S1 = mu + (log S1)',

  with explicit state-dependent coefficients (`riccK.py`: residual 1e-10; `<S2/S1> = mu` to 1e-12). With N = K/S1:
  `N' = kG1/S1 - (1 - mu + 2 S2/S1) N + p N^2`, `lambda_c = -<pN>`.
  Scalar comparison (sub/super-solutions of (R_K) give pointwise bounds on K) is only SKETCHED: the Möbius wrap-around
  through K = infinity has to be excluded.
  Target C in this form: the neutral-mean explicit trial (S2/S1, S3/S1) has `d_R = eta - S2/S1 = -pK/S1`,
  `e_R = zeta - S3/S1 = -p'K/S1`, i.e. a fully explicit sign pattern (-sign p, -sign p').
- *Explicit almost-solutions of the compound system* `w' = Cw` and their defects `w' - Cw`:
  `omega^Y = (P,y,x)`: `(0,0,x(2p+A-lambda))`;  `omega^c = (1,mu,mu)`: `(0,0,k_H - x)`;
  `omega^R = f^r_E = (S1,S2,S3)`: `k G1(U) (0,p,p')`;  `f^e_U = (-p',y',0)`: `x(0,p,p')`;
  chemical radial `f^(1-e^{-U},p,H)`: `(1-e^{-U}) p (p',-y',0)`. In general `(f^v)' = C(f^v) + f^[f,v]`.
- *Markov form.* In ratio coordinates `w = (P r1, y r2, x r3)`: `r1' = (y/P)(r2-r1)`, `r2' = (x/y)(r3-r2)`,
  `r3' = P(r1-r3) + (lambda - A - 2p) r3`: a cyclic three-state chain with rates (y/P, x/y, P) plus the potential
  `lambda - A - 2p` on state 3; lambda_c is its Feynman–Kac exponent.
- *Integral form.* `K_A = (D+A)^{-1}(D+1)^{-1}(D+mu)^{-1} x` (causal filters) has spectral radius < 1 iff the cycle is
  attracting; x a^ is its adjoint eigenfunction; scalar criterion: Phi > 0 with `Lc[Phi] >= 0` suffices (positivity of
  the intermediate (D+mu)Phi, (D+1)(D+mu)Phi is then automatic).

### Numerical-only / conditional
- S1 > 0 and S4 > 0 on all W^u(E) spirals tested (16 parameter pairs x 2 phases, `attractor_check.py`). The proof of
  Theorem R goes through along any bounded entire orbit carrying the positive adjoint direction with K>0, W>0 that crosses
  p = 0 before and after the time considered — i.e. conditional on brief 3(b). Then every ray from E is a global section of
  the planar flow on the attractor and `B = 1/S1` is a smooth positive weight on D\{E} with `div(F/S1) = -pK/S1^2`.
- `<p Lam_E/S1> > <p K/S1> > 0` on all 55 cycles (ratio 1.13 … 870); `S3 >= 0` on all 55 cycles; `N - 2<w> > 0` (Section 3). Unproved.

---------------------------------------------------------------------------------------------------------------------
## 2. Target A — numerics (engine `cvx.py`; drivers `famA2.py`, `famA_joint.py`, `famA_joint2.py`, `famA_pointwise.py`, `explicit_trials.py`; receipt `table_A_out.txt`)

Score = S/(-lambda_c) in (-inf, 1]; a family is alive on a cycle iff its (concave-programme) optimum is > 0.

| family (direction eta~, zeta~) | result |
|---|---|
| `omega^Y` = (y/P, x/P) | S = -<(lambda-A-2p)_+> < 0 on all 55 cycles (-0.03 … -35). Dead (first-order inexact). |
| `omega^c` = (mu, mu) | S = -<(x-k_H)_+>/mu < 0 on all 55. Dead. |
| radial T_R = (S2/S1, S3/S1) | mean of eta~ exactly mu, residuals kG1 (p,p')/S1; S <= 0 by construction; S2 changes sign. Dead. |
| T_K: K <- Lam_E on the constraint line (first-order exact, parameter-free) | S < 0 wherever it is a positive direction: S/(-lambda_c) -> -2.2 (mu=1), -2.8 (0.5), -7.5 (0.2), -24 (0.1), -100 (0.05) as rho -> 1; loses positivity for rho >~ 1.03–1.3. Dead (wrong second-order balance). |
| exp-linear: log eta~, log zeta~ affine in (U, p, (y-k)/k) | dead on all 55 cycles (optimum -0.24 … -1.87). |
| monomials x^a Z^b Y^c (features U, log y, log P), per-cycle optimum | **dead for mu <= 0.2 near Hopf**: optimum -0.41 (mu=.05), -0.40 (.1), -0.10 (.2) at rho=1.001, negative up to rho ~ 1.03; alive per cycle for mu = 0.5 (0.34 … 0.74) and mu = 1 (0.27 … 0.82). Cross-checked with Nelder–Mead/Powell (same optimum to 3 digits). |
| + Lam_E exponent (U, y, P, Lam_E) | negative at mu = 0.05 (rho <= 1.01: -0.12 … -0.08) and mu = 0.1 (rho <= 1.003: -0.04, -0.03); positive elsewhere (0.002 … 0.88). Dead near Hopf at small mu. |
| monomials in X, Z, Y, Z* (features U, log y, log P, log P*), `P* = P + P' = (1-mu)P + y` (the other filter order) | **alive per cycle on all 55** (0.585 … 0.985 for mu < 1, -> 1 as rho -> 1 because 5 parameters per component match the whole second-order jet; degenerate at mu = 1 where P* = y: 0.27 … 0.82). |
| same, exponents and prefactors common to all rho at fixed mu | none found (best min-ratio -0.96 … -20; joint soft-min solver, values are lower bounds of a concave optimum). |
| same, common exponents, free k-dependent prefactors, rho in [1.001,10] | none found (best min-ratio -0.46 (mu=.5) … -15 (mu=.05)). |
| same, rho in [1.001, 1.1] only (common exponents, free prefactors) | alive for mu = 0.2, 0.5, 1 (+0.34, +0.62, +0.27); none found for mu = 0.05, 0.1 (-0.55, -0.08). |
| same, rho in [2, 10] only (everything common) | alive for mu <= 0.5 (+0.27 … +0.37); mu = 1 marginal (-0.005). |
| pointwise on-cycle version (min over t and i, scalar factor phi in the same family), X,Z,Y,Z* | best found (soft-min solver; a positive value is a genuine on-cycle certificate, a negative one is a lower bound, not a refutation): positive only near Hopf (rho <= 1.03 for mu <= 0.1, rho <= 1.1 for mu = 0.2, 0.5; never for mu = 1); negative for rho >= 1.3 (-0.2 … -3.6). |
| pointwise, features + Lam_E, S1 | solver lower bounds only: positive up to rho = 1.1 (mu <= 0.1), up to rho = 2 (mu = 0.2, 0.5), for 1.1 <= rho <= 2 (mu = 1); negative for rho >= 3. The low-temperature soft-min is ill-conditioned (one value is below that of a sub-family), so these are not refutations. Runs with more singular features (S4, Q_E) did not finish in the time box. |

Reading: the averaged criterion can be met cycle by cycle once the family has >= 4 exponents per component, but
(i) the exponents must vary with k (they must track the 3-jet of W^u(E) near Hopf); (ii) an averaged inequality along an
unknown cycle is not provable without the cycle — the provable form is the pointwise one on a region known to contain all
cycles, and already on the cycle itself the pointwise monomial certificate disappears for rho >~ 1.3.

---------------------------------------------------------------------------------------------------------------------
## 3. Target B — numerics (`targetB_identity.py`, `shapes.py`; receipts `targetB_out.txt`, `shapes_out.txt`)

- Identity verified (Proposition B). `chi(gamma) = <x a^>/<a^>` to 6 digits, and < k_H on all 55 cycles.
- Shapes: a^, beta, kappa are unimodal in every resolved case (2 turning points; the extra turning points reported for
  small mu and large rho are 1e-12-level noise where beta -> 1+ in the trough). a^ peaks on the rising arc
  (phase 0.63–0.98 after the crest) and is minimal just after the crest (phase 0.001–0.18); max a^/min a^ = 1.3 … 1e3.
- Equivalent forms: `(gamma - A)(N - 2<w>) = <w (lambda + A - 2 beta)>`, w = x a^, N = <c^P + b^y + a^x>; `N - 2<w> > 0`
  on all 55 cycles, so stability <=> `<w beta>/<w> > (lambda+A)/2` <=> `<w p>/<w> > (lambda + A)/2 - gamma`
  (the left side exceeds (lambda+A)/2 by only 0.0004 … 0.27).
- No rearrangement/Chebyshev argument found: the sign is decided by the phase lag between a^ and x (a^ is in the middle of
  its range at the crest of x), a second-order effect near Hopf.

## 4. Target C — numerics (`shapes.py`, `prof_print.py`)

- `<d> = lambda_c` (d = eta - y/P) verified; d and e each change sign exactly twice per period (resolved cases); d > 0 on
  50–95 % of the period with a short deep negative excursion before the crest; `<e> < 0` on all 55 cycles.
- The (d,e) system is a forced focus (`d' = e - a d`, `e' = -zeta~ d - b e - zeta~ V~`), not cooperative; no
  sign-definite quadratic balance: `(zeta~ d^2 + e^2)/2` has derivative `(zeta~/2)(p + mu - eta~ - 2a) d^2 - b e^2 - zeta~ V~ e`, indefinite.
- Better form of the same question (Theorem R): with the radial trial instead of the Y-chain,
  `(d_R, e_R) = -(K/S1)(p, p')` exactly, so Target C becomes `<pK/S1> > 0` with K > 0 solving (R_K).

---------------------------------------------------------------------------------------------------------------------
## 5. Obstructions (sharpest form)

1. **First-order rigidity** (Lemma A1): any state-function certificate must reproduce the 2-jet of W^u(E) at Hopf
   exactly and has only O(eps^2) slack; every natural explicit direction (chains, radial, constants) fails at order eps,
   and the parameter-free first-order-exact one (K <- Lam_E) has the wrong second-order balance.
2. **Averaging needs the cycle.** `lambda_c <= -<min g>` is an average along an unknown orbit; the only exact averages
   available are of total derivatives and `min` is not linear, so an averaged certificate is not a proof. A proof needs
   pointwise inequalities `g_i(z) >= 0` on an explicit region containing all cycles, with omega singular at E.
   Near E this is easy at fixed k (omega ~ Q_E^{1+delta} w_E gives g ~ delta*alpha); the difficulty is global, and on the
   cycles themselves pointwise monomial certificates stop existing at rho ~ 1.3.
3. **Sign of p.** Every radial/scaling identity has the form `(log s)' = eta - mu + p N/s` with N > 0 (N = K, W, K + c zeta,
   or -Delta for a second solution). On the rising arc this is an explicit upper bound for eta; on the falling arc it is a
   lower bound. A proof along these lines needs one explicit upper bound for eta (equivalently for K, or K/S1) on {p < 0}
   that is first-order exact at E. Proposition T is the two-cycle face of the same obstruction: radial ratios grow on one
   half-plane and shrink on the other.
4. Targets B and C produce only identities; the sign is carried by a phase lag, not by a monotone or rearrangement structure.

## 6. What I would try next

- (R_K) on the falling arc only: explicit sub/super-solutions `K~ = Lam_E (1 + a U + b p)`-type with a, b from the 3-jet of
  W^u(E); make the scalar comparison rigorous.
- Hand Theorem R to the Dulac lane: with `B = S1^{-a} S4^{-(1-a)} K^{c1} W^{c2} zeta^{c3} eta^{c4} e^{Psi(z)}`, c1+…+c4 = 0,
  `Dlog B + eta - mu = -p[a K/S1 + (1-a) W/S4] + c1 kG1/K + c2 x p^2/W + c3 x/zeta + c4 (zeta/eta + A) + D Psi`; whether a one-signed
  combination exists along W^u(E) spirals + the cycle is a linear programme (not run here: it uses eta, zeta).
- Exploit Proposition T: look for a second monotonicity (in the hodograph plane (p,p'), where S4 > 0 gives star-shape, or
  for the Liénard variable s) with the opposite sense; a Zhang-Zhifen-type comparison of `\oint (K/S1) dU` on nested cycles.
- Finish the pointwise programme with an exponential-cone LP solver instead of soft-min, on a 3D neighbourhood of the cycle.

---------------------------------------------------------------------------------------------------------------------
## 7. Files and commands (python = /Users/admin/.venvs/analysis-py313/bin/python -W ignore, run inside this directory)

    base.py                 build cycles.pkl (55 cycles; log-variable integration; eta,zeta forward, beta,kappa backward)   ~25 s
    targetB_identity.py     Proposition B identities                       -> targetB_out.txt
    shapes.py               Target B / C shape tables                      -> shapes_out.txt
    prof_print.py MU RHO N  sampled profile of one cycle;  prof_K.py MU RHO N  the same with K, S1, pK/S1
    radial_check.py         Theorem R: S1, S4 > 0, exponent formulas, ODE residuals, beta > 1   -> radial_check_out.txt
    attractor_check.py      S1 = Q_E + U Lam_E, <(p-alpha U)Lam_E/Q_E> = alpha, S1,S4 > 0 on W^u(E)  -> attractor_check_out.txt
    polar_check.py          polar orbit equation of Proposition T
    riccK.py                scalar reduction (R_K), K/Lam_E range, <S2/S1> = mu     -> riccK_out.txt
    bounds_line.py          constraint-line bounds on eta (0 violations) and their first-order looseness
    posfuncs.py             positivity scan of candidate explicit functions (S1 > 0, S3 >= 0, S2 indefinite)
    cvx.py                  concave-programming engine
    famA2.py FEATS [q]      per-cycle optimum (FEATS e.g. U,y,P,Pstar)  -> famA2_*.pkl ;  table_A.py -> table_A_out.txt
    famA_joint.py / famA_joint2.py FEATS [rmin rmax]   common exponents (and prefactors / free prefactors)
    famA_pointwise.py FEATS pointwise on-cycle optimum (pw_a.txt = run with U,y,P,Pstar,LamE,star)
    explicit_trials.py      parameter-free trials                          -> explicit_trials_out.txt
    famA.py                 first (Nelder–Mead/Powell) version, kept as the independent cross-check of cvx.py
