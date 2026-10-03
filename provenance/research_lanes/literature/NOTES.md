# Lane `literature` — NOTES (2026-10-03)

Scope: prior-art and applicability search for the Wilhelm–Heinrich (WH) uniqueness conjecture
(`U''' + A U'' + mu U' + k(e^U - 1) = 0`, `k > mu(1+mu)`), questions 1–6 of the lane prompt, plus the
numerical tests required by the honesty rule for every lemma this lane proposed.

Everything is in this directory. `raw/` holds fetched open texts and API dumps (arXiv, OpenAlex, Semantic Scholar,
one author-hosted survey). No file outside the lane directory was written.

---------------------------------------------------------------------------------------------------------------

## 0. Bottom line

1. **No resolution found.** Nothing located up to 2026-10-03 resolves or partially resolves the WH cycle-count
   question beyond Smith (2012) and Boros–Hofbauer (2023, Thm 13 and the paragraph after it). Status: open.
2. **No published theorem applies off the shelf.** Every uniqueness theorem inspected either (a) needs a separable
   planar Liénard structure (`f/g` monotone, Zhang type) that the reduced WH flow does not have, (b) needs a
   Massera/star-shaped hypothesis that was tested here and fails on the cycles themselves, (c) needs a concave
   monotone return map (H. L. Smith 1986 / Krasnosel'skii), tested here and refuted in all power charts, or (d) is
   proved only in a singular regime (steep/step nonlinearity, relaxation limit, near Hopf) or by computer.
3. **The closest solved problem is Jones' conjecture** for Wright's equation (same nonlinearity `e^U - 1`, delay
   operator instead of the two-stage filter). It was not solved by a structural theorem. It was solved by
   "every periodic solution is stable => unique" (Xie) + large-parameter asymptotics (Xie, `alpha >= 5.67`) +
   computer-assisted a-priori bounds and Floquet estimates (`[1.9, 6.0]`) + a separate computer-assisted
   Fourier-space proof near Hopf (`(pi/2, 1.9]`). This is the only route with a track record.
4. **By-product proved in this lane (new relative to the v2.0 manuscript and the brief):** every projected cycle in
   the `(U,p)` plane is strictly star-shaped about `E` and strictly convex, and
   `lambda_c = -<p K / S*> = -<p W / V>` with `S* = p^2 - U p'`, `V = p'^2 - p p''` (Section 3.1). These are exact;
   they do not decide the sign of `lambda_c`.
5. Lane status: **obstruction** (nothing in the literature resolves the conjecture or yields a new region of
   uniqueness without new work).

---------------------------------------------------------------------------------------------------------------

## 1. Q1 — has the WH uniqueness question been resolved?

Evidence gathered (all negative):

| check | tool | result |
|---|---|---|
| works citing Smith 2012 (DOI 10.1007/s10910-011-9946-9) | OpenAlex `cites:W1980509781`; Semantic Scholar graph API | 8 (+1) citing works: multistationarity, bursting, integrability, network-inheritance papers; none on cycle count |
| works citing Boros–Hofbauer 2023 (DOI 10.1016/j.nonrwa.2023.103839) | OpenAlex `cites:W4317940786`; Semantic Scholar | 6–8 citing works (Banaji–Boros–Hofbauer network papers, Craciun–Erban planar cycles, one CRN-computation preprint); none on WH uniqueness |
| works citing WH 1995 / WH 1996 since 2012 | OpenAlex | 45 + 12; none on uniqueness |
| arXiv API | `all:"Wilhelm-Heinrich"`, `abs:"smallest chemical reaction system"`, author listings for Boros / Hofbauer | only arXiv:2202.11034 and arXiv:1109.5204 |
| OpenAlex keyword search | several phrasings, 2012– | nothing relevant |
| WebSearch | 4 phrasings (model names + uniqueness / proof / 2024–2026) | nothing relevant |
| MATCH 94 (2025) 863 (Gao–Jin–Marwan et al., cites Smith 2012; "Hopf and Bautin bifurcations") | full text read | studies the **Willamowski–Rössler** system, not WH. Irrelevant; no conflict with the always-negative WH Lyapunov coefficient |

Exact published status (full text, arXiv:2202.11034v2, after Theorem 13): "As remarked by Hal Smith, the main open
problem about the dynamics of the mass-action system (6) is the number of limit cycles in case
kappa1 - kappa4 > kappa3 + kappa5. He says, numerical simulations suggest there is only one."
Smith (arXiv:1109.5204, after Theorem 1.1): "The main open problem concerning the dynamics of (1.2) is the number of
periodic orbits in case (c). Numerical simulations suggest there is only one."

Limits of this check: Google Scholar itself was not queried (no API); OpenAlex/Semantic Scholar citation graphs are
incomplete for very recent preprints; non-English and unindexed venues are not covered. A negative bounded search is
not a proof of absence.

---------------------------------------------------------------------------------------------------------------

## 2. Source table

Access levels: **FT** = full text read in this lane; **1P** = first page / abstract text only; **AB** = abstract or
publisher record only; **SN** = search-engine snippet only; **CIT** = citation located, content not accessed
(Elsevier open-archive PDFs returned HTTP 403 to automated access and were not circumvented).
Statements marked (mem) are standard statements recalled from the literature and **not** re-verified against the
full text in this lane.

### 2.1 Q1/Q2 — the model itself and the distributed-delay logistic equation

| source | exact statement | access | verdict for WH |
|---|---|---|---|
| H. L. Smith, J. Math. Chem. 50 (2012) 989; arXiv:1109.5204 | Thm 1.1(c): for `k > k3 + k5` every solution off `S` and off the 1-D stable manifold of `E` converges to a nontrivial periodic orbit; at least one is orbitally asymptotically stable; finitely many. Cycle number stated open. Part (b) is proved with the Li–Muldowney second-compound criterion, weights `diag((1-2e)/x, (1-e)/x, -1/k5)`, max-norm: `mu(B) = max{-e k3/(1-e), -e k5, x/(1-2e) - (k3+k5)}`, then `<x> = k`. | FT | Baseline. Smith's own below-Hopf proof is exactly brief item (e) with `phi = 1` (`Lc[1] = k_H - x`); above Hopf the same weights give a positive average, so a state-dependent weight is the literature-consistent open step. |
| Boros & Hofbauer, NARWA 72 (2023) 103839; arXiv:2202.11034 | Thm 13 (i)–(vi): unique positive equilibrium, competitive, permanent, supercritical Hopf at `kappa1-kappa4 = kappa3+kappa5`, GAS for `0 < kappa1-kappa4 <= kappa3+kappa5`, stable limit cycle and finitely many periodic orbits above. Number of cycles: open. | FT | Baseline; open. |
| S. Ruan, "Delay differential equations in single species dynamics" (NATO ASI survey, 2006), §9 | Weak kernel `G = a e^{-as}`: equilibrium globally stable (Thm 7, Lyapunov function). Strong kernel `G = a^2 s e^{-as}`: linear-chain system (9.11) `x' = r x(1 - z/K), y' = a(x - y), z' = a(y - z)`; equilibrium stable iff mean delay `T = 2/a < 4/r`; Hopf at `T = 4/r`, bifurcating solutions orbitally stable (Thm 8). General kernel: local stability if `int s G < 1/r` (Thm 9). No global or uniqueness statement for periodic solutions. | FT | (9.11) **is** WH with `mu = 1` (`lambda = r/a`, Hopf at `lambda = 2 = A`). Confirms brief (i). Only local Hopf results; survey cites Cushing 1977, MacDonald 1978, Fargue 1973, Landman 1980, Simpson 1980 for the same local theory. |
| Cushing (1977, LN Biomath. 20); MacDonald (1978, LN Biomath. 27); Wörz-Busekros, SIAM J. Appl. Math. 35 (1978); Dai, Math. Biosci. 53 (1981) | Linear chain trick; local Hopf bifurcation and stability for gamma kernels; global stability of equilibria for Lotka–Volterra systems with continuous delay under matrix conditions (Wörz-Busekros); existence of nonconstant periodic solutions in predator–prey systems with continuous delay (Dai). (mem) | CIT | No uniqueness or global orbital-stability theorem for the periodic solutions of the gamma/Erlang-kernel logistic equation was found in any of these lines. |
| A. López-Nieto, arXiv:2506.07660 (2025) | For scalar DDEs `x' = r f(x(t), x(t-1))` monotone in the delayed argument, each connected component of periodic orbits (over `r`) is an annulus; `x(t), x(t-1)` solve an integrable planar "predator–prey" ODE; periodic solutions characterised by two time maps. | AB | Discrete delay only. Structural analogue (planar integrable description of all periodic orbits); no distributed-delay version found. |

### 2.2 Q6 — Jones' conjecture (Wright's equation `x'(t) = -alpha (e^{x(t-1)} - 1)`)

| source | exact statement | access | what transfers |
|---|---|---|---|
| X. Xie, thesis (Rutgers 1991); J. Differential Equations 103 (1993) 350 | (i) If `alpha > pi/2` and every slowly oscillating periodic solution (SOPS) is asymptotically stable then the SOPS is unique up to translation. (ii) For `alpha >= 5.67` every SOPS is asymptotically stable (asymptotic description of SOPS for large `alpha`, after Nussbaum, then a Floquet-multiplier bound via an eigenvalue problem with a phase-condition boundary term). Xie (thesis p. 97): the argument "can not be used ... when alpha is close to pi/2". | statements via FT of arXiv:1705.02432 and arXiv:0909.4107; originals CIT | (i) is brief item (a). (ii) is the template for a large-`rho` analytic regime. |
| Jaquette, Lessard, Mischaikow, J. Differential Equations 263 (2017); arXiv:1705.02432 | Unique SOPS for `alpha in [1.9, 6.0]`. Method: reduce every SOPS to three numbers `kappa(x) = (q, qbar, x(1))`; branch-and-prune over a compact region of `R^3` with interval "bounding functions" `l_K <= x <= u_K` valid for **all** SOPS with `kappa(x) in K`; then bound the Floquet multipliers of all such SOPS (Xie's formulation); rigorous DDE integrator. | FT | For the WH ODE the analogue is far cheaper: v2.0 already gives one crest, an explicit box, nestedness; "all cycles" are fixed points of a planar crest map, so bounding functions become interval enclosures of a 2-D return map and its derivative (CAPD-type). |
| van den Berg & Jaquette, J. Differential Equations 264 (2018); arXiv:1704.00029 | Wright's conjecture; and: the Hopf branch has no folds for `alpha in (pi/2, pi/2 + 6.830e-3]`; for `alpha in (pi/2, pi/2 + 5.53e-3]` at most one SOPS (Thm 4.10), via an explicit, scaling-independent neighbourhood of the bifurcation point in Fourier space and a Newton–Kantorovich argument. | FT (abstract, intro, Cor. 4.8, Thm 4.10 located) | Template for an **explicit** Hopf collar. The v2.0 collar (Lemma 7.1) is existential. |
| J. Jaquette, J. Differential Equations 266 (2019); arXiv:1801.09806 | Unique SOPS for `alpha in (pi/2, 1.9]`, completing Jones' conjecture. Stability-based arguments are "infeasible" near Hopf ("dynamics becomes center-like"); instead periodic orbits = zeros of a functional on Fourier coefficients; a-priori bounds confine all zeros to a compact set; branch-and-bound with a Krawczyk operator proves exactly one zero. Also: no isolas. Open: monotonicity of the period in `alpha`. | FT | Near Hopf, switch from "all cycles are stable" to "zero-counting in a compact set of Fourier coefficients". |
| J.-P. Lessard, J. Differential Equations 248 (2010); arXiv:0909.4107 | Rigorous continuation: no folds on the Hopf branch for `alpha in [pi/2 + eps, 2.3]`. | FT (intro) | Branch regularity is weaker than uniqueness (isolas remain); same distinction as "neutral cycle" vs "count" in v2.0. |
| van den Berg, Lessard, Queirolo, SIAM J. Appl. Dyn. Syst. 20 (2021) 573; arXiv:2006.13373 | General method to validate Hopf bifurcations and folds of cycles in ODE systems: blow-up (desingularisation) near the Hopf point turns the branch into a regular continuation problem; validated continuation of cycles from onset. | AB | Direct ODE tool for an explicit collar at fixed `mu` (and, with `mu` as second continuation parameter, on compact `mu`-intervals). |
| Church, Dai, Hénot, Lappicy, SIAM J. Appl. Dyn. Syst. (2026); arXiv:2504.03058 | Computer-assisted **global families** of stable periodic orbits in a 3-D predator ecosystem model, between two transcritical bifurcations (Butler–Waltman problem). | AB | Shows current validated-continuation technology handles a full parameter interval for a 3-D ODE population model; existence/stability of a family, not exclusion of other cycles. |
| Kolesov & Rozov, Sb. Math. 202 (2011) 829; Glyzin, Kolesov, Rozov, Math. Notes 101 (2017) 71 | Hutchinson's equation, large parameter: complete asymptotic expansion of **a** stable relaxation cycle. 3-D repressilator ODE model: existence, asymptotics and stability of the relaxation cycle. | AB | Asymptotic machinery for the large-`rho` regime (relay-type limit in `U = log x`), but these papers construct one stable cycle in a prescribed class; exclusion of other cycles needs an a-priori confinement step as in Xie. |

What a WH analogue would need (2-parameter, non-compact parameter set): (R1) explicit collar `1 < rho <= 1 + eps(mu)`;
(R2) analytic large-`rho` stability of every cycle, uniform in `mu`; (R3) compact middle region by validated
numerics on the crest map; (R4) the `mu -> 0` end (by reciprocity also `mu -> infinity`), which has no counterpart
in the one-parameter Wright problem and must be handled analytically.

### 2.3 Q3 — convex / planar structure of attractors of 3-D competitive and cyclic feedback systems

| source | exact statement | access | verdict |
|---|---|---|---|
| M. W. Hirsch, SIAM J. Math. Anal. 13 (1982) 167 | Limit sets of competitive/cooperative systems: orthogonal projection along any positive direction maps a **limit set** homeomorphically and equivariantly onto an invariant set of a Lipschitz vector field in a hyperplane; in dimension 3 a compact limit set without equilibria is a closed orbit or a cylinder of closed orbits. | AB | Covers limit sets only, not the whole attractor; no convexity. |
| Mallet-Paret & Smith, J. Dynam. Differential Equations 2 (1990) 367 | Poincaré–Bendixson for monotone cyclic feedback systems; planar projection of omega-limit sets is injective (mem). | AB | Same: limit sets; no statement on injectivity for the full attractor, none on convexity, none on uniqueness. |
| H. L. Smith, J. Differential Equations 65 (1986) 361 | Existence, location and stability type of periodic orbits of competitive/cooperative systems; invariant manifolds of periodic orbits; PB for 3-D irreducible systems. | AB | Prior art for the real positive transverse multipliers in v2.0 Prop. 2.1 (not checked line by line). No uniqueness. |
| Zhu & Smith, J. Differential Equations 110 (1994) 143 | Dissipative, 3-D, competitive, irreducible, unique equilibrium with negative Jacobian determinant: equilibrium stable or there is an orbitally stable periodic orbit; if analytic and the equilibrium is unstable there is an orbitally asymptotically stable one. | SN | Existence of a stable cycle only (this is what Smith 2012 used). |
| Gedeon & Mischaikow, J. Dynam. Differential Equations 7 (1995) 141 | Global attractor of cyclic feedback systems admits a Morse decomposition (by the integer Lyapunov function); gradient-like structure between Morse sets; example with shift dynamics inside the lowest Morse set for non-monotone cyclic feedback. | SN | Says nothing about the number of cycles inside one Morse set, nor about convexity. |
| H.-O. Walther, Mem. AMS 113 no. 544 (1995) | For `x'(t) = -mu x(t) + f(x(t-1))`, `f` monotone bounded smooth negative feedback: the attractor in the lowest oscillation class is a **two-dimensional graph, with possibly several nested periodic solutions**. | SN (publisher description) | Closest prior art for brief (b) (attractor = 2-D graph bounded by a periodic orbit). Explicitly compatible with several nested cycles: the graph property alone does not give uniqueness. |
| Ivanov & Lani-Wayda, Commun. Pure Appl. Anal. (2025); arXiv:2408.06548 | (N+1)-dimensional delay systems with cyclic monotone negative feedback: a two-dimensional invariant manifold (strong unstable manifold of 0) on which orbits spiral out to a bounding periodic orbit, assuming essentially only instability of 0. | AB | Same picture as `closure(W^u(E))` in brief (b); no uniqueness. |
| S. Baigent, Proc. Edinb. Math. Soc. 55 (2012) 53; Nonlinearity 26 (2013) 1001; DCDS-B 24 (2019) 1697; Tineo, Appl. Math. Comput. 123 (2001); Zeeman & Zeeman, Trans. AMS 355 (2003) 713 | Carrying simplices of totally competitive Lotka–Volterra (and Kolmogorov map) systems are convex/concave when the flow (map) renders convex all planes whose normals lie in a fixed invariant cone; method: follow the curvature of an evolving surface, Gauss image confined to an invariant cone. Zeeman–Zeeman: if the carrying simplex lies to one side of its tangent plane at the interior fixed point, a split Lyapunov function excludes nontrivial recurrence. | AB | Mechanism identical to brief (b)/(j) (normals `(zeta,-eta,-1)` in an invariant cone + convex nonlinearity). Not a prior statement for WH (not Kolmogorov-totally-competitive; attractor is a disc, not a carrying simplex). The Zeeman–Zeeman conclusion is false here (cycles exist): the analogue `Lambda_E > 0` is true but `Lambda_E'` is not signed. |
| H. L. Smith, Nonlinear Anal. 10 (1986) 1037, as quoted in Farcot & Gouzé, INRIA RR-6875 (arXiv:0903.1969), Thm 1 | `T: [0,p] -> [0,p]` continuous, `C^1` in `(0,p)`, `DT(0)` exists. (M) `DT(x) > 0` for `0 < x < p`. (C) `DT(y) < DT(x)` if `0 < x < y < p`. `Tp < p`. If `T0 = 0` and `spr DT(0) > 1` there is a unique nonzero fixed point `q`, and `T^n x -> q` for every `x != 0`. | FT (of the quoting paper) | **Tested here, refuted** (Section 3.3): the crest map satisfies (M) on the tested rectangles but (C) fails in every power chart tried. |
| Glass & Pasternack, J. Math. Biol. 6 (1978) 207; Farcot & Gouzé, Acta Biotheor. 57 (2009) 429 | Piecewise-linear (step) negative feedback loops: existence, **uniqueness and stability** of the limit cycle (dimension >= 3), via a monotone concave first-return map and the theorem above. | AB / FT of RR-6875 | Works because the return map of the switching system is explicitly fractional-linear. No smooth analogue. |

### 2.4 Q4 — planar uniqueness theorems (reduced system `u'' + mu u' + h(u,u') = k`; Liénard form `U' = P(U,s), s' = -g(U)`)

| source | exact statement | access | verdict |
|---|---|---|---|
| Zhang Zhifen (1958; Appl. Anal. 23 (1986) 63); Zhang et al., "Qualitative Theory of Differential Equations", AMS 1992 | `x' = -phi(y) - F(x), y' = g(x)` (equivalently `x' = h(y) - F(x), y' = -g(x)`): `x g(x) > 0`, `phi` increasing with `y phi(y) > 0`, `f = F'`; if `f(x)/g(x)` is non-decreasing on `(a,0)` and on `(0,b)` (not constant near 0) there is at most one limit cycle in the strip, stable if it exists. (mem + SN) | SN | Needs a **separable** right-hand side. The exact reduced field `P(U,s)` is not separable. Not applicable as stated. |
| Zeng, Zhang, Gao, Bull. London Math. Soc. 26 (1994) 213; Xiao & Zhang, Nonlinearity 16 (2003) 1185; Kuang & Freedman, Math. Biosci. 88 (1988) 67; Hwang, J. Math. Anal. Appl. 238 (1999) 179 | Survey of Zhang-type theorems for `x' = phi(y) - F(x), y' = -g(x)`; a new such theorem and its use for general predator–prey systems after a change of variables; Gause-type uniqueness criterion `d/dx[(x g'(x) + g(x) - x g(x) p'(x)/p(x))/(-gamma + q(x))] <= 0` (mem). | AB | All for separable generalized Liénard systems. No theorem for jointly convex non-separable `P(U,s)` was found. |
| Levinson–Smith 1942; Sansone 1949; Massera 1954; Carletti & Villari, J. Math. Anal. Appl. 307 (2005), arXiv:math/0307372; Sabatini & Villari, Le Matematiche 65 (2010), arXiv:1101.2761 | Thm 1.1 (Carletti–Villari): `F` with exactly three transversal zeros `x2 < 0 < x1`, monotone outside `[x2,x1]`, `G(x1) = G(x2)`, growth condition (E): unique stable cycle. Thm 6 (Sabatini–Villari), system `x' = P, y' = Q` with unique equilibrium at 0: if `alpha := P(xQ_x + yQ_y) - Q(xP_x + yP_y) >= 0` (radial angular monotonicity) and a non-degeneracy condition, at most one limit cycle. | FT | Separable ones: not applicable. RAM: for the reduced field `alpha = -v K` exactly (Section 3.2), which changes sign with `v` on every cycle. **Fails.** |
| M. Sabatini, arXiv:1006.2020 (Chaos Solitons Fractals 2010) | `x'' + x' Phi(x,x') + g(x) = 0`, `g' > 0`, `x g > 0`: if `S := (sgn(x) sqrt(2G)/g)[2G(Phi_x g - Phi g')/g^2 + Phi] + y Phi_y != 0` (strict star-shapedness after the Conti–Filippov transformation) there is at most one limit cycle. | FT | Covers inseparable mixed terms, hence the exact form of the reduced WH equation. **Tested, fails:** `S` takes both signs on the cycles (Section 3.2). Structural reason: `grad Phi(E) != 0`. |
| Gasull & Giacomini, J. Differential Equations 185 (2002) 54; arXiv:1305.3402; arXiv:2101.03874 (J. Differential Equations 2021) | Extended Bendixson–Dulac: if `M_s := <grad V, X> - s V div X` does not change sign on a region, limit cycles do not cross `{V = 0}`, at most `l(U_i)` cycles in each component, all hyperbolic, stability = sign of `-V M_s`. Thm 2.3 (2021): with `V = Q Pdot - P Qdot` (orbit-curvature numerator) and `s = 1`, `M_1 = P^2 Q(P_xx - 2Q_xy) + P Q^2(2P_xy - Q_yy) + Q^3 P_yy - P^3 Q_xx`. | FT | This is the published form of brief item (d). The curvature choice was evaluated: `M_1 = p W` (Section 3.1), sign-changing, so Thm 2.3 does **not** apply; but it yields the new exact formula and the convexity lemma. |
| Li & Zhang, J. Differential Equations 124 (1996) 407 | Criterion for monotonicity of the ratio of two Abelian integrals over ovals of `y^2/2 + Psi(x) = h`. | AB | Tool for any future near-Hamiltonian wedge (`rho - 1 = O(mu)`); see Section 3.4 for why the naive fixed-`rho` reduction is void. |
| Duff 1953; Perko, "Differential Equations and Dynamical Systems" §4.6 | Rotated vector fields: in a one-parameter family with `P Q_lambda - Q P_lambda` of one sign, limit cycles expand/contract monotonically, semistable cycles split or vanish. (mem) | CIT | Not applicable as stated: for the reduced field in the parameter `k` the determinant is `v(1 - d_k h)`, which changes sign with `v`; and a rotated family still allows semistable cycles to appear. |

### 2.5 Q5 — 3-D cyclic negative feedback, third-order equations, compound matrices

| source | exact statement | access | verdict |
|---|---|---|---|
| S. P. Hastings, Rocky Mountain J. Math. 7 (1977) 513 | Goodwin-type 3-D system with parameters `alpha, beta, gamma, rho`: if `rho` (steepness) is sufficiently large and a parameter combination is sufficiently close to a critical value, there is a unique non-constant periodic solution attracting all trajectories except those near the equilibrium or on a one-dimensional curve; structurally stable. | 1P | Singular (near-step) regime only. Method (explicit return map near the switching limit) does not use any structure WH lacks, but WH's nonlinearity is not near a step for moderate `rho`. |
| Bingxi Li, Nonlinear Anal. 5 (1981) 13–19, "Uniqueness and stability of a limit cycle for a third order dynamical system arising in neuron modelling" (DOI 10.1016/0362-546X(81)90065-1) | This is the "1981 neuron-model paper" of the bundle's SOURCES.md. Model: Stein et al. (Kybernetik 1974) neuron model, for which an der Heiden (Biol. Cybern. 1976) proved existence. Reference list (Crossref): Stein; an der Heiden; Oberhettinger (transform tables); Miller and MacCamy (Volterra integral equations). | CIT (hypotheses **not** obtained) | The model is a third-order Lur'e loop with three real lags and a sigmoid, i.e. the bounded-nonlinearity cousin of WH. Hypotheses unknown; the 7-page length and the Volterra-equation references suggest a restricted (steep or special) regime. Worth obtaining; cannot be relied on. |
| Bingxi Li, Nonlinear Anal. 5 (1981) 931–958 | "Periodic orbits of autonomous ODEs: theory and applications" (torus principle; cited in 2025 as a survey containing uniqueness conditions). | CIT | Not accessed. |
| L.-S. Dai, J. Differential Equations 31 (1979) 392; L. A. Sánchez, J. Math. Anal. Appl. 317 (2006) | Existence, uniqueness and global asymptotic stability of the periodic solution of the modified Michaelis–Menten mechanism (3-D); later dynamics paper on the same system. | CIT | Hypotheses not obtained. Dai follows Hastings (parameter-restricted). |
| R. A. Smith, Proc. Roy. Soc. Edinburgh A 85 (1980); J. Differential Equations 69 (1987) 265; Burkin, Differ. Equ. 51 (2015) 1717 (survey); Miranda-Villatoro, Forni, Sepulchre, arXiv:1710.01645 | Quadratic-cone (frequency-domain, shifted circle criterion) hypothesis gives a 2-D reduction and Poincaré–Bendixson; with an unstable equilibrium, existence of an orbitally stable closed trajectory. The Leonov–Burkin school constructs systems with a unique equilibrium and an arbitrary number of stable cycles. (mem + SN) | SN / CIT | Existence only; needs slope bounds on `e^U` over the attractor, which degrade as `rho` grows. No uniqueness theorem. |
| Muldowney, Rocky Mountain J. Math. 20 (1990); Li & Muldowney, J. Differential Equations 124 (1996) 425; SIAM J. Math. Anal. 27 (1996) 1070; Math. Biosci. 125 (1995) 155 (SEIR) | A periodic orbit is orbitally asymptotically stable with asymptotic phase if the second compound equation `z' = (df/dx)^{[2]}(p(t)) z` is asymptotically stable; an omega-limit set without equilibria is a stable hyperbolic periodic orbit iff the second compound equation is uniformly asymptotically stable; Bendixson criterion `q2 < 0` with a matrix weight `A(x)` and Lozinskii measure. SEIR proofs use state-dependent weights such as `E/I` and the vanishing of averaged log-derivatives. | SN + FT of Smith 2012 (which applies it) | This is brief item (e). It is a criterion, not a theorem whose hypotheses WH is known to satisfy. No published instance was found in which it proves uniqueness of a cycle around an unstable equilibrium in a 3-D competitive system. |
| Katz, Giordano, Margaliot, arXiv:2407.00461 (2025); J. Differential Equations (2025), arXiv:2503.09155; Taghvafard, Proskurnikov, Cao (2016, Goodwin–Smith) | Strongly 2-cooperative systems (include negative cyclic feedback): strong PB property; explicit sets of initial conditions converging to periodic orbits. Conclusion of arXiv:2503.09155: "An interesting topic for further research is to find conditions that guarantee uniqueness of the periodic orbit". Goodwin–Smith: uniqueness of periodic orbits "seems to be an open problem" (snippet). | FT (2503.09155 conclusion) / SN | Independent confirmation that, as of 2025, no general uniqueness theorem exists for this class. |
| Hofbauer & So, Appl. Math. Lett. 7 (1994); later 3-D competitive Lotka–Volterra examples with 2–4 limit cycles | 3-D competitive systems can have several limit cycles. (mem) | CIT | Competitiveness + dimension 3 cannot by themselves give uniqueness. |
| Sepulchre & Stan, Systems Control Lett. 54 (2005) 809; Stan & Sepulchre, IEEE TAC 52 (2007) 256 | Lur'e feedback of a **passive** linear system with a static sector nonlinearity: (almost) global convergence to a unique limit cycle for the gain slightly above the Hopf (or pitchfork) value. | SN | Passivity fails (relative degree 3). Conceptual prior art for a global Hopf collar only. |
| Llibre, Ponce and co-authors (3-D piecewise-linear systems); Megretski / Gonçalves–Megretski–Dahleh (relay feedback) | Uniqueness/global stability of limit cycles for 3-D piecewise-linear and relay systems under structural conditions (e.g. "Hopf-scroll" sewing systems; surface Lyapunov functions via LMIs). | SN | Piecewise-linear only. |

---------------------------------------------------------------------------------------------------------------

## 3. Statements formulated and tested in this lane

Python: `/Users/admin/.venvs/analysis-py313/bin/python <script>` from this directory. All runs are single-process,
floating point (DOP853 via `whlib`), each under 2 minutes. Output of each script is saved next to it as `*.out`.

### 3.1 PROVED — star-shaped and convex projections; two exact exponent formulas

Notation: on a cycle, `L2 := D^2 + (mu - eta) D + zeta` (the second factor of manuscript eq. (2.9)); `L2[p] = 0`
(eq. (2.8)); `K > 0` (Cor. 3.2); `W > 0` (Lemma 3.3).

**Lemma A (Wronskian gauge identity).** For every `C^2` function `phi(t)` along a cycle put
`S_phi := p phi' - p' phi`. Then
`S_phi' = (eta - mu) S_phi + p L2[phi]`.

Proof. `S_phi' = p phi'' - p'' phi`. Substitute `p'' = (eta - mu) p' - zeta p`:
`S_phi' = p phi'' - (eta-mu) p' phi + zeta p phi = (eta-mu)(p phi' - p' phi) + p[phi'' - (eta-mu) phi' + zeta phi]`. ∎

Two instances:
* `phi = U`: `S_U = p^2 - U p' =: S*`, and `L2[U] = p' + (mu - eta) p + zeta U = k - y - eta p + zeta U = K`.
  Hence `S*' = (eta - mu) S* + p K`.
* `phi = -p'`: `S = p'^2 - p p'' =: V`, and differentiating `L2[p] = 0` gives
  `L2[p'] = eta' p' - zeta' p = -W`, so `L2[-p'] = W`. Hence `V' = (eta - mu) V + p W`.

**Theorem B.** On every nonconstant positive cycle (any `mu > 0`, `k > 0` at which one exists), for all `t`:
(i) `S* = p^2 - U p' > 0` (the projected cycle in the `(U,p)` plane is strictly star-shaped about `E = (0,0)`,
    clockwise);
(ii) `V = p'^2 - p p'' > 0` (the projected cycle has nowhere-vanishing curvature; with the Jordan property of
    Thm 4.1 it is a strictly convex curve).

Proof. By Thm 4.1 `p` has exactly two zeros per period, both simple; they split the period into one open arc with
`p < 0` (crest to trough) and one with `p > 0`. At a zero of `p`: `V = p'^2 > 0`; and `S* = -U p' = U (y - k) > 0`,
because by eq. (4.1) `x - k = (1 + eta)(y - k)`, so `U` and `y - k` have the same nonzero sign. Let `F` be `S*` or `V`
and `G` be `K` or `W` (positive). At any zero of `F`, `F' = p G`. On the closed `p < 0` arc `F` is positive at both
ends; if it had a zero, take the last one `t*`: `F > 0` just after, so `F'(t*) >= 0`, but `F'(t*) = p G < 0`. On the
`p > 0` arc take the first zero: `F > 0` just before, so `F'(t*) <= 0`, but `F'(t*) = p G > 0`. Contradiction in both
cases. ∎

**Corollary C (exact exponent formulas).** For all `a, b >= 0`, not both zero,
`lambda_c = <eta> - mu = - < p (a K + b W) / (a S* + b V) >`,
in particular `lambda_c = -<p K / S*> = -<p W / V>`; equivalently, with `dU = p dt`,
`lambda_c T = - oint (K/S*) dU = - oint (W/V) dU`.

Proof. `S := a S* + b V > 0` is periodic and `(log S)' = (eta - mu) + p (aK + bW)/S`; average over a period. ∎

**Corollary D (half-cycle bounds).** Let `I_+`, `I_-` be the integrals of `eta - mu` over the rising (`p > 0`) and
falling (`p < 0`) arcs, and `r := (y_crest - k)/(k - y_trough)`. Then `I_+ < 2 log r` and `I_- > -2 log r`
(from `(V e^{-int(eta-mu)})' = p W e^{-int(eta-mu)}` and `V = (y - k)^2` at the two zeros of `p`). The analogous
bounds with `S*` use `U (y - k)` at crest and trough.

Status: **proved**, conditional only on the v2.0 results quoted (Prop. 2.2, Cor. 3.2, Lemma 3.3, Thm 4.1), which
the brief instructs to take as given. Relation to the literature: `V` is the Gasull–Giacomini orbit-curvature
function and `S*` the angular function of the planar reduced field; their Dulac defects are `M_1 = pW` and `pK`.
Neither formula decides the sign of `lambda_c`, because of the factor `p`.

Numerical checks (all pass):
* `curvature_dulac_test.py` — 21 cycles (`mu in {0.05,0.2,1}`, `rho` 1.05–8): `V > 0` on all; residual of
  `V' = (eta-mu)V + pW` at most `6e-10` (relative); `-<pW/V>` equals `<eta> - mu` to all printed digits.
* `curvature_dulac_test_mu_gt1.out` — `mu in {3,10}`, `rho in {1.1,2,5}`: same.
* `wronskian_gauge_test.py` — 20 cycles (`mu in {0.05,0.3,1,4}`, `rho` 1.05–8): `S* > 0`, `V > 0`; both
  identities to `<= 1.3e-9`; four formulas (`(a,b) = (1,0),(0,1),(1,1),(1,10)`) agree with `<eta> - mu`
  (worst gap `1e-5` is against the monodromy value at `m_t ~ 1e-6`, a known cancellation in `whlib`).
* `branch_compare_test.py` — Corollary D holds on 18 cycles. Observed: `I_+ < 0 < I_-` on every tested cycle
  (the rising arc contracts, the falling arc expands, the sum is negative).

### 3.2 REFUTED — Massera-type and curvature-Dulac hypotheses for the reduced planar field

Reduced field (manuscript Prop. 5.2; brief (b)): `X = (P,Q) = (v, k - mu v - h(u,v))`, with `h_u = zeta`,
`h_v = -eta` on cycles.

* **Radial angular monotonicity (Sabatini–Villari Thm 6).**
  `alpha = P(uQ_u + vQ_v) - Q(uP_u + vP_v) = v[h - k - u h_u - v h_v] = -v K`. Since `K > 0`, `alpha` has the sign of
  `-v`: it changes sign on every cycle. Not applicable. (`starshaped_test.py` confirms numerically.)
* **Strict star-shapedness after Conti–Filippov (Sabatini, arXiv:1006.2020, Thm 3)** with `g(u) = h(u,0) - k`,
  `Phi(u,v) = mu + (h(u,v) - h(u,0))/v`. `starshaped_test.py` (interpolating `h(u,0)`, `zeta(u,0)` from 134–400
  crossings of `W^u(E)` with `p = 0`): on the cycle `S` ranges over `[-0.25, 0.54]` (`mu=1, rho=1.1`),
  `[-0.36, 1.11]` (`mu=1, rho=1.3`), `[-0.21, 0.58]` (`mu=0.3, rho=1.1`), negative throughout the quadrant
  `u>0, v<0`, positive throughout `u<0, v>0`. Gross sign change, far above interpolation error. Not applicable.
  Structural reason: near `E`, `S ~ u Phi_u/sqrt(zeta_E) + v Phi_v`, a nonzero linear form, because
  `Phi_v(E) = h_vv(E)/2 > 0` (at Hopf `c02 = mu/D_mu`) and `Phi_u(E) = -eta_u != 0`. Massera-type theorems require
  the damping to have a critical point at the equilibrium; here it is monotone through it.
* **Gasull–Giacomini curvature Dulac function (arXiv:2101.03874, Thm 2.3).** With `P = v`,
  `M_1 = v[h_uu v^2 + 2 h_uv v Q + h_vv Q^2] = v * (Phi^T D^2h Phi) = p W`. Changes sign with `p`. Not applicable.
* **Pointwise branch comparison.** `lambda_c T = -int [ (W/V)_rising(U) - (W/V)_falling(U) ] dU`. The integrand is
  **not** pointwise positive: negative on 46%–95% of the `U`-range on 18 cycles (`branch_compare_test.py`).
  A pointwise-in-`U` comparison of the two branches cannot work.

### 3.3 REFUTED (hypothesis C) / SUPPORTED (hypothesis M) — H. L. Smith's monotone-concave theorem for the crest map

Crest map `P` on `{p = 0, p' < 0}`, charts `xi = ((e^{aU}-1)/a, ((y/k)^b - 1)/b)` centred at `E`.

* `concavity_grid_test.py`, 9x9 grids over `[0.1,1.5] x [0.1,1.5]` times the crest coordinates, at
  `(mu,rho) = (1,2), (1,1.1), (0.2,1.5)`: all 243 return derivatives are entrywise positive and all returns stay in
  `{U>0, y>k}` (hypothesis (M) supported; consistent with the compound-system argument of Prop. 4.4 whenever
  `x >= y > k`). Hypothesis (C) (entries of `DP` non-increasing along ordered pairs) is violated in 96–281 of 832
  comparisons for every chart tried (`(a,b) = (1,1), (a*,1), (1/2,1/2), (0,1), (0,0)`), with relative increases up
  to 3x, concentrated away from the diagonal. Along the single ray from `E` through the crest, (C) is nearly but
  not robustly satisfied (0–10 violations of 68). **The theorem's hypotheses fail on any order interval.**
* `sublinear_chart_test.py`: the necessary fixed-point consequence `M v <= v`, `v_i = phi_i/phi_i'`, holds on all
  24 computed cycles (`mu in {0.2,1,5}`, `rho` 1.02–8) in the chart `(a,b) = (0.75,0.75)` and in `(1,1)`, but this
  is only because `rho >= 1.02`: from the manuscript's Hopf expansions the `O(h)` direction mismatch of a power
  chart is
  `r_v - r_e = h [ (1-b)/(2A^2) + a/(2A) - (mu^2+4mu+1)/(2 A D_mu) ] + O(h^2)`, `D_mu = mu^2+6mu+1`,
  (reproduces the manuscript's `mu(1-mu)(mu+3)/(4A^2 D_mu)` at `a=b=1/2`), so near Hopf `M v <= v` fails unless
  `a = (mu^2+4mu+1)/D_mu - (1-b)/A`. Charts must be `mu`-dependent even for the necessary condition.

### 3.4 REFUTED — Zhang's theorem on the formal slow-manifold Liénard reduction

Idea tested: for `mu -> 0` at fixed `rho`, slave `z = y/k` to `e^U` and obtain, in `tau = sqrt(k) t`,
`U_tautau + (mu/eps - eps e^U) U_tau + e^U - 1 = O(eps^2)`, `eps = sqrt(k)`. This truncated Liénard equation
satisfies Zhang's hypothesis for every `rho > 1` (`d/dU[(1/rho' - e^U)/(e^U - 1)] = e^U(1 - 1/rho')/(e^U-1)^2 > 0`).
Outcome (`toda_ratio.py`):
* On level sets of `H = q^2/2 + e^U - U - 1` one has `oint (e^U - 1) q^2 dtau = -oint q^2 dq = 0` identically, so
  the first-order Melnikov function is `eps (1 - mu/k) oint q^2 dtau > 0` on every level: the truncated equation
  has **no** limit cycle when `rho(1+mu) > 1` ("at most one" is satisfied vacuously).
* Actual WH cycles at `mu in {0.2, 0.05, 0.02}`, `rho in {1.5,3,8}` are far from level sets (`spread(H)/mean(H)`
  between 1.1 and 7; `mean H` of order `1/mu` or larger). The slaving assumption fails on the fast pulses.
The exact identity `<e^U p^2>/<p^2> = 1/rho` cannot hold on near-level curves (where the ratio is 1). A
near-Hamiltonian description can only be valid in a wedge `rho - 1 = O(mu)`; that regime, not fixed `rho`, is where
Abelian-integral monotonicity (Li–Zhang 1996) could be used.

---------------------------------------------------------------------------------------------------------------

## 4. Assessment of brief items against the literature

* (a) pointwise criterion — same logic as Xie's Theorem (stable => unique) and Sansone's argument; sound given
  Smith 2012 Thm 1.1(c) (the complement of the 1-D stable manifold in the open octant is connected).
* (b) attractor as a graph — prior art: Walther (DDE, graph with possibly several nested cycles), Ivanov–Lani-Wayda
  (2-D invariant manifold spiralling to a bounding cycle), Hirsch/MPS (limit sets only). Convexity: Baigent-type
  mechanism; no prior statement for WH found. None of these implies uniqueness.
* (d) Dulac formulation — published form: Gasull–Giacomini extended Bendixson–Dulac (`M_s` criterion); 3-D
  invariant-surface form: Busenberg–van den Driessche (J. Math. Anal. Appl. 172 (1993) 463, (mem)). The two most
  natural `V` (angular `S*`, curvature `V`) give defects `pK`, `pW`: positive gauges times the sign-changing `p`.
* (e) supersolution/compound criterion — it is the Li–Muldowney method; Smith 2012 used it with `phi = 1` below
  Hopf. SEIR-type proofs (state-dependent weights + vanishing averaged log-derivatives) are the methodological
  template for brief (c).
* (i) distributed-delay form — confirmed (Ruan eq. (9.11) for `mu = 1`).
* (k) Liénard plane — no published Zhang-type theorem for non-separable `P(U,s)` was found.
* No brief item was contradicted by the literature.

---------------------------------------------------------------------------------------------------------------

## 5. Sharpest obstruction

The nonlinearity `e^U - 1` is sublinear for `U < 0` and superlinear for `U > 0`, and the effective damping
`mu - eta` is monotone through the equilibrium rather than extremal there. Every structural uniqueness theorem
located needs one of: two-sided sublinearity (Cao 1996 for delay equations; Krasnosel'skii/Smith concavity), a
separable Liénard form (Zhang), a damping critical point at the equilibrium (Massera, Sabatini), a step-like
nonlinearity (Hastings, Glass–Pasternack, Farcot–Gouzé), or a small parameter. For exactly this nonlinearity in
the delay setting (Wright's equation) the uniqueness question stayed open for 57 years and was closed only by
asymptotics plus computer-assisted proof; Xie's thesis records that the stability argument cannot reach the Hopf
point. In this lane the same obstruction appears concretely as the factor `p` in every Dulac defect computed
(`alpha = -vK`, `M_1 = pW`), and as the failure of hypothesis (C) for the crest map.

---------------------------------------------------------------------------------------------------------------

## 6. What to try next (ordered)

1. Use Corollary C as a design space: find a periodic `phi` built from state functions with `L2[phi] > 0` and the
   right boundary signs, for which `<p L2[phi]/S_phi>` can be signed (for instance through a second gauge whose
   defect carries the opposite half-cycle weighting). The observed pattern `I_+ < 0 < I_-` suggests comparing the
   rising arc against the falling arc at equal `y` or equal `p'` rather than equal `U`.
2. Li–Muldowney weights: search for a matrix weight (not diagonal) on the compound system, singular at `E`, using
   the zero-mean log-derivatives of brief (c) together with `log S*` and `log V` (two new zero-mean
   log-derivatives: `<(log S*)'> = <(log V)'> = 0`).
3. Jones-type programme: explicit collar by the van den Berg–Lessard–Queirolo desingularisation; validated crest
   map (CAPD) on compact `(mu, rho)` boxes using `M v <= v` in a `mu`-tuned chart or direct enclosure of `spr M`;
   Xie/Kolesov–Rozov-type asymptotics for `rho -> infinity`; separate analysis of `mu -> 0`.
4. Obtain the two papers of Bingxi Li (1981) and Dai (1979), whose hypotheses could not be read here.

---------------------------------------------------------------------------------------------------------------

## 7. Search log and access limits

* WebSearch calls used: 26 of 45. WebFetch: 8. OpenAlex, Semantic Scholar graph, Crossref and arXiv APIs via curl
  or urllib (no API key, no personal data sent). Semantic Scholar `snippet_search` was not used.
* Full texts read: arXiv 1109.5204, 2202.11034, 1801.09806, 1705.02432, 1704.00029, 0909.4107 (intro), 1101.2761,
  math/0307372, 1006.2020, 1305.3402, 2101.03874, 0903.1969 (theorem section), 2503.09155 (abstract, conclusion),
  2407.00461 (abstract); Ruan's survey (author-hosted PDF, §9); MATCH 94 (2025) 863 (to identify the system).
* Not accessible (HTTP 403 / paywall; not circumvented): Hastings 1977 beyond the first page; Li 1981 (both);
  Dai 1979; R. A. Smith 1987; H. L. Smith 1986 (both); Zhu–Smith 1994; Xie 1991/1993; Li–Muldowney 1996;
  Mallet-Paret–Smith 1990; Gedeon–Mischaikow 1995; Walther 1995; Zeng–Zhang–Gao 1994; Xiao–Zhang 2003;
  Kuang–Freedman 1988. Statements for these are from abstracts, snippets, quoting papers, or memory, as marked.
* Scripts in this directory: `toda_ratio.py`, `sublinear_chart_test.py`, `concavity_grid_test.py`,
  `starshaped_test.py`, `curvature_dulac_test.py`, `branch_compare_test.py`, `wronskian_gauge_test.py`
  (+ `arx.py`, `oa.py`, `oaw.py`, `oaabs.py` API helpers). Outputs: matching `*.out` files.
