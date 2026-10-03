# s07_attractor.tex — audit report (2026-10-03)

The original writer was interrupted before self-checks and wrote no report. This report is the auditing
editor's. It covers the fidelity audit against the verified source, the disposition of the 14 defects of the
cross-model review, the re-checks of every displayed identity, and the compile test.

Abbreviations: SRC = `ROOT/work/lanes/planar/attractor_theorem.tex` (line numbers); PV = `ROOT/work/lanes/planar/VERIFY.md`;
DV = `ROOT/work/lanes/dulac/VERIFY.md`; SOL = `ROOT/work/sol/review_s07_response.md`.

## 0. Outcome

- No theorem, proposition or corollary had to be weakened. The statements of `thm:attractor` and `cor:criterion` are
  unchanged. No `\label` was renamed or removed (49 labels, 7 exported).
- **One status sentence was wrong by omission and was corrected.** The section said only that the two hypotheses of
  Proposition `prop:crestaxis`(iii) are unproved and pointed to floating-point evidence "bearing on them". The
  evidence is negative at some parameters: in `ROOT/work/lanes/cartography/stests_results.json` (62 parameter pairs)
  `iota` is not decreasing at 19 pairs and `R_1(u)/u` is not decreasing at 9 pairs (all with very small or very large
  `mu`; the 9 are among the 19; both hold at the other 43). The editor recounted these from the JSON; they are rows
  A5, A6 of the table in `s12_negative.tex`. Part (iii) is a correct implication, but its hypotheses are refuted in
  floating point on part of the parameter domain. The text now says so, before the proposition and after its proof,
  and points to Section `sec:negative`.
- Every proof step survived being written out. The 14 cross-model defects were each checked against the text:
  12 led to a repair (a corrected statement of an imported result, or an added argument); 2 (defects 8 and 10)
  were not gaps in the text as it stood, and there the wording was made explicit. Defect 1 (the transversal lemma)
  is a real defect of statement: the term "transversal" was undefined; with the classical definition the lemma is
  the standard one and both uses are valid.
- All displayed identities were re-derived (53 exact sympy checks) and spot-checked in floating point (54 checks on
  three parameter pairs); all pass.
- Compiles with the wrapper; only cross-section references and citations are undefined.

Limits, stated plainly: (a) the imported results (I2)–(I4) are cited by chapter of Hartman's book; the book was not
available to the editor, so theorem numbers and the printed wording of the transversal lemma were not re-opened
(PV says the same). The definition of a transversal now written in (I2)(b) is the classical one supplied in the
editing brief. (b) The remark for `k <= k_H` imports Smith (2012, Thm 1.1(b)) and Boros–Hofbauer (2023, Thm 13);
these were cross-checked against `ROOT/work/lanes/literature/NOTES.md` only. (c) The numerical spot checks are
ordinary floating point at `rho <= 1.5` and are not part of any proof. (d) The repairs were written and checked by
the same editor; the repaired text has not been read again by the second model or by an independent referee.

## 1. Statements included and their sources

| Statement in s07 | Label | Source |
|---|---|---|
| Setting; Riccati system; normal equation | `att:eq:ric`, `att:eq:normal` | SRC §1, l. 36–58 |
| Imported results (I1)–(I5) | — | SRC §1.1, l. 60–91 (source item (I6), dependence on v2.0, became the sentence after the list) |
| Theorem (A1)–(A8) | `thm:attractor` | SRC Thm 1.1, l. 94–172 |
| Remark for `k <= k_H` | — | SRC l. 173–174 |
| Lemma, explicit absorbing box | `lem:absorbing` | SRC Lemma 2.1, l. 183–227 |
| Proof of (A1) | — | SRC l. 229–260 |
| Linear system, Birkhoff constants | `att:eq:lin`, `att:eq:birkhoff` | SRC l. 262–278 |
| Lemma, uniform positivity and contraction | `att:lem:cone` | SRC Lemma 3.1, l. 280–302 |
| Lemma, suprema of bounded C^1 functions | `att:lem:sup` | SRC Lemma 3.2, l. 304–312 |
| Proposition, unique positive entire Riccati solution | `prop:riccatientire` | SRC Prop. 3.3, l. 314–366 |
| Remark on the size of the rate | — | PV §1 "(A2) Riccati", remark |
| Proof of (A2) | — | SRC l. 368–386 |
| Theorem, two-orbit comparison | `thm:twoorbit` | SRC Thm 4.1, l. 389–409 |
| Proof of (A3) | — | SRC l. 410–413 |
| Lemma, gauges on the attractor | `att:lem:gauges` | SRC Lemma 4.2, l. 415–434 (residual `-eta' O` from PV §1) |
| Proofs of (A4), (A5) | `att:eq:jet` | SRC l. 437–459 |
| Lemma, reduced flow on the u-axis | `att:lem:axis` | SRC Remark 5.1, l. 463–466 (promoted to a lemma) |
| Remark, gauges in terms of h | `att:rem:convexgauges` | SRC Remark 5.1, l. 460–463 |
| Proposition, radial scaling | `att:prop:scaling` | planar `NOTES.md` §2(b), l. 137–140, in the corrected form of PV §3 |
| Linearisation at E, `E^u = n_E^perp` | `att:eq:Eu` | SRC l. 469–481 |
| Proof of (A6) | — | SRC l. 483–514 |
| Lemma, geometry of projected cycles | `att:lem:cyclegeom` | SRC Lemma 6.1, l. 516–535 |
| Lemma, compactness of the crest set | `att:lem:crests` | SRC Lemma 6.2, l. 538–551 |
| Proof of (A7) | — | SRC l. 556–576 |
| Multipliers, `n(a_0)^perp` | `att:eq:nperp` | SRC l. 578–589 |
| Proof of (A8) | `att:eq:dulacmean` | SRC l. 591–641; step (d) with the detail of PV §1 "(A8)" |
| Corollary, parameter-wise criterion | `cor:criterion` | (i), (iii): SRC (A8)(c),(d). (ii): new packaging; the attraction statement is proved in the section from Lemma `lem:absorbing`, (A1), (A6) and (I3) |
| Remark, what is and is not proved | `att:rem:status` | SRC l. 643–650 |
| Proposition, crest return map | `prop:crestaxis` | SRC Prop. 8.1, l. 655–686; continuity of `T_1` at the end point is added and proved |
| Corollary, two constraints on Dulac functions | `att:cor:dulac` | dulac `NOTES.md` Prop. B1, l. 126–136, in the repaired form of DV §3.4 |
| Remark, divergence in other planes | `att:rem:planes` | SRC Remark 8.2, l. 691–709 |

Not ported, by design: SRC §9 (numerical consistency checks) and the numerical sentence after Prop. 8.1; the section
points to Section `sec:evidence` instead. The clause "B must blow up at least like r^{-2}" of the dulac lane's
Prop. B1(ii) is not in this section (DV §3.4 item 3 found it under-argued; Section `sec:negative` treats it).

Renaming used: source `\cD` → `\Dom`; `D` → `\Dsc`; `q = lambda - p` → `\ell`; `psi_t` → `\chi_t`; `J`, `J_c` (curves) →
`\mathcal J`, `\mathcal J_\xi`; crest `c` → `\xi`; `G` → `\Omega`; `delta_x, M_x, ...` → `\underline x, \overline x, ...`; `c(t)` → `\beta(t)`;
`m_0, M_0, Delta_0, tau_0` → `\underline m, \overline m, d_0, \nu_0`; `Sigma, C(t,s)` → `\triangle, \triangle(t,s)`;
`Psi` (adjoint solution) → `\tilde n`; `Sigma_3, P_3` → `\mathcal T, \Pi`; `tau(u)` → `T_1(u)`; `g~, G~, R~` → `g_*, G_*, R_*`;
`m, s` (Remark 8.2) → `w_1, w_2`; points `a, b, z` → bold.

## 2. VERIFY corrections applied

Planar VERIFY (PV):
1. Claim 5 (scaling identity): the sign statement holds for `0 <= theta <= 1` only. Proposition `att:prop:scaling`
   states `b_theta >= 0` for `0 <= theta <= 1`, the lower bound `(1-theta) K(theta q)` for `0 < theta < 1`,
   `theta q` in D \ {0}, and the reversed sign `b_theta <= (1-theta) K(q) < 0` for `theta > 1`. It uses the convex
   extension `H`, so `theta q` in D is assumed only where it is needed. The geometric sentence is PV's corrected one.
   PV writes the lower bound with strict inequality; the proof gives `>=` and the section states `>=`.
2. (A8)(d), return-time detail: present in full (a nearby cycle meets `{p = 0, U > 0}` only at its own crest, by
   Lemma `att:lem:cyclegeom` and injectivity of `pi`, so the return selected by the implicit function theorem is
   the return to its own crest).
3. (A2), remark that the rate (pull-back estimate) is quantitatively weak: stated after Proposition `prop:riccatientire`.
4. Prop. 8.1(i): monotonicity of `R_1` is derived from (ii) (`R_1' > 0`), as PV recommends. Prop. 8.1(iii): only
   positivity of `R_*` is used (PV: `R~(0+) >= 0` suffices); the source's parenthetical on the unstable focus is
   not needed and is not ported.
5. Lemma gauges (iii): the residual `-eta' O` off the constraint is stated.
6. "Convex C^1 graph" is stated as PV requires: `h = H|_D` with `H` convex, `h` in `C^1(int D)`, `(h, g)` satisfying
   the Whitney hypothesis; D is not claimed convex.

Dulac VERIFY §3.4 (DV):
1. B1(i): the cycle-average argument needs a hyperbolic cycle. Applied twice: (A8)(f) claims the impossibility of the
   reverse inequality only in the form the average identity proves (on `Gamma_in` or `Gamma_out`, unless it vanishes
   a.e.), and says that `= 0` a.e. gives `lambda_c = 0`, which is not excluded; Corollary `att:cor:dulac`(ii) is the
   flux repair proposed by DV, in return-map form, valid whether or not the cycle is hyperbolic.
2. B1(ii): the dependence on the then unproved continuity of `eta` at E is removed, because (A2), (A4) are proved in
   this section; Corollary `att:cor:dulac`(i).
3. The under-argued clause "at least like r^{-2}" is not stated here.

## 3. Cross-model review (SOL): disposition of the 14 defects

| # | Defect | Checked against the text | Disposition |
|---|---|---|---|
| 1 | (I2)(b) false for an arbitrary "transversal segment" | Correct as a defect of statement: the term was undefined, and the counterexample (a segment through the equilibrium of a rotation) fits the loose reading. With the classical definition the lemma is the standard one | **Repaired.** (I2)(b) now defines a transversal (closed segment on which F neither vanishes nor is parallel to it, hence crossed in one direction), states the monotonicity lemma for three intersection times, and says why the loose reading fails. Both uses now verify the definition (5 and 11) |
| 2 | (A1)(a): induction covers only starting points in B' | Correct | **Repaired.** Induction written out with times `s_j`; forward invariance of `B_1` derived from it for points `phi_r(b)` |
| 3 | (A1)(a): finite cover alone does not give uniform absorption | Correct | **Repaired.** `t_Q` is the largest of the finitely many entry times; each point stays in `B_1` after its own entry |
| 4 | (A1)(b): nestedness not used explicitly | Correct (minor) | **Repaired.** Monotonicity of `t -> phi_t(B_1)` proved and used |
| 5 | (A6): transversal through `q_*` | Correct | **Repaired.** `I` is a short closed segment with midpoint `q_*`, perpendicular to `F(q_*)`, shown to be a transversal; isolation of crossings and crossings near `q_*` are derived from the sign of the normal component. The step "the orbit meets `I` on one side of `q_*`" was dropped: the conclusion (a strictly monotone bi-infinite sequence cannot have the same limit at both ends) does not need it |
| 6 | Lemma cyclegeom: crossing count | Correct (minor) | **Repaired.** Finiteness (transverse crossings are isolated, period compact) and the crossing-number formula (branch of the angle on the slit plane) are written out |
| 7 | (A7): null-homotopy fact not imported | Correct (minor) | **Repaired.** Added to (I4) with its one-line reason; cited in the proof of (A7) |
| 8 | (A8)(d): `Pi(a_j)` in `{U > 0}` | The text already shrank `T'` so that `U > 0` on `Pi(T')` | **Wording made explicit** (continuity of `Pi`, `Pi(a_xi) = a_xi`, `xi > 0`). No mathematical change |
| 9 | (iii): strict decrease of `iota` versus strict concavity of `R_*` | Correct: the converse was asserted | **Repaired.** Both directions proved (mean value theorem; a concave differentiable function has nonincreasing derivative, and equal values at two points force affinity between them) |
| 10 | (iii): strict decrease of `R_*(eps)/eps` | The implication in the text is valid: a nonincreasing function with equal values at two points is constant between them, so `R_* = c eps` there | **No change of mathematics needed; argument written in full**, with the three-point concavity inequality displayed |
| 11 | Cor. dulac(ii): (I2b) on the whole half-axis | Correct as stated before the repair of defect 1 | **Repaired.** By Lemma `att:lem:axis`, `F = (0, -g_*(s))`, `g_*(s) > 0` on `0 < s <= xi_out`, so every closed subsegment is a transversal; the monotone sequence accumulates at `xi_in`, hence increases |
| 12 | Cor. dulac(ii): citation for continuity at `xi_in` | Correct (minor) | **Repaired.** Two cases: `xi_in < xi_out` (the C^1 statement of Prop. crestaxis(i)) and `xi_in = xi_out` (the end-point statement) |
| 13 | Remark planes: existence of the implicit functions | Correct (minor) | **Repaired.** The two coercivity inequalities for `H` are derived from the affine functions in its definition; the four maps are bijections of R; C^1 on the images of int D by the implicit function theorem |
| 14 | Opening sentence ambiguous | Correct (minor) | **Repaired.** "(A1)–(A5) import only the Birkhoff–Hopf theorem (I1); (A6)–(A8) use in addition (I2)–(I5)" |

Verdict lines of SOL: Theorem `thm:attractor` "GAP" and Corollary `att:cor:dulac` "GAP" rested on defect 1 and its
consequences; Proposition `prop:crestaxis` "GAP" on defects 9–10. All are closed by the repairs above without
changing a statement.

## 4. Fidelity audit: discrepancies found and fixed

Statement by statement and step by step against SRC. No sign or index error from the renaming was found. Changes:

1. Proposition `prop:riccatientire`, Bounds: the source clause on the derivatives of `eta, zeta` had been dropped;
   restored as "their derivatives are bounded" (the source wording was ambiguous; the clause is not used later).
2. (A8)(b): the points `b_j` lie on the local stable manifold "for all large j" (as (I3b) gives), not for all j.
3. (A8)(c): `W^u_loc(Gamma_out)` is taken inside the open set `Dom` before Lemma `lem:absorbing` is applied; (A1)
   is cited for the conclusion.
4. (A8)(a): the conclusion for `Gamma_out` is stated (the multiplier of `J_out` for the extended field is `m_t`,
   whatever the extension), as PV notes.
5. Lemma `lem:absorbing`: the statement orders the constants strictly; neither SRC nor the section proved it. One
   sentence added (`beta_1 < 1/mu`, so `overline x > 2k > k/4 > underline x`, and the other two follow).
6. (A7): joint continuity of `(t, q) -> Theta_t(q)` stated, since it is what makes the family of loops a homotopy.
7. Proposition `prop:crestaxis`(i): the neighbourhoods `N_0, N_1` are taken inside `{u > 0}`, so the zero found is
   a return to the half-axis.
8. Corollary `att:cor:dulac`(ii): the segment on which `P_Lambda` is defined is closed, so the integral is defined
   from its left end.
9. Introduction of the section: "divergence `eta - mu`" qualified by "on the interior of D".
10. Notation (PLAN rule on global symbols): `delta` is the global `k - k_H`; its four local uses were renamed
    (small time → `vartheta`; small energy → `eps'`; tangent vector in `att:eq:Eu` and variational solution → bold `v`,
    added to the local-symbols paragraph). `D` is the global operator d/dt; `D Pi` → `d Pi` (as `d pi` already in the
    section) and `D^2 h` → `nabla^2 h`.
11. Three displays wider than the text block (`beta(t)`, `att:eq:returnderivative`, the four `D log` identities) were
    broken into two lines; content and label unchanged. "note that" removed (no second person).
12. Status of Proposition `prop:crestaxis`(iii): see §0. The source (SRC l. 687–689) reported `iota` decreasing at 14
    pairs with `rho <= 3`; the later cartography computations found failures at extreme `mu`. The sentence before
    the proposition and the paragraph after its proof were rewritten accordingly; the reference `sec:negative` was
    added.

Note for the owners of other sections (not edited here): `s14_status.tex` cites Proposition `prop:crestaxis` for the
statement that a single sign change of the one-turn integral suffices; that statement is the proposition on one sign
change in `s12_negative.tex` (internal label `neg:prop:signchange`), which builds on `prop:crestaxis`(i),(ii).

Checked and left unchanged: (A8)(f) differs from SRC ("the reverse inequality is impossible") and is the form that
its proof and DV §3.4 support; the remark `att:rem:status` replaces SRC's "equivalent up to the degenerate case" by
the two one-way statements; Lemma `att:lem:crests`, the degree argument of (A7), and step (b) of (A8) follow SRC
step for step.

## 5. Formulas re-checked

Script `ROOT/work/framework_checks/s07_attractor.check.py`, output `s07_attractor.check.out` (single process,
about 5 s; `python = /Users/admin/.venvs/analysis-py313/bin/python`). Result: 107 checks, all pass.

Part A, exact (sympy), 53 identities:
- `J^T n`, and the three components of `n' + J^T n + (1+eta) n` (`att:eq:ric` ⇔ `att:eq:normal`);
- `att:eq:xyl` (`x' = (lambda - ell) x`, `ell' = y - mu ell`);
- `H' = B H` ⇔ Riccati, both directions; `(e^{int(1+eta)} n)' = -J^T(...)`;
- `nu_0 = tanh(d_0/4) = (M - m)/(M + m)`;
- `att:eq:deltaode`; `O' = -(1+eta) O`; `K' + (1+eta) K = k psi(U)` (identically); `x R(-U) = k psi(U)`;
  `W' + (1+eta) W - x p^2 = -eta' O`; the displayed `zeta''`, `eta''`, `p''`;
- `att:eq:axis` from `att:eq:pde` at `v = 0`; `K = k - [h - u h_u - v h_v]`;
- `att:eq:scaling`, and the algebra of both bounds in Proposition `att:prop:scaling`;
- characteristic polynomial of `J_E`, its root `-gamma_E`, value `k_H` of the cubic at `A`, `J_E e_E = -gamma_E e_E`,
  `J_E^T n_E = -gamma_E n_E`, the constant Riccati solution at E, `eta_E - mu = gamma_E - A`, exponent of `m_t`;
- wedge transport `w' = (div F) w` and the two boundary evaluations of `att:eq:returnderivative`;
- Remark `att:rem:planes`: `w_1'`, `w_2'`, four Jacobians, four divergences (by implicit differentiation), four
  `D log` identities, and divergence `= (eta - mu) + D log(Jacobian)` in each plane;
- Lemma `lem:absorbing`: the change of variables in `beta(t)` and the convolution identity `s e^{-mu_0 s}`.

Checked by hand in addition (not in the script): every step of Lemma `lem:absorbing`; the factor count in Lemma
`att:lem:cone`(ii); the order of the four bounds in Proposition `prop:riccatientire`; the two-sided inequality
`att:eq:whitney`; convexity and concavity of the four implicit functions (epigraph and hypograph as sublevel sets);
the orientation statement after Proposition `att:prop:scaling`; the stable eigenvector `e_s + r f` of the return map.

Part B, floating point (whlib, DOP853), at `(mu, rho) = (1, 1.5), (0.3, 1.5), (3, 1.3)`, 18 checks each:
`O = 0`, `K, W, Lambda_E > 0` and cone exclusion on the cycle; the bounds `att:eq:constants`; `lambda_c` from the five
cycle averages of Remark `att:rem:planes` against the monodromy (agreement 1e-12); `M^T n = m_s n`, `n.f = 0`,
`n.e_s != 0`, multipliers `{1, m_s, m_t}` and `det M = e^{-AT}`; the crest identity; the upper constants of the
absorbing box; the Birkhoff constants and the pull-back estimate from data (0,0), (50,0), (0,50), (50,50);
on an orbit of `W^u(E)`: the axis identity at every crest (1e-12), increasing crests (`att:eq:outward`),
supporting planes between crest points of different turns, both scaling inequalities on the crest half-axis for
`theta < 1` and `theta > 1` (586 to 1748 ordered pairs), and `att:eq:returnderivative` against central differences
between neighbouring orbits (relative error at most 2e-7).

## 6. Compile test and bibliography

- `pdflatex -interaction=nonstopmode -halt-on-error test_s07_attractor.tex`, run twice in `BUNDLE/paper`: exit 0,
  18 pages, no LaTeX error, no overfull box. Undefined references are the cross-section labels `eq:logsystem`,
  `eq:adjointmatrices`, `eq:projective`, `eq:ratios`, `eq:cones`, `eq:crestidentities`, `lem:projectivebounds`,
  `prop:multipliers`, `prop:crestmap`, `thm:support`, `thm:geometry`, `sec:adjoint`, `sec:geometry`, `sec:negative`,
  `sec:evidence` (all exported by their sections in the PLAN) and the citations (the wrapper has no bibliography).
  Auxiliary files and the PDF were deleted; the wrapper is kept.
- Citation keys used: `Birkhoff57 Bushell73 Seneta81 Hartman02 Thomassen92 Whitney34` (all in
  `BUNDLE/paper/bib/s07_attractor.tex`) and the global keys `Smith12`, `BH23` (present in `BUNDLE/paper/bibliography.tex`).
  The six DOIs of the fragment were resolved through the Crossref API with `ROOT/work/bibcheck/bibcheck.py`;
  authors, year, title, journal, volume and first page agree. The header comment of the fragment records this.
- Cross-section consistency checked by reading each citing sentence: the statements cited from this section by
  `s01`, `s03`, `s06`, `s09`, `s10`, `s11`, `s12`, `s14` (`cor:criterion`(i)–(iii), `thm:attractor`(A2), (A4)–(A6),
  (A8)(f), `prop:riccatientire`(i), `prop:crestaxis`, `att:cor:dulac`(i), (ii)) read as those sections use them, with
  the one imprecision in `s14` noted at the end of §4. The quantities referred to in
  `s02`–`s04` (`eq:adjointmatrices`, `eq:ratios`, `eq:projective`, `lem:projectivebounds`, `eq:cones`,
  `eq:crestidentities`, `thm:support`, `prop:multipliers`, `thm:geometry`, `prop:crestmap`) are the ones this section
  says they are.
