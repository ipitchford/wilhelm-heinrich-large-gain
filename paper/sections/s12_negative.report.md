# Report for `s12_negative.tex` (Section "Certificates that fail, and what a proof must do")

Date: 2026-10-03. ROOT = `/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03`, BUNDLE = `ROOT/v3/wilhelm_heinrich_v3_0`.

## Files written

| file | purpose |
|---|---|
| `BUNDLE/paper/sections/s12_negative.tex` | the section (526 source lines, 19 pages in the wrapper) |
| `BUNDLE/paper/bib/s12_negative.tex` | three new `\bibitem`s: `FG10`, `Sabatini10`, `Duff53` |
| `BUNDLE/paper/test_s12_negative.tex` | compile wrapper (kept; auxiliary files deleted) |
| `ROOT/v3/checks/s12_negative_check.py`, `.out` | formula re-checks: 43 exact (sympy) + 113 floating-point (whlib), all pass |
| `ROOT/v3/checks/s12_iota_check.py`, `.out` | size of the non-monotonicity of the one-turn integral at 7 parameter pairs (recomputed with the cartography sampler, read-only import); source of the example quoted after Table `neg:tab:candidates` |
| this report | |

Compile test: `pdflatex test_s12_negative.tex` (twice) ends with exit status 0, **no LaTeX errors**, no overfull boxes.
The only warnings are undefined cross-section references and citations, as expected for a single-section wrapper.
A scratch assembly (outside the bundle) of all sections that exist today (s01–s12, s14, appA, three bib fragments)
resolves every reference and citation made by this section except `sec:evidence` (not yet written), and reports no
multiply defined label. Part numbering of the cited statements (`prop:kappa` (i)–(vi), `prop:strongidentity`
(i)–(v), `thm:gauges` (i)–(vi), `prop:crestaxis` (i)–(iii), `att:cor:dulac` (i)–(ii), `prop:crestmonotone` (a)–(d),
`thm:attractor` (A1)–(A8)) was re-checked against the current files at 13:55.

## Outcome reported

No statement in the section claims or implies that the uniqueness conjecture is solved. Every refuted candidate is
labelled "refuted (fp)" = refuted by floating-point computation unless a proof of failure is given; survivors are
labelled "survives (fp); unproved".

## Statements, labels and sources

Exported labels: `sec:negative thm:failed eq:radial eq:radialfailure prop:nodulac prop:nomonomial lem:rigidity`.
All other labels carry the prefix `neg:`.

| statement (label) | status in paper | source |
|---|---|---|
| Theorem `thm:failed`, eqs `eq:radial`, `eq:radialfailure` (square-root radial certificate fails near Hopf, 0<mu<1) | proved; proof ported verbatim | v2.0 `manuscript.tex` Section 8, lines 527–580 |
| Proposition `neg:prop:charts` (exact two-component identities; first-order test for any vector; power charts (a,b); second component fails for mu>1; linear chart fails for every mu) | proved (new) | chart coefficient: `lanes/literature/NOTES.md` §3.3; second-component formula: `lanes/freeform/NOTES.md` §3.3 and `lanes/monotone/NOTES.md` X2; identities and proof written here |
| paragraph "Relation to concave monotone maps" (H. L. Smith's theorem as quoted in Farcot–Gouzé; its consequence M v <= v fails) | discussion + proved consequence | `lanes/literature/NOTES.md` §2.3 (last rows), §3.3; `raw/0903.1969.txt` (Theorem 1 re-read) |
| Proposition `prop:nodulac` (i) forced sign, (ii) no Dulac function regular at E (extended to C^1 functions of (U,p,y,eta,zeta)), (iii) forced singularity e^Lambda >= c|q|^{-theta} for every theta<2, (iv) local singular weights Q_E^{-c}, c>1, exist | proved; (i) and the planar case of (ii) by citation of `att:cor:dulac` | `lanes/dulac/NOTES.md` Prop. B1, §§2–4; `lanes/dulac/VERIFY.md` §3.4; `lanes/planar/NOTES.md` §1.3 item 3; s07 `att:cor:dulac` |
| Proposition `prop:nomonomial` (no supersolution M(t)(c1 l, c2 y, c3 x); quantitative bound min_i <g_i> <= -(gamma_E - A); equality case) | proved; first part ported, bound new | `lanes/freeform/NOTES.md` L5 |
| one-line exclusions of constant pairs and of (eta_0, zeta_0) | proved | `lanes/supersolution/NOTES.md` §2 table, `VERIFY.md` Claim 4 last bullet |
| Lemma `lem:rigidity` (i) expansion of cycles near Hopf, (ii) expansion of the adjoint ratios, (iii) rigidity of trial pairs, (iv) state functions, eq. `neg:eq:jet` | proved; (iii)–(iv) ported, (i)–(ii) new | `lanes/supersolution/NOTES.md` Lemma A1; `VERIFY.md` Claim 4 |
| Remark `neg:rem:jet`, Remark `neg:rem:selfcontained` | proved remarks | new |
| list "Floating-point results for explicit families" | floating-point evidence | `lanes/supersolution/NOTES.md` §2, `table_A_out.txt` |
| Proposition `neg:prop:planar` (i) F∧(DF q) = -pK, (ii) curvature/angular Dulac defects pW, pK, (iii) Liénard plane not separable near Hopf | proved; (iii) new | `lanes/literature/NOTES.md` §§2.4, 3.2; `lanes/planar/NOTES.md` §2(b); s05 Remark (c); s07 `att:prop:scaling` |
| paragraphs on Massera/Sabatini–Villari, Sabatini, Gasull–Giacomini, Zhang, rotated fields | discussion | `lanes/literature/NOTES.md` §§2.4, 3.2–3.4; `lanes/planar/NOTES.md` §2, `VERIFY.md` §3 |
| Table `neg:tab:candidates` (32 rows: A1–A8, B1–B14, C1–C4, D1–D6) | status per row | `lanes/cartography/NOTES.md` (S1–S6), planar, dulac, literature, monotone, supersolution notes; every count re-read from the output file named in the row |
| Proposition `neg:prop:signchange` (one sign change of the one-turn integral iota, or strict monotonicity of the one-turn mean, implies uniqueness; hyperbolic attraction when iota(xi_out)<0) | proved, conditional (hypotheses unproved) | new; uses `prop:crestaxis`, `cor:criterion`, (A2), (A4) |
| Section `neg:sec:requirements` (what a proof must accomplish; Wright's equation and Jones' conjecture) | discussion | `lanes/literature/NOTES.md` §2.2, §5; s05 last subsection; lanes' "sharpest obstruction" paragraphs |

Results already in other sections are cited, not re-proved: `att:prop:scaling`, `att:cor:dulac`, `prop:crestaxis`,
`cor:criterion`, `thm:attractor` (A2),(A4)–(A6),(A8), `prop:kappa`, `prop:strongidentity`, `prop:lvchain`,
`thm:gauges`, `prop:crestmonotone`, `prop:crestmap`, `prop:riccatientire`, `lem:projectivebounds`, `lem:collar`.

## VERIFY corrections applied

1. `dulac/VERIFY.md` §3.4 item 1 (the cycle-average argument for "only R <= 0 is possible" needs a hyperbolic
   attracting cycle). `prop:nodulac`(i) rests on `att:cor:dulac`(ii) (flux/return-map argument on the inner side of
   the innermost cycle, no hyperbolicity). A paragraph after the proposition states why the average identity alone
   is insufficient.
2. §3.4 item 2 (dependence on continuity of eta on the attractor at E). Now a proved input: `thm:attractor`
   (A2), (A4), (A5), (A6) are cited at the point of use.
3. §3.4 item 3 ("B must blow up at least like r^{-2}" was argued only for pure powers). Replaced by the proved
   statement (iii): for every theta<2, e^Lambda >= c|q|^{-theta} near E. The exponent 2 itself needs
   eta - eta_E = O(|q|), i.e. C^2 regularity of the graph at E, which s07 does not prove; the text says so.
   Part (iv) records that weights Q_E^{-c}, c>1, do work locally (planar notes §1.3 item 3).
4. §3.4 item 4 (reporting detail, "-2.000 sigma exactly"): the only LP margin quoted is that of the gauge family
   (`lp_G.out`, -2.0000/-2.0001), described as "to four digits".
5. `planar/VERIFY.md` §3 (bracket sign valid only for 0 <= theta <= 1). The scaling identity is cited from
   `att:prop:scaling`, which already has the corrected two-sided statement; this section uses only the
   infinitesimal form F∧(DF q) = -pK and the corrected geometric sentence.
6. `supersolution/VERIFY.md` §4 item 4 ("reduces to an upper bound on the falling arc" is a sufficient route, not
   an equivalence): the last subsection says "a proof by way of these formulas … needs", followed by the explicit
   sentence that this is a sufficient route and not a statement about every conceivable proof.
7. `supersolution/VERIFY.md` Claim 4 scope remark: the lemma is claimed only for trial pairs with constant
   order-zero limits (plus j-dependent constants sigma_j, tau_j, which the lemma shows to be o(xi_j)), and row C3 /
   the last subsection say "C^1 functions of the state". The extension to phase-dependent order-zero terms
   sketched by the referee is **not** included.
8. `monotone/VERIFY.md` §3.2 (one-directional relation between lambda_c<0 on every cycle and uniqueness): the
   section only ever says that uniqueness follows from lambda_c<0 on every cycle.

## Formulas re-checked after renaming (recorded in `ROOT/v3/checks/s12_negative_check.out`)

Renamings: crest amplitude h -> xi (h is the graph function in v3.0); H=y-k written out; M -> M_{(U,z)};
k_H -> \kH; lane symbols V,Q,S1,S4 -> \Sg,\Xg; q=lambda-p -> \ell; sigma=Re lambda_u -> alpha_E; S -> <m>.

Exact (sympy), all zero residual: C_mu = A^2+mu; invariance of the quadratic centre graph (5 coefficients);
r_e; r_v - r_e for the square-root chart and for the general chart (a,b); limits a=b=1/2 and (0,1); stable vector;
M_0 derived from the linear system by projection, trace 1+q_f, det q_f, eigenvector (1,1/A); coefficient of
`eq:radialfailure` = (M_0)_12 × mismatch; second-component coefficient; a_*, b_* solve the linearised Riccati
system and equal the gradient of the quadratic centre graph; first-order profiles solve the linear system at Hopf;
Riccati perturbation equations; first-order defects B_1,B_2,B_3 and the vanishing positive combination;
D Q_E = 2 alpha_E Q_E; F∧(DF q) = -pK; defects of M(c1 l, c2 y, c3 x); (A+theta)(mu+theta)(1+theta) = chi(A+theta);
row B13 identity (Filippov-plane radial function); lambda_c = -A mu^2 xi^2/(2 C_mu D_mu) from the l_1 of Appendix A.

Floating point (whlib): near-Hopf cycles at mu = 0.3, 1, 3 with rho-1 = 1e-5…4e-4: the two exact component
identities; (r_v-r_e)/xi; both leading coefficients; M_{(U,z)} -> M_0 at rate O(xi); r_e = zeta(crest)/k;
Perron slope; least-squares first-order coefficients of eta, zeta against (U,p) (agree with `neg:eq:astar` to
1e-4); lambda_c/xi^2; k<e^U p^2> = k_H<p^2> to 1e-12. Signs observed: first component >0 at mu=0.3 for
xi = 0.014, 0.044 and <0 at xi = 0.086; second component >0 at mu=3; both <0 at mu=1 (consistent with the text).
Lotka–Volterra means and the monomial bound on 3 cycles with 2000 random (c1,c2,c3) each; limits
eta-mu -> 2 alpha_E and D log Q_E -> 2 alpha_E on W^u(E).

## New bibliography entries (verified 2026-10-03)

`FG10` (Crossref DOI 10.1080/00207720903144552; preprint arXiv:0903.1969 verified through the arXiv API; the
theorem number refers to the preprint, which is what was read), `Sabatini10` (Crossref DOI
10.1016/j.chaos.2010.07.002; arXiv:1006.2020), `Duff53` (Crossref DOI 10.2307/1969724). Keys used from
`bib/s01_intro.tex`: Smith12, SmithHL86, SV10, GG21, Massera54, Zhang86, Xie93, JLM17, vdBJ18, Jaquette19.

## Evidence files named in Table `neg:tab:candidates` (for the evidence section)

Relative to `provenance/research_lanes/`: `cartography/stests_results.json`, `cartography/stests_summary.csv`,
`cartography/s3a_cycle_scan.json`, `cartography/cycle_tests.csv`; `planar/job2_a.json`, `job2_b.json`,
`job2_c.json`, `job2_d.json`, `job2_massera.json`; `dulac/lp_G.out`, `lp_run2.out`, `lp_shared.out`,
`explore2.out`, `explore4.out`, `quadrant_facts.out`; `literature/concavity_grid_test.out`,
`sublinear_chart_test.out`, `starshaped_test.out`, `branch_compare_test.out`, `toda_ratio.out`;
`monotone/retmap.out`, `twelve.out`; `supersolution/table_A_out.txt`.
Relative to `provenance/`: `v2_0_bundle/evidence/manuscript_summary.json` (row A1).

## Differences between lane notes and their output files (the output files were followed)

- planar notes: "beta positive on 21–53 % in every chart" — `job2_massera.json` gives 28 %–52 % (`frac_beta_pos`).
- planar notes: LP optimum "+0.07…+0.15" — `job2_c.json` gives 0.031…0.152 (still positive in 8 of 8 runs).
- dulac notes: quadrant sign facts "on 30 cycles" — `quadrant_facts.out` lists 20 cycles (row D4 says 20).
- monotone notes X4 "false for rho >= 1.3" — `twelve.out`: fails on 10 of 15 cycles (row D2).
- cartography S1: the one-turn *integral* iota fails to be monotone at 19 of 62 parameter pairs
  (`stests_results.json`, field `integral_monotone_decreasing_in_U_inside`), so the second hypothesis of
  `prop:crestaxis`(iii) is refuted by floating-point computation (row A6); the first hypothesis at 9 of 62 (row A5).
  What survives is the one-turn *mean* (row A8).

## Not ported, and points that deserve a second reading

1. **Theorem `thm:failed`** is ported verbatim and, as in v2.0, relies on the classical smooth expansion of the
   Hopf branch (k = k_H + O(xi^2), m_t = 1 + O(xi^2), M_{(U,z)} = M_0 + O(xi)), which the paper does not prove.
   Remark `neg:rem:selfcontained` derives the sign statements with o(xi) remainders from `lem:rigidity`(i),(ii);
   the remark's convergence M_{(U,z)} -> M_0 is argued in three lines (continuity of the variational flow and of
   the normalised flow direction at the crest) and should be read critically.
2. **`lem:rigidity`(i),(ii) are new proofs** (the lane assumed the expansions). (i) uses Arzelà–Ascoli on rescaled
   cycles, the zero count of `thm:geometry` and the new identity k<e^U p^2> = k_H<p^2>; (ii) uses
   `prop:riccatientire`(i) for the limit and a variation-of-constants estimate. They passed the numerical checks
   above but have had no independent adversarial reading.
3. The hypothesis of `lem:rigidity` is "cycles converging to E as k_j -> k_H". The sentence that every sequence of
   cycles with k_j -> k_H has this property refers to the proof of `lem:collar` (s09) and is marked "not used".
4. `prop:nodulac`(iii) gives every exponent theta<2, not theta=2 (see correction 3).
5. `neg:prop:signchange` and `neg:prop:charts`, `neg:prop:planar`(iii) and the quantitative part of
   `prop:nomonomial` are new statements written for this section; proofs are short and complete in the text.
6. Massera54, Zhang86 and Duff53 were not obtained; the text says so and describes only the general form of their
   hypotheses. H. L. Smith's theorem is quoted from the Farcot–Gouzé preprint.
7. References by description rather than by label (internal labels of other sections avoided): "the last remark of
   Section `sec:attractor`" (currently `att:rem:planes`) and "the distributed-delay form of Section `sec:model`"
   (currently `eq:delay`). The two internal labels named in the assignment, `att:prop:scaling` and `att:cor:dulac`,
   are referenced directly.
8. Lane items deliberately left out of the table: the bound kappa < 0.6 (refuted, side remark), "<eta(lambda-p)> <= k"
   and "lambda_c <= c_H <eta p>" (their refutation rests on unproved near-Hopf asymptotics), the LP leads of
   `planar/job2_c5/c6` (ill-conditioned, not validated), and the fold-direction programme of the monotone lane.
9. Observation for the frame editor (not my file): `s14_status.tex` item (a) writes log m_t ≈ -2π ε^3 r with
   ε = sqrt(mu)/(1+mu), r = (rho-1)(1+mu)^2/mu. The normal form with Appendix A gives
   -2πσ(rho-1)/(1+σ^2) = -2π ε (rho-1)/(1+ε^2), so the s14 expression is its leading term for small ε only.

## Cross-model review (work/sol/review_s12_response.md)

Verdicts returned: HOLDS for `thm:failed`, `neg:prop:charts`, `prop:nodulac`, `prop:nomonomial`, `neg:prop:planar`, `neg:prop:signchange`; GAP for `lem:rigidity` and for several prose overclaims. The editing pass was cut off before writing this record; the orchestrator completed the last repairs and checked the rest against the current text.

| Point | Disposition |
|---|---|
| Opening claim on all comparison-type certificates | Repaired: restricted to the classes covered by `lem:rigidity`. |
| "failure confined to a neighbourhood of Hopf" | Repaired: the proposition proves failure near Hopf; no confinement claimed. |
| Smith's concave-map theorem quoted second-hand | Repaired: the quoting source is named. |
| Sign of a Dulac function | Repaired: restated for a function defined near the innermost cycle (forced sign at interior points). |
| Closed ellipse and exit time in `prop:nodulac`(iii) | Repaired. |
| Equality case in `prop:nomonomial` | Repaired: identities stated. |
| `lem:rigidity`(i) period convergence | Repaired: two consecutive negative-slope zeros on a fixed interval, C^1 convergence, one crest per period. |
| `lem:rigidity`(ii) uniqueness of the entire Riccati solution | Repaired: Proposition `prop:riccatientire` cited; its hypotheses cover the constant driving function. Locally uniform versus global bounds distinguished. |
| `lem:rigidity`(iii) means over [0, T_j] | Repaired: explicit estimate. |
| Remark on the Hopf-amplitude formula | Repaired by restatement: the formula is now "expected" from the normal form, the conversion is not carried out, it is not used in any proof, and its agreement with two computed cycles is reported (floating point). |
| Self-contained remark recovers sign and leading coefficient only | Repaired in wording. |
| `neg:prop:planar`(iii) interval statement; phase notation | Repaired. |
| Table row A2 coordinates; row B8 qualifier | Repaired. |
| `neg:prop:signchange`(ii) endpoint case; (iii) order of the proof | Repaired. |
| Return time and shrinking turn as u -> 0 | Repaired: argument from the focus (angular speed bounded below near the origin, Gronwall bound on the radius). |
| "Every exact formula ..." | Repaired: the gauge formulas are named. |
| "A proof needs an upper bound" | Repaired: "would need", for the routes examined. |

The section compiles without errors.
