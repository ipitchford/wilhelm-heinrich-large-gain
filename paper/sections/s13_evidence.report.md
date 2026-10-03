# s13_evidence.tex — report of the evidence engineer (3 October 2026)

Section label `sec:evidence`; internal labels `num:sec:executed`, `num:sec:classes`, `num:sec:results`,
`num:sec:cartography`, `num:sec:review`, `num:sec:failure`, `num:sec:missing`, `num:tab:counts`,
`num:tab:lgconstants`, `num:tab:cartography`. Compile test: `paper/test_s13_evidence.tex` (no LaTeX error, no overfull
box; undefined references are the cross-section ones). All 56 labels referenced exist in the section files as of
this report.

## 1. Where every number comes from

No number in the section is typed from memory. Three sources:

1. Macros `\Vthree...` (115) and the nine framework macros, defined in `paper/generated/v3_macros.tex` and
   `evidence_macros.tex`, written by `code/make_tables.py` from the receipts. `evidence/v3/v3_summary.json` lists, for each
   macro, its TeX value and the receipt field it is taken from. `code/check_consistency.py` rebuilds every generated
   file from the receipts and compares it with the file on disk, and checks that every macro used in any section is
   defined. The section loads the macro files itself when compiled alone (`\ifdefined` guard), so it works both in the
   wrapper and in `manuscript.tex`.
2. Design parameters of the checks (parameter lists, tolerances, counts of grid points), which are constants of the
   scripts under `code/v3/` and are quoted from them.
3. Numbers of the three retained failures in `num:sec:failure` (49 of 50; 65 of 66, 1.8e-9, 775 k, 11 236 of 137 303;
   36 of 37, 3.26e9, 3e-13, 7e-12). They are static receipts (`evidence/exact_checks_initial_failure.json`,
   `evidence/v3/v3_gauges_initial_failure.json`, `evidence/v3/v3_largegain_lemmas_small_mu_supplement.json` and its log);
   `check_consistency.py` asserts the counts.

Tables: `generated/v3_counts.tex`, `generated/v3_largegain_constants.tex` (notation of Table `lg:tab:constants` of
s10 as it stood when this was written: k_0, L_*, varphi_*, chi_ii, tilde R, vartheta_1, chi_0, hat R, hat vartheta_N,
hat w, c_log), `generated/v3_cartography.tex`.

## 2. Statements and their sources

| Statement | Source |
|---|---|
| counts table, runtime, memory | receipts; `evidence/run_receipt.json` |
| framework block (31 cycles, audits, cross-method, reciprocal, envelope) | v2.0 manuscript §9 (condensed), receipts `evidence/*.json` re-executed in this bundle |
| list of what the 148 exact identities cover; the two cofactors | `evidence/v3/v3_exact_identities.json` (groups, `cofactor_of_O`) |
| high-precision constants | `v3_largegain_constants.json`; wrapped scripts `lib/lg_constants_author.py`, `lib/lg_constants_referee.py` |
| gauges, weak adjoint, attractor, large gain, wedge, refutations | the receipts of the same names |
| cartography statements (i)–(iv), recomputation | `v3_cartography.json`, shipped tables, lane `cartography/NOTES.md` |
| adversarial-review record | `work/lanes/*/VERIFY.md`, `work/attack_results.json`, `work/sol/README.md` and responses, section report files |
| failure semantics | `code/run_all.py`, `v3_controls.json`, the retained receipts |
| missingness | this report, `provenance/research_lanes/SKIPPED.md` |

## 3. What was changed relative to the state found (the previous engineer was interrupted)

Found: v2.0 scripts copied into `code/`, `code/v3/common.py`, `cyclekit.py`, `v3_exact_identities.py`, `v3_gauges.py`,
15 library copies, two receipts. The v2.0 suite and `v3_exact_identities` ran and passed as found. `v3_gauges` had one
failing record (65/66).

- `v3_gauges`: the failing predicate was `cycle_located` at (mu, rho) = (0.02, 20). The crest return is reproduced to
  2.5e-12 relative; the gate compared the y-residual divided by k with 1e-9, and y_crest = 775 k there. The residual is
  now relative to the crest state. The failing receipt is kept (`v3_gauges_initial_failure.json`). Every mathematical
  predicate passed in both versions.
- `cyclekit`: for mu > 1 the gauge W is evaluated in the cancellation-free variables eps = eta − (mu−1),
  r = y − (mu−1)(lambda−p) (identity checked by sympy at each run). Before, the sign of W was unresolved at 11 236 of
  137 303 samples (all mu ≥ 3); now at none. Optional extra integrands were added to the profile integrator.
- New: `v3_largegain_constants`, `v3_largegain_lemmas`, `v3_attractor`, `v3_weak_adjoint`, `v3_wedge`,
  `v3_cartography`, `v3_negative`, `v3_controls`; `run_all.py`, `make_tables.py`, `check_consistency.py`;
  library copies of the referee's large-gain scripts, of the mu > 1 supporting-plane script, of the period-identity
  script and of the cartography code.

## 4. Discrepancies with other sections found while building the evidence (for the editor)

1. **s05_gauges.tex**, sentence after the proof of Theorem `thm:gauges`: "eta_1 takes negative values on part of every
   cycle examined". On the 66 cycles of `v3_gauges` eta_1 is negative somewhere on 37 and positive throughout on 29
   (the small-amplitude cycles; eta_1 → mu > 0 at the Hopf point). eta_0 is negative somewhere on 10, all with
   mu ≤ 0.1, which agrees with the sentence. s13 states the receipt values; the sentence in s05 should be changed to
   "on part of every cycle of sufficiently large amplitude examined" or to the counts.
2. **s10_largegain.tex**, Remark `lg:rem:scope`(a): log m_t at (mu, k) = (1, 300) is printed as −375.7. The computed
   value is −375.648, which rounds to −375.6 (the research output printed −375.65). s13 prints −375.6 from the receipt.
   The check in `v3_largegain_lemmas` accepts the printed value only because of its 20 per cent allowance on half a
   unit of the last digit.
3. **s11/s12 cross-model reviews**: `work/sol/review_s11_response.md` and `review_s12_response.md` are dated after the
   report files of s11 and s12. s13 therefore says of these two reviews only that their findings are in the review
   directory. If the dispositions are applied and recorded, the two sentences in `num:sec:review` can say so.
4. The digest's "within 4e-4" for the selection rule at amplitude ≤ 1, mu = 1e-4 is 4.03e-4 at amplitude 1; s13 prints
   4.0e-4 from the receipt and does not use the phrase.

No discrepancy was found between receipts and: Table `lg:tab:constants` (s10), Tables `wed:tab:selection` and
`wed:tab:fixedrho` (s11; all rows recomputed in `v3_wedge`), rows B13 and D1 of Table `neg:tab:candidates` (126 of 300
cycles; 300 of 300 cycles), the refutation parameters (0.0027, 5.02) and (1, 5) of s12.

## 5. Formulas re-checked

The section introduces no new formula. Two identities are used by the new code and are checked by sympy inside it:
the 8-monomial form of W in the variables (x, q, zeta, eta, eps, r) (`cyclekit.w_identity_check`, run at the start of
`v3_gauges`), and the three quotient identities behind Proposition `prop:nomonomial` (`v3_negative`, exact records).
The expression used for eta_1, eta_1 = (x Sigma + k psi(U) p')/Xi − 1, follows from `gau:Xistate`.

## 6. Not done

- No interval arithmetic.
- The structural tables `stests_summary.csv`, `stests_results.json`, `cycle_tests.csv`, `limit_map.csv`,
  `s3a_cycle_scan.json`, `dipscan.json` are shipped and checked for the statements read off them; they are not recomputed.
- Research outputs cited by s11 and s12 other than those wrapped (`fixed_rho` is recomputed; `v4_theoremC`,
  `relax_map`, `lambda_of_q`, `variance_route`, the LP searches of the Dulac line, `job2_*`, `retmap`, `twelve`,
  `table_A`, the literature tests) are shipped under `provenance/research_lanes/` and not re-executed.
- The author's lemma script at (mu, k) = (0.02, 1000) is outside the default run (432 s, one unresolved comparison;
  receipt kept).
