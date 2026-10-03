# Report — `s05_gauges.tex` (Section 5: convex star-shaped projections, exact exponent formulas, scalar reduction, further identities)

Written by the auditing editor on 3 October 2026. The section text was drafted by a writer who was interrupted before
self-checks and before a report; the text was then refereed by a second model (`ROOT/work/sol/review_s05_s06_response.md`)
and repaired by the orchestrator. This report records a statement-by-statement and step-by-step comparison of the final text
with the lane sources, the repairs made at the audit, and the checks run.

ROOT = `/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03`. Lane files are `ROOT/work/lanes/<lane>/NOTES.md` and `VERIFY.md`.

## 1. Statements included and their sources

| Label | Content | Source | Change relative to the source |
|---|---|---|---|
| `gau:wronskian` (Lemma) | `V_{f,g}' = (eta-mu)V_{f,g} + f L2[g] - g L2[f]`; `L2[1]=zeta, L2[U]=K, L2[p]=0, L2[-p']=W` | literature §3.1 Lemma A (l. 151–162); dulac §1 "General principle" (l. 108–111) | stated for a general pair `(f,g)`; the lanes' `L2[y-k]=W` is the same value because `-p' = mu p + (y-k)` and `L2[p]=0` |
| `thm:gauges` (i) | `Sigma>0`, `Xi>0`, `eq:gaugeodes`, state forms | dulac Thm D(a),(b) (l. 37–42, 65–81); supersolution Thm R 1–2 (l. 29–32, 40–51); literature Thm B | proof = integrating-factor two-arc argument of dulac/supersolution |
| `thm:gauges` (ii) | star shape, polar angle decreases by exactly 2 pi | dulac D(a); supersolution R1 | the count "exactly 2 pi" is written out (the lanes quoted the winding number) |
| `thm:gauges` (iii) | curvature, tangent angle, strict supporting lines `gau:supportline`, hull | dulac D(a) proof (l. 82–85) | the supporting-line inequality and the hull statement are proved here; the lanes cited the standard fact for curves of nonvanishing curvature |
| `thm:gauges` (iv) | `eq:exactexponent`, mixed formula `gau:mixedexponent` | dulac D(c); supersolution R3; literature Cor. C | weights `a,b>=0` of Cor. C normalised to `vartheta in [0,1]` |
| `thm:gauges` (v) | shifted centres `Sigma_c` | supersolution R4 (l. 38, 53–54) | none |
| `thm:gauges` (vi) | `eta_0`, `eta_1`, means `mu`, `gau:etadifference` | dulac D(b),(c) (l. 41–50, 87–92) | none |
| `gau:halfcycle` (Cor.) | half-cycle bounds | literature §3.1 Cor. D | both gauges written out |
| `gau:remarks` | log-concavity; general Wronskian formula; relation to the Bendixson–Dulac criterion of [GG21, Thm 2.3] | supersolution "Corollaries" first bullet; dulac "General principle"; literature §3.1 last paragraph and §2.4 | none |
| `prop:kappa` (i),(ii) | affine parametrisation by `kappa`, Riccati equation `eq:kappariccati` | dulac D(c) (l. 50–52, 90–92); supersolution "Reductions/Constraint line" (l. 169–179) | lane `N=K/S1` is `kappa` |
| `prop:kappa` (iii) | compound form `gau:compoundsource`, defect form | dulac D(e) (l. 58–62, 105–106) | none |
| `prop:kappa` (iv) | branch form of the exponent | dulac §6 item 1; freeform §3.9; literature §3.2 last bullet | identity only; the refuted pointwise sign is not asserted |
| `prop:kappa` (v),(vi) | certificates by state functions; supporting planes bound `kappa` from below | dulac §4 O1, O2 (l. 222–241) | stated on cycles only (VERIFY §3.5) |
| `gau:dynamicbounds` (Remark) | `zeta<y`; `eta<ell` for `mu<=1/2`; resulting upper bounds for `kappa` | dulac O2 (l. 236–239); supersolution "Corollaries" | none |
| `prop:polar` | polar system, orbit equation, polar exponent formula | supersolution R3 and proof of Prop. T (l. 91–94); freeform L6 | none |
| `prop:nestedradii` | radial comparison of nested cycles, shifted centre | supersolution Prop. T (l. 76–99) | shifted centre proved directly with `K+c zeta>0` (supersolution VERIFY, Claim 2) instead of the convex envelope |
| `prop:monotoneintegrals` | `<(alpha p+beta p')(alpha eta+beta zeta)> < 0`, orderings | dulac D(d) (l. 54–56, 94–103) | the sign of the chord parameter is obtained from `gau:supportline` instead of the lane's continuity argument |
| `gau:convexcriterion` (Cor.) | convex criterion | dulac Prop. CD (l. 113–124) | wording of the sufficient conditions corrected (VERIFY §3.3) |
| `prop:holder` | Hölder sandwich | freeform L1 (l. 16–36) | AM–GM in place of the generalised Hölder inequality; equality case written out |
| `prop:lvchain` | Lotka–Volterra chain; Cole–Hopf identity; second-order forms; chain pair | freeform L2, L3 (l. 38–66); supersolution "explicit almost-solutions" (`omega^Y`) | none |
| `prop:strongidentity` (i)–(iv) | pairing with `nu=e^{At}(c,b,a)`; two instances; scalar forms; `frak r>1`, `frak s>(1-mu) frak r` | supersolution Prop. B (l. 151–165) | (iv) is stated for all `mu>0`, the lane for `mu<=1`; see §5 |
| `prop:strongidentity` (v) | comparison bound, independence of `phi`, convexity | supersolution Lemma A0 (l. 106–127) | identification of the compound spectrum written out |
| `prop:gaussmap` | injectivity of the Gauss map, orientation | freeform L4 (l. 78–87) | orientation argument written out |
| closing subsection | what the identities do not decide | freeform §3.6; supersolution §5.3; dulac O1–O2 | none |

Not ported, with reasons (none is used by another section): the identity `Sigma = Q_E + U Lambda_E` and the mean identity
`<(p-alpha_E U)Lambda_E/Q_E> = alpha_E` (supersolution "Corollaries", last bullet; exact and verified in the lane, omitted by
the writer); the general wedge rule `f^v` and the defects of `f^e_U` and of the chemical radial field; the Markov form and
the integral (next-generation) form of supersolution "Reductions", and the equivalence stated after freeform L3: these rest on
an imported threshold theorem that the lanes did not reprove. The scalar comparison for the Riccati equation is marked "only
sketched" in the lane and is correctly absent. Lemma A1, Prop. B1 and freeform L5 belong to Section `sec:negative`.

## 2. VERIFY corrections and their status

- dulac VERIFY §3.3 (Prop. CD: the monotonicity conditions need a set containing the coordinate-wise paths): applied in the paragraph after `gau:convexcriterion`.
- dulac VERIFY §3.5 (O1–O2 proved on cycles only): applied; `prop:kappa`(v),(vi) are statements on a cycle.
- dulac VERIFY §1.1 (`eta<lambda-p` fails for `mu>1/2`; `K/V<0.6` false): the restriction `mu<=1/2` is in `gau:dynamicbounds`; the false bound is not in the text. Reconfirmed at the audit: `eta<ell` fails at `(mu,rho)=(1,2)` and `(3,1.5)`.
- supersolution VERIFY, Claim 5 (positivity of `(a,b,c)` was assumed): the text cites Lemma `lem:strongstable`.
- supersolution VERIFY, Claim 2 (shifted centre): direct proof used.
- supersolution VERIFY §4.4 ("reduces to an upper bound" is a sufficient route, not an equivalence): the closing subsection says "a proof by these means would need".
- dulac VERIFY §3.4 (Prop. B1) and supersolution VERIFY Claim 4 (Lemma A1) concern Section `sec:negative`, not this file.

## 3. Repairs made after the second-model review (checked at the audit)

The review found only wording defects in this file. The orchestrator's repairs were checked and are correct and complete:
`prop:polar` (the sign statement now refers to the fraction under the integral sign); `prop:holder` (`rho=k/kH` is named with a
reference to `eq:parameters`); `prop:lvchain`(ii) (only the second and third rows of the adjoint system are asserted; the
first row is stated as not assumed); `prop:gaussmap` (the orientation argument states that the curve lies below its
horizontal tangent at the top point and that the interior is on the right of the velocity).

## 4. Discrepancies found at the audit and fixed

1. **Notation clash.** `eta_1` is the comparison function of `thm:gauges`(vi); `prop:gaussmap` used `(eta_1,zeta_1)`, `(eta_2,zeta_2)` for the adjoint ratios at two points. Rewritten as `(eta,zeta)(z_1)`, `(eta,zeta)(z_2)`.
2. **Sentence after `thm:gauges`.** It said that `eta_0` and `eta_1` need not be positive but gave evidence for `eta_0` only. Checked numerically: `eta_1<0` on part of each of the six cycles of the check script (minima −0.06 to −0.98), `eta_0<0` at `(0.2,2)` and `(0.1,1.1)`. The sentence now says this. Section `sec:evidence` should carry these numbers (they are in `s05_gauges.check.out`).
3. **Remark `gau:remarks`(b).** The formula needs `f,g` periodic; added.
4. **Paragraph after `prop:nestedradii`.** The sentence that a radially scaled copy of a cycle "is crossed inwards/outwards" presupposed a planar flow, which is constructed only in Section `sec:attractor`. Replaced by what `gau:radialratio` proves (the ratio increases with theta on one half-plane and decreases on the other), and the infinitesimal statement is now marked as formal and not used.
5. **Paragraph after `prop:lvchain`.** "three explicit pairs, ..., `(eta_1, .)`" replaced by "three explicit functions of the state".
6. **Paragraph after `prop:strongidentity`.** "The two formulas (gau:strongformulas) ... express the same fact" was accurate for the second formula only; corrected.
7. **Closing subsection.** "the linear equations (eq:projective)": the projective equations are not linear; reworded ("in which x enters as a coefficient").
8. **Proof of `thm:gauges`(iii).** The last step (interior of the hull = bounded complementary component) was implicit; one sentence added.
9. **Opening paragraph.** The list of results used now names Section `sec:adjoint` and the compound system `eq:compound`.
10. **Layout.** Two displays overfull by 50 pt and 35 pt (`prop:kappa`(v) and the proof of (iii)) were split; no overfull box remains.

No sign or index error was found in any displayed formula after the renaming `V,S*,S1 -> Sigma`, `Q,V,S4 -> Xi`,
`V3,S2 -> Sigma eta_0`, `beta,kappa(lane) -> frak r, frak s`, `q,P -> ell`.

## 5. Statements weakened or strengthened

- Weakened: none.
- Stronger than the lane source: `prop:strongidentity`(iv) holds for every `mu>0` (lane: `mu<=1`). The proof was re-derived line by line: `d=c-(1-mu)b` satisfies `d'=-mu d-xa`, and the periodicity argument needs only `gamma>mu`, which follows from `eq:etaminbasic`. Confirmed numerically at `mu=3`. The same bounds are proved independently in Lemma `lem:strongstable`(b).

## 6. Formulas re-checked and how

- **By hand**, line by line: every display of the section (Lemma `gau:wronskian`; all six parts of `thm:gauges` and each step of its proof; `gau:halfcycle`; Cramer's rule, the `kappa`-coefficient of `W`, the Riccati equation, the three rows of the compound form, the branch integral, the two-point formula `gau:deltakappa`; the two dynamic bounds; the polar system and the substitution `dt=-dtheta/sigma`; `g(s_2)-g(s_1)`; the chord-sign computation and both pairs of supporting inequalities in `prop:monotoneintegrals`; the AM–GM and Jensen steps and the convexity ranges of `f_mu`, `chi`, `G_mu`; the Lotka–Volterra and Cole–Hopf computations; `nu'=-C^T nu`, the two instances, the scalar identities, the defects `g_i`, `g_2-g_1=R_1/eta~`, `g_3-g_1=R_2/zeta~`, the variation-of-constants inequality; the Liénard form).
- **Exactly (sympy)**: `ROOT/work/framework_checks/s05_gauges.check.py`, Part A, 94 identities, 0 failures. Blocks A1–A6 are the writer's; block A7 (22 identities) was added at the audit: the general Wronskian lemma with arbitrary `f,g,eta,zeta`; the mixed gauge; `(p/U)'`, `(p'/p)'`; both identities of `gau:dynamicbounds`; the intermediate steps of the proof of `gau:Xistate`; `(Sigma eta_0)'`; the `kappa`-coefficient of `W`; `f_mu''`, `G_mu'`, `chi''`; the factor chain `(D+A)(D+1)(D+mu)ell`; `a''`, `a'''`; the polar identities for a shifted centre; `pK/Sigma=K sin/(r sigma)`.
- **Floating point** (same script, Part B; whlib cycles at `(mu,rho)` = (0.3,1.5), (1,1.2), (1,2), (0.2,2), (3,1.5), (0.1,1.1)): 0 failures. Writer's checks: gauge positivity and gauge equations; four exponent formulas against `<eta>-mu` and the monodromy; means of `eta_0`, `eta_1`; half-cycle bounds; `kappa`-Riccati; branch integral; polar form; monotone integrals in 200 random directions; Hölder sandwich; Lotka–Volterra means; strong-stable formulas; comparison bound on 300 random trial pairs; orientation of the Gauss curve. Added at the audit (block (j)): half-cycle bounds with `Sigma`; minima of `eta_1`; `zeta<y`, `eta<ell` and the three upper bounds for `kappa`; `gau:supportline` and the supporting planes on 240 x 24001 ordered pairs per cycle; the orderings at 400 equal-`U` and 400 equal-`p` levels; `gamma<gamma_E`; the compound monodromy (entrywise nonnegative, spectral radius `m_t` to 9 digits).
- Proposition `prop:nestedradii` cannot be tested on data (no parameter with two cycles is known); its algebra is in blocks A4 and A7.

Output of the last run: `ROOT/work/framework_checks/s05_gauges.check.out` ("exact identities failed: 0", "numerical checks failed: 0").

## 7. Check script

`s05_gauges.check.py` was moved from `BUNDLE/paper/sections/` to `ROOT/work/framework_checks/` and extended there (block A7 and block (j)). Run:
`/Users/admin/.venvs/analysis-py313/bin/python -W ignore s05_gauges.check.py > s05_gauges.check.out` (single process, a few minutes).

## 8. Compile test and bibliography

`pdflatex -interaction=nonstopmode -halt-on-error test_s05_gauges.tex` in `BUNDLE/paper`: no error, 15 pages, no overfull box; the undefined references are cross-section labels and the citation key. The auxiliary files (`.aux`, `.log`, `.out`, `.pdf`) were removed from `BUNDLE/paper`; the wrapper is kept. The only citation, `GG21`, is defined in `BUNDLE/paper/bib/s01_intro.tex`; no bibliography fragment was added.
