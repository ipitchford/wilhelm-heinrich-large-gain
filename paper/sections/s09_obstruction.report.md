# Report: s09_obstruction.tex

Written by the orchestrator from the v2.0 manuscript, Section 7 (Lemma 7.1, Theorem 7.2, Corollary 7.3), with the referee's request 4.3 applied.

## Statements

| Label | Content | Source |
|---|---|---|
| `lem:collar` | Exactly one cycle, hyperbolic and attracting, for k in a right-hand neighbourhood of k_H at fixed mu. Proof in six steps with two nested neighbourhoods. | v2.0 Lemma 7.1, rewritten |
| `thm:neutral` | Neutral set relatively closed; count locally constant off it; a second cycle at k_1 needs a neutral cycle at some k_* in (k_0, k_1]. | v2.0 Theorem 7.2, proof made explicit |
| `obs:cor:uniqueset` | The set of parameters at which every cycle has lambda_c < 0 is open, its relative boundary lies in the neutral set, it contains the collars and the large-gain region. | new; follows from `thm:compact`, `thm:neutral`, `cor:criterion` |
| `obs:rem:twoways` | Relation between the parameter-wise criterion and the continuation statement. | new remark |
| `cor:neutralconditions` | Equivalent neutrality tests; variance window. | v2.0 Corollary 7.3 |

## Changes relative to v2.0

- Lemma `lem:collar`: nested neighbourhoods V inside U; common box from `lem:statebox` in physical variables; finite-cover time tau; uniform continuity in k on the compact box; invariance of each cycle under the time-tau map.
- Step 2 (critical global stability) now cites the published statement and also gives a derivation from this paper's own results: Jensen on f_mu makes every cycle at 0 < k <= k_H attracting; by `thm:attractor` (A1), (A4), (A5) the attractor is a connected set on which the flow is planar; Poincaré–Bendixson in reversed time, Lyapunov stability of E and orbital stability of cycles leave only E. The v2.0 route through Smith's permanence and Poincaré–Bendixson theorems is no longer needed in this step.
- Expansion of f_mu(mu + q) written out in the proof of the variance window.

## Cross-model review (work/sol/review_s08_s09_response.md)

Verdicts returned: HOLDS for `thm:neutral`, `obs:cor:uniqueset`, `cor:neutralconditions`; GAP for `lem:collar`.

| Point | Disposition |
|---|---|
| Step 1 wording "remains in a neighbourhood for all times" | Repaired: complete orbits, hence periodic orbits, lie on the centre manifold. |
| Step 1 crossing orientation | Repaired: the appendix gives a positive crossing speed; the sentence now says so. |
| Step 2 hypotheses of imported permanence and Poincaré–Bendixson results unchecked | Repaired by replacement: published statement cited for the exact system, plus a self-contained derivation (above). |
| "eta is nonconstant" | Repaired: reference to the remark after `eq:gainprimitive`. |
| Persistence of the transverse unit multiplier in the limit | Repaired: argued through the crest-return derivatives. |
| Implicit-function step | Repaired: transversality of the crest section and invertibility of I - M stated. |
| Strictness in the variance window | Repaired: explicit argument for both inequalities. |
| `\dd` without a variable | No change: the macro expands to dt. |

## Checks

- Expansion f_mu(mu+q) = k_H + (A^2+mu) q + 2A q^2 + q^3 re-derived by hand.
- Wrapper `paper/test_s09_obstruction.tex` compiles without errors.
