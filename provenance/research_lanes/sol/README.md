# Cross-model review record (3 October 2026)

A second model (`gpt-5.6-sol`, OpenAI Responses API) was used as an adversarial reviewer. Its replies are hypotheses; each was checked against the text before any change was made. Files:

| File | Content |
|---|---|
| `attack_prompt.txt`, `attack_response.md` | Request to refute the claim that the remaining sign statement (weak exponent negative on every cycle) is out of reach of the listed exact structure. Verdict returned: PARTIAL. It proved only what the paper already has and named the same missing monotonicity inequality as the research lanes. |
| `review_s05_s06_prompt.txt`, `review_s05_s06_response.md` | Adversarial refereeing of the gauge and weak-adjoint sections as drafted. Every statement returned HOLDS except Proposition `prop:crestmonotone` (two repairable gaps, Steps 3 and 4), plus six wording points. All were repaired in the section files. |
| `review_s07_prompt.txt`, `review_s07_response.md` | Adversarial refereeing of the attractor section as drafted: 14 points (one imported lemma stated without the definition of a transversal that makes it true; steps to write out in the proof of (A1); strictness claims in the crest-axis proposition; citations). Dispositions are recorded in `paper/sections/s07_attractor.report.md`. |
| `review_s10_prompt.txt`, `review_s10_response.md` | Adversarial refereeing of the large-gain section as ported. Verdict HOLDS for the theorem, both corollaries and all thirteen lemmas, with the key computations and numerical margins re-derived; eight minor points (a derivative to display, one strict inequality, one citation, the notation of one constant, two steps in the second uniqueness proof, the normalised variational equation). All eight were applied to the section file. |
| `review_s08_s09_prompt.txt`, `review_s08_s09_response.md` | Compactness and obstruction sections. HOLDS for the state box, equilibrium separation, speed and compactness statements, the neutral-locus theorem and both corollaries; GAP for the Hopf-collar lemma (crossing orientation, wording of the centre-manifold step, unchecked hypotheses of an imported step). Repaired; see `paper/sections/s09_obstruction.report.md`. |
| `review_s11_prompt.txt`, `review_s11_response.md` | Wedge section. HOLDS for the identities, the Melnikov proposition, the kurtosis theorem, the root corollary and the normal-form lemma; GAP for the conditional wedge theorem (details to write out; status of the Fenichel import). Dispositions in `paper/sections/s11_wedge.report.md`. |
| `review_s12_prompt.txt`, `review_s12_response.md` | Negative-results section. HOLDS for six statements; GAP for the rigidity lemma (two steps to write out) and for several prose overclaims. Dispositions in `paper/sections/s12_negative.report.md`. |
| `call_sol.sh` | The call wrapper. `*.json` are raw API responses; `*.json.req` the requests. |

## Erratum in `attack_prompt.txt`

Item (6) of that prompt wrote the polar form of the exponent as an integral of `(K/r) sin(theta)` against `|d theta|`. The correct weight has an extra factor `1/sigma^2`, with `sigma = Sigma/r^2 = |theta'|`:

    lambda_c T = - integral over one turn of  K sin(theta) / (r sigma^2)  d theta.

The same item also gave the correct equivalent forms (`-<pK/Sigma>` and the integral of `kappa_up - kappa_down` against `dU`), and the verdict did not rely on the mis-stated one. The paper's Proposition `prop:polar` has the correct formula, which the second review re-derived.
