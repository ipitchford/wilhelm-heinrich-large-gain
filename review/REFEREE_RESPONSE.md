# Response to the referee report on v2.0

The report on v2.0 (reproduced verbatim in `review/v2_0_referee_report.md`) recommended minor revisions and found no fatal mathematical defect. This note records how each point was handled in v3.0. Section and label names refer to `paper/manuscript.pdf`.

The report assesses v2.0. It is not an assessment of the new material in v3.0, and none of its conclusions is transferred to that material.

## Major comments

| Report item | Request | Action in v3.0 |
|---|---|---|
| 4.1 Novelty | Make the novelty claim theorem-specific; separate established results, direct proofs of general phenomena, and results claimed as new; keep the access limitation for unobtained full texts. | Section 1.4 and Table 1 classify every main result as imported, known, direct proof, or new, with the nearest prior work named for each. Analogues found in a wider search are cited (planar projection of limit sets of competitive systems; two-dimensional attractors of scalar delay equations; convex carrying simplices; the proof of Jones' conjecture). The statement that the 1995 and 1996 originals were not obtained is retained, and `SOURCES.md` records the access level of every source. |
| 4.2 Compactness | Give the uniform compactness argument prominence; extract the equilibrium-separation and projected-speed arguments into named lemmas; state that the constants depend only on the compact parameter set. | Section 8 now consists of Lemma 8.1 (uniform physical state box), Lemma 8.2 (uniform separation from the equilibrium), Lemma 8.3 (uniform projected speed) and Theorem 8.4, preceded by a dependency summary. The three inequalities named in the report are displayed with their constants, and each lemma states that its constants depend only on the parameter set. The remainder estimate in the separation argument is written out. |
| 4.3 Hopf collar | Make the local-to-global passage explicit with two nested neighbourhoods. | Lemma 9.1 is rewritten in six steps: the local Hopf neighbourhood, global stability at the critical value, the common box in physical variables, a compact neighbourhood inside the Hopf neighbourhood reached by the critical flow in a fixed time (finite cover), continuity in the parameter on the compact box, and invariance of each cycle under the fixed-time map. |

## Minor comments

| Request | Action |
|---|---|
| Distinguish return matrices in (U,y) and (U,z) | Proposition 4.4 states that M is in section coordinates (U,y) and gives the similarity to the matrix in (U,z); Section 12 uses the subscripted matrix. |
| Captions, numbers and time units for tables | All tables are numbered and captioned; periods are labelled by their time variable; the amplitude-fixed sampling design is stated in the caption of the cycle table. |
| Notation table | Table 2 in Section 1.5. |
| Periodic orbit, numerical approximation, validated enclosure | The three terms are kept distinct. Nothing in the paper is an interval enclosure, and Section 13 says so for each class of computation. |
| Schematic of the projected cycle | Figure 1 shows the projected attractor and cycle and the six ordered extrema at one parameter. It is an illustration and is labelled as such. |

## Suggestions for strengthening

| Suggestion | Outcome |
|---|---|
| A transferable theorem for a broader class of convex negative-feedback systems | Not attempted as a general theorem. Two results are stated for the model but their proofs use only the convex competitive structure: the attractor theorem (Section 7) and the kurtosis bound for strictly convex potentials (Section 11), which is stated for any strictly convex potential. |
| Rigorous exclusion of neutral cycles on a nontrivial parameter region | Done for large gain: Theorem 10.1 shows that every cycle is hyperbolic and attracting when k >= 1000 (mu <= 1), with the reciprocal statement for mu >= 1, and uniqueness follows. |
| Validating one continued branch would not exclude other branches | Agreed. The parameter-wise criterion (Corollary 7.13) avoids continuation: a bound on the weak exponent of every cycle at a parameter gives uniqueness at that parameter. |

## What the report said would change its decision

The report named two events: a failure of the strict-support-to-equilibrium-separation argument, or an older general theorem already supplying the compactness and envelope results. The separation argument was written out in full for Lemma 8.2 and survived. The wider literature search found analogues but no theorem covering these results; the search was bounded and is described in `SOURCES.md`.

## Scope

The report noted that v2.0 should not be described as having effectively solved the uniqueness problem. That remains true of v3.0: uniqueness is proved for large gain and near each Hopf point, and is open in between.
