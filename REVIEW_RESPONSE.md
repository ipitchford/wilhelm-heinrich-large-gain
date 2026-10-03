# Response to the supplied v3 review

3 October 2026. Revised candidate 3.1.0-candidate. This is an internal revision record, not an external acceptance decision.

## Dispositions

1. **Attractor scope — corrected.** The abstract and opening scope now distinguish the positive-interior attractor from the extinction equilibrium and connections in a nonnegative-phase-space attractor. The disk domain is asserted only above Hopf; at and below Hopf the positive attractor is the equilibrium. “Convex graph” means graph of a convex function, not a convex three-dimensional subset.
2. **R derivative — corrected and regression-tested.** Section 10 now gives R'(k)=(2/3)k^(-2/3)(3k^(1/3)-2). The positivity argument survives. The revision controls reject the old derivative.
3. **W interpretation — corrected.** Section 7 describes the Hessian quadratic form, rather than the full material second derivative. The missing acceleration/gradient term is explicitly displayed and covered by a counterexample control.
4. **The 0.83 contraction statement — strengthened checking.** The author-side comparison now tests 0.83, not 0.86. A supplied exact-rational interval implementation was adapted with explicit root/domain guards and portable paths. It encloses 55 endpoint comparisons plus a repeated explicit theorem-bound comparison. The total is 56, not 56 independent proof obligations. At k=1000, q lies near 0.07717738297297608456 and the margin below 0.83 exceeds 0.0028226170270239. This does not certify sampled monotonicity grids, orbit integrations or all analytic reductions.
5. **Assumption (F) — discharged locally.** The revised Section 11 includes a full cutoff construction on a complete Euclidean base and a hypothesis-by-hypothesis application of Eldering Theorem 3.2. Joint regularity follows by adjoining a frozen scaled parameter. The formerly conditional local wedge theorem is now a written-proof result. Source mapping and a targeted internal check are retained under review/. The checker is the agent that proposed the construction, so this is not independent external review. Neither global wedge uniqueness, uniformity as r tends to zero, nor an explicit threshold follows.
6. **Carrying-simplex antecedents — compared at hypothesis level.** The full-text portions of Baigent (2019) and Mierczynski (2020) relevant to their hypotheses were read. Their face-preserving competitive-map settings do not directly include the Wilhelm–Heinrich flow in physical coordinates. The revised positioning explains the strict entire-trajectory supports, logarithmic graph and divergence formula. It does not assert that no other transformation or prior theorem can cover the result.
7. **Priority — qualified.** Smith and Boros–Hofbauer were checked against their primary preprints. Li (1981) and other full-text gaps remain listed. “Not found” is bounded-search evidence, not a novelty certificate.
8. **Reproduction versus environment — separated.** The supplied review audit reports scientific checks passing but an aggregate failed receipt due to dependency-version mismatch. Its raw failed receipt is preserved in the private input audit and not rewritten. The publication replay uses the pinned reference environment. The unresolved historical small-mu comparison remains failed; another integrator's success does not change that receipt.
9. **Reading route — added.** A short main-result route precedes the literature discussion. The small-rate-ratio extension, failed-certificate catalogue, computation section and future programme are identified as supplementary, without discarding the complete research record.
10. **Reusable mathematics — highlighted.** The convex-potential kurtosis theorem is identified as a reusable result. A general-feedback classification or an all-parameter oscillator theorem is not claimed; either would require new research, not an editorial strengthening.

## Review and input identities

- Original research ZIP SHA-256: dbf75f19278da9310c3e08f3c72e57335fc5bd3d888162554409e6cdc306f0cc.
- Supplied audit ZIP SHA-256: 2098cee1f5e7c8fb82ed349344d373b6a88e6c438dd28c697306c1d1499c6e99.
- The original v3 manifest is preserved in provenance/v3_0_manifest.sha256. Historical v2 and research-lane materials retain their original status.
- The supplied review's prospective upper-three-star assessment is not an official REF judgement and is not a public assurance badge.

## Assurance ceiling

Unrefereed mathematical candidate. Written proofs and source-mapped theorem applications have internal checking, selected algebraic/interval calculations have executable verification, and finite floating-point experiments are diagnostics. No specialist review, unaffiliated reproduction, formal verification, empirical validation or exhaustive priority assessment is established.
