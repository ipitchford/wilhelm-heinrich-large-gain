# Final targeted internal check of the inserted slow-graph lemma

Date: 2026-10-03. Disposition: **no material mathematical gap found in the inserted derivation of (F)**. This is an internal adversarial application check by the agent that proposed the construction, not an independent external review or verification of the entire wedge theorem.

## Exact material inspected

- `paper/sections/s11_slowgraph.tex`, all 65 lines, SHA-256 `31899758625b64d36c58afc12c596ed87e2f618551c2474d01e433bb5e257161`.
- Adjacent property statement, invocation, theorem statement, and initial proof steps in `paper/sections/s11_wedge.tex`; file SHA-256 at check time `bf032836352eb7eb63090f5085c3ad439cb67b5fd3a017bf2e2d07379e391469`.
- Eldering, author-hosted May 1, 2018 pre-proof PDF: <https://www.jaapeldering.nl/files/NHIM-noncompact-book.pdf>, Theorem 3.2 and Remark 3.3(iii), printed p.62/PDF page 70; theorem statement reopened during this final check. Definitions 1.8/1.11 and §4.2 were inspected in the preceding check. Author's errata: <https://www.jaapeldering.nl/files/errata.pdf>.

Theorem 3.2 explicitly allows vector fields on a uniformly sized neighborhood of the zero section. It returns a graph over the entire base with the requested finite regularity, under bounded/uniformly continuous derivatives, sufficiently small C^1 perturbation, and an attracting normally hyperbolic initial graph. Remark 3.3(iii) permits a constant family when the initial graph is zero. The lemma uses exactly this setting, not a boundary persistence result.

## Adversarial checks

1. **Coordinate calculation.** With z=Q+g(U), differentiation gives z_s=-z+epsilon(g'-1)P-r epsilon^2 g. Since g'-1=g, the displayed transformed equation is correct. The transformation is epsilon-independent and smooth.

2. **Initial graph and spectral gap.** On X=R^4 the zero section is invariant and complete. Its tangent derivative is the identity and its stable normal derivative is exp(-s). The splitting is constant, its projections have norm one, and there is no unstable bundle. In the source's rate notation one may take rho_M=0 and rho_-= -1, so every finite regularity gap is satisfied. The auxiliary graph-family condition is satisfied with h_sigma=0 identically. Taking regularity max(N,2) also covers N=1.

3. **Uniform tube and mixed derivatives.** Fix d>0. On |z|<d, the cutoff field and unperturbed field are bounded. Each mixed derivative of their difference through any fixed finite order is O(delta_*)+O(delta_*^2), uniformly in (U,P,r,a,z). Compact support of chi controls g and polynomial factors; compact support of b and bounds on its derivatives control a. Differentiation in a produces delta_* b^(m)(a), not a reciprocal power of delta_*. Terms involving z are uniformly bounded on the fixed tube. Bounded derivatives of one higher order give uniform continuity of all required derivatives. Thus the cited theorem's norm hypotheses genuinely hold.

4. **Flow issue.** The theorem does not require the full ambient tube to be invariant. Its graph is chosen with sufficiently small height inside the tube. Even the displayed global extension has at most linear growth in z: its z equation is -z plus uniformly bounded forcing, and its other nonzero components are bounded by C(1+|z|). There is no finite-time escape obstruction to using its flow. All tangent/normal exponential estimates required for the initial graph are exact and uniform.

5. **Parameter regularity and epsilon=0.** The parameter a is part of the base of a single C^j graph, so this is not an unsupported smooth selection among separately existing manifolds. At a=0 every base component vanishes. Tangency of the graph then forces its normal component -h to vanish. Hence h(U,P,r,0)=0 exactly. For fixed positive delta_*, composition with a=epsilon/delta_* preserves C^N. Derivatives may grow as delta_* decreases, but (F) requires no uniform derivative bound as this construction parameter tends to zero.

6. **Boundary and original equation.** Smoothness holds on all base coordinates, thus on a neighborhood of the closed target set and signed epsilon values around zero. There is no loss of derivatives at epsilon=0 or at the spatial boundary. On the target neighborhood chi=1 and b(epsilon/delta_*)=epsilon/delta_*; substituting z=w+g into graph tangency yields the exact stated invariance equation. At epsilon=0 this reduces to w=-g. Local, rather than global, invariance for the original system is the appropriate conclusion.

## Minor clarification and remaining limits

- Optional proof-explicitness improvement: before the theorem application, choose an allowed graph-height bound eta<d/2, then choose delta_* small enough for that eta. This is already available from the theorem's quantifiers; its omission is not a material gap.
- The adjacent wedge proof still contains the phrase “the assumption (F)” after (F) has become a proved property. This is a stale editorial label, not a logical issue. No manuscript edits were made during this check.
- The revised local theorem statement retains a fixed compact r interval, a fixed energy ball, and a nonexplicit threshold. Discharging (F) does not prove global uniqueness, statement (N), uniformity as r tends to zero, or the relaxation-limit claims.
- This check does not independently certify Eldering's full proof or the remainder of the manuscript. It confirms a direct, hypothesis-matched use of the stated primary theorem. It is sufficient to support the local (F) lemma as written, while publication-level review and all other evidence gates remain separate.
