# Targeted check of assumption (F)

Date: 2026-10-03. Scope: bounded internal mathematical/source check, not an independent referee report or a review of the whole wedge proof. No manuscript edits.

## Finding

A precise freely accessible primary persistence theorem supports a short derivation of (F), using an explicit cutoff construction below. This avoids the manuscript's unchecked overflowing-boundary argument. The result is a proposed, source-mapped proof for review, not an authorization to silently remove all conditional labels from v3.0. Retain the existing conditional status in this release unless the derivation is incorporated and independently checked. The original Fenichel/Jones theorem numbering remains unverified.

## Sources actually checked

1. Jaap Eldering, *Normally hyperbolic invariant manifolds — the noncompact case*, author-hosted pre-proof PDF, title-page date May 1, 2018: <https://www.jaapeldering.nl/files/NHIM-noncompact-book.pdf>. Theorem 3.2 and Remark 3.3(iii), printed p.62 (PDF page 70); Definitions 1.8/1.11, printed pp.17–18; §4.2, printed p.114 (PDF page 122). The author's publication page identifies the PDF as a pre-proof version of the 2013 monograph: <https://www.jaapeldering.nl/research.php>. Published-book DOI: <https://doi.org/10.2991/978-94-6239-003-4>.

   Theorem 3.2 gives a globally defined smooth invariant graph for a sufficiently small perturbation of an attracting normally hyperbolic zero section over a complete bounded-geometry base. Remark 3.3(iii) permits a constant auxiliary family when the initial graph is zero. Section 4.2 explains the frozen, scaled parameter construction. These portions, not the entire monograph proof, were inspected.

2. Author's errata, dated July 1, 2021: <https://www.jaapeldering.nl/files/errata.pdf>, pp.1–2. Corrections include a theorem notation fix, a sharper norm in Corollary 3.5, clarification of boundary hypotheses, and a replacement norm in Appendix C. The errata do not withdraw Theorem 3.2. The boundary extension is not used below.

3. Jones publisher record: <https://doi.org/10.1007/BFb0095239>. Bibliographic metadata accessible; full chapter not verified. The CiteSeerX full-text route returned 403. Fenichel publisher route <https://www.sciencedirect.com/science/article/pii/0022039679901529> failed to fetch. These access outcomes are not evidence against either theorem.

The package's `SOURCES.md`, sections 5–6, accurately records the earlier Fenichel/Jones access gap. It does not yet contain Eldering or the present argument. A targeted Scholar query did not supply a usable theorem text; the author-hosted text is the substantive evidence.

## Explicit hypothesis mapping and derivation

Let K be the prescribed compact subset of the variables x=(U,P,r). Choose a bounded open O containing K, with compact closure, and a smooth compactly supported cutoff chi(x) equal to one on a neighborhood of that closure. Set g(U)=exp(U)-1 and z=Q+g(U). In fast time s=tau/epsilon, direct differentiation gives

    U_s = epsilon P,
    P_s = epsilon (z-g),
    r_s = 0,
    z_s = -z + epsilon gP - r epsilon^2 g.

The last identity uses g'=g+1. This coordinate change is smooth and invertible without an epsilon-dependent singularity.

Choose a smooth compactly supported scalar function b(a) that equals a for |a|<=1. For a small fixed positive lambda, put e(a)=lambda b(a), and adjoin the frozen coordinate a_s=0. On the trivial bundle X times Y, with X=R^4 having coordinates (U,P,r,a) and Y=R having coordinate z, define

    U_s = chi e(a) P,
    P_s = chi e(a) (z-g),
    r_s = a_s = 0,
    z_s = -z + chi [e(a) gP - r e(a)^2 g].

The unperturbed field is v0=(0,0,0,0,-z). Its zero section M=X times {0} consists of equilibria. The tangent flow is the identity; the stable normal flow is exp(-s). Thus it is k-normally attracting for every finite k, with constant splitting and bounded projections. X is smooth, connected, complete, Euclidean, and of bounded geometry. The initial graph h=0 meets the theorem's small-graph condition identically.

Work on any fixed uniformly sized tube |z|<d around M, as allowed explicitly in Theorem 3.2. Both fields and all required finite derivatives are bounded and uniformly continuous there. Compact support of chi controls the exponential g and factors P,r; compact support of b controls a. Every derivative of e is O(lambda). Consequently the perturbation tends to zero in C^k, in particular C^1, as lambda tends to zero. There is no need for global bounds on exp(U) before cutoff.

Apply Theorem 3.2 with regularity k=max{N,2} and Hölder index zero, taking its auxiliary family constant as permitted by Remark 3.3(iii). For sufficiently small lambda it gives an invariant graph z=h(U,P,r,a) of class C^k over all of X. Joint parameter regularity is obtained here by putting a in X itself, rather than assuming that separately chosen invariant manifolds depend smoothly on epsilon.

For 0<=epsilon<=epsilon0=lambda/2 define

    w(U,P,r,epsilon) = h(U,P,r,epsilon/lambda)-g(U).

This is jointly C^N on a neighborhood of closure(O) times [0,epsilon0]. At a=0 the modified vector field has zero base velocity and normal velocity -z. Tangency of the invariant graph therefore forces h(x,0)=0, giving w(x,0)=-g(U) exactly. For x in O and a=epsilon/lambda, chi=1 and e(a)=epsilon; hence the modified system is exactly the original transformed fast system there. Graph tangency yields

    epsilon (w_U P + w_P w) = -w-epsilon P-(1+r epsilon^2)g(U),

the manuscript's invariance equation. It also holds at epsilon=0 by the identity just established. This is precisely (F).

The construction bypasses boundary hypotheses by extending the cutoff system over a boundaryless complete base, then restricting a smooth graph to O. It does not assert that a chosen slow manifold is canonical or globally invariant for the unmodified system outside O.

## Dependency and limits

- `paper/sections/s11_wedge.tex:363` states (F). At line 400 the proof requests only N=6 and K=Omega_(xi_b+1) times [r1,r2]. Step 0 needs the joint epsilon Taylor expansion with C^2 remainder. The subsequent return-map differentiability arguments use the same joint regularity.
- This construction permits epsilon0 and the selected graph to depend on N and K. It does not give one C-infinity graph on a fixed epsilon interval, analyticity, explicit thresholds, or bounds uniform over growing energy balls.
- It does not establish statement (N), global cycle uniqueness, uniformity as r tends to zero, the fixed-rho relaxation limit, or correctness of the rest of the wedge proof.
- No primary-source PDFs were added to the publication bundle, and no local checksum/archive of the new source was produced in this bounded check. A publication revision should archive/version the accessed source outside redistributed materials and review this proof before promoting the theorem's status.
- Recommendation: the (F) gap looks locally repairable by established theory, not a new open research problem. Keep v3.0 conditional as written until a deliberate, reviewed revision adds this lemma and its source mapping. Do not relabel the existing unchecked Fenichel/Jones sketch as verified.
