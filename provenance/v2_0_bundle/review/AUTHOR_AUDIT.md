# Author-side mathematical audit

This is a recorded self-audit, not an independent referee report or formal proof-assistant verification.

## Critical dependency chain

1. Positive adjoint evolution gives positive periodic ratios without presuming cycle stability. Liouville and Perron dominance identify both transverse multipliers, including a possible additional +1.
2. The normal equation and exact exponential remainder give strict support between bounded entire solutions. Strictness is justified by ODE uniqueness, not by an unjustified pointwise positive remainder at every instant.
3. Support in both directions proves injective projection and the two different equilibrium cone exclusions. The cycle-to-equilibrium inequality alone does not supply both signs.
4. Tangent orthogonality makes all axis crossings transverse and oriented. Jordan winding counts one crest and trough. The second-order tangent equation and the linear-filter derivative then order the six nondegenerate extrema.
5. The periodic physical state box is uniform on compact positive parameter sets. Cone exclusion bounds the angle from the stable eigendirection. A positive quadratic form on the unstable plane keeps cycles away from the equilibrium, uniformly on compact strictly supercritical sets.
6. Explicit period-independent bounds on the projective ratios and tangent identities give a lower projected speed. Unimodality bounds arclength, hence least period. Crest transversality prevents a limit period from becoming a multiple traversal.
7. Properness plus hyperbolic continuation makes cycle count locally constant away from the relatively closed neutral locus. The global Hopf starting collar uses a common **physical** positive box, avoiding the parameter-dependent log-coordinate-domain pitfall.
8. A finite fixed-rate path from that collar to any finite multiple-cycle point must meet a neutral cycle. The neutral point may be the endpoint. Nothing excludes neutrality itself.

## New envelope argument

For a finite union of cycles and the equilibrium, all support slopes are bounded and strictly signed. The supremum of their affine planes is a finite convex Lipschitz function. Strict pairwise support identifies its unique contact point at each selected orbit point. The curvature identity gives a nondegenerate phase maximum and the explicit rank-one Hessian.

Forward epigraph invariance follows by evolving every base phase and using the finite-time support equation; the family of base phases is invariant. The graph defect is nonpositive and is not generally zero. Therefore this construction must not be called a global invariant manifold or used to assert equivalence of all planar and spatial trajectories. The envelope depends on actual selected cycles and is not an algorithm for finding them.

## All-positive-rate audit

No filter-order swap is used in the final theorem. The projective lower bound is eta > max(0, mu-1), not merely eta > 0. This is what makes the cubic gain polynomial increasing and strictly convex for mu > 1 as well. The period inequality and all compactness arguments are valid for arbitrary positive mu. Reciprocal-rate numerical comparisons check the coordinate symmetry without changing the labels in the six-extrema theorem.

## Neutrality and degeneracy

The +1 eigenvalue of the two-dimensional crest return is simple because its other eigenvalue is exp(-(1+mu)T)<1. The full monodromy also has its autonomous tangent multiplier +1 and may have a Jordan block; no second linearly independent periodic variational solution is assumed. A neutral event is necessary for a count change but does not by itself establish a fold or multiple cycles.

The variance window follows from the exact cubic gain identity, mean neutrality, q > -min(mu,1) and an explicit upper projective bound. Its inequalities do not contradict each other. No unweighted mean sign follows from a weighted identity or curvature alone.

## Computational audit

The active exact checker tests residual equality to zero; it does not count literal pass assertions as checks. The initial off-constraint residual mistake is preserved and explained. The full active code was subsequently rerun, including regeneration of the 31 BVPs, before release.

All numerical diagnostics are floating point. Finite grids do not prove strict signs between samples, and failed historical neutral solves do not exclude neutral cycles. The one-crest proof, not the sampled roots, establishes exhaustive phase structure. Tiny determinants are not treated as relatively accurate floating-point enclosures.

## Remaining uncertainties

The analysis still requires independent specialist scrutiny, particularly the strict-support-to-properness chain and its new envelope consequences. No fatal gap was identified in this self-audit, but self-audit is not independent validation. Historical priority remains open because the literature search is bounded and original full texts remain missing.
