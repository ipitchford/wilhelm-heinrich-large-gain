# Changes from the attached v1.2 working paper

## Mathematical changes

- Remove the restriction 0 < mu <= 1. Use fast time equal to the original second filter rate and allow every positive physical rate ratio. All statements retain the original intermediate filter Z.
- Prove explicit positive upper and lower bounds on both projective ratios independently of cycle period. This replaces the earlier qualitative positive-evolution bound in the compactness proof.
- Extend the one-crest result to one nondegenerate maximum and minimum in each concentration, with the strict order X-min, Z-min, Y-min, X-max, Z-max, Y-max.
- Quantify the ordering of crest heights and their supporting slopes for nested cycles.
- Construct a common convex epigraph from a finite selected family, prove its forward invariance and strict entry, and compute the local rank-one Hessian and exact graph defect. Explicitly reject the stronger, unsupported claim of a global invariant manifold.
- State relative closedness of the neutral parameter locus, local constancy of the count off that locus, a trace test and a necessary neutral-variance window.
- Make the fixed-rate obstruction endpoint-inclusive: the neutral event can occur at the proposed counterexample parameter itself.
- In the Hopf-collar proof use a common box in physical positive variables; a box assembled in parameter-dependent log coordinates need not lie in the critical chemical domain.
- Retain the analytically false square-root radial certificate as a documented negative result, not a cycle counterexample.

## Evidence and availability

Six v1.2 files were available: its manuscript PDF and LaTeX and four receipts. No source code or meshes for its described continuation, multistart, supporting-geometry or neutral-search experiments were available in those attachments. Accordingly no reproduction of those experiments is claimed.

The active code was written and executed afresh. It provides 31 complete BVP polynomial meshes and 31 audit profiles; four integration-method repetitions; 12 coordinate-symmetry comparisons; four envelope checks; and four deliberate-error controls. The two preserved historical neutral-search receipts are treated as failed searches only.

The earlier v0.3 review is preserved as provenance. Its favourable assessment is not attributed to this paper. The full original 1995 and 1996 source papers remain an access gap; the current mathematics does not silently import unseen results from them.

## Current recorded failures

The first exact run failed one of 50 tests because the off-constraint residual used the wrong multiple of a constraint that vanishes on cycles. The original failed receipt and script are retained; the corrected exact suite passes all 50 equalities.

The intentionally false exact-sign check must fail and exit with status 1. Its receipt is retained with passed=false. The harness reports success for detecting this failure. The square-root contraction test fails on nine sampled return matrices, consistent with its analytic obstruction. These failures are neither deleted nor relabelled as successful certificates.
