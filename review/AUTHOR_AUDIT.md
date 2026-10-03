# Author-side audit of v3.0

This is a recorded self-audit with adversarial internal review. It is not an independent referee report and not a formal verification.

## Dependency chain of the headline result (uniqueness at large gain)

1. Positive adjoint ratios along every cycle (Section 2) and the explicit linear gauge Lambda_E > 0 (Corollary 3.2). Neither presumes stability.
2. One crest and one trough (Theorem 4.1), the period bound T > 2 pi / sqrt(mu) (Proposition 4.4) and the state box (Lemma 8.1).
3. The positive lobe {U > 0} of every cycle lasts less than pi / omega, by a Sturm comparison that uses only Lambda_E > 0 (Lemma 10.5). With the period bound this forces a long negative lobe.
4. A four-duration pigeonhole on the negative lobe gives a first lower bound on the slope at the up-crossing; a spike-mass inequality and a quasi-convexity bootstrap raise it to the relaxation scale: T >= (11/60) lambda/mu and slope >= 0.819 lambda.
5. The crest interval {x > c_0} is short, with an explicit bound.
6. A linear estimate for the compound system (Lemma 10.16), with one Lyapunov function on the trough and another on the crest, gives <eta> <= 0.83 mu on every cycle.
7. The parameter-wise criterion (Corollary 7.13) turns "every cycle is attracting" into uniqueness and global attraction. It rests on the attractor theorem, not on continuation. A second proof through Smith's convergence theorem is given.
8. The reciprocal-rate correspondence transports cycles and multipliers, not positive domains, to mu >= 1.

Steps 3 to 6 are statements about every cycle. Each inequality was tested on computed cycles and on non-periodic solution segments satisfying the same hypotheses, to guard against use of a property of the particular cycle that numerics finds. Every numerical constant is recomputed at 30 or more digits by a script in the bundle.

## Dependency chain of the attractor theorem

1. Dissipativity with an explicit absorbing box, proved directly (no permanence theorem imported).
2. Along every bounded entire orbit the positive adjoint direction exists and is unique, by contraction of Hilbert's projective metric; its ratios are bounded and continuous on the attractor.
3. The exact comparison identity gives strict supporting planes between any two bounded entire orbits. Hence the projection to the (U,p) plane is injective, and the attractor is the graph of the restriction of a convex function; continuity of the normals gives C^1.
4. Above the Hopf value the two cone exclusions make the equilibrium a source inside the attractor. Planar Poincaré–Bendixson theory on the reduced flow, the ordering of crests and a winding-number argument identify the projection with the closed Jordan domain of the outermost cycle.
5. Multipliers of the planar flow equal the weak multipliers; innermost and outermost cycles cannot be repelling; the parameter-wise criterion follows.

The imported classical results are listed in Section 7.1 with the hypotheses checked at each use. One of them (monotone intersections with a transversal) was first stated without the definition of a transversal that makes it true; the cross-model review caught this and the definition is now explicit.

## Points where the paper was corrected during internal review

- A first-order Melnikov premise for the small-rate-ratio limit was false: that function has no zero. The Hamiltonian regime is a thinner wedge selected at third order; Section 11 is written on the corrected basis.
- Return-map concavity and sublinearity on the crest half-axis, which would give uniqueness, fail in floating point at extreme rate ratios. They are reported as refuted. What survives every test is monotonicity of the one-turn mean of eta - mu; it is unproved.
- A proposition on Dulac functions needed a hyperbolicity hypothesis in one step; it was restated and re-proved.
- The crest return map was first shown monotone without proving that the return exists on the whole region; the existence and regularity proof was added.
- A scaling identity was stated without the range 0 <= theta <= 1 on which its sign holds.
- The polar formula for the weak exponent was mis-stated in one working note (not in the paper); see `provenance/research_lanes/sol/README.md`.

## What the computations do and do not show

Exact symbolic checks test algebraic identities only. Floating-point checks test inequalities at finitely many parameters and phases. High-precision checks confirm the numerical constants and scalar side conditions of the large-gain proof. None of these formalises the non-algebraic steps: Perron–Frobenius and Hilbert-metric contraction, planar topology, compactness, Sturm comparison, or the bootstrap logic.

## Remaining uncertainties

- The large-gain proof is long. It was re-derived lemma by lemma by an independent referee agent and checked by a second model, with no error found beyond wording; that is strong internal evidence and not specialist review.
- The attractor theorem imports classical planar and invariant-manifold results in their textbook form; the textbook was not re-read.
- The wedge theorem takes a property from Fenichel's theorem whose original statement could not be accessed.
- Priority is open; see `SOURCES.md`.
