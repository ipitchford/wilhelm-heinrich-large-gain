# Changes from v2.0, and what is missing

## What v3.0 adds

| Area | v2.0 | v3.0 |
|---|---|---|
| Uniqueness | Not proved anywhere away from the Hopf collar. | Proved for large gain with an explicit threshold: k >= 1000 for mu <= 1, and k >= 1000 mu^3 for mu >= 1 (Section 10). Every cycle there is hyperbolic and attracting, with explicit bounds on its period, up-crossing slope, crest duration and weak exponent. |
| Criterion for uniqueness | Continuation from the Hopf collar: non-uniqueness needs a neutral cycle somewhere on the path. | In addition, a parameter-wise criterion needing no continuation: if every cycle at a parameter has lambda_c < 0, the cycle is unique and attracts every positive orbit off the stable manifold of the equilibrium (Corollary 7.13). |
| Supporting planes | Base a cycle or the equilibrium. | Any two bounded entire orbits (Theorem 7.6). |
| Invariant geometry | A convex epigraph built from the cycles, explicitly not an invariant manifold. | The global attractor itself is the graph of a convex C^1 function over a planar Jordan domain; the flow on it is planar with divergence eta - mu (Theorem 7.1). The v2.0 envelope construction is superseded and is not repeated. |
| Shape of cycles | Injective projection, one crest, six ordered extrema. | Also strictly convex and star-shaped in projection, with exact formulas for the weak exponent, a scalar Riccati reduction and a polar form (Section 5). |
| Linear structure | Positive adjoint; exact exponent. | Also the strong-stable vector, the weak adjoint with a fixed order of twelve events, strong-stable supporting planes, a rank-one identity, and a monotone crest return map on an explicit region (Section 6). |
| Small rate ratio | Not treated. | Exact identities in frequency-one time; a kurtosis bound for convex potentials; a local uniqueness theorem in the near-Hamiltonian wedge, conditional on a property taken from Fenichel's theorem (Section 11). |
| Negative results | One failed radial certificate. | Also: constraints on Dulac functions, no monomial supersolution, first-order rigidity of certificates at Hopf, and a table of candidate structural statements with their status (Section 12). |
| Referee's requests | — | Compactness split into named lemmas with explicit constants; Hopf collar with nested neighbourhoods; positioning table; notation table; coordinate labels of return matrices; a schematic figure. See `REFEREE_RESPONSE.md`. |

The v2.0 statements on the adjoint, supporting planes, geometry, compactness, the neutral-cycle obstruction and the failed radial certificate are retained with unchanged mathematics.

## What is not established

- Uniqueness at every oscillatory parameter. The region kH < k < 1000 max(1, mu^3) outside the Hopf collars is open.
- Any statement excluding a neutral cycle in that region, or exhibiting one.
- Global uniqueness in the small-rate-ratio wedge: the wedge theorem is local in phase space, its threshold in mu is not explicit, and the exclusion of large non-attracting cycles is reduced to unproved estimates.
- An explicit Hopf collar.
- Priority of the results marked new.

## Evidence and availability

- All computations are ordinary floating point, exact symbolic algebra, or high-precision evaluation of constants. Nothing is an interval enclosure.
- The floating-point cartography covers a finite set of parameters. It cannot detect two cycles closer than the scan spacing, and it says nothing about parameters outside the ranges listed in Section 13.
- The v2.0 historical gaps stand: the code and meshes of the v1.2 experiments were never available and are not reproduced; the two v1.2 neutral-search receipts remain failed searches.
- Full texts not accessed are listed in `SOURCES.md`, Section 6. One of them (a 1981 uniqueness theorem for a third-order neuron model) could bear on priority or method and has not been seen.

## Review status

- Each new theorem was attacked by an independent referee agent, which re-derived the proofs and wrote its own numerical tests; the records are in `provenance/research_lanes/*/VERIFY.md`.
- The assembled sections were then refereed adversarially by a second model; the prompts, replies and dispositions are in `provenance/research_lanes/sol/` and in the section reports under `paper/sections/*.report.md`.
- Gaps found at either stage and their repairs are listed in Section 13 of the paper. No theorem was withdrawn. Several statements were corrected or narrowed before inclusion, and those corrections are recorded in the same places.
- None of this is review by specialists in the field.

## Recorded failures kept in the bundle

- The first exact-identity run of v2.0 (49 of 50) and its explanation, as in v2.0.
- Deliberate-error controls, which must fail, for both the v2.0 and the v3.0 suites.
- Initial failures of v3.0 check scripts that were traced to the scripts and not to the mathematics are kept as separate receipts with an `_initial_failure` suffix.
- Candidate structural statements refuted by computation are reported as refuted, with parameters, in Section 12.
