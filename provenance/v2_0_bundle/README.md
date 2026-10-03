# Wilhelm–Heinrich: finite-parameter obstruction and cycle geometry

**Version 2.0 — 2 October 2026. Research candidate for consideration by Evidence Press.**

Start with `paper/manuscript.pdf`. The paper is self-contained for its stated direct finite-parameter results, with standard theory and the published model-specific inputs explicitly identified. It does not require the earlier singular-limit drafts.

## Exact claim boundary

For arbitrary positive filter rates, every positive nonconstant periodic orbit has a simple clockwise logarithmic phase-plane projection and precisely one crest and trough in each concentration. Their six extrema have a fixed order. Distinct cycles are nested, their crests are strictly ordered, and their supporting planes form a convex forward-invariant epigraph. That envelope is **not** asserted to be an invariant manifold.

All cycles form a compact family, after fixing crest phase, on every compact strictly oscillatory parameter set. The cycle count can change along a finite parameter path only if a nonconstant cycle acquires a second Floquet multiplier +1. Starting from the global Hopf collar, any failure of uniqueness therefore requires such a neutral cycle at a finite parameter, possibly at the target endpoint.

**Not proved:** that no neutral cycle exists; uniqueness at every oscillatory parameter; a numerical exclusion over a parameter rectangle; or a validated multicycle counterexample. A neutral cycle alone would not refute uniqueness. The new analytic arguments have not received independent mathematical referee review, and priority is not established.

## Contents

- `paper/`: final PDF, editable LaTeX and evidence-generated tables.
- `code/`: fresh exact checks, periodic boundary-value solves, projective/linear-adjoint and return-map audits, envelope checks, negative controls, table generator, run driver and consistency checks.
- `evidence/cycles/`: 31 complete saved collocation polynomials and meshes, not just figures.
- `evidence/profiles/`: the 31 projective, state, support and curvature profiles used by the current audit.
- `evidence/`: all current JSON receipts, including the deliberately failed claim and the original checker failure.
- `provenance/`: eight byte-identical input snapshots (six v1.2 attachments and two earlier-review files), plus the first failed exact-check script.
- `review/`: author-side proof audit, dependencies, missingness and changes from v1.2.
- `SOURCES.md`, `CLAIMS.json`, `STATUS.json`, `ENVIRONMENT.json`, `MANIFEST.sha256`: evidence and release boundaries.

## Reproduce

The reference run used Python 3.13.5, NumPy 2.3.5, SciPy 1.17.0 and SymPy 1.14.0. Create an isolated environment, install `requirements.txt`, and work in a copy of the bundle if you need to preserve the original manifest:

```sh
python -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
python code/run_all.py --regenerate
python code/check_consistency.py
```

`--regenerate` recomputes all 31 periodic approximations. Omitting it reuses the saved collocation polynomials but re-executes all algebra, projective/adjoint propagation, nonlinear first-return integrations, phase geometry, cross-method, reciprocal-rate, envelope and deliberate-error checks. No online service or API key is needed. The driver exits with non-zero status on a required failure. Reproduction updates the current receipts and their timing fields; new byte hashes may therefore differ even when the scientific values reproduce.

The scientific suite checks **50 exact equalities**, **31 periodic approximations**, **31 cycle audits**, **4 cross-method repetitions**, **12 reciprocal-rate comparisons**, **4 envelope checks** and **4 deliberate-error controls**. These are not independent statistical replicates. Counts and empirical discrepancies are generated from receipts, not maintained by hand in the numerical tables.

The BVP generator may legitimately return rejected attempts; the full run driver requires the stated 31-case reference outcome when reproducing this release. Every attempted solve is retained. Exponentials are never clipped in the active BVP solver: invalid guarded iterates abort. Numerical outputs are ordinary floating-point approximations, not interval certificates.

## Interpret the failures correctly

`evidence/exact_checks_initial_failure.json` records a 49/50 checker run. An off-constraint residual used an incorrect multiple of the tangent constraint O. The correct term is pO, not (1+eta)pO. The on-cycle identity remained unchanged. The initial implementation is preserved in `provenance/exact_checks_initial.py` and is **not** the active checker.

`evidence/deliberately_wrong_exact_receipt.json` intentionally has `passed=false`: the correct negative-control outcome is rejection of a deliberately false sign identity and exit status 1. Reversed phase order and wrong winding are tested with the same predicates used in the active numerical audit. A false radial contraction certificate is correctly rejected even though the cycle approximation is attracting.

The two v1.2 neutral-search receipts record failed historical solves. Their source code and meshes were not among the attached v1.2 material, so they were **not rerun**. They supply no exclusion of neutral cycles. The supplied review concerns v0.3, not this paper.

## Build and verify

See `BUILD.md` for PDF compilation. Run `python code/verify_manifest.py` before modifying the release to verify byte integrity. The manifest does not prove the analysis. Rendered-page checks and release consistency results are in `review/` and `evidence/release_checks.json`.
