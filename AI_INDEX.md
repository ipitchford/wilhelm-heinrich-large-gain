# AI index: Wilhelm–Heinrich large-gain uniqueness

## Start here

Anonymous, version 3.1.0-candidate. DOI https://doi.org/10.5281/zenodo.23126164.
Read the [manuscript](paper/manuscript.pdf), [claim ledger](CLAIMS.json),
[status record](STATUS.json) and [review response](REVIEW_RESPONSE.md).
The ledger's “proved” means a written proof in this unrefereed candidate;
it does not mean external acceptance or machine-checked proof.

## Model and claims

For positive physical rates a0,b0,r0, the candidate proves one attracting
positive periodic orbit when r0 >= 1000 max(a0,b0)^2/min(a0,b0).
In normalized variables this is k >= 1000 max(1,mu^3).
It attracts positive solutions outside the equilibrium's stable manifold.
The main argument is [Section 10](paper/sections/s10_largegain.tex), using
the [adjoint multipliers](paper/sections/s02_adjoint.tex) and the
[attractor uniqueness criterion](paper/sections/s07_attractor.tex).

The positive-interior attractor is a graph of a convex function. Only above
Hopf is its base a closed Jordan disk; below/at Hopf it is the equilibrium.
This does not assert convexity of the graph as a subset of three-dimensional space.
The exact multiplier product alone does not prove both multipliers are below one.

The [local slow-graph lemma](paper/sections/s11_slowgraph.tex) discharges
the old assumption (F) using Eldering's theorem. The resulting
[wedge theorem](paper/sections/s11_wedge.tex) remains local to a fixed energy
ball and nonexplicit in its small-parameter threshold. The same section contains
a reusable kurtosis inequality for strictly convex potentials.

## Open and negative results

All-parameter cycle uniqueness remains open between the Hopf collars and
the explicit large-gain region. No second cycle was found; finite searches
cannot exclude one. Global small-rate-ratio uniqueness and r-to-zero uniformity
are not established. See [remaining programme](paper/sections/s14_status.tex)
and [failed certificates](paper/sections/s12_negative.tex).

## Evidence and reproduction

Use the pinned [environment](ENVIRONMENT.json) and [requirements](requirements.txt).
The [driver](code/run_all.py) with --regenerate recomputes the main suite,
including 31 saved cycles. It changes receipts and generated tables, so replay
in a fresh copy. The [run receipt](evidence/run_receipt.json) and
[evidence discussion](paper/sections/s13_evidence.tex) define coverage.
Some cartography tables and historical research scripts are not regenerated.

Run [rational_interval_audit.py](code/rational_interval_audit.py) separately
for outward rational enclosures of 56 selected endpoint comparisons;
its [receipt](evidence/rational_interval_results.json) is not an orbit enclosure.
Run [revision_controls.py](code/revision_controls.py) normally and with -O;
wrong derivatives, missing chain-rule terms, invalid domains and a false bound
must be rejected. Original semantic controls remain in the main suite.
Use [MANIFEST.sha256](MANIFEST.sha256) for byte integrity and
[BUILD.md](BUILD.md) for compilation. Integrity is not mathematical truth.

## Sources, provenance and assurance

[SOURCES.md](SOURCES.md) records access depth and unresolved Li (1981) and
other priority gaps. [RESEARCH_GATES.json](RESEARCH_GATES.json) separates
retrospective intake checking from the prospective local repair.
[F source mapping](review/F_SOURCE_MAPPING.md) and
[internal check](review/F_INTERNAL_CHECK.md) delimit the added lemma's audit.
Historical failed receipts remain failures, including the unresolved
small-mu quadrature comparison; successful alternative integration does not erase them.
Original prose/data use CC0-1.0 and original code MIT, with exclusions in
[LICENSES.md](LICENSES.md). No specialist peer review, independent external
reproduction, formal verification or historical priority is asserted.
