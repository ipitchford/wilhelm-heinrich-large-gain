# Wilhelm–Heinrich oscillator: convex attractor geometry, exact stability formulas and large-gain uniqueness

**Version 3.1.0-candidate — 3 October 2026. Anonymous. Unrefereed research candidate.**

DOI: https://doi.org/10.5281/zenodo.23126164

See [AI_INDEX.md](AI_INDEX.md) for claim and replay boundaries.

Start with `paper/manuscript.pdf` (compiled from the accompanying LaTeX). The paper is stand-alone: it restates the model, proves every statement it uses or names the imported result, and does not require v2.0.

## Exact claim boundary

**Proved.**

- *Uniqueness at large gain.* If r0 >= 1000 max(a0,b0)^2 / min(a0,b0) in the original rates (k >= 1000 for mu <= 1; k >= 1000 mu^3 for mu >= 1), every cycle is hyperbolic and attracting, with explicit bounds on its period, its slope at the up-crossing, the duration of its crest and its weak exponent. The cycle is unique and attracts every positive solution off the one-dimensional stable manifold of the equilibrium.
- *The attractor is a convex planar graph.* For all positive parameters the global attractor is the graph of a convex function over a planar set, C^1 in the interior. Above the Hopf value that set is the closed Jordan domain bounded by the outermost projected cycle, and the flow on the attractor is planar with divergence eta − mu. At any parameter where every cycle has ⟨eta⟩ < mu the cycle is unique.
- *Shape of every cycle.* The projection to the (U,p) plane is strictly convex and star-shaped; the weak exponent equals −⟨pK/Sigma⟩ = −⟨pW/Xi⟩ for two positive convexity gauges; the weak adjoint has a fixed pattern of twelve events per period.
- *The v2.0 framework*, with the referee's requested revisions: positive adjoint and exact exponent, supporting planes, six ordered extrema, uniform compactness of cycle families, global uniqueness in a non-explicit collar of each Hopf point, and the neutral-cycle obstruction.
- *Small rate ratio.* Exact identities; a kurtosis bound for strictly convex potentials; a local uniqueness theorem in the near-Hamiltonian wedge, with the formerly assumed slow-graph property (F) now derived from Eldering's precisely cited persistence theorem. This remains local in a fixed energy ball and has no explicit threshold.
- *Negative results.* Several natural certificates cannot work; each is stated with a proof or with the parameters of its floating-point refutation.

**Not proved.** Uniqueness at every oscillatory parameter. The region between the Hopf collar and the large-gain threshold is open: no neutral cycle is excluded there and none is exhibited. No multiple-cycle example was found. The statement that would close the question — ⟨eta⟩ < mu on every cycle — is unproved, as is the monotonicity statement that survived every numerical test. Global uniqueness for small rate ratio is not proved. Priority is not established, and the new arguments have not been reviewed by specialists.

## Contents

| Path | Content |
|---|---|
| `paper/` | `manuscript.pdf`; editable LaTeX (`manuscript.tex`, `preamble.tex`, `sections/`, `bibliography.tex`, `bib/`); generated tables (`generated/`); Figure 1 and its script (`figures/`); a report per section (`sections/*.report.md`) recording sources, corrections and checks. |
| `code/` | The v2.0 suite, unchanged in behaviour, and the v3.0 checks (`code/v3/`); the run driver, table generator and consistency checker. |
| `evidence/` | JSON receipts of every check; the 31 saved v2.0 cycle approximations and profiles; v3.0 receipts and the cartography tables (`evidence/v3/`). |
| `logs/` | Logs of the reference run. |
| `provenance/v2_0_bundle/` | The v2.0 bundle, byte for byte, with its own provenance. |
| `provenance/research_lanes/` | Working records of the nine research lanes behind v3.0: notes, proofs as first written, scripts and outputs, the adversarial referees' records (`VERIFY.md`), the cross-model review prompts and replies (`sol/`), and the checks made while assembling the paper (`framework_checks/`). |
| `review/` | The referee report on v2.0 (verbatim), the response to it, the author-side audit, and the list of changes and missing items. |
| `CLAIMS.json`, `STATUS.json`, `SOURCES.md`, `BUILD.md`, `ENVIRONMENT.json`, `requirements.txt`, `MANIFEST.sha256` | Claim ledger with numbered statements, status flags, source ledger with access levels, build and environment records, byte manifest. |

## Reproduce

Python 3.13 with the packages in `requirements.txt` (exact versions in `ENVIRONMENT.json`). Work in a copy of the bundle if the manifest must stay valid: a run rewrites receipts, logs and generated tables.

```sh
python -m venv .venv && . .venv/bin/activate
python -m pip install -r requirements.txt
python code/run_all.py              # every check, on the 31 saved cycles
python code/run_all.py --regenerate # also recompute the 31 saved cycles
python code/run_all.py --quick      # skip the cartography recomputation
python code/check_consistency.py    # receipts, tables, macros and references agree
python code/verify_manifest.py      # byte integrity against MANIFEST.sha256
```

No network access or API key is needed. Building the PDF is described in `BUILD.md`.

## Evidence at a glance

| Check | Kind | Passed / attempted |
|---|---|---|
| Framework identities | exact | 50 / 50 |
| Periodic boundary-value solves | floating point | 31 / 31 |
| Cycle geometry and adjoint audits | floating point | 31 / 31 |
| DOP853–Radau repetitions | floating point | 4 / 4 |
| Reciprocal-rate comparisons | floating point | 12 / 12 |
| Envelope Hessian checks | floating point | 4 / 4 |
| Deliberate-error controls (framework) | control | 4 / 4 |
| New algebraic identities | exact | 148 / 148 |
| Gauges and exponent formulas on cycles | floating point | 66 / 66 |
| Constants of the large-gain proof | high precision | 141 / 141 |
| Large-gain lemmas on cycles and on non-periodic segments | floating point | 29,561 / 29,561 |
| Attractor: supporting planes, gradient, crest-axis formula | floating point | 156 / 156 |
| Weak adjoint: event order, dual planes, rank-one identity | floating point | 180 / 180 |
| Kurtosis bound, period identities, third-order formula | high precision and floating point | 211 / 211 |
| Cartography: shipped tables and recomputed subsample | floating point | 140 / 140 |
| Refutations reproduced | floating point | 123 / 123 |
| Deliberate-error controls | control | 6 / 6 |

The reference run of `code/run_all.py` passed in 7.4 minutes; its receipt is `evidence/run_receipt.json`. A supplementary run of the large-gain lemma script at mu = 0.02, k = 1000, outside the default run, is kept with one comparison recorded as failed: two numbers of size about 3.26e9 that differ by about 1e-3, below the resolution of the quadrature; a separate integrator covers the same pair and passes; this does not turn the unresolved comparison into a pass (`evidence/v3/v3_largegain_lemmas_small_mu_supplement.json`).

Exact checks test algebraic identities. High-precision checks confirm the numerical constants of the large-gain proof. Floating-point checks test inequalities and formulas at finitely many parameters. The revision adds outward rational interval enclosures for 56 selected scalar comparisons. These do not certify trajectories, the analytic monotonicity reductions, or the non-algebraic proof steps. Counts are not independent replications and are not converted into a probability of correctness.

## How the results were checked

Each new theorem was first proved in a research lane and then attacked by a separately instructed referee agent that re-derived the proof and wrote its own numerical tests. The assembled sections were then refereed adversarially by a second model. Gaps found at either stage were repaired in the text; the dispositions are recorded in `paper/sections/*.report.md` and summarised in Section 13 of the paper. This is internal adversarial review. It is not review by specialists in the field.

## Interpret the failures correctly

- Receipts with the suffix `_initial_failure` record a first run that failed because of an error in the checking script, not in the mathematics; the corrected run is the receipt without the suffix. Both are kept.
- Deliberate-error controls must fail; a control receipt with `passed: false` inside a passing control suite is the correct outcome.
- The v2.0 records keep their meaning: the 49-of-50 first exact run, the two failed historical neutral searches, and the failed radial certificate are preserved and are not evidence for or against uniqueness.
- "Refuted (floating point)" in Section 12 of the paper means that a candidate statement fails at stated parameters in ordinary floating-point computation.

## Open decisions

| Decision | Options | Note |
|---|---|---|
| Specialist review of the large-gain proof and the attractor theorem | Seek it before any public release; or release as a candidate with the present caveats | Both proofs survived two internal adversarial passes; neither has been read by a specialist. |
| The remaining region | Commission the programme of Section 14 (explicit Hopf collar, lower large-gain threshold, computer-assisted core, small-rate-ratio analysis); or leave the question open | The programme is a multi-stage project with one analytic piece (small rate ratio) of uncertain difficulty. |
| Literature access | Obtain the full texts listed in `SOURCES.md`, Section 6, and revise the positioning table | One 1981 paper on a third-order neuron model has not been seen and could bear on priority or method. |
| Scope of a first publication | One paper as assembled; or two (framework and attractor; large-gain theorem) | The assembled paper is long. |
