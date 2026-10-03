# Wilhelm–Heinrich oscillator: convex attractor geometry, exact stability formulas and large-gain uniqueness

**Version 3.0 — 3 October 2026. Research candidate for consideration by Evidence Press.**

Start with `paper/manuscript.pdf` ({{PAGES}} pages). The paper is stand-alone: it restates the model, proves every statement it uses or names the imported result, and does not require v2.0.

## Exact claim boundary

**Proved.**

- *Uniqueness at large gain.* If r0 >= 1000 max(a0,b0)^2 / min(a0,b0) in the original rates (k >= 1000 for mu <= 1; k >= 1000 mu^3 for mu >= 1), every cycle is hyperbolic and attracting, with explicit bounds on its period, its slope at the up-crossing, the duration of its crest and its weak exponent. The cycle is unique and attracts every positive solution off the one-dimensional stable manifold of the equilibrium.
- *The attractor is a convex planar graph.* For all positive parameters the global attractor is the graph of a convex function over a planar set, C^1 in the interior. Above the Hopf value that set is the closed Jordan domain bounded by the outermost projected cycle, and the flow on the attractor is planar with divergence eta − mu. At any parameter where every cycle has ⟨eta⟩ < mu the cycle is unique.
- *Shape of every cycle.* The projection to the (U,p) plane is strictly convex and star-shaped; the weak exponent equals −⟨pK/Sigma⟩ = −⟨pW/Xi⟩ for two positive convexity gauges; the weak adjoint has a fixed pattern of twelve events per period.
- *The v2.0 framework*, with the referee's requested revisions: positive adjoint and exact exponent, supporting planes, six ordered extrema, uniform compactness of cycle families, global uniqueness in a non-explicit collar of each Hopf point, and the neutral-cycle obstruction.
- *Small rate ratio.* Exact identities; a kurtosis bound for strictly convex potentials; a local uniqueness theorem in the near-Hamiltonian wedge, conditional on one property taken from Fenichel's theorem.
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

{{REPRODUCE}}

## Evidence at a glance

{{CHECKS}}

Exact checks test algebraic identities. High-precision checks confirm the numerical constants of the large-gain proof. Floating-point checks test inequalities and formulas at finitely many parameters. None is an interval enclosure, and none formalises the non-algebraic steps of the proofs. Counts are not independent replications and are not converted into a probability of correctness.

## How the results were checked

Each new theorem was first proved in a research lane and then attacked by an independent referee agent that re-derived the proof and wrote its own numerical tests. The assembled sections were then refereed adversarially by a second model. Gaps found at either stage were repaired in the text; the dispositions are recorded in `paper/sections/*.report.md` and summarised in Section {{EVIDENCE_SECTION}} of the paper. This is internal adversarial review. It is not review by specialists in the field.

## Interpret the failures correctly

- Receipts with the suffix `_initial_failure` record a first run that failed because of an error in the checking script, not in the mathematics; the corrected run is the receipt without the suffix. Both are kept.
- Deliberate-error controls must fail; a control receipt with `passed: false` inside a passing control suite is the correct outcome.
- The v2.0 records keep their meaning: the 49-of-50 first exact run, the two failed historical neutral searches, and the failed radial certificate are preserved and are not evidence for or against uniqueness.
- "Refuted (floating point)" in Section 12 of the paper means that a candidate statement fails at stated parameters in ordinary floating-point computation.

## Open decisions

| Decision | Options | Note |
|---|---|---|
| Specialist review of the large-gain proof and the attractor theorem | Seek it before any public release; or release as a candidate with the present caveats | Both proofs survived two internal adversarial passes; neither has been read by a specialist. |
| The remaining region | Commission the programme of Section {{STATUS_SECTION}} (explicit Hopf collar, lower large-gain threshold, computer-assisted core, small-rate-ratio analysis); or leave the question open | The programme is a multi-stage project with one analytic piece (small rate ratio) of uncertain difficulty. |
| Literature access | Obtain the full texts listed in `SOURCES.md`, Section 6, and revise the positioning table | One 1981 paper on a third-order neuron model has not been seen and could bear on priority or method. |
| Scope of a first publication | One paper as assembled; or two (framework and attractor; large-gain theorem) | The assembled paper is long. |
