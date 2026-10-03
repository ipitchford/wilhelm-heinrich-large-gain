<!-- Referee report on v2.0, supplied with the request that led to v3.0 (3 October 2026). Reproduced verbatim as provenance; it assesses v2.0, not this version. Links of the form sandbox:/mnt/data/... point to the referee's own workspace and are not part of this bundle. -->

# Review: *Finite-parameter obstructions and cycle geometry in the Wilhelm–Heinrich oscillator*, v2.0

**Editorial recommendation: Minor revisions.**  
**Prospective REF assessment: 3* — internationally excellent.**  
**Principal finding: I found no fatal mathematical defect in the stated obstruction theorem. The strongest contribution is the compactness theorem that makes the obstruction argument possible, rather than the continuation argument by itself.**

## 1. Executive assessment

This is a substantial theoretical paper with a carefully delimited claim. It does not prove uniqueness of the Wilhelm–Heinrich cycle throughout the oscillatory parameter region. Instead, it establishes detailed restrictions on every possible cycle and proves that failure of uniqueness at a finite parameter must be accompanied, somewhere along a finite fixed-rate path, by a nonconstant cycle with an additional transverse multiplier \(+1\).

The strongest result is Theorem 6.1: cycles cannot escape towards the chemical boundary, approach the equilibrium, or acquire unbounded least periods while parameters remain in a compact strictly oscillatory region. Its proof combines the supporting-plane inequalities, equilibrium cone exclusion, a uniform lower bound on projected speed, and the one-crest geometry. I checked this chain closely and found it mathematically convincing.

The computational evidence is also strong within its stated scope. The full supplied workflow reproduced, including regeneration of all 31 cycle approximations. A separately written check in the original concentration coordinates corroborated four return matrices and the failed radial certificate.

The remaining concerns concern precise novelty attribution, presentation of the uniformity arguments, and the extent of the contribution’s significance. These are not reasons to reject the paper or demand a uniqueness proof before publication. My recommendation is **minor revisions**, with a **3*** prospective REF rating and greater confidence in the mathematical assessment than in definitive historical priority.

## 2. Scope and evidence limits

I reviewed the manuscript.pdf[17-page manuscript](sandbox:/mnt/data/wh_review/wilhelm_heinrich_finite_parameter_v2_0/paper/manuscript.pdf), its LaTeX source, the active computational scripts, saved cycle meshes, numerical and symbolic receipts, source ledger, claim statements, and author-side audit and missingness documentation. The theorem and page references below refer to this v2.0 manuscript.

The assessment standard is that of a specialist research paper in dynamical systems or mathematical biology. No particular journal’s additional criteria were supplied.

Three kinds of evidence must remain distinct:

- **Analytic evidence:** the proofs in the manuscript, assessed by checking their equations, hypotheses, deductions and limiting arguments.
- **Computational evidence:** execution of the supplied workflow and additional independently implemented numerical checks.
- **External evidence:** primary literature and official REF guidance inspected through live searches on **2 October 2026**.

This review is not proof-assistant verification or an interval-arithmetic certificate. My independent numerical implementation reused supplied cycle seeds; it was not an independent discovery of those cycles. The historical referee material in the archive was not treated as an endorsement of v2.0.

Full texts of the original Wilhelm–Heinrich papers and the foundational Mallet-Paret–Smith paper were not obtained. Their publisher records and abstracts were checked; the accessible Smith and Boros–Hofbauer preprints supplied the model-specific comparison.

## 3. Contribution and positioning

### Strongest defensible thesis

The paper establishes that **finite-parameter multiplicity cannot arise through several otherwise plausible global degenerations: any failure of uniqueness must encounter a nonconstant neutral cycle with one transverse multiplier equal to \(+1\), while the other remains strictly contracting**.

In the original physical variables, the required transverse multipliers are

\[
1,\qquad \exp\!\left[-(a_0+b_0)T_{\mathrm{phys}}\right].
\]

This is a necessary obstruction, not a sufficient condition for multiplicity. The manuscript correctly allows the neutral event at the endpoint and does not assume a generic saddle-node bifurcation or two independent ambient unit eigenvectors.

The contribution separates into several layers:

| Result | Evidence in manuscript | Assessment |
|---|---|---|
| Positive adjoint and exact transverse exponent | Propositions 2.1–2.2 | Sound, useful formulation; general-theory antecedents need careful attribution. |
| Strict supporting-plane comparison | Theorem 3.1 and its consequences | A central structural ingredient with substantial mathematical value. |
| Simple projection, six-extrema order and ordered crests | Theorems 4.1–4.2; Corollary 4.3 | Convincing direct derivations; some underlying phenomena belong to established monotone-feedback theory. |
| Convex forward-invariant epigraph and rank-one Hessian | Theorem 5.1; Proposition 5.2 | Interesting additional structure, correctly distinguished from an invariant manifold. |
| Uniform compactness of all cycles | Theorem 6.1 | The strongest and most consequential part of the paper. |
| Closed neutral locus and preservation of cycle count | Lemma 7.1; Theorem 7.2 | Credible consequence of compactness, the global Hopf collar and hyperbolic continuation. |
| Failure of the square-root radial certificate | Theorem 8.1 | A useful analytic negative result, not a counterexample to uniqueness. |

The established literature already supplies permanence, competitive dynamics, existence of an attracting cycle above Hopf, and finiteness of the number of cycles. Smith explicitly left their number unresolved; Boros and Hofbauer retained that question in their later account. The present paper should therefore be credited for its additional geometry and parameter-uniform control—not for rediscovering periodic behaviour or excluding chaos. [arXiv](https://arxiv.org/html/1109.5204v1)

## 4. Major comments

**I identified no fatal defect or major mathematical repair required.** The following are the principal substantive comments; the requested revisions are bounded.

### 4.1. Make the novelty claim theorem-specific

**Issue.** The manuscript acknowledges familiar monotone-feedback geometry and incomplete priority verification, but the reader still needs a sharper account of precisely which statements are being advanced as new.

**Manuscript evidence.** Section 1.3 explicitly distinguishes established inputs and acknowledges that direct derivation does not establish historical priority. The source ledger is similarly candid.

**Why it matters.** Planar projection, ordering phenomena and Floquet restrictions have a substantial theoretical background. Mallet-Paret and Smith’s foundational paper establishes a Poincaré–Bendixson theory for monotone cyclic feedback systems. Its existence does not establish that it already contains the present results, but it makes a theorem-level comparison necessary before attaching strong originality claims to the whole package. [Springer](https://link.springer.com/article/10.1007/BF01054041)

**Required revision.** Add a compact positioning table separating:

1. established model-specific results;
2. direct specialisations or new proofs of general phenomena;
3. the supporting-plane, envelope and uniform compactness statements claimed as new.

Where an older full text remains inaccessible, retain that limitation rather than implying that a negative search establishes priority. This is the main qualification on my originality rating.

### 4.2. Give the uniform compactness argument greater prominence

**Issue.** The most important proof occupies a relatively compressed section among numerous structural results.

**Manuscript evidence.** Theorem 6.1, pp. 9–10, proves the state box, equilibrium separation, speed bound, period bound and preservation of least period in five steps.

**Why it matters.** A bounded family of trajectories is not automatically a compact family of nonconstant periodic orbits with bounded least periods. The manuscript supplies the additional arguments; those arguments deserve to be unmistakable.

**Required revision.** Extract the equilibrium-separation argument and the projected-speed argument into named intermediate lemmas, or give them a short dependency summary. In particular, collect the constants establishing

\[
|r|\le C_Q|P_ur|,
\qquad
V'\ge c_QV>0,
\qquad
|r|\le C_Q'\sqrt{p^2+(p')^2},
\]

and explicitly state that their constants depend only on the compact parameter set \(Q\), not on the particular cycle or its period.

The necessary reasoning is already present. This is a request to expose the proof’s central mechanism, not an allegation that uniformity is missing.

### 4.3. Make the global Hopf-collar localisation completely explicit

**Issue.** Local Hopf uniqueness alone would not justify the starting cycle count used in the obstruction theorem.

**Manuscript evidence.** Lemma 7.1, p. 11, correctly goes further: it uses a common positive box in physical coordinates, critical global attraction, finite-time parameter dependence and invariance of any periodic orbit.

**Why it matters.** This is the step that excludes a remote cycle coexisting with the small Hopf cycle near threshold. Without it, the final continuation argument would start from an unjustified global count.

**Required revision.** Introduce two nested local neighbourhoods explicitly: a smaller one containing the image of the common box at the critical parameter, and a larger one in which the local Hopf classification applies. Then spell out how parameter continuity and cycle invariance place every nearby cycle in the larger neighbourhood.

I consider the existing argument sound. The requested addition would make an important local-to-global passage easier to verify independently.

## 5. Rigour, results and inference

### 5.1. Model, scaling and parameter coverage

The transformations preserve the original filter order. In particular, the intermediate filter is \(Z\), the feedback filter is \(Y\), and the final theorems do not silently restrict \(\mu=a_0/b_0\) to \(\mu\le1\).

This matters because reciprocal-rate equivalence of a scalar formulation is not the same as preserving the labels of the chemical concentrations.

I checked the reduction to

\[
U'''+(1+\mu)U''+\mu U'+k(e^U-1)=0
\]

and the conversion back to physical time in the neutral multiplier. I found no discrepancy. The original concentration scaling agrees with Smith’s model formulation. [arXiv](https://arxiv.org/html/1109.5204v1)

### 5.2. Adjoint and Floquet identities

The positive adjoint construction does not assume cycle stability. That is essential: otherwise its later use to constrain potentially unstable cycles would be circular.

The multiplier formulas

\[
m_s=e^{-(1+\langle\eta\rangle)T},
\qquad
m_t=e^{(\langle\eta\rangle-\mu)T}
\]

are consistent with

\[
m_sm_t=e^{-(1+\mu)T}.
\]

The identification of the potentially unstable transverse exponent as

\[
\lambda_c=\langle\eta\rangle-\mu
\]

is useful and correctly preserves the unresolved sign.

I also checked the role of the strengthened bound

\[
\eta>\max(0,\mu-1).
\]

It is important for the all-positive-rate argument and for convexity of the cubic gain polynomial. The paper does not incorrectly infer stability above Hopf from that convexity.

### 5.3. Supporting planes: a genuinely important step

Theorem 3.1 reduces the comparison of two bounded entire solutions to

\[
\Delta'+(1+\eta)\Delta
=
-x_1\left(e^{U_2-U_1}-1-(U_2-U_1)\right).
\]

The right-hand side is nonpositive, and integration from the remote past gives the supporting inequality. The strictness argument uses uniqueness of solutions when the exponential remainder vanishes identically; it does not require that remainder to be positive at every instant.

That distinction is handled correctly.

The two equilibrium comparisons are also separately needed. A one-sided comparison alone would not justify excluding both cones around the stable eigendirection. The manuscript supplies both directions and subsequently uses them in the compactness proof.

### 5.4. Geometry and return-map restrictions

The projection argument is stronger than a numerical phase portrait. Two-way support excludes self-intersections and intersections between distinct projected cycles. The crossing identities and winding argument then give a unique crest and trough.

The six-extrema order concerns the original labelled concentrations:

\[
X_{\min}
<
Z_{\min}
<
Y_{\min}
<
X_{\max}
<
Z_{\max}
<
Y_{\max}
<
X_{\min}+T,
\]

where the inequalities order event times.

The positive crest-return matrix and the period estimate yield

\[
\det M=e^{-(1+\mu)T}<e^{-4\pi}
\approx 3.4873\times10^{-6}.
\]

This is strong **area contraction**, not a proof that both transverse directions contract. The manuscript correctly avoids that inference. It also distinguishes the positivity restrictions from an unproved global contraction criterion.

### 5.5. The convex envelope is not overinterpreted

The envelope is constructed from the supporting planes of selected cycles and the equilibrium. The maximum-of-affine-functions construction gives convexity and the relevant one-sided invariance.

The local Hessian formula is credible: uniqueness of the supporting phase and the negative second phase derivative permit the implicit-function argument, producing a positive semidefinite rank-one Hessian.

Two limitations are correctly retained:

**First**, the envelope depends on actual selected cycles; it is not an independent algorithm for locating all cycles.

**Second**, its graph has a nonzero invariance defect away from the selected trajectories. Consequently, the associated planar field does not establish a global equivalence between arbitrary planar and spatial trajectories.

These qualifications are mathematically substantive, not merely defensive wording.

### 5.6. Compactness and the obstruction theorem

Theorem 6.1 survives the principal objections I tested.

The positive state box is uniform on compact positive parameter sets. The explicit lower bound can be extremely small, but it remains positive; practical sharpness is not required for the theorem.

Equilibrium separation is not inferred merely from instability. Cone exclusion prevents cycle points from aligning with the stable eigendirection. The quadratic form on the unstable plane then excludes the minimum of that form from a sufficiently small equilibrium neighbourhood.

The tangent identities convert this separation into a positive lower bound on projected speed. The one-maximum/one-minimum property bounds total variation and hence projected length. Dividing length by the speed lower bound gives the period bound.

Finally, transverse crest crossings prevent a limiting period from becoming a multiple traversal of a smaller cycle.

These arguments establish the compactness needed for the neutral locus to be relatively closed and for hyperbolic cycles to persist without unseen additional branches appearing nearby. The final obstruction follows. It remains a necessary obstruction: neither the mean criterion nor the variance window excludes neutrality.

### 5.7. Reproducibility audit

The original archive passed its byte-integrity check: **121 of 121 manifest entries matched**.

I then ran the complete workflow in a working copy, including fresh regeneration:

| Check | Observed outcome |
|---|---:|
| Symbolic equalities | 50/50 passed |
| Regenerated cycle approximations | 31/31 accepted |
| Cycle geometry and adjoint audits | 31/31 passed |
| Cross-method checks | 4/4 passed |
| Reciprocal-rate comparisons | 12/12 passed |
| Envelope checks | 4/4 passed |
| Negative controls | 4/4 passed |
| Final manuscript/evidence consistency checks | 107/107 passed |

The reproduced maximum discrepancy between the two adjoint-exponent formulations was approximately \(3.53\times10^{-12}\). The largest reported return-closure error was approximately \(2.09\times10^{-11}\).

The deliberately wrong symbolic check failed as intended. The preserved initial failed check is explained rather than suppressed. Those are useful features of the audit trail.

Nevertheless, the 31 accepted solutions come from **amplitude-fixed boundary-value solves**. They are not an exhaustive parameter-space survey, and the reciprocal-rate comparisons are symmetry checks rather than independent discoveries.

### 5.8. Independently implemented checks

I wrote an additional implementation using

\[
X'=\lambda X(1-Y),\qquad
Y'=\mu(Z-Y),\qquad
Z'=X-Z,
\]

without importing the supplied code. It integrated the concentration equations and their variational system, applied the event-time correction at the first crest return, and checked the resulting derivative against nonlinear finite differences.

Four cases were tested, including small amplitude, larger amplitude, equal rates and rate ratios on both sides of one.

The largest absolute difference from the supplied return matrices was **\(3.22\times10^{-9}\)**. At the common finite-difference step \(5\times10^{-6}\), the largest difference from the independently integrated variational derivative was **\(2.76\times10^{-8}\)**.

The small-amplitude example also reproduced

\[
\frac{(Mv)_1}{v_1}\approx1.0000581423>1,
\qquad
m_t\approx0.9999882842<1.
\]

Thus the proposed radial certificate fails even though that cycle is attracting. This corroborates the interpretation of Theorem 8.1.

The script, inputs, results, workflow receipts and logs are available in the wilhelm_heinrich_v2_referee_checks_2026-10-02.zip[referee computational-check bundle](sandbox:/mnt/data/wilhelm_heinrich_v2_referee_checks_2026-10-02.zip). These remain floating-point checks, not rigorous numerical enclosures.

## 6. External literature check

### Search record and limitations

The live search used web discovery, arXiv author preprints, Springer publisher records, ScienceDirect records, the AIMS journal site and official REF publications. Searches combined the model names with “periodic-orbit uniqueness”, “limit cycles”, “supporting planes”, “convexity” and “neutral cycles”, alongside exact-title searches and searches of monotone cyclic feedback theory.

I prioritised the original equations and directly relevant theorems. Results about modified, homogenised, forced or mass-conserving systems were not treated as results about the original oscillator. This was a bounded search, not a systematic priority certificate.

### Decision-relevant comparisons

**Wilhelm and Heinrich, 1995 and 1996.** The publisher abstracts confirm the model’s introduction and the local bifurcation analysis, including a stable cycle near Hopf. They do not provide a basis for asserting that the present global compactness theorem is already known. Full-text priority comparison remains incomplete. [Springer](https://link.springer.com/article/10.1007/BF01165134)

**Smith, 2012.** The accessible preprint establishes the global dynamical framework and explicitly leaves the above-Hopf cycle count open. Crucially, convergence to some cycle is not uniqueness of that cycle. The manuscript’s use of Smith respects this distinction. [arXiv](https://arxiv.org/html/1109.5204v1)

**Boros and Hofbauer, 2023.** Theorem 13 restates the relevant original-system conclusions, including critical global equilibrium stability, and again identifies the cycle count as unresolved. Its subsequent multicycle phenomena concern a different mass-conserving extension; they do not contradict the present paper. [arXiv](https://arxiv.org/html/2202.11034v2)

**Mallet-Paret and Smith, 1990.** This is the essential general-theory comparator. The publisher abstract verifies the scope of its monotone cyclic feedback theory, but I did not obtain sufficient full text to settle the priority of each geometric refinement here. [Springer](https://link.springer.com/article/10.1007/BF01054041)

**Van den Berg and Queirolo, 2021.** Their validated-continuation framework illustrates the distinction between numerically following periodic solutions and rigorously validating a branch. It provides methodological context for possible future computational work, not evidence that the present diagnostic grid is an exhaustive exclusion. [AIMS Sciences](https://www.aimsciences.org/article/doi/10.3934/jcd.2021004)

**Search conclusion:** I found no inspected source establishing an all-parameter uniqueness theorem for the original equations or directly superseding the complete support–compactness argument. That supports the paper’s relevance, but does not settle historical novelty.

## 7. Minor comments

**Coordinate labels for return matrices.** Distinguish the matrices represented in \((U,y)\) and \((U,z)\), with \(z=y/k\), by subscripts or an explicit similarity transformation. Their eigenvalues agree, but their entries and componentwise test vectors are coordinate-dependent.

**Tables.** Give the numerical tables captions and numbers, and identify whether each period is in fast time, rescaled time or physical time. State the amplitude-fixed sampling design in the relevant caption.

**Notation.** The coexistence of \(\lambda\) as a model parameter and \(\lambda_c\) as an exponent is manageable, but a compact notation table would reduce unnecessary checking.

**Terminology.** Preserve the distinction between a periodic orbit, a numerical approximation to it, and a validated numerical enclosure. The existing text generally does this well.

**Presentation.** A single schematic showing the projected cycle, its supporting plane and the one-sided envelope would help readers understand the geometry. This is optional and should not replace any argument.

## 8. Prioritised revision plan

### Must fix before final publication

Make the originality statement theorem-specific and consistent with what the inspected literature establishes. Preserve the present scope: a finite-parameter obstruction and structural cycle theory, not a solution of the uniqueness conjecture. I identified no central mathematical repair that must precede publication.

### Should fix to strengthen the paper

Promote the equilibrium-separation and speed-bound arguments into clearly identifiable lemmas or a dependency summary. Make the nested-neighbourhood step in the global Hopf collar explicit. Clarify return-map coordinates and numerical table captions.

### Could improve presentation or future work

A transferable theorem for a broader class of convex negative-feedback systems would materially strengthen the significance case. Rigorous exclusion of neutral cycles on a nontrivial parameter region would also be valuable. Validating one continued branch, however, would still require a separate argument to exclude other branches.

These are research extensions, not conditions for accepting the present scoped results.

## 9. Editorial recommendation

**Minor revisions.**

The paper merits publication as a structural result about the original oscillator. Its main arguments are coherent, its use of established inputs is explicit, and its computations reproduce. The negative results and missing historical evidence are reported responsibly.

My confidence is **moderate to high** concerning the validity of the stated theorem chain, and **moderate** concerning originality.

The decision would change if a specialist identified a failure of the strict-support-to-equilibrium-separation argument, or if an older general theorem already supplied the claimed new compactness and envelope results under the same hypotheses. Neither issue was established in this review.

I would not require the authors to solve the all-parameter uniqueness problem before publishing this paper. Conversely, the present work should not be described as having effectively solved it: the admissibility or exclusion of neutral cycles remains substantive.

## 10. Provisional REF calibration

The most appropriate assessment area is **UoA 10: Mathematical Sciences**, within Main Panel B. The primary contribution is mathematical, despite the chemical model. [REF 2029](https://2029.ref.ac.uk/panels/main-and-sub-panels/main-panel-b/)

The current REF 2029 guidance states that detailed panel criteria and the quality-rating system will be published subsequently. Accordingly, this calibration uses the established REF 2021 output definitions and Main Panel B descriptors as an explicitly provisional benchmark—not as final REF 2029 rules. [REF 2029](https://2029.ref.ac.uk/guidance/section-1-overview/)

| Dimension | Provisional rating | Rationale | Confidence |
|---|---|---|---|
| **Originality** | **3*** | The support–envelope–compactness combination appears to offer important additional structure. Exact separation from older general theory remains incomplete. | Moderate |
| **Significance** | **3*** | It removes several global failure mechanisms and precisely identifies the remaining finite-parameter obstruction. The resulting influence is plausible and substantial, but not yet convincingly world-leading. | Moderate |
| **Rigour** | **3*** | Careful hypotheses, attention to degeneracies and limiting periods, convincing analytic arguments, reproducible diagnostics and effective negative controls. | Moderate–high |
| **Overall output quality** | **3*** | An internationally excellent specialist contribution, judged holistically rather than by averaging the dimensions. | Moderate |

### Why not 4* at present?

Main Panel B’s established four-star descriptors look for such features as exceptional novelty, field-leading results or the capacity to exert major influence on a research theme. The three-star descriptors encompass important knowledge and techniques with a likely lasting influence without requiring fundamental new concepts. [REF 2021](https://2021.ref.ac.uk/media/1450/ref-2019_02-panel-criteria-and-working-methods.pdf)

This paper currently fits the latter description more securely.

The obstacle is **not simply that the original conjecture remains open**. A partial result can be world-leading. Here, however, the evidence does not yet establish that the compactness and convex-support method has the novelty or wider consequences needed for that stronger judgement. The final continuation step is comparatively standard once the difficult compactness result has been secured.

A credible four-star case could arise from a demonstrably reusable theorem with major consequences beyond this example, or from a decisive classification or exclusion of the neutral mechanism. More numerical examples and more polished packaging alone would not establish that case.

**My present assessment is a strong 3* candidate, not a defensible 4* claim.** This is an indicative output-quality judgement, not an official panel decision or a determination of submission eligibility.

## 11. References

### Model and global dynamics

Wilhelm, T., & Heinrich, R. (1995). Smallest chemical reaction system with Hopf bifurcation. *Journal of Mathematical Chemistry, 17*, 1–14. DOI: 10.1007/BF01165134. Publisher abstract and metadata consulted. [Springer](https://link.springer.com/article/10.1007/BF01165134)

Wilhelm, T., & Heinrich, R. (1996). Mathematical analysis of the smallest chemical reaction system with Hopf bifurcation. *Journal of Mathematical Chemistry, 19*, 111–130. DOI: 10.1007/BF01165179. Publisher abstract and metadata consulted. [Springer](https://link.springer.com/article/10.1007/BF01165179)

Smith, H. L. (2012). Global dynamics of the smallest chemical reaction system with Hopf bifurcation. *Journal of Mathematical Chemistry, 50*, 989–995. DOI: 10.1007/s10910-011-9946-9. Author preprint consulted. [arXiv](https://arxiv.org/html/1109.5204v1)

Boros, B., & Hofbauer, J. (2023). Some minimal bimolecular mass-action systems with limit cycles. *Nonlinear Analysis: Real World Applications, 72*, Article 103839. DOI: 10.1016/j.nonrwa.2023.103839. Author preprint consulted, especially Theorem 13 and the subsequent discussion. [arXiv](https://arxiv.org/html/2202.11034v2)

### General theory and numerical validation

Mallet-Paret, J., & Smith, H. L. (1990). The Poincaré–Bendixson theorem for monotone cyclic feedback systems. *Journal of Dynamics and Differential Equations, 2*, 367–421. DOI: 10.1007/BF01054041. Publisher abstract and metadata consulted. [Springer](https://link.springer.com/article/10.1007/BF01054041)

Van den Berg, J. B., & Queirolo, E. (2021). A general framework for validated continuation of periodic orbits in systems of polynomial ODEs. *Journal of Computational Dynamics, 8*(1), 59–97. DOI: 10.3934/jcd.2021004. Publisher abstract and article information consulted. [AIMS Sciences](https://www.aimsciences.org/article/doi/10.3934/jcd.2021004)

### Assessment framework

Research Excellence Framework. (2019). *Panel criteria and working methods* (REF 2019/02), especially paragraphs 190–202 and the UoA 10 descriptor. [REF 2021](https://2021.ref.ac.uk/media/1450/ref-2019_02-panel-criteria-and-working-methods.pdf)

Research Excellence Framework. (2022). *Guidance on REF 2021 results*. [REF 2021](https://2021.ref.ac.uk/guidance-on-results/guidance-on-ref-2021-results/index.html)

Research Excellence Framework. (n.d.). *Section 1—Overview: REF 2029 guidance*. Retrieved October 2, 2026. [REF 2029](https://2029.ref.ac.uk/guidance/section-1-overview/)
