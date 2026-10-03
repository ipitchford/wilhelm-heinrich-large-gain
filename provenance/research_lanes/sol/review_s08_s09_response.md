### Defects

- Lemma \(\ref{lem:collar}\), Step 1, “Every solution that remains in a sufficiently small neighbourhood of \(E_{\rm ph}\) for all times lies on the local centre manifold” | **MINOR** | This is false if “for all times” means only all forward times: points on stable fibres can remain near the equilibrium and converge to the centre manifold without lying on it. The local centre manifold contains the maximal invariant set in a sufficiently small neighbourhood, hence in particular every periodic orbit contained there. | Replace the sentence by: “Every complete orbit contained in a sufficiently small neighbourhood, and hence every periodic orbit contained there, lies on the local centre manifold.”

- Lemma \(\ref{lem:collar}\), Step 1, use of the Hopf calculation | **GAP** | A negative first Lyapunov coefficient and merely “transversal crossing” do not by themselves identify on which side \(k>\kH\) the attracting branch lies. The asserted supercritical orientation requires the sign of \(\frac{d}{dk}\Re z(k)\) at \(k=\kH\). | Add the computation
  \[
  \frac{dz}{dk}=-\frac1{\chi_E'(z)},\qquad
  \Re\frac{dz}{dk}\Big|_{z=i\sqrt\mu,k=\kH}
  =\frac{1}{2\bigl(\mu+(1+\mu)^2\bigr)}>0.
  \]
  Together with \(l_1<0\), this puts the stable branch on \(k>\kH\).

- Lemma \(\ref{lem:collar}\), Step 2, “permanence and the competitive Poincaré–Bendixson property \cite[Theorems 1.3 and 1.4]{Smith12} give…” | **GAP** | The hypotheses of the imported results are not checked. In particular, the displayed physical system has a negative-feedback cyclic sign pattern, not literally the standard competitive sign pattern in the displayed coordinates. The proof must identify the relevant coordinate/sign transformation or invoke the precise monotone cyclic-feedback theorem, and it must establish the required precompactness/permanence and irreducibility assumptions. As written, the deduction about all \(\omega\)-limit sets is unsupported. The separate citation to BH23 repairs this only if that cited theorem really states global convergence for this exact system at \(k=\kH\). | State BH23’s theorem precisely and verify that its variables and hypotheses match, or supply the missing transformation and hypothesis checks for the cited Poincaré–Bendixson theorem.

- Lemma \(\ref{lem:collar}\), Step 2, “On a nonconstant cycle \(\eta\) is nonconstant” | **MINOR** | This is needed for strict Jensen but is asserted without proof. | If \(\eta\) were constant, the first projective equation would make \(\zeta\) constant; then the second would make \(x\) constant. Hence \(U,p,y\) are constant, contradicting that the orbit is a cycle.

- Theorem \(\ref{thm:neutral}\), “smooth dependence of the variational equation preserves the extra unit multiplier” | **MINOR** | Smooth convergence of monodromy matrices preserves an eigenvalue at \(1\), but every autonomous periodic orbit already has the trivial flow multiplier \(1\). The sentence must explain why the *transverse* unit multiplier persists rather than merely rediscovering the trivial one. The least-period convergence supplied by Theorem \(\ref{thm:compact}\) prevents period folding, but the transverse return map or algebraic multiplicity must be mentioned. | Use the crest return maps: crest transversality persists, their derivatives converge, and each neutral derivative has eigenvalue \(1\); therefore the limiting two-dimensional return derivative has eigenvalue \(1\). Alternatively, say that the full monodromy has algebraic multiplicity at least two at \(1\), and continuity of the transverse multipliers gives the conclusion.

- Theorem \(\ref{thm:neutral}\), “Each persists uniquely by the implicit-function theorem for the local return map” | **MINOR** | The common transverse section and persistence of the return map are left implicit. This is harmless but needed for a fully self-contained local-count argument. | Choose a local section through each crest. Since \(p'<0\) at a crest, the section remains transverse for nearby parameters; the return map is then \(C^1\), and \(I-DP\) is invertible because neither transverse multiplier equals \(1\).

- Corollary \(\ref{cor:neutralconditions}\), proof of the strict upper inequality in \(\eqref{eq:variancewindow}\), “a nonconstant continuous \(q\) cannot attain an extremal constant on every interval where \(q^2>0\)” | **MINOR** | This is imprecise and does not directly characterize equality in the averaged inequality. For the upper bound one uses \(q^3\le bq^2\), \(b=\eta_+-\mu\), and must rule out equality almost everywhere. | Equality would imply \(q(t)\in\{0,b\}\) whenever \(q(t)^2>0\). Continuity on the connected time circle then forces \(q\) to be constant unless it passes through values strictly between \(0\) and \(b\), where the inequality is strict. Mean zero and nonconstancy therefore force strictness. For the lower bound, \(q>-\min(\mu,1)\) gives pointwise strictness wherever \(q\ne0\).

- Equations \(\eqref{cpt:eq:V}\) and surrounding displays, use of `\(\dd\)` without an explicit variable | **MINOR** | The integration variable is typographically incomplete unless `\dd` is a macro that already includes \(dt\), which would be unusual given later use of `\(\mathrm d\sigma\)`. | Write \(\mathrm dt\).

### Verdicts

- Lemma \(\ref{lem:statebox}\): **HOLDS** — the \(R,S\) equations, extrema arguments, two-filter tail estimate, averaging identities, and logarithmic lower bound all re-derive correctly.
- Lemma \(\ref{lem:eqseparation}\): **HOLDS** — the spectral decomposition, \(V\)-bounds, identity \(DV(r)J_Er=|P_ur|^2\), cone estimate, nonlinear remainder bound, and minimum-of-\(V\) argument all check.
- Lemma \(\ref{lem:speed}\): **HOLDS** — both identities in \(\eqref{eq:speedrelations}\) and the resulting norm estimate are correct.
- Theorem \(\ref{thm:compact}\): **HOLDS** — the total-variation length estimate, uniform period bound, nonconstant periodic limit, crest preservation, and exclusion of a multiple limiting period all check.
- Lemma \(\ref{lem:collar}\): **GAP** — the conclusion depends on an unchecked imported global-stability/Poincaré–Bendixson assertion and omits the Hopf crossing orientation.
- Theorem \(\ref{thm:neutral}\): **HOLDS** — subject to the minor clarification about convergence of the transverse return multiplier.
- Corollary \(\ref{obs:cor:uniqueset}\): **HOLDS** — compactness, multiplier continuity, local count constancy, and the boundary argument are valid.
- Corollary \(\ref{cor:neutralconditions}\): **HOLDS** — the multiplier equivalences, trace factorization, cubic expansion, and variance bounds are correct, with only the strictness argument needing cleaner wording.
