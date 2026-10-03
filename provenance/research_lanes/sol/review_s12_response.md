- **Opening paragraph, “Every certificate of comparison type must agree with the true adjoint ratios to first order at the Hopf point.”** | **GAP** | Lemma \(\ref{lem:rigidity}\) proves this only for trial pairs admitting the specific first-order expansion in part (iii), and for state-dependent pairs with the regularity assumed in part (iv). It does not cover arbitrary “comparison type” certificates, for example nonsmooth or singular pairs, pairs depending nonlocally on the orbit, or pairs without an asymptotic expansion. | Restrict the claim to the classes covered by Lemma \(\ref{lem:rigidity}\).

- **After Proposition \(\ref{neg:prop:charts}\), “For a fixed chart the failure is confined to a neighbourhood of the Hopf point.”** | **GAP** | The proposition proves that failure occurs sufficiently near Hopf when \(c(a,b)\ne0\). It does not prove that failure cannot also occur away from Hopf. The floating-point tests cited afterward do not establish confinement either. | Replace by “For a fixed chart, the proposition proves failure in a neighbourhood of the Hopf point”; any assertion that failure is confined there must be removed or proved globally.

- **Relation to concave monotone maps, quoted Smith theorem** | **MINOR** | The theorem is quoted second-hand from a preprint, while the original source was not checked. Although it is explicitly not used in a proof, the exact hypotheses and conclusion therefore have not been independently verified here. | Quote the precise accessible source as the source of the stated formulation, or verify the original theorem.

- **Proposition \(\ref{prop:nodulac}\)(i), “a Dulac function for the reduced flow can only have the sign \(R_\Lambda\le0\)”** | **MINOR** | What is proved is only that no \(C^1\) \(\Lambda\) on an open set containing all of \(\Dsc\setminus\{0\}\) can satisfy \(R_\Lambda\ge0\) there. The unqualified wording could be read as excluding positive-sign Dulac functions on smaller subdomains, which the argument does not do. | Say “a global Dulac function on the punctured attractor domain can only have the nonpositive sign.”

- **Proposition \(\ref{prop:nodulac}\)(iii), choice of \(N_1=\{Q_E\le r_1^2\}\)** | **MINOR** | Strictly, \(N_1\) is a closed ellipse rather than a neighbourhood in the usual open-set sense, and the proof needs the whole closed ellipse to lie in \(N\cap\inter\Dsc\). This is clearly obtainable but should be stated. | Choose \(r_1>0\) so that \(\{Q_E\le r_1^2\}\subset N\cap\inter\Dsc\), and refer to it as a compact sublevel set.

- **Proposition \(\ref{prop:nodulac}\)(iii), existence of the first exit time \(t_q\)** | **MINOR** | The proof states that the lower exponential bound for \(Q_E(\chi_tq)\) implies a first time at which \(Q_E=r_1^2\), but does not explicitly rule out finite-time loss of existence. This is harmless because the orbit remains in a compact subset of \(\inter\Dsc\) until exit. | Add that the reduced flow exists as long as it remains in the compact ellipse, so the exponential growth forces finite exit.

- **Proposition \(\ref{prop:nomonomial}\), “with equality if and only if \((c_1,c_2,c_3)\) is proportional to \((\mu,\eta_E,\zeta_E)\)”** | **MINOR** | The proof establishes this, but the converse uses both identities
  \[
  \frac{\eta_E}{\mu}=\frac{\gamma_E-1}{\mu},\qquad
  \frac{\zeta_E}{\eta_E}=\gamma_E-\mu,
  \]
  which are only implicit. | State these identities explicitly when checking the equality case.

- **Lemma \(\ref{lem:rigidity}\)(i), Step 1, “the solutions … bounded on \(\mathbb R\) are exactly \(c\,r_1(\cdot+\phi)\) with \(c\ge0\)”** | **MINOR** | The zero solution is not represented by a phase in a meaningful way, and a general nonzero center-plane solution can indeed be written in that form only after absorbing the sign into a phase shift. The subsequent limit is nonzero, so this causes no damage. | Say “the nonzero bounded solutions are \(c\,r_1(\cdot+\phi)\), \(c>0\), together with the zero solution.”

- **Lemma \(\ref{lem:rigidity}\)(i), Step 2, passage from subsequential convergence to \(T_j\to T_0\)** | **GAP** | The argument says that comparison of consecutive crest zeros gives \(T_j\to T_0\), but initially no upper bound on \(T_j\) has been established. One must first use convergence on a fixed interval containing at least two consecutive negative-slope zeros of the limiting \(p_1\) to force the approximating function to have corresponding zeros; because the cycle has exactly one negative-slope zero per period, their distance is \(T_j\). The intended repair is local but currently omitted. | Explicitly choose a fixed interval containing two consecutive negative-slope zeros of the limit, use \(C^1\) convergence and simplicity to obtain corresponding zeros of \(\hat p_j\), and invoke the “one crest per period” theorem to identify their separation with \(T_j\).

- **Lemma \(\ref{lem:rigidity}\)(ii), Step 1, use of Proposition \(\ref{prop:riccatientire}(i)** | **GAP** | The uniqueness of the nonnegative entire Riccati solution for constant \(x=\kH\) is essential to the proof, but it is not among the supplied established results and no proof is given here. Without that result the convergence \((\eta_j,\zeta_j)\to(\mu,\mu)\) is not established. | Either cite the exact earlier proposition with its hypotheses and verify that “nonnegative entire” includes the limiting solution, or prove the constant-\(x\) uniqueness directly.

- **Lemma \(\ref{lem:rigidity}\)(ii), Step 2, “\(\mathbf N_j=\xi_j(0,\kH U_1)^{\mathsf T}+o(\xi_j)\)”** | **MINOR** | Part (i) defines \(o(\xi_j)\) only uniformly on compact time intervals, whereas the preceding estimates use suprema over all of \(\mathbb R\). The later convolution argument only needs the first-order forcing expansion on bounded intervals and the global \(O(\xi_j)\) bound for the tail, but the text does not distinguish these two modes of control. | State explicitly that the displayed first-order expansion is locally uniform, while \(\sup_{\mathbb R}|\mathbf N_j|=O(\xi_j)\) is used for the remote tail.

- **Lemma \(\ref{lem:rigidity}\)(iii), order-one averaging** | **MINOR** | The statement that the means of \(a',b'\) over \([0,T_j]\) tend to zero uses \(T_j\to T_0\), not exact \(T_j\)-periodicity of \(a,b\). This is correct, but the size needed is \(o(1)\), and it should be made explicit because \(T_j\ne T_0\). | Note that, for example,
  \[
  \frac1{T_j}\int_0^{T_j}a'(t)\,dt
  =\frac{a(T_j)-a(0)}{T_j}=o(1)
  \]
  by \(T_j\to T_0\) and \(T_0\)-periodicity.

- **Remark \(\ref{neg:rem:jet}\)(b), imported formula**
  \[
  \lambda_c=-\frac{A\mu^2}{2C_\mu D_\mu}\xi^2+o(\xi^2)
  \]
  | **GAP** | The displayed first Lyapunov coefficient in the supplied context does not by itself give this coefficient in the crest-amplitude variable \(\xi\). One also needs the normalization of the Hopf eigenvector, the relation between normal-form radius and \(\xi\), and the relation between the radial Floquet exponent and \(\lambda_c\). None is shown here. Calling the value “imported” does not identify a proved source for the conversion. | Supply the amplitude and time-normalization calculation, or cite a precise earlier result that already states this formula in the variable \(\xi\).

- **Remark \(\ref{neg:rem:selfcontained}\), final sentence** | **MINOR** | The argument only obtains remainders \(o(\xi_j)\), not the \(O(\xi^2)\) remainder stated in Theorem \(\ref{thm:failed}\) and Proposition \(\ref{neg:prop:charts}\)(ii). It does recover the sign, as claimed, but it does not independently recover the displayed asymptotic precision. | Explicitly say that the self-contained argument proves the leading coefficient and sign only, while the \(O(\xi^2)\) remainder uses smooth Hopf expansions.

- **Proposition \(\ref{neg:prop:planar}\)(iii), “This interval contains every value lying strictly between the minimum and maximum of \(w_2\) on \(\Dsc\)”** | **MINOR** | The image \(w_2(\inter\Dsc)\) is an interval because \(\inter\Dsc\) is connected and \(w_2\) is continuous, but the relation between its endpoints and extrema on the closed domain requires that interior points approximate boundary extrema. This is true for a Jordan domain but is not stated. | Replace the general assertion by the only fact needed: the common value on the small cycle lies in \(w_2(\inter\Dsc)\), since the cycle itself is contained in \(\inter\Dsc\).

- **Proposition \(\ref{neg:prop:planar}\)(iii), notation \(w_2(3T_0/8)\) on \(\Gamma_j\)** | **MINOR** | The cycle has period \(T_j\), not \(T_0\); evaluating at the fixed time \(3T_0/8\) is legitimate, but the notation may suggest a phase normalized to the cycle period. | Clarify that the cycles are parametrized in the original time with the crest at \(0\), and fixed times are allowed because convergence is local in time.

- **Table A2, “fails near the Hopf point: proved … with \(a\to0,b=1\)”** | **MINOR** | Proposition \(\ref{neg:prop:charts}\) is formulated in \((U,z)\), while row A2 refers to \((U,y-k)\). The two are related by the constant scaling \(y-k=k(z-1)\) at fixed parameter, so the conclusion is correct, but the coordinate conversion is omitted. | State that constant scaling of the second chart coordinate conjugates the return derivatives and preserves the componentwise radial inequality.

- **Table B8, status “refuted (fp)”** | **MINOR** | The reported finite-sample linear-program result does not refute the mathematical statement; the paper’s introductory convention acknowledges this weaker use of “refuted by floating-point computation,” but the table’s terse wording can still be mistaken for a proof. | Keep the explicit qualifier “on the sampled discretization” in the status cell.

- **Proposition \(\ref{neg:prop:signchange}\)(ii), case \(u_0=\xi_{\rm out}\)** | **MINOR** | The proof derives \(d>0\) on the whole interior while \(d\) is strictly increasing there and \(d(\varepsilon_{\max})=0\), which is impossible by continuity. Thus that endpoint case cannot actually occur, but the contradiction is not mentioned. The uniqueness conclusion remains valid because the hypothesis is then inconsistent. | Note explicitly that \(u_0=\xi_{\rm out}\) is impossible under \(d(\varepsilon_{\max})=0\).

- **Proposition \(\ref{neg:prop:signchange}\)(iii), definition of \(u_0\)** | **MINOR** | The proof does not separately discuss the possibility that the positivity set is all of \((0,\xi_{\rm out})\). It later excludes precisely that case when proving the endpoint value is negative, so the logic is recoverable but the definition “with some \(u_0\in[0,\xi_{\rm out})\)” is established only after the endpoint contradiction. | Reorder the proof: first exclude \(\bar\iota(\xi_{\rm out})=0\), then define the unique sign-change point.

- **After Proposition \(\ref{neg:prop:signchange}\), “the return time \(T_1(u)\) stays bounded and the arc of one turn shrinks to the origin”** | **GAP** | Bounded return time alone does not imply that the whole return arc shrinks to the equilibrium. This follows from local continuous dependence and the hyperbolic-focus linearization, but that argument is not supplied or cited. | Invoke a local focus/Poincaré-map theorem giving \(T_1(u)\to2\pi/\omega_E\) and uniform shrinking of the return arc, or prove it from the \(C^1\) linearization estimates.

- **Final subsection, “Every exact formula for the exponent obtained in Section … has the form \(\lambda_c=-\langle pN/S\rangle\) with \(N,S>0\)”** | **GAP** | This is a universal claim about every exact formula in another section, but only the two formulas involving \(K/\Sigma\) and \(W/\Xi\) are supplied in the context. The statement is stronger than what is documented here. | Enumerate the formulas meant and verify the claim for each, or replace “every” by “the gauge formulas used here.”

- **Final subsection, “A proof … needs an upper bound …”** | **MINOR** | Even with the later qualification that this is only a sufficient route, “needs” is too strong: alternative integral, topological, variational, or return-map arguments need not furnish such a pointwise upper bound. | Replace by “A proof following these gauge/Dulac formulas would need…”

### Verdicts

- Theorem \(\ref{thm:failed}\): **HOLDS** — the center-graph crest expansion, mismatch coefficient, limiting matrix entries, and first-component asymptotic were re-derived.
- Proposition \(\ref{neg:prop:charts}\): **HOLDS** — both component identities, the two leading coefficients, and \(c(a,b)\) were re-derived.
- Proposition \(\ref{prop:nodulac}\): **HOLDS** — the factorization at \(E\), \(Q_E'=2\alpha_EQ_E\), singularity estimate, and local model \(-c\log Q_E\) were re-derived.
- Proposition \(\ref{prop:nomonomial}\): **HOLDS** — all three defect means, the product inequality, monotonicity range of \(\chi\), and equality case were checked.
- Lemma \(\ref{lem:rigidity}\): **GAP** — the period-convergence step is abbreviated and the essential Riccati uniqueness theorem is imported without enough information here.
- Proposition \(\ref{neg:prop:planar}\): **HOLDS** — \(-pK\), the two gauge defects, the coordinate derivative \(1/(1+\eta)\), and the Hopf-phase comparison were re-derived.
- Proposition \(\ref{neg:prop:signchange}\): **HOLDS** — the energy-coordinate derivative, unimodal-displacement argument, and endpoint contradiction were checked.
