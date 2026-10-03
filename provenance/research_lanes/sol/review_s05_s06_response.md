### Defects

- **Proposition \(\ref{prop:polar}\), sentence “The integrand is positive for \(\theta\in(0,\pi)\)”** | **MINOR** | The displayed formula includes an overall minus sign:
  \[
  \lambda_cT=-\int_0^{2\pi}\frac{K\sin\theta}{r\sigma^2}\,d\theta.
  \]
  Thus the fraction \(K\sin\theta/(r\sigma^2)\) is positive on \((0,\pi)\), but the signed integrand of the displayed integral is negative there. | Say “the fraction under the integral sign” or absorb the minus sign into the integrand and reverse the stated signs.

- **Proposition \(\ref{prop:holder}\), phrase “\(\avg{e^{U/3}}<\rho^{-1/3}\)”** | **MINOR** | The supplied context and local notation do not define \(\rho\). The proof implicitly uses
  \[
  \rho=\frac{k}{k_H},
  \]
  since \(k\avg{e^{U/3}}^3<k_H\). | Define \(\rho:=k/k_H\) before the proposition or replace \(\rho^{-1/3}\) by \((k_H/k)^{1/3}\).

- **Proposition \(\ref{prop:lvchain}\)(ii), phrase “These are the second and third rows of \(H'=BH\)”** | **MINOR** | For arbitrary positive \(\varphi\), the constructed \(H\) does not satisfy the full equation \(H'=BH\): the first row would require \(H_1'=xH_3\), which is precisely not assumed. Only the identities \(H_3'=H_2+H_3\) and \(H_2'=H_1+\mu H_2\) hold by construction. The subsequent computation itself does not use the false first row. | Replace by “The definitions enforce the second and third row identities of \(H'=BH\); the first-row defect is \(H_1'/H_3\).”

- **Lemma \(\ref{wk:lem:line}\), use of “inequality \(\eqref{eq:equilibriumreverse}\)”** | **MINOR** | That imported label and statement are not among the supplied context. The needed inequality is nevertheless directly derivable from \(\Lambda_E>0\):
  \[
  (y-k)-\zeta_EU+\eta_Ep>0
  \]
  and \(y-k=-p'-\mu p\), giving
  \[
  -p'+(\eta_E-\mu)p-\zeta_EU>0.
  \]
  | Cite Corollary \(\ref{cor:cones}\) directly and include this one-line substitution.

- **Lemma \(\ref{lem:strongstable}\)(d), sentence “By the Sturm separation theorem every nontrivial real solution of this equation vanishes somewhere”** | **MINOR** | As written, this overstates Sturm separation: the conclusion uses the additional fact that \(p\) has consecutive zeros on the whole real line, obtained by periodic repetition of its two simple zeros. Without explicitly invoking those repeated consecutive zeros, the quoted general statement is false for an arbitrary second-order equation with a solution having only an isolated zero pair. | State: “Since the periodic function \(p\) has an infinite sequence of consecutive simple zeros, any linearly independent solution has a zero between each consecutive pair; a dependent solution vanishes at the zeros of \(p\).”

- **Proposition \(\ref{prop:crestmonotone}\), Step 3, phrase “The same dichotomy applies on \((t_a,\infty)\)”** | **GAP** | The second application is compressed too far. For \(p>0\), one must separately establish: if \(U\to+\infty\), then the monotonicity of \(U\) implies \(x=ke^U\to\infty\), the convolution formula for \(y'+y=x\) implies \(y\to\infty\), and then \(p'\le k-y\) forces \(p\) to cross zero in finite time. If \(U\) is bounded, precompactness and the invariance argument must be repeated. These facts are plausible and enough, but are not actually supplied by the phrase “same dichotomy.” | Add those arguments explicitly. This is a local repair.

- **Proposition \(\ref{prop:crestmonotone}\), Step 4, claim that nearby solutions have “the same pattern of zeros”** | **GAP** | Transversality at \(t_a,t_b\) and negativity on a short initial interval do not alone exclude additional zeros for nearby solutions elsewhere in \((0,t_b)\). One also needs compact sign bounds for \(p\) away from small neighborhoods of \(t_a,t_b\). | Choose disjoint neighborhoods of \(t_a,t_b\); on the compact complement, \(p\) is bounded away from zero. Continuous dependence then excludes extra zeros, and the implicit-function theorem tracks the two transversal zeros.

- **Proposition \(\ref{prop:gaussmap}\), final orientation argument** | **MINOR** | The phrase “has its interior on the right of the velocity” is asserted rather than justified. At a top point of a regular Jordan curve whose tangent points right, the local interior is indeed below the tangent and hence to the right, but this uses the Jordan separation/local-side property. | Add one sentence invoking the local Jordan-curve side theorem, or define orientation via the signed area/tangent convention.

### Checked identities that do not produce defects

I re-derived the Wronskian equation, all four values of \(L_2\), both gauge ODEs, the state formula for \(\Xi\), the formulas for \(\eta_0,\eta_1\), the \(\kappa\)-Riccati equation, the compound-source equation, the polar system and exponent substitution, the nested-radius derivative, the quadratic form in Proposition \(\ref{prop:monotoneintegrals}\), the Hölder/Jensen inequalities, the Lotka–Volterra system, the Cole–Hopf identity, the strong-pairing formulas, the scalar equations for \(\mathfrak r\), the twelve-event zero-crossing derivatives, and the rank-one boundary term. No sign or algebra error was found in those computations beyond the wording issues above.

### Verdicts

- Lemma \(\ref{gau:wronskian}\): **HOLDS**
- Theorem \(\ref{thm:gauges}\): **HOLDS**
- Corollary \(\ref{gau:halfcycle}\): **HOLDS**
- Proposition \(\ref{prop:kappa}\): **HOLDS**
- Proposition \(\ref{prop:polar}\): **HOLDS**
- Proposition \(\ref{prop:nestedradii}\): **HOLDS**
- Proposition \(\ref{prop:monotoneintegrals}\): **HOLDS**
- Corollary \(\ref{gau:convexcriterion}\): **HOLDS**
- Proposition \(\ref{prop:holder}\): **HOLDS**
- Proposition \(\ref{prop:lvchain}\): **HOLDS**
- Proposition \(\ref{prop:strongidentity}\): **HOLDS**
- Proposition \(\ref{prop:gaussmap}\): **HOLDS**
- Lemma \(\ref{lem:strongstable}\): **HOLDS**
- Corollary \(\ref{wk:cor:reform}\): **HOLDS**
- Lemma \(\ref{wk:lem:weakadjoint}\): **HOLDS**
- Theorem \(\ref{thm:twelve}\): **HOLDS**
- Lemma \(\ref{wk:lem:line}\): **HOLDS**
- Theorem \(\ref{thm:dualsupport}\): **HOLDS**
- Proposition \(\ref{prop:crestmonotone}\): **GAP**
- Proposition \(\ref{prop:rankone}\): **HOLDS**
