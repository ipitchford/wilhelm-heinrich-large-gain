- Lemma \(\ref{lg:lem:id}\)(a), “\(w_k<0\) by Lemma \(\ref{lg:lem:signs}(a)\)” | MINOR | This is a forward dependency. It is not circular because the proof of Lemma \(\ref{lg:lem:signs}\) uses only the cone exclusion, not Lemma \(\ref{lg:lem:id}\), but the dependency order is obscured. | Move Lemma \(\ref{lg:lem:signs}\)(a) before Lemma \(\ref{lg:lem:id}\), or insert its one-line cone argument here.

- Lemma \(\ref{lem:lgshortlobe}\), “Consequently \(T_->T_0(k)\), and \(T_0\) is increasing in \(k\)” | MINOR | The monotonicity is asserted without differentiation. | Since
  \[
  \omega_0(k)=\frac12\sqrt{3k^{2/3}-4k^{1/3}-1},
  \]
  its radicand has derivative
  \[
  2k^{-2/3}(k^{1/3}-1)>0\qquad(k>1),
  \]
  so \(\omega_0\) and hence \(2\pi-\pi/\omega_0\) are increasing for the range used, \(k\ge1000\).

- Lemma \(\ref{lg:lem:rise}\)(a), “for \(0<L\le p_k^2/(2k)\) the time \(t_L\in(s_0,T_-)\)” | MINOR | The displayed non-strict estimate \(|U_{\min}|\ge p_k^2/(2k)\) alone does not ensure that when \(L=p_k^2/(2k)\), the rising passage lies strictly after \(s_0\); equality could nominally put the level at the trough. | Record the strict estimate. On the rising phase,
  \[
  p'=-\mu p+|w|<k,
  \]
  because \(p>0\) away from \(s_0\) and \(|w|<k\). Thus \(|U_{\min}|>p_k^2/(2k)\), which gives \(t_L\in(s_0,T_-)\), including the endpoint value of \(L\).

- Lemma \(\ref{lg:lem:P1}\), “Also \(P_1<\lambda\)” | MINOR | This conclusion is not proved in the lemma’s proof. The argument only proves \(P_1\le p_k+k\bar s\). | Cite (F5): the whole cycle lies in the chemical domain, where \(p(t)<\lambda\); hence its maximum \(P_1<\lambda\).

- Lemma \(\ref{lg:lem:boot}\), “\(c_T=\frac{0.055}{0.3}=0.1833\ldots\)” versus Theorem \(\ref{thm:largegain}\), “\(c_T=0.1833\)” | MINOR | These are not literally the same constant:
  \[
  \frac{0.055}{0.3}=0.183333\ldots>0.1833.
  \]
  The proof establishes the theorem’s weaker decimal bound, but the notation identifies two unequal numbers. | Define \(c_T:=0.1833\), noting that the proof gives the stronger constant \(11/60=0.183333\ldots\), or use \(c_T=11/60\) throughout.

- Proof of Corollary \(\ref{cor:lgunique}\), second proof, “\(\Dom\setminus W^s(E)\) is the union of the basins; in particular there is at least one cycle” | MINOR | The conclusion also requires \(\Dom\setminus W^s(E)\neq\varnothing\), which is evident because a one-dimensional immersed stable manifold cannot fill the three-dimensional open domain, but it is not stated. | Add that \(W^s(E)\cap\Dom\) has empty interior, so the complement is nonempty.

- Proof of Corollary \(\ref{cor:lgunique}\), second proof, “cannot be a union of two or more disjoint nonempty open sets” | MINOR | Connectedness rules out a separation into two relatively open sets. For a possibly infinite family of basins, one should explicitly group one basin against the union of all the others. | If two basins exist, choose one \(B\); then \(B\) and the union of all other basins are disjoint, nonempty, relatively open sets whose union is \(\Dom\setminus W^s(E)\), contradicting connectedness.

- Proof of Corollary \(\ref{cor:lgreciprocal}\), “the scalar variational equations coincide in the time \(\tau\)” | MINOR | The claim is correct, but the coefficient \(x=ke^U\) also changes and the invariance is not displayed; merely saying the equations coincide is too compressed for the multiplier conclusion. | Write the normalized variational equation
  \[
  V_{\tau\tau\tau}+\sigma V_{\tau\tau}+V_\tau+\sigma\rho e^U V=0.
  \]
  Both \(\sigma\) and \(\rho\) are invariant under the reciprocal map, so the monodromy operators over the corresponding \(\tau\)-period are identical.

No defect was found in the following checked computations:

- In Lemma \(\ref{lg:lem:id}\)(a),
  \[
  \int_{T_-}^T w=P_d+p_k,\qquad
  M_+=w_d-w_k+P_d+p_k.
  \]
- In Lemma \(\ref{lg:lem:id}\)(b), multiplication by \(U''=p'\) gives
  \[
  A\int p'^2=\int xp^2=\mu A\int p^2,
  \]
  hence \(\int(x-k_H)p^2=0\).
- In Lemma \(\ref{lem:lgshortlobe}\),
  \[
  4\omega^2=3\gamma_E^2-2A\gamma_E-(1-\mu)^2,
  \]
  and the stated lower bound for \(\omega\) follows for \(k\ge1\).
- The formulas defining \(c_{\rm i}\), \(c_{\rm ii}\), \(c_{\rm b}\), the estimates in Phases (i) and (ii), and the \(\Phi_\mu\) bounds were re-integrated.
- In Lemma \(\ref{lg:lem:mass}\), multiplication leading to \(Q\) gives exactly
  \[
  Q=\frac14p_k^2\Bigl(\frac12p_k^2-2\mu P_1R\Bigr)-\frac1{\sqrt2}p_kk.
  \]
- In Lemma \(\ref{lg:lem:boot}\),
  \[
  f'(z)\ \text{has the sign of}\ 
  (\mu z+1)(\mu z^2-z-6.3),
  \]
  so the claimed sign reduction is correct.
- The numerical margins \(T_0(1000)\approx5.89277\), \(\vartheta_1(1000)\approx0.19402\), \(T_1(1000)\approx43.898\), \(z_{10}\approx10.169\), and the final \(q(1000)<0.0772\) are consistent with the displayed formulas.
- In Lemma \(\ref{lem:lglinear}\), the transformed compound system, the weighted \(\mathcal L\) identity, the Dini estimate for \(\mathcal G\), and the matching factor \(1+\alpha s+\beta s^2\) were re-derived.
- Under the reciprocal map,
  \[
  \lambda'=\frac{\lambda}{\mu},\qquad \tilde T=\mu T,\qquad
  \rho'=\rho,
  \]
  and applying the map twice returns \((\mu,k)\). The domain mismatch calculation is also correct.

Verdicts:

- Lemma \(\ref{lg:lem:id}\): **HOLDS**
- Lemma \(\ref{lem:lgshortlobe}\): **HOLDS**
- Lemma \(\ref{lg:lem:signs}\): **HOLDS**
- Lemma \(\ref{lg:lem:i}\): **HOLDS**
- Lemma \(\ref{lg:lem:ii}\): **HOLDS**
- Lemma \(\ref{lg:lem:rise}\): **HOLDS**
- Lemma \(\ref{lg:lem:stage1}\): **HOLDS**
- Lemma \(\ref{lg:lem:P1}\): **HOLDS**
- Lemma \(\ref{lg:lem:mass}\): **HOLDS**
- Lemma \(\ref{lg:lem:stage2}\): **HOLDS**
- Lemma \(\ref{lg:lem:boot}\): **HOLDS**
- Lemma \(\ref{lg:lem:crest}\): **HOLDS**
- Lemma \(\ref{lem:lglinear}\): **HOLDS**
- Theorem \(\ref{thm:largegain}\): **HOLDS**
- Corollary \(\ref{cor:lgunique}\): **HOLDS**
- Corollary \(\ref{cor:lgreciprocal}\): **HOLDS**
