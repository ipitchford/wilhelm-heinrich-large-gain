1. **MINOR — Lemma \(\ref{lem:id}\) (“Identities”).** No error found. The integrations give
   \[
   \int_{T_-}^T w=P_d+p_k,\qquad
   w_d-w_k=-\int_{T_-}^T w+M_+,
   \]
   hence \(M_+=w_d-w_k+P_d+p_k\). Part (b) also follows from
   \[
   A\int p'^2=\int xp^2=A\mu\int p^2=k_H\int p^2.
   \]

2. **MINOR — Lemma \(\ref{lem:short}\) (“The positive lobe is short”).** No error found. In particular,
   \[
   4\omega^2=3\gamma_E^2-2(1+\mu)\gamma_E-(1-\mu)^2,
   \]
   and substituting \(\gamma_E>k^{1/3}\), \(1+\mu\le2\), \((1-\mu)^2\le1\) gives the stated \(\omega_0\). The Sturm integration has the correct boundary sign.

3. **MINOR — Lemma \(\ref{lem:signs}\).** No error found. The zero-crossing uniqueness in part (b) follows because at any zero of \(w\) in the negative lobe,
   \[
   w'=-k\varphi<0.
   \]

4. **MINOR — Lemma \(\ref{lem:i}\) (“Phase (i)”).** No error found. The sign convention is correct:
   \[
   |p|'=w-\mu|p|
   \]
   while \(p<0,w>0\), and \(\mu\le1\) indeed gives
   \(e^{-\mu(s-r)}\ge e^{-(s-r)}\).

5. **MINOR — Lemma \(\ref{lem:ii}\) (“Phase (ii)”).** No error found. The backward terminal-value formula is
   \[
   |p|(s)=\int_s^{s_0}e^{\mu(r-s)}|w(r)|\,dr,
   \]
   not a forward stable convolution; the proof uses this correct formula. Also
   \[
   \int_0^u r(1-e^{-r})\,dr=\frac{u^2}{2}-1+(1+u)e^{-u}.
   \]

6. **MINOR — Lemma \(\ref{lem:rise}\) (“Rising phase”).** No error found. In part (c),
   \[
   p(s_\ell)=\int_{s_0}^{s_\ell}e^{-\mu(s_\ell-r)}|w(r)|\,dr
   \]
   gives the stated \(\Psi_\mu\), and
   \[
   p_k\ge p(s_\ell)e^{-\mu\tau_b}
       \ge p(s_\ell)-\mu\ell
   \]
   follows from
   \(\ell\ge p(s_\ell)(1-e^{-\mu\tau_b})/\mu\).

7. **MINOR — Lemma \(\ref{lem:stage1}\) (“Stage 1”).** No error found. The four-duration argument is valid: if simultaneously
   \(s_w<a_1,\ u<a_2,\ r_\ell<a_3,\ \tau_b<a_4\), then
   \[
   T_-<a_1+a_2+a_3+a_4=5.89,
   \]
   contradicting \(T_->5.8927\ldots\). The quoted lower coefficient
   \[
   \varphi_\ell\Psi_1(1.48)-2/1000\approx0.3745
   \]
   is consistent.

8. **MINOR — Lemma \(\ref{lem:P1}\).** No error found. The apparent issue that \(\bar s\) might not lie in \([0,s_+]\) is harmless: the argument is by contradiction under \(s_+\ge\bar s\); otherwise \(s_+<\bar s\) already. The bound \(P_1<\lambda\) is simply (F5).

9. **MINOR — Lemma \(\ref{lem:mass}\) (“Mass from slope”).** No error found. The change of variables gives
   \[
   m(s)=\int_0^{U(s)}\frac{k(e^V-1)}{p}\,dV
   \le \frac{k(e^{U(s)}-1-U(s))}{\theta p_k}.
   \]
   Expanding \(Q\) yields
   \[
   Q=\frac{p_k^4}{8}\left(
   1-\frac{4\mu P_1H}{p_k^2}-\frac{4\sqrt2\,k}{p_k^3}
   \right),
   \]
   so the displayed \(\delta_Q\) is correct after \(P_1\le(1+\varepsilon_1)p_k\).

10. **MINOR — Lemma \(\ref{lem:stage2}\) (“Stage 2”).** No error found. The lower bound
    \[
    r_\ell\ge0.4708T_--0.9416\log T_--0.6987
    \]
    exceeds \(0.3T_-\) for \(T_-\ge22\), and
    \[
    \Psi_\mu(s)\ge0.95\frac{1-e^{-\mu(s-3)}}{\mu}
    \]
    follows by restricting its defining integral to \(r\in[3,s]\).

11. **GAP — Lemma \(\ref{lem:boot}\), Step 3, interval placement.** The proof invokes endpoint quasi-convexity on
    \([z_1(k),z^\sharp]\) without explicitly proving that this interval is nonempty, i.e.
    \(z_1(k)\le z^\sharp\), uniformly in \(\mu\). The missing computation is easy but required:
    \[
    z_1(k)
    =\frac{0.3c_1^2\sqrt{1-\delta_1(K_0)}}{\sqrt8}\,k-3.15
    <0.0149k-3.15,
    \]
    whereas
    \[
    z^\sharp=\frac{0.055k}{\mu^2}\ge0.055k.
    \]
    Thus \(z_1(k)<z^\sharp\), but this is not stated where the interval argument depends on it.

12. **MINOR — Lemma \(\ref{lem:boot}\), quasi-convexity itself.** No error found once the preceding interval placement is supplied. Direct differentiation gives
    \[
    \operatorname{sgn}h'(z)
      =\operatorname{sgn}(\mu z^2-z-6.3),
    \]
    so \(h\) decreases and then increases, implying its maximum on a compact interval occurs at an endpoint. The implication
    \(h(z)>0.3\kappa_0k\Rightarrow z>z^\sharp\) is then valid.

13. **GAP — Lemma \(\ref{lem:crest}\), opening claim about numerical side conditions.** The sentence
    > “All majorants below have the form \((c+d\log\lambda)\lambda^{-a}\) or are monotone for obvious reasons, so each numerical side condition needs to be checked only at \(\lambda=K_0\)”
    
    is not a proof and is literally inaccurate. Quantities such as
    \[
    \delta'
      =\frac{4(1+\varepsilon_1)\tilde H_\lambda}{p_k}
       +\frac{4k}{\theta p_k^3},
    \quad
    \hat W=\frac{c_p^2\lambda^2e^{-\bar\tau_1}}{A_1},
    \quad
    v_f=b\sqrt{\frac{2(\hat H+\hat L)}{\hat W}}
    \]
    are nested products and quotients, not of the asserted simple form. To reduce checks such as
    \[
    b^2\!\left[1-\frac23\left(1+\frac{\lambda}{\hat W}\right)v_f\right]\ge1,
    \qquad
    v_f\le\log(1+\hat W/\lambda)
    \]
    to \(\lambda=1000\), one must prove that \((1+\lambda/\hat W)v_f\) decreases and that
    \(\log(1+\hat W/\lambda)-v_f\) increases. No such derivative or comparison is provided.

14. **GAP — Lemma \(\ref{lem:crest}\), post-crest convolution.** The displayed estimate for \(|p|\) uses the kernel \(e^{-(v-r)}\), whereas the exact equation is
    \[
    |p|'=w-\mu|p|,\qquad
    |p|(v)=\int_0^v e^{-\mu(v-r)}w(r)\,dr.
    \]
    The stated lower bound is valid only after explicitly using
    \[
    e^{-\mu(v-r)}\ge e^{-(v-r)}
    \]
    and only on the interval where the lower majorant for \(w\) is nonnegative. The proof has arranged \(v\le v_+\), so this can be repaired, but the necessary sign/kernel argument is omitted.

15. **MINOR — Lemma \(\ref{lem:crest}\), remaining algebra.** No further algebraic error found. In particular,
    \[
    X_\theta\ge\frac{\theta p_k^3}{4}(1-\delta')
    \]
    correctly includes both losses \(4\mu P_1H/p_k^2\) and \(4k/(\theta p_k^3)\), and the quadratic-root bound for \(v_c\) is consistent.

16. **MINOR — Lemma \(\ref{lem:lyap}\).** No error found. The transformed system is
    \[
    a'=-\mu a+b,\qquad b'=-b+c,\qquad c'=-Ac+xa,
    \]
    with multiplier \(m_t\). For the max norm,
    \[
    D^+\mathcal N
    \le\max\{s-\mu,\ s-1,\ x/s^2-A\}\mathcal N
    \le(s-\mu+x/s^2)\mathcal N,
    \]
    using \(\mu\le1\). The matching factor is exactly
    \[
    s^2\bigl(s^{-2}+\alpha s^{-1}+\beta\bigr)
    =1+\alpha s+\beta s^2.
    \]

17. **MINOR — Proof of Theorem \(\ref{thm:main}\).** No error found conditional on the gaps in Lemma \(\ref{lem:crest}\). The numerical estimate
    \[
    q(1000)\le\frac{5.544+8.606}{183.33}<0.0772
    \]
    gives
    \(\avg\eta\le0.8272\mu<0.83\mu\).

18. **GAP — Corollary \(\ref{cor:unique}\), topology of \(W^s(E)\).** The proof introduces facts not contained in (F1)–(F7) and not derived from the equations:
    > “\(W^s(E)\cap\mathcal D\) is a countable union of compact \(C^1\) arcs, a set of topological dimension one, which cannot disconnect a connected \(3\)-manifold.”
    
    At minimum this requires the stable-manifold theorem, a proof that the global stable manifold is an injectively immersed one-manifold in \(\mathcal D\), the countable-closed-sum theorem for covering dimension, and a separation theorem for subsets of dimension at most one in a \(3\)-manifold. None is imported or proved. Under the document’s rule that everything besides (F1)–(F7) must follow from those facts and the equations, the connectedness of
    \(\mathcal D\setminus W^s(E)\) is unjustified.

19. **FATAL — Corollary \(\ref{cor:recip}\), claimed conjugacy of the positive domains.** The assertion
    > “the same change of variables conjugates the two flows on their positive domains”
    
    is false. Under
    \[
    \mu'=\frac1\mu,\qquad k'=\frac{k}{\mu^3},\qquad
    \tilde U(\tilde t)=U(\tilde t/\mu),
    \]
    one has
    \[
    \tilde p=\frac p\mu,\qquad
    \tilde y
      =k'-\mu'\tilde p-\tilde p'
      =\frac{k-\mu p-\mu p'}{\mu^3}.
    \]
    This is not a positive multiple of
    \[
    y=k-\mu p-p'.
    \]
    Concrete counter-computation: take \(\mu=2\), \(k=8000\), \(p=0\), and choose an original-domain state with \(y=1\). Then
    \[
    p'=k-y=7999,
    \]
    while the transformed coordinate is
    \[
    \tilde y=\frac{8000-2(7999)}{8}=-999.75<0.
    \]
    Thus a point of the original positive domain maps outside the reciprocal positive domain. Consequently the claimed correspondence of basins and stable manifolds does not follow from this transformation.

20. **GAP — Corollary \(\ref{cor:recip}\), possible repair not given.** The false positive-domain conjugacy is unnecessary, but the proof does not supply the correct replacement. One can transfer cycles and their multipliers through the scalar equation, obtain uniqueness of cycles from uniqueness at \((1/\mu,k/\mu^3)\), and then apply (F7) directly at the original parameters to get convergence in the original \(\mathcal D\). That argument would still inherit the topology gap in Corollary \(\ref{cor:unique}\), but it avoids the false domain-conjugacy claim.
