## (a) Verdict: PARTIAL

The listed structure yields several rigorous partial conclusions, but I do not see a valid deduction of  
\[
\langle\eta\rangle<\mu
\]
for every periodic orbit. In particular, positivity of \(K\), convexity of \(h\), and the Riccati equation for \(\kappa\) do not by themselves determine the sign of the signed line integral
\[
\oint_\Gamma \kappa\,dU.
\]

What can be proved is:

1. Every periodic orbit which is innermost around \(E\) satisfies
   \[
   \lambda_c=\langle\eta\rangle-\mu\le 0.
   \]
2. The strict inequality holds for the regimes already listed in (10): sufficiently small Hopf cycles, and \(k\ge1000,\ \mu\le1\).
3. A precise global monotonicity inequality for the one-turn divergence integral would close the proof, imply strict negativity, and imply uniqueness. It is stated below.

---

## (b) Rigorous partial argument

### 1. Reduction to the planar flow

On the invariant graph \(y=h(U,p)\), the dynamics is
\[
U'=p,\qquad p'=F(U,p):=k-\mu p-h(U,p).
\]
Since
\[
h_U=\zeta,\qquad h_p=-\eta,
\]
the planar divergence is
\[
\operatorname{div}(p,F)=F_p=-\mu-h_p=\eta-\mu.
\]
Consequently, for a periodic orbit \(\Gamma\),
\[
\int_0^T\operatorname{div}(p,F)\,dt
=(\langle\eta\rangle-\mu)T
=\lambda_cT.
\]
Thus the nontrivial planar Floquet multiplier is
\[
m_t=e^{\lambda_cT}.
\]

The desired statement is exactly that every planar periodic orbit is attracting:
\[
m_t<1.
\]

### 2. An innermost periodic orbit cannot have \(\lambda_c>0\)

Let \(\Gamma\) be a periodic orbit whose interior contains no other periodic orbit. Its interior contains the unique equilibrium \(E\), which is an unstable focus for \(k>k_H\).

Choose a point \(q\) in the interior of \(\Gamma\), sufficiently close to \(E\), with \(q\ne E\). Because \(E\) is a hyperbolic unstable focus in the planar flow, \(q\) leaves a sufficiently small neighborhood of \(E\) in forward time. Its forward orbit remains in the compact closed region bounded by \(\Gamma\).

Its \(\omega\)-limit set cannot be \(E\): a hyperbolic source has no nontrivial forward orbit converging to it. By the planar Poincaré–Bendixson theorem, the \(\omega\)-limit set therefore contains a periodic orbit. By the assumed innermost property of \(\Gamma\), that periodic orbit must be \(\Gamma\). Hence some points from the interior approach \(\Gamma\) in forward time.

If \(\lambda_c>0\), then
\[
m_t=e^{\lambda_cT}>1,
\]
so the local Poincaré return map is repelling on both sides of \(\Gamma\). No point sufficiently close to \(\Gamma\), other than a point on \(\Gamma\), can then converge to \(\Gamma\) in forward time. This contradicts the conclusion above. Therefore
\[
\boxed{\lambda_c\le0}
\]
for every innermost periodic orbit.

This argument does not exclude the nonhyperbolic possibility \(\lambda_c=0\): a periodic orbit with multiplier \(1\) can be semistable and attract from one side.

### 3. Why positivity of \(K\) does not directly settle the sign

The exact formula is
\[
\lambda_cT=-\oint_\Gamma \kappa\,dU,
\qquad
\kappa=\frac K\Sigma>0.
\]
Because \(dU=p\,dt\) changes sign between the rising and falling arcs, the pointwise inequality \(\kappa>0\) gives no sign for the integral.

A natural stronger conjecture would be
\[
\kappa_{\rm up}(U)>\kappa_{\rm down}(U)
\quad\text{for every }U,
\]
but this is not plausible and is already incompatible with the leading Hopf expansion.

Indeed, at the Hopf point \(k=k_H=\mu(1+\mu)\), write \(\omega=\sqrt\mu\) and use the leading harmonic orbit
\[
U=a\cos\phi,\qquad p=-a\omega\sin\phi,\qquad \phi=\omega t.
\]
Then
\[
\Sigma=p^2-U p'=\mu a^2+O(a^3).
\]
Since
\[
K'+(1+\eta)K=k\psi(U),
\qquad
\psi(U)=\frac{U^2}{2}+O(U^3),
\]
and \(\eta=\mu+O(a)\), the quadratic part \(K_2\) satisfies
\[
\omega\frac{dK_2}{d\phi}+(1+\mu)K_2
=\frac{k_Ha^2}{4}(1+\cos2\phi).
\]
Thus
\[
K_2=a^2\bigl(C_0+\alpha\cos2\phi+\beta\sin2\phi\bigr),
\qquad \beta>0.
\]
At the two points with the same \(U=a\cos\phi\), the \(\sin 2\phi\) term changes sign. Hence
\[
K_{\rm up}-K_{\rm down}
\]
has, to leading order, opposite signs for \(U<0\) and \(U>0\). The same is true for \(\kappa_{\rm up}-\kappa_{\rm down}\), because \(\Sigma=\mu a^2+O(a^3)\).

Therefore a pointwise branch ordering is false even infinitesimally near Hopf. The desired sign comes from a higher-order integrated imbalance, consistent with
\[
\lambda_c\sim-c(\mu)a^2.
\]

### 4. Established strict cases

From the supplied results, one rigorously has
\[
\lambda_c<0
\]
in the following cases:

- sufficiently small periodic orbits on the supercritical Hopf branch, because
  \[
  \lambda_c\sim-c(\mu)\,a^2,\qquad c(\mu)>0;
  \]
- every periodic orbit when
  \[
  \mu\le1,\qquad k\ge1000.
  \]

These do not exclude additional, larger cycles for parameters merely close to the Hopf curve unless an independent global continuation or uniqueness argument is supplied.

---

## (c) The most promising missing inequality

Use the crest section
\[
\mathcal S=\{(U,p):p=0,\ U>0,\ F(U,0)<0\}.
\]
Let \(P(a)\) be the next-crest return map starting from \((a,0)\), and let \(\tau(a)\) be the corresponding return time. Define the one-turn divergence integral
\[
I(a):=\int_0^{\tau(a)}
\bigl(\eta(U(t),p(t))-\mu\bigr)\,dt.
\]

At a fixed point \(P(a)=a\), corresponding to a periodic orbit,
\[
I(a)=\lambda_cT.
\]

Set
\[
V(a):=-F(a,0)=h(a,0)-k>0
\]
on the crest section, and introduce the normalized section coordinate
\[
q=q(a),\qquad \frac{dq}{da}=\frac1{V(a)}.
\]
Let \(\widehat P\) denote the return map in the \(q\)-coordinate.

The standard planar Poincaré-map derivative formula gives
\[
P'(a)
=\frac{F(P(a),0)}{F(a,0)}e^{I(a)}
=\frac{V(P(a))}{V(a)}e^{I(a)}.
\]
Therefore
\[
\widehat P'(q(a))
=
\frac{q'(P(a))P'(a)}{q'(a)}
=e^{I(a)}.
\]

Hence the following single inequality would suffice:

\[
\boxed{\frac{d}{dq}I(a(q))<0}
\tag{M}
\]

throughout the domain of the crest return map.

Equivalently,
\[
\boxed{\widehat P''(q)<0.}
\]

### Why (M) closes the proof

Near the unstable focus, the return map moves points outward, so
\[
\widehat P(q)-q>0
\]
at the small-amplitude end of the section. Under (M), \(\widehat P'\) is strictly decreasing, so
\[
D(q):=\widehat P(q)-q
\]
is strictly concave. Starting positive near the source end, such a displacement can have at most one zero. At any zero \(q_*\), strict concavity forces
\[
D'(q_*)<0,
\]
and therefore
\[
\widehat P'(q_*)<1.
\]
Using \(\widehat P'(q_*)=e^{I(a_*)}\), this gives
\[
I(a_*)<0,
\]
hence
\[
\lambda_c<0.
\]
It also gives uniqueness of the periodic orbit.

### Plausibility assessment

I regard (M) as plausible but genuinely unproved from the listed identities:

- It is exactly consistent with the reported numerical observation that the one-turn mean of \(\eta-\mu\) decreases with crest amplitude.
- The negative first Lyapunov coefficient verifies its local consequence near Hopf.
- The normalization \(dq/da=1/V(a)\) removes the endpoint velocity factor from the Poincaré derivative, making \(I\) the logarithmic derivative of the normalized return map.

On the other hand, neither convexity of \(h\) nor positivity of \(K\) appears to imply (M) pointwise. The Hopf calculation above shows that simpler branchwise comparisons fail at leading order. Thus proving (M), or an integrated variant strong enough to imply strict concavity of \(\widehat P\), is the substantive missing step.
