### Defects

1. **Imported result (I2)(b), “Successive intersections of an orbit with a transversal segment are monotone along the segment, and a periodic orbit meets a transversal segment in at most one point.”** | **FATAL as stated** | This is false for an arbitrary transversal segment. For the circular flow \(q'=(-y,x)\), the vertical segment \(\{x=0,-2<y<2\}\) is transverse at both intersections with the unit circle, and the same periodic orbit meets it twice. Successive intersections with such a segment need not be monotone. Several later arguments cite this overbroad version. | Replace it by the corresponding statement for a sufficiently small local transversal contained in a flow box, or formulate an order-preserving local Poincaré map. In the half-axis arguments, prove monotonicity directly from the positive derivative of the return map and the absence of intervening fixed points.

2. **Proof of (A1)(a), “since \(t_{\mathbf z_i}\ge1\), induction gives \(\varphi_t(\mathbf z)\in\mathcal B_1\) for all \(t\ge0\). Hence \(\mathcal B_1\) is forward invariant.”** | **GAP** | The induction described directly treats initial points in \(\mathcal B'\), not arbitrary initial points in \(\mathcal B_1=\varphi([0,t_*]\times\mathcal B')\). The final sentence does not follow literally from the preceding induction. | First prove that every orbit starting in \(\mathcal B'\) remains in \(\mathcal B_1\): at each return to \(\mathcal B'\), select a new covering neighbourhood and repeat. Then, if \(z=\varphi_r(b)\in\mathcal B_1\), \(b\in\mathcal B'\), \(0\le r\le t_*\), one has \(\varphi_t(z)=\varphi_{t+r}(b)\in\mathcal B_1\). This repairs the argument locally.

3. **Proof of (A1)(a), “Any compact \(Q\subset\Dom\) is covered by finitely many \(V_{\mathbf z}\), so \(\varphi_t(Q)\subset\mathcal B_1\) for all \(t\ge t_Q\).”** | **GAP** | A finite cover alone does not give this conclusion unless one also uses the preceding “once a point reaches \(\mathcal B'\), its future remains in \(\mathcal B_1\)” argument. The entry times differ among the covering sets. | Let \(t_Q\) be the maximum of the finitely many selected entry times. Each point has entered \(\mathcal B'\) by its own entry time and thereafter remains in \(\mathcal B_1\); then the claimed uniform conclusion follows.

4. **Proof of (A1)(b), converse invariance, “by injectivity of \(\varphi_s\) the point \(\mathbf w=\varphi_{j-s}(\mathbf b_j)\) does not depend on \(j\), lies in \(\varphi_{j-s}(\mathcal B_1)\) for all \(j\ge s\), hence in \(\cA\).”** | **MINOR** | Membership at the discrete times \(j-s\to\infty\) does not immediately equal membership in the intersection over every real \(t\ge0\); the needed use of nestedness is omitted. | For arbitrary \(t\ge0\), choose \(j\) with \(j-s\ge t\). Since the sets are decreasing, \(\varphi_{j-s}(\mathcal B_1)\subset\varphi_t(\mathcal B_1)\), hence \(\mathbf w\in\varphi_t(\mathcal B_1)\).

5. **Proof of (A6), “Take a transversal segment \(I\) through \(q_*\in\mathcal J\). … By (I2b) these points … form a strictly monotone bi-infinite sequence.”** | **GAP** | This invokes the false global version of (I2)(b). An arbitrary transversal segment can meet the periodic orbit more than once and does not carry the asserted global ordering. | Choose \(I\) inside a sufficiently small flow box around \(q_*\), so that \(\mathcal J\cap I=\{q_*\}\), and use the local order-preserving return correspondence. With that repair, a strictly monotone sequence confined to one side of \(q_*\) cannot tend to \(q_*\) in both time directions.

6. **Lemma \(\ref{att:lem:cyclegeom}\), “The winding number … is the signed number of its crossings of the ray.”** | **MINOR** | The proof omits the finiteness and regularity justification needed for this crossing-count formula. Transversality makes each crossing isolated, but the argument should state that compactness then yields finitely many crossings. | Note that every crossing is transverse by Lemma \(\ref{att:lem:axis}\); the crossing set is therefore discrete and compact, hence finite. The standard ray-crossing formula then applies.

7. **Proof of (A7), winding-number argument, “The latter loop extends to the map … into \(\mathbb R^2\setminus\{q_0\}\) and so has winding number \(0\).”** | **MINOR** | This uses the additional standard fact that a loop extending continuously over a disc in the punctured plane is null-homotopic and has zero winding number. It is not explicitly among the imported statements, though it is elementary. | Add this fact to (I4), or give the one-line null-homotopy argument.

8. **Proof of (A8)(d), return map near a crest, “The point \(\Pi(\mathbf a_j)\) lies on \(\Gamma_j\) and in \(\{p=0,U>0\}\). … Therefore \(\Pi(\mathbf a_j)=\mathbf a_j\).”** | **MINOR** | The local return map is constructed using a return time near \(T\), but the proof does not explicitly ensure that \(\Pi(\mathbf a_j)\) remains in the neighbourhood where \(U>0\). The statement follows only after shrinking both the domain and image of the return map. | Shrink \(\mathcal T'\) so that \(\Pi(\mathcal T')\subset\{U>0\}\cap\mathcal T\). Then the unique positive-axis intersection property forces \(\Pi(\mathbf a_j)=\mathbf a_j\).

9. **Proposition \(\ref{prop:crestaxis}\)(iii), “\(\iota\) is strictly decreasing … equivalently … \(R_*\) is strictly concave.”** | **GAP** | The displayed derivative gives
   \[
   R_*'(\varepsilon)=e^{\iota(G_*^{-1}(\varepsilon))}.
   \]
   Strict decrease of \(\iota\) indeed implies strict concavity. The converse needs the additional theorem that a differentiable strictly concave function has strictly decreasing derivative. This is true only after supplying the standard argument: concavity gives nonincrease, and equality of derivatives at two points forces affinity on the intervening interval, contradicting strict concavity. It is asserted rather than proved. | Add this argument explicitly.

10. **Proposition \(\ref{prop:crestaxis}\)(iii), “It is strictly decreasing, since equality at two points would make \(R_*\) linear between them.”** | **GAP** | The preceding argument only establishes that \(R_*(\varepsilon)/\varepsilon\) is nonincreasing. Equality of this ratio at two points does not, by itself and without tracking the omitted \(R_*(\delta)\)-term or an extension to \(\varepsilon=0\), immediately imply linearity. | Use strict concavity directly. For \(0<\delta<\varepsilon_1<\varepsilon_2\),
   \[
   R_*(\varepsilon_1)>
   \frac{\varepsilon_2-\varepsilon_1}{\varepsilon_2-\delta}R_*(\delta)
   +\frac{\varepsilon_1-\delta}{\varepsilon_2-\delta}R_*(\varepsilon_2),
   \]
   unless the function is affine on \([\delta,\varepsilon_2]\). Since \(R_*(\delta)>0\), passing to \(\delta\downarrow0\) yields the required strict ratio inequality, or alternatively first establish a concave extension at \(0\).

11. **Corollary \(\ref{att:cor:dulac}\)(ii), “By (I2b) these points … are strictly monotone.”** | **GAP** | Again this cites the false global transversal assertion for the whole positive half-axis. | Use Proposition \(\ref{prop:crestaxis}\)(ii), which gives \(R_1'>0\). Since there is no fixed point in \((0,\xi_{\rm in})\), continuity implies \(R_1(u)-u\) has constant sign there; convergence of the orbit to \(\mathcal J_{\rm in}\) forces the sign to be positive. Thus \(R_1^j(u)\uparrow\xi_{\rm in}\).

12. **Corollary \(\ref{att:cor:dulac}\)(ii), citation “By Proposition \(\ref{prop:crestaxis}\)(i), \(R_1\) and \(T_1\) are continuous at \(\xi_{\rm in}\).”** | **MINOR** | The statement of Proposition \(\ref{prop:crestaxis}\)(i) explicitly asserts endpoint extension only at \(\xi_{\rm out}\), not continuity at an arbitrary interior fixed point such as \(\xi_{\rm in}\). Its proof does establish local \(C^1\) regularity at every \(u_0<\xi_{\rm out}\), so the needed fact is available, but the citation is imprecise. | Cite the local \(C^1\) assertion on \((0,\xi_{\rm out})\), or state separately that \(\xi_{\rm in}<\xi_{\rm out}\) unless the cycle is unique, in which case the endpoint clause applies.

13. **Remark \(\ref{att:rem:planes}\), “the solutions are unique because \(H\) has slope at least \(\zeta_->0\) in \(u\) and at most \(-\eta_-<0\) in \(v\).”** | **MINOR** | Uniqueness follows from strong monotonicity, but global existence of a solution for every right-hand side also requires the limits to be \(\pm\infty\). Those limits are not stated. | Derive, for \(u_2>u_1\),
   \[
   H(u_2,v)-H(u_1,v)\ge \zeta_-(u_2-u_1),
   \]
   and analogously in \(v\). These inequalities imply the required coercive limits and hence global solvability.

14. **Opening claim, “Their proofs import one result, the Birkhoff--Hopf contraction theorem.”** | **MINOR** | The proofs of (A6)–(A8) explicitly import Poincaré–Bendixson theory, invariant-manifold theory, orbital stability, Jordan/Schoenflies topology, winding-number facts, and Whitney extension. The sentence is only correct if “their” refers exclusively to (A1)–(A5), but its antecedent is ambiguous after the preceding paragraph discusses the whole theorem. | Replace by “The proofs of (A1)–(A5) import only the Birkhoff–Hopf theorem; (A6)–(A8) additionally use (I2)–(I5).”

### Verdicts

- **Lemma \(\ref{lem:absorbing}\) (Explicit absorbing box): HOLDS** — the variation-of-constants estimates, logistic upper bound, and threshold-interval lower bound check out.
- **Theorem \(\ref{thm:attractor}\) (Attractor graph): GAP** — the principal statements appear repairable, but the proof contains the false imported transversal theorem and several consequent gaps.
- **Lemma \(\ref{att:lem:cone}\) (Uniform positivity and contraction): HOLDS** — the Peano–Baker comparison, irreducibility, projective diameter bound, and factor count were re-derived.
- **Lemma \(\ref{att:lem:sup}\) (Suprema of bounded \(C^1\) functions): HOLDS** — the penalisation calculation and derivative estimate are correct.
- **Proposition \(\ref{prop:riccatientire}\) (Unique positive entire Riccati solution): HOLDS** — nested projective images, uniqueness, explicit bounds, and logarithmic contraction were checked.
- **Theorem \(\ref{thm:twoorbit}\) (Two-orbit comparison): HOLDS** — the nonlinear remainder identity, adjoint cancellation, integral formula, and strictness argument are correct.
- **Lemma \(\ref{att:lem:gauges}\) (Gauges on the attractor): HOLDS** — the identities for \(O,K,W,\Lambda_E\), including the cancellation in the \(W'\) computation, were re-derived.
- **Lemma \(\ref{att:lem:axis}\) (Reduced flow on the \(u\)-axis): HOLDS** — the axis identity and equilibrium conclusion follow correctly from the PDE.
- **Proposition \(\ref{att:prop:scaling}\) (Radial scaling): HOLDS** — the wedge identity and both convexity inequalities have the stated signs.
- **Lemma \(\ref{att:lem:cyclegeom}\) (Geometry of projected cycles): HOLDS** — after adding the omitted finiteness justification for transverse ray crossings.
- **Lemma \(\ref{att:lem:crests}\): HOLDS** — the compactness and separating-cycle argument are valid.
- **Corollary \(\ref{cor:criterion}\) (Parameter-wise criterion): HOLDS conditional on repairing Theorem \(\ref{thm:attractor}\)’s transversal gap.**
- **Proposition \(\ref{prop:crestaxis}\) (Crest return map): GAP** — parts (i)–(ii) hold; part (iii) needs a corrected strict-concavity/ratio argument.
- **Corollary \(\ref{att:cor:dulac}\) (Two constraints on Dulac functions): GAP** — the conclusion is repairable, but the proof invokes the false global transversal statement.
