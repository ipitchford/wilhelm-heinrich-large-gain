# Lane `monotone` — notes (2026-10-03)

Target of the lane: use the cooperative/convex backward flow (brief 3(j)) and the Mallet-Paret–Smith
sign structure to pin down the sign of the weak exponent lambda_c = <eta> − mu on every cycle.

**Outcome in one paragraph.** The sign of lambda_c was NOT obtained. What was obtained is a new,
fully proved structural theorem about every cycle (file `ds_theorem.tex` / `ds_theorem.pdf`):
the weak adjoint has an explicit formula, its three components have exactly two zeros each, interlaced
with the six extrema ("twelve-event order"), and the strong-stable plane span{z'(t), e_s(t)} at every
point of a cycle is a strict supporting plane of the whole cycle, of E, and of everything nested inside
("dual supporting planes", the strong-stable counterpart of v2.0 Thm 3.1). The same analysis yields a
precise obstruction: every identity obtained by pairing the weak adjoint with a periodic test function
is the *same* rank-one identity, carrying the single unknown sign of e^{−lambda_c T} − 1; and the
first-order consequences of convexity of the backward flow are swamped by the strong direction.
No counter-evidence to the conjecture was found (all 23 distinct computed cycles are attracting).

Notation: z=(U,p,y), z'=(p,p',y'), x=ke^U, A=1+mu, J the Jacobian, L3 = D^3+AD^2+mu D+x.
Six extrema (v2.0 Thm 4.2), in cyclic order: t_tr (trough of U) < t_Zmin < t_Ymin (= max of p) < t_cr (crest)
< t_Zmax < t_Ymax (= min of p) < t_tr+T. Arcs C=(t_Ymin,t_cr): p>0,p'<0,y'>0; F=(t_Ymax,t_tr+T): p<0,p'>0,y'<0.
Strong-stable Floquet vector e_s=(a,−b,−c), ratios r=b/a, s=c/a. psi_p := p' + r p. theta := −p'/p on C∪F.

---------------------------------------------------------------------------------------------------

## 1. Proved (complete proofs written out in `ds_theorem.tex`; every step also checked numerically)

**R1 (strong-stable vector).** e_s=(a,−b,−c) with a,b,c>0 for all t (Perron–Frobenius for the backward
cooperative irreducible linear system). r' = r^2 − mu r − s, s' = s r − s − x, and r > max(1,mu), s>0.
Equivalent scalar form: r'' + (A−3r) r' + r(r−1)(r−mu) = x, <r> = gamma = 1+<eta>.
(Mirror image of the adjoint equation g''+(3g−A)g'+g(g−1)(g−mu)=x for g=1+eta; so the conjecture is
"<r> < A" for the positive periodic solution r of the time-reversed equation.)

**R2 (weak adjoint formula).** Psi := e^{At} z' × e_s = e^{At} a ( r y' − s p',  y' + s p,  −(p' + r p) ) solves
Psi' = −J^T Psi, Psi·z'=Psi·e_s=0, Psi(t+T)=e^{−lambda_c T}Psi(t). Psi·v = e^{At} det(v, z', e_s).
Psi_3 = e^{At}(a'p − a p') is (minus) the Wronskian of the tangent p and the strong-stable solution a.

**R3 (twelve-event order).** Each Psi_j has exactly two simple zeros per period, interlacing the six extrema:

    t_tr < [Psi_2 ↑] < t_Zmin < [Psi_1 ↑] < t_Ymin < [Psi_3 ↑ = t_2] < t_cr
         < [Psi_2 ↓] < t_Zmax < [Psi_1 ↓] < t_Ymax < [Psi_3 ↓ = t_1] < t_tr + T.

So Psi>0 at the crest, Psi<0 at the trough; Psi_3>0 exactly on the "crest arc" (t_2,t_1) and <0 on the
"trough arc" (t_1,t_2+T); t_1 in F, t_2 in C are the phases where the (U,p)-projection of e_s is tangent
to the projected cycle (p'/p = −r). Proof idea: at a zero of theta−r, (theta−r)' = y'/p + s > 0; at a zero
of y'/p+s the derivative is s psi_p/p > 0; at a zero of y'/p' − s/r it is x psi_p/(r p') > 0.
This answers lane item 1 for Psi: the sign pattern IS forced (not by N-theory alone — in 3D, N=2 for the
adjoint only says "not in ± the cone" — but by N-theory-style Riccati crossing arguments plus the six-extrema order).

**R4 (one planar supporting line).** The tangent line of the projected cycle at phase t_1 strictly supports the
projected cycle: L(t) := p'(t_1)(U(t)−U(t_1)) − p(t_1)(p(t)−p(t_1)) > 0 for t ≠ t_1 (and L>0 inside the curve,
in particular at E). Proof: L' vanishes only where theta = r(t_1) > mu, exactly once in C and once in F
(theta' = theta^2 − mu theta + y'/p > 0 there), and L''(t_1) = p^2 r(r−mu) + p y' > 0. The same proof gives:
the tangent line at ANY phase of C∪F with −p'/p >= mu is a supporting line.

**R5 (dual supporting planes, "DS").** For every cycle, all t, all s ≠ 0 mod T:

    det( z(t+s) − z(t), z'(t), e_s(t) ) < 0,   det( E − z(t), z'(t), e_s(t) ) < 0,   det( z''(t), z'(t), e_s(t) ) < 0,

and the first inequality also holds for every point of every cycle nested inside (more generally every entire
solution whose projection is on the inner side of the tangent line of R4 at the times t_1+jT).
Proof: B := Psi·(z_2 − z) has B' = Psi_3 · x · (e^{d}−1−d) (exact chord identity), so on each period window B is
maximal at t_1, where Psi(t_1) = kappa (p', −p, 0), kappa<0, and B(t_1) = kappa · L < 0 by R4. For z'':
(Psi·z'')' = Psi_3 x p^2 and (Psi·z'')(t_1) = kappa L''(t_1) < 0.
Consequences: Wr(p',p,a) never vanishes; in the Floquet frame (z', e_t, e_s) at any point of the cycle, every chord
to another point of the cycle and the vector to E have negative e_t-coordinate (DS) and negative e_s-coordinate
(v2.0 Thm 3.1). The cycle is "doubly convex": supported by both the centre-unstable planes and the strong-stable planes.

**R6 (rank-one window identity = the exact obstruction).** For every T-periodic phi:

    ∫_{t_1}^{t_1+T} Psi_3 L3[phi] dt = (e^{−lambda_c T} − 1) |Psi_2(t_1)| (phi' + r phi)(t_1).

Hence lambda_c<0  <=>  ∫Psi_3 x > 0  <=>  ∫Psi_3 x p^2 < 0  <=>  ∫Psi_3 (1−e^U+Ue^U) < 0  <=> ∫Psi_3 x c(U(·+s)−U) < 0
for any/every s (window integrals, unnormalised Psi). All test functions give the same single sign.
Exact by-product: ∫_{t_1}^{t_1+T} Psi_3 x (p^2 + alpha_*) dt = 0 with alpha_* = [ |p|(r−mu) + |y'|/r ](t_1).
Sturm–Liouville form: with varpi = e^{At}a (varpi(t+T) = e^{−lambda_c T} varpi(t)) and
L2' := D^2 + (A−r)D + beta, beta=(r−1)(r−mu)−2r', one has L3[phi] = L2'[phi'+r phi], varpi L2'[psi] = (varpi psi')' + varpi beta psi,
L2'[psi_p]=0, L2'[r] = x > 0. The conjecture is equivalent to: there is a periodic psi with psi(t_1)>0 and
psi_p · L2'[psi] <= 0 (L2'[psi] must change sign exactly at t_1 and t_2).

**R7 (crest return map off the cycles).** On the crest section {p=0, y>k} with coordinates (U,y), the region
R = {x >= y > k} is mapped into its interior by the first return map P, and DP is entrywise positive at every
point of R (not only at cycles). Proof: for a section variation (dU,0,dy) >= 0 the wedge with z' is
((y−k)dU,(x−y)dU,(y−k)dy) >= 0, it evolves by the cooperative irreducible compound system w'=Cw along ANY orbit,
and at the return it equals ((y_1−k)dU_1,(x_1−y_1)dU_1,(y_1−k)dy_1) > 0. At a cycle crest the Perron eigenvector of
M=DP is exactly (1, zeta(0)) (leading bundle = n^⊥, n=(zeta,−eta,−1)); so any Krasnoselskii-type certificate vector v
(Mv<=v, v>0) must have slope within O(1−m_t) of zeta(crest).

**R8 (what convexity of the backward flow gives; all first order).** In coordinates (U,−p,−(y−k)) the backward
time-t map F_t is monotone and convex (cooperative field with convex components). Therefore:
(i) (DF_T(z_0) − I)(z_0 − E) >= 0, (DF_T − I)(z(s) − z_0) <= 0, (DF_T − I) z''(0) <= 0 in the cone order
[= beta(1/m_t − 1)e_t + delta(e^{gamma T} − 1)e_s with the Floquet coordinates of R5]; numerically all three vectors are
parallel to e_s to within the factor (1/m_t−1)/(e^{gamma T}−1): the weak direction is invisible at this order.
(ii) (proof sketch, not written out to referee standard) the forward flow maps closed convex lower sets (closed under
U↓, p↑, y↑) to convex lower sets; the convex epigraph of v2.0 Thm 5.1 is such a set. (iii) the scaled curve
E+theta(z−E) has residual k[theta(e^U−1)−(e^{theta U}−1)] e_3 (>=0 for theta<1, <=0 for theta>1), i.e. the flow crosses the
cone over the cycle with apex E downward inside and upward outside; for nested cycles this only reproduces
(y−k)/|(U,p)| increasing along rays, i.e. convexity of the attractor graph. No stability information.

---------------------------------------------------------------------------------------------------

## 2. Numerically supported, not proved

- **N1.** J0 := <a^ psi_p> = −<Psi^_3> > 0 on 15 cycles (ds_check.py). Equivalent: <Psi^_2> > 0 (periodic
  normalisation), i.e. d(return defect)/dk has a definite weak component, so k would be a regular parameter along
  the cycle family even at a neutral cycle. Formula: <Psi^_2> = (gamma−mu) <a^(2r−gamma)p>,
  <Psi^_1> = (mu+lambda_c)<Psi^_2>, <x Psi^_3> = −lambda_c(lambda_c+1)(lambda_c+mu) <Psi^_3>.
- **N2.** <Psi^_3 x nu^2> > 0 on 15 cycles (nu = U-component of weak Floquet vector, periodic normalisation,
  nu(crest)>0). (On non-neutral cycles this is not a coordinate-free quantity.)
- **N3.** The (U,p)-projection is strictly convex (p'^2 − p p'' > 0) and star-shaped about 0 (p^2 − U p' > 0) on
  15 cycles (twelve.py). Proved only partially: star-shapedness where U>=0 (Lambda_E>0 gives
  p^2 − U p' > p^2 − (eta_E−mu) p U + zeta_E U^2, positive definite because E has a nonreal unstable pair);
  convexity on [t_Zmin,t_Ymin] and [t_Zmax,t_Ymax] (p, p', y' of one sign: p'^2 + mu p p' + p y' > 0 termwise) and on
  the parts of C∪F where −p'/p >= mu. Open on (t_tr,t_Zmin), (t_cr,t_Zmax) and where −p'/p < mu.
- **N4.** The weak Floquet vector e_t=(nu,nu',dy) has the analogous interlacing: zeros of nu in (t_cr,t_Zmax) and
  (t_tr,t_Zmin); of dy in (t_Zmax,t_Ymax),(t_Zmin,t_Ymin); of nu' in (t_Ymax,t_tr+T),(t_Ymin,t_cr) (15 cycles).
- **N5.** DP entrywise positive at 245/245 grid points of the section including points with x<y (retmap.py);
  the U-component of P is radially concave about E (D^2 P_U[rad,rad] < 0 at 245/245 points). If the latter were
  proved it would give only the bound m_t <= 1 + K(0) M_12 / U(0), not the sign.
- **N6.** All computed cycles have lambda_c<0 (23 distinct parameter pairs, mu in [0.05,6], rho in [1.01,5]).

## 3. Refuted / shown not to work

- **X1.** Componentwise positivity of the wedges (z−E)∧z' and z'∧z'' (which would have made DS trivial): the
  (U,y)-component is negative on 6–45% of the period for rho >= 1.3 (wedges.py). DS holds anyway.
- **X2.** Radial sublinearity of P about E: violated at 5/98–95/98 test pairs (retmap.py); analytically, near Hopf
  the radial vector has slope 1/A + c h/2 against the Perron slope 1/A + c h, c=(mu^2+4mu+1)/(A D_mu), so the second
  component of Mv<=v fails at O(h) (complementing v2.0 Thm 8.1 where the first component fails for the sqrt direction).
  Any fixed scaling family must match zeta(crest) to O(h^2); only the attractor itself does.
- **X3.** Order-concavity/convexity of P in (U,y): D^2P[e_i,e_i] has both signs (retmap.py).
- **X4.** (eta−mu)^2 < 4 zeta along cycles (frozen leading block complex): false for rho >= 1.3 (twelve.py).
- **X5.** Sign of lambda_c from N-theory alone: for the adjoint system in 3D the lap number only says Psi is not in
  ± the cone; the full pattern (R3) needs the Riccati crossing argument and still leaves the sign of lambda_c free (R6).
- **X6 (argued, not computed).** Weights of the form x(beta p^2 − alpha + delta(e^{−U}−1+U)) cannot have the sign
  pattern of Psi_3 on large cycles: that needs rising speed > falling speed at equal U on [U(t_1),U(t_2)], false for rho>~2.

---------------------------------------------------------------------------------------------------

## 4. Identities verified (exact = sympy `sym_checks.py`, 25/25 pass; numeric = `identities.py`, residual <= 3.4e-6)

Exact: L3[1]=x; L3[U]=k(1−e^U+Ue^U); L3[p]=0; L3[p']=−x p^2; L3[y]=x(k+p^2); (D+1)(D+mu)(lambda−p)=x;
L3[a]=0 with a'=−ra; r''+(A−3r)r'+r(r−1)(r−mu)=x; L3[phi]=L2'[(D+r)phi]; L2'[r]=x; Psi'=−J^T Psi and Psi·z'=0 for
the formula in R2; theta'=theta^2−mu theta+y'/p; the three zero-crossing derivatives of R3; p'^2−pp''=p'^2+mu pp'+py';
compound equation w'=Cw for w=v∧u=(v2u1−v1u2, v1u3−v3u1, v2u3−v3u2); (v×u)'=(trJ−J^T)(v×u); chord remainder
F(z2)−F(z1)−J(z1)(z2−z1) = x1(e^d−1−d)e_3; scaled-curve residual.
Numeric: <r>=1+<eta>; <x a^>=chi(gamma)<a^>, <b^>=gamma<a^>, <c^>=gamma(gamma−mu)<a^>; <chi(r)>=k;
brief (f): lambda_c<c^ q+b^ y+a^ x> = <a^ x(lambda−A−2p)>; <Psi^_2>=−(1+lambda_c)<Psi^_3>; <Psi^_1>=(mu+lambda_c)<Psi^_2>;
lambda_c<Psi^_1>=<x Psi^_3>; −<Psi^_3>=<a^(2r−gamma)p>; <xPsi^_3>=−<a^ x p(2r−gamma−p)>; brief (j):
lambda_c<Psi^·z''>=−<Psi^_3 x p^2>, lambda_c<Psi^·(E−z)>=−k<Psi^_3 g>... (sign conventions as in identities.py);
e_t∧z' ∝ (1,eta,zeta); M(1,zeta(0))=m_t(1,zeta(0)); the window identities of R6 for phi=1, p', U and the alpha_* orthogonality.
Brief items checked: (f) correct; (j) correct (both identities, and convexity/monotonicity consequences C1–C3).

## 5. Lane items, answered

1. *Sign patterns.* Psi: fully determined and proved (R3). e_s: in the open cone, r>max(1,mu) (R1). e_t: analogous
   interlacing observed (N4), not proved. sigma=Psi·z'' has no zeros (R5). Crest/trough signs of Psi are forced.
2. *Second-order identity.* sigma'=Psi_3 x p^2 with sigma<0 everywhere and extremal at t_1,t_2; same for Psi·(E−z)
   and for all chords. Integrating between phases gives exactly R6: the boundary terms at t_1 have known signs,
   but the integrand changes sign at t_1,t_2, so e^{−lambda_c T}>1 is NOT forced. Equivalent reformulations:
   lambda_c<0 <=> growth of |sigma| on the trough arc exceeds its decay on the crest arc.
3. *Krasnoselskii.* P is monotone on R and P(R)⊂R (R7), but not sublinear/order-concave about E (X2, X3); the
   obstruction near Hopf is analytic. Only the U-component is radially concave (N5).
4. *Scaling comparisons.* Give the convex-graph structure and the cone inequalities of R8 only; the compound
   system sees the forcing (0,p,p')·(e_3-defect), which has the right sign only on the arcs where p and p' have the
   same sign; on C and F every e_3-defect has the wrong sign in one component. DS is the 1-D substitute.

## 6. Sharpest obstruction

(O1) By R6, pairing the weak adjoint with any periodic test function yields one and the same scalar identity; the
sign of lambda_c is the sign of a difference between the crest-arc and trough-arc integrals of |Psi_3| against a
nonnegative weight, and for different admissible weights (x versus x p^2, g(U), chord weights) the dominant arc is
different. A sign proof by this route needs a test function whose L3-image changes sign exactly at the two
tangency phases t_1,t_2, which are defined through the non-explicit strong-stable ratio r.
(O2) The standard cone (in which the backward flow is monotone and convex) has e_s as Perron direction; all its
first-order consequences are dominated by the e_s-component. The weak direction is the Perron direction of the
compound (= positive adjoint) system, for which no convexity/sublinearity in a fixed scaling family is available:
near Hopf a certificate direction must agree with zeta(crest) to O(h^2).
(O3) Near Hopf, <Psi^_3 x p^2> is O(h^4) with cancellation at O(h^3): any inequality proof must be exact through
the order that encodes the first Lyapunov coefficient.

## 7. What I would try next

1. Fold-direction programme (sketch, NOT proved): at the first neutral parameter k* on a fixed-mu path the neutral
   cycle is either degenerate or a birth fold, which forces sign(u·D^2P[v,v]) = −sign(u·∂_kP) (u, v the left/right
   Perron vectors of M). So the conjecture follows from "u·∂_kP and u·D^2P[v,v] have the same sign at every neutral
   cycle". Numerically the non-neutral analogues J0 (N1) and <Psi^_3 x nu^2> (N2) are both positive on all computed
   cycles. Flow formulas at a neutral cycle need care because the unit multiplier generically has a Jordan block
   (nu(t+T)=nu(t)+tau' p(t), and then ∫Psi_3 x p nu = tau' sigma ≠ 0).
2. Prove N1 (J0>0): it would make k a regular parameter along the whole cycle family for each mu
   (stable <=> crest moves outward with k), a clean global-continuation statement.
3. Dual certificate problem suggested by R1: find a positive periodic r̄ with r̄''+(A−3r̄)r̄'+chi(r̄) >= x,
   r̄^2−mu r̄−r̄' > 0 and <r̄> <= A (time-reversed analogue of brief (f); Collatz–Wielandt for the backward system).
   Equivalent renewal form: lambda_c<0 <=> spectral radius of phi ↦ g3*(x phi) on T-periodic functions is < 1,
   g3 = density of Exp(A)+Exp(1)+Exp(mu) scaled by 1/k_H ("R0<1 in a periodic environment", while <x>/k_H = rho > 1).
4. Use DS quantitatively: Psi_3 = −e^{At}a^2 (p/a)' and the crest-to-crest formula
   ∫Psi_3 x p^n dt = (n+1)^{-1} ∫ e^{At} a x p^{n+1} (A + p − (n+2) r) dt reduce everything to a^-weighted moments.

## 8. Reproduction (each run < 30 s, single process)

    PY=/Users/admin/.venvs/analysis-py313/bin/python ; cd ROOT/work/lanes/monotone
    $PY mono_lib.py          # library self-test (Floquet data: e_s ratios by backward Riccati, e_t, Psi)
    $PY signs.py             # zero phases of Psi, e_t, sigma, Psi.(E-z) on 15 cycles           -> signs.out
    $PY wedges.py            # X1                                                                -> wedges.out
    $PY retmap.py            # item 3: monotonicity / sublinearity / second differences of P     -> retmap.out
    $PY stableplane.py       # DS for all chords (12 cycles) + chord identities                  -> stableplane.out
    $PY ds_check.py          # every step of the DS proof + J0, I2 on 15 cycles                  -> ds_check.out
    $PY twelve.py            # twelve-event order, e_t analogue, r>max(1,mu), q1,q2,q3           -> twelve.out
    $PY ds_independent.py    # DS via monodromy eigenvectors (8 further cycles), P on R          -> ds_independent.out
    $PY convex_cone.py       # R8(i) cone inequalities                                           -> convex_cone.out
    $PY sym_checks.py        # 25 exact identities                                               -> sym_checks.out
    $PY identities.py        # averaged / window identities                                      -> identities.out
    pdflatex ds_theorem.tex  # the proved theorem (5 pages)

Numerical caveats: floating point only. For rho >~ 3 and mu <= 0.05 the trough is so deep that signs of
quantities multiplied by x (e.g. y' near the trough, dy_t) are roundoff noise (visible in signs.out for mu=0.05);
return-map derivatives after very long excursions are below roundoff and are reported as "unresolved" in
ds_independent.out. Web searches used: 0.
