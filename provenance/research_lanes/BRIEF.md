# Wilhelm–Heinrich uniqueness conjecture — shared brief (2026-10-03)

Read this whole file before doing anything. It is the only context you share with the other lanes.

## 0. Fixed paths and rules

- ROOT = `/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03`
- v2.0 paper (LaTeX, 686 lines, read it): `ROOT/v2_0_source/wilhelm_heinrich_finite_parameter_v2_0/paper/manuscript.tex`
- Shared numerics: `ROOT/work/whlib.py` (import with `sys.path.insert(0, 'ROOT/work')`). Do **not** edit it; copy if you need a variant.
- Python: `/Users/admin/.venvs/analysis-py313/bin/python` (numpy 2.4, scipy 1.17, sympy 1.14, mpmath, python-flint 0.9).
- Your lane directory: `ROOT/work/lanes/<lane>/`. Put every script, data file and your `NOTES.md` there. Write nothing outside it.
- Machine is shared and has ~3 GB free RAM. Single-process numerics only, each run under ~5 minutes and ~1 GB. Anything longer must be wrapped: `/Users/admin/.local/bin/research-run run --cpu 1 --memory-gb 1 --label wh-<lane> -- <command>`.
- Honesty rule (hard): a statement is "proved" only if every step is written out and you would defend it to a hostile referee. Test every proposed lemma numerically **before** trying to prove it, and report refutations plainly. A precise obstruction is a valuable result; a hopeful sketch labelled as a proof is a harmful one. Never weaken the target silently.
- Never use the semantic-scholar `snippet_search` tool.

## 1. The problem

Scaled system (paper eq. 1.4–1.5), parameters mu > 0, k > 0, A = 1 + mu, lambda = k/mu:

    x' = p x,   p' = k − mu p − y,   y' = x − y,        x = k e^U,  U' = p.

Chemical variables: X = x/k, Z = y/k, Y = 1 − p/lambda (so Y > 0 means p < lambda), with
X' = lambda X(1−Y), Y' = mu(Z−Y), Z' = X − Z. Equilibrium E: (U,p,y) = (0,0,k).
Scalar form: U''' + A U'' + mu U' + k(e^U − 1) = 0. Hopf at k_H = mu(1+mu); write rho = k/k_H.
Oscillatory region P = {k > k_H}. Reciprocal symmetry mu -> 1/mu (frequency-one time) preserves Floquet multipliers, so mu in (0,1] suffices for stability questions.

**Conjecture (open since Smith 2012, restated Boros–Hofbauer 2023):** for every (mu,k) in P there is exactly one positive nonconstant periodic orbit.

## 2. What v2.0 already proves (take as given; proofs are in the manuscript)

Along any cycle there are positive T-periodic adjoint ratios eta, zeta with

    eta' = zeta + (mu−1) eta − eta^2,    zeta' = x − (1+eta) zeta,   eta > max(0, mu−1).

- Transverse multipliers: m_s = exp(−(1+<eta>)T) < m_t = exp((<eta> − mu)T). Weak exponent **lambda_c = <eta> − mu**. Attracting iff <eta> < mu.
- Tangent orthogonality: zeta p − eta p' = x − y (= y'). Tangent equation p'' = (eta − mu) p' − zeta p.
- Factorisation D^3 + A D^2 + mu D + x = (D + 1 + eta)(D^2 + (mu − eta) D + zeta).
- Moment identity: (zeta + eta + eta^2/2)' = x − f(eta), f(s) = s(s+1)(s+1−mu); hence k = <f(eta)>, and Jensen gives f(<eta>) < k. Note f(mu) = k_H.
- Supporting planes (Thm 3.1): for a cycle or E as base z1 with its (eta,zeta), and any other log-bounded entire solution z2,
  Delta = zeta(U2−U1) − eta(p2−p1) − (y2−y1) satisfies Delta' + (1+eta)Delta = −x1[e^{U2−U1} − 1 − (U2−U1)] and Delta < 0.
- Consequences: K := zeta U − eta p − y + k > 0 with K' + (1+eta)K = k(U e^U − e^U + 1);
  W := zeta' p − eta' p' > 0 with W' + (1+eta) W = x p^2;
  with gamma_E the root > A of gamma(gamma−1)(gamma−mu) = k, eta_E = gamma_E − 1, zeta_E = eta_E(eta_E+1−mu):
  Lambda_E := (y−k) − zeta_E U + eta_E p > 0 on every cycle (explicit linear function!), and Lambda_E' = −gamma_E Lambda_E + k(e^U − 1 − U).
- Geometry: every cycle has an injective clockwise projection to the (U,p) plane, one crest/trough, cycles are nested, crests ordered; crest-return derivative M (section p=0, coords (U,y)) is entrywise positive with det M = exp(−AT) < exp(−4 pi); T > 2 pi/sqrt(mu).
- Compactness of all cycles on compact subsets of P; cycle count locally constant off the neutral set N = {some cycle has <eta> = mu}; count = 1 near Hopf. Hence non-uniqueness needs a neutral cycle (m_t = 1) somewhere on a fixed-mu path from Hopf.
- Moment identities on any cycle: <x> = k, <p> = 0, <p'^2> = mu <p^2>, <x p^2> = A <p'^2> = mu A <p^2>.

**Remaining obligation = the conjecture:** show <eta> < mu on every cycle (or otherwise exclude non-uniqueness).

## 3. New observations from this session (derived by the orchestrator; NOT independently verified — check before relying on them)

(a) *Pointwise criterion.* For k > k_H every positive orbit off the 1-D stable manifold of E converges to a cycle (Smith 2012). If at a given (mu,k) every cycle is hyperbolic attracting, the open basins partition a connected set, so the cycle is unique. So it suffices, parameter by parameter, to show every cycle has <eta> < mu. No continuation needed.

(b) *Attractor is a convex graph (believed provable).* The positive adjoint direction exists along every bounded entire orbit (cone contraction / pull-back), so Thm 3.1 extends to any two distinct bounded entire orbits, in both directions. Hence the global attractor 𝒜 (all log-bounded entire orbits) projects injectively to the (U,p) plane; 𝒜 = graph of a convex C^1 function y = h(U,p) over the closed disc D bounded by the outermost projected cycle, with grad h = (zeta, −eta). The flow on 𝒜 is the planar C^1 flow
    u' = v,  v' = k − mu v − h(u,v),   divergence = eta − mu,
  and h satisfies h_u v + h_v (k − mu v − h) = k e^u − h. E is a source of the planar flow (unstable focus). K = k − [h − u h_u − v h_v] >= 0 and W = Phi^T D^2h Phi >= 0 are convexity gauges.

(c) *Exact formulas for the exponent on a cycle* (log-derivatives 𝔻 log F along the flow):
    𝔻log zeta = x/zeta − (1+eta);  𝔻log eta = zeta/eta + mu − (1+eta);  𝔻log K = k(Ue^U − e^U + 1)/K − (1+eta);  𝔻log W = x p^2/W − (1+eta);
    𝔻log x = p;  𝔻log y = x/y − 1;  𝔻log(lambda − p) = y/(lambda − p) − mu;  𝔻log Lambda_E = k(e^U−1−U)/Lambda_E − gamma_E.
  Hence lambda_c = <x/zeta> − A = <zeta/eta> − 1 = <k(Ue^U−e^U+1)/K> − A = <x p^2/W> − A, and <x/y> = 1, <y/(lambda−p)> = mu.
  With q = eta − mu: k − k_H = (A^2+mu)<q> + 2A<q^2> + <q^3>, so stability <=> <q^2(2A+q)> > k − k_H.

(d) *Dulac formulation.* By (b), uniqueness follows if there is a function Lambda on D \ {E} (may be singular at E) with 𝔻Lambda + eta − mu of one fixed sign on D \ {E}. Candidates: Lambda = a log W + b log zeta + c log eta + d log K + (functions of state), a+b+c+d = 1 cancels eta. Near E any such B = e^Lambda must behave like r^{-2} (K and W both vanish quadratically at E).

(e) *Supersolution criterion in 3D* (no reduction needed). Compound system w' = C w, C = [[−mu,1,0],[0,−1,1],[x,0,−A]], Perron exponent = lambda_c. If omega(z) > 0 is a C^1 vector function of the state on a set containing all cycles (singular at E allowed) with 𝔻omega >= C omega along cycles, strict somewhere, then lambda_c < 0 for every cycle. Scalar form: phi > 0, (D+mu)phi > 0, (D+1)(D+mu)phi > 0 and Lc[phi] := (D+A)(D+1)(D+mu)phi − x phi >= 0.
  Explicit values: Lc[1] = k_H − x;  Lc[lambda − p] = x(2p + A − lambda);  Lc[p] = A(k − x) − 2 x p.
  Any regular omega fails pointwise at E (there the inequality reads k <= k_H), so a certificate must use the average or be singular at E. The inequality is tight to leading order at Hopf (it must encode the sign of the first Lyapunov coefficient).

(f) *Explicit trial with the neutral mean.* (eta~, zeta~) := (mu Z/Y, mu X/Y) = (y/(lambda−p), x/(lambda−p)) satisfies the first Riccati equation exactly, has <eta~> = mu exactly, and misses the second by zeta~(2p + A − lambda). So lambda_c = <eta − eta~>. Equivalent exact identity: with a(t) > 0 the U-component of the strong-stable Floquet solution (a' = −b, b' = −mu b − c, c' = −x a − c, all > 0, decay e^{−gamma t}, gamma = 1+<eta>) and a^ = e^{gamma t} a periodic,
    sign(lambda_c) = sign < a^ x (lambda − A − 2p) >,   and   gamma(gamma−1)(gamma−mu) = <x a^>/<a^>  (so stable <=> a^-weighted mean of x is below k_H).

(g) *Necessary condition.* Jensen on logs gives gamma(gamma−1)(gamma−mu) >= k e^{<U>}; so every attracting cycle has <U> < −log rho.

(h) *Cross-parameter comparison.* In V = log x coordinates the vector field depends on k only through the constant in p'. For k2 >= k1, every bounded entire solution at k2 lies strictly above every supporting plane of every cycle at k1 (one direction only): Delta' + (1+eta1)Delta = −eta1(k2−k1) − x1[...] < 0.

(i) *Distributed-delay logistic form.* (D+1)(D+mu)U' = −k(e^U − 1), i.e. x' = x(lambda − G*x) with G the two-stage kernel (convolution of exponentials with rates 1 and mu). The conjecture is the two-stage linear-chain analogue of Jones' conjecture for Wright's equation (proved by Jaquette 2019 with computer assistance plus Xie's large-parameter analysis).

(j) *Monotone convex structure.* In coordinates (a,b,c) = (U, −p, −(y−k)) the backward flow is cooperative and its vector field is componentwise convex, so backward time-t maps are monotone and convex for the standard cone. Second-order identity: for the weak adjoint Floquet solution Psi (Psi' = −J^T Psi, Psi(t+T) = e^{−lambda_c T} Psi(t), Psi·z' = 0), d/dt(Psi·z'') = Psi_3 x p^2, and for any adjoint solution Psi and any second solution z2: d/dt[Psi·(z2 − z1)] = Psi_3 x1 [e^{U2−U1} − 1 − (U2−U1)].

(k) *Liénard plane.* s := p − (y−k) + mu U has s' = −k(e^U − 1). On 𝒜, p = P(U,s) is jointly convex with P_s = 1/(1+eta), P_U = (zeta − mu)/(1+eta). Radial (Massera) scaling in the Filippov plane fails at U = 0: the relevant quantity is (1+eta)E = K − (zeta − mu) chi(U), chi(U) = U − 2G/g >= 0.

(l) *Numerics (floating point).* For mu in {0.05, 0.2, 1} and rho from 1.02 to 20 the computed cycle has m_t decreasing monotonically in rho: e.g. mu=1: rho=1.1 -> 0.78, 2 -> 0.145, 5 -> 5.6e-4, 8 -> 8.5e-7; troughs become astronomically deep (U_min ~ −285 at rho=8, −5e4 at rho=100). All 31 bundle cycles are attracting. The hard regime is near Hopf, where lambda_c -> 0−. `find_cycle` in whlib is unreliable for mu <= 0.01 near Hopf (warm-up lands on E); use amplitude-fixed BVP there (see the bundle's `code/cycles.py`).

## 4. whlib quick reference

    from whlib import *
    c = find_cycle(mu, k)              # dict: U0,y0 (crest), T, M (2x2), m_t, m_s, lambda_c, residual
    pr = cycle_profile(mu, k, c)       # arrays t,U,p,y,x,eta,zeta,K,W ; pr['lambda_c_riccati']
    orb = unstable_manifold_orbit(mu, k)   # spiral on W^u(E) with eta,zeta,K,W (fills the attractor disc)
    cs = continue_in_rho(mu, [1.05,1.2,...])
    equilibrium_data(mu,k)             # gammaE, etaE, zetaE, unstable eigenpair

## 5. What to return

Write `ROOT/work/lanes/<lane>/NOTES.md` containing: exact statements attempted; for each, status (proved with full proof / sketch / numerically supported / refuted) and the numerical tests run with commands; the sharpest obstruction you hit; what you would try next. Then return the structured summary requested in your prompt. Keep the returned summary short; the detail lives in NOTES.md.
