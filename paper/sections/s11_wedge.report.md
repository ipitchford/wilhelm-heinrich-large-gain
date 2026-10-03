# s11_wedge.tex — writer's report (2026-10-03)

Section: "Small rate ratio: the near-Hamiltonian wedge and the relaxation limit" (`sec:wedge`).
Outcome reported in the section: exact identities, a kurtosis theorem and a **local, conditional** wedge theorem are
proved; the relaxation regime is floating-point evidence plus formal asymptotics; global uniqueness for small `mu` is
**not** proved, and the section says so in its first paragraph, in the theorem statement, in a remark and in a
closing subsection.

ROOT = `/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03`, BUNDLE = `ROOT/v3/wilhelm_heinrich_v3_0`.

## Files written

| file | content |
|---|---|
| `BUNDLE/paper/sections/s11_wedge.tex` | the section (711 lines; 19 pages in the stand-alone wrapper) |
| `BUNDLE/paper/bib/s11_wedge.tex` | `\bibitem{Fenichel79}`, `\bibitem{Jones95}` (Crossref-verified metadata) |
| `BUNDLE/paper/test_s11_wedge.tex` | compile wrapper (auxiliary files and PDF deleted after the test) |
| `BUNDLE/paper/sections/s11_wedge.check.py`, `.check.out` | re-check of every renamed formula (42 checks, all pass, 3 s) |
| `BUNDLE/paper/sections/s11_wedge.report.md` | this file |

Compile test: `pdflatex` on the wrapper, two passes: **0 LaTeX errors, 0 overfull boxes**. Undefined references are
only cross-section ones, all in the exported lists of PLAN.md: `app:hopf app:reproduce cor:cones cor:criterion
eq:adjointmatrices eq:gainprimitive eq:logsystem eq:orthogonality eq:projective eq:scalar eq:statebox eq:system
eq:variancewindow lem:collar lem:projectivebounds lem:statebox prop:crestmap prop:multipliers sec:evidence`, plus the
citations `Fenichel79`, `Jones95`. No internal (`wed:`) label is undefined or multiply defined.

## Statements included, with sources

Source lane: `ROOT/work/lanes/hamiltonian/NOTES.md` (H), referee report `ROOT/work/lanes/hamiltonian/VERIFY.md` (V),
`ROOT/work/lanes/cartography/NOTES.md` (C).

| label | statement | status in the paper | source |
|---|---|---|---|
| `lem:tauidentities` (exported) | exact identities (a)–(f) on every cycle, plus (b),(c),(f) in original variables | proved, full proof | H §1 Lemma 1; V §2 |
| `wed:prop:melnikov` | `M_1(xi) = (1-rho) oint U'^2`, no zero; ratio `oint U'^2/oint(e^U-1)^2 = rho` | proved | H §2 Prop. 2; C §4.2 |
| `thm:kurtosis` (exported) | (i) `<P^2> >= Pm^2/3`; (ii) `<P^4> <= (6/5)Pm^2<P^2> - Pm^4/5 <= (9/5)<P^2>^2`; (iii) `I0'=T0`, `I1'=3I0`, `Rsel' >= 4/5`; (iv) `Rsel` bijection | proved, full proof | H §3.2 Theorem A; V §1 |
| `wed:cor:root` | unique simple root `xi*(r)` of `M_3 = I0 (r - Rsel)`, `3r/4 <= xi* <= 5r/4`, slope `= rT0 - 2I0 = -I0 Rsel' <= -(4/5) I0` | proved | H §3.2 (last paragraph of Theorem A) |
| `wed:rem:sharp` | strictness, equality case `|U|`, floating-point kurtosis values | remark | H §3.2 Remarks, `toda.out` |
| `wed:lem:normalform` | slow-manifold coefficients `w0..w3`, modified energy, `phi_3`, `oint phi_3 = M_3` | proved (exact algebra) | H §3.1; V §3 "Algebra" |
| `wed:rem:exactselection` | bracket of (f) on the centre equals `(2/3)<P^4>` | remark (computation shown) | H §1 Remark; V §2 "Extra (X2)" |
| `wed:rem:hopf` | consistency with the Hopf coefficients of Appendix A | remark | H §3.1 "Sanity" |
| statement **(F)** | the exact property of the Fenichel slow manifold that is used, with hypothesis check | imported, cited | H §4.1 Step 0; V §3 "Step 0" |
| `thm:wedge` (exported) | local uniqueness in `Omega_{xi_b}`, crest energy `xi*(r)+O(eps)`, `log m_t = -eps^3 I0 Rsel' + O(eps^4)` | **conditional on Fenichel's theorem**, local in phase space, `mu_0` not explicit, explicitly not global | H §4.1 Theorem C; V §3 |
| `wed:rem:scope` | (a) local; (b) not uniform as `r -> 0`, `r -> infinity`; (c) divergence route to the multiplier formula | remark | H §4.2, §4.3; V §3 "Independent derivation" |
| §11.6 (`wed:sec:relax`) | floating-point tables: test of the third-order predictions; fixed-`rho` cycles; formal limit map and its fixed points; `log m_t` convergence; `Lambda` monotone on grids; three regimes and formal matching | **floating point / formal**, separated from proofs | H §2.2, §3.3, §5; V §3 (P1); C §2.1 |
| §11.7 (`wed:sec:missing`) | statement (N); a-priori window; variance ceiling `wed:eq:ceiling` (proved); exact identities `wed:eq:energyidentity`, `wed:eq:spike` (proved); the two a-priori estimates (P-T), (P-E) and the unproved spike step | open, stated as open | H §4.2, §4.2b |

## Renaming to the global notation (lane symbol -> paper symbol)

`sigma=(1+mu)/sqrt(mu)` -> not used (`1/varepsilon`); `eps` -> `\varepsilon`; `'` (d/dtau) -> dot; `G(U)` -> `g(U)=\Rem'(U)`;
`Phi` (Toda potential) -> `\Rem` (global), `Phi` kept for the general potential of `thm:kurtosis`; `E=e^U` -> written out
(E is the equilibrium); `H0` -> `\mathcal E_0`; energy level `h` -> `\xi` (h is the attractor graph); `h*(r)`, `h_a`, `h_b` ->
`\xi_*(r)`, `\xi_a`, `\xi_b`; `D_b`, `K` -> `\Omega_{\xi_b}`, `\Omega_{\xi_b+1}` (K is a gauge); `J0`, `J1` -> `I_0`, `I_1`
(J is the Jacobian); `Rsel` -> `R_{\mathrm{sel}}`; `M(h)` -> `\mathcal M_1`, `M3` -> `\mathcal M_3` (M is the crest-return
derivative); `M_eps` -> "the graph"; `X_eps` -> `\mathcal V_\varepsilon` (X is a chemical variable); `Psi` -> `\Theta`
(Psi is the weak adjoint); `sigma(U,h)` -> `\Upsilon`; `h1`, `h2`, `Ht` -> `\mathcal E_1`, `\mathcal E_2`, `\tilde{\mathcal E}`;
`F3`, `R4` -> `\phi_3`, `\phi_4` (F is the reduced planar field); `k_s` -> `N`; `D=Q-w` -> `\varpi` (D is d/dt);
crest amplitude `a` -> `\alpha` (a is a strong-stable component); section coordinate `c` and `Dsp` -> `\alpha` and `\beta`;
`g(a,eps,r)`, `g3`, `g4` -> `\chi`, `\chi_3`, `\chi_4`; `S` (section) -> `\mathcal L_+`; measure `pi`, atom variable `a`,
ratio `t` -> `\nu_1`, `s`, `\omega`; pulse map `S(q)` (H) / `G(p)` (C) -> `\mathcal G(q)`; `Sigma` (H) / `S*` (C) -> `\vartheta`
and, at the fixed point, `S(\rho)` (as in `s14_status.tex`); `F_lambda` / `Pi` -> `\Pi^{\rm rel}_\lambda`; `Lambda(q)` kept;
`q=eta-mu`, `qh` -> `\hat\eta`; `s`, `V`, `E1`, `L` -> `\mathsf s`, `\mathsf V`, `\mathsf E`, `\mathsf L`.

## VERIFY.md corrections applied

1. **C^1 -> C^2 in Step 2 of Theorem C** (V §3 "Step 2"). Applied: in the paper's Step 3 the remainder `\phi_4` and the
   displacement `\beta` are stated to be of class `C^{N-4} = C^2` (N = 6), `L` bounds first *and second* derivatives of
   `\beta`, and the text marks the single place where the second derivatives are used (the `O(eps^4)` remainder of the
   return-map slope `wed:eq:returnslope`).
2. **eps-dependent section coordinate** (V §3 "A second cosmetic point"). Resolved by parametrising the section by the
   crest abscissa `\alpha` (independent of eps) and using the modified energy only through the displacement identity
   `wed:eq:displacement`; the monotonicity of `U -> \tilde{\mathcal E}(U,0)` on a fixed larger interval is stated.
3. **Two unstated facts in Step 4** (V §3 "Step 4": one crossing of a transversal; transversality down to the origin;
   uniqueness of the planar equilibrium). Made explicit: the polar angle is strictly decreasing on the punctured ball
   for small eps (`wed:eq:angular`), a closed orbit is a Jordan curve of winding number -1 and meets the half-axis once;
   the only planar equilibrium is the origin (lift to the 3-D equilibrium).
4. V §1 remark on `J0' = T0` ("differentiates across an integrable endpoint singularity; standard"): replaced by a
   complete Tonelli argument for both `I0' = T0` and `I1' = 3 I0`, with continuity of `T0`.
5. V §4 "things a reader should not over-read": each of the four cautions is in the text (Theorem statement, Remark
   `wed:rem:scope`, Remark `wed:rem:sharp` last sentence, last sentence after `lem:tauidentities`).

No statement that VERIFY left unrefereed (H §4.2a, §4.2b, §5, §6) is presented as a theorem, with one exception that I
re-derived myself and machine-checked: the variance ceiling `wed:eq:ceiling` and the two exact identities
`wed:eq:energyidentity`, `wed:eq:spike` of H §4.2b (checks A7, A8a–A8d, B2 below). The *use* of these identities
towards (N) is presented as a formal sketch that "has not been made rigorous".

## Deviations from the lane text (content preserved, presentation changed)

- Theorem C's proof is reorganised into Steps 0–5: the polar-coordinate return map (lane Step 3) is constructed first
  (Step 2) on the whole ball, then the annulus (Step 3, lane Step 2), then the neighbourhood of the origin (Step 4).
  No step was deleted; the compactness argument for the unique zero of `\beta` and the mean-value step in Step 4 are
  written out in more detail than in the lane.
- The imported theorem is isolated as an explicit statement **(F)**: existence of `w in C^N` on a neighbourhood of a
  compact set, `w = -g` at eps = 0, satisfying the invariance equation. The paper explains how (F) follows from the
  cited theorem (normal hyperbolicity, graph over a ball, local invariance => invariance equation) and adds one
  sentence on formulations that require an overflowing-invariant base domain (collar modification outside the region
  of interest). `thm:wedge` begins "Assume (F)".
- H §4.2a (Claim C', Hopf end `r -> 0`) is **not ported** as a statement: the lane marks it "detailed outline, not
  counted as proved". The paper records it as missing item (c) in §11.7 and in Remark `wed:rem:scope`(b).
- H §4.3 (uniformity table) is condensed into Remark `wed:rem:scope`(b) and the paragraph "Three regimes".
- H §6 (Dulac hint, amplitude formulation) and H §8 are not ported (programme material; belongs to s12/s14 if wanted).
- The energy bound for the crest in §11.7 uses `eq:statebox` (`e^U = X <= M_0`) as in the lane; the lane's sentence
  "no bound on |p| at negative p" is made precise: `p >= -lambda (M_0 - 1)` from `Y <= M_0`, of order `1/mu`.

## Formulas re-checked after renaming

`s11_wedge.check.py` (output `s11_wedge.check.out`; `/Users/admin/.venvs/analysis-py313/bin/python`, 3 s, single
process). **All 42 checks pass.**

- Part A, exact (sympy): recursion and closed forms of `w0..w3`; invariance equation through `eps^3`; partials of
  `\Upsilon`; `\Theta' = g^2`; cohomological identity; the full identity `wed:eq:drift` with `\phi_3`, `\phi_4`;
  `m_2`, `m_4` and the two polynomial identities of `thm:kurtosis`(ii); the exact `\varpi` equation; energy balance;
  (f) in original variables; moment identity; `(X1) = -(orthogonality residual)`; energy identity on
  `{orthogonality = 0}`; `\zeta = \mu(1 + A L \hat\eta)` and the spike identity; the explicit pulse; quiescent-phase
  formulas and the small-`q_-` expansion; Hopf-coefficient algebra; `\partial_P w_3`.
- Part B, floating point (not interval-certified): identities of `wed:prop:melnikov` for `rho = 1, 2.5, 0.4`
  (defect <= 1.2e-12); averages `wed:eq:todaaverages` and `oint \phi_3 = r I0 - (2/3) I1` (<= 2.3e-12);
  `thm:kurtosis`(i),(ii) slacks positive; `I0' = T0`, `I1' = 3 I0`, `Rsel'` formula by central differences (<= 3.9e-7);
  Remark `wed:rem:exactselection` (6e-13); divergence route (3e-12); roots `xi*(r)` for `r = 0.3, 1, 3, 12`
  (0.301252, 1.013914, 3.121936, 13.402232, in agreement with `verify/v4_theoremC.out`) and the slope identity;
  `Rsel/xi -> 1`, `T0 -> 2 pi`; the two elementary bounds used in §11.7; Lemma `lem:tauidentities` (a)–(f), the
  original-variable forms, the moment identity and the energy identity on three computed cycles
  (`(mu,rho) = (0.3,1.5), (1,2), (0.05,2)`, worst defect 1.7e-12 for the identities, 9.2e-12 for the energy identity);
  one row of Table `wed:tab:selection` reproduced (0.99599, 0.99261); `\mathcal G(1) = 2.000000000`; small-`q` values
  of `\mathcal G`; fixed point of the limit map at `rho = 1.5` (`q* = 1.150216`, `m_0 = 0.679827`, `S = 2.43869`,
  `lim mu U_min = -1.03364`); matching constant 2.10819.

Numbers in the four tables were copied from the lane outputs and compared digit by digit with
`fixed_rho.out`, `ampcycle.out`, `relax_map.out`, `lambda_of_q.out` (hamiltonian lane), `verify/v4_theoremC.out`,
and `limit_map_stdout.txt`, `limit_map.csv` (cartography lane). Rounding: `mu T` in Table `wed:tab:fixedrho` is given
to four decimals exactly as printed (the lane's NOTES rounded 5.9425 to 5.943).

## Evidence files named in the section (for the evidence engineer)

The section names scripts and outputs by base name and points to `sec:evidence` / `app:reproduce` for their location.
Current locations: `ROOT/work/lanes/hamiltonian/` — `ampcycle.py/.out`, `toda.py/.out`, `fixed_rho.py/.out`,
`relax_map.py/.out`, `lambda_of_q.py/.out`, `variance_route.py/.out`, `verify/v4_theoremC.py/.out`;
`ROOT/work/lanes/cartography/` — `limit_map.py`, `limit_map.csv`, `limit_map_stdout.txt`.
Only `ampcycle.py` and `toda.py` are in `BUNDLE/code/v3/lib/` at the time of writing. If the bundle renames any of
these files, the `\texttt{...}` names in §11.6–11.7 must be changed to match. `s11_wedge.check.py` may be moved to
`BUNDLE/code/v3/`; it locates `whlib.py`, `toda.py`, `ampcycle.py` in `BUNDLE/code/v3/lib` first and falls back to
the lane directories.

## Bibliography

`Fenichel79` and `Jones95` were verified through the Crossref API on 2026-10-03 (queries recorded in the bib
fragment): Fenichel, J. Differential Equations 31(1), 53–98 (1979), DOI 10.1016/0022-0396(79)90152-9; Jones,
"Geometric singular perturbation theory", pp. 44–118 in *Dynamical Systems* (Lecture Notes in Mathematics, ed.
R. Johnson), Springer 1995, DOI 10.1007/BFb0095239. Crossref gives no volume number for the Jones chapter; none is
printed.

## Doubts and things I could not do

1. **Fenichel's theorem was not read in the original.** The publisher page for Fenichel (1979) returned HTTP 403 and
   an online copy of Jones's notes returned HTTP 401. The pinpoint citations "Theorem 9.1" (Fenichel) and
   "Theorems 1 and 2" (Jones) are those of the lane and agree with my recollection, but were **not checked against the
   texts**. The *form* of the statement (smooth fast–slow system; compact normally hyperbolic critical manifold that
   is a graph over a compact simply connected domain with smooth boundary; locally invariant graph, `C^r` for every
   finite `r`, `O(eps)`-close) was checked against one secondary restatement that I did read (J.-M. Ginoux,
   arXiv:2012.06770, Section II.B, which follows Jones and Kaper). That restatement lists overflowing invariance of
   the base domain as a hypothesis; this is why the paper contains the collar-modification sentence. The joint `C^N`
   dependence on `(U,P,r,eps)` **down to eps = 0** is used essentially (Steps 0, 3, 4) and rests on the cited theorem.
   Consequence for the frame: `s01_intro.tex` says "No theorem is imported from a text that was not read"; for
   Fenichel's theorem this sentence is not currently true and should be amended, or the two sources should be read.
2. `s14_status.tex` item (d) says global uniqueness in the wedge "would follow from Theorem `thm:wedge` together with
   two a-priori estimates". §11.7 is more cautious, and I believe correctly so: the two estimates give a bounded
   oscillation of `log E` on non-attracting cycles, but the final "spike" step from there to (N) is only sketched in
   the lane and is listed as missing item (b). The s14 sentence should add "and a step that is only sketched".
   Also, uniqueness for all `r in (0, r_2]` additionally needs the Hopf-end version (missing item (c)).
3. Symbol `S(\rho)` (limiting scaled period) follows `s14_status.tex`; `S` is also the signature matrix of
   `eq:adjointmatrices`. The local-symbols paragraph says so. If the frame prefers, rename to `S_*(\rho)` in both files.
4. Step 2 of the proof uses standard smooth dependence of solutions of the polar equation on initial values and
   parameters, including the initial value `\varsigma = 0`; Step 1 uses that Floquet multipliers are invariant under a
   `C^N` change of coordinates. Both are standard and are stated, not re-proved.
5. The sharpness sentence in Remark `wed:rem:sharp` claims only that equality forces a constant density (potential
   `|U|`, outside the class); the lane's stronger wording "the constant 4/5 is sharp: `Rsel' -> 4/5` as `h -> infinity`"
   for the Toda potential is supported by floating point only (kurtosis 1.7988 at `xi = 10^3`; the referee reports
   1.7999988 at `10^6`) and is reported as such.
6. The relaxation subsection is formal and floating point by design. The statement that each `q` has exactly one
   `Lambda(q)` relies on a monotonicity checked on a 500-point grid (lane H §5); the text says so.
7. Not checked by me: the referee's fixed-`r` numbers quoted in §11.6 (1.16, 4.37, 19.4; 0.9997, 0.9993, 0.9983) were
   copied from `verify/v4_theoremC.out` and not recomputed, apart from the roots `xi*(r)`.

## Cross-model review (work/sol/review_s11_response.md)

Verdicts returned: HOLDS for `lem:tauidentities`, `wed:prop:melnikov`, `thm:kurtosis`, `wed:cor:root`, `wed:lem:normalform`; GAP for `thm:wedge`. The editing pass that applied the repairs was cut off before writing this record; the dispositions below were checked by the orchestrator against the current text.

| Point | Disposition |
|---|---|
| 1 strictness in Remark `wed:rem:sharp` | Repaired: positivity of the Stieltjes measure away from 0 is now derived from strict increase of f. |
| 2 constant density and |U| | Repaired: restated for V-shaped potentials, outside the class of the theorem. |
| 3 derivation of (F) from Fenichel's theorem | Repaired by restatement: (F) is now an explicit standing assumption of Theorem `thm:wedge`; the derivation from Fenichel's theorem is labelled an unchecked sketch, and every dependent statement says it is conditional on (F). |
| 4 equilibrium as a bounded entire solution | Repaired in wording. |
| 5 convexity and radial single-valuedness of the energy ball | Repaired. |
| 6 winding number | Repaired: strict decrease of the polar angle gives winding number -m, m >= 1, and the Jordan property gives m = 1. |
| 7 monotonicity inputs | Repaired. |
| 8 return point inside the monotonicity interval | Repaired: the annulus bound on the return point is stated. |
| 9 epsilon_4 when c_4 = 0 | Repaired: epsilon_4 = min{epsilon_3, c_3/(1+c_4)}. |
| 10 crest is a maximum | Repaired: strict local maximum of U at the return point. |
| 11 "reduces to" in Remark `wed:rem:exactselection` | Repaired: described as a formal limit. |
| 12 whole cycle in the energy ball | Repaired: uniform energy distance d_b from the boundary and the energy-drift bound give containment in a smaller ball. |
| 13 hyperbolicity | Repaired: both transverse multipliers in (0,1), flow multiplier simple. |
| 14 averages in t and tau | Repaired in wording. |
| 15 differential in a display | No change needed where `\dd` is used (it expands to dt). |
| 16 "every strictly convex potential" | Repaired: the hypotheses of Theorem `thm:kurtosis` are named. |

The check script and its output were moved to `work/framework_checks/` (copied into `provenance/research_lanes/framework_checks/`). The section compiles without errors.
