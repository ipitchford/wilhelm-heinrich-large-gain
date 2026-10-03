# v3.0 paper and bundle — binding plan for all writers (2026-10-03)

## Outcome being written up (do not overstate)

The all-parameter uniqueness conjecture was **not** proved. The paper is a stand-alone account of what is proved:

1. the v2.0 framework (positive adjoint, supporting planes, cycle geometry, compactness, neutral-cycle obstruction);
2. the **attractor theorem**: the global attractor is a convex C^1 graph over a planar Jordan domain; exact planar reduction; uniqueness holds at a parameter whenever every cycle there has lambda_c < 0;
3. **exact stability formulas**: every projected cycle is strictly convex and star-shaped; lambda_c = −<pK/Sigma> = −<pW/Xi>; the scalar Riccati reduction; weak-adjoint twelve-event order; dual supporting planes;
4. **large-gain uniqueness**: for 0 < mu <= 1 and k >= 1000 (and k >= 1000 mu^3 for mu >= 1) every cycle is hyperbolic attracting, hence unique;
5. a local small-rate-ratio wedge theorem (conditional on Fenichel's theorem), negative results about certificates, floating-point cartography, and the remaining programme.

Status words to use: "proved" (complete proof in the paper), "conditional on [cited theorem]", "floating-point evidence", "refuted", "open". Never write or imply that the conjecture is solved.

## Fixed paths

- ROOT = `/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03`
- BUNDLE = `ROOT/v3/wilhelm_heinrich_v3_0`
- Paper sections: `BUNDLE/paper/sections/` (one file per section, **no preamble, no \begin{document}**)
- Shared preamble (read it; never redefine its macros): `BUNDLE/paper/preamble.tex`
- v2.0 source of the framework sections: `ROOT/v2_0_source/wilhelm_heinrich_finite_parameter_v2_0/paper/manuscript.tex`
- Lane material: `ROOT/work/lanes/<lane>/` with `NOTES.md`, proofs (`*.tex`), scripts, and the adversarial referee's `VERIFY.md` (apply every correction it lists)
- Lane result digest: `ROOT/work/attack_results.json`
- Shared numerics: `ROOT/work/whlib.py`; Python `/Users/admin/.venvs/analysis-py313/bin/python`
- Machine is short of RAM: single-process numerics, runs under ~5 minutes, otherwise wrap with `/Users/admin/.local/bin/research-run run --cpu 1 --memory-gb 1 --label wh-<name> -- <cmd>`.

## Section map, owners and exported labels

Other sections may reference exactly these labels. Internal labels must carry the section prefix shown.

| File | Content | Prefix | Exported labels |
|---|---|---|---|
| `s01_intro.tex` (frame, later) | model, coordinates, main results, positioning, notation | — | `eq:physical eq:parameters eq:filters eq:change eq:system eq:logsystem eq:scalar sec:model thm:main` |
| `s02_adjoint.tex` | positive adjoint, exact exponent (v2 §2) | `adj:` | `sec:adjoint eq:adjointmatrices eq:ratios eq:projective eq:etaminbasic prop:multipliers eq:multipliers eq:exponent prop:factor eq:orthogonality eq:tangent eq:factorisation lem:projectivebounds eq:gainprimitive` |
| `s03_support.tex` | supporting planes for cycles (v2 §3) | `sup:` | `sec:support thm:support eq:delta eq:supportidentity cor:cones eq:positiveK eq:equilibriumreverse eq:LambdaE eq:cones lem:curvature eq:curvature` |
| `s04_geometry.tex` | one crest, six extrema, crest map (v2 §4) | `geo:` | `sec:geometry thm:geometry eq:crestidentities thm:six eq:six cor:orderedcrests eq:crestorder prop:crestmap eq:compound eq:determinant` |
| `s05_gauges.tex` | star/convex projection, exact exponent formulas, scalar reduction, further identities | `gau:` | `sec:gauges thm:gauges eq:gaugeodes eq:exactexponent prop:kappa eq:kappariccati prop:polar prop:nestedradii prop:monotoneintegrals prop:holder prop:lvchain prop:strongidentity prop:gaussmap` |
| `s06_weakadjoint.tex` | strong-stable vector, weak adjoint, twelve-event order, dual supporting planes, rank-one identity, monotone crest map | `wk:` | `sec:weakadjoint lem:strongstable thm:twelve thm:dualsupport prop:rankone prop:crestmonotone` |
| `s07_attractor.tex` | attractor theorem and planar reduction | `att:` | `sec:attractor thm:attractor thm:twoorbit lem:absorbing prop:riccatientire prop:crestaxis cor:criterion` |
| `s08_compactness.tex` | uniform compactness of cycle families (v2 §6, restructured) | `cpt:` | `sec:compact lem:statebox eq:statebox eq:twosidedbox lem:eqseparation eq:anglebound lem:speed eq:speedrelations thm:compact` |
| `s09_obstruction.tex` | Hopf collar, neutral locus, count preservation (v2 §7) | `obs:` | `sec:obstruction lem:collar eq:neutralset thm:neutral cor:neutralconditions eq:neutralcriteria eq:variancewindow` |
| `s10_largegain.tex` | large-gain uniqueness theorem with full proof and constants | `lg:` | `sec:largegain thm:largegain cor:lgunique cor:lgreciprocal lem:lgshortlobe lem:lglinear` |
| `s11_wedge.tex` | exact tau-time identities, kurtosis bound, small-rate-ratio wedge, relaxation limit (numerical) | `wed:` | `sec:wedge lem:tauidentities thm:kurtosis thm:wedge` |
| `s12_negative.tex` | failed certificates and proved obstructions (v2 §8 plus new) | `neg:` | `sec:negative thm:failed eq:radial eq:radialfailure prop:nodulac prop:nomonomial lem:rigidity` |
| `s13_evidence.tex` (evidence engineer) | executed computations, counts, failure semantics | `num:` | `sec:evidence` |
| `s14_status.tex` (frame, later) | status, remaining region, programme | — | `sec:status` |
| `appA_hopf.tex` | Hopf coefficient (v2 App. A) | `hopf:` | `app:hopf` |
| `appB_reproduce.tex` (evidence engineer) | reproduction boundary | — | `app:reproduce` |

Citation keys available (bibliography is merged later; add any new entry as a `\bibitem` line with verified metadata in `BUNDLE/paper/bib/<section>.tex`, using only sources recorded in `ROOT/work/lanes/literature/NOTES.md` or the v2.0 bibliography, or verified through the Crossref API): `WH95 WH96 Smith12 BH23 MPS90`.

## Global notation (binding)

- Parameters: `\mu`, `k`, `A=1+\mu`, `\lambda=k/\mu`, `\kH=\mu(1+\mu)`, `\rho=k/\kH`, `\delta=k-\kH`. Oscillatory domain `\cP`.
- State `(U,p,y)`, `x=ke^U`, equilibrium `E=(0,0,k)`; chemical `X=x/k, Z=y/k, Y=1-p/\lambda`; `\ell:=\lambda-p`.
- Time `t` is the fast time of eq. (1.3); `\tau=\sqrt\mu\,t` is frequency-one time. Primes are d/dt unless a section says otherwise in its first paragraph. `D` is the operator d/dt in operator identities; `\Dt` is the derivative of a state function along the flow.
- A *cycle* is a nonconstant positive periodic orbit with least period `T`; `\avg{f}` is the time average.
- Adjoint ratios `\eta,\zeta`; `\gamma=1+\avg\eta`; weak exponent `\lambda_c=\avg\eta-\mu`; multipliers `m_s<m_t`.
- Convexity gauges: `K=\zeta U-\eta p-y+k`, `W=\zeta'p-\eta'p'`, `\LamE=(y-k)-\zeta_EU+\eta_Ep`.
- `\psiU(U)=Ue^U-e^U+1`, `\Rem(d)=e^d-1-d`.
- **Star gauge** `\Sg:=p^2-Up'` and **convexity gauge** `\Xg:=p'^2-pp''` (lanes called these V, S*, S1 and Q, V, S4: rename). `\kappa:=K/\Sg`. `\etaz:=[U(x-y)-p(y-k)]/\Sg`.
- Strong-stable Floquet vector `e_s=(a,-b,-c)`, `a,b,c>0`; its ratios are written `\mathfrak r=b/a`, `\mathfrak s=c/a` (not r, s). Weak adjoint `\Psi`.
- Attractor `\cA`, chemical domain `\Dom`, projection `\pi(U,p,y)=(U,p)`, `\Dsc=\pi(\cA)`, graph function `h` (convex extension `H`), reduced planar field `F(u,v)=(v,k-\mu v-H(u,v))`.
- Crest-return derivative `M` is in section coordinates `(U,y)`; in coordinates `(U,z)`, `z=y/k`, write `M_{(U,z)}=\diag(1,1/k)\,M\,\diag(1,k)`.
- Compound matrix `C` as in (4.x) `eq:compound`. Cones `C_\pm` as in `eq:cones`.
- A section may introduce local symbols, listed in a short paragraph at its start; a local symbol must not reuse a global one with a different meaning.

## Writing rules

- Stand-alone paper. No references to "lanes", "the brief", "agents", or to v2.0 as something the reader must consult. Cross-reference by label.
- UK spelling. Academic "we". No second person. No promotional adjectives. No rhetorical questions, no "it is worth noting", no "crucially/importantly/notably" openers, no summary sentence that restates the paragraph. State what is proved, by what, under which hypotheses.
- Every theorem carries a complete proof or an explicit statement of the imported result with a citation and a check that its hypotheses hold. Port proofs from the lane files faithfully; do not shorten a proof by deleting steps. Apply every correction recorded in the lane's VERIFY.md and say so in your report file.
- Floating-point evidence is never mixed into a proof. Refer to it as "floating-point evidence" and point to Section `sec:evidence`.
- After any renaming of symbols in a formula, re-check the formula (sympy or a numeric spot check with whlib) and record the check.
- Compile-test your section with a wrapper `BUNDLE/paper/test_<file>.tex` (`\documentclass[11pt,a4paper]{article}\input{preamble}\begin{document}\input{sections/<file>}\end{document}`); undefined cross-section references are expected, LaTeX errors are not. Delete the wrapper's auxiliary files afterwards, keep the wrapper.
- Write `BUNDLE/paper/sections/<file>.report.md`: statements included, source file and location for each, VERIFY corrections applied, formulas re-checked, anything you could not port and why.
- Do not edit files owned by another section.
