# Report: s10_largegain.tex (Section 10, "Uniqueness at large gain")

ROOT = `/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03`, BUNDLE = `ROOT/v3/wilhelm_heinrich_v3_0`.

Source: `ROOT/work/lanes/largegain/theorem.tex`, lines 33-584 (everything from `\section{Statement}` to the last
remark). The lane's `VERIFY.md` re-derived every lemma and found no error; it lists two cosmetic remarks, both applied.

The port is mechanical and reproducible:

    section = HEADER + INTRO + T(S(body))

- `ROOT/v3/checks/s10_largegain_port.py` regenerates the section file from the lane source. `S` is a list of
  exact-string structural edits (each must match exactly once); `T` is the list of token renames.
- `ROOT/v3/checks/s10_largegain_check.py` (output `s10_largegain_check.out`, 218 checks, all pass) verifies the port.
- No proof was retyped, shortened or paraphrased. 487 of the 553 lines of the token-normalised source `T(body)` are
  carried over verbatim; the remaining differences are the 34 hunks printed in the check output and listed in §5.

## 1. Statements included

| Label | Content | Source lines |
|---|---|---|
| `lg:eq:sys` | system in `(U,p,w)`, `w=y-k` | 34-41 |
| `thm:largegain` | (i) `T >= 0.1833 lambda/mu`; (ii) `p >= 0.819 lambda` at the up-crossing; (iii) `{x>c_0}` is one interval mod `T` of length `<= (71.79+12.28 log(lambda/1000))/(0.819 lambda)`; (iv) `<eta> <= 0.83 mu`, so `0<m_s<m_t<1` | 43-54; proof 506-522 |
| `cor:lgunique` | uniqueness and attraction off `W^s(E)` for `0<mu<=1`, `k>=1000` | 56-60; proof 524-538 |
| `cor:lgreciprocal` | (a) reciprocal map, stated explicitly; (b) `mu>=1`, `k>=1000 mu^3`; `rho_0(mu)` | 62-71; proof 540-557 |
| `lg:sec:imports` | imported facts (F1)-(F8), events, fixed functions, `lg:eq:v2` | 81-123 |
| `lg:lem:id` | mass identity; `<(x-k_H)p^2>=0` | 126-139 |
| `lem:lgshortlobe` | `T_+ < pi/omega`; `gamma_E-mu < k^{1/3}+1`; `omega >= omega_0(k)`; `T_- > T_0(k)` | 141-164 |
| `lg:lem:signs` | signs of `w`, `p`, `phi` on the two lobes | 166-180 |
| `lg:lem:i`, `lg:lem:ii`, `lg:lem:rise`, `lg:lem:stage1` | phases (i), (ii), rising phase, Stage 1 (`p_k >= 0.374 k`) | 183-275 |
| `lg:eq:wm`, `lg:lem:P1`, `lg:lem:mass` (`lg:eq:Etheta`, `lg:eq:Q`, `lg:eq:massfinal`) | spike: `P_1` bound and mass from slope | 278-328 |
| `lg:lem:stage2` (`lg:eq:pkT`), `lg:lem:boot` (`lg:eq:Tm`, `lg:eq:hz`) | Stage 2 and the bootstrap to `T_- > 0.1833 lambda/mu`, `p_k >= 0.819 lambda` | 331-395 |
| `lg:lem:crest` | crest duration, including the monotonicity paragraph | 398-476 |
| `lem:lglinear` | linear estimate for the weak exponent | 479-504 |
| `lg:tab:constants` (subsection `lg:sec:constants`) | table of constants | 560-574 |
| `lg:rem:scope`, `lg:rem:mu` (subsection `lg:sec:scope`) | closing remarks on the threshold | 73-79 and 576-584 |

Counts checked by the script: 1 theorem, 13 lemmas, 2 corollaries, 2 remarks, 16 proofs in the source and in the section.
All 28 displayed formulas of the source are present (26 verbatim after renaming; 2 reformatted, see §5). Of the 722
inline formulas, 696 are verbatim and 26 belong to the rewritten passages listed in §5 (the script asserts this exact
set and that each replacement is present). All 78 distinct decimal constants of the source proof are present; the only
numeric literals that disappear are the v2.0 location numbers (1.5, 2.1, 3.2, 4.1, 4.4, 6.1, 2012).

Exported labels (PLAN.md): `sec:largegain thm:largegain cor:lgunique cor:lgreciprocal lem:lgshortlobe lem:lglinear`,
all present. Every other label has the prefix `lg:`. No reference is dangling: 28 labels of other sections are
referenced, all in PLAN.md's exported lists and all defined in the written sections except `sec:evidence` (s13, not
yet written).

## 2. Macros of the source and their mapping

| Source | Mapping |
|---|---|
| `\avg{.}` | preamble `\avg` (identical definition) |
| `\dd` = `\,\mathrm d` | the preamble `\dd` is `\,\mathrm dt`, so the source macro is expanded inline: `\,\mathrm d` (23 uses); `\mathrm d` without thin space in `dC_tau/d log lambda` (2 uses) |
| `theorem`, `lemma`, `corollary`, `remark` environments | preamble environments (numbered within the section) |
| `\setlist`, packages, title, author, abstract-like status paragraph, `\begin{document}` | removed |
| `k_H`, `\Lambda_E`, `\mathcal D`, `\mathbb R` (written out in the source) | `\kH` (10), `\LamE` (4), `\Dom` (15), `\R` (4) |

## 3. Label mapping

`thm:main -> thm:largegain`, `cor:unique -> cor:lgunique`, `cor:recip -> cor:lgreciprocal`,
`lem:short -> lem:lgshortlobe`, `lem:lyap -> lem:lglinear`; every other source label `X -> lg:X`
(`eq:sys eq:v2 eq:wm eq:Etheta eq:Q eq:massfinal eq:pkT eq:Tm eq:hz lem:id lem:signs lem:i lem:ii lem:rise
lem:stage1 lem:P1 lem:mass lem:stage2 lem:boot lem:crest rem:scope rem:mu`). New labels: `sec:largegain`,
`lg:sec:statement imports identities negativelobe spike bootstrap crest linear constants scope`, `lg:tab:constants`.

## 4. Symbols renamed (T) and symbols kept

Renamed because the source symbol has another meaning in the global notation of PLAN.md / Section 1:

| Source | Paper | Reason |
|---|---|---|
| `\ell=2`, `\varphi_\ell`, `s_\ell`, `r_\ell` | `L_\ast`, `\varphi_\ast`, `s_\ast`, `r_\ast` | `\ell=\lambda-p` |
| `D` (trough depth), `U(s_0)=-D` | `\|U_{\min}\|`, `U(s_0)=U_{\min}` | `D=d/dt`; `U_min`, `U_max` as in s05/s14 |
| `h` (crest value of `U`) | `U_{\max}` | graph function `h` |
| `h(z)` (function of the bootstrap) | `f(z)` | graph function `h` |
| `H=\log M_0`, `\hat H`, `\tilde H`, `\tilde H_\lambda` | `R`, `\hat R`, `\tilde R`, `\tilde R_\lambda` | convex extension `H`; adjoint solution `\mathbf H` (kept, same meaning) |
| `K_0` | `k_0` | gauge `K` |
| `W_c`, `\hat W` | `w_c`, `\hat w` | gauge `W` (they are `w` at the crest and its lower bound) |
| `X_\theta` | `N_\theta` | chemical variable `X=x/k` |
| `\Psi_\mu`, `\Psi_1` | `\Phi_\mu`, `\Phi_1` | weak adjoint `\Psi` (s06 uses `\Psi_1` for its component) |
| `\kappa_{\rm ii}`, `\kappa_0` | `\chi_{\rm ii}`, `\chi_0` | `\kappa=K/\Sigma` |
| `\delta_Q`, `\delta_1`, `\delta'`, `\hat\delta'` | `\vartheta_Q`, `\vartheta_1`, `\vartheta_N`, `\hat\vartheta_N` | `\delta=k-k_H`; a prime would read as d/dt |
| `(\delta U,\delta p,\delta w)` | `(V_U,V_p,V_w)` | `\delta=k-k_H` |
| `(a,b,c)` of the linear estimate | `(G_1,G_2,G_3)` | components of `e_s` (the triple is a solution of `eq:compound`, a different object) |
| `\mathcal N` (max functional) | `\mathcal G` | `\cN` is the neutral set |
| generic threshold `a`; `b=1.02` | `\xi`; `b_f` | components of `e_s`; `a_0`, `b_0` are the physical rates |
| points `a,b,c` in the corollary | `\mathbf a,\mathbf b,\mathbf c` | convention of s07 |
| `(s,t)`, `\Gamma_n`, `\gamma_n` in the same proof | `(r,v)`, `\Gamma_j`, `\mathbf g_j` | `t` is time, `n` the normal, `\gamma=1+<eta>` |
| dummy `\gamma` in the short-lobe lemma | instantiated at `\gamma_E` | `\gamma=1+<eta>` |
| dummy `T` in the Stage-2 lemma | instantiated at `T_-` | `T` is the least period |
| generic constants `c,d` (monotonicity criterion) | `q_0,q_1` | `c` component of `e_s`; `d` argument of `\Rem` (practice of s04) |
| `c_b` (log-coefficient, italic subscript) | `c_{\log}` | the source also has the function `c_{\rm b}` |
| `m=\min(\mu,1/\mu)` | `\underline\mu` | the source also uses `m(s)` for the mass function |

Kept as local symbols and listed in the paragraph "Local symbols" at the start of the section: `w`, the times
`s_w,s_0,T_\pm,\tau_1,\tau_b,\tau_c,\bar\tau_1` (stated to be measured in `t`, not in `\tau=\sqrt\mu t`), `u`, `v`
and its relatives, `z,z_1,z_{10},z^\sharp`, `r`, `s`, `\varsigma`, `\sigma` (two meanings, both stated), `\omega`,
`\omega_0`, `T_0(k),T_1(k),T'`, `P_d,P_w,P_1,p_k,w_d,w_k`, `M_+`, `m(s)`, `M_0` (the number of `lem:statebox`), `E_\theta`,
`Q`, `L,L_0,\hat L`, `a_1..a_4,a_c,A_1`, `c_0,c_1,c_T,c_p,c_\sharp`, `c_{\rm i},c_{\rm ii},c_{\rm b}`, `C_\tau,\bar C_\tau`,
`\nu,\alpha,\beta,\theta`, `\varphi,\varphi_2`, `\varepsilon_1` family, `q,g`, `V`, `\mathcal L`, `D^+` (Dini derivative),
`W^s(E)` (as `W^u(E)` in s07). The paragraph says explicitly that the subscripted `M_+`, `C_\tau`, `E_\theta`, `A_1`,
`T_0`, `T_1`, `T'`, `a_j`, `c_j` are numbers unrelated to `M`, `C`, `C_\pm`, `E`, `e_s`, `a_0`, `b_0`, that `u,v,z` are
not the coordinates of s04/s07, and what a prime means on each non-time quantity. This is a judgement call: these are
subscripted scalars that cannot be read as the global object (same practice as `M_0`, `C_Q` in s08 and `\kappa_0` in
s09). If a stricter rule is wanted they are simple token renames in `s10_largegain_port.py`.

## 5. Structural edits (S); everything else is the source text

1. Nine `\section` headings became `\subsection`s under `\section{Uniqueness at large gain}\label{sec:largegain}`; a
   short introduction (what is proved, five steps of the proof) and the paragraph "Local symbols" were added.
2. Statement: the opening lines refer to `eq:logsystem`, `eq:parameters`, `sec:model`, `eq:ratios`,
   `prop:multipliers` instead of v2.0.
3. Theorem (iii): **VERIFY.md remark 1** applied: the set is `{t: x(t)>c_0}`, "one interval per period, that is, a single
   interval modulo T" (with the stated time origin the set `{t in [0,T)}` wraps around `t=0`).
4. Fixed functions: **VERIFY.md remark 2** applied: one clause proving that `\Phi_\mu(u)` is increasing in `u`,
   `d/du \Phi_\mu = (1-e^{-u}) - \mu\Phi_\mu >= (1-e^{-u})e^{-\mu u} > 0`.
5. Imported facts: see §6. The lane's "v2.0" locations are replaced by labels; (F6) has a one-line derivation; (F7)
   has the hypothesis check; (F8) is new.
6. Uniqueness corollary: see §7. Title "Uniqueness; uses (F7)" -> "Uniqueness"; the statement adds "the cycle is the
   omega-limit set of every such orbit".
7. Reciprocal corollary: part (a) states the map `(mu,k) -> (mu',k') = (1/mu, k/mu^3)` explicitly (involution,
   preserves `rho`, `lambda'=lambda/mu`, `\tilde U(\tilde t)=U(\tilde t/\mu)`, `\tilde T=\mu T`, bijection of cycles, equal
   multipliers, positive domains not mapped onto each other unless `mu=1`); part (b) is the source statement. The
   source proof is kept in full, including the sentence that only cycles and multipliers are transported. Added to the
   proof: the map is an involution of `\cP`; "this proves (a) except for its last sentence"; the explicit computation
   behind "the positive domains do not correspond" (`{U'<\lambda, U''<\mu(\lambda-U')}` is mapped onto
   `{\tilde U_t<\lambda', \tilde U_{tt}<\lambda'-\tilde U_t}`, the domain at `(mu',k')` has the factor `\mu'`); the
   check `(mu,k) in \cP`; and the translation of `k>=1000`, `k>=1000 mu^3` into `rho>=rho_0(mu)`.
8. Lemma (identities) and the reciprocal proof: the scalar equation is cited as `eq:scalar` (formula kept).
9. Lemma (linear estimate): one sentence "This is the linear system `eq:compound`" after the system for
   `(G_1,G_2,G_3)`; the one-line display is broken into two lines (`gather*`), since with the longer symbols it was
   57 pt overfull.
10. Lemma (Stage 1): the four thresholds are displayed (the inline formula was 35 pt overfull).
11. Lemma (crest): the two references to `constants.py` became references to `lg:sec:constants`.
12. Constants: the table is a float with caption in its own subsection, with the sentence that every constant and
    scalar side condition is recomputed in high precision by a script in the bundle with receipt
    `evidence/v3/v3_largegain_constants.json`. The row `\hat H, \hat L, \delta', a_c, A_1` was split in two (5 pt
    overfull); `\hat\varepsilon_1(k_0)=0.06637` was added (the proof quotes it as 0.0664) and the entry `\delta'` is
    labelled `\hat\vartheta_N`, since the number 0.07772 is the value of the majorant.
13. Closing remark (`lg:sec:scope`): the source remark "Honest scope" was moved to the end and retitled; (a) now
    says "floating-point" and points to `sec:evidence`; (b) states `rho_0(mu)=1000/(mu(1+mu))` and `rho_0 ~ 1000/mu`;
    (c) unchanged. The source remark on where mu-uniformity is lost follows it unchanged.

## 6. Imported facts (F1)-(F7) of the source -> labels of this paper

| Fact | Content | Reference in the section |
|---|---|---|
| F1 | positive adjoint Floquet solution `\mathbf H'=B\mathbf H`; `m_t=exp((<eta>-mu)T)`, `0<m_s<m_t`, `m_s m_t=e^{-AT}` | `sec:adjoint`, `eq:adjointmatrices`, `eq:ratios`, `prop:multipliers`, `eq:multipliers`, `eq:exponent` |
| F2 | `\LamE>0`; exclusion of the cones `C_\pm` | `cor:cones`, `eq:LambdaE`, `eq:cones` |
| F3 | one crest and one trough, both simple zeros of `p` | `thm:geometry` |
| F4 | `T > 2 pi/sqrt(mu)` | `prop:crestmap`, `eq:determinant` |
| F5 | `x <= k M_0`, `M_0=(lambda+mu)(lambda+1)/mu`; `y>0`, `p<lambda` | `lem:statebox`, `eq:statebox` |
| F6 | `<x>=k`, `<p>=0`, `<p'^2>=mu<p^2>` | derived in one line from `lg:eq:sys`; last identity as in the proof of `prop:crestmap` (`eq:scalar`) |
| F7 | Smith's convergence theorem | `\cite[Theorem~1.1(c)]{Smith12}`, with the check `r_0>a_0+b_0 <=> k>k_H` and octant = `\Dom` (`eq:physical`, `eq:change`) |
| F8 (new) | parameter-wise criterion | `cor:criterion`(ii) |

Not used and therefore not listed: `eq:projective`, `eq:orthogonality`, `eq:crestidentities` (they enter only through
the proofs of the imported statements). `eq:compound` is cited in the proof of `lem:lglinear` for identification only.

## 7. Logic of the uniqueness corollary

`cor:criterion`(ii), as written in s07, states: if every cycle at `(mu,k) in \cP` has `lambda_c<0`, there is exactly one
cycle, it is hyperbolic and orbitally asymptotically stable, **and** it is the omega-limit set of every solution in
`\Dom` that does not converge to `E` (i.e. off the one-dimensional stable manifold of `E`). Its proof in s07 uses
`thm:attractor` (A1), (A6), (A8)(d), `lem:absorbing` and classical invariant-manifold theory; it does not use Smith's
convergence theorem. Hence `cor:criterion`(ii) suffices on its own for uniqueness and for the attraction of all
positive orbits off `W^s(E)`; Smith's theorem is not needed.

The section therefore proves `cor:lgunique` as follows: `k_H <= 2 < 1000 <= k` gives `(mu,k) in \cP`;
`thm:largegain`(iv) gives `lambda_c <= -0.17 mu < 0` on every cycle; `cor:criterion`(ii) gives all assertions. The
lane's own argument (basins of hyperbolic attracting cycles are open and disjoint, their union is `\Dom \ W^s(E)` by
Smith's theorem, and `\Dom \ W^s(E)` is connected by the measure-zero cone argument) is kept verbatim as a "second
proof, from (F7) in place of (F8)", which is independent of s07. For `mu>=1` both proofs are applied at `(mu,k)`
itself (`k_H <= 2 mu^2 < 1000 mu^3 <= k`).

## 8. Checks

`ROOT/v3/checks/s10_largegain_check.py` -> `s10_largegain_check.out`: 218 checks, 0 failures.

- A. The section file equals the output of the port script.
- B. `T` is a bijective relabelling: `T^{-1}(T(body)) == body` (whitespace-insensitive; 257 of 553 source lines carry
  a renamed token), `T` is idempotent on the ported text, the label map is injective, and no clashing token is left
  in the mathematics (bare `a b c d h H D \gamma`, `\ell \delta \kappa \Psi K_0 W_c X_\theta \mathcal N`, ...). This
  is the "diff against the source after normalising the renamed tokens": after `T` the only differences are the 34
  hunks of §5, printed in full in the output.
- C. Coverage of environments, displays, inline formulas and numeric literals (figures in §1).
- D. Labels and references.
- E. Every constant of the table and every numerical claim and scalar side condition of the proofs, recomputed with
  mpmath at 50 digits from the formulas **as written in the section's notation** (independent of the lane script):
  closed forms of `\Phi_1`, `c_i`, `c_ii` against their integrals; `omega_0`, `T_0`; Stage 1 and Stage 2 numbers;
  `\vartheta_1`, `T_1`, `z_{10}`, `\chi_0`, both endpoint values and a brute-force check that `f(z) > 0.3\chi_0 k` is
  impossible on `[z_1(k), z^\sharp]`; the crest majorants at `lambda=k_0` and all side conditions, the monotonicity of
  the majorants and `C_\tau <= 71.79+12.28 log(lambda/k_0)` on 601 values of `lambda` in `[1e3,1e15]`; `q(k_0)<0.0772`;
  the final margin; `omega >= omega_0(k)` and `gamma_E-mu<k^{1/3}+1` against true roots on a grid. Each printed number
  is checked to occur in the section and to be the rounding (or the stated one-sided bound) of the recomputed value.
- F. sympy: the system for `(G_1,G_2,G_3)` is `eq:compound`; `\mathcal L'=-\nu\mathcal L+(\beta x-\mu+\nu)G_1`;
  `c_0=(\mu-\nu)/\beta`; `U''-2\sigma U'+\zeta_E U=-\LamE`; `4\omega^2=3\gamma_E^2-2A\gamma_E-(1-\mu)^2`; factorisation
  of the characteristic polynomial; `M_0` in its three forms; `k/c_0`; the reciprocal map (involution, `rho'=rho`,
  `lambda'=lambda/mu`, `sigma'=sigma`, transformation of the scalar equation, the images of the positive domains,
  `rho_0`); the variational equation; `d/du \Phi_\mu`; the sign of `f'(z)`; the monotonicity criterion.

Lane scripts re-run **unchanged** after the renaming:
- `constants.py`: "ALL CONSTANT CHECKS PASSED for K0 = 1000.0", 55 checks (`ROOT/v3/checks/s10_largegain_lane_constants.out`;
  also executed inside the check script).
- `check_lemmas.py`: "ALL LEMMA CHECKS PASSED", 11 cycles, 376 inequality lines OK, 0 FAIL
  (`ROOT/v3/checks/s10_largegain_lane_check_lemmas.out`; run through `research-run`, identical to the lane's
  `check_lemmas.out` apart from the queue messages).
The lane directory was only read.

## 9. Compile test

`BUNDLE/paper/test_s10_largegain.tex`, pdflatex twice with the output directory in a scratch folder (no auxiliary
file was written in the bundle; the wrapper is kept): no LaTeX error, no overfull or underfull box, 11 pages.
Undefined references: only the 28 labels of other sections listed in the check output, and the citation `Smith12`
(the wrapper has no bibliography). No new bibliography entry is needed.

## 10. Not ported, and points for the frame and the evidence engineer

- Not ported: the title block and the "Status" paragraph of the lane file (they describe lane scripts, the
  cross-model review and the refereeing status; the paper states status in s01/s14).
- `lg:rem:scope`(a) quotes floating-point values: `log m_t = -13.98, -67.5, -375.7` at `mu=1`, `k=16, 60, 300`, and
  `T mu/lambda ~ 1.305`, `p_k/lambda ~ 1` at large gain. They come from the lane's `check_lemmas.out` (re-run here:
  -13.982, -67.506, -375.65; `T mu/lambda` = 1.3070 ... 1.3045 for `k>=1000`). Section 13 should carry a receipt for
  them, or the sentence should be reduced to the qualitative statement.
- The sentence in `lg:sec:constants` cites the receipt `evidence/v3/v3_largegain_constants.json`; the bundle script
  that writes it does not exist yet. Part E of `s10_largegain_check.py` is a complete recomputation in the paper's
  notation and can be used for it.
- The paragraph "Local symbols" is long (the proof uses more than a hundred local symbols); it is complete rather
  than short.
- The check script and the port script live in `ROOT/v3/checks/` (next to the s07 check), not in the bundle. If the
  section file is edited by hand, check A of `s10_largegain_check.py` will report the difference; the port script
  would overwrite hand edits and should then not be re-run.
- `L_0` here is `log(k/c_0)`; s08 uses a local `L_0` with another meaning inside `lem:statebox`. `M_0` is the same
  number in both sections.
