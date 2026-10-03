# Report — `s06_weakadjoint.tex` (Section 6: strong-stable vector, weak adjoint, twelve-event order, dual supporting planes, monotone crest map, rank-one identity)

Written by the auditing editor on 3 October 2026. The section text was drafted by a writer who was interrupted before
self-checks and before a report; it was then refereed by a second model (`ROOT/work/sol/review_s05_s06_response.md`) and
repaired by the orchestrator. This report records the comparison of the final text with the sources, the repairs made at the
audit, and the checks run.

ROOT = `/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03`. Sources: `ROOT/work/lanes/monotone/ds_theorem.tex` (cited
by label and line), `monotone/NOTES.md` (R1–R7) and `monotone/VERIFY.md`.

## 1. Statements included and their sources

| Label | Content | Source | Change relative to the source |
|---|---|---|---|
| §`wk:sec:signs` | event times `wk:eq:events`, sign table `wk:eq:signs`, arcs `I_+`, `I_-`, `theta=-p'/p` and `wk:eq:theta`, decay of `n.v` (`wk:eq:ndecay`) | ds_theorem §1, (F3), arcs `C`, `F`, eq. (theta) (l. 42–65) | arcs renamed `C,F -> I_+,I_-`; `wk:eq:ndecay` is new here (from the proof of Prop. `prop:factor`) and is used in `lem:strongstable`(d) and `prop:crestmonotone`(d) |
| `lem:strongstable` (a),(b) | `e_s=(a,-b,-c)`, `a,b,c>0`; Riccati system; `frak r>max(1,mu)`, `frak s>max(0,(1-mu)frak r)` | ds_theorem `lem:es` (l. 68–92); NOTES R1 | the conjugacy argument for `frak r>1` is written out as in VERIFY §1 ("Lemma 1") |
| `lem:strongstable` (c) | scalar equation `wk:eq:rscalar`, factorisation `wk:eq:rfactor`, `<frak r>=gamma` | NOTES R1 (scalar form, mean); R6 (`L3[phi]=L2'[phi'+r phi]`) | none |
| `lem:strongstable` (d) | `frak r` is the only periodic solution of `wk:eq:rscalar` | not in the lane | new statement with a complete proof (re-derived at the audit, see §5) |
| `wk:cor:reform` and the paragraph after it | `lambda_c<0 <=> <frak r> < A` | NOTES R1 | wording corrected per VERIFY §3.2 |
| mirror equation for `g=1+eta` | | NOTES R1 | none |
| `wk:lem:weakadjoint` | `Psi=e^{At} Gamma' x e_s`, adjoint equation, quasi-periodicity, chord identities `wk:eq:chord` | ds_theorem `lem:psi` (l. 94–113); NOTES R2 | none |
| `thm:twelve` | two simple zeros of each `Psi_j`, order `wk:eq:twelve`, sign of `q=p'+frak r p` | ds_theorem `thm:twelve` (l. 116–164); NOTES R3 | zeros renamed `tau_j^+-, t_2, t_1 -> t_j^+-, t_3^+, t_3^-`; the full derivative identities for `theta-frak r`, `g_2`, `g_1` are displayed (VERIFY §1, "Theorem 4") |
| `wk:lem:line` | supporting line at `t_*=t_3^-`; `wk:eq:linedata` | ds_theorem `lem:line` (l. 167–192); NOTES R4 | the three facts of `wk:eq:linedata` are collected in the statement |
| `thm:dualsupport` | three determinant inequalities; extension to entire solutions and nested cycles | ds_theorem `thm:DS` (l. 194–226); NOTES R5 | nestedness cited from Theorem `thm:geometry` (VERIFY §3.3) |
| `wk:rem:double` | two families of supporting planes; Wronskian statements | ds_theorem `cor:frame` (l. 228–235) | the weak Floquet vector is used only if `lambda_c != 0`; otherwise any vector of `n^perp` with `Psi.e_t>0` |
| `prop:crestmonotone` | first return on `Q={k e^U >= y > k}` exists, is `C^1`, maps `Q` into its interior, `DP>0`; Perron vector at a crest | NOTES R7 (l. 79–85); VERIFY §3.1 and §4 (i)–(iv) | see §2 |
| `prop:rankone` | rank-one window identity; four instances | ds_theorem `prop:rank` (l. 241–265); NOTES R6 (fourth instance, chord weight) | none |
| `wk:rem:secondorder` | Sturm–Liouville and Lagrange form, `wk:eq:ranktwo` | NOTES R6 | none |
| closing paragraph | what the identity gives and does not give; `alpha_*` orthogonality | ds_theorem Remark (l. 266–278); NOTES R6, §6 (O1) | none |

Not ported: NOTES R8 (first-order consequences of convexity of the backward flow; part (ii) is a sketch not written to
referee standard) and the second half of the remark in ds_theorem that rests on it; NOTES R4's side remark on tangent lines
at phases with `-p'/p >= mu` (superseded by Theorem `thm:gauges`); the numerical items N1–N6 and the refuted items X1–X6
(evidence and Section `sec:negative`).

## 2. VERIFY corrections and their status

- **VERIFY §3.1 and §4 (existence and `C^1` regularity of the first return on the whole region).** Applied in the proof of `prop:crestmonotone`: Step 1 = (i) (no zero of `p` on a solution converging to `E`), Step 2 = (ii) (no tangency), Step 3 = (iii) (existence of `t_a`, `t_b`), Step 4 = (iv) (transversality and regularity), Step 5 = the lane's wedge computation. Statement (a) asserts existence and `C^1` regularity on an open subset of the section containing `Q`.
- **VERIFY §3.2 (wording).** The text says that `<frak r> < A` is equivalent to `lambda_c<0` for the cycle under consideration, that it implies uniqueness if it holds on every cycle (Corollary `cor:criterion`), and that the converse is not claimed. The word "conjecture" does not occur in the section.
- **VERIFY §3.3 (citation of nestedness).** Cited from Theorem `thm:geometry`.
- VERIFY §3.4 concerns the lane's numerics, not the text.

## 3. Repairs made after the second-model review (checked at the audit)

- `lem:strongstable`(d), Sturm separation: the proof now uses the infinite sequence of consecutive simple zeros of the periodic function `p` and treats dependent solutions separately. Correct.
- `prop:crestmonotone`, Step 3, second half: both alternatives are written out. If `U -> +infinity`, then `x` increases to infinity, `y(t) >= (1-e^{-1}) x(t-1)` by the variation-of-constants formula and the monotonicity of `x`, and `p' <= k-y` forces a zero of `p` in finite time. If `U` is bounded, the orbit is bounded, the omega-limit set is `{E}`, and Step 1 gives a contradiction with `p(0)=0`. Each inequality was re-derived; correct and complete.
- `prop:crestmonotone`, Step 4: the zero pattern of nearby solutions is controlled on `(0,eps]` by the integral form of `p`, on the two compact intervals where the reference `p` is bounded away from zero, and on the two short intervals where `p'` has a fixed sign. This excludes additional zeros; the implicit function theorem then gives the `C^1` return time. Correct and complete.
- `wk:lem:line`: the inequality used for `U(t_*)<0` is now derived in the text from `Lambda_E>0` (Corollary `cor:cones`) by the substitution `y-k=-p'-mu p`. Correct.

## 4. Discrepancies found at the audit and fixed

1. **Sign in the Wronskian sentence after `wk:lem:weakadjoint`.** The text called `Psi_3=e^{At}(a'p-ap')` "the negative of the Wronskian of `p` and `a`". With the convention used in Remark `wk:rem:double` of the same section, `Wr(p,a)=pa'-p'a`, the two are equal: `Psi_3=e^{At}Wr(p,a)`. Corrected.
2. **Sentence before the proof of `wk:lem:line`.** "It needs only Lemma `lem:strongstable` and the sign table" omitted the two other inputs of the proof (the Jordan-curve statement of Theorem `thm:geometry` and Corollary `cor:cones`). Corrected.
3. **Closing paragraph.** "`L_3[phi]` must change sign exactly at the two tangency phases" claimed more than follows; replaced by "is nonnegative on the crest arc and nonpositive on the trough arc, so it can change sign only at the two tangency phases".
4. **Layout.** Three displays overfull by 50 pt, 12 pt and 88 pt (the arcs `I_+-`, `wk:eq:linedata`, `wk:eq:rankinstances`) were split; no overfull box remains.

No sign or index error was found in the renamed formulas (`r,s -> frak r, frak s`; `psi_p -> q`; `z -> Gamma`; `c(u) -> Rem`; `kappa -> beta`; `ell -> L_Gamma`).

## 5. Statements weakened or added

- Weakened: none.
- Added relative to the lane: `lem:strongstable`(d). Proof checked step by step: for a periodic solution `R`, `a~=exp(-int R)` gives a positive Floquet solution `v` of the variational equation; if its multiplier differed from `m_s`, `wk:eq:ndecay` and the periodicity of `n` would give `n.v=0`, so `a~` would solve the second-order equation `eq:tangent`, which `p` also solves; by Sturm separation `a~` would vanish somewhere. The multiplier `m_s` is simple because `m_s<min(1,m_t)`.

## 6. Formulas re-checked and how

- **By hand**, line by line: the sign table against `eq:six`; `wk:eq:theta`; `wk:eq:ndecay`; `wk:eq:abc`, the conjugate system for `c~`, `wk:eq:ricc`, `wk:eq:rscalar`, the expansion of `wk:eq:rfactor`, the computation in (d); the cross product `wk:eq:psi`, `wk:eq:psicomp`, `wk:eq:quasi`, `wk:eq:chord`; the three crossing-derivative identities and every interval of the proof of `thm:twelve`; `L_Gamma'`, `L_Gamma''(t_*)`, `L(0,0)`; `beta<0` and the two values of `B(t_*+jT)`; all five steps of `prop:crestmonotone` and part (d); `wk:eq:L3values`; the boundary term of `wk:eq:rankone` and the four instances; the Lagrange identity and `varpi q'=-Psi_2` at `t_*`; `alpha_*`.
- **Exactly (sympy)**: `ROOT/work/framework_checks/s06_weakadjoint_check.py`, Part A, items A01–A43, all pass. A01–A35 are the writer's. Added at the audit: A36 (Sturm–Liouville form), A37 (Lagrange identity), A38 (`Psi_3` as a Wronskian, `Wr(p,a)=-aq`), A39 (`det(v_f1,v_f2,v_f3)=-Wr(f1,f2,f3)`), A40 (`n.e_s`), A41 (`alpha_*`), A42 (the cross-product identity for a general matrix), A43 (sum and product of the two unstable roots).
- **Floating point** (same script, Part B; cycles at `(mu,rho)` = (1,1.5), (0.3,2), (3,1.3)): all pass. Writer's items B01–B10: `<frak r>-A=<eta>-mu`; `frak r>max(1,mu)`; the twelve-event order; signs of `Psi` at crest and trough; the supporting line; the three determinant inequalities on 200 x 200 pairs; quasi-periodicity of `varpi`; the rank-one identity for `phi=1, p', -U`; the Perron vector `(1,zeta(crest))`; `DP>0` and `x_1>y_1>k` at 12 points of `Q` including `x_0=y_0`. Added at the audit: B11 (fourth instance, chord shift `T/3`: window integral equals `(e^{-lambda_c T}-1) beta L` to 8 digits, `beta<0<L`) and B12 (`int Psi_3 x (p^2+alpha_*) = 0`, residual below 1e-9 of `int |Psi_3| x p^2`).

Output of the last run: `ROOT/work/framework_checks/s06_weakadjoint_check.out` ("ALL PASS", 85 items).
Part (d) of `lem:strongstable` and the existence part of `prop:crestmonotone` are not numerical statements; they rest on the proofs.

## 7. Check script

`s06_weakadjoint_check.py` was moved from `BUNDLE/paper/sections/` to `ROOT/work/framework_checks/` and extended there (A36–A43, B11, B12). Run:
`/Users/admin/.venvs/analysis-py313/bin/python -W ignore s06_weakadjoint_check.py > s06_weakadjoint_check.out` (single process, a few minutes).

## 8. Compile test and bibliography

`pdflatex -interaction=nonstopmode -halt-on-error test_s06_weakadjoint.tex` in `BUNDLE/paper`: no error, 11 pages, no overfull box; the undefined references are cross-section labels and the citation key. The auxiliary files (`.aux`, `.log`, `.out`, `.pdf`) were removed from `BUNDLE/paper`; the wrapper is kept. The only citation, `MPS90`, is defined in `BUNDLE/paper/bib/s01_intro.tex`; no bibliography fragment was added.
