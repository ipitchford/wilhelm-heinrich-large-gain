# Report: s02_adjoint.tex

Source: `ROOT/v2_0_source/wilhelm_heinrich_finite_parameter_v2_0/paper/manuscript.tex`, Section 2, lines 115-198. Ported verbatim apart from the changes listed below. No proof step was altered or removed.

## Statements included

| Label | Content | v2.0 location |
|---|---|---|
| `eq:adjointmatrices` | J, B = -S J^T S, S = diag(1,-1,-1) | l. 118-122 |
| (text) | strict positivity of H' = BH, Perron eigenvalue, positive Floquet solution | l. 123 |
| `eq:ratios` | eta = H_2/H_3 > 0, zeta = H_1/H_3 > 0 | l. 124-126 |
| `eq:projective` | eta' = zeta + (mu-1) eta - eta^2, zeta' = x - (1+eta) zeta | l. 128-131 |
| `eq:etaminbasic` | eta > max(0, mu-1) | l. 133-135 |
| (text) | gamma = 1 + <eta> > 1 | l. 136 |
| `prop:multipliers`, `eq:multipliers`, `eq:exponent` | 0 < m_s < m_t, m_s < 1, m_s m_t = e^{-AT}; m_s = e^{-(1+<eta>)T}, m_t = e^{(<eta>-mu)T}, lambda_c = <eta> - mu; with proof | l. 138-153 |
| (remark) | <eta> = mu leaves nonlinear stability undecided; ambient Jordan block not excluded | l. 154 |
| `prop:factor`, `eq:orthogonality`, `eq:tangent`, `eq:factorisation` | zeta p - eta p' - (x-y) = 0; p'' = (eta-mu)p' - zeta p; operator factorisation; with proof | l. 157-177 |
| `lem:projectivebounds` | eta_-, eta_+, zeta_-, zeta_+ from x_- <= x <= x_+; with proof | l. 179-191 |
| `eq:gainprimitive` | (zeta + eta + eta^2/2)' = x - f_mu(eta); k = <f_mu(eta)>; eta nonconstant | l. 193-198 |

All 14 exported labels of PLAN.md are present. Internal labels added: `adj:sec:projective`, `adj:sec:factor` (subsection anchors).

## Changes relative to v2.0

1. The positive adjoint solution is written `\mathbf H=(H_1,H_2,H_3)^T` (bold) instead of `H`. Reason: PLAN.md reserves the scalar `H` for the convex extension of the graph function `h`. Components `H_1,H_2,H_3` are unchanged, which matches the usage of the large-gain proof (`\mathbf H`, `H_i`).
2. A two-sentence opening paragraph lists the symbols introduced here (J, B, S, bold H, n, O).
3. "A useful exact primitive is" became "An exact primitive is" (no evaluative adjective).
4. Last sentence: "does not yet decide the sign" became "by itself it does not determine the sign" (mathematical caveat kept, status wording removed).
5. No scope sentence occurred in this section; nothing else was removed.

## VERIFY.md corrections

None of the six lane VERIFY.md files records a correction against v2.0 Section 2. `largegain/VERIFY.md` (F1) and `monotone/VERIFY.md` (section 1) re-derived Proposition 2.1 and Proposition 2.2 and confirmed them; `dulac/VERIFY.md` (section 3.1) re-derived `n' + J^T n + (1+eta) n = 0`. No change was required.

## Formulas re-checked

Script `ROOT/work/framework_checks/check_s02_s04_appA.py`, output `check_s02_s04_appA.out` (run with `/Users/admin/.venvs/analysis-py313/bin/python`, exit status 0, "ALL CHECKS PASSED"). Section 2 items, all exact sympy identities (12 PASS):

- B = -S J^T S; eta', zeta' and H_3' = (1+eta) H_3 from H' = BH (after the rename to bold H);
- n' + J^T n + (1+eta) n = 0; O' + (1+eta) O = 0 along the flow; p'' = -mu p' - x + y; p'' - [(eta-mu) p' - zeta p] = O;
- the factorisation (D+1+eta)[D^2+(mu-eta)D+zeta] = D^3 + A D^2 + mu D + x on a test function;
- the gain primitive; the two quadratic roots eta_+, eta_- of Lemma `lem:projectivebounds`.

Floating-point spot check with `whlib` at (mu, rho) = (0.1, 1.5), (1, 2), (3, 1.3): m_t = exp((<eta>-mu)T) to relative 1e-6; m_s = exp(-(1+<eta>)T) at the two cases where exp(-AT) is above rounding (one case skipped and reported as SKIP); <x> = k = <f_mu(eta)>; orthogonality residual below 1e-6 relative.

## Not ported

Nothing was left out.

## Compile test

`BUNDLE/paper/test_s02_adjoint.tex`: pdflatex, two passes, no LaTeX errors, no overfull or underfull boxes, 2 pages; only `eq:logsystem` undefined (owned by s01). Auxiliary files were written to a scratch output directory, none in the bundle.
