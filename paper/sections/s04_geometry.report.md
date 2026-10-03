# Report: s04_geometry.tex

Source: `ROOT/v2_0_source/wilhelm_heinrich_finite_parameter_v2_0/paper/manuscript.tex`, Section 4, lines 258-324. Ported verbatim apart from the changes listed below. No proof step was altered or removed.

## Statements included

| Label | Content | v2.0 location |
|---|---|---|
| `thm:geometry`, `eq:crestidentities` | regular injective clockwise projection to (U,p) enclosing the origin; one nondegenerate maximum and minimum of U and of p; projections of distinct cycles disjoint and nested; with proof | l. 259-274 |
| `thm:six`, `eq:six` | each of X, Z, Y has one maximum and one minimum per period, all nondegenerate, in the cyclic order X_min < Z_min < Y_min < X_max < Z_max < Y_max; with proof | l. 276-292 |
| `cor:orderedcrests`, `eq:crestorder` | zeta_i (U_o-U_i) < y_o-y_i < zeta_o (U_o-U_i) at crests of nested cycles; trough statement; with proof | l. 294-303 |
| `prop:crestmap`, `eq:determinant`, `eq:compound` | crest-return derivative M in (U,y) is entrywise positive, eigenvalues m_s, m_t, det M = e^{-AT} < e^{-4 pi}, T > 2 pi / sqrt(mu); similarity to M_{(U,z)}; with proof | l. 305-323 |
| (remark) | small determinant does not exclude m_t >= 1 | l. 324 |

All 10 exported labels of PLAN.md are present. Internal label added: `geo:eq:similarity`.

## Changes relative to v2.0

1. Proposition `prop:crestmap` now says "expressed in the section coordinates (U,y)" and adds the similarity `M_{(U,z)} = diag(1,1/k) M diag(1,k) = [[M_11, k M_12],[M_21/k, M_22]]` for z = y/k, with the statement that eigenvalues, determinant and sign pattern agree, while entries and componentwise inequalities on test vectors depend on the coordinates (a test vector xi for M corresponds to diag(1,1/k) xi). One paragraph was appended to the proof for this added claim.
2. Section variations are written `(\xi_U, 0, \xi_y)` instead of `(\delta U, 0, \delta y)`, because PLAN.md fixes `\delta = k - k_H`; `(y-k)\delta U` would be ambiguous.
3. In the Wirtinger equality case the constants `c, d` were renamed `c_0, c_1` (PLAN.md fixes a, b, c as components of the strong-stable vector, and d is the argument of `\Rem`).
4. Title of Theorem `thm:six`: "Chemical phase order, without swapping filters" became "Chemical phase order in the original labelling". The statement is unchanged.
5. A short opening paragraph lists the local symbols (subscripts o and i, xi_U, xi_y, w, z).
6. No scope sentence occurred in this section; nothing else was removed.

## VERIFY.md corrections

No lane VERIFY.md records a correction against v2.0 Section 4. `largegain/VERIFY.md` (F3, F4) read Theorem 4.1 in full and re-derived the Wirtinger step ("No gap found"); `monotone/VERIFY.md` (section 1) audited Theorems 4.1, 4.2 and Proposition 4.4 as imported; `dulac/VERIFY.md` (section 3.2(e)) confirmed the wedge sign convention of `eq:compound`. `monotone/VERIFY.md` (statement 5) notes that Proposition 4.4 is local to a cycle crest; the ported statement keeps the word "local", so no change was required.

## Formulas re-checked

Script `ROOT/work/framework_checks/check_s02_s04_appA.py`, output `check_s02_s04_appA.out` (exit status 0). Section 4 items, exact sympy identities (11 PASS):

- crest identities x - k = (1+eta)(y-k), x - y = eta (y-k) from O = 0 at p = 0;
- y'' = x p - y'; Y'' - zeta p/lambda = -[(eta-mu) p' + O]/lambda; y' = zeta p - eta p' - O;
- w' = C w for the sign-conjugated wedge of the tangent and a variation, and the initial wedge ((y-k) xi_U, (x-y) xi_U, (y-k) xi_y) at a crest (after the rename to xi);
- entries, determinant and trace of diag(1,1/k) M diag(1,k);
- 2 pi (sqrt(mu) + 1/sqrt(mu)) = 2 pi A / sqrt(mu); U''' + mu U' = 0 for U = c_0 + c_1 cos(sqrt(mu) t + theta) (after the rename).

Floating-point spot check with `whlib` at (mu, rho) = (0.1, 1.5), (1, 2), (3, 1.3): M and M_{(U,z)} entrywise positive; eigenvalues of M_{(U,z)} equal those of M; det M = exp(-AT) within an absolute tolerance (the determinant is below rounding at mu = 0.1, as the script states); exp(-AT) < exp(-4 pi); T > 2 pi / sqrt(mu); <p'^2> = mu <p^2>.

## Not ported

Nothing was left out.

## Compile test

`BUNDLE/paper/test_s04_geometry.tex`: pdflatex, two passes, no LaTeX errors, no overfull or underfull boxes, 3 pages; undefined references only to labels owned by s01, s02 and s03. Auxiliary files were written to a scratch output directory, none in the bundle.
