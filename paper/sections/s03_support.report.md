# Report: s03_support.tex

Source: `ROOT/v2_0_source/wilhelm_heinrich_finite_parameter_v2_0/paper/manuscript.tex`, Section 3, lines 200-256. Ported verbatim apart from the changes listed below. No proof step was altered or removed.

## Statements included

| Label | Content | v2.0 location |
|---|---|---|
| `thm:support`, `eq:delta`, `eq:supportidentity` | base = a cycle or E, second solution = distinct log-bounded entire solution; Delta' + (1+eta) Delta = -x_1 Rem(U_2-U_1); Delta < 0; with proof | l. 201-221 |
| (remark) | hypothesis holds for cycles and E, not assumed for forward trajectories; pointer to `thm:twoorbit` / `sec:attractor` | l. 222 plus one new sentence |
| `cor:cones`, `eq:positiveK`, `eq:LambdaE`, `eq:equilibriumreverse`, `eq:cones` | K > 0 with K' + (1+eta)K = xU - x + k = k psi(U); gamma_E, eta_E, zeta_E; Lambda_E > 0 with Lambda_E' = -gamma_E Lambda_E + k Rem(U); exclusion of the closed cones C_+, C_-; with proof | l. 224-245 |
| `lem:curvature`, `eq:curvature` | W = zeta' p - eta' p', W' + (1+eta) W = x p^2, W > 0; with proof | l. 247-256 |

All 11 exported labels of PLAN.md are present. No internal label was needed.

## Changes relative to v2.0

1. `\Rem(d) := e^d - 1 - d` is defined in the statement of Theorem `thm:support`, and `eq:supportidentity` is written with it. The proof's remainder vector is written both in the v2.0 form and as `(0,0,x_1 Rem(U_2-U_1))`.
2. In the integral representation the dummy variable `\tau` was renamed `s` (PLAN.md reserves `\tau` for the frequency-one time). "Strict convexity makes the bracket nonnegative. If it vanished ..." became "Strict convexity makes Rem nonnegative, with Rem(d) = 0 only at d = 0. If the integrand vanished ...": same step, reworded because the bracket is now `\Rem`.
3. `\psiU(U) := U e^U - e^U + 1` is defined after `eq:positiveK`, which now reads `K' + (1+eta)K = xU - x + k = k psi(U)`.
4. Corollary `cor:cones`: the linear gauge `\LamE := (y-k) - zeta_E U + eta_E p` is introduced with label `eq:LambdaE`, together with `\LamE > 0` and `\LamE' = -gamma_E \LamE + k Rem(U)` along every cycle. `eq:equilibriumreverse` keeps the v2.0 inequality `zeta_E U - eta_E p - y + k < 0`, now stated as the equivalent form. The proof keeps every v2.0 sentence and adds: Delta = -K and x Rem(-U) = k psi(U) for the first comparison; the explicit eigenvector equation `B (zeta_E, eta_E, 1)^T = gamma_E (zeta_E, eta_E, 1)^T` at x = k; Delta = -Lambda_E, x_1 = k, 1 + eta_E = gamma_E for the second comparison.
5. The remark after Theorem `thm:support`: "It is not silently assumed for every forward trajectory" became "It is not assumed for arbitrary forward trajectories", and one sentence was added pointing to Theorem `thm:twoorbit` in Section `sec:attractor` for arbitrary pairs of distinct log-bounded entire orbits. The theorem itself is stated exactly as in v2.0 (base = a cycle or E).
6. Lemma `lem:curvature`, last proof sentence: the v2.0 reference to "the exact-algebra check" (a script) was replaced by the mathematical statement that the residual of the identity off O = 0 is -eta' O, so the identity holds along cycles and not for independent values of (U,p,y,eta,zeta).
7. A short opening paragraph lists the local symbols (z_1, z_2, x_j, J_1, Delta, r). The solution names `z_1, z_2` are kept from v2.0; `z = y/k` occurs only inside Proposition `prop:crestmap` of Section `sec:geometry`.

## VERIFY.md corrections

No lane VERIFY.md records a correction against v2.0 Section 3. `largegain/VERIFY.md` (F2) re-derived `Delta' + gamma_E Delta = -k(e^U - 1 - U)` and `Lambda_E = -Delta > 0`; `dulac/VERIFY.md` (section 3.1) checked the hypotheses of Theorem 3.1 and confirmed that its proof uses neither k > k_H nor stability. These agree with the text as ported.

## Formulas re-checked

Script `ROOT/work/framework_checks/check_s02_s04_appA.py`, output `check_s02_s04_appA.out` (exit status 0). Section 3 items, exact sympy identities (12 PASS):

- Delta' + (1+eta) Delta + x_1 Rem(U_2-U_1) = 0 along two solutions; x_2 - x_1 - x_1 (U_2-U_1) = x_1 Rem(U_2-U_1);
- K' + (1+eta) K = xU - x + k; xU - x + k = k psi(U); x Rem(-U) = k psi(U); K = -Delta with second solution E;
- equilibrium eigenvector equation with k = gamma_E (gamma_E-1)(gamma_E-mu); A(A-1)(A-mu) = mu(1+mu) (so gamma_E > A exactly when k > k_H);
- Lambda_E' = -gamma_E Lambda_E + k Rem(U); Lambda_E = -(zeta_E U - eta_E p - y + k); Lambda_E = -Delta with base E;
- W' + (1+eta) W - x p^2 = -eta' O.

Floating-point spot check with `whlib` at (mu, rho) = (0.1, 1.5), (1, 2), (3, 1.3): K > 0, W > 0, Lambda_E > 0 at all 4001 samples of each cycle; no sample in C_+ or C_-.

## Not ported

Nothing was left out.

## Compile test

`BUNDLE/paper/test_s03_support.tex`: pdflatex, two passes, no LaTeX errors, no overfull or underfull boxes, 2 pages; undefined references only to labels owned by s01, s02 and s07. Auxiliary files were written to a scratch output directory, none in the bundle.
