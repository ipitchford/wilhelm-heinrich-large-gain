# Report: appA_hopf.tex

Source: `ROOT/v2_0_source/wilhelm_heinrich_finite_parameter_v2_0/paper/manuscript.tex`, Appendix A, lines 633-670. Ported verbatim apart from the changes listed below. No computation step was altered or removed.

## Content included

| Item | Content | v2.0 location |
|---|---|---|
| `app:hopf` | section anchor | l. 633 |
| `hopf:eq:scaled` | U_ttt + sigma U_tt + U_t + sigma rho (e^U - 1) = 0 in tau = sqrt(mu) t, sigma = (1+mu)/sqrt(mu) >= 2, rho = k/k_H | l. 634-642 |
| (display) | characteristic polynomial (nu+sigma)(nu^2+1) at rho = 1; d Re nu / d rho = sigma / (2(1+sigma^2)) > 0 | l. 643-647 |
| (display) | J_0, q = (1,i,-1)^T, left eigenvector (i sigma, sigma+i, 1)/(-2+2i sigma) | l. 648-653 |
| (display) | bilinear and trilinear forms; formula for G_21 | l. 654-663 |
| `hopf:eq:l1` | G_21 = sigma(sigma+3i)/[3(sigma+2i)(-1+i sigma)], l_1 = Re G_21 / 2 = -sigma/[(sigma^2+1)(sigma^2+4)] < 0 | l. 664-669 |
| (text) | supercritical branch; local asymptotic stability of the critical equilibrium | l. 670 |

Exported label `app:hopf` present. Internal labels added: `hopf:eq:scaled`, `hopf:eq:l1`. The file does not issue `\appendix`; `manuscript.tex` does so before the input.

## Changes relative to v2.0

Symbol renames forced by the global notation of PLAN.md (each formula re-checked, see below):

1. Left eigenvector `\ell` became `\chi` (PLAN.md fixes `\ell = \lambda - p`).
2. Multilinear forms `B`, `C` became `\mathcal B`, `\mathcal C` (B is the adjoint matrix of `eq:adjointmatrices`, C the compound matrix of `eq:compound`).
3. Spectral variable `z` became `\nu`, and the generic arguments `(v,w,z)` of the forms became `(\xi,\vartheta,\varsigma)` (z = y/k and the planar coordinates (u,v) are global).
4. `\rho = k/(\mu(1+\mu))` is written `\rho = k/\kH = k/(\mu(1+\mu))`.

Other changes:

5. A short opening paragraph lists the local symbols and states that subscripts tau are derivatives in tau = sqrt(mu) t.
6. "For completeness, use frequency-one time ..." became "Use the time tau ...". "The first-order Jacobian is" became "In the variables (U, U_tau, U_tautau) the first-order Jacobian at rho = 1 is" (clarification of what v2.0 computes).
7. Scope phrase removed: the last sentence "This recovers the established local result in a fixed convention, not a new claim that Hopf criticality was previously unknown" became "This recovers the established local result (see [WH95, WH96] and [BH23, Theorem 13]) in a fixed convention". The citations use existing keys only; `literature/NOTES.md` records Theorem 13 of Boros and Hofbauer (full text read) as stating the supercritical Hopf bifurcation and stability through the critical value, and v2.0 Section 1.3 attributes the local Hopf analysis to Wilhelm and Heinrich.

## VERIFY.md corrections

No lane VERIFY.md records a correction against v2.0 Appendix A. No change was required.

## Formulas re-checked

Script `ROOT/work/framework_checks/check_s02_s04_appA.py`, output `check_s02_s04_appA.out` (exit status 0). Appendix items, exact sympy identities with symbolic sigma > 0 (12 PASS):

- sigma rho = k / mu^{3/2} and the coefficient of U_tautau after rescaling `eq:scalar`; sigma - 2 = (1 - sqrt(mu))^2 / sqrt(mu);
- characteristic polynomial at rho = 1; d Re nu / d rho at nu = i, rho = 1;
- J_0 q = i q; chi J_0 = i chi; chi q = 1;
- J_0^{-1} B(q, qbar) = (1,0,0)^T; first component of (2iI - J_0)^{-1} B(q,q) = sigma/[3(sigma+2i)];
- G_21 and l_1 = Re G_21 / 2 = -sigma/[(sigma^2+1)(sigma^2+4)].

## Not ported

Nothing was left out.

## Compile test

`BUNDLE/paper/test_appA_hopf.tex`: pdflatex, two passes, no LaTeX errors, no overfull or underfull boxes, 1 page; undefined: `eq:scalar` (owned by s01) and the three citation keys (bibliography merged later). Auxiliary files were written to a scratch output directory, none in the bundle.
