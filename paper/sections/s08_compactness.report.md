# Report — `s08_compactness.tex` (Section 8: uniform compactness of cycle families)

Written by the auditing editor on 3 October 2026. The section text was drafted by a writer who was interrupted before
self-checks and before a report. This report records the comparison of the text with its source, the checks run and the
(small) repairs made.

ROOT = `/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03`. Source: the v2.0 manuscript
`ROOT/v2_0_source/wilhelm_heinrich_finite_parameter_v2_0/paper/manuscript.tex`, Section 6 (Theorem 6.1 and its proof, Steps 1–5,
lines 391–461). The orchestrator had already re-derived Lemma `lem:eqseparation` and Lemma `lem:speed` line by line; the audit
concentrated on Lemma `lem:statebox` and Theorem `thm:compact` and re-read the other two.

## 1. Statements included and their sources

| Label | Content | Source (v2.0 proof of Thm 6.1) | Change relative to the source |
|---|---|---|---|
| `lem:statebox`, `eq:statebox`, `eq:twosidedbox` | `X,Y,Z in [beta_0,M_0]` for every positive periodic solution of `eq:filters`, all `mu,k>0`; `M_Q`, `beta_Q` | Step 1 (l. 397–423) | every inequality written out (see §3); constants renamed `m,L,d -> mu_wedge, L_0, beta_0`; the tail bound is proved by an explicit integral estimate instead of the coupling remark; the upper bounds for `Y,Z` come from the maximum principle instead of the convolution representation |
| `cpt:eq:logbox` | the box in log coordinates | Step 1, last sentence | written as a display |
| preliminaries of §8.2 | spectrum of `J_E`, `v_E`, `n_E`, projection `P_u`, formula `cpt:eq:expformula`, quadratic form `V`, bounds `cpt:eq:Vbounds`, identity `cpt:eq:Vlinear`, constants `cpt:eq:constants` | Step 2 (l. 425–441) | all explicit; the v2.0 stable eigenvector `(1,-gamma_E,k/(1-gamma_E))` equals `v_E` because `k=gamma_E(gamma_E-1)(gamma_E-mu)` |
| `lem:eqseparation`, `eq:anglebound` | (i) `|r| <= C_Q |P_u r|`; (ii) `DV >= c_Q V` when `|P_u r| <= eps_Q`; (iii) `|r| >= d_Q` | Step 2 | v2.0 obtained the constant in `eq:anglebound` by a compactness-and-contradiction argument; here `C_Q=max_Q(1+|v_E|)` is explicit and follows from both cone exclusions of Corollary `cor:cones` |
| `lem:speed`, `eq:speedrelations`, `cpt:eq:speed` | uniform lower bound for the projected speed | Step 3 (l. 443–451) | constant `C_Q'` explicit |
| `thm:compact` (a)–(d) | state box, separation, period bounds, sequential compactness with least limiting period | statement of Thm 6.1 and Steps 4–5 (l. 453–459) | `T_Q` explicit in terms of `beta_Q, M_Q, lambda_Q, sigma_Q`; the proof that the limit starts at its crest and that its period is least is written out |
| closing paragraph | what is and is not claimed | l. 461 | the constants are now explicit in principle, not existential |

Nothing in the source was dropped. No VERIFY file exists for this section (it is a port of v2.0, not of a lane).

## 2. Discrepancies found and fixed

- No mathematical discrepancy was found. Every hypothesis of the v2.0 theorem is kept (`Q` compact in `cP` for the theorem; `lem:statebox` alone holds for all `mu,k>0`, as in v2.0).
- Layout: the display of `M_0, mu_wedge, L_0, beta_0` (overfull by 21 pt) was split into two lines, and the formula for `DV(r)h` (inline, overfull by 35 pt) was displayed; the sentence on the lower bound for `V` now states the integral it comes from. No overfull box remains.

Statements weakened: none.

## 3. Lemma `lem:statebox`: the inequalities, as checked

1. `R=Z/X`, `S=Y/X`: `R'=1-R(1+lambda(1-Y))`, `S'=mu R-S(mu+lambda(1-Y))` (sympy A01, A02). At a minimum of `R`, `R=1/(1+lambda(1-Y))>1/(lambda+1)` because `0<1+lambda(1-Y)<1+lambda`. At a minimum of `S`, `S=mu R/(mu+lambda(1-Y))>mu/((lambda+1)(lambda+mu))=1/M_0`.
2. At a maximum of `X`, `Y=1`, so `max X=1/S<M_0`; at extrema of `Z` and `Y`, `Z=X` and `Y=Z`, which gives the upper bounds and `min Y>=min Z>=min X`.
3. Tail: `I(s)=(e^{-s}-e^{-mu s})/(mu-1)` (`s e^{-s}` for `mu=1`) (A04, A05); the two expressions for `theta_mu` agree (A06); `theta_mu` is the second-filter response and the survival function of `Exp(1)+Exp(mu)` (A07, A08); `I(s)<=s e^{-mu_wedge s}`; `theta_mu<=(1+mu_wedge s)e^{-mu_wedge s}<=2e^{-mu_wedge s/2}` (A09 and the grid test B01, seven values of `mu` including `mu` near 1).
4. Excursion below `1/2`: `<X>=1` gives `max X>=1`; on the component `(t_0,t_1)` of `{X<1/2}` containing the minimiser, `Z-1/2<=(M_0-1/2)e^{-s}` and `Y-1/2<=(Y(t_0)-1/2)e^{-mu s}+mu(M_0-1/2)I(s)<=(M_0-1/2)theta_mu(s)` (A11; variation of constants with a positive kernel); for `s>=L_0`, `theta_mu<=1/(4M_0)` (A10), so `Y<3/4` and `X'>0`; hence `t_m-t_0<L_0`.
5. Lower bound: `(log X)'>=-lambda(M_0-1)` gives `log X(t_m)>=log(1/2)-lambda(M_0-1)L_0=log beta_0`; the case `min X>=1/2` is covered because `beta_0<=1/2`.

Theorem `thm:compact`: the arclength bound (`sqrt(p^2+p'^2)<=|U'|+|p'|`, total variation of a function with one maximum and one minimum per period), the period bound, the passage to the limit (uniform box, continuous dependence, `z_*(T_*)=z_*^0`, nonconstancy from (b)), the identification of the limit's crest (`U_*(0)>=0`, `U != 0` at zeros of `p`) and the least-period argument (a second downward zero of `p_*` would persist under `C^1` convergence) were each re-derived; no gap.

## 4. Formulas re-checked and how

- **By hand**: all of the above, the preliminaries of §8.2 (`chi_E(-gamma)=k-f_mu(gamma-1)`; nonreality of the unstable pair; `v_E`, `n_E`, `n_E.v_E=f_mu'(eta_E)`; `P_u`; `cpt:eq:expformula`; the three bounds on `V`; `DV(r)J_E r=|P_u r|^2`), parts (i)–(iii) of `lem:eqseparation`, and `lem:speed`.
- **Exactly (sympy)**: `ROOT/work/framework_checks/s08_compactness_check.py`, Part A, items A01–A25, all pass (ratio equations, tail identities, change of variables, characteristic polynomial and its factorisation, eigenvector relations, projection identities, `(J_E-alpha_E)^2P_u=-omega_E^2P_u`, the exponential formula as the solution of `M'=-J_E M`, `r'=J_E r+N(r)`, the integral form of `Rem`, the speed relations, the crest identity).
- **Floating point** (same script, Part B; whlib cycles at `(mu,rho)` = (0.3,1.5), (1,2), (3,1.5), (0.1,1.1)), all pass:
  ratio bounds (e.g. at (1,2): `min Z/X=0.2124>0.2000`, `min Y/X=0.0522>0.0400`); upper box and the orderings of extrema; means equal to 1; the three excursion inequalities pointwise and `t_m-t_0<L_0` (e.g. 2.62 < 10.60); `log min X>=log beta_0`; spectrum of `J_E`; `cpt:eq:expformula` against the matrix exponential; `V` as the solution of the Lyapunov equation `J_E^T G+GJ_E=P_u^T P_u` against the defining integral; the bounds `cpt:eq:Vbounds`, the identity `cpt:eq:Vlinear` and the remainder bound of part (ii) on 2000 random vectors; `Lambda_E>0`, both cone exclusions, `c_r<|P_u r|` and `eq:anglebound` on every cycle point (largest `|r|/|P_u r|` between 1.06 and 1.36, against `1+|v_E|` between 2.9 and 8.0); `cpt:eq:Vdot` against finite differences; the separation constants for the one-point set `Q={(mu,k)}`; the speed relations (residual below 2e-11); the arclength inequality and the period bounds; one downward zero of `p` per period.

Output: `ROOT/work/framework_checks/s08_compactness_check.out` ("ALL PASS", 78 items).

The numerical checks also show how far the constants are from sharp, as the closing paragraph of the section says: for the
one-point sets, `log beta_0` is between −1018 and −3684 while `log min X` is between −13 and −2.8; `d` is between 4e-9 and
2e-5 while the smallest `|r|` on the cycle is between 0.53 and 3.6; the constant `C_Q'^2` is of order 10^887 to 10^3203. No
cycle point had `|P_u r|<=eps`, so part (ii) of `lem:eqseparation` was exercised only through its general bounds on random
vectors.

## 5. Check script

The section came without a check script. `s08_compactness_check.py` was written at the audit and placed in
`ROOT/work/framework_checks/` (not in the paper's sections directory). Run:
`/Users/admin/.venvs/analysis-py313/bin/python -W ignore s08_compactness_check.py > s08_compactness_check.out` (single process, about one minute).

## 6. Compile test and bibliography

`pdflatex -interaction=nonstopmode -halt-on-error test_s08_compactness.tex` in `BUNDLE/paper`: no error, 7 pages, no overfull box; the undefined references are cross-section labels. The auxiliary files (`.aux`, `.log`, `.out`, `.pdf`) were removed from `BUNDLE/paper`; the wrapper is kept. The section has no citation; no bibliography fragment was added.
