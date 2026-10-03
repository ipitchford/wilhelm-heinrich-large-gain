# Cartography tables: files and columns

These tables were computed during the research phase with a compiled integrator
(`code/v3/lib/cart.py`: numba DOP853 with coefficients imported from scipy, default relative
tolerance 1e-12, events located by Newton iteration on the step length). They are shipped
unchanged. Everything is ordinary floating point; nothing here is an interval enclosure, and a
search that finds no second cycle does not exclude one.

`code/v3/v3_cartography.py` checks the statements read off these tables, recomputes 51 grid
cycles with the bundle's own integrator (`code/v3/cyclekit.py`, scipy, logarithmic variables) and
re-runs 24 of the 240 uniqueness hunts. Its receipt is `evidence/v3/v3_cartography.json`.

Notation: `rho = k/(mu(1+mu))`, `U = log(x/k)`, crest = point of a cycle with `p = 0`, `U > 0`;
`m_t` is the weak transverse multiplier, `lambda_c = <eta> - mu`, `log m_t = lambda_c T`.

## `grid_table.csv` (300 rows) and `branch_table.csv` (827 rows)

Grid: `mu` in {0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.85, 1} times
`rho - 1 = logspace(1e-3, 199, 25)`. Each grid cycle was located by amplitude-fixed shooting
(crest amplitude prescribed, unknowns `y` at the crest and `k`), continued in amplitude from the
Hopf point, followed by a secant iteration on the amplitude to reach the target `rho`.
`branch_table.csv` holds the amplitude branches themselves (same columns, no `secant_iters`).

| column | meaning |
|---|---|
| `mu`, `rho`, `k` | parameters |
| `U_max`, `U_min` | crest and trough values of `U` on the cycle |
| `y0_over_k` | `y/k` at the crest |
| `T`, `T_scaled` | least period in the time of the scaled system; `T sqrt(mu)/(2 pi)` |
| `log_mt` | `log m_t = T(<eta> - mu)` from the quadrature of the periodic adjoint ratio |
| `log_mt_M` | logarithm of the larger eigenvalue of the 2x2 crest-return derivative (blank when `m_t < 1e-7`) |
| `lambda_c`, `lambda_c_T` | `<eta> - mu`; `lambda_c T` (equal to `log_mt`) |
| `mean_eta`, `mean_U`, `meanU_plus_log_rho` | `<eta>`, `<U>`, `<U> + log rho` |
| `q2`, `q3` | `<(eta-mu)^2>`, `<(eta-mu)^3>` |
| `moment_identity_residual` | `[(A^2+mu)<q> + 2A<q^2> + <q^3> - (k-k_H)]/(k-k_H)`, `q = eta - mu`, `A = 1+mu` |
| `mean_x_rel_err` | `<x>/k - 1` |
| `newton_residual` | crest-map residual of the located cycle (`U` absolute, `y` relative to `k`) |
| `closure` | mismatch after integrating one period with the profile system |
| `ratio_defect`, `settle_periods` | periodicity defect of `(eta, zeta, K, W)`; number of periods used to reach the periodic ratios |
| `eta_crest`, `zeta_crest` | adjoint ratios at the crest |
| `secant_iters` | secant steps on the amplitude (grid only) |

## `hunt_results.jsonl` (240 lines), `hunt_results_supplement.jsonl`, `hunt_results_pairs.npz`, `hunt_summary.txt`

One JSON object per random parameter point: `mu` log-uniform in [1e-3, 1e3], `rho - 1`
log-uniform in [1e-3, 299]. The random stream is consumed inside time-boxed loops, so the point
sequence is not reproducible from the seed; the recomputation reads the points from this file.
Line 110 is a failure of the integrator in direct coordinates (stiff trough at `mu = 141.65`,
`rho = 195.66`); `hunt_results_supplement.jsonl` holds the same point computed at the reciprocal
rate `1/mu`. `hunt_summary.txt` gives aggregate counts.

Keys: `mu`, `rho`, `k`, `U_star`, `y_star_over_k`, `T`, `log_mt`, `lambda_c`, `m_t`,
`branch_monotone` (`rho` increasing along the amplitude branch up to the target),
`primary_residual`.
Displacement scan of the crest return map on the attractor: `scan_seeds`, `scan_pairs`,
`scan_wu_pairs` (pairs on the unstable manifold of `E`), `scan_fail`, `scan_wrong_sign` (pairs
with `sign(U_{n+1} - U_n)` different from `sign(U_star - U_n)`), `scan_U_over_Ustar_range`,
`scan_r_monotone_violations`, `scan_r_max_increase`, `scan_r_slope_max`, `scan_r_at_E`, where
`r(U) = log(phi(U)/U)` and `phi` is the crest return map in the coordinate `U`.
Long runs: `long` (list of initial condition tag, number of returns, final gap to the primary
crest, status), `long_same`, `long_different`, `long_other`, `long_heading_bad`.
Newton searches from perturbed guesses: `newton` = {`same`, `equilibrium`, `failed`, `new`}; `new`
lists fixed points distinct from the primary one and is empty at every point.
`hunt_results_pairs.npz`: per point index `i`, array `p{i}` with rows
`(U_n, U_{n+1}, y_n, source)`; source 0 or 1 = section seeds launched from two different `y`,
2 = orbit on the unstable manifold of `E`.

## `stests_summary.csv` (62 rows) and `stests_results.json`

Structural tests on the attractor disc at `mu` in {0.002, 0.01, 0.05, 0.2, 1, 5, 50} times
`rho - 1` in {0.003, 0.03, 0.3, 1.5, 2.5, 4, 8, 20}, and at (0.0027, 5.02), (0.002, 13.92),
(500, 34), (0.5, 5), (2, 5), (0.5, 12).
`S1_*`: the one-turn mean of `eta - mu` along the spiral of the unstable manifold (turn by turn,
and as a function of crest amplitude); `S1_a_at_E = eta_E - mu`, `S1_a_at_cycle = lambda_c`.
`S2_*[c]`: concavity of the crest return map in the coordinate `c` in {`U`, `y-k`, `e^U-1`,
`sqrtK`, `sqrtW`, `Lambda_E`}: `S2_concave_inside`, the largest finite-difference slope of the
derivative (`> 0` means convex somewhere), the extent of the convex window in `U/U_star`, the
largest elasticity (`< 1` everywhere means `phi(c)/c` decreasing), and the same outside the cycle.
`S3a_*`, `S3L_*`: sign of `d/dU[(mu - zeta)/(k(e^U - 1))]`, and of the analogue with
`(mu - zeta)/(1 + eta)`, along the section `p = 0` of the attractor (crest side `U > 0`, trough
side `U < 0`); `*_worst_normalised` lies in [-1, 1], negative means failure.
`S3b_disc_*`: the same derivative at fixed Lienard variable over the disc.
`S4_*`: smallest and median ratio of the eigenvalues of the Hessian of the graph function `h`
over reliably computed points; `hessian_reliable_fraction`.
`S5_*`: sign statistics of `K - (zeta - mu) chi(U)`, `chi(U) = U - 2(e^U - U - 1)/(e^U - 1)`.
`S6_*`: `zeta/eta` and the four ratios `x/zeta`, `zeta/eta + mu`, `k(Ue^U - e^U + 1)/K`,
`x p^2/W`, each minus `A`.
`smoothness_r = log m_s/log m_t` at the cycle.

## `cycle_tests.csv` (300 rows, the grid cycles)

`convex_projection`, `inflections`, `time_fraction_nonconvex`, `nonconvex_*`: sign changes of
`p p'' - p'^2` along the cycle. `min_zeta_over_eta`, `max_zeta_over_eta`,
`time_fraction_zeta_over_eta_le_1`. `R1..R4_minus_A_min/max`: pointwise ranges of
`R1 = x/zeta`, `R2 = zeta/eta + mu`, `R3 = k(Ue^U - e^U + 1)/K`, `R4 = x p^2/W`, minus `A`.
`Ehat_*`: `K - (zeta - mu) chi(U)`. `min_K_over_k`, `min_W`, `LambdaE_min_over_k`,
`max_p_over_lambda`: positivity gauges. `mean_x_over_y_minus_1`,
`mean_y_over_lam_minus_p_minus_mu`: the two Lotka-Volterra means by the trapezoid rule on the
recorded mesh; meaningful only for `rho - 1 < 0.5`, because at larger `rho` the integrands are
spikes that this mesh does not resolve.

## `s3a_cycle_scan.json`, `dipscan.json`

`s3a_cycle_scan.json`: the S3a quantity at the cycle crest for `mu` in {0.3, ..., 3} times `rho`
in {2.5, ..., 10}: `Qn_at_cycle`, `Qn_min_section`, `r` (the smoothness exponent).
`dipscan.json`: scan of the interior dip of `r(U) = log(phi(U)/U)` at small `mu`.

## `limit_map.csv` (85 rows), `limit_map_stdout.txt`

The formal `mu -> 0` limit at fixed `rho` (matched asymptotics, not a theorem): `p_in` (value of
`p` before the pulse), `G` (minus the value of `p` after the pulse of `x' = px`, `p' = -y`,
`y' = x - y`), `pulse_log_x_max`, `b_star`, `Lambda=rho_limit` (the `rho` whose limit cycle has
this `p_in`), `m0_limit_multiplier`, `log_m0`, `dG_dp`. `limit_map_stdout.txt` compares the
limit map with cycles of the full system at `mu` = 0.01, 0.002, 0.0005.

## Logs

`grid_log.txt`, `hunt_stdout.txt`, `stests_stdout.txt`, `tol_check_stdout.txt` (the weakest
exponents, the non-monotone return-map ratio at `mu = 0.0027`, `rho = 5.02` and the slope of the
section ratio at `mu = 1`, `rho = 5`, each at integrator tolerances 1e-10, 1e-12, 1e-13).
