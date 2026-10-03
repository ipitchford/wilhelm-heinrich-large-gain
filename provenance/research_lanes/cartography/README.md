# Cartography lane — files and columns

Everything here is ordinary floating point (numba DOP853, coefficients imported from scipy; default
relative tolerance 1e-12). Nothing is an interval certificate. Narrative, verdicts and commands: `NOTES.md`.

Python: `/Users/admin/.venvs/analysis-py313/bin/python`. All scripts are single-process.

## Code

| file | purpose |
|---|---|
| `cart.py` | compiled integrator (modes: state; state+eta,zeta,K,W+quadratures; state+variational+d/dk; state+pull-back+one transversal variation), crest return map, Newton at fixed k, amplitude-fixed Newton (unknowns y0,k), pulled-back periodic profile |
| `test_core.py` | self-test against the four v2.0 bundle receipts |
| `grid.py` | Part 1: amplitude continuation from Hopf, secant onto target rho, tables |
| `hunt.py` | Part 2: random-parameter hunt for a second cycle |
| `attractor.py` | sampling of the attractor disc with pulled-back eta, zeta, K, W and the Hessian of h |
| `stests.py` | Part 3: S1–S6 on the attractor disc |
| `cycle_tests.py` | S4b/S5/S6 on the 300 grid cycles |
| `s3scan.py` | S3a at the cycle crest over a (mu, rho) window |
| `portraits.py` | S6 phase portraits of the four source ratios |
| `limit_map.py` | formal mu→0 limit: pulse scattering map G, limit return map, comparison with the full system |
| `tol_check.py` | tolerance-tightening of the decisive numbers |
| `dipscan.py` | near-tangency scan of r(U) = log(phi(U)/U) at small mu (`dipscan.json`, `dipscan_stdout.txt`) |
| `melnikov_check.py` | numerical check of the Hamiltonian-limit identity ∮U_ss² = ∮U_s² |
| `scratch/`, `superseded/` | debugging scripts; the first hunt run (before the crest-event sign filter) |

## Tables

### `grid_table.csv` (300 rows = 12 mu × 25 rho) and `branch_table.csv` (827 rows, amplitude branches)
mu ∈ {0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.85, 1}; rho − 1 = logspace(1e-3 … 199, 25).
`branch_table.csv` has the same columns (no `secant_iters`) at fixed crest amplitudes h = U_max from 0.01 up to rho ≈ 260.

| column | meaning |
|---|---|
| `mu`, `rho`, `k` | parameters, rho = k/(mu(1+mu)) |
| `U_max`, `U_min` | crest and trough of U = log(x/k) on the cycle |
| `y0_over_k` | y/k at the crest (crest is (U_max, p=0, y0)) |
| `T`, `T_scaled` | least period in fast time t; T·sqrt(mu)/(2π) |
| `log_mt` | log m_t = T(⟨eta⟩ − mu) from the pulled-back Riccati quadrature (reliable at any size) |
| `log_mt_M` | log of the larger eigenvalue of the 2×2 crest-return derivative (blank when m_t < 1e-7) |
| `lambda_c`, `lambda_c_T` | ⟨eta⟩ − mu; and lambda_c·T (= `log_mt`) |
| `mean_eta`, `mean_U`, `meanU_plus_log_rho` | ⟨eta⟩, ⟨U⟩, ⟨U⟩ + log rho (brief item (g): must be < 0) |
| `q2`, `q3` | ⟨(eta−mu)²⟩, ⟨(eta−mu)³⟩ |
| `moment_identity_residual` | [(A²+mu)⟨q⟩ + 2A⟨q²⟩ + ⟨q³⟩ − (k−k_H)]/(k−k_H) (brief item (c)) |
| `mean_x_rel_err` | ⟨x⟩/k − 1 |
| `newton_residual` | crest-map residual of the located cycle (U in absolute, y relative to k) |
| `closure` | mismatch after integrating exactly one period with the profile system |
| `ratio_defect`, `settle_periods` | periodicity defect of (eta, zeta, K, W) and pull-back periods used |
| `eta_crest`, `zeta_crest` | pulled-back adjoint ratios at the crest |
| `secant_iters` | secant steps on the amplitude to hit the target rho |

### `cycle_tests.csv` (300 rows, the grid cycles)
`convex_projection`, `inflections`, `time_fraction_nonconvex`, `nonconvex_*`: sign changes of p·p'' − p'² along the cycle (S4b).
`min_zeta_over_eta`, `max_zeta_over_eta`, `time_fraction_zeta_over_eta_le_1` (S6).
`R1..R4_minus_A_min/max`: pointwise ranges of R1 = x/zeta, R2 = zeta/eta + mu, R3 = k(Ue^U − e^U + 1)/K, R4 = x p²/W, minus A = 1 + mu.
`Ehat_*`: (1+eta)E = K − (zeta − mu)·chi(U), chi(U) = U − 2(e^U − U − 1)/(e^U − 1) (S5).
`min_K_over_k`, `min_W`, `LambdaE_min_over_k`, `max_p_over_lambda`: positivity gauges.
`mean_x_over_y_minus_1`, `mean_y_over_lam_minus_p_minus_mu`: brief identities ⟨x/y⟩ = 1 and ⟨y/(lambda−p)⟩ = mu by trapezoid on the recorded mesh; **only meaningful for rho − 1 < 0.5** (errors ≤ 1.6e-8 and 8e-5); at larger rho the integrands are spikes and the trapezoid rule on this mesh is not accurate.

### `hunt_results.jsonl` (240 lines), `hunt_results_pairs.npz`, `hunt_summary.txt`, `hunt_results_supplement.jsonl`
`hunt_summary.txt`: aggregate counts. `hunt_results_supplement.jsonl`: the one point whose direct-coordinate run failed (index 110), redone at the reciprocal rate.
One JSON object per random point (seed 20261003; mu log-uniform in [1e-3, 1e3]; rho − 1 log-uniform in [1e-3, 299]).
Keys: `mu`, `rho`, `k`, `U_star`, `y_star_over_k`, `T`, `log_mt`, `lambda_c`, `m_t`, `branch_monotone` (rho(h) increasing along the amplitude branch up to the target), `primary_residual`;
scan: `scan_seeds`, `scan_pairs`, `scan_wu_pairs`, `scan_fail`, `scan_wrong_sign` (pairs with sign(U_{n+1} − U_n) ≠ sign(U* − U_n)), `scan_U_over_Ustar_range`, `scan_r_monotone_violations`, `scan_r_max_increase`, `scan_r_slope_max`, `scan_r_at_E` (r(U) = log(phi(U)/U));
long runs: `long` (list: initial condition tag, number of returns, final gap to the primary crest, status), `long_same`, `long_different`, `long_other`, `long_heading_bad`;
Newton: `newton` = {`same`, `equilibrium`, `failed`, `new` (list of distinct fixed points; empty everywhere)}.
The `.npz` holds, per point index i, array `p{i}` with rows (U_n, U_{n+1}, y_n, source) — source 0/1 = section seeds launched from two different y0, 2 = W^u(E) orbit.

### `stests_results.json` (full detail) and `stests_summary.csv` (62 rows, one per parameter point)
Parameter points: mu ∈ {0.002, 0.01, 0.05, 0.2, 1, 5, 50} × rho − 1 ∈ {0.003, 0.03, 0.3, 1.5, 2.5, 4, 8, 20}, plus (0.0027, 5.02), (0.002, 13.92), (500, 34), (0.5, 5), (2, 5), (0.5, 12).
Summary columns: `S1_*` (turn-by-turn and amplitude-monotone tests of the one-turn mean of eta − mu; `S1_a_at_E` = eta_E − mu, `S1_a_at_cycle` = lambda_c);
`S2_concave_inside[c]`, `S2_max_dphi1_dU_inside[c]` (largest finite-difference slope of phi_c' in U; > 0 means convex somewhere), `S2_convex_from/to_U/U*[c]` (extent of the convex window), `S2_max_elasticity_inside_U>0.01U*[c]` (elasticity = d log phi_c / d log c; < 1 everywhere ⇔ phi_c(c)/c decreasing), `S2_concave_outside[c]`, for c ∈ {U, y−k, e^U−1, sqrtK, sqrtW, Lambda_E};
`S3a_*`, `S3L_*`: sign of d/dU[(mu − zeta)/(k(e^U − 1))] and of d/dU[f/g], f = (mu − zeta)/(1+eta), along the attractor section p = 0 (crest side U > 0, trough side U < 0); `*_worst_normalised` ∈ [−1, 1] (negative = fails);
`S3b_disc_*`: the same derivative at fixed Liénard variable s over the whole disc; `min_(zeta-mu)/mu_on_U=0`;
`S4_*`: smallest/median ratio of the eigenvalues of D²h (in (U, p/sqrt(mu)) scaling) over reliably computed points; `hessian_reliable_fraction`;
`S5_*`: sign statistics of K − (zeta − mu)·chi(U); `S6_*`: zeta/eta and the four ratios minus A; `b_min_*`: positivity of K and W on the disc;
`smoothness_r` = log m_s / log m_t = A·T/|log m_t| − 1 at the cycle.

### `s3a_cycle_scan.json`
S3a quantity at the cycle crest for mu ∈ {0.3 … 3} × rho ∈ {2.5 … 10}: `Qn_at_cycle`, `Qn_min_section`, `r`.

### `limit_map.csv` (85 rows)
`p_in` (pre-pulse p), `G` (= −p after the pulse of x' = px, p' = −y, y' = x − y), `pulse_log_x_max`, `b_star`, `Lambda=rho_limit` (the rho whose mu→0 limit cycle has this p_in), `m0_limit_multiplier`, `log_m0`, `dG_dp`.

### `portrait_mu<mu>_rho<rho>.npz/.png` (three sets: (1, 2), (0.1, 1.5), (0.01, 3))
Arrays `t, U, p, y, x, eta, zeta, K, W, R1, R2, R3, R4` on one period from the crest, scalars `mu, rho, k, T, A, lambda_c`, `means` (the four cycle means of R_i, each equal to 1 + ⟨eta⟩).

### Logs
`grid_log.txt`, `grid_stdout.txt`, `hunt_stdout.txt`, `stests_stdout.txt`, `limit_map_stdout.txt`, `tol_check_stdout.txt`, `dipscan_stdout.txt`, `melnikov_check_stdout.txt`.
