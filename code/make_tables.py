#!/usr/bin/env python3
"""Generate every number that the paper takes from the evidence, directly from the strict JSON receipts.

Writes, under paper/generated/:
  evidence_macros.tex, evidence_counts.tex, cycle_examples.tex, phase_example.tex   (framework suite, as in v2.0)
  v3_counts.tex               table: check, kind, passed, attempted
  v3_macros.tex               headline numbers and worst residuals as macros \\VthreeXxx (letters only)
  v3_cartography.tex          log m_t on a subgrid of mu x rho of the shipped grid table
  v3_largegain_constants.tex  the constants of the large-gain proof as recomputed at 40 digits
and evidence/manuscript_summary.json, evidence/v3/v3_summary.json (macro values with the SHA-256 of every source receipt).

Macros whose value contains \\times or a minus sign are meant for math mode.
`build_v3()` returns the generated texts without writing; check_consistency.py uses it to verify that the files on disk
are exactly what the receipts give.
"""
from __future__ import annotations
import csv, hashlib, json, sys
from pathlib import Path
import numpy as np
from cycles import ROOT, strict_write

GEN = ROOT / 'paper' / 'generated'
EV3 = ROOT / 'evidence' / 'v3'
V3_RECEIPTS = ['v3_exact_identities', 'v3_gauges', 'v3_largegain_constants', 'v3_largegain_lemmas', 'v3_attractor',
               'v3_weak_adjoint', 'v3_wedge', 'v3_cartography', 'v3_negative', 'v3_controls']
LABELS = {
    'v3_exact_identities': ('New algebraic identities', 'exact', 'residual tested equal to zero'),
    'v3_gauges': ('Gauges and exponent formulas on cycles', 'floating point', 'cycle'),
    'v3_largegain_constants': ('Constants of the large-gain proof', 'high precision', 'scalar comparison'),
    'v3_largegain_lemmas': ('Large-gain lemmas on cycles and on non-periodic segments', 'floating point', 'inequality evaluated'),
    'v3_attractor': ('Attractor: supporting planes, gradient, crest-axis formula', 'floating point', 'parameter pair and predicate'),
    'v3_weak_adjoint': ('Weak adjoint: event order, dual planes, rank-one identity', 'floating point', 'cycle and predicate'),
    'v3_wedge': ('Kurtosis bound, period identities, third-order formula', 'high precision and floating point', 'level, identity, cycle or table row'),
    'v3_cartography': ('Cartography: shipped tables and recomputed subsample', 'floating point', 'table statement, cycle or search'),
    'v3_negative': ('Refutations reproduced', 'floating point', 'refutation, identity or candidate'),
    'v3_controls': ('Deliberate-error controls', 'control', 'wrong assertion rejected'),
}


def read(n):
    return json.loads((ROOT / 'evidence' / n).read_text())


def read3(n):
    def bad(x):
        raise ValueError('non-standard JSON constant ' + x)
    return json.loads((EV3 / f'{n}.json').read_text(), parse_constant=bad)


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def sci(v, digits=3):
    if v == 0:
        return '0'
    s = f'{v:.{digits}e}'
    a, b = s.split('e')
    return a + r'\times10^{' + str(int(b)) + '}'


def grp(n):
    """integer with thin-space grouping (valid in text and in math mode)"""
    s = f'{abs(int(n)):,}'.replace(',', r'\,')
    return ('-' if n < 0 else '') + s


def fix(v, d):
    return f'{v:.{d}f}'


def sig3(v):
    """three significant digits, trailing zeros kept, no exponent (for |v| >= 0.1)"""
    s = np.format_float_positional(float(v), precision=3, unique=False, fractional=False, trim='k')
    return s[:-1] if s.endswith('.') else s


# ----------------------------------------------------------------------------------------------- framework suite (v2.0)
def run_v2():
    exact = read('exact_checks.json'); trial = read('cycle_trials.json'); num = read('numerical_audit.json'); cross = read('cross_method.json')
    rec = read('reciprocal_rates.json'); env = read('envelope_checks.json'); neg = read('negative_controls.json')
    rows = num['records']; by = {r['id']: r for r in rows}
    summary = {'schema_version': '1.0', 'claim_boundary': 'Finite numerical diagnostics, not interval certificates and not exhaustive cycle exclusion.',
               'counts': {'exact_equalities': [exact['number_passed'], exact['number_checked']],
                          'cycle_approximations': [trial['number_accepted'], trial['number_attempted']],
                          'cycle_geometry_and_adjoint_audits': [num['number_passed'], num['number_attempted']],
                          'cross_method': [cross['number_passed'], cross['number_checked']],
                          'reciprocal_rates': [rec['number_passed'], rec['number_checked']],
                          'envelope': [env['number_passed'], env['number_checked']],
                          'negative_controls': [neg['number_passed'], neg['number_checked']]},
               'maxima': {'exponent_cross_formulation': max(r['exponent_difference'] for r in rows),
                          'orthogonality_scaled': max(r['orthogonality_max_scaled'] for r in rows),
                          'state_return_closure': max(r['return_closure_error'] for r in rows),
                          'moment_relative_residual': max(v for r in rows for v in r['moment_relative_residuals'].values()),
                          'cross_method_matrix_difference': max(r['return_matrix_max_difference'] for r in cross['records']),
                          'cross_method_exponent_difference': max(r['exponent_difference'] for r in cross['records']),
                          'reciprocal_all_comparisons': max(v for r in rec['records'] for key, v in r.items() if key.startswith('max_difference') or key.startswith('difference_')),
                          'envelope_hessian_scaled_error': max(r['relative_scaled_error'] for r in env['records']),
                          'envelope_defect_residual': max(r['defect_identity_max_residual'] for r in env['records'])},
               'square_root_certificate_failures': sum(not all(v < 1 for v in r['square_root_test_ratios']) for r in rows),
               'number_negative_computed_exponents': sum(r['lambda_c_projective'] < 0 for r in rows),
               'phase_example': by['mu_1_h_1']['phase_fractions'],
               'false_radial_example': {k: by['mu_0.3_h_0.01'][k] for k in ['mu', 'amplitude', 'square_root_test_ratios', 'weak_multiplier']}}
    src = ['exact_checks.json', 'cycle_trials.json', 'numerical_audit.json', 'cross_method.json', 'reciprocal_rates.json', 'envelope_checks.json', 'negative_controls.json']
    summary['source_sha256'] = {n: sha(ROOT / 'evidence' / n) for n in src}
    summary['passed'] = all(exact['passed'] and v[0] == v[1] for v in summary['counts'].values())
    out = {}
    names = {'exponent_cross_formulation': 'ExponentAgreement', 'orthogonality_scaled': 'OrthogonalityError', 'state_return_closure': 'ReturnClosure',
             'moment_relative_residual': 'MomentError', 'cross_method_matrix_difference': 'CrossMatrixError',
             'cross_method_exponent_difference': 'CrossExponentError', 'reciprocal_all_comparisons': 'ReciprocalError',
             'envelope_hessian_scaled_error': 'EnvelopeHessianError', 'envelope_defect_residual': 'EnvelopeDefectError'}
    macros = ['\\newcommand{\\' + n + '}{' + sci(summary['maxima'][k]) + '}' for k, n in names.items()]
    out['evidence_macros.tex'] = '% Automatically generated; do not edit.\n' + '\n'.join(macros) + '\n'
    labels = {'exact_equalities': 'Exact symbolic equalities', 'cycle_approximations': 'Fresh periodic boundary-value solves',
              'cycle_geometry_and_adjoint_audits': 'Cycle geometry and adjoint audits', 'cross_method': 'DOP853--Radau repetitions',
              'reciprocal_rates': 'Reciprocal-rate coordinate comparisons', 'envelope': 'Local envelope Hessian checks',
              'negative_controls': 'Deliberate-error controls'}
    lines = [r'\begin{tabular}{lrr}', r'\toprule', r'Check & Accepted/passed & Attempted \\', r'\midrule']
    for k, v in summary['counts'].items():
        lines.append(f'{labels[k]} & {v[0]} & {v[1]} ' + r'\\')
    lines.extend([r'\bottomrule', r'\end{tabular}'])
    out['evidence_counts.tex'] = '\n'.join(lines) + '\n'
    lines = [r'\begin{tabular}{rrrrr}', r'\toprule', r'$\mu$ & $h$ & $T$ (fast time) & $\lambda_c$ & $m_t$ \\', r'\midrule']
    for ident in ['mu_0.1_h_3', 'mu_0.3_h_0.01', 'mu_1_h_1', 'mu_10_h_3']:
        r = by[ident]
        lines.append(f"{r['mu']:.3g} & {r['amplitude']:.3g} & {r['period_fast']:.10f} & ${sci(r['lambda_c_projective'])}$ & {r['weak_multiplier']:.10f} " + r'\\')
    lines.extend([r'\bottomrule', r'\end{tabular}'])
    out['cycle_examples.tex'] = '\n'.join(lines) + '\n'
    lines = [r'\begin{tabular}{lrrrrrr}', r'\toprule', r'Event & $X_{\max}$ & $Z_{\max}$ & $Y_{\max}$ & $X_{\min}$ & $Z_{\min}$ & $Y_{\min}$ \\', r'\midrule']
    lines.append('Fraction of period & ' + ' & '.join(f'{x:.6f}' for x in summary['phase_example'].values()) + r'\\')
    lines.extend([r'\bottomrule', r'\end{tabular}'])
    out['phase_example.tex'] = '\n'.join(lines) + '\n'
    return summary, out


# ----------------------------------------------------------------------------------------------- v3
def build_v3():
    R = {n: read3(n) for n in V3_RECEIPTS}
    v2 = {n: read(n + '.json') for n in ['exact_checks', 'cycle_trials', 'numerical_audit', 'cross_method', 'reciprocal_rates', 'envelope_checks', 'negative_controls']}
    out = {}
    M = {}                                  # macro name -> (tex value, source description)

    def put(name, tex, source):
        if not name.isalpha():
            raise ValueError('macro names use letters only: ' + name)
        if name in M:
            raise ValueError('duplicate macro ' + name)
        M[name] = (tex, source)

    # counts table
    P = r'>{\raggedright\arraybackslash}p'
    lines = [r'\begin{tabular}{@{}' + P + '{5.3cm}' + P + '{2.7cm}rr' + P + r'{3.5cm}@{}}', r'\toprule',
             r'Check & Kind & Passed & Attempted & Unit counted \\', r'\midrule']
    v2rows = [('Framework identities', 'exact', v2['exact_checks']['number_passed'], v2['exact_checks']['number_checked'], 'residual tested equal to zero'),
              ('Periodic boundary-value solves', 'floating point', v2['cycle_trials']['number_accepted'], v2['cycle_trials']['number_attempted'], 'cycle'),
              ('Cycle geometry and adjoint audits', 'floating point', v2['numerical_audit']['number_passed'], v2['numerical_audit']['number_attempted'], 'cycle'),
              ('DOP853--Radau repetitions', 'floating point', v2['cross_method']['number_passed'], v2['cross_method']['number_checked'], 'cycle'),
              ('Reciprocal-rate comparisons', 'floating point', v2['reciprocal_rates']['number_passed'], v2['reciprocal_rates']['number_checked'], 'pair of cycles'),
              ('Envelope Hessian checks', 'floating point', v2['envelope_checks']['number_passed'], v2['envelope_checks']['number_checked'], 'phase'),
              ('Deliberate-error controls (framework)', 'control', v2['negative_controls']['number_passed'], v2['negative_controls']['number_checked'], 'wrong assertion rejected')]
    for lab, kind, p_, a_, unit in v2rows:
        lines.append(f'{lab} & {kind} & {grp(p_)} & {grp(a_)} & {unit} ' + r'\\')
    lines.append(r'\midrule')
    tot_p = tot_a = 0
    for n in V3_RECEIPTS:
        lab, kind, unit = LABELS[n]
        c = R[n]['counts']
        tot_p += c['passed']; tot_a += c['attempted']
        lines.append(f'{lab} & {kind} & {grp(c["passed"])} & {grp(c["attempted"])} & {unit} ' + r'\\')
    lines.extend([r'\bottomrule', r'\end{tabular}'])
    out['v3_counts.tex'] = '% Automatically generated by code/make_tables.py; do not edit.\n' + '\n'.join(lines) + '\n'

    # macros
    ex, ga, lc, ll, at, wk, wd, ca, ng, ct = [R[n] for n in V3_RECEIPTS]
    w = ga['worst_residuals']
    put('VthreeExactPassed', grp(ex['counts']['passed']), 'v3_exact_identities counts.passed')
    put('VthreeExactAttempted', grp(ex['counts']['attempted']), 'v3_exact_identities counts.attempted')
    put('VthreeFrameworkExact', grp(v2['exact_checks']['number_checked']), 'exact_checks number_checked')
    put('VthreeGaugeCycles', grp(ga['counts']['attempted']), 'v3_gauges counts.attempted')
    put('VthreeGaugeCyclesPassed', grp(ga['counts']['passed']), 'v3_gauges counts.passed')
    put('VthreeGaugeSaved', grp(ga['number_saved_cycles']), 'v3_gauges number_saved_cycles')
    put('VthreeGaugeExtended', grp(ga['number_extended_cycles']), 'v3_gauges number_extended_cycles')
    put('VthreeGaugeSamples', grp(w['total_samples']), 'v3_gauges worst.total_samples')
    put('VthreeGaugeUnresolved', grp(w['unresolved_sign_samples']), 'v3_gauges worst.unresolved_sign_samples')
    put('VthreeGaugeExponentGap', sci(w['max_exponent_gap_relative']), 'v3_gauges worst.max_exponent_gap_relative')
    put('VthreeGaugeEtaZeroGap', sci(w['max_eta0_mean_relative']), 'v3_gauges worst.max_eta0_mean_relative')
    put('VthreeGaugePolarError', sci(w['max_polar_angle_error']), 'v3_gauges worst.max_polar_angle_error')
    put('VthreeGaugeMinSigma', sci(w['min_star_gauge_ratio']), 'v3_gauges worst.min_star_gauge_ratio')
    put('VthreeGaugeMinXi', sci(w['min_convexity_gauge_ratio']), 'v3_gauges worst.min_convexity_gauge_ratio')
    put('VthreeGaugeMinK', sci(w['min_K_ratio']), 'v3_gauges worst.min_K_ratio')
    put('VthreeGaugeMinW', sci(w['min_W_ratio']), 'v3_gauges worst.min_W_ratio')
    put('VthreeGaugeMaxLambda', sci(w['max_lambda_c']), 'v3_gauges worst.max_lambda_c')
    put('VthreeGaugeMaxFormEigenvalue', sci(w['max_monotone_form_eigenvalue']), 'v3_gauges worst.max_monotone_form_eigenvalue')
    put('VthreeGaugeUmin', grp(round(ga['ranges']['U_min'])), 'v3_gauges ranges.U_min (rounded)')
    put('VthreeGaugeHolderLower', fix(w['max_holder_lower_over_middle'], 6), 'v3_gauges worst.max_holder_lower_over_middle')
    put('VthreeGaugeHolderUpper', fix(w['max_holder_middle_over_upper'], 6), 'v3_gauges worst.max_holder_middle_over_upper')
    put('VthreeGaugeEtaZeroNegative', grp(ga['cycles_with_eta0_negative_somewhere']), 'v3_gauges cycles_with_eta0_negative_somewhere')
    put('VthreeGaugeEtaOneNegative', grp(ga['cycles_with_eta1_negative_somewhere']), 'v3_gauges cycles_with_eta1_negative_somewhere')
    put('VthreeGaugeAttracting', grp(ga['number_attracting']), 'v3_gauges number_attracting')

    put('VthreeLgConstChecks', grp(lc['counts']['attempted']), 'v3_largegain_constants counts.attempted')
    put('VthreeLgConstAuthor', grp(lc['number_author_comparisons']), 'v3_largegain_constants number_author_comparisons')
    put('VthreeLgConstReferee', grp(lc['number_referee_checks']), 'v3_largegain_constants number_referee_checks')
    put('VthreeLgConstMinMargin', sci(lc['worst_residuals']['min_relative_margin_author_script'], 1), 'v3_largegain_constants worst.min_relative_margin_author_script')
    g = ll['groups']
    put('VthreeLgLemmaChecks', grp(ll['counts']['attempted']), 'v3_largegain_lemmas counts.attempted')
    put('VthreeLgCycleChecks', grp(g['A referee cycles']['attempted'] + g['B author cycles']['attempted']), 'v3_largegain_lemmas groups A + B attempted')
    put('VthreeLgBelowChecks', grp(g['B2 author cycles below the threshold']['attempted']), 'v3_largegain_lemmas group B2 attempted')
    put('VthreeLgSegmentChecks', grp(sum(g[k]['attempted'] for k in g if k[0] in 'CDEF')), 'v3_largegain_lemmas groups C-F attempted')
    lw = ll['worst_residuals']
    put('VthreeLgMaxEtaRatio', fix(lw['max_mean_eta_over_mu_on_cycles'], 3), 'v3_largegain_lemmas worst.max_mean_eta_over_mu_on_cycles')
    put('VthreeLgMinPeriodRatio', fix(lw['min_T_mu_over_lambda_on_cycles'], 3), 'v3_largegain_lemmas worst.min_T_mu_over_lambda_on_cycles')
    put('VthreeLgMinStageRatio', fix(lw['min_pk_over_k_stage1_nonperiodic'], 3), 'v3_largegain_lemmas worst.min_pk_over_k_stage1_nonperiodic')
    lb = lw['log_mt_below_threshold']
    put('VthreeLgLogMtSixteen', fix(lb['mu_1_k_16'], 2), 'v3_largegain_lemmas worst.log_mt_below_threshold.mu_1_k_16')
    put('VthreeLgLogMtSixty', fix(lb['mu_1_k_60'], 1), 'v3_largegain_lemmas worst.log_mt_below_threshold.mu_1_k_60')
    put('VthreeLgLogMtThreeHundred', fix(lb['mu_1_k_300'], 1), 'v3_largegain_lemmas worst.log_mt_below_threshold.mu_1_k_300')

    aw = at['worst_residuals']
    put('VthreeAttPairs', grp(at['parameter_pairs'] + at['cancellation_free_pairs']), 'v3_attractor parameter_pairs + cancellation_free_pairs')
    put('VthreeAttBatteryPairs', grp(at['parameter_pairs']), 'v3_attractor parameter_pairs')
    put('VthreeAttChecks', grp(at['counts']['attempted']), 'v3_attractor counts.attempted')
    put('VthreeAttOrderedPairs', sci(aw['ordered_pairs_total'], 2), 'v3_attractor worst.ordered_pairs_total')
    put('VthreeAttNonnegPairs', grp(aw['pairs_nonnegative_total']), 'v3_attractor worst.pairs_nonnegative_total')
    put('VthreeAttBeyondRounding', grp(aw['pairs_beyond_rounding_total']), 'v3_attractor worst.pairs_beyond_rounding_total')
    put('VthreeAttWorstSeparated', sci(aw['max_support_value_separated_pairs'], 1), 'v3_attractor worst.max_support_value_separated_pairs')
    put('VthreeAttGradPairs', grp(aw['gradient_fd_pairs_applicable']), 'v3_attractor worst.gradient_fd_pairs_applicable')
    put('VthreeAttGradSamples', grp(aw['gradient_fd_samples_total']), 'v3_attractor worst.gradient_fd_samples_total')
    put('VthreeAttGradError', sci(aw['max_gradient_fd_error'], 1), 'v3_attractor worst.max_gradient_fd_error')
    put('VthreeAttGradMedian', sci(aw['max_gradient_fd_median_error'], 1), 'v3_attractor worst.max_gradient_fd_median_error')
    put('VthreeAttReturnError', sci(aw['max_return_derivative_error'], 1), 'v3_attractor worst.max_return_derivative_error')
    put('VthreeAttReturnFdError', sci(aw['max_return_derivative_fd_error'], 1), 'v3_attractor worst.max_return_derivative_fd_error')
    put('VthreeAttExponentGap', sci(aw['max_exponent_difference_riccati_vs_qr'], 1), 'v3_attractor worst.max_exponent_difference_riccati_vs_qr')
    put('VthreeAttCrestValueError', sci(aw['max_crest_axis_value_error'], 1), 'v3_attractor worst.max_crest_axis_value_error')

    ww = wk['worst_residuals']
    put('VthreeWeakCycles', grp(wk['cycles']), 'v3_weak_adjoint cycles')
    put('VthreeWeakChecks', grp(wk['counts']['attempted']), 'v3_weak_adjoint counts.attempted')
    put('VthreeWeakChordPairs', grp(ww['chord_pairs_total']), 'v3_weak_adjoint worst.chord_pairs_total')
    put('VthreeWeakChordViolations', grp(ww['chord_violations_total']), 'v3_weak_adjoint worst.chord_violations_total')
    put('VthreeWeakWindowError', sci(ww['max_window_error_counted_test_functions'], 1), 'v3_weak_adjoint worst.max_window_error_counted_test_functions')
    put('VthreeWeakTies', grp(ww['cycles_with_tied_events']), 'v3_weak_adjoint worst.cycles_with_tied_events')
    put('VthreeWeakPhiYGood', grp(ww['cycles_with_phi_y_error_below_1e-5']), 'v3_weak_adjoint worst.cycles_with_phi_y_error_below_1e-5')
    put('VthreeWeakExponentSpread', sci(ww['max_exponent_spread'], 1), 'v3_weak_adjoint worst.max_exponent_spread')

    dw = wd['worst_residuals']
    put('VthreeWedgeChecks', grp(wd['counts']['attempted']), 'v3_wedge counts.attempted')
    put('VthreeWedgeKurtosisToda', fix(dw['max_kurtosis_toda'], 7), 'v3_wedge worst.max_kurtosis_toda')
    put('VthreeWedgeKurtosisTodaMargin', sci(1.8 - dw['max_kurtosis_toda'], 1), '9/5 - v3_wedge worst.max_kurtosis_toda')
    put('VthreeWedgeKurtosisOtherMargin', sci(1.8 - dw['max_kurtosis_other_convex'], 1), '9/5 - v3_wedge worst.max_kurtosis_other_convex')
    put('VthreeWedgeMinRselPrime', fix(dw['min_Rsel_prime'], 7), 'v3_wedge worst.min_Rsel_prime')
    put('VthreeWedgeIdentityDefect', sci(dw['max_identity_defect'], 1), 'v3_wedge worst.max_identity_defect')
    put('VthreeWedgeIdentityCycles', grp(dw['cycles']), 'v3_wedge worst.cycles')
    put('VthreeWedgeSelectionError', sci(dw['max_selection_error_smallest_mu_amplitude_le_1'], 1), 'v3_wedge worst.max_selection_error_smallest_mu_amplitude_le_1')
    put('VthreeWedgeExponentError', sci(dw['max_exponent_error_smallest_mu_amplitude_le_1'], 1), 'v3_wedge worst.max_exponent_error_smallest_mu_amplitude_le_1')
    ratios = [r['ratio'] for r in wd['records'] if 'ratio' in r]
    put('VthreeWedgeRatioMin', fix(min(ratios), 1), 'v3_wedge records: smallest error ratio between mu = 1e-3 and mu = 1e-4')
    put('VthreeWedgeRatioMax', fix(max(ratios), 1), 'v3_wedge records: largest error ratio between mu = 1e-3 and mu = 1e-4')
    put('VthreeWedgeRatioPredicted', fix(wd['records'][[i for i, r in enumerate(wd['records']) if 'ratio' in r][0]]['predicted_for_second_order'], 2),
        'v3_wedge records: predicted_for_second_order')
    put('VthreeWedgeSelectionErrorAll', sci(dw['max_selection_error_smallest_mu'], 1), 'v3_wedge worst.max_selection_error_smallest_mu')
    put('VthreeWedgeExponentErrorAll', sci(dw['max_exponent_error_smallest_mu'], 1), 'v3_wedge worst.max_exponent_error_smallest_mu')

    cw = ca['worst_residuals']
    put('VthreeCartChecks', grp(ca['counts']['attempted']), 'v3_cartography counts.attempted')
    put('VthreeCartGridCycles', grp(cw['grid_cycles']), 'v3_cartography worst.grid_cycles')
    put('VthreeCartBranchCycles', grp(cw['branch_cycles']), 'v3_cartography worst.branch_cycles')
    put('VthreeCartHuntPoints', grp(cw['hunt_points']), 'v3_cartography worst.hunt_points')
    put('VthreeCartHuntDirect', grp(cw['hunt_points_computed_directly']), 'v3_cartography worst.hunt_points_computed_directly')
    put('VthreeCartMaxLambda', sci(cw['grid_max_lambda_c'], 2), 'v3_cartography worst.grid_max_lambda_c')
    put('VthreeCartMaxLogMt', sci(cw['grid_max_log_mt'], 2), 'v3_cartography worst.grid_max_log_mt')
    put('VthreeCartMinLogMt', fix(cw['grid_min_log_mt'], 1), 'v3_cartography worst.grid_min_log_mt')
    put('VthreeCartBranchMaxLambda', sci(cw['branch_max_lambda_c'], 1), 'v3_cartography worst.branch_max_lambda_c')
    put('VthreeCartHuntMaxLogMt', sci(cw['hunt_max_log_mt'], 2), 'v3_cartography worst.hunt_max_log_mt')
    put('VthreeCartHuntPairs', grp(cw['hunt_pairs']), 'v3_cartography worst.hunt_pairs')
    put('VthreeCartHuntLongRuns', grp(cw['hunt_long_runs']), 'v3_cartography worst.hunt_long_runs')
    put('VthreeCartHuntNotFinished', grp(cw['hunt_long_not_finished']), 'v3_cartography worst.hunt_long_not_finished')
    put('VthreeCartHuntNewton', grp(cw['hunt_newton_starts']), 'v3_cartography worst.hunt_newton_starts')
    put('VthreeCartRecomputedCycles', grp(cw['recomputed_grid_cycles']), 'v3_cartography worst.recomputed_grid_cycles')
    put('VthreeCartRecomputedHunts', grp(cw['recomputed_hunts']), 'v3_cartography worst.recomputed_hunts')
    put('VthreeCartRecomputedDiffT', sci(cw['recomputed_grid_max_rel_diff_T'], 1), 'v3_cartography worst.recomputed_grid_max_rel_diff_T')
    put('VthreeCartRecomputedDiffLogMt', sci(cw['recomputed_grid_max_diff_log_mt'], 1), 'v3_cartography worst.recomputed_grid_max_diff_log_mt')
    put('VthreeCartRecomputedHuntPairs', grp(cw['recomputed_hunt_pairs']), 'v3_cartography worst.recomputed_hunt_pairs')
    put('VthreeCartRecomputedLongRuns', grp(cw['recomputed_hunt_long_runs']), 'v3_cartography worst.recomputed_hunt_long_runs')
    put('VthreeCartRecomputedNotFinished', grp(cw['recomputed_hunt_long_runs_not_finished']), 'v3_cartography worst.recomputed_hunt_long_runs_not_finished')
    put('VthreeCartRecomputedNewton', grp(cw['recomputed_hunt_newton_starts']), 'v3_cartography worst.recomputed_hunt_newton_starts')
    put('VthreeCartRecomputedNewPoints', grp(cw['recomputed_hunt_new_fixed_points']), 'v3_cartography worst.recomputed_hunt_new_fixed_points')
    put('VthreeCartIndependentChecks', grp(cw['independent_primary_checks']), 'v3_cartography worst.independent_primary_checks')
    put('VthreeCartIndependentDiffLogMt', sci(cw['independent_max_diff_log_mt'], 1), 'v3_cartography worst.independent_max_diff_log_mt')
    put('VthreeCartNegativeEhatCycles', grp(cw['cycle_tests_cycles_with_negative_Ehat']), 'v3_cartography worst.cycle_tests_cycles_with_negative_Ehat')

    nw = ng['worst_residuals']
    put('VthreeNegChecks', grp(ng['counts']['attempted']), 'v3_negative counts.attempted')
    put('VthreeNegRise', fix(nw['return_ratio_rise'], 4), 'v3_negative worst.return_ratio_rise')
    put('VthreeNegRatioAtE', fix(nw['return_ratio_at_E'], 4), 'v3_negative worst.return_ratio_at_E')
    put('VthreeNegSlope', fix(nw['section_ratio_slope_near_cycle'], 3), 'v3_negative worst.section_ratio_slope_near_cycle')
    put('VthreeNegSlopeMid', fix(nw['section_ratio_slope_mid'], 3), 'v3_negative worst.section_ratio_slope_mid')
    put('VthreeNegSignChanges', grp(nw['candidate_sign_changes']), 'v3_negative worst.candidate_sign_changes')
    put('VthreeControls', grp(ct['counts']['attempted']), 'v3_controls counts.attempted')
    put('VthreeControlsPassed', grp(ct['counts']['passed']), 'v3_controls counts.passed')
    put('VthreeTotalAttempted', grp(tot_a), 'sum of counts.attempted over the ten v3 receipts')
    put('VthreeTotalPassed', grp(tot_p), 'sum of counts.passed over the ten v3 receipts')

    rr_path = ROOT / 'evidence' / 'run_receipt.json'
    if rr_path.exists():
        rr = json.loads(rr_path.read_text())
        put('VthreeRuntimeMinutes', fix(rr['computation_seconds'] / 60.0, 1), 'run_receipt computation_seconds / 60')
        put('VthreePeakMemoryMb', grp(round(rr['peak_rss_mb'])), 'run_receipt peak_rss_mb (rounded)')
    # (the macro file is assembled below, after the macros taken from the grid table)

    # cartography table from the shipped grid
    with (EV3 / 'cartography' / 'grid_table.csv').open() as f:
        grid = [{k: float(v) if v else None for k, v in r.items()} for r in csv.DictReader(f)]
    mus = [0.002, 0.01, 0.05, 0.2, 0.5, 1.0]
    idx = [0, 5, 10, 14, 18, 24]
    lines = [r'\begin{tabular}{@{}l' + 'r' * len(mus) + r'@{}}', r'\toprule',
             r'$\rho-1$ & ' + ' & '.join(f'$\\mu={m:g}$' for m in mus) + r' \\', r'\midrule']
    cart_values = {}
    for i in idx:
        row = []
        rho = None
        for m in mus:
            col = sorted((r for r in grid if r['mu'] == m), key=lambda r: r['rho'])
            rho = col[i]['rho']
            v = col[i]['log_mt']
            cart_values[f'mu={m:g},rho_index={i}'] = v
            row.append(f'${sci(v, 2)}$' if abs(v) < 0.1 else f'${sig3(v)}$')
        lines.append(f'${sig3(rho - 1)}$ & ' + ' & '.join(row) + r' \\')
    lines.extend([r'\bottomrule', r'\end{tabular}'])
    top = [sorted((r for r in grid if r['mu'] == m), key=lambda r: r['rho'])[-1]['log_mt'] for m in sorted({r['mu'] for r in grid})]
    put('VthreeCartTopRhoMaxLogMt', fix(max(top), 1), 'evidence/v3/cartography/grid_table.csv: largest log_mt at the largest rho')
    put('VthreeCartTopRhoMinLogMt', fix(min(top), 1), 'evidence/v3/cartography/grid_table.csv: smallest log_mt at the largest rho')
    out['v3_macros.tex'] = ('% Automatically generated by code/make_tables.py; do not edit.\n'
                            '% Values containing \\times or a minus sign are for math mode.\n'
                            + ''.join('\\newcommand{\\' + n + '}{' + t + '}\n' for n, (t, _) in M.items()))
    out['v3_cartography.tex'] = '% Automatically generated by code/make_tables.py from evidence/v3/cartography/grid_table.csv; do not edit.\n' + '\n'.join(lines) + '\n'

    # large-gain constants, in the notation of the large-gain section
    h = lc['headline_constants_at_K0']
    val = lambda k_, d=5: f"{float(h[k_]['value']):.{d}g}"
    rowsc = [
        (r'$\omega_0(k_0)$, $T_0(k_0)$', f"${val('om0')}$, ${val('T0', 6)}$", r'$T_0(k_0)\ge a_1+a_2+a_3+a_4=' + val('window_sum', 3) + '$'),
        (r'$k_0c_{\rm i}(a_1)$, $\chi_{\rm ii}(a_2)$, $k_0c_{\rm b}(a_4)$', f"${val('k0_ci_a1')}$, ${val('chi_ii_a2', 4)}$, ${val('k0_cb_a4')}$",
         r'$\ge2$, $\ge2$, $\ge L_\ast/\varphi_\ast=' + val('strip_ratio') + '$'),
        (r'bound for $c_1$', f"${val('c1_bound')}$", r'$\ge c_1=0.374$'),
        (r'$\tilde R(k_0)$, $\bar\varepsilon_1$, $\vartheta_1(k_0)$, $T_1(k_0)$', f"${val('Ht')}$, ${val('eps1', 4)}$, ${val('d1')}$, ${val('T1')}$", r'$\vartheta_1<1$, $T_1\ge22$'),
        (r'$z_{10}$, $\chi_0$', f"${val('z10')}$, ${val('kap0')}$", ''),
        (r'bound for $c_p$; $c_T$', f"${val('cp_bound')}$; ${val('cT')}$", r'$\ge c_p=0.819$'),
        (r'$\hat R$, $\hat L$, $\hat\varepsilon_1$ at $\lambda=k_0$', f"${val('Hhat', 6)}$, ${val('Lhat')}$, ${val('e1', 4)}$", ''),
        (r'$\hat\vartheta_N$, $a_c$, $A_1$ at $\lambda=k_0$', f"${val('dp', 4)}$, ${val('ac')}$, ${val('A1', 6)}$", r'$\hat\vartheta_N<1$'),
        (r'$\bar\tau_1$, $\hat w$, $v_f$ at $\lambda=k_0$', f"${val('tau1')}$, ${val('W')}$, ${val('vf')}$", ''),
        (r'$C_\tau(k_0)$, $c_{\log}$', f"${val('Ctau')}$, ${val('cb_log')}$", r'$\le71.79$, $\le12.28$'),
        (r'$q(k_0)$', f"${val('q')}$", r'$<0.0772$'),
        (r'final margin', f"${val('margin_term_1', 4)}-{val('margin_term_2', 3)}-{val('margin_term_3', 3)}={val('Phi_', 4)}$", r'$>0$'),
    ]
    lines = [r'\begin{tabular}{@{}lll@{}}', r'\toprule', r'Quantity & Recomputed value & Compared with \\', r'\midrule']
    lines += [f'{a} & {b} & {c} ' + r'\\' for a, b, c in rowsc]
    lines.extend([r'\bottomrule', r'\end{tabular}'])
    out['v3_largegain_constants.tex'] = '% Automatically generated by code/make_tables.py from evidence/v3/v3_largegain_constants.json; do not edit.\n' + '\n'.join(lines) + '\n'

    summary = {'schema_version': '3.0', 'macros': {n: {'tex': t, 'source': s_} for n, (t, s_) in M.items()},
               'cartography_table_values': cart_values,
               'totals': {'attempted': tot_a, 'passed': tot_p},
               'receipts': {n: {'kind': R[n]['kind'], 'passed': R[n]['passed'], 'attempted': R[n]['counts']['attempted'],
                                'n_passed': R[n]['counts']['passed']} for n in V3_RECEIPTS},
               'source_sha256': {f'evidence/v3/{n}.json': sha(EV3 / f'{n}.json') for n in V3_RECEIPTS}
                                | {'evidence/v3/cartography/grid_table.csv': sha(EV3 / 'cartography' / 'grid_table.csv')},
               'passed': all(R[n]['passed'] for n in V3_RECEIPTS),
               'claim_boundary': 'Finite computations. Exact symbolic equalities are exact; everything else is floating point or '
                                 'high-precision arithmetic without interval enclosure.'}
    return summary, out


def run():
    GEN.mkdir(parents=True, exist_ok=True)
    s2, o2 = run_v2()
    strict_write(ROOT / 'evidence/manuscript_summary.json', s2)
    s3, o3 = build_v3()
    (EV3 / 'v3_summary.json').write_text(json.dumps(s3, indent=1, allow_nan=False) + '\n')
    for name, text in {**o2, **o3}.items():
        (GEN / name).write_text(text)
    print(json.dumps({'framework_counts': s2['counts'], 'v3_receipts': s3['receipts'], 'macros': len(s3['macros'])}, indent=1))
    return bool(s2['passed'] and s3['passed'])


if __name__ == '__main__':
    sys.exit(0 if run() else 1)
