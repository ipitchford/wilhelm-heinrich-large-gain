#!/usr/bin/env python3
"""Tests of the kurtosis bound, of the tau-time period identities, and of the third-order selection formula.

A  kurtosis bound      lib/theoremA.py (the referee's script; mpmath, 30 digits, tanh-sinh quadrature in U, no ODE):
                       on 14 levels h in [1e-6, 1e6] of the potential e^U - U - 1: <P^4> <= (9/5) <P^2>^2, the first
                       inequality of the mixture argument, <P^2> >= P_max^2/3 and 4h/5 <= Rsel <= 4h/3; Rsel' from central
                       differences against (2/3)(3 - kurtosis), J_0' = T_0, J_1' = 3 J_0; Rsel' >= 4/5; Rsel increasing on 121
                       levels; ten further convex potentials (the bound must hold); two non-convex potentials (the bound
                       must be violated, which shows that the test can fail); monotonicity of the velocity density.
B  period identities   lib/wedge_identities.py: the seven identities (a)-(g) on 14 cycles computed with an independent
                       integrator, mu in {0.02, 0.3, 1} x rho in {1.01, 1.5, 4, 12}, (3, 2), (0.02, 1.003); each defect is
                       normalised by the largest term entering it; pass: defect <= 1e-8.
C  third-order formula amplitude-fixed cycles (lib/ampcycle.py) at crest amplitudes 0.3, 1, 2, 3 and
                       mu in {3e-2, 1e-2, 3e-3, 1e-3, 3e-4, 1e-4}: (rho-1)/eps^2 against Rsel(h) and log m_t/eps^3 against
                       T_0 (2/3)(<P^4> - 3<P^2>^2)/<P^2> = -J_0 Rsel'.  Pass, fixed in advance: every cycle is located
                       (residual <= 1e-9) and attracting; the two relative errors decrease along the mu sequence; the
                       error ratio between mu = 1e-3 and mu = 1e-4 lies in [5, 20] (second order in eps predicts 9.98);
                       at mu = 1e-4 the errors are below 5e-4 and 1e-3 for amplitudes up to 1, and below 1e-2 otherwise.
D  fixed rho           the six cycles at rho in {1.5, 3}, mu in {1e-1, 1e-2, 1e-3} with the bundle integrator: m_t in (0, 1),
                       the identity <(rho e^U - 1) p^2> = 0 to 1e-6 relative, and the values U_max, U_min, eps^2 R(U_max),
                       mu U_min, mu T, m_t against the numbers printed in the table of the wedge section (half a unit of the
                       last printed digit, with a 20 per cent allowance).
The rows of the selection table printed in the wedge section are compared in the same way with the ratios computed in C.
Part A is a high-precision quadrature, parts B, C and D are floating point.  Nothing here is an interval enclosure, and the
agreement in C is evidence for the asymptotic formula, not a proof of the wedge theorem (which cites Fenichel's theorem).
"""
from __future__ import annotations
import sys, time, warnings
import numpy as np
from common import LIB, Clock, write_receipt, exec_script

AMPS = [0.3, 1.0, 2.0, 3.0]
MUS = [3e-2, 1e-2, 3e-3, 1e-3, 3e-4, 1e-4]
# rows printed in the wedge section: (amplitude, mu) -> (selection ratio, exponent ratio), five decimals
PRINTED_SELECTION = {(0.3, 1e-2): (0.96373, 0.95298), (0.3, 1e-3): (0.99619, 0.99501), (0.3, 1e-4): (0.99962, 0.99950),
                     (1.0, 1e-3): (0.99599, 0.99261), (1.0, 1e-4): (0.99960, 0.99925), (2.0, 1e-3): (0.99302, 0.97829),
                     (2.0, 1e-4): (0.99929, 0.99775), (3.0, 1e-3): (0.98318, 0.93877), (3.0, 1e-4): (0.99821, 0.99329)}
# (rho, mu) -> printed (U_max, U_min, eps^2 R(U_max), mu U_min, mu T, m_t) with the number of printed decimals
PRINTED_FIXED_RHO = {(1.5, 1e-1): (2.553, -13.5, 0.768, -1.352, 3.1928, 0.5692), (1.5, 1e-2): (4.480, -104.5, 0.811, -1.045, 2.5209, 0.6669),
                     (1.5, 1e-3): (6.737, -1032.5, 0.833, -1.032, 2.4470, 0.6785), (3.0, 1e-1): (4.093, -96.4, 4.531, -9.639, 5.9425, 0.1379),
                     (3.0, 1e-2): (6.195, -801.7, 4.736, -8.017, 5.3313, 0.1810), (3.0, 1e-3): (8.476, -7882.3, 4.779, -7.882, 5.2685, 0.1861)}
FIXED_RHO_DECIMALS = (3, 1, 3, 3, 4, 4)


def part_A(recs):
    t0 = time.time()
    ns, text, st = exec_script(LIB / 'theoremA.py', log_name='v3_wedge__theoremA.out')
    fail = list(ns['FAIL'])
    def add(ident, ok, **kw):
        recs.append(dict(group='A kurtosis bound', id=ident, passed=bool(ok), **kw))
    add('polynomial identities of the mixture argument', 'T1' not in fail)
    for row in ns['TODA_ROWS']:
        h = row['h']
        add(f'Toda level h={h}: kurtosis <= 9/5 and companions', f'T2 h={h}' not in fail, kurtosis=row['kurtosis'], Rsel_over_h=row['Rsel_over_h'])
        add(f'Toda level h={h}: derivative formulas', f'T2-derivative h={h}' not in fail)
    add("Rsel' >= 4/5 on the level grid", "T2 Rsel' < 4/5" not in fail, min_Rsel_prime=str(ns['minRp']))
    add('Rsel strictly increasing on 121 levels', 'T2 monotone' not in fail)
    for row in ns['OTHER_ROWS']:
        add(f"convex potential {row['potential']}: kurtosis <= 9/5", ('T3 ' + row['potential']) not in fail, max_kurtosis=row['max_kurtosis'])
    for row in ns['NONCONVEX_ROWS']:
        add(f"non-convex potential {row['potential']}: bound violated (power check)", row['exceeds'], max_kurtosis=row['max_kurtosis'])
    add('velocity density increasing (three levels, 199 grid points each)', 'T5' not in fail)
    kurt = [float(r['kurtosis']) for r in ns['TODA_ROWS']]
    other = [float(r['max_kurtosis']) for r in ns['OTHER_ROWS']]
    print('A kurtosis:', 'failures', fail, f'{time.time()-t0:.1f}s', flush=True)
    return dict(max_kurtosis_toda=max(kurt), max_kurtosis_other_convex=max(other), min_Rsel_prime=float(ns['minRp']),
                script_failures=fail, exit_status=st)


def part_B(recs):
    t0 = time.time()
    ns, text, st = exec_script(LIB / 'wedge_identities.py', log_name='v3_wedge__identities.out')
    worst = 0.0
    for row in ns['ROWS']:
        ident = f"cycle mu={row['mu']:g} rho={row['rho']:g}"
        if 'error' in row:
            recs.append(dict(group='B period identities', id=ident, passed=False, error=row['error']))
            continue
        for name, v in row['defects'].items():
            recs.append(dict(group='B period identities', id=f'{ident}: identity ({name})', passed=bool(v <= 1e-8), defect=v))
            worst = max(worst, v)
    print('B identities:', len(ns['ROWS']), 'cycles, worst defect', worst, f'{time.time()-t0:.1f}s', flush=True)
    return dict(max_identity_defect=worst, cycles=len(ns['ROWS']),
                min_U=min((r['U_min'] for r in ns['ROWS'] if 'U_min' in r), default=None))


def part_C(recs):
    t0 = time.time()
    from ampcycle import amp_cycle, eps_of_mu
    from toda import level_data, Phi
    out = {}
    for a in AMPS:
        h = float(Phi(a))
        d = level_data(h)
        r_pred = (2.0 / 3.0) * d['J1'] / d['I']
        m_pred = d['T0'] * (2.0 / 3.0) * (d['p4'] - 3 * d['p2'] ** 2) / d['p2']
        er, em = [], []
        for mu in MUS:
            eps = eps_of_mu(mu)
            ident = f'amplitude {a:g}, mu {mu:g}'
            try:
                c = amp_cycle(mu, a)
            except Exception as exc:
                recs.append(dict(group='C third-order formula', id=ident + ': cycle located', passed=False, error=repr(exc)))
                er.append(np.nan); em.append(np.nan)
                continue
            r_num = (c['rho'] - 1) / eps ** 2
            m_num = c['logmt'] / eps ** 3
            er.append(abs(r_num / r_pred - 1)); em.append(abs(m_num / m_pred - 1))
            if (a, mu) in PRINTED_SELECTION:
                ps, pm = PRINTED_SELECTION[(a, mu)]
                recs.append(dict(group='C third-order formula', id=ident + ': printed table row reproduced',
                                 passed=bool(abs(r_num / r_pred - ps) <= 6e-6 and abs(m_num / m_pred - pm) <= 6e-6),
                                 computed=[r_num / r_pred, m_num / m_pred], printed=[ps, pm]))
            recs.append(dict(group='C third-order formula', id=ident + ': cycle located and attracting',
                             passed=bool(c['residual'] <= 1e-9 and c['logmt'] < 0), rho=c['rho'], residual=c['residual'],
                             selection_numeric=r_num, selection_predicted=r_pred, exponent_numeric=m_num, exponent_predicted=m_pred,
                             selection_relative_error=er[-1], exponent_relative_error=em[-1], log_mt=c['logmt']))
        er, em = np.array(er), np.array(em)
        cap_r, cap_m = (5e-4, 1e-3) if a <= 1.0 else (1e-2, 1e-2)
        for name, e, cap in (('selection', er, cap_r), ('exponent', em, cap_m)):
            recs.append(dict(group='C third-order formula', id=f'amplitude {a:g}: {name} error decreasing in mu',
                             passed=bool(np.all(np.isfinite(e)) and np.all(np.diff(e) < 0)), errors=e.tolist()))
            ratio = float(e[3] / e[5])
            recs.append(dict(group='C third-order formula', id=f'amplitude {a:g}: {name} error ratio mu=1e-3 over mu=1e-4 in [5, 20]',
                             passed=bool(5.0 <= ratio <= 20.0), ratio=ratio, predicted_for_second_order=(eps_of_mu(1e-3) / eps_of_mu(1e-4)) ** 2))
            recs.append(dict(group='C third-order formula', id=f'amplitude {a:g}: {name} error at mu=1e-4 below {cap:g}',
                             passed=bool(e[5] <= cap), error=float(e[5])))
        out[a] = (er, em)
    print('C third order:', {a: (float(v[0][-1]), float(v[1][-1])) for a, v in out.items()}, f'{time.time()-t0:.1f}s', flush=True)
    return {'max_selection_error_smallest_mu_amplitude_le_1': float(max(out[a][0][-1] for a in AMPS if a <= 1)),
            'max_exponent_error_smallest_mu_amplitude_le_1': float(max(out[a][1][-1] for a in AMPS if a <= 1)),
            'max_selection_error_smallest_mu': float(max(out[a][0][-1] for a in AMPS)),
            'max_exponent_error_smallest_mu': float(max(out[a][1][-1] for a in AMPS))}


def part_D(recs):
    t0 = time.time()
    import math
    import cyclekit as ck
    worst = 0.0
    for (rho, mu), printed in PRINTED_FIXED_RHO.items():
        ident = f'fixed rho={rho:g}, mu={mu:g}'
        try:
            k = rho * ck.kH(mu); kh = ck.kH(mu)
            cyc = ck.find(mu, rho)
            pr = ck.profile(mu, k, cyc, extras=lambda d: [(d['x'] - kh) * d['p']**2, d['p']**2], n_extras=2)
        except Exception as exc:
            recs.append(dict(group='D fixed rho', id=ident, passed=False, error=repr(exc)))
            continue
        T, lam = pr['T'], float(pr['means'][0])
        Umax, Umin = cyc['U0'], float(pr['g']['U'].min())
        eps2 = mu / (1 + mu)**2
        vals = (Umax, Umin, eps2 * (math.expm1(Umax) - Umax), mu * Umin, mu * T, math.exp(lam * T))
        defect = abs(float(pr['extra_means'][0])) / (kh * float(pr['extra_means'][1]))
        worst = max(worst, defect)
        ok_print = all(abs(v - q) <= 0.6 * 10.0**(-n) for v, q, n in zip(vals, printed, FIXED_RHO_DECIMALS))
        recs.append(dict(group='D fixed rho', id=ident + ': 0 < m_t < 1 and the period identity', m_t=vals[5], identity_defect=defect,
                         passed=bool(0 < vals[5] < 1 and defect <= 1e-6 and cyc['res_rel'] <= 1e-9)))
        recs.append(dict(group='D fixed rho', id=ident + ': printed table row reproduced', computed=list(vals), printed=list(printed),
                         passed=bool(ok_print)))
    print('D fixed rho: worst identity defect', worst, f'{time.time()-t0:.1f}s', flush=True)
    return {'max_fixed_rho_identity_defect': worst}


def main():
    warnings.filterwarnings('ignore'); np.seterr(all='ignore')
    clock = Clock()
    recs = []
    worst = {}
    for part in (part_A, part_B, part_C, part_D):
        try:
            worst.update(part(recs))
        except Exception as exc:
            recs.append(dict(group=part.__name__, id='exception', passed=False, error=repr(exc)))
    groups = {}
    for r in recs:
        g = groups.setdefault(r['group'], {'attempted': 0, 'passed': 0})
        g['attempted'] += 1; g['passed'] += int(r['passed'])
    extra = {'groups': groups, 'precision': {'A': 'mpmath, 30 significant digits', 'B': 'double precision', 'C': 'double precision', 'D': 'double precision'},
             'counting_rule': 'A: one count per level, potential or structural test; B: one per (cycle, identity); C: one per cycle, per convergence predicate and per printed table row; D: two per cycle'}
    ok = write_receipt('v3_wedge', 'float', recs, worst=worst, clock=clock, extra=extra,
                       parameters={'toda_levels': '14 values of h in [1e-6, 1e6]', 'identity_cycles': 14, 'amplitudes': AMPS, 'mu': MUS},
                       scope='Quadratures on level curves and computed cycles. The asymptotic comparison is evidence, not a proof.')
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
