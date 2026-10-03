#!/usr/bin/env python3
"""Floating-point test of the lemma inequalities of the large-gain proof.

A  referee cycles      lib/lg_referee_cycles.py: cycles at k in {1000, 3000}, mu in {1, 0.3, 0.02}, located with the
                       referee's hybrid integrator (numerical spike, closed-form trough); every lemma of the chain and the
                       four conclusions of the theorem are evaluated on each cycle.
B  author cycles       lib/lg_check_lemmas.py: the same chain evaluated with the author's integrator on cycles at
                       (mu, k) in {(1,1000), (1,3000), (0.5,1000), (0.3,1000), (0.3,3000), (0.1,1000)}.
B2 below the threshold the same script at mu = 1, k in {16, 60, 300}: the lemmas that hold without k >= 1000, and log m_t
                       against the values printed in the paper (-13.98, -67.5, -375.7).
C  author segments     lib/lg_segments.py: the derivation steps on 300 random non-periodic orbit segments.
D  author stage test   lib/lg_stage.py: the Stage-1 and Stage-2 implications on 400 random non-periodic negative lobes.
E  referee segments    lib/lg_referee_segments.py: 300 negative lobes, 300 spikes, and the linear estimate for 100 linear
                       cooperative systems with an arbitrary periodic coefficient x(t).
F  referee stage test  lib/lg_referee_stage.py: 1500 non-periodic deep lobes with k >= 1000.

Counting rule: one count per inequality evaluated.  A comparison passes when it holds within the tolerance stated in the
wrapped script (1e-7 relative for the segment tests; none for the cycle tests of A and B).
A cycle is one orbit; the tests cannot probe the quantifier "every cycle".  The non-periodic tests C-F exercise the
derivation steps in regimes that the computed cycle does not visit.  Nothing here is an interval enclosure.

--small-mu-supplement runs the author's script B at (mu, k) = (0.02, 1000) (about seven minutes) and writes a separate
receipt.  It is not part of the default run: one of its 37 comparisons is below the resolution of the script's
trapezoid quadrature (see the receipt), and the same parameter pair is covered by the referee's integrator in A.
"""
from __future__ import annotations
import argparse, sys, time, warnings
import numpy as np
from common import LIB, Clock, write_receipt, exec_script, capture

AUTHOR_CYCLES = [(1.0, 1000.0), (1.0, 3000.0), (0.5, 1000.0), (0.3, 1000.0), (0.3, 3000.0), (0.1, 1000.0)]
# below the threshold of the theorem: only the lemmas that do not need k >= 1000 are evaluated by the wrapped script;
# the third entry is the value of log m_t printed in the remark on the scope of the threshold
BELOW_THRESHOLD = [(1.0, 16.0, -13.98, 0.006), (1.0, 60.0, -67.5, 0.06), (1.0, 300.0, -375.7, 0.06)]
MAIN = {'__name__': '__main__'}


def author_cycle_records(cases, log_name, group='B author cycles'):
    import lg_check_lemmas as LC
    recs = []
    for i, (mu, k) in enumerate(cases):
        t0 = time.time()
        try:
            (ok, rows, info), _ = capture(LC.analyse, mu, k, verbose=True, return_rows=True, log_name=log_name, append=i > 0)
            failed = [r for r in rows if not r['passed']]
            recs.append(dict(group=group, id=f'mu_{mu:g}_k_{k:g}', mu=mu, k=k, attempted=len(rows),
                             n_passed=len(rows) - len(failed), failed_comparisons=failed, passed=not failed,
                             T=info['T'], T_mu_over_lambda=info['T'] * mu / info['lam'], pk_over_lambda=info['pk'] / info['lam'],
                             log_mt=info['log_mt'], mean_eta_over_mu=1.0 + info['log_mt'] / (info['T'] * mu),
                             lambda_tau_c=info['lam'] * info['tau_c'], seconds=time.time() - t0))
        except Exception as exc:
            recs.append(dict(group=group, id=f'mu_{mu:g}_k_{k:g}', mu=mu, k=k, attempted=1, n_passed=0,
                             passed=False, error=repr(exc), seconds=time.time() - t0))
        print(recs[-1]['id'], 'PASS' if recs[-1]['passed'] else 'FAIL', f"{recs[-1]['seconds']:.1f}s", flush=True)
    return recs


def supplement():
    clock = Clock()
    recs = author_cycle_records([(0.02, 1000.0)], 'v3_largegain_lemmas_small_mu_supplement__author.out')
    r = recs[0]
    f = r['failed_comparisons'][0] if r.get('failed_comparisons') else None
    note = ('Author script at mu = 0.02, k = 1000 (lambda = 50000, T about 3.26e6).  The comparison "C: M+ < k T-" is '
            'equivalent to the positivity of the integral of x over the negative lobe, which is of order k/p_k; M+ and k T- '
            'are of order 3e9 and M+ is obtained by trapezoid quadrature on 400001 nodes, whose relative error (measured by '
            'the mass-identity comparison of the same run, about 7e-12) exceeds the relative gap of the two sides.  The '
            'comparison is below the resolution of the script; it is reported as failed, not relabelled.  The referee '
            'integrator covers this parameter pair in v3_largegain_lemmas (group A).')
    extra = {'note': note, 'part_of_default_run': False}
    if f:
        extra['unresolved_comparison'] = dict(name=f['name'], lhs=f['lhs'], rhs=f['rhs'],
                                              relative_gap=abs(f['lhs'] - f['rhs']) / max(abs(f['lhs']), abs(f['rhs'])))
    return write_receipt('v3_largegain_lemmas_small_mu_supplement', 'float', recs, worst={}, clock=clock,
                         attempted=r['attempted'], n_passed=r['n_passed'], extra=extra,
                         parameters={'mu': 0.02, 'k': 1000.0},
                         scope='One computed cycle; supplementary, retained with its unresolved comparison.')


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--small-mu-supplement', action='store_true'); a = ap.parse_args()
    warnings.filterwarnings('ignore'); np.seterr(all='ignore')
    if a.small_mu_supplement:
        supplement()
        return 0       # the supplement records an unresolved comparison by design; its exit status is not a verdict
    clock = Clock()
    recs = []
    # ---- A
    t0 = time.time()
    ns, _, st = exec_script(LIB / 'lg_referee_cycles.py', init=MAIN, log_name='v3_largegain_lemmas__referee_cycles.out')
    by = {}
    for r in ns['RECORDS']:
        by.setdefault(r['tag'], []).append(r)
    for tag, rows in by.items():
        failed = [r for r in rows if not r['passed']]
        recs.append(dict(group='A referee cycles', id=tag, attempted=len(rows), n_passed=len(rows) - len(failed),
                         failed_comparisons=failed, passed=not failed))
    print('A referee cycles:', sum(r['attempted'] for r in recs), 'checks', f'{time.time()-t0:.1f}s', flush=True)
    # ---- B
    recs += author_cycle_records(AUTHOR_CYCLES, 'v3_largegain_lemmas__author_cycles.out')
    # ---- B2: cycles below the threshold (parameter-free lemmas only), and the printed values of log m_t
    below = author_cycle_records([(m, k_) for m, k_, _, _ in BELOW_THRESHOLD], 'v3_largegain_lemmas__author_cycles_below_threshold.out',
                                 group='B2 author cycles below the threshold')
    for r, (m, k_, printed, tol) in zip(below, BELOW_THRESHOLD):
        if 'log_mt' in r:
            ok = abs(r['log_mt'] - printed) <= tol
            r.update(printed_log_mt=printed, printed_value_reproduced=bool(ok), attempted=r['attempted'] + 1,
                     n_passed=r['n_passed'] + int(ok), passed=bool(r['passed'] and ok))
    recs += below
    # ---- C
    t0 = time.time()
    ns, _, st = exec_script(LIB / 'lg_segments.py', init=MAIN, log_name='v3_largegain_lemmas__author_segments.out')
    slack_c = []
    for name, (cnt, nf, slack) in ns['stats'].items():
        recs.append(dict(group='C author segments', id=name, attempted=int(cnt), n_passed=int(cnt - nf),
                         min_relative_slack=float(slack), passed=nf == 0))
        slack_c.append(float(slack))
    seg_outcomes = dict(ns['outcomes'])
    print('C author segments:', seg_outcomes, f'{time.time()-t0:.1f}s', flush=True)
    # ---- D
    t0 = time.time()
    ns, _, st = exec_script(LIB / 'lg_stage.py', init=MAIN, log_name='v3_largegain_lemmas__author_stage.out')
    recs.append(dict(group='D author stage test', id='Stage 1: T_- >= 5.89 implies p_k >= 0.374 k', attempted=int(ns['n1']),
                     n_passed=int(ns['n1'] - ns['f1']), min_pk_over_k=float(ns['worst1']), passed=ns['f1'] == 0 and ns['n1'] > 0))
    recs.append(dict(group='D author stage test', id='Stage 2: T_- >= 22 implies r_ell >= 0.3 T_-', attempted=int(ns['n2']),
                     n_passed=int(ns['n2'] - ns['f2']), min_rl_over_Tminus=float(ns['worst2']), passed=ns['f2'] == 0 and ns['n2'] > 0))
    recs.append(dict(group='D author stage test', id='strip time tau_b < 0.36 on all segments', attempted=1,
                     n_passed=int(ns['maxtb'] < 0.36), max_tau_b=float(ns['maxtb']), passed=bool(ns['maxtb'] < 0.36)))
    stage_outcomes = dict(ns['out'])
    print('D author stage:', stage_outcomes, f'{time.time()-t0:.1f}s', flush=True)
    # ---- E
    t0 = time.time()
    ns, text, st = exec_script(LIB / 'lg_referee_segments.py', init=MAIN, log_name='v3_largegain_lemmas__referee_segments.out')
    worst_e = {n: float(v[0]) for n, v in ns['worst'].items()}
    recs.append(dict(group='E referee segments', id='all inequality families', attempted=int(ns['ntest']),
                     n_passed=int(ns['ntest'] - ns['nfail']), tightest_relative_slack_by_family=worst_e,
                     passed=ns['nfail'] == 0 and ns['ntest'] > 0))
    print('E referee segments:', ns['ntest'], 'tests', ns['nfail'], 'failures', f'{time.time()-t0:.1f}s', flush=True)
    # ---- F
    t0 = time.time()
    ns, _, st = exec_script(LIB / 'lg_referee_stage.py', argv=['1500'], init=MAIN, log_name='v3_largegain_lemmas__referee_stage.out')
    s_, m_ = ns['stats'], ns['mins']
    nF = int(s_['lobes'] + s_['hyp1'] + s_['hyp2'])
    recs.append(dict(group='F referee stage test', id='tau_b < 0.36; Stage-1 and Stage-2 conclusions under their hypotheses',
                     attempted=nF, n_passed=nF - int(s_['fail']), lobes=int(s_['lobes']), stage1_hypotheses_met=int(s_['hyp1']),
                     stage2_hypotheses_met=int(s_['hyp2']), min_pk_over_k=float(m_['pk_over_k']), min_rl=float(m_['rl']),
                     min_rl_over_Tminus=float(m_['rl_over_T']), max_tau_b=float(m_['taub']), passed=s_['fail'] == 0 and s_['hyp1'] > 0))
    print('F referee stage:', dict(s_), f'{time.time()-t0:.1f}s', flush=True)

    attempted = sum(r['attempted'] for r in recs)
    n_passed = sum(r['n_passed'] for r in recs)
    cyc = [r for r in recs if r['group'] == 'B author cycles' and 'T' in r]
    worst = {'max_mean_eta_over_mu_on_cycles': max(r['mean_eta_over_mu'] for r in cyc) if cyc else None,
             'min_T_mu_over_lambda_on_cycles': min(r['T_mu_over_lambda'] for r in cyc) if cyc else None,
             'min_pk_over_lambda_on_cycles': min(r['pk_over_lambda'] for r in cyc) if cyc else None,
             'min_relative_slack_author_segments': min(slack_c) if slack_c else None,
             'min_relative_slack_referee_segments': min(worst_e.values()),
             'log_mt_below_threshold': {r['id']: r.get('log_mt') for r in below},
             'referee_segment_note': 'negative slacks of the referee segment tests are equality cases within the 1e-7 relative tolerance of the wrapped script',
             'min_pk_over_k_stage1_nonperiodic': min(float(ns['mins']['pk_over_k']), [r for r in recs if r['id'].startswith('Stage 1')][0]['min_pk_over_k'])}
    groups = {}
    for r in recs:
        g = groups.setdefault(r['group'], {'attempted': 0, 'passed': 0})
        g['attempted'] += r['attempted']; g['passed'] += r['n_passed']
    extra = {'groups': groups, 'author_segment_outcomes': seg_outcomes, 'author_stage_outcomes': stage_outcomes,
             'counting_rule': 'one count per inequality evaluated',
             'excluded_from_default_run': 'author script at (mu, k) = (0.02, 1000): see v3_largegain_lemmas_small_mu_supplement.json'}
    ok = write_receipt('v3_largegain_lemmas', 'float', recs, worst=worst, clock=clock, extra=extra, attempted=attempted,
                       n_passed=n_passed, passed=bool(recs) and all(r['passed'] for r in recs),
                       parameters={'referee_cycles': 'k in {1000, 3000} x mu in {1, 0.3, 0.02}', 'author_cycles': AUTHOR_CYCLES,
                                   'segments': 'random non-periodic data; seeds fixed in the wrapped scripts'},
                       scope='Inequalities evaluated on computed orbits. Not an enclosure; one cycle per parameter pair.')
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
