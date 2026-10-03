#!/usr/bin/env python3
"""Cartography: consistency of the shipped tables, and recomputation of a stated subsample.

The tables under evidence/v3/cartography/ were produced during the research phase by a compiled integrator
(lib/cart.py, numba) and are shipped unchanged; README.md there lists files and columns.  This script does three things.

T   statements read off the shipped tables (always run; no integration):
      grid_table.csv     300 cycles = 12 mu x 25 rho; every lambda_c < 0; log m_t and lambda_c strictly decreasing in rho at
                         each mu; <U> + log rho < 0; Newton residuals
      branch_table.csv   827 amplitude-branch cycles; every lambda_c < 0; rho strictly increasing and log m_t strictly
                         decreasing in the crest amplitude on each of the 12 branches (no fold)
      hunt_results.jsonl 240 random parameter points; 239 computed, one failure in direct coordinates redone at the reciprocal
                         rate (hunt_results_supplement.jsonl); every lambda_c < 0; no return-map pair with the wrong sign of
                         the displacement; no long run converging to a different crest; no new Newton fixed point
      limit_map.csv      the formal mu -> 0 limit map: Lambda strictly increasing, 0 < m_0 < 1
R1  recomputation of grid cycles with the bundle's own integrator (cyclekit: scipy DOP853, logarithmic variables, Newton on the
    crest map with variational equations), started from the shipped crest.  Subsample, fixed in advance: rho indices
    0, 6, 12, 18 for every mu (48 cycles, rho - 1 from 1e-3 to 9.4) and rho index 24 (rho = 200) for mu in {0.2, 0.5, 1}.
    Compared with the shipped row: T (1e-7 relative), crest U (1e-7), trough U (1e-6 relative), log m_t (1e-6 (1+|log m_t|)).
R2  recomputation of uniqueness hunts at the shipped points with index 0, 10, ..., 230 (24 points), with the research code
    (lib/cart_hunt.py; time-boxed, so the number of pairs varies from run to run): primary cycle against the shipped line;
    lambda_c < 0; no wrong-sign displacement pair; no long run ending elsewhere; no new fixed point; rho increasing along
    the amplitude branch.  The primary cycle is also recomputed with the bundle's own integrator when rho <= 100, at the
    reciprocal rate when mu > 100 and rho > 5 (T' = mu T, same log m_t).  A shipped failure line is rerun at the reciprocal
    rate, as in the research phase, and compared with the supplement.

--tables-only runs T alone and writes v3_cartography_tables_only.json (used by run_all.py --quick).
Floating point throughout.  A finite search that finds no second cycle does not exclude one.
"""
from __future__ import annotations
import argparse, csv, json, math, sys, time, warnings
import numpy as np
from common import EVID, Clock, write_receipt, capture

CART = EVID / 'cartography'
GRID_RHO_IDX = (0, 6, 12, 18)
GRID_TOP_MU = (0.2, 0.5, 1.0)
HUNT_IDX = tuple(range(0, 240, 10))
HUNT_BUDGET = 8.0


def _conv(v):
    if v in ('', None):
        return None
    if v in ('True', 'False'):
        return v == 'True'
    try:
        return float(v)
    except ValueError:
        return v


def rows_csv(name):
    with (CART / name).open() as f:
        return [{k: _conv(v) for k, v in r.items()} for r in csv.DictReader(f)]


def rows_jsonl(name):
    return [json.loads(l) for l in (CART / name).read_text().splitlines() if l.strip()]


def by_mu(rows):
    d = {}
    for r in rows:
        d.setdefault(r['mu'], []).append(r)
    return d


def part_T(recs):
    def add(ident, ok, **kw):
        recs.append(dict(group='T shipped tables', id=ident, passed=bool(ok), **kw))
    out = {}
    grid = rows_csv('grid_table.csv')
    g = by_mu(grid)
    add('grid: 300 rows on 12 mu x 25 rho', len(grid) == 300 and len(g) == 12 and all(len(v) == 25 for v in g.values()))
    add('grid: lambda_c < 0 on every cycle', all(r['lambda_c'] < 0 for r in grid), cycles=len(grid))
    worst = max(grid, key=lambda r: r['lambda_c'])
    out.update(grid_cycles=len(grid), grid_max_lambda_c=worst['lambda_c'], grid_max_lambda_c_at={'mu': worst['mu'], 'rho': worst['rho'], 'log_mt': worst['log_mt']},
               grid_max_log_mt=max(r['log_mt'] for r in grid), grid_min_log_mt=min(r['log_mt'] for r in grid),
               grid_mu_range=[min(g), max(g)], grid_rho_range=[min(r['rho'] for r in grid), max(r['rho'] for r in grid)],
               grid_min_U=min(r['U_min'] for r in grid), grid_max_newton_residual=max(r['newton_residual'] for r in grid))
    for mu, v in sorted(g.items()):
        v = sorted(v, key=lambda r: r['rho'])
        add(f'grid mu={mu:g}: log m_t strictly decreasing in rho', all(b['log_mt'] < a['log_mt'] for a, b in zip(v, v[1:])))
        add(f'grid mu={mu:g}: lambda_c strictly decreasing in rho', all(b['lambda_c'] < a['lambda_c'] for a, b in zip(v, v[1:])))
    add('grid: <U> + log rho < 0 on every cycle', all(r['meanU_plus_log_rho'] < 0 for r in grid),
        max_value=max(r['meanU_plus_log_rho'] for r in grid))
    add('grid: Newton residual <= 1e-9 on every cycle', out['grid_max_newton_residual'] <= 1e-9)
    branch = rows_csv('branch_table.csv')
    b = by_mu(branch)
    add('branch: 827 rows on 12 branches', len(branch) == 827 and len(b) == 12)
    add('branch: lambda_c < 0 on every cycle', all(r['lambda_c'] < 0 for r in branch), cycles=len(branch))
    out.update(branch_cycles=len(branch), branch_max_lambda_c=max(r['lambda_c'] for r in branch), branch_max_rho=max(r['rho'] for r in branch))
    for mu, v in sorted(b.items()):
        v = sorted(v, key=lambda r: r['U_max'])
        add(f'branch mu={mu:g}: rho strictly increasing in crest amplitude (no fold)', all(q['rho'] > p_['rho'] for p_, q in zip(v, v[1:])))
        add(f'branch mu={mu:g}: log m_t strictly decreasing in crest amplitude', all(q['log_mt'] < p_['log_mt'] for p_, q in zip(v, v[1:])))
    hunt = rows_jsonl('hunt_results.jsonl')
    good = [h for h in hunt if 'error' not in h]
    supp = rows_jsonl('hunt_results_supplement.jsonl')
    add('hunt: 240 points, 239 computed, 1 failure line', len(hunt) == 240 and len(good) == 239)
    add('hunt: the failure is covered by the reciprocal-rate supplement', len(supp) == 1 and 'error' not in supp[0]
        and all(abs(supp[0]['rho'] - h['rho']) < 1e-9 * h['rho'] and abs(supp[0]['mu'] * h['mu'] - 1) < 1e-9 for h in hunt if 'error' in h))
    allh = good + supp
    add('hunt: lambda_c < 0 at every computed point', all(h['lambda_c'] < 0 for h in allh), points=len(allh))
    add('hunt: no wrong-sign displacement pair', sum(h['scan_wrong_sign'] for h in allh) == 0,
        pairs=int(sum(h['scan_pairs'] for h in good)), pairs_with_supplement=int(sum(h['scan_pairs'] for h in allh)))
    add('hunt: no long run converging to a different crest', sum(h['long_different'] for h in allh) == 0 and sum(h['long_heading_bad'] for h in allh) == 0,
        long_runs=int(sum(len(h['long']) for h in good)), reached_primary=int(sum(h['long_same'] for h in good)),
        not_finished=int(sum(h['long_other'] for h in good)))
    add('hunt: no new Newton fixed point', sum(len(h['newton']['new']) for h in allh) == 0,
        newton_starts=int(sum(h['newton']['same'] + h['newton']['equilibrium'] + h['newton']['failed'] + len(h['newton']['new']) for h in good)))
    add('hunt: rho increasing along every amplitude branch', all(h['branch_monotone'] for h in allh))
    out.update(hunt_points=len(hunt), hunt_points_computed_directly=len(good), hunt_max_log_mt=max(h['log_mt'] for h in allh), hunt_mu_range=[min(h['mu'] for h in hunt), max(h['mu'] for h in hunt)],
               hunt_rho_range=[min(h['rho'] for h in hunt), max(h['rho'] for h in hunt)], hunt_pairs=int(sum(h['scan_pairs'] for h in good)),
               hunt_long_runs=int(sum(len(h['long']) for h in good)), hunt_long_not_finished=int(sum(h['long_other'] for h in good)),
               hunt_newton_starts=int(sum(h['newton']['same'] + h['newton']['equilibrium'] + h['newton']['failed'] for h in good)))
    ct = rows_csv('cycle_tests.csv')
    names = ('R1', 'R2', 'R3', 'R4')
    both = [all(r[f'{n}_minus_A_min'] < 0 < r[f'{n}_minus_A_max'] for n in names) for r in ct]
    add('cycle tests: each of the four pointwise candidates takes both signs on every grid cycle', len(ct) == 300 and all(both),
        cycles=int(sum(both)))
    add('cycle tests: projected cycle convex on every grid cycle', len(ct) == 300 and all(r['convex_projection'] is True for r in ct))
    st = rows_csv('stests_summary.csv')
    out.update(cycle_tests_cycles_with_negative_Ehat=int(sum(r['Ehat_min_over_k'] < 0 for r in ct)),
               cycle_tests_cycles_with_zeta_over_eta_on_both_sides_of_1=int(sum(r['min_zeta_over_eta'] < 1 < r['max_zeta_over_eta'] for r in ct)),
               stests_parameter_pairs=len(st))
    lim = sorted(rows_csv('limit_map.csv'), key=lambda r: r['p_in'])
    add('limit map: Lambda strictly increasing', all(q['Lambda=rho_limit'] > p_['Lambda=rho_limit'] for p_, q in zip(lim, lim[1:])), rows=len(lim))
    add('limit map: 0 < m_0 < 1', all(0 < r['m0_limit_multiplier'] < 1 for r in lim))
    out.update(limit_map_rho_range=[lim[0]['Lambda=rho_limit'], lim[-1]['Lambda=rho_limit']])
    return out


def part_R1(recs):
    import cyclekit as ck
    grid = by_mu(rows_csv('grid_table.csv'))
    cases = []
    for mu, v in sorted(grid.items()):
        v = sorted(v, key=lambda r: r['rho'])
        idx = list(GRID_RHO_IDX) + ([24] if mu in GRID_TOP_MU else [])
        cases += [(mu, i, v[i]) for i in idx]
    worst = dict(T=0.0, U_max=0.0, U_min=0.0, log_mt=0.0)
    lam_max = -math.inf
    for mu, i, row in cases:
        t0 = time.time()
        ident = f'grid mu={mu:g} rho index {i} (rho-1={row["rho"]-1:.4g})'
        try:
            k = row['k']
            cyc = ck.from_crest(mu, k, row['U_max'], row['y0_over_k'] * k)
            pr = ck.profile(mu, k, cyc)
            T, lam = pr['T'], float(pr['means'][0])
            d = dict(T=abs(T - row['T']) / row['T'], U_max=abs(cyc['U0'] - row['U_max']) / (1 + abs(row['U_max'])),
                     U_min=abs(float(pr['g']['U'].min()) - row['U_min']) / (1 + abs(row['U_min'])),
                     log_mt=abs(lam * T - row['log_mt']) / (1 + abs(row['log_mt'])))
            ok = d['T'] <= 1e-7 and d['U_max'] <= 1e-7 and d['U_min'] <= 1e-6 and d['log_mt'] <= 1e-6 and lam < 0 and cyc['res_rel'] <= 1e-9
            for n in worst:
                worst[n] = max(worst[n], d[n])
            lam_max = max(lam_max, lam)
            recs.append(dict(group='R1 grid recomputation', id=ident, mu=mu, rho=row['rho'], passed=bool(ok), differences=d,
                             lambda_c=lam, log_mt=lam * T, T=T, shipped_log_mt=row['log_mt'], location_residual=cyc['res_rel'],
                             seconds=time.time() - t0))
        except Exception as exc:
            recs.append(dict(group='R1 grid recomputation', id=ident, mu=mu, rho=row['rho'], passed=False, error=repr(exc), seconds=time.time() - t0))
        print(ident, 'PASS' if recs[-1]['passed'] else 'FAIL', recs[-1].get('differences', recs[-1].get('error')), f"{recs[-1]['seconds']:.1f}s", flush=True)
    return dict(recomputed_grid_cycles=len(cases), recomputed_grid_max_rel_diff_T=worst['T'], recomputed_grid_max_diff_U_max=worst['U_max'],
                recomputed_grid_max_rel_diff_U_min=worst['U_min'], recomputed_grid_max_diff_log_mt=worst['log_mt'],
                recomputed_grid_max_lambda_c=lam_max)


def independent_primary(mu, rho, line):
    """primary cycle with the bundle's own integrator; returns (status, differences)"""
    import cyclekit as ck
    if rho > 100:
        return 'not attempted (rule: rho > 100)', None
    if mu > 100 and rho > 5:
        m2 = 1.0 / mu
        cyc = ck.find(m2, rho)
        pr = ck.profile(m2, rho * ck.kH(m2), cyc)
        T, lam = pr['T'] / mu, float(pr['means'][0]) * mu          # T' = mu T,  lambda_c' T' = lambda_c T
        how = 'reciprocal rate'
    else:
        k = line['k']
        cyc = ck.from_crest(mu, k, line['U_star'], line['y_star_over_k'] * k)
        pr = ck.profile(mu, k, cyc)
        T, lam = pr['T'], float(pr['means'][0])
        how = 'direct'
    return how, dict(T=abs(T - line['T']) / line['T'], log_mt=abs(lam * T - line['log_mt']) / (1 + abs(line['log_mt'])))


def part_R2(recs):
    hunt = rows_jsonl('hunt_results.jsonl')
    supp = rows_jsonl('hunt_results_supplement.jsonl')
    try:
        import cart_hunt as H
    except Exception as exc:
        recs.append(dict(group='R2 hunt recomputation', id='import of the compiled research code (needs numba)', passed=False, error=repr(exc)))
        return {}
    tot = dict(pairs=0, wrong=0, long=0, long_same=0, long_other=0, different=0, newton=0, new=0)
    worst = dict(U_star=0.0, T=0.0, log_mt=0.0, ind_T=0.0, ind_log_mt=0.0)
    n_ind = 0
    lam_max = -math.inf
    for i in HUNT_IDX:
        t0 = time.time()
        line = hunt[i]
        via = 'direct coordinates'
        if 'error' in line:
            line = supp[0]
            via = 'reciprocal rate (shipped line is a failure in direct coordinates)'
        mu, rho = line['mu'], line['rho']
        ident = f'hunt point {i}: mu={hunt[i]["mu"]:.5g} rho={rho:.6g}'
        try:
            r, _ = capture(H.hunt_point, mu, rho, np.random.default_rng(20261003 + i), budget_s=HUNT_BUDGET)
            d = dict(U_star=abs(r['U_star'] - line['U_star']) / (1 + abs(line['U_star'])), T=abs(r['T'] - line['T']) / line['T'],
                     log_mt=abs(r['log_mt'] - line['log_mt']) / (1 + abs(line['log_mt'])))
            how, ind = independent_primary(mu, rho, line)
            checks = {
                'primary_cycle_matches_shipped': d['U_star'] <= 1e-8 and d['T'] <= 1e-8 and d['log_mt'] <= 1e-7,
                'attracting': r['lambda_c'] < 0,
                'no_wrong_sign_pair': r['scan_pairs'] > 0 and r['scan_wrong_sign'] == 0,
                'no_long_run_elsewhere': r['long_different'] == 0 and r['long_heading_bad'] == 0,
                'no_new_fixed_point': len(r['newton']['new']) == 0,
                'branch_monotone': bool(r['branch_monotone']),
            }
            if ind is not None:
                checks['independent_integrator'] = ind['T'] <= 1e-7 and ind['log_mt'] <= 1e-6
                worst['ind_T'] = max(worst['ind_T'], ind['T']); worst['ind_log_mt'] = max(worst['ind_log_mt'], ind['log_mt'])
                n_ind += 1
            for n in ('U_star', 'T', 'log_mt'):
                worst[n] = max(worst[n], d[n])
            nn = r['newton']
            tot['pairs'] += r['scan_pairs']; tot['wrong'] += r['scan_wrong_sign']; tot['long'] += len(r['long'])
            tot['long_same'] += r['long_same']; tot['long_other'] += r['long_other']; tot['different'] += r['long_different']
            tot['newton'] += nn['same'] + nn['equilibrium'] + nn['failed'] + len(nn['new']); tot['new'] += len(nn['new'])
            lam_max = max(lam_max, r['lambda_c'])
            recs.append(dict(group='R2 hunt recomputation', id=ident, index=i, mu=mu, rho=rho, via=via, passed=all(checks.values()),
                             checks=checks, differences_from_shipped=d, independent_integrator=how, independent_differences=ind,
                             lambda_c=r['lambda_c'], log_mt=r['log_mt'], scan_pairs=r['scan_pairs'], scan_wrong_sign=r['scan_wrong_sign'],
                             long_runs=[{k_: v for k_, v in q.items() if k_ in ('ic', 'nev', 'status', 'gap_after_iter')} for q in r['long']],
                             newton={'same': nn['same'], 'equilibrium': nn['equilibrium'], 'failed': nn['failed'], 'new': len(nn['new'])},
                             seconds=time.time() - t0))
        except Exception as exc:
            recs.append(dict(group='R2 hunt recomputation', id=ident, index=i, mu=mu, rho=rho, via=via, passed=False, error=repr(exc),
                             seconds=time.time() - t0))
        q = recs[-1]
        print(ident, 'PASS' if q['passed'] else 'FAIL', {k_: v for k_, v in q.get('checks', {}).items() if not v}, q.get('error', ''),
              q.get('independent_integrator'), f"{q['seconds']:.1f}s", flush=True)
    return dict(recomputed_hunts=len(HUNT_IDX), recomputed_hunt_pairs=tot['pairs'], recomputed_hunt_wrong_sign_pairs=tot['wrong'],
                recomputed_hunt_long_runs=tot['long'], recomputed_hunt_long_runs_reaching_primary=tot['long_same'],
                recomputed_hunt_long_runs_not_finished=tot['long_other'], recomputed_hunt_long_runs_elsewhere=tot['different'],
                recomputed_hunt_newton_starts=tot['newton'], recomputed_hunt_new_fixed_points=tot['new'],
                recomputed_hunt_max_lambda_c=lam_max, recomputed_hunt_max_diff_U_star=worst['U_star'],
                recomputed_hunt_max_rel_diff_T=worst['T'], recomputed_hunt_max_diff_log_mt=worst['log_mt'],
                independent_primary_checks=n_ind, independent_max_rel_diff_T=worst['ind_T'], independent_max_diff_log_mt=worst['ind_log_mt'])


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--tables-only', action='store_true'); a = ap.parse_args()
    warnings.filterwarnings('ignore'); np.seterr(all='ignore')
    clock = Clock()
    recs = []
    worst = {}
    parts = [part_T] if a.tables_only else [part_T, part_R1, part_R2]
    for part in parts:
        try:
            worst.update(part(recs))
        except Exception as exc:
            recs.append(dict(group=part.__name__, id='exception', passed=False, error=repr(exc)))
        print(part.__name__, 'done', f'{clock.seconds():.0f}s', flush=True)
    groups = {}
    for r in recs:
        g = groups.setdefault(r['group'], {'attempted': 0, 'passed': 0})
        g['attempted'] += 1; g['passed'] += int(r['passed'])
    extra = {'groups': groups, 'recomputation_included': not a.tables_only,
             'second_cycle_found': bool(worst.get('recomputed_hunt_new_fixed_points', 0) or worst.get('recomputed_hunt_long_runs_elsewhere', 0)),
             'counting_rule': 'T: one count per table statement; R1: one per recomputed cycle; R2: one per recomputed hunt point',
             'subsample': {'grid_rho_indices': list(GRID_RHO_IDX), 'grid_top_rho_for_mu': list(GRID_TOP_MU), 'hunt_indices': list(HUNT_IDX),
                           'hunt_time_budget_seconds': HUNT_BUDGET}}
    name = 'v3_cartography_tables_only' if a.tables_only else 'v3_cartography'
    ok = write_receipt(name, 'table', recs, worst=worst, clock=clock, extra=extra,
                       parameters={'tables': sorted(p.name for p in CART.iterdir() if p.is_file())},
                       scope='Shipped floating-point tables and a recomputed subsample. No interval enclosure; a search that finds no '
                             'second cycle does not exclude one.')
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
