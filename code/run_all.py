#!/usr/bin/env python3
"""Re-run every computation the paper relies on as evidence.  Ordinary numerics; never a formal proof.

    python code/run_all.py               framework suite on the 31 saved cycles, then all v3 checks (about ten minutes)
    python code/run_all.py --regenerate  also recompute the 31 amplitude-held boundary-value problems of the framework suite
    python code/run_all.py --quick       skip only the cartography recomputation (the shipped tables are still checked)

No network access is used.  Every step runs as one single-threaded process (BLAS and numba threads set to 1), one step
at a time.  Each step writes its own strict-JSON receipt under evidence/ and its output under logs/.  All steps are run
even if an earlier one fails, so that one failure does not hide another; the exit status is non-zero if any required
step, or any required receipt, failed.  The deliberately wrong assertions of the controls must fail: a control receipt
that passes is itself a failure.

Run on a copy if the SHA-256 manifest of the release must stay unchanged: receipts, logs and generated tables are rewritten.
"""
from __future__ import annotations
import argparse, hashlib, json, os, platform, resource, subprocess, sys, time
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CODE = ROOT / 'code'
EV = ROOT / 'evidence'
LOGS = ROOT / 'logs'

FRAMEWORK_RECEIPTS = ['exact_checks', 'numerical_audit', 'cross_method', 'reciprocal_rates', 'envelope_checks', 'negative_controls']
V3_RECEIPTS = ['v3_exact_identities', 'v3_gauges', 'v3_largegain_constants', 'v3_largegain_lemmas', 'v3_attractor',
               'v3_weak_adjoint', 'v3_wedge', 'v3_cartography', 'v3_negative', 'v3_controls']
MUST_FAIL = ['deliberately_wrong_exact_receipt.json', 'v3/v3_control_exact_wrong_sign.json',
             'v3/v3_control_largegain_wrong_c1.json', 'v3/v3_control_largegain_wrong_cp.json']


def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def checked_json(path):
    def bad(value):
        raise ValueError(f'Non-standard JSON value: {value}')
    return json.loads(path.read_text(), parse_constant=bad)


def peak_rss_mb():
    r = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    return r / (1024.0 * 1024.0) if sys.platform == 'darwin' else r / 1024.0


def versions():
    out = {}
    for name in ('numpy', 'scipy', 'sympy', 'mpmath', 'numba', 'llvmlite'):
        try:
            out[name] = __import__(name).__version__
        except Exception as exc:        # recorded, not hidden; numba is needed only by the cartography recomputation
            out[name] = 'not importable: ' + repr(exc)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--regenerate', action='store_true', help='recompute the 31 saved cycles of the framework suite')
    ap.add_argument('--quick', action='store_true', help='skip the cartography recomputation (tables are still checked)')
    args = ap.parse_args()
    LOGS.mkdir(exist_ok=True)
    (LOGS / 'v3').mkdir(exist_ok=True)
    code_files = sorted(CODE.glob('*.py')) + sorted((CODE / 'v3').glob('*.py')) + sorted((CODE / 'v3' / 'lib').glob('*.py'))
    receipt = {'schema_version': '3.0', 'started_utc': datetime.now(timezone.utc).isoformat(), 'python': sys.version,
               'platform': platform.platform(), 'packages': versions(), 'regenerated_cycles': args.regenerate, 'quick': args.quick,
               'threads': 1, 'network_used': False,
               'code_sha256': {str(p.relative_to(ROOT)): digest(p) for p in code_files}, 'steps': [], 'passed': False}
    steps = [('exact_checks', 'exact_checks.py', [])]
    if args.regenerate:
        steps.append(('cycle_regeneration', 'cycles.py', []))
    elif not (EV / 'cycle_trials.json').exists():
        raise FileNotFoundError('Missing cycle trials: use --regenerate.')
    steps += [('numerical_audit', 'numerical_audit.py', []), ('additional_checks', 'additional_checks.py', [])]
    for name in V3_RECEIPTS:
        extra = ['--tables-only'] if (name == 'v3_cartography' and args.quick) else []
        steps.append((name, f'v3/{name}.py', extra))
    env = os.environ.copy()
    env.update(OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', MKL_NUM_THREADS='1', VECLIB_MAXIMUM_THREADS='1',
               NUMEXPR_NUM_THREADS='1', NUMBA_NUM_THREADS='1', PYTHONWARNINGS='ignore')

    def run_step(name, script, extra):
        t = time.monotonic()
        log = LOGS / f'{name}.log'
        with log.open('w') as f:
            result = subprocess.run([sys.executable, str(CODE / script)] + extra, cwd=ROOT, env=env, stdout=f, stderr=subprocess.STDOUT)
        row = {'name': name, 'script': 'code/' + script, 'arguments': extra, 'exit_status': result.returncode,
               'elapsed_seconds': time.monotonic() - t, 'log': str(log.relative_to(ROOT)), 'peak_rss_mb_so_far': peak_rss_mb()}
        receipt['steps'].append(row)
        print(f"{name:28s} exit {result.returncode}   {row['elapsed_seconds']:7.1f} s", flush=True)
        return result.returncode

    problems = []
    for name, script, extra in steps:
        if run_step(name, script, extra):
            problems.append(f'step {name} exited with a non-zero status; see logs/{name}.log')
    receipt['computation_seconds'] = sum(s['elapsed_seconds'] for s in receipt['steps'])
    receipt['peak_rss_mb'] = peak_rss_mb()

    # scientific verdicts are read from the receipts, not from exit codes alone
    summary = {}
    try:
        trial = checked_json(EV / 'cycle_trials.json')
        expected = {f'mu_{m:.8g}_h_{h:.8g}' for m in [.1, .3, 1., 10 / 3, 10.] for h in [.05, .2, .5, 1., 2., 3.]} | {'mu_0.3_h_0.01'}
        if {r['id'] for r in trial['attempts']} != expected or len(trial['attempts']) != 31 or not all(r['accepted'] for r in trial['attempts']):
            problems.append('the 31 accepted target cycles of the framework suite were not reproduced')
        for name in FRAMEWORK_RECEIPTS:
            data = checked_json(EV / f'{name}.json')
            n_att = data.get('number_checked', data.get('number_attempted'))
            summary[name] = {'kind': data.get('kind'), 'passed': bool(data.get('passed')), 'attempted': n_att, 'n_passed': data.get('number_passed')}
            if not data.get('passed', False):
                problems.append(f'receipt {name} did not pass')
        for name in V3_RECEIPTS:
            fname = 'v3_cartography_tables_only' if (name == 'v3_cartography' and args.quick) else name
            data = checked_json(EV / 'v3' / f'{fname}.json')
            summary[fname] = {'kind': data['kind'], 'passed': bool(data['passed']), 'attempted': data['counts']['attempted'],
                              'n_passed': data['counts']['passed']}
            if not data['passed']:
                problems.append(f'receipt {fname} did not pass')
        for rel in MUST_FAIL:
            if checked_json(EV / rel)['passed']:
                problems.append(f'control receipt {rel} passed a deliberately wrong assertion')
    except Exception as exc:
        problems.append('receipts could not be read: ' + repr(exc))
    receipt['receipts'] = summary
    receipt['passed'] = not problems
    receipt['problems'] = problems
    receipt['completed_utc'] = datetime.now(timezone.utc).isoformat()
    receipt['scope'] = ('All specified finite checks were executed. Exact symbolic equalities are exact; all other checks are '
                        'floating-point or high-precision computations without interval enclosure. This is not a proof assistant.')
    (EV / 'run_receipt.json').write_text(json.dumps(receipt, indent=2, allow_nan=False) + '\n')

    # generated tables, then bookkeeping consistency (both read the run receipt written above)
    for name, script in (('make_tables', 'make_tables.py'), ('check_consistency', 'check_consistency.py')):
        if run_step(name, script, []):
            problems.append(f'step {name} exited with a non-zero status; see logs/{name}.log')
    receipt['passed'] = not problems
    receipt['problems'] = problems
    receipt['completed_utc'] = datetime.now(timezone.utc).isoformat()
    receipt['total_seconds'] = sum(s['elapsed_seconds'] for s in receipt['steps'])
    (EV / 'run_receipt.json').write_text(json.dumps(receipt, indent=2, allow_nan=False) + '\n')
    print()
    print('computation: %.1f min; peak memory of a step: %.0f MB' % (receipt['computation_seconds'] / 60.0, receipt['peak_rss_mb']))
    if problems:
        print('FAILED:')
        for p in problems:
            print('  -', p)
    else:
        print('ALL REQUIRED CHECKS PASSED (finite computations; see the scope field of each receipt)')
    return 0 if receipt['passed'] else 1


if __name__ == '__main__':
    sys.exit(main())
