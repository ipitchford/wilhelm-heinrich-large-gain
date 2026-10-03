#!/usr/bin/env python3
"""Bookkeeping consistency of receipts, generated tables, macros, references and provenance.

This verifies that the numbers printed from the evidence are exactly what the receipts contain, that every receipt and
file the paper refers to exists, that the receipts were produced by the shipped code, and that failures which are meant
to be kept are still there.  It does not verify any proof and it is not a numerical enclosure.
Run after code/run_all.py (which calls it as its last step); writes evidence/release_checks.json.
"""
from __future__ import annotations
import hashlib, json, re, sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
EV = ROOT / 'evidence'
EV3 = EV / 'v3'
GEN = ROOT / 'paper' / 'generated'
SECTIONS = ROOT / 'paper' / 'sections'
LANES = ROOT / 'provenance' / 'research_lanes'
V3_RECEIPTS = ['v3_exact_identities', 'v3_gauges', 'v3_largegain_constants', 'v3_largegain_lemmas', 'v3_attractor',
               'v3_weak_adjoint', 'v3_wedge', 'v3_cartography', 'v3_negative', 'v3_controls']
KIND_EXACT = 'exact_symbolic_equality'
FRAMEWORK_MACROS = ['ExponentAgreement', 'OrthogonalityError', 'ReturnClosure', 'MomentError', 'CrossMatrixError',
                    'CrossExponentError', 'ReciprocalError', 'EnvelopeHessianError', 'EnvelopeDefectError']
ZERO = ('0', 'Matrix([[0, 0, 0], [0, 0, 0], [0, 0, 0]])', 'Matrix([[0], [0], [0]])', 'Matrix([[0, 0, 0]])', 'Matrix([[0], [0]])')


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def read(p):
    def bad(x):
        raise ValueError('Non-standard JSON constant ' + x)
    return json.loads(Path(p).read_text(), parse_constant=bad)


def run():
    checks = []

    def ck(name, ok, detail=None):
        checks.append({'id': name, 'passed': bool(ok), 'detail': detail})

    def guarded(name, fn):
        try:
            fn()
        except Exception as exc:
            ck(name + '_exception', False, repr(exc))

    # ------------------------------------------------------------------ framework suite (as in v2.0)
    def framework():
        exact = read(EV / 'exact_checks.json')
        ck('exact_count', len(exact['checks']) == exact['number_checked'] == 50)
        ck('exact_actual_equalities', all(row['passed'] and row['residual'] in ZERO for row in exact['checks']),
           'Residual strings are checked as well as pass flags.')
        for name, key, n in [('numerical_audit', 'records', 31), ('cross_method', 'records', 4), ('reciprocal_rates', 'records', 12),
                             ('envelope_checks', 'records', 4), ('negative_controls', 'controls', 4)]:
            d = read(EV / f'{name}.json')
            ck(name + '_actual_count', len(d[key]) == n and d.get('number_checked', d.get('number_attempted')) == n)
            ck(name + '_all_predicates', d['passed'] and all(x['passed'] for x in d[key]) and d['number_passed'] == n)
        trials = read(EV / 'cycle_trials.json'); rows = trials['attempts']
        ck('BVP_count', len(rows) == 31 and len({q['id'] for q in rows}) == 31 and trials['number_attempted'] == 31)
        ck('BVP_acceptance', trials['number_accepted'] == 31 and all(q['solver_success'] and q['accepted'] and all(q['acceptance_checks'].values())
                                                                     and not q['guard_triggered'] for q in rows))
        bad = []
        for q in rows:
            with np.load(ROOT / q['mesh'], allow_pickle=False) as a:
                if not all(np.all(np.isfinite(a[k])) for k in a.files):
                    bad.append(q['id'])
        ck('saved_cycle_meshes_finite', not bad, bad)
        bad = []
        for q in read(EV / 'numerical_audit.json')['records']:
            with np.load(ROOT / q['profile'], allow_pickle=False) as a:
                if not (all(np.all(np.isfinite(a[k])) for k in a.files) and all(q['checks'].values())):
                    bad.append(q['id'])
        ck('saved_profiles_finite_and_all_predicates', not bad, bad)
        summary = read(EV / 'manuscript_summary.json')
        ck('framework_table_sources_unchanged', all(sha(EV / n) == h for n, h in summary['source_sha256'].items()))
        old = read(EV / 'exact_checks_initial_failure.json')
        ck('framework_real_checker_failure_retained', not old['passed'] and old['number_passed'] == 49 and old['number_checked'] == 50)
        wrong = read(EV / 'deliberately_wrong_exact_receipt.json')
        ck('framework_deliberate_wrong_claim_rejected', not wrong['passed'] and wrong['number_passed'] < wrong['number_checked'])
    guarded('framework', framework)

    # ------------------------------------------------------------------ v3 receipts
    R = {}

    def receipts():
        for n in V3_RECEIPTS:
            d = read(EV3 / f'{n}.json')
            R[n] = d
            fields = all(k in d for k in ('name', 'kind', 'passed', 'counts', 'worst_residuals', 'parameters', 'seconds', 'scope', 'records'))
            ck(n + '_fields', fields and d['name'] == n)
            ck(n + '_certification_flag', (d['kind'] == KIND_EXACT) != bool(d.get('not_interval_certified', False)),
               'every receipt that is not an exact symbolic equality carries not_interval_certified: true')
            ck(n + '_passed', d['passed'] is True and d['counts']['failed'] == 0 and d['counts']['passed'] == d['counts']['attempted'] > 0)
            recs = d['records']
            if all('attempted' in r for r in recs):
                tot, npass = sum(r['attempted'] for r in recs), sum(r['n_passed'] for r in recs)
            else:
                tot, npass = len(recs), sum(bool(r['passed']) for r in recs)
            ck(n + '_counts_match_records', tot == d['counts']['attempted'] and npass == d['counts']['passed'], [tot, npass])
        ex = R['v3_exact_identities']
        ck('v3_exact_residuals_are_zero', ex['kind'] == KIND_EXACT and all(r['passed'] and r['residual'] == '0' for r in ex['records']),
           'residual strings are checked as well as pass flags')
        ck('v3_gauges_all_signs_resolved', R['v3_gauges']['worst_residuals']['unresolved_sign_samples'] == 0)
        ca = R['v3_cartography']
        ck('v3_cartography_includes_recomputation', ca.get('recomputation_included') is True and not ca.get('second_cycle_found'))
    guarded('v3_receipts', receipts)

    # ------------------------------------------------------------------ failures that must be there
    def kept():
        w = read(EV3 / 'v3_control_exact_wrong_sign.json')
        failed = [r['id'] for r in w['records'] if not r['passed']]
        ck('control_exact_wrong_sign_fails', not w['passed'] and failed == ['injected_wrong_sign_star_gauge'])
        for tag in ('c1', 'cp'):
            w = read(EV3 / f'v3_control_largegain_wrong_{tag}.json')
            ck(f'control_largegain_wrong_{tag}_fails', not w['passed'] and any(not r['passed'] for r in w['records']))
        g0 = read(EV3 / 'v3_gauges_initial_failure.json')
        bad = [r for r in g0['records'] if not r['passed']]
        ck('retained_gauge_checker_failure', not g0['passed'] and g0['counts'] == {'attempted': 66, 'passed': 65, 'failed': 1}
           and len(bad) == 1 and bad[0]['id'] == 'ext_mu_0.02_rho_20' and [k for k, v in bad[0]['checks'].items() if not v] == ['cycle_located'])
        sm = read(EV3 / 'v3_largegain_lemmas_small_mu_supplement.json')
        fc = sm['records'][0]['failed_comparisons']
        ck('retained_small_mu_supplement', not sm['passed'] and sm['counts'] == {'attempted': 37, 'passed': 36, 'failed': 1}
           and [f['name'] for f in fc] == ['C: M+ < k T-'] and sm.get('part_of_default_run') is False)
    guarded('kept_failures', kept)

    # ------------------------------------------------------------------ generated tables and macros
    macros = {}

    def generated():
        sys.path.insert(0, str(ROOT / 'code'))
        import make_tables as mt
        s2, o2 = mt.run_v2()
        s3, o3 = mt.build_v3()
        for name, text in {**o2, **o3}.items():
            p = GEN / name
            ck('generated_' + name, p.exists() and p.read_text() == text, 'file on disk equals the text rebuilt from the receipts')
        stored = read(EV3 / 'v3_summary.json')
        ck('v3_summary_sources_unchanged', all(sha(ROOT / n) == h for n, h in stored['source_sha256'].items()))
        ck('v3_summary_macros_current', stored['macros'] == s3['macros'])
        for line in (GEN / 'v3_macros.tex').read_text().splitlines() + (GEN / 'evidence_macros.tex').read_text().splitlines():
            m = re.match(r'\\newcommand\{\\([A-Za-z]+)\}\{(.*)\}$', line)
            if m:
                macros[m.group(1)] = m.group(2)
        ck('macro_names_letters_only', all(n.isalpha() for n in macros) and all(n in macros for n in FRAMEWORK_MACROS))
        ck('every_macro_matches_its_receipt', all(macros.get(n) == v['tex'] for n, v in s3['macros'].items()),
           f"{len(s3['macros'])} v3 macros compared with values recomputed from the receipts")
    guarded('generated', generated)

    # ------------------------------------------------------------------ what the paper refers to
    def references():
        used, inputs, paths, files, receipts_named = set(), set(), {}, {}, set()
        tex = sorted(SECTIONS.glob('*.tex')) + [p for p in (ROOT / 'paper').glob('*.tex') if not p.name.startswith('test_')]
        for f in tex:
            s = f.read_text()
            used |= set(re.findall(r'\\(Vthree[A-Za-z]+)', s))
            used |= {m for m in FRAMEWORK_MACROS if re.search(r'\\' + m + r'(?![A-Za-z])', s)}
            inputs |= set(re.findall(r'\\input\{(generated/[^}]+)\}', s))
            for m in re.findall(r'\\path\{([^}]*)\}', s):
                paths.setdefault(m, f.name)
            for m in re.findall(r'\\texttt\{([^}]*)\}', s):
                t = m.replace('\\_', '_').replace('\\allowbreak', '').strip()
                if re.search(r'\.(py|out|json|jsonl|csv|txt|md|jl|sh)$', t) and ' ' not in t:
                    files.setdefault(t, f.name)
        ck('every_used_macro_is_defined', used <= set(macros), sorted(used - set(macros)))
        ck('every_generated_input_exists', all((ROOT / 'paper' / (i if i.endswith('.tex') else i + '.tex')).exists() for i in inputs), sorted(inputs))
        names = {}
        for base in (LANES, ROOT / 'code', ROOT / 'evidence', ROOT / 'provenance' / 'v2_0_bundle'):
            for p in base.rglob('*'):
                if p.is_file() and '__pycache__' not in p.parts:
                    names.setdefault(p.name, []).append(p)
        roots = (LANES, ROOT / 'provenance', ROOT, EV, EV3, ROOT / 'code', ROOT / 'code' / 'v3', ROOT / 'code' / 'v3' / 'lib')
        missing = []
        for ref, where in list(paths.items()) + list(files.items()):
            if ref.endswith('/'):
                ok = any((r / ref).is_dir() for r in roots)
            else:
                ok = any((r / ref).is_file() for r in roots) or ('/' not in ref and ref in names)
            if not ok:
                missing.append([ref, where])
        ck('every_referenced_file_exists', not missing, missing or f'{len(paths) + len(files)} file references resolved inside the bundle')
        named = [ref for ref in list(paths) + list(files) if ref.startswith('evidence/')]
        missing_named = [ref for ref in named if not ((ROOT / ref).is_dir() if ref.endswith('/') else (ROOT / ref).is_file())]
        ck('every_referenced_receipt_exists', not missing_named, missing_named or f'{len(named)} evidence references resolved')
    guarded('references', references)

    # ------------------------------------------------------------------ run record and code identity
    def runrecord():
        rr = read(EV / 'run_receipt.json')
        ck('run_record_passed', rr['passed'] and all(x['exit_status'] == 0 for x in rr['steps']),
           {'regenerated_cycles': rr['regenerated_cycles'], 'quick': rr['quick']})
        ck('run_record_single_thread_no_network', rr.get('threads') == 1 and rr.get('network_used') is False)
        changed = [n for n, h in rr['code_sha256'].items() if not (ROOT / n).exists() or sha(ROOT / n) != h]
        ck('executed_code_unchanged', not changed, changed)
        current = sorted(str(p.relative_to(ROOT)) for pat in ('code/*.py', 'code/v3/*.py', 'code/v3/lib/*.py') for p in ROOT.glob(pat))
        ck('no_unrecorded_code', current == sorted(rr['code_sha256']), sorted(set(current) ^ set(rr['code_sha256'])))
        ck('peak_memory_below_1500_MB', rr['peak_rss_mb'] < 1500, rr['peak_rss_mb'])
        envp = ROOT / 'ENVIRONMENT.json'
        if envp.exists():
            e = read(envp)
            ck('environment_matches_run', all(e.get(k) == v for k, v in rr['packages'].items()), rr['packages'])
        if not rr['regenerated_cycles']:
            src = ROOT / 'provenance' / 'v2_0_bundle' / 'evidence' / 'cycles'
            diff = [p.name for p in sorted((EV / 'cycles').glob('*.npz')) if not (src / p.name).exists() or sha(p) != sha(src / p.name)]
            ck('saved_cycles_identical_to_v2_0', not diff and len(list((EV / 'cycles').glob('*.npz'))) == 31, diff)
    guarded('run_record', runrecord)

    # ------------------------------------------------------------------ provenance
    def provenance():
        man = ROOT / 'provenance' / 'v2_0_bundle' / 'MANIFEST.sha256'
        bad, n = [], 0
        for line in man.read_text().splitlines():
            if not line or line.startswith('#'):
                continue
            h, name = line.split('  ', 1)
            n += 1
            p = man.parent / name
            if not p.is_file() or sha(p) != h:
                bad.append(name)
        ck('v2_0_bundle_matches_its_manifest', n > 0 and not bad, bad or f'{n} files')
        bad = []
        for p in sorted((ROOT / 'code' / 'v3' / 'lib').glob('*.py')):
            head = p.read_text().splitlines()[:3]
            m = re.search(r'research file `([^`]+)`', head[1]) if len(head) > 1 else None
            if not (m and head[2].startswith('# --- Changes from the original:')):
                bad.append(p.name + ': no provenance header')
                continue
            rel = m.group(1)
            rel = rel[len('lanes/'):] if rel.startswith('lanes/') else rel[len('work/'):] if rel.startswith('work/') else rel
            if not (LANES / rel).is_file():
                bad.append(p.name + ': original not in provenance: ' + rel)
        ck('library_copies_have_provenance', not bad, bad)
        cart = EV3 / 'cartography'
        diff = [p.name for p in sorted(cart.iterdir()) if p.is_file() and p.name != 'README.md' and p.suffix != '.npz'
                and not ((LANES / 'cartography' / p.name).is_file() and sha(p) == sha(LANES / 'cartography' / p.name))]
        ck('shipped_cartography_tables_equal_research_outputs', not diff, diff or 'every shipped table and log is byte-identical to the research output')
        ck('skipped_list_present', (LANES / 'SKIPPED.md').is_file())
        ck('no_font_files', not any(p.suffix.lower() in ('.ttf', '.otf', '.woff', '.woff2') for p in ROOT.rglob('*') if p.is_file()))
    guarded('provenance', provenance)

    result = {'kind': 'release_bookkeeping_and_finite_receipt_consistency', 'checks': checks, 'number_checked': len(checks),
              'number_passed': sum(x['passed'] for x in checks), 'passed': all(x['passed'] for x in checks),
              'scope': 'Not proof verification and not a numerical enclosure.'}
    (EV / 'release_checks.json').write_text(json.dumps(result, indent=2, allow_nan=False, default=str) + '\n')
    print(f"{result['number_passed']}/{result['number_checked']} consistency checks passed")
    for q in checks:
        if not q['passed']:
            print(q)
    return result['passed']


if __name__ == '__main__':
    sys.exit(0 if run() else 1)
