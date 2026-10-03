#!/usr/bin/env python3
"""Fill README.md from the template, write STATUS.json, copy frame scripts into provenance."""
import json, re, shutil, subprocess
from pathlib import Path
ROOT = Path('/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03')
B = ROOT/'v3'/'wilhelm_heinrich_v3_0'
F = ROOT/'work'/'frame'

log = (B/'paper'/'manuscript.log').read_text(errors='replace')
pages = int(re.search(r'Output written on manuscript\.pdf \((\d+) pages', log).group(1))
aux = (B/'paper'/'manuscript.aux').read_text()
lab = {m.group(1): m.group(2) for m in re.finditer(r'\\newlabel\{([^}]+)\}\{\{([^}]*)\}', aux)}

# checks table from the generated counts
counts = (B/'paper'/'generated'/'v3_counts.tex').read_text()
rows = []
for line in counts.splitlines():
    if '&' in line and not line.lstrip().startswith('Check'):
        cells = [c.strip().rstrip('\\').strip() for c in line.split('&')]
        if len(cells) >= 4:
            name = re.sub(r'--', '–', cells[0]); passed = cells[2].replace('\\,', ','); att = cells[3].replace('\\,', ',')
            rows.append(f'| {name} | {cells[1]} | {passed} / {att} |')
checks = '| Check | Kind | Passed / attempted |\n|---|---|---|\n' + '\n'.join(rows)
rr = json.loads((B/'evidence'/'run_receipt.json').read_text())
minutes = None
try:
    from datetime import datetime
    a = datetime.fromisoformat(rr['started_utc']); b = datetime.fromisoformat(rr['completed_utc'])
    minutes = round((b - a).total_seconds()/60, 1)
except Exception:
    pass
checks += (f"\n\nThe reference run of `code/run_all.py` {'passed' if rr.get('passed') else 'FAILED'}"
           + (f" in {minutes} minutes" if minutes else '') + '; its receipt is `evidence/run_receipt.json`.'
           " A supplementary run of the large-gain lemma script at mu = 0.02, k = 1000, outside the default run, is kept with one comparison"
           " recorded as failed: two numbers of size about 3.26e9 that differ by about 1e-3, below the resolution of the quadrature;"
           " the independent integrator covers the same parameter pair and passes (`evidence/v3/v3_largegain_lemmas_small_mu_supplement.json`).")

reproduce = """Python 3.13 with the packages in `requirements.txt` (exact versions in `ENVIRONMENT.json`). Work in a copy of the bundle if the manifest must stay valid: a run rewrites receipts, logs and generated tables.

```sh
python -m venv .venv && . .venv/bin/activate
python -m pip install -r requirements.txt
python code/run_all.py              # every check, on the 31 saved cycles
python code/run_all.py --regenerate # also recompute the 31 saved cycles
python code/run_all.py --quick      # skip the cartography recomputation
python code/check_consistency.py    # receipts, tables, macros and references agree
python code/verify_manifest.py      # byte integrity against MANIFEST.sha256
```

No network access or API key is needed. Building the PDF is described in `BUILD.md`."""

t = (F/'README_template.md').read_text()
t = (t.replace('{{PAGES}}', str(pages)).replace('{{REPRODUCE}}', reproduce).replace('{{CHECKS}}', checks)
      .replace('{{EVIDENCE_SECTION}}', lab['sec:evidence']).replace('{{STATUS_SECTION}}', lab['sec:status']))
assert '{{' not in t
(B/'README.md').write_text(t)

status = {
 'version': '3.0', 'date': '2026-10-03', 'status': 'research_candidate',
 'title': 'Convex attractor geometry, exact stability formulas and large-gain uniqueness for the Wilhelm–Heinrich oscillator',
 'original_all_parameter_conjecture_solved': False,
 'uniqueness_proved_regions': ['k >= 1000 for 0 < mu <= 1', 'k >= 1000 mu^3 for mu >= 1', 'non-explicit right-hand neighbourhood of each Hopf point'],
 'open_region': 'k_H < k < 1000 max(1, mu^3) outside the Hopf collars',
 'certified_multicycle_counterexample': False,
 'neutral_cycle_excluded_in_open_region': False,
 'conditional_results': {'thm:wedge': 'assumption (F), a consequence expected from Fenichel\'s theorem; not verified against the original source'},
 'self_contained_for_stated_claims': True,
 'independent_specialist_review': False,
 'internal_adversarial_review': True,
 'cross_model_review_of_all_new_sections': True,
 'priority_established': False,
 'interval_certified_computations': False,
 'default_driver_passed': bool(rr.get('passed')),
 'default_driver_minutes': minutes,
 'pdf_pages': pages,
 'claims_file': 'CLAIMS.json',
}
(B/'STATUS.json').write_text(json.dumps(status, indent=2, ensure_ascii=False) + '\n')

dst = B/'provenance'/'research_lanes'/'frame'
dst.mkdir(parents=True, exist_ok=True)
for f in ['assemble.py', 'make_claims.py', 'finalise.py', 'README_template.md']:
    shutil.copy2(F/f, dst/f)
shutil.copy2(ROOT/'v3'/'PLAN.md', dst/'PLAN.md')
for d in ['sol', 'framework_checks']:
    src = ROOT/'work'/d; tgt = B/'provenance'/'research_lanes'/d
    if src.exists():
        shutil.copytree(src, tgt, dirs_exist_ok=True, ignore=shutil.ignore_patterns('__pycache__', '*.aux', '*.log'))
print('README, STATUS written; pages', pages, '; runtime', minutes)
