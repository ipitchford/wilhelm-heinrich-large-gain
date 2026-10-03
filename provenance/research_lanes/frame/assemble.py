#!/usr/bin/env python3
"""Merge bibliography fragments, then report label/citation problems for the assembled manuscript."""
import re, glob, sys, collections
from pathlib import Path
P = Path('/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03/v3/wilhelm_heinrich_v3_0/paper')
man = (P/'manuscript.tex').read_text()
order = re.findall(r'\\input\{sections/([^}]+)\}', man)
secs = {s: (P/'sections'/f'{s}.tex').read_text() if (P/'sections'/f'{s}.tex').exists() else None for s in order}
missing = [s for s, t in secs.items() if t is None]
# --- bibliography
entries = {}
dups = []
for f in sorted(glob.glob(str(P/'bib'/'*.tex'))):
    for m in re.finditer(r'\\bibitem\{([^}]+)\}(.*?)(?=\\bibitem|\Z)', Path(f).read_text(), re.S):
        key, body = m.group(1), ' '.join(l for l in m.group(2).splitlines() if not l.strip().startswith('%')).strip()
        if key in entries:
            if ' '.join(entries[key].split()) != ' '.join(body.split()):
                dups.append((key, Path(f).name))
            continue
        entries[key] = body
cited = collections.OrderedDict()
for s in order:
    if secs[s] is None: continue
    for m in re.finditer(r'\\cite(?:\[[^\]]*\])?\{([^}]+)\}', secs[s]):
        for k in m.group(1).split(','):
            cited.setdefault(k.strip(), s)
for m in re.finditer(r'\\cite(?:\[[^\]]*\])?\{([^}]+)\}', man):
    for k in m.group(1).split(','): cited.setdefault(k.strip(), 'manuscript')
undefined = [k for k in cited if k not in entries]
unused = [k for k in entries if k not in cited]
def sortkey(k):
    b = re.sub(r'\\[a-zA-Z]+|[{}\\\'"`^~]', '', entries[k])
    return b.lower()
keys = sorted([k for k in entries if k in cited], key=sortkey)
out = ['\\begin{thebibliography}{99}'] + [f'\\bibitem{{{k}}} {entries[k]}' for k in keys] + ['\\end{thebibliography}', '']
(P/'bibliography.tex').write_text('\n'.join(out))
# --- labels
labels = collections.Counter(); where = collections.defaultdict(list)
for s in order:
    if secs[s] is None: continue
    for l in re.findall(r'\\label\{([^}]+)\}', secs[s]): labels[l] += 1; where[l].append(s)
refs = collections.defaultdict(set)
for s in order:
    if secs[s] is None: continue
    for l in re.findall(r'\\(?:eq)?ref\{([^}]+)\}', secs[s]): refs[l].add(s)
for l in re.findall(r'\\(?:eq)?ref\{([^}]+)\}', man): refs[l].add('manuscript')
print('sections missing:', missing)
print('bib entries used:', len(keys), '| cited but undefined:', undefined, '| defined but unused (dropped):', unused, '| conflicting duplicates:', dups)
print('duplicate labels:', {l: where[l] for l, c in labels.items() if c > 1})
print('undefined refs:', {l: sorted(v) for l, v in refs.items() if l not in labels})
bad = collections.defaultdict(list)
for s in order:
    if secs[s] is None: continue
    for pat in (r'\blanes?\b', r'\bthe brief\b', r'\bagents?\b', r'v2\.0', r'\byou\b', r'\byour\b', r'it is worth noting', r'[Cc]rucially', r'[Ii]mportantly', r'[Nn]otably'):
        n = len(re.findall(pat, secs[s]))
        if n: bad[s].append(f'{pat}:{n}')
print('style flags:', dict(bad))
