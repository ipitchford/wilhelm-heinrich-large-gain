#!/usr/bin/env python3
"""Verify release bytes against MANIFEST.sha256 (standard library only)."""
from __future__ import annotations
import hashlib,json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]

def main():
    entries=[];failures=[];manifest=ROOT/'MANIFEST.sha256'
    if not manifest.exists():print('MANIFEST.sha256 is missing',file=sys.stderr);return 1
    for lineno,line in enumerate(manifest.read_text().splitlines(),1):
        if not line or line.startswith('#'):continue
        try:
            expected,name=line.split('  ',1);p=(ROOT/name).resolve()
            if not p.is_relative_to(ROOT) or len(expected)!=64:raise ValueError('unsafe path or invalid digest')
            got=hashlib.sha256(p.read_bytes()).hexdigest();ok=got==expected
            entries.append({'path':name,'matched':ok})
            if not ok:failures.append({'path':name,'expected':expected,'actual':got})
        except Exception as exc:failures.append({'line':lineno,'error':repr(exc)})
    passed=bool(entries) and not failures
    print(json.dumps({'number_checked':len(entries),'number_matched':sum(x['matched'] for x in entries),'passed':passed,'failures':failures,'scope':'Byte integrity only, not proof verification.'},indent=2))
    return 0 if passed else 1
if __name__=='__main__':sys.exit(main())
