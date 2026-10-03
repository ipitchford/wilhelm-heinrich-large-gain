#!/usr/bin/env python3
"""Shared helpers for the v3 checks: paths, strict JSON receipts, timing.

Every receipt written through `write_receipt` is strict JSON (no NaN or Infinity) and carries
name, kind, passed, counts, worst_residuals, parameters, seconds and, for anything that is not an
exact symbolic equality, "not_interval_certified": true.
"""
from __future__ import annotations
import json, math, sys, time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
V3 = Path(__file__).resolve().parent
LIB = V3 / 'lib'
EVID = ROOT / 'evidence' / 'v3'
if str(LIB) not in sys.path:
    sys.path.insert(0, str(LIB))
if str(V3) not in sys.path:
    sys.path.insert(0, str(V3))

KINDS = {
    'exact': 'exact_symbolic_equality',
    'float': 'floating_point_diagnostic',
    'mp': 'high_precision_constant_check',
    'control': 'deliberate_error_regression',
    'table': 'floating_point_table_recomputation',
}


def clean(o):
    """Convert numpy / mpmath objects to plain JSON types; non-finite floats become None."""
    try:
        import numpy as np
    except Exception:  # pragma: no cover
        np = None
    if isinstance(o, dict):
        return {str(k): clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple, set)):
        return [clean(v) for v in o]
    if np is not None:
        if isinstance(o, np.ndarray):
            return clean(o.tolist())
        if isinstance(o, np.bool_):
            return bool(o)
        if isinstance(o, np.integer):
            return int(o)
        if isinstance(o, np.floating):
            o = float(o)
    if isinstance(o, bool) or o is None or isinstance(o, (int, str)):
        return o
    if isinstance(o, float):
        return o if math.isfinite(o) else None
    if isinstance(o, complex):
        return [clean(o.real), clean(o.imag)]
    return str(o)


def strict_write(path: Path, data: dict):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(clean(data), indent=1, allow_nan=False) + '\n')


def read_json(path: Path):
    def bad(x):
        raise ValueError('non-standard JSON constant ' + x)
    return json.loads(Path(path).read_text(), parse_constant=bad)


class Clock:
    def __init__(self):
        self.t0 = time.monotonic()

    def seconds(self):
        return time.monotonic() - self.t0


def write_receipt(name, kind, records, *, worst, parameters, clock, scope, extra=None, passed=None,
                  attempted=None, n_passed=None, out=None):
    """records: list of dicts each with a boolean 'passed'.  Returns the receipt's 'passed' flag."""
    if attempted is None:
        attempted = len(records)
    if n_passed is None:
        n_passed = sum(bool(r.get('passed')) for r in records)
    if passed is None:
        passed = attempted > 0 and n_passed == attempted
    receipt = {
        'schema_version': '3.0',
        'name': name,
        'kind': KINDS[kind],
        'passed': bool(passed),
        'counts': {'attempted': int(attempted), 'passed': int(n_passed), 'failed': int(attempted - n_passed)},
        'worst_residuals': worst,
        'parameters': parameters,
        'seconds': clock.seconds(),
        'scope': scope,
    }
    if kind != 'exact':
        receipt['not_interval_certified'] = True
    if extra:
        receipt.update(extra)
    receipt['records'] = records
    strict_write(out or (EVID / f'{name}.json'), receipt)
    print(f"{name}: {n_passed}/{attempted} passed -> {'PASS' if passed else 'FAIL'}  [{receipt['seconds']:.1f}s]", flush=True)
    for r in records:
        if not r.get('passed'):
            print('   FAILED:', json.dumps(clean({k: v for k, v in r.items() if k not in ('detail_rows',)}))[:600], flush=True)
    return bool(passed)


LOGS = ROOT / 'logs' / 'v3'


def exec_script(path, *, argv=None, init=None, log_name=None):
    """Execute a wrapped research script in its own namespace, with stdout captured.

    Returns (namespace, captured_text, exit_status).  A SystemExit raised by the script is caught and its
    status recorded; any other exception propagates.  The captured text is written to logs/v3/<log_name>.
    Used for the lane scripts under lib/ that run at import time and report through print().
    """
    import contextlib, io
    path = Path(path)
    ns = {'__name__': '__wrapped__', '__file__': str(path)}
    if init:
        ns.update(init)
    buf = io.StringIO()
    old_argv = sys.argv
    sys.argv = [str(path)] + list(argv or [])
    status = 0
    try:
        with contextlib.redirect_stdout(buf):
            try:
                exec(compile(path.read_text(), str(path), 'exec'), ns)
            except SystemExit as exc:
                status = exc.code if isinstance(exc.code, int) else (0 if exc.code is None else 1)
    finally:
        sys.argv = old_argv
    text = buf.getvalue()
    if log_name:
        LOGS.mkdir(parents=True, exist_ok=True)
        (LOGS / log_name).write_text(text)
    return ns, text, status


def capture(fn, *args, log_name=None, append=False, **kwargs):
    """Call fn with stdout captured; returns (result, text).  Text is written/appended to logs/v3/<log_name>."""
    import contextlib, io
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        result = fn(*args, **kwargs)
    text = buf.getvalue()
    if log_name:
        LOGS.mkdir(parents=True, exist_ok=True)
        with (LOGS / log_name).open('a' if append else 'w') as f:
            f.write(text)
    return result, text
