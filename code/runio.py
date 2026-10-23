"""Loading of recorded runs: OUT.json or its gzip-compressed form OUT.json.gz (as published)."""
import gzip
import json
import os


def resolve(fn):
    """fn if it exists, else fn + '.gz' if that exists, else fn."""
    if os.path.exists(fn) and os.path.exists(fn + '.gz'):
        raise SystemExit('both %s and %s.gz exist; remove one (audit R7-9)' % (fn, fn))
    if os.path.exists(fn):
        return fn
    if os.path.exists(fn + '.gz'):
        return fn + '.gz'
    return fn


def exists(fn):
    return os.path.exists(resolve(fn))


def load_json(fn):
    fn = resolve(fn)
    if fn.endswith('.gz'):
        with gzip.open(fn, 'rt') as f:
            return json.load(f)
    with open(fn) as f:
        return json.load(f)
