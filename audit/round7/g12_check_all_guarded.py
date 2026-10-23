"""AUDIT-7 g12: run the project's code/check_all.py unchanged, but with every write redirected into audit/round7/final/
(check_all writes ../RESULTS.md, which is outside the auditor's area).  Python bytecode and numba caches are not
written into code/.  usage: python g12_check_all_guarded.py OUTDIR"""
import builtins, os, sys, runpy
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.abspath(sys.argv[1]) if len(sys.argv) > 1 else os.path.join(HERE, 'final')
os.makedirs(OUT, exist_ok=True)
os.environ.setdefault('NUMBA_CACHE_DIR', os.path.join(HERE, '..', 'numba_cache'))
sys.dont_write_bytecode = True
CODE = os.path.normpath(os.path.join(HERE, '..', '..', 'code'))
_open = builtins.open


def guarded(file, mode='r', *a, **k):
    if isinstance(file, str) and any(c in mode for c in 'wax+'):
        target = os.path.join(OUT, os.path.basename(file))
        print('[guard] write redirected:', file, '->', target, file=sys.stderr)
        file = target
    return _open(file, mode, *a, **k)


builtins.open = guarded
os.chdir(CODE)
sys.path.insert(0, CODE)
sys.argv = ['check_all.py']
runpy.run_path(os.path.join(CODE, 'check_all.py'), run_name='__main__')
