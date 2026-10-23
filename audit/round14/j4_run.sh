#!/bin/sh
PY="python"
cd "$(dirname "$0")"
NUMBA_CACHE_DIR="$(cd ../numba_cache && pwd)"; export NUMBA_CACHE_DIR
R="<repo>/code/runs"
{ "$PY" -B j4_disc_leaves.py "$R/KQF.r2.json" 30 4 4 41; "$PY" -B j4_disc_leaves.py "$R/KQF.r3.json" 20 4 4 42; } > j4_chain1.out 2>&1 &
{ "$PY" -B j4_disc_leaves.py "$R/KQF2.r3.json" 30 4 4 43; "$PY" -B j4_disc_leaves.py "$R/KQF2.r4.json" 30 4 4 44; } > j4_chain2.out 2>&1 &
wait
echo "j4 done $(date +%H:%M:%S)"
