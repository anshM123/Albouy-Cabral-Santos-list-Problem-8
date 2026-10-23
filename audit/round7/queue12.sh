#!/bin/sh
# AUDIT-12 core-1 queue: after snapshot part A, the final C3 / face charts (g4), then snapshot part B.
PY="${PY:-python}"
HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"
until grep -q "== snapshot A done" snap12/chainA.out 2>/dev/null; do sleep 10; done
NUMBA_CACHE_DIR="$(cd ../numba_cache && pwd)"; export NUMBA_CACHE_DIR
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
C=../../code/c3/runs
c3files=""; for c in CK CE Cm PhiP PhiM L2aE L2am L2aYp L2aYm L2KE L2Km L2KYp L2KYm FE FM FIP FIM; do c3files="$c3files $C/$c.json"; done
{
  date "+== %H:%M:%S g4: final C3 charts + face charts: replay, 20 spot leaves x 2 points, 100 leaves re-certified per file"
  "$PY" -B g4_c3_spot.py 20 2 23 100 $c3files
  date "+== g4 final done %H:%M:%S"
} > final12/chainD.out 2>&1
sh snap12b.sh
echo "queue12 finished $(date +%H:%M:%S)"
