#!/bin/sh
# AUDIT-12 snapshot pass, part C (the other core, after the final-family g1): g4 on Ca up to r7.
PY="${PY:-python}"
HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"
OUT="$HERE/snap12"; mkdir -p "$OUT"
NUMBA_CACHE_DIR="$(cd ../numba_cache && pwd)"; export NUMBA_CACHE_DIR
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
C=../../code/c3/runs
cafiles="$C/Ca.json $C/Ca.r1.json $C/Ca.r2.json $C/Ca.r3.json $C/Ca.r4.json $C/Ca.r5.json $C/Ca.r6.json $C/Ca.r7.json"
{
  date "+== snapshot C start %H:%M:%S: g4 on Ca + r1..r7: replay, 20 spot leaves x 2 points, 40 leaves re-certified per file"
  "$PY" -B g4_c3_spot.py 20 2 29 40 $cafiles
  date "+== snapshot C done %H:%M:%S"
} > "$OUT/chainC.out" 2>&1
