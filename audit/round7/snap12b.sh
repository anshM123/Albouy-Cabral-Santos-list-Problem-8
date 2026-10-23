#!/bin/sh
# AUDIT-12 snapshot pass, part B (one core): g1 on the KQ/KQF snapshot files, g4 on Ca up to r7.
PY="${PY:-python}"
HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"
OUT="$HERE/snap12"; mkdir -p "$OUT"
NUMBA_CACHE_DIR="$(cd ../numba_cache && pwd)"; export NUMBA_CACHE_DIR
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
R=../../code/runs; C=../../code/c3/runs
NAMES="KQ0T KQ0T.r1 KQ0T.r2 KQ0T.r3 KQ2 KQ2.r1 KQ2.r2 KQ2.r3 KQ2.r4 KQB KQB.r1 KQB.r2 KQB.r3 KQB.r4 KQB.r5 KQB.r6 KQK1 KQK1.r1 KQK1.r2 KQK1.r3 KQK1.r4 KQK1.r5 KQK1.r6 KQF KQF.r1 KQF2 KQF2.r1 KQF2.r2 KQS KQS.r1"
files=""; for n in $NAMES; do files="$files $R/$n.json"; done
cafiles="$C/Ca.json $C/Ca.r1.json $C/Ca.r2.json $C/Ca.r3.json $C/Ca.r4.json $C/Ca.r5.json $C/Ca.r6.json $C/Ca.r7.json"
{
  date "+== snapshot B start %H:%M:%S"
  echo "== g1: replay + 10 count leaves per file, 5 spot leaves x 2 points at 100 digits (KQ and KQF engines)"
  "$PY" -B g1_kq_wind.py 0 10 5 2 29 $files
  echo "== (g4 on Ca: snap12c.sh, run on the other core)"
  date "+== snapshot B done %H:%M:%S"
} > "$OUT/chainB.out" 2>&1
