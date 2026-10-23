#!/bin/sh
# AUDIT-12 snapshot pass, part A (one core): the frozen 15:00 snapshot only (levels up to the closed ones).
PY="${PY:-python}"
HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"
OUT="$HERE/snap12"; mkdir -p "$OUT"
NUMBA_CACHE_DIR="$(cd ../numba_cache && pwd)"; export NUMBA_CACHE_DIR
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
export AUDIT_MAXLEV="KQ0T:3;KQ2:4;KQB:6;KQK1:6;KQF:1;KQF2:2;KQS:1;C3 Ca:7"
NAMES="KQ0T KQ0T.r1 KQ0T.r2 KQ0T.r3 KQ2 KQ2.r1 KQ2.r2 KQ2.r3 KQ2.r4 KQB KQB.r1 KQB.r2 KQB.r3 KQB.r4 KQB.r5 KQB.r6 KQK1 KQK1.r1 KQK1.r2 KQK1.r3 KQK1.r4 KQK1.r5 KQK1.r6 KQF KQF.r1 KQF2 KQF2.r1 KQF2.r2 KQS KQS.r1"
{
  date "+== snapshot A start %H:%M:%S"
  echo "== g2: replay of the snapshot families up to the closed levels (+ face charts for the p-dependency)"
  "$PY" -B g2_replay_all.py --prune-samples 4 --only "KQ0T,KQ2,KQB,KQK1,KQF,KQF2,KQS,C3 Ca,face FE,face FM,face FIP,face FIM"
  date "+== %H:%M:%S g9: 200 regular/deep leaves per file"
  "$PY" -B g9_recert_all.py 200 29 $NAMES --kinds reg,deep,tw,disc
  date "+== %H:%M:%S g9: 300 winding leaves per file (nb_cert, then nb_cert2)"
  "$PY" -B g9_recert_all.py 300 29 $NAMES --kinds wind
  date "+== snapshot A done %H:%M:%S"
} > "$OUT/chainA.out" 2>&1
