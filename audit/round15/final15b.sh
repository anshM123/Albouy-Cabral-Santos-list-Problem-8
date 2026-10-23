#!/bin/sh
# AUDIT-15, chains B-E, G and the g4 step of H again, with relative paths (the absolute paths contain spaces).
PY="python"
HERE="$(cd "$(dirname "$0")" && pwd)"; cd "$HERE"
OUT="$HERE/out"
NUMBA_CACHE_DIR="$(cd ../numba_cache && pwd)"; export NUMBA_CACHE_DIR
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
CR=../../code/runs; CC=../../code/c3/runs
NEWKQ="KQ0T.r4 KQ0T.r5 KQ2.r5 KQ2.r6 KQ2.r7 KQ2.r8 KQB.r7 KQB.r8 KQB.r9 KQB.r10 KQF.r2 KQF.r3 KQF2.r3 KQF2.r4 KQK1.r7 KQK1.r8 KQK1.r9 KQK1.r10 KQK1.r11 KQK1.r12 KQS.r2 KQS.r3 KQS.r4 KQS.r5"
kqrel=""; for n in $NEWKQ; do kqrel="$kqrel $CR/$n.json"; done
( { date "+== B start %H:%M:%S"; "$PY" -B j6_recert_all_tw_disc.py $CR/KQ0T.r4.json $CR/KQ0T.r5.json $CR/KQ2.r5.json $CR/KQ2.r6.json $CR/KQ2.r7.json $CR/KQ2.r8.json $CR/KQF.r2.json $CR/KQF.r3.json; date "+== B done %H:%M:%S"; } > "$OUT/B_tw_disc.out" 2>&1 ) &
( { date "+== C start %H:%M:%S"; "$PY" -B j6_recert_all_tw_disc.py $CR/KQB.r7.json $CR/KQB.r8.json $CR/KQB.r9.json $CR/KQB.r10.json $CR/KQS.r2.json $CR/KQS.r3.json $CR/KQS.r4.json $CR/KQS.r5.json; date "+== C done %H:%M:%S"; } > "$OUT/C_tw_disc.out" 2>&1 ) &
( { date "+== D start %H:%M:%S"; "$PY" -B j6_recert_all_tw_disc.py $CR/KQK1.r7.json $CR/KQK1.r8.json $CR/KQK1.r9.json $CR/KQK1.r10.json $CR/KQK1.r11.json $CR/KQK1.r12.json $CR/KQF2.r3.json $CR/KQF2.r4.json; date "+== D done %H:%M:%S"; } > "$OUT/D_tw_disc.out" 2>&1 ) &
( { date "+== E start %H:%M:%S"; "$PY" -B j6_recert_all_tw_disc.py $CC/Ca.r8.json $CC/Ca.r9.json $CC/Ca.r10.json; date "+== E done %H:%M:%S";
    cd ../round7 && "$PY" -B g4_c3_spot.py 30 3 47 150 ../../code/c3/runs/Ca.r8.json ../../code/c3/runs/Ca.r9.json ../../code/c3/runs/Ca.r10.json; date "+== E g4 done %H:%M:%S"; } > "$OUT/E_tw_c3.out" 2>&1 ) &
( cd ../round7 && { date "+== G start %H:%M:%S"; "$PY" -B g1_kq_wind.py 400 100 15 2 47 $(for n in $NEWKQ; do printf ' ../../code/runs/%s.json' $n; done); date "+== G done %H:%M:%S"; } > "$OUT/G_g1.out" 2>&1 ) &
wait
echo "final15b done $(date +%H:%M:%S)"
