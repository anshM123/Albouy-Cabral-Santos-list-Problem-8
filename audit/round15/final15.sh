#!/bin/sh
# AUDIT-15 final pass over the complete run set: 8 serial chains (8 cores), outputs in audit/round15/out/.
PY="python"
HERE="$(cd "$(dirname "$0")" && pwd)"; cd "$HERE"
OUT="$HERE/out"; mkdir -p "$OUT"
NUMBA_CACHE_DIR="$(cd ../numba_cache && pwd)"; export NUMBA_CACHE_DIR
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
R7="$HERE/../round7"; R13="$HERE/../round13"
CR="$HERE/../../code/runs"; CC="$HERE/../../code/c3/runs"
NEWKQ="KQ0T.r4 KQ0T.r5 KQ2.r5 KQ2.r6 KQ2.r7 KQ2.r8 KQB.r7 KQB.r8 KQB.r9 KQB.r10 KQF.r2 KQF.r3 KQF2.r3 KQF2.r4 KQK1.r7 KQK1.r8 KQK1.r9 KQK1.r10 KQK1.r11 KQK1.r12 KQS.r2 KQS.r3 KQS.r4 KQS.r5"
kqpaths=""; for n in $NEWKQ; do kqpaths="$kqpaths $CR/$n.json"; done
capaths="$CC/Ca.r8.json $CC/Ca.r9.json $CC/Ca.r10.json"
# A: the project checker (unchanged, output redirected) and my own replay of every family
( cd "$R7" && { date "+== A start %H:%M:%S"; "$PY" -B g12_check_all_guarded.py "$OUT"; date "+== %H:%M:%S g2 full replay";
  "$PY" -B g2_replay_all.py --prune-samples 4; date "+== A done %H:%M:%S"; } > "$OUT/A_checkall_g2.out" 2>&1 ) &
# B, C, D, E: EVERY tau-winding and disc-winding leaf of the new files re-certified
( { date "+== B start %H:%M:%S"; "$PY" -B j6_recert_all_tw_disc.py $CR/KQ0T.r4.json $CR/KQ0T.r5.json $CR/KQ2.r5.json $CR/KQ2.r6.json $CR/KQ2.r7.json $CR/KQ2.r8.json $CR/KQF.r2.json $CR/KQF.r3.json; date "+== B done %H:%M:%S"; } > "$OUT/B_tw_disc.out" 2>&1 ) &
( { date "+== C start %H:%M:%S"; "$PY" -B j6_recert_all_tw_disc.py $CR/KQB.r7.json $CR/KQB.r8.json $CR/KQB.r9.json $CR/KQB.r10.json $CR/KQS.r2.json $CR/KQS.r3.json $CR/KQS.r4.json $CR/KQS.r5.json; date "+== C done %H:%M:%S"; } > "$OUT/C_tw_disc.out" 2>&1 ) &
( { date "+== D start %H:%M:%S"; "$PY" -B j6_recert_all_tw_disc.py $CR/KQK1.r7.json $CR/KQK1.r8.json $CR/KQK1.r9.json $CR/KQK1.r10.json $CR/KQK1.r11.json $CR/KQK1.r12.json $CR/KQF2.r3.json $CR/KQF2.r4.json; date "+== D done %H:%M:%S"; } > "$OUT/D_tw_disc.out" 2>&1 ) &
( { date "+== E start %H:%M:%S"; "$PY" -B j6_recert_all_tw_disc.py $capaths; date "+== E done %H:%M:%S"; } > "$OUT/E_tw_c3.out" 2>&1 ) &
# F: sampled re-certification of regular, deep and winding leaves of every new KQ file (g9, final engines)
( cd "$R7" && { date "+== F start %H:%M:%S"; "$PY" -B g9_recert_all.py 400 47 $NEWKQ --kinds reg,deep,wind; date "+== F done %H:%M:%S"; } > "$OUT/F_g9.out" 2>&1 ) &
# G: g1 (replay, 400 winding + 100 count + 100 tw + 100 disc leaves re-certified, 100-digit spots incl. tw/disc)
( cd "$R7" && { date "+== G start %H:%M:%S"; "$PY" -B g1_kq_wind.py 400 100 15 2 47 $kqpaths; date "+== G done %H:%M:%S"; } > "$OUT/G_g1.out" 2>&1 ) &
# H: 100-digit spot checks of the new files (g3), C3 charts r8-r10 (g4), own zeros inside TW contours + enclosure (j1)
( { date "+== H start %H:%M:%S"; cd "$R7" && "$PY" -B g3_spot_all.py 60 3 47 $NEWKQ;
  "$PY" -B g4_c3_spot.py 30 3 47 150 $capaths; cd "$R13";
  for n in KQ2.r5 KQB.r8 KQK1.r7 KQS.r3 KQ0T.r4 KQF.r2; do "$PY" -B j1_tw_leaves.py "$CR/$n.json" kq @ 0 0 0 0 15 3 3 51; done;
  "$PY" -B j1_tw_leaves.py "$CC/Ca.r8.json" Ca @ 0 0 0 0 15 3 3 52;
  date "+== H done %H:%M:%S"; } > "$OUT/H_spots.out" 2>&1 ) &
wait
grep -h -E "VIOLATION|MISMATCH|PROBLEM|FAILURE|NOT OK|NOT COMPLETE|BOUNDS: NO|not present|problems [1-9]|violations [1-9]|mismatches [1-9]|failures [1-9]|odd N_mult [1-9]|Traceback|Error" "$OUT"/*.out > "$OUT/ALERTS.txt"
echo "final15 done $(date +%H:%M:%S); alerts: $(wc -l < "$OUT/ALERTS.txt")"
