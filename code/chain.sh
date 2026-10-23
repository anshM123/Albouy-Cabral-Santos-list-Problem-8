#!/bin/sh
# Continuation chain of one run family: wait for level FROM to finish (its JSON is written when every box has
# returned), list its open boxes (next_level.py), run the next level on them with a per-box time cap, and repeat
# until no box is open.
# usage: sh chain.sh KIND NAME FROM WORKERS TMAX S0 S1 S2 S3        (run from code/; KIND kq or c3)
PY="${PY:-python}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1   # one BLAS thread per process (memory)
EXTRA="${EXTRA:-}"                 # extra kq_run flags, e.g. EXTRA=--facesplit
KIND=$1; NAME=$2; LEV=$3; W=$4; TMAX=$5; S0=$6; S1=$7; S2=$8; S3=$9
if [ "$KIND" = kq ]; then DIR=runs; else DIR=c3/runs; fi
while true; do
  if [ "$LEV" = 0 ]; then F=$DIR/$NAME.json; else F=$DIR/$NAME.r$LEV.json; fi
  while [ ! -f "$F" ]; do sleep 20; done
  # the driver writes the file in one json.dump at the end: wait until its size is stable
  SZ0=-1; SZ1=$(stat -c %s "$F"); while [ "$SZ0" != "$SZ1" ]; do sleep 10; SZ0=$SZ1; SZ1=$(stat -c %s "$F"); done
  OUT=$("$PY" next_level.py "$KIND" "$NAME" "$LEV") || { echo "$(date +%H:%M) next_level failed, retrying"; sleep 30; continue; }
  echo "$(date +%H:%M) $OUT"
  N=$(echo "$OUT" | sed -n 's/.*open boxes \([0-9][0-9]*\) ->.*/\1/p')   # the summary line only (not "dropped ...")
  if [ "$N" = 0 ]; then echo "$(date +%H:%M) $NAME: no open box after level $LEV"; break; fi
  NEXT=$((LEV + 1))
  if [ "$KIND" = kq ]; then
    "$PY" kq_run.py --boxes "$DIR/${NAME}_r${NEXT}_boxes.json" --scale $S0 $S1 $S2 $S3 "$TMAX" "$W" "$DIR/$NAME.r$NEXT.json" $EXTRA > "$DIR/$NAME.r$NEXT.log" 2>&1
  else
    (cd c3 && "$PY" c3_run.py "$NAME" --boxes "runs/${NAME}_r${NEXT}_boxes.json" --scale $S0 $S1 $S2 $S3 "$TMAX" "$W" "runs/$NAME.r$NEXT.json" $EXTRA > "runs/$NAME.r$NEXT.log" 2>&1)
  fi
  LEV=$NEXT
done
