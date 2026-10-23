#!/bin/sh
# Resume an interrupted continuation level (records kept by --resume), then continue as chain.sh from that level.
# usage: sh resume_chain.sh KIND NAME LEV WORKERS TMAX S0 S1 S2 S3      (EXTRA: the level's own driver flags)
PY="${PY:-python}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
EXTRA="${EXTRA:-}"
KIND=$1; NAME=$2; LEV=$3; W=$4; TMAX=$5; S0=$6; S1=$7; S2=$8; S3=$9
if [ "$KIND" = kq ]; then
  "$PY" kq_run.py --boxes "runs/${NAME}_r${LEV}_boxes.json" --scale $S0 $S1 $S2 $S3 "$TMAX" "$W" "runs/$NAME.r$LEV.json" $EXTRA --resume >> "runs/$NAME.r$LEV.log" 2>&1
else
  (cd c3 && "$PY" c3_run.py "$NAME" --boxes "runs/${NAME}_r${LEV}_boxes.json" --scale $S0 $S1 $S2 $S3 "$TMAX" "$W" "runs/$NAME.r$LEV.json" $EXTRA --resume >> "runs/$NAME.r$LEV.log" 2>&1)
fi
exec sh chain.sh "$KIND" "$NAME" "$LEV" "$W" "$TMAX" "$S0" "$S1" "$S2" "$S3"
