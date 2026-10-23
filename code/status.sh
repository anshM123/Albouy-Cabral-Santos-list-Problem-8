#!/bin/sh
date +%H:%M:%S
for f in runs/*.log c3/runs/*.log; do
  case "$f" in *smallchain*|*chain*) continue;; esac
  echo "$(basename $f .log): $(tail -n 1 $f 2>/dev/null | cut -c1-120)"
done
for f in runs/*.fails.jsonl; do s=$(wc -c < $f); [ "$s" -gt 0 ] && echo "FAILS/INCOMPLETE logged: $f ($s bytes)"; done
