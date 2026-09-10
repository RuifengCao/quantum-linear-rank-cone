#!/bin/sh
# usage: ./s4b_run_lrs.sh clr5.ine [estimate|run]
# PRIMARY s4 engine since R5.  Reverse search: memory-safe (~MB), single-threaded.
# R5 measured on the real cone: estimate ~1.4M bases / ~1 h; actual run confirmed.
# 'estimate' probes the tree first (seconds..minutes); 'run' enumerates for real.
IN=${1:?input .ine file}
MODE=${2:-run}
command -v lrs >/dev/null 2>&1 || { echo "lrs not found: apt-get install -y lrslib  (or conda -c conda-forge lrslib)" >&2; exit 1; }
if [ "$MODE" = estimate ]; then
  cp "$IN" "${IN%.ine}.est.ine"
  printf 'maxdepth 2\nestimates 10\n' >> "${IN%.ine}.est.ine"
  lrs "${IN%.ine}.est.ine" | grep -E "Estimat|Totals" 
  exit 0
fi
setsid nohup lrs "$IN" "${IN%.ine}.lrs.out" > "${IN%.ine}.lrs.log" 2>&1 &
echo "lrs pid $!   output ${IN%.ine}.lrs.out"
echo "progress:  wc -l ${IN%.ine}.lrs.out    (rays stream out; CLR_5 expects ~7,943 + headers)"
echo "kill:      pkill -x lrs"
