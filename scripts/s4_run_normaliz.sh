#!/bin/sh
# usage: ./s4_run_normaliz.sh clr5.in [threads] [primal|dual]
# R5 NOTE: on the real CLR_5 the default primal path EXPLODED (91.5M intermediate
# hyperplanes at gen 62/1905, ~1.31x/gen).  Use scripts/s4b_run_lrs.sh as the
# primary engine; this script remains for cross-checks (try 'dual').
# Extreme-ray enumeration is the memory-hungry step: run it ALONE on the node.
# Binary resolution order: $NORMALIZ  ->  normaliz on PATH  ->  bundled kit/bin/normaliz.
IN=${1:?input .in file}
T=${2:-$(nproc)}
ALGO=${3:-primal}
AFLAG=""
[ "$ALGO" = dual ] && AFLAG="-d"
HERE=$(cd "$(dirname "$0")" && pwd)
NMZ=${NORMALIZ:-}
[ -z "$NMZ" ] && command -v normaliz >/dev/null 2>&1 && NMZ=normaliz
[ -z "$NMZ" ] && [ -x "$HERE/../bin/normaliz" ] && NMZ="$HERE/../bin/normaliz"
if [ -z "$NMZ" ]; then
  echo "normaliz not found. Options:" >&2
  echo "  * use the bundled static binary: chmod +x $HERE/../bin/normaliz" >&2
  echo "  * or: conda install -c conda-forge normaliz" >&2
  exit 1
fi
echo "using: $NMZ ($($NMZ --version 2>/dev/null | head -1))"
setsid nohup "$NMZ" $AFLAG -c -x=$T "$IN" > "${IN%.in}.log" 2>&1 &
echo "normaliz pid $!  threads $T  log ${IN%.in}.log"
echo "watch:  tail -f ${IN%.in}.log      kill:  pkill -x normaliz"
