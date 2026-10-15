#!/bin/bash
# adversarial searches for Theorem 12, run one after the other (sandbox, 2 cores shared)
cd "$(dirname "$0")/../.."   # the repository root
run() { echo "### $*"; nice -n 10 python3 scripts/s34_six_party_survey.py "$@"; }
run --transfer-anneal "1,1,1,1,1,1,1" --seconds 60 --seed 4201
run --transfer-anneal "2,2,2,2,2,2,2" --seconds 120 --seed 4202
run --transfer-anneal "2,2,2,2,2,2,2" --seconds 120 --seed 4203
run --transfer-anneal "3,3,3,3,3,3,3" --seconds 240 --seed 4204
run --transfer-anneal "2,2,2,1,1,1,2" --seconds 120 --seed 4205
run --transfer-anneal "1,2,1,2,1,2,2" --seconds 120 --seed 4206
run --transfer-anneal "3,2,3,2,3,2,3" --seconds 120 --seed 4207
run --transfer-anneal "2,2,2,2,2,2,2" --seconds 120 --seed 4208 --pull-nf 1.0
run --transfer-anneal "3,3,3,3,3,3,3" --seconds 240 --seed 4209 --pull-nf 1.0
echo "### done"
