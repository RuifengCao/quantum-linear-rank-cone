#!/bin/sh
# R45 sampling runs (as run, with the scratch binary ./sample7q = cc/s39_sample7q.c; two queues, one per core).
# Build: gcc -O3 -march=native -o sample7q cc/s39_sample7q.c
# Then: python3 scripts/s36_one_qudit_fields.py --build s_gf --out table.npz
# (the orbits and counts are reproduced; the certificates may differ if the compiler rounds the density differently).
./sample7q 2 2000000 201 s_gf2
./sample7q 3 4000000 301 s_gf3
./sample7q 4 20000000 401 s_gf4 > s_gf4.log
( ./sample7q 5 100000000 501 s_gf5 > s_gf5.log; ./sample7q 9 100000000 901 s_gf9 > s_gf9.log; ./sample7q 11 50000000 1101 s_gf11 > s_gf11.log; ./sample7q 25 100000000 2501 s_gf25 > s_gf25.log ) &
( ./sample7q 7 100000000 701 s_gf7 > s_gf7.log; ./sample7q 13 50000000 1301 s_gf13 > s_gf13.log; ./sample7q 27 100000000 2701 s_gf27 > s_gf27.log; ./sample7q 49 50000000 4901 s_gf49 > s_gf49.log ) &
wait
./sample7q 8 100000000 801 s_gf8 > s_gf8.log
./sample7q 16 100000000 1601 s_gf16 > s_gf16.log
