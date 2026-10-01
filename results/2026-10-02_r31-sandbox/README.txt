R31 sandbox runs (2-core sandbox, 2026-10-01/02).  Row numbers refer to
results/2026-10-02_a2-status/undecided_reps.npy unless stated otherwise.

two_ci_{a,b,c,d}.json/.log  two simultaneous common informations (prototype of scripts/s25_two_ci.py):
                            every unordered pair of disjoint pairs with |X|, |Y| <= 2, Shannon on the
                            8-element extension.  a, b: the 10 smallest undecided orbits (rows 273 396 434
                            484 517 807 970 973 981 995); c, d: 6 larger ones (rows 52 193 362 475 699 904).
                            All feasible, about 90,000 LPs, no solver trouble (floating-point verdicts).
s25_981.json/.log           the repository script on row 981 (R27 frontier ray 294): same result as the
                            prototype (106 pairs, 5,565 pairs of pairs, feasible).
gf4_small10.json            GF(4) graph states at lambda = 1 (s24 kernel, 290 s each) on the 10 smallest
                            undecided orbits: none realised; best residual per row.
gf4_sample16.json           the same, 55 s each, on 16 random larger undecided orbits: none; residual 16-74.
exhaustive4.log             cc/s24_exhaustive6.c on rows 43 44 195 225 of data/stab5_extreme_reps.npy
                            (S_P = 1 for every party): no graph state with one qubit (all 2^15) or one qutrit
                            (all 3^15) per party realises them at lambda = 1.
exhaustive_controls.log     the same program on rows 3 and 45 (realised at lambda = 1): hits found over
                            GF(2) and GF(3), as it must.
