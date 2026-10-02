R33 sandbox runs (2026-10-02, 2-core container; no server time used).

s28_t2.log           scripts/s28_generalized_cubes.py --t 2 --fields 2 3 5 --exhaustive01
                     t = 2 (the cube states of R32): self-duality weights, h_S = 2r, contraction = configuration
                     J - I; Lemma 4 (LB) = -2 on the GF(2) state, (LA) = -2 on the GF(3), GF(5) states;
                     the GF(2) and GF(q) vectors differ in 1 coordinate. Lemma 4 on all 2,097,152 {0,1}-point
                     configurations of GF(2)^3, GF(3)^3, GF(5)^3: the valid form never < 0, the
                     wrong-characteristic form reaches -1 (sensitivity). Pena-Sarria Example 6 (a)/(b) agree.
s28_t3.log           scripts/s28_generalized_cubes.py --t 3 --fields 3 5 7 11 --random 300 --exhaustive01
                     t = 3 (Theorem 2' for p = 3 against q = 5, 7, 11): (LB) = -2 on the GF(3) state, (LA) = -2
                     on the GF(q) states; the GF(3) and GF(q) vectors differ in 1 coordinate (of 511).
                     Valid forms on 300 random arrangements and 2.2 M random {0,1}-point configurations of
                     GF(q)^4 per field: never < 0. (Random sampling of GF(q)^4 rarely hits the J - I
                     configuration, so the sensitivity lines are not informative here except over GF(2).)
s28_t5.log           scripts/s28_generalized_cubes.py --t 5 --fields 3 5 7 11
                     t = 5 (Theorem 2' for p = 5 against q = 7, 11): (LB) = -2 on the GF(5) state, (LA) = -2 on
                     the GF(7), GF(11) states, 1 differing coordinate (of 8,191). GF(3) is included as an
                     extra check (t = 5 is self-dual over GF(3) too, (LA) = -2 there); the pair (3, 5) itself
                     is covered by t = 3.
dfz1401_misprint.log scripts/s28_generalized_cubes.py --dfz1401
                     DFZ arXiv:1401.2507v1 Theorem 3.1 as printed: -4 on C = X = a line, all else zero;
                     125 / 254 violating arrangements among all subspace arrangements of GF(2)^2 / GF(3)^2.
selftest_full.log    python3 selftest.py --full: 36 gate lines, ALL PASS (275 s).
