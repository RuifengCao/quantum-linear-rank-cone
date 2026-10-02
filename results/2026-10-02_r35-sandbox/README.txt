R35 sandbox runs (2026-10-02, 2-core container; no server time used).

s30_profile.log          scripts/s30_dimension_profile.py --n4 --one-cut 2 3 4 5 6 7 --six
                         Theorem 7: 46/46 extreme rays of h^-1(Shannon + Ingleton) on five elements have integer
                         graph states with unimodular cut blocks (stored lifts; re-verified over GF(p), p <= 13).
                         Proposition 6: rank changes of the generalised cube points over GF(l) vs Q exactly as
                         predicted for t = 2..7, l <= 13, l not dividing t - 1 (30 cases; single cut when l = t).
                         Proposition 8: the six stored certificates den * (f + f o D) = sum num_i e_i are exact.
s30_probe6_n5pilot.log   scripts/s30_dimension_profile.py --probe6 --n5-pilot 60
                         Encoding check of Proposition 8: f(r + r*) >= 0 on all {0,1}-point configurations of
                         GF(p)^3 and 2.2 M random ones of GF(p)^4, while f itself reaches -1 there.
                         Five-party pilot (not certified): of the first 60 known extreme rays of Stab_5, 21 lift to
                         all fields with all signs +1, 19 after a sign search, 20 not lifted.
selftest_full.log        python3 selftest.py --full: 40 gate lines, ALL PASS (369 s).

An independent referee (separate agent, own code) re-derived Theorem 7 (own double description and lrs: the same 46
rays), Proposition 6 (30 cases) and Proposition 8 (certificates re-found; 1,750 pure states x 5,040 relabellings x
6 forms, minimum 0); its scripts are not part of the repository.
