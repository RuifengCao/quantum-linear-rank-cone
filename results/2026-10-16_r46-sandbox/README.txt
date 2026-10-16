R46 sandbox runs (2026-10-16; 2-core container; no server time used).  Sections 4.13, 4.14 and 4.5 of
docs/R32_reduction_and_local_dimension.md: Lemma 4* (Lemma 4 with the derived subspace C = (A_1 + B_1) ∩ ... ∩ (A_M + B_M)),
Theorem 17 (for primes p < q the stabilizer cones are incomparable from n = 2p + 2 on), Corollary 19 (Theorem 5 from
2 max(d, d') + 2 parties), Proposition 18 (the T8 and T8' vectors lie in every cone), and DFZ arXiv:1401.2507
Theorems 3.1 and 4.1.  The statements are proofs; this folder holds the computations that check them.

s37_check_full.log     python3 scripts/s37_derived_c.py --check --full: gates G42 and G42b.  Values of (LB*) and (LA*)
                       on T(t) over eight fields (M = 4, 5, 6) and over Z_N; all 81 choices of j(.) at M = 4; the
                       balancing coefficients; the Theorem 17 witnesses (doubled codes, entropies from the [G | G]
                       code); Proposition 18; DFZ 4.1 (two transcriptions, Theorem 4.2 values, balanced minima over
                       the 8! relabellings on h_T8', h_T8, h_AG) and DFZ 3.1 (re-expanded minus printed, the
                       counterexample, T8); then all 8-tuples of subspaces of GF(2)^2, GF(3)^2, GF(5)^2, coordinate
                       lines with every line of GF(2)^4 and GF(3)^4, 1.6 M random arrangements, M = 8 (168 s).
selftest_full.log      python3.12 selftest.py --full on a clean clone of commit 9902a3d (R45) with the R46 changes
                       copied in: 57 gate lines, all pass (1036 s; nothing else running).  The gate lines of R45 are
                       unchanged.  The same with python3.13: 57 lines, all pass (1024 s; log not kept).

The exploratory probes of this round were run outside the repository: DFZ 4.1 on h_T8' in all relabellings; the face
of h_T8'; the GF(9) / GF(4) / GF(8) constructions behind Proposition 18; searches for two-copy families in a normal form
(for h_T8: exhaustive over GF(2), none; over GF(5) and GF(7), 1,760 and 13,517 normalised families, among them the
reductions of the integer family of Proposition 18; for h_T8': over GF(3), found at once); the first tests of Lemma 4*.
Everything they established that the note uses is recomputed by scripts/s37_derived_c.py.  The independent referee's
code is not part of the repository either (see the note, section 4.13, "Independent referee (R46)").
