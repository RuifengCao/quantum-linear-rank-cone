R36 sandbox runs (2026-10-02, 2-core container; no server time used).

s31_six_party.log        scripts/s31_six_party.py --derive --witness --invariant --face --lemmas 200 --stress 400
                         --exhaustive
                         Theorem 9: the transfer form f9 (67 terms) recomputed from Lemmas 0, A-E and Lemma 4 (LB)_2
                         equals data/s31_six_party.npz; the seven-qubit simplex-code state has h_S = r_F + r_F^*
                         (state vector) and f9(h_S) = -3; (f9 + f9 o D)(h_F) = -6; invariant form
                         (1/3) f9 = 210 v_P + 30 alpha + 21 nu - 421 mu (DFZ (65) variant: 154, 26, 5, -309);
                         exact Shannon faces at h_F and h_NF are cone(r, r^*) (ranks 125 / 126); lemma-by-lemma checks
                         with the actual subspaces on 200 families over GF(2), GF(3), GF(5) (all lemma slacks >= 0;
                         chain >= 0 in odd characteristic, -3 over GF(2)); f9 >= 0 on 400 random arrangements over
                         GF(3), GF(5), GF(7) and on all {0,1}-point configurations of GF(3)^3 and GF(5)^3, while it
                         reaches -3 over GF(2)^3.
s30_six_bkl.log          scripts/s30_dimension_profile.py --six --bkl
                         Proposition 8 with eight forms (BKL (6.8) and (6.17) corrected added): all certificates
                         exact.  BKL transcription checks: Fano / non-Fano values, the counterexample to the printed
                         (6.17) (-3 over GF(2), GF(3), GF(5)), exhaustive minima 0, 0 and -3 (printed (6.17)).
selftest_full.log        python3 selftest.py --full: 44 gate lines, ALL PASS (357 s).

An independent referee (separate agent, own code) re-derived Lemmas 0, A-E and the composition, recomputed f9 exactly
(no difference), confirmed the values from r_F + r_F^* and from the state vector, re-derived the odd-order group case
(Lemma 4(G), Theorem 1(d)), and stress-tested f9 >= 0: about 8.2 M arrangements each over GF(3), GF(5), GF(7), all
105,413,504 seven-tuples of points of PG(2,3) or zero (minimum 0), about 13,000 subgroup families of odd-order abelian
groups; over GF(2), f9 < 0 exactly on the 168 Fano frames.  A second agent explored the reverse direction at n = 6
(open; exact face reduction recorded in section 4.8 of the note).  Their scripts are not part of the repository.
