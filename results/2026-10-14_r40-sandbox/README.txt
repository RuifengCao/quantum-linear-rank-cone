R40 sandbox runs (2026-10-14 09:05 to about 13:00; 2-core container; no server time used).
Two roadmap items: S5(a), finished (docs/R32_reduction_and_local_dimension.md section 4.9), and M3, steps 1-2 (six
parties, the other direction; section 4.10).

S5(a): the five rows that R39 left open (587, 589, 612, 630, 635) are covered.  All 641 known extreme rays of Stab_5
are now certified for every local dimension d: 566 by unimodular realisations (data/s32_stab5_unimodular.npz,
unchanged) and 75 by covering certificates (data/s32_stab5_cover.npz: 14 new pieces, 9 rational realisations and 5
GF(2) realisations).

s5a_campaign.log        the steps after the R39 cutoff, in order (the same scratch orchestration around
                        cc/s32_unimodular.c as in R39; the log continues results/2026-10-14_r39-sandbox/campaign.log):
                          push 2   rational searches, 1200 s per row, on the five rows: from the qubit certificate
                                   (start temperature 1.0) and from R39's best near miss (start temperature 0.3)
                          push 3   step A (as in R39: rational and right mod the uncovered primes, from the row's
                                   rational realisation, 60 s), then rounds of 600-s rational searches from the best
                                   near miss so far (start temperature 0.3), each followed by step A
                          push 4   row 612: rounds of two 600-s searches from its two best near misses
                          push 5   row 612: rounds of two 300-s searches from its rational realisation, one right mod
                                   all uncovered primes (3, 5, 73), one right mod one of them
                        Every reported matrix had its bad primes recomputed exactly before it entered the state.
s5a_raw_kernel_output.txt   the raw kernel output of every call after the R39 cutoff ("#" lines: step, seconds,
                        seed, active primes, ring, alphabet, start and end temperature, sweeps, neutral-move
                        probability, extra checks, dump flag, number of jobs, wall time; FOUND lines carry the
                        matrices; BEST lines the best residual and, with "dump 1", the best matrix).  The near-miss
                        starts are BEST matrices of earlier calls (in this file or in the R39 file).
s5a_summary.json        the five rows: pieces kept in the data (minimal covering subsets), bad primes, provenance
s5a_check_cover.log     python3 scripts/s32_unimodular.py --check-cover (199 pieces, 75 rows, ALL OK)

M3, steps 1-2 (section 4.10): states of seven systems with at most one qudit per system.

m3_enumeration.log      cc/s35_enum7.c for p = 2 and p = 3 (the lists in data/s34_seven_qudit_vectors.npz)
m3_survey.log           python3 scripts/s34_six_party_survey.py --survey: the 53 orbits of qutrit-only vectors and the
                        qubit searches at multiplicity 2 (60 s per orbit; raw kernel lines at the end)
m3_nonfano.log          the non-Fano state at multiplicities 1-3 (runs of 120 s and 600 s) and the Steane controls
                        (cc/s34_anneal7.c, driven by a scratch script; run 3 with s34_six_party_survey.run_anneal)
m3_fields.log           python3 scripts/s34_six_party_survey.py --fields 2,3,4,5,7,8,9 --controls 2
                        (cc/s36_target7.c: graphs over GF(q) with cut-rank vector S_F or S_NF; Proposition 10), and
                        the counts of cc/s35_enum7.c for the same vectors
m3_check_full.log       python3 scripts/s34_six_party_survey.py --check --full
m3_selftest.log         python3 scripts/s34_six_party_survey.py --self-test

Both: selftest_full.log  python3 selftest.py --full (49 gate lines)
