R45 sandbox runs (2026-10-15; 2-core container; no server time used).  Section 4.10 of
docs/R32_reduction_and_local_dimension.md, "Other fields": the entropy vectors of GF(q)-linear stabilizer states of
seven systems with one qudit each (the cut-rank vectors of weighted graphs over GF(q) on seven vertices), over 13
fields; Propositions 15 and 16.  Propositions 15 and 16 are proofs; this folder holds the computations and checks.

sampling_runs.sh         the sampler runs behind data/s36_one_qudit_fields.npz (cc/s39_sample7q.c, run under its
                         scratch name ./sample7q, identical apart from the usage string): field, number of graphs,
                         seed.  The order in the file is not the order in which they ran.
sampling.log             per field: graphs sampled, seed, distinct vectors sampled, and orbits / vectors after closing
                         under S_7 (equal to the table); the sampler's own output lines.
table_build.log          python3 scripts/s36_one_qudit_fields.py --build: the table from the sampler outputs; its arrays
                         are those of data/s36_one_qudit_fields.npz.
enumeration.log          cc/s40_enum7q.c: GF(2) (with and without the relabelling reduction), GF(3) (both), GF(4) and
                         GF(5), three slices of GF(3), and slice 0 of 400 over GF(7) (the heaviest); counts, times and
                         the comparison of the closed outputs with the stored lists and the table.
enumeration_counts.py    the numbers of prefixes and graphs cc/s40_enum7q.c examines, counted from support patterns
enumeration_counts.log   alone (q = 2, 3, 4, 5, 7); the GF(7) count (815,738,564,501 graphs) gives the cost estimate.
s36_check_full.log       python3 scripts/s36_one_qudit_fields.py --check --full: the checks of gates G41 and G41b
                         (table, exhaustive GF(2)-GF(5), target searches for Proposition 15, a sampling run,
                         Proposition 16).
prop16.log               python3 scripts/s36_one_qudit_fields.py --prop16: the 70 minors of [I_4 | J_4 - I_4], the
                         brute-force representability counts of T8, T8' and AG(3,2) over GF(2) .. GF(9), and S_T8
                         against S_T8'.
selftest_full.log        python3.12 selftest.py --full on a clean clone of commit 5ff0736 (R44) with the R45 changes
                         copied in: 55 gate lines, all pass (921 s; nothing else running).  The gate lines of R44
                         are unchanged.  README.md and the R32 note were edited afterwards (wording only).

Not kept: the raw sampler outputs s_gf<q>.bin / .w (about 30 MB; sampling_runs.sh regenerates them) and the raw
outputs of cc/s40_enum7q.c (G41b regenerates them).
