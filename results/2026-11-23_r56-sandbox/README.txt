Gate 1 of the second paper (R56), sandbox: undecided orbits excluded by two common informations, and more realised by
a fixed-structure CSS search.  Inputs: results/2026-11-23_r54-sandbox/undecided_reps.npy (the 793 orbits undecided
after R54; 'the R54 list'), results/2026-11-23_ci2-ingleton/ (job ci2-ingleton, run on a 32-thread laptop on
2026-11-23, largest N first; its 6 h deadline cut it off after 433 of the 793 rows: all 401 with N >= 36 and 32 of
the 128 with N = 34; the job ci2-ingleton2 resumes it), the sandbox searches of 2026-11-23 and data/.
Nothing in data/ changed and gates G0-G49 are as before; gate G50 checks this folder and the job folder.  Counts
(status.json): 1,915,913 excluded orbits of the R30 list (1,915,910 + 3); 1,003 + 5 = 1,008 known orbits of extreme
rays of Stab5; 793 - 5 - 3 = 785 undecided.

certificates_ci2.json  one exact certificate for each of the 3 rows that ci2-ingleton reported infeasible (R30 rows
                       810, 962 and 1146; N = 58, 42 and 36), made by scripts/s43_certify_ci2.py.  Fields: index (row
                       of the R54 list), ray, P1 = (X1, Y1) and P2 = (X2, Y2) (bit masks of A..F; bit 6 is Z1, so P2
                       may contain Z1), terms (descriptor and rational multiplier: 'H' and 'M' Shannon rows, 'I' a
                       set-valued conditional Ingleton row, 'C' one of the six common-information equalities; the
                       multipliers of the inequality rows are nonnegative, those of the equalities have either sign),
                       F_S (the inequality on the 31 entropies of five parties that the certificate proves, integer
                       coefficients) and F_S_at_ray (< 0).  The combination cancels on every set that meets Z1 = W_X1
                       cap W_Y1 or Z2 = W_X2 cap W_Y2; the form F left on A..F satisfies F(h_norm(S)) >= 0 for every
                       stabilizer state S of every prime dimension (the arrangement extended by the two intersections
                       is again a subspace arrangement), while F(h_norm(r)) < 0, so r is not in the stabilizer cone of
                       any prime dimension.  Checked by epr1kit.stabcert.verify_exclusion_ci2, which rebuilds every
                       row from its descriptor without the row-building code of s42 and s43; s43 runs it before it
                       writes a certificate and in --check, and gate G50 runs it again.  Failed certifications: 0.
                       The three inequalities F_S lie in three different S6 orbits and each takes the value -2 on its
                       ray; every S6 image of each holds on all 1,008 known orbits and on the two qutrit rows below,
                       as it must, and cuts none of the 785 undecided rows (the inequality of R30 row 962 also cuts
                       row 810).
realisations_5.npz     one certificate for each of the 5 rows realised at lambda = 2 by the fixed-structure search
                       (R30 rows 73, 418, 518, 662 and 793; N = 30, 30, 28, 32 and 27 qubits at lambda = 1), in the
                       format of results/2026-11-23_r54-sandbox/realisations_8.npz with row54 (index in the R54 list)
                       in place of row53: row, row54, rep, lam (2), sizes, rows (adjacency rows of a bipartite graph
                       state), cols and m (the matrix M), source ('css-fixed').  Each row has three parties whose
                       subspaces the rank data force into a normal form (pairwise intersections one common subspace T,
                       sum T plus a direct sum; proved for any field); the search fixed them as coordinate subspaces
                       and annealed the other three (scripts/s44_fixed_search.py, kernel cc/s44_fixed_anneal.c; rows
                       73 and 418 were found first by a sandbox prototype driver of the same kernel and again by s44;
                       provenance in status.json).  Re-checked: s41's verify (column-group ranks, cut ranks of the
                       graph, stabcert) and stabcert on the stored graph; each row is an extreme ray of QLR5 (tight
                       rank 30) in an S6 orbit that is none of the 1,003 known before.
realisations_gf3.json  realisations over GF(3) at lambda = 2 of R30 rows 807 and 973 (rows 535 and 651 of the R54
                       list), the two undecided orbits with the fewest qubits: a matrix M over GF(3) whose column
                       groups have the rank function h_S (stored by columns, p:q bit masks of the entries 1 and 2) and
                       the weighted bipartite qutrit graph state of M row-reduced to [I | A].  Checked: the 63
                       column-group ranks over GF(3) (stabcert.gfp_rank) and the entropy vector of the graph
                       (stabcert.verify_realisation_gfp).  Both rows are extreme rays of QLR5, which contains the
                       qutrit stabilizer cone, so they are extreme rays of that cone.  Whether qubits realise them is
                       open; they stay in undecided_reps.npy.
undecided_reps.npy     the other 785 rows of the R54 list, in the same order (N <= 21: 2, N = 22-32: 257, N = 33-40:
                       372, N = 41-64: 154); 183 of them have a forced triple.
status.json            the counts above, by N band; the excluded rows with the pairs of their certificates; the run of
                       ci2-ingleton (rows tested, LPs, solver trouble); the provenance of each realisation (driver,
                       seed, seconds to the hit, the forced triple); the qutrit rows; the status for every local
                       dimension carried over from R54 with the new rows open; the frontier.
frontier_1008.npz      the frontier of the cone generated by the S6 images of the 1,008 known orbits (438,122 rays).
                       14 of the 340 classes of results/2026-11-23_r54-sandbox/frontier_1003.npz are violated by a row
                       realised here, so they do not hold on Stab5; 326 classes remain; for the undecided rows that
                       none of them cuts, the linear program of the S6 sandbox gave a facet of the cone violated by
                       the row (8 new classes).  Fields: classes, index_1003 (-1 for a new class), tight, cut (334 x
                       785), rows (R30 rows of the undecided rows), greedy_cover (182 classes cutting all 785 rows).
                       Every class was re-checked in exact integer arithmetic to be valid on all rays with a tight set
                       of rank 30.

Re-check: python3 selftest.py (gate G50).
The scripts of the sandbox steps are not in the repository (sandbox bundle EPR1_r56_sandbox_2026-11-23.zip,
sha256 b019d011f5c025079d80360e64b0c7093d364e6e0966c0732cdf4d03eb98b28f, kept by the authors).
