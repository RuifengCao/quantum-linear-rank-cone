E2 at maxsize 2 completed (R57), sandbox.  Inputs: results/2026-11-24_ci2-ingleton2/ (job ci2-ingleton2, 4.7 h on a
32-thread laptop on 2026-11-24; it resumed ci2-ingleton and tested the 360 rows that run had not reached, so all 793
rows of results/2026-11-23_r54-sandbox/undecided_reps.npy (the R54 list) are now tested at maxsize 2 with nested
pairs: 14,362,340 LPs in all, none with solver trouble), results/2026-11-23_r56-sandbox/ and data/.  Nothing in data/
changed and gates G0-G50 are as before; gate G51 checks this folder and the job folder.  Counts (status.json):
1,915,914 excluded orbits of the R30 list (1,915,913 + 1); 1,008 known orbits of extreme rays of Stab5; 785 - 1 = 784
undecided.

certificates_ci2.json  one exact certificate, for the one row that ci2-ingleton2 reported infeasible and R56 had not
                       certified: row 591 of the R54 list (R30 row 885, N = 32; P1 = (A, BC), P2 = (B, AF)), in the
                       format of results/2026-11-23_r56-sandbox/certificates_ci2.json (69 terms; F_S takes the value
                       -2 on the ray).  Checked by epr1kit.stabcert.verify_exclusion_ci2, which rebuilds every row
                       from its descriptor without the row-building code of s42 and s43; s43 runs it before it writes
                       a certificate and in --check, and gate G51 runs it again.  The other three rows the job
                       reported infeasible, 537, 643 and 785, are the three of R56, with the pairs of their R56
                       certificates.  The new inequality lies in an S6 orbit different from the three of R56; every S6
                       image of it holds on all 1,008 known orbits and on the two qutrit rows of R56, as it must, and
                       none cuts any of the 784 undecided rows.
undecided_reps.npy     the other 784 rows of results/2026-11-23_r56-sandbox/undecided_reps.npy, in the same order
                       (N <= 21: 2, N = 22-32: 256, N = 33-40: 372, N = 41-64: 154).
status.json            the counts; the excluded row with the pairs of its certificate; the run of ci2-ingleton2 (rows,
                       LPs, time per LP, the seven rows realised in R56, all feasible as they must be); the R30 rows
                       of the undecided list; the frontier.  Frontier: no new realisation, so the 334 classes of
                       results/2026-11-23_r56-sandbox/frontier_1008.npz stand; restricted to the 784 remaining rows,
                       every row is still cut and the greedy cover of R56 (182 classes) still covers them all (a new
                       greedy cover needs 181).

Re-check: python3 selftest.py (gate G51).
The build script is not in the repository (sandbox bundle EPR1_r57_sandbox_2026-11-24.zip,
sha256 08c0f38c75e0d4592e9fe469787658c9acf786530189e0df420bed0e4c2d446e, kept by the authors).
