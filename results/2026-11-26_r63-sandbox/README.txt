Laptop jobs of 2026-11-25 and 2026-11-26 (R63), sandbox: the first piece of E2 at maxsize 3, the last two pieces of
css-sweep3 and the four pieces of s32-new5.  Inputs: results/2026-11-25_ci3/ (job ci3-8h, 7.9 h on a 32-thread laptop:
scripts/s42_ci2_ingleton.py --maxsize 3 on the 777 rows of results/2026-11-24_r60-sandbox/undecided_reps.npy, first
piece of a resumable run), results/2026-11-26_css-sweep3d/ and css-sweep3e/ (1.0 h each: scripts/s41_css_search.py, one
search of 3,600 s on each of the last 46 of the 142 rows with N <= 30 of the R57 list), results/2026-11-26_s32-new5a/
... s32-new5d/ (1.0 h each: scripts/s32_unimodular.py, the unimodular search for every local dimension on the rows
realised after submission), results/2026-11-24_r60-sandbox/ and data/.  Nothing was realised or excluded: 1,915,914
excluded / 1,015 known / 777 undecided orbits of the R30 list, as in R60; the undecided list and the frontier of R60
stand.  Nothing in data/ changed and gates G0-G53 are as before; gate G54 checks this folder and the seven job folders.

status.json            the counts (unchanged); every_d: 114 of the 374 orbits realised after submission now have a
                       certificate for every local dimension d (111 from s32-new4, plus R30 rows 588 from s32-new5a and
                       313 and 361 from s32-new5b; s32-new5c and s32-new5d found none), 260 are open; ci3: the first
                       piece did 283 of the 777 rows (N = 38-64, largest first; 8,193,065 LPs, none with solver
                       trouble), every one feasible, so no candidate exclusion yet, and 494 rows (N = 16-38) are left
                       for the next pieces; css_sweep3_de: no row realised (best residuals per piece, the comparison
                       with css-sweep2, R30 row 979 of piece d already realised in R59); the job folders with their run
                       times.
ci3 (its output is in results/2026-11-25_ci3): re-checked here row by row: done = combos = the count of s42's own
enumeration at maxsize 3 without nesting, no solver trouble, no infeasible row, the summary equal to the rows.  The
verdicts are floating point; a candidate exclusion would still need an exact certificate (scripts/s43_certify_ci2.py).
The run continues on the laptop from ci3.json (jobs ci3-3h ... ci3-10h).
s32-new5 (its output is in the four job folders): each file holds the 111 certificates of s32-new4 and the new ones of
its piece; all four re-checked here with the jobs' own --check --full (the catalogue extended by the six realisation
files of R49-R54).  The three new certificates were also checked with independent code: W symmetric with zero diagonal,
party sizes kappa times the singles, and for each of the 31 sets X of parties the cut block has Smith invariants in
{0,1} (sympy, over ZZ) and rational rank kappa r(X), so W mod d is a graph state of Z_d qudits with entropy vector kappa
r for every d >= 2.
css-sweep3 (output in the job folders d and e): 46 rows, none realised; against css-sweep2 the best residual is lower on
11, the same on 24 and higher on 11 [empirical].  With pieces a-c (R60) all 142 rows with N <= 30 of the R57 list were
searched once more: 4 realised.

Re-check: python3 selftest.py (gate G54).
The build script is not in the repository (sandbox bundle EPR1_r63_sandbox_2026-11-26.zip, sha256
80d7a2e746c9ead936acf7f26f566a82c0d878b1692d6f438b099e1f80417ada, kept by the authors).
