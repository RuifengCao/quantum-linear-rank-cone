Laptop job of the night of 2026-11-26/27 (R68), sandbox: the second piece of E2 at maxsize 3.  Input:
results/2026-11-27_ci3/ (job ci3-8h, 7.9 h on a 32-thread laptop, run on the R67 commit: scripts/s42_ci2_ingleton.py
--maxsize 3 --hours 7.9 --resume on results/2026-11-24_r60-sandbox/undecided_reps.npy, continuing ci3.json of
results/2026-11-25_ci3/), results/2026-11-26_r66-sandbox/ and data/.  Nothing in data/ changed and gates G0-G56 are as
before; gate G57 checks this folder and the job folder.  Nothing is realised or excluded here: 1,915,914 excluded /
1,019 known / 773 undecided orbits of the R30 list, as in R66.

status.json            the counts of R66; every_d of R66; the ci3 block: the first piece (R63: 283 rows, N = 38-64,
                       8,193,065 LPs) and this one (308 rows, N = 32-38, 8,819,416 LPs in 28,442 s), all feasible with
                       no solver trouble; 591 of 777 rows done (17,012,481 LPs), 186 left (N = 16-32, 5,273,482 LPs,
                       about 4.7 h at the rate of this piece [estimate]); the rows of the R60 list realised in R64 and
                       R66: row 80 (R30 row 150) came out feasible in this piece, as a realised row has to, and rows 25,
                       334 and 397 (R30 rows 54, 539, 644) are still to do.
Re-checked: the job folder (manifest, exit code, no warnings, the log resumes from the 283 rows of the first piece and
ends INCOMPLETE); every row of the first piece carried over unchanged; every row of ci3.json feasible with done = combos
= s42's own count of LPs at maxsize 3 without nesting and N = h(all); the summary equal to the rows.  The verdicts are
floating point: a feasible row is not excluded by these LPs, which says nothing about whether it is realisable.

Re-check: python3 selftest.py (gate G57).
The build script is not in the repository (sandbox bundle EPR1_r68_sandbox_2026-11-27.zip, sha256
cd1c7486ce5c78be14809ad3eeba13c2672594832f0a62161e35778a8078b926, kept by the authors).
