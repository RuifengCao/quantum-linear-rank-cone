Laptop job of the night of 2026-11-27/28 (R73), sandbox: the third and last piece of E2 at maxsize 3.  Input:
results/2026-11-28_ci3/ (job ci3-6h, 4.7 h on a 32-thread laptop, run on the R71 commit: scripts/s42_ci2_ingleton.py
--maxsize 3 --hours 5.9 --resume on results/2026-11-24_r60-sandbox/undecided_reps.npy, continuing ci3.json of
results/2026-11-27_ci3/), results/2026-11-27_r71-sandbox/ and data/.  Nothing in data/ changed and gates G0-G59 are as
before; gate G60 checks this folder and the job folder.  Nothing is realised or excluded here: 1,915,914 excluded /
1,024 known / 768 undecided orbits of the R30 list, as in R71.

status.json            the counts of R71; every_d of R71; the ci3 block with the three pieces (283 rows, N = 38-64,
                       8,193,065 LPs; 308 rows, N = 32-38, 8,819,416 LPs; 186 rows, N = 16-32, 5,273,482 LPs in 16,834
                       s): all 777 rows of the R60 list done, 22,285,963 LPs, every row feasible, no solver trouble,
                       complete.  The 768 rows undecided in R71 are among them; the 9 rows of the R60 list realised in
                       R64, R66, R69 and R71 came out feasible, as realised rows have to.
So E2 at maxsize 3 without nesting (s42: the Shannon and set-valued Ingleton rows on the common-information extensions
of up to three variables) excludes none of the 768 undecided rows [computed; floating-point LP verdicts].  A feasible
row is not excluded by these LPs, which says nothing about whether it is realisable.
Re-checked: the job folder (manifest, exit code, no warnings, the log resumes from the 591 rows of the first two pieces
and ends COMPLETE); every row of the earlier pieces carried over unchanged; all 777 rows feasible with done = combos =
s42's own count of LPs at maxsize 3 without nesting and N = h(all); the summary equal to the rows and complete.

Re-check: python3 selftest.py (gate G60).
The build script is not in the repository (sandbox bundle EPR1_r73_sandbox_2026-11-28.zip, sha256
38d15c28dfbcfae2f5c44bd48bd40789221345b451d46fb6818fab74a96e7108, kept by the authors).
