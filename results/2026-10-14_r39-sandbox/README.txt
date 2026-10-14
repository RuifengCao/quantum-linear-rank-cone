R39 sandbox runs (2026-10-13 23:00 to 2026-10-14 03:43; 2-core container, 7 GB; no server time used). Roadmap item
S5(a), continued: the 64 known extreme rays of Stab_5 that R38 left open (docs/R32_reduction_and_local_dimension.md
section 4.9).

Result: 636 of 641 rows certified for every local dimension d: 566 by unimodular realisations
(data/s32_stab5_unimodular.npz) and 70 by covering certificates (data/s32_stab5_cover.npz).  R39 closed
59 of R38's 64 open rows (2 unimodular, 57 covered).  Open: 5 rows (587, 589, 612, 630, 635).

campaign.log            the steps in order (scratch orchestration around cc/s32_unimodular.c; the same kinds of
                        searches are available as scripts/s32_unimodular.py --cover, R39 order).  Step names:
                          A    rational searches that must also be right mod every uncovered prime of the row (at
                               most four), started from the row's rational realisation with the fewest uncovered
                               primes; 40 s per row and round
                          B    rational searches (residual mod 2^31 - 1 only), plain, from the qubit certificate at
                               multiplicity 1 (Bl: the same, longer, with the best matrix dumped)
                          Bw, Bi   the same in ring mode, Z[omega] and Z[i] (kappa 2); Ba2: entries in {-2, ..., 2}
                          Bd   rational searches at multiplicity 2 from two copies of the qubit certificate (row 200:
                               from its qubit certificate at multiplicity 2)
                          final push   Bl at 900 s on the last six rows without a rational realisation, beside
                               searches right mod 3 and 5 for row 576 (300 s each); then A at 60 s for row 594
                        Every reported matrix had its bad primes recomputed exactly before it entered the state; a
                        report with wrong rational ranks never occurred.  Before the campaign: the state was built from
                        every kernel report of R38 (rows 539 and 631: unimodular; row 574: covered) and a calibration
                        (probes/calibration_p3.txt) covered rows 488, 493, 496, 613.  Row 564 was finished by side runs
                        (probes/row564.txt).
raw_kernel_output.txt   the raw kernel output of every call ("#" lines: step, seconds, seed, active primes, ring, alphabet,
                        sweeps, neutral-move probability, extra checks, number of jobs, wall time; FOUND lines carry the
                        matrices; BEST lines the best residual, and with "dump 1" the best matrix)
open_rows.json          for every row still open: N, party sizes, the best residuals of its rational searches, and for the
                        dumped best matrices the cuts whose exact rational rank differs from the target (+: too high)
summary.json            counts: unimodular and covering rows, kappa, pieces, how each of the 64 rows was closed
check_full.log          python3 scripts/s32_unimodular.py --check --full
check_cover.log         python3 scripts/s32_unimodular.py --check-cover
s32_selftest.log        python3 scripts/s32_unimodular.py --self-test
selftest_full.log       python3 selftest.py --full (47 gate lines)
smoke_s32_unimodular.log  the smoke variant of the server job s32-unimodular (jobs.json): passed
ci_all_pairs.json/.log  python3 scripts/s33_ci_oddchar.py --rows <the 22 rows open at 00:20> --pairs all
                        (exploratory; floating-point verdicts; no exclusion)
probes/calibration_p3.txt   rows 488, 493, 496, 613 (every known rational realisation bad at 3), 30 s per job:
                        GF(3) at multiplicity 1 and GF(9) through Z[i] from random starts versus a rational search that
                        is also right mod 3, started from the known rational realisation
probes/row564.txt       row 564: four more rational searches (three found, all bad at 3), GF(3) and GF(9) searches
                        (none), then rational searches right mod 3 from its four rational realisations, 300 s (three found)
