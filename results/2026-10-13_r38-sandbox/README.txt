R38 sandbox runs (2026-10-13, about 4 h; 2-core container, 7 GB; no server time used). Roadmap item S5(a): unimodular
realisations of the 641 known extreme rays of Stab_5 (docs/R32_reduction_and_local_dimension.md section 4.9).

Result: 577 of 641 rows certified for every local dimension d: 564 by unimodular realisations
(kappa = 1: 513, 2: 50, 3: 1; data/s32_stab5_unimodular.npz) and 13 by covering certificates
(data/s32_stab5_cover.npz).  Open: 64 rows (200, 488, 492, 493, 495, 496, 497, 501, 502, 503, 504, 506, 508, 509, 525, 539, 549, 556, 557, 559, 560, 564, 569, 572, 573, 574, 575, 576, 578, 579, 580, 582, 583, 584, 585, 586, 587, 588, 589, 590, 591, 593, 594, 595, 596, 597, 598, 605, 606, 607, 609, 610, 611, 612, 613, 630, 631, 633, 634, 635, 636, 637, 638, 640).

campaign.log            the passes in order (scratch orchestration around cc/s32_unimodular.c; the same strategies
                        are available as scripts/s32_unimodular.py --search): jobs, new certificates, wall time.
                        Pass 1 (plain, from the qubit certificates, 20 s) and pass 2 (plain, 15 s) covered the 629
                        rows with qubit certificates; H1 (neutral moves, 10 s) all open rows; H2 (Eisenstein ring,
                        15 s, stopped early); qutrit-start and gf3-search-start runs on the rows that qubits realise
                        only at multiplicity 2 and on the 12 ledger rows; H3 (Eisenstein from qutrit certificates,
                        300 s; neutral, 45 s) and H4 (neutral, 90 s; stopped after 34 of 87 rows to free the cores
                        for the covering run) on the rest.  "CERT-FAIL 607" is the report with Smith invariant 127
                        that the exact check rejected (the kernel then checked primes up to 97 only; now below 1000).
summary.json            counts for the unimodular file: kappa, strategy, the 12 ledger rows, rows by size N.
cover_run.log           python3 scripts/s32_unimodular.py --cover --seconds 15 --seed 2001 (rounds, rows covered)
cover_raw.txt           its raw kernel output (FOUND lines carry the rational realisations)
cover_open.json         for every row still open: its stored GF primes, the bad-prime sets of the rational
                        realisations found (covering run and probes/badprimes_example.txt), the primes still uncovered.
check_cover.log         python3 scripts/s32_unimodular.py --check-cover
check_full.log          python3 scripts/s32_unimodular.py --check --full (M1, M2, GF(p) p <= 19, Z_4 / Z_6 state vectors)
selftest_full.log       python3 selftest.py --full (47 gate lines)
s32_selftest.log        python3 scripts/s32_unimodular.py --self-test
ci_oddchar.json/.log    python3 scripts/s33_ci_oddchar.py --rows <the rows open before the covering run> --pairs disjoint
                        (exploratory; floating-point verdicts)
probes/strategy_sample20.txt   success at 10 s on a random sample of 20 hard rows, per search setting
probes/badprimes_example.txt   three rational realisations of rows 456, 488, 636 and their exact bad primes (the
                        observation behind the covering certificates)
probes/gf3_lambda1_sample6.txt GF(3)-only searches at multiplicity 1, 60 s, on six rows (456, 488, 525, 575, 591,
                        636): none found, although row 636 is realised mod 3 by its rational realisation above
