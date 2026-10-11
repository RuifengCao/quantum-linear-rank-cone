R37 runs (2026-10-11, 2-core container; no server time used).

The verification scripts of the paper, run from checks/ of the R37 tree without EPR1_REPO (so they find data/ and
results/ through the default path), and the full self-test. Python 3.11.15, numpy 2.4.4, scipy 1.17.1, sympy 1.14.0,
lrslib 0.71b, gcc 13.3.0, Linux x86_64. Times are wall-clock seconds on this container.

check_sec3.log                      python3 check_sec3.py                          ALL OK            22 s
check_sec3_lift.log                 python3 check_sec3_lift.py                     ALL OK             6 s
check_sec4a.log                     python3 check_sec4a.py                         ALL OK            41 s
check_sec4b.log                     python3 check_sec4b.py                         ALL OK            27 s
check_sec4b_facet.log               python3 check_sec4b_facet.py                   1,414 tight, rank 30   34 s
check_dfz6.log                      python3 check_dfz6.py                          ALL OK            70 s
check_facets119.log                 python3 check_facets119.py                     117 facets        40 s
check_facets119_with_literature.log python3 check_facets119.py --with-literature   119 facets        40 s
check_realizations.log              python3 check_realizations.py                  ALL OK             4 s
check_sec5.log                      python3 check_sec5.py                          ALL OK            27 s
check_sec6.log                      python3 check_sec6.py                          ALL OK            22 s
selftest_full.log                   python3 selftest.py --full                     44 gate lines, ALL PASS (332 s)

checks/README.md lists which statements of the paper each script covers.
