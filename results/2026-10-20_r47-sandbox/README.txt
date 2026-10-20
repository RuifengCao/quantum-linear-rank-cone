R47 sandbox runs (2026-10-20; 2-core container; no server time used).  Sections 4.15 and 4.16 of
docs/R32_reduction_and_local_dimension.md: stabilizer quantum secret sharing (Lemma 20, Proposition 21) and the
distances of the witnesses from the cones.  The statements are proofs and exact rational certificates; this folder
holds the computations that check them.

s38_solve.log          python3 scripts/s38_secret_sharing.py --solve: gate G43b.  The two linear programs for Gamma_F
                       with Lemma 4 (LB), t = 2, on all minors (pure; mixed with an environment) re-solved with HiGHS,
                       optimum 43/41 with an exactly feasible rounded point; the two controls ((LA) in place of (LB);
                       Shannon inequalities only) with optimum 1; the programs with the Jafari-Khazaei
                       inequalities on contractions, optimum 4/3 with exact feasible points, with and without
                       self-duality; then the exact checks of gate G43 (contractions, access structures,
                       certificates, minima over all minors, margins).
selftest_full.log      python3.12 selftest.py --full on a clean copy of the working tree (R46, commit 467f6c3, with the
                       R47 changes; lrs installed): 59 gate lines, all pass (1036 s on two cores).  The 57 gate
                       lines that the R46 log also has are identical to it, except G29b, which here ran the exact
                       lrs enumeration (46 rays; skipped in the R46 log); G43, G43b are new.

The certificates in data/s38_qss_certificates.json were generated with --certify (HiGHS duals repaired to exact
rational multipliers with SymPy).  The probes that preceded them (floating-point linear programs from 2026-10-16 on,
with DFZ (65), DFZ (91) and f9; superseded, see the CHANGELOG) were run outside the repository; everything the note
uses is recomputed by scripts/s38_secret_sharing.py.  The independent referee's checker is checks/check_rev_qss.py.
