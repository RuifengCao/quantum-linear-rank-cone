R44 sandbox runs (2026-10-15; 2-core container; no server time used).  Theorem 14 of
docs/R32_reduction_and_local_dimension.md section 4.12: f^(a) = f12 + (1 - 1/a) U_Z >= 0 on the Z_{2^a} stabilizer
cone at six parties, f^(a)(S_NF) = -1/a.  Theorem 14 is a proof; this folder holds only checks.  All values are
integers: entropies in bits and a f^(a) = a f12 + (a - 1) U_Z.

s35_check_full.log       python3 scripts/s35_two_power_transfer.py --check --full --trials 400: the checks of gates
                         G40 and G40b with 400 random states over Z_4 and 400 over Z_8 (every step of the proof,
                         chain: 0 failures, minimum a f^(a) = 0, attained 39 and 27 times; minimum f12 = -2 and -4),
                         the entropy formula against state vectors (largest deviation 1.3e-15), Lemma 4(G') tests
                         (valid forms minimum 0; the form for a smaller exponent -1), and the two exhaustive families
                         of the gate.
lifts3.log               r44_runs.py lifts3: cc/s38_minform_zn.c in lifts mode, sixteen families of 2^21 codes
                         B + 2^s P (3 x 7, P a 0/1 matrix, one qudit per system): all 0/1 codes over Z_4, Z_8,
                         Z_16; the lifts seven-point + 2P (Z_4, Z_8), + 4P (Z_8), + 8P (Z_16); the Z_4 seven-point
                         code inside Z_8 + 4P; and eight random 0/1 bases + 2P over Z_4.  Per family: the minimum of
                         a f^(a) (0, or positive) and of four forms the proof does not give (U_Z coefficient lowered
                         by 1/(2a), e-terms dropped, 3 S(Z) in U_Z, U_W instead of 4 U_W), each minimiser rebuilt and
                         recomputed in Python, and every step of the proof at the minimiser of a f^(a).
lifts4.log               r44_runs.py lifts4: all 2^28 codes 4 x 7 with 0/1 entries over Z_4.
anneal.log               r44_runs.py anneal: cc/s38_minform_zn.c on CSS codes over Z_4 and Z_8, five runs of two
                         processes, 120 s each: one qudit per system unpulled; pulled towards full single-system
                         entropies; pulled towards S^(2^a); two qudits per system pulled towards 2 S^(2^a); unequal
                         numbers of qudits per system over Z_4, pulled towards full single-system entropies.  6.5e7
                         moves in all; minimum 0; every minimiser passed every step of the proof.
selftest_full.log        python3 selftest.py --full on a clean clone of commit 88279f7 (R42) with the R43 and R44
                         changes copied in, Python 3.12.
r44_runs.py              the script of lifts3.log, lifts4.log and anneal.log.
