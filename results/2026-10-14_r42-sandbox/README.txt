R42 sandbox runs (2026-10-14 evening; 2-core container; no server time used).  Roadmap item M3, finished: six
parties, the other direction for the closed cones (docs/R32_reduction_and_local_dimension.md section 4.11,
Theorem 12).  Theorem 12 is a proof; this folder holds only checks and an adversarial search.

transfer_check_400.log  python3 scripts/s34_six_party_survey.py --transfer --full --trials 400: f12(S_NF) = -1,
                        f12(S_F) = 0, 31 terms; on the GF(2) list (10,773 vectors) minimum 0, attained 38 times; on
                        the GF(3) list (27,334 vectors) one negative value, at S_NF.  Every step of the proof
                        (transfer_chain: dim L_T = k(T); the restrictions of the line pieces, their dimensions and
                        their orthogonality; Lemmas T1-T4; (LA)_2 on the arrangement of Lemma T5; f12 >= 0) on 400
                        states: 200 random compressed qubit graph states (1-3 qubits per system), 100 compressed CSS
                        states of random GF(2) families pushed towards the non-Fano dependencies, 100 copies (one or
                        two) of the Steane state after random local symplectic maps.  No step failed; minimum f12 0,
                        attained by 101 states.  The dual arrangement of Lemma T5 built explicitly for 8 states, and
                        the compression (symplectic reduction) checked against 8 uncompressed states.
transfer_anneal.log     nine runs of cc/s37_minform7.c (commands: run_anneal.sh): simulated annealing minimising f12
                        over qubit graph states with the given numbers of qubits per system; the last two runs are
                        pulled towards 2 S_NF and 3 S_NF (weight 1 per unit of L1 distance).  For each run the wrapper
                        rebuilds the graph with the smallest f12, recomputes its entropies exactly, reduces it to S(P)
                        qubits per system and runs every step of the proof on it.  Empirical only.
run_anneal.sh           the commands of transfer_anneal.log (run from this folder or from the repository root).
selftest_full.log       python3 selftest.py --full (51 gate lines, G39 and G39b new), run with Python 3.12 on a clean
                        clone of commit e2f74fd (R41) with this round's files copied in.
