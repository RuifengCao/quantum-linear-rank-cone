R32 sandbox runs (2026-10-02, 2-core container; no server time used).

s26_self_test.log   scripts/s26_reduction.py --self-test --trials 24
                    CSS entropy formula vs state vectors, and the CSS state of a graph state's Lagrangian
                    code vs 2S: max deviations 4e-16 / 6e-16.
s26_n4.log          scripts/s26_reduction.py --n4 --lrs --realise --out data/s26_stab4_certs.npz
                    n = 4: 820 pulled-back DFZ (1)-(24) rows, 0 not implied by Shannon + Ingleton (LP);
                    lrs: 46 extreme rays of h^-1(Shannon + Ingleton); all 46 realised by qubit graph states
                    (41 at lambda = 1, 5 at lambda = 2) and re-verified.
s27_cube_states.log scripts/s27_cube_states.py --state-vector --random 400 --out data/s27_cube_witnesses.npz
                    cube states: identically self-dual, contraction = DFZ Fano / non-Fano configuration,
                    pulled-back (65) = -2 on the qubit state, pulled-back (91) = -2 on the qutrit state,
                    differ only on ABCZ (2 vs 4), state vectors agree; doubled Fano / non-Fano / T8 codes
                    (Proposition 3): h_S = M + |X| on the visible parties, balanced (65)/(91) = -1;
                    random arrangements / graph states: no violation of the valid forms.
s27_probe6.log      scripts/s27_cube_states.py --probe6 240   (exploratory six-party checks, negative)
                    symmetrised Fano / non-Fano: minima 19 / 16 over all relabellings; annealing best 4.
selftest_full.log   python3 selftest.py --full: 34 gate lines, ALL PASS (266 s).
