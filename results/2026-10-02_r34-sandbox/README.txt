R34 sandbox runs (2026-10-02, 2-core container; no server time used).

s29_zd.log         scripts/s29_zd_stabilizer.py --cases 2:4 2:2 2:3 4:8 4:4 3:9 3:3 6:4 6:2 6:3 --state-vector
                   --lemma-tests 2200000
                   Generalised cube states over Z_N (units of log p): every case free, unit weights, identically
                   self-dual, h_S = 2r, contraction = J - I configuration over Z_N.
                   Pulled-back Lemma 4(G) values (LA, LB):
                     Z_4 t=2: (-2, -2)   [outside the qubit cone and every odd-p cone]
                     Z_2 t=2: (0, -2)    Z_3 t=2: (-2, 0)    [the R32 cube states]
                     Z_8 t=4: (-2, -4)   Z_4 t=4: (0, -4)    [(LA)_4 separates Z_8 from Z_4]
                     Z_9 t=3: (-2, -2)   Z_3 t=3: (0, -2)    [(LA)_3 separates Z_9 from Z_3]
                     Z_4 t=6: (-2, -2)   Z_2 t=6: (0, -2)    Z_3 t=6: (0, -2)   [(LA)_6 separates Z_4 from Z_6]
                   Z_4 cube state against its state vector (4^8 amplitudes, 254 cuts): deviation 8.9e-16 bits;
                   Z_4 cube (bits) = qubit cube (bits) + qutrit cube (trits): True.
                   Lemma 4(G) on point configurations: (LA)_2, (LB)_2 over Z_4^3 (all 2,097,152): min -1 (not
                   claimed there); (LB)_3 over Z_4^4, (LA)_4 over Z_4^5, (LB)_2 over Z_9^3 (2.2 M random each):
                   min 1, 0, 0 (claimed valid); (LA)_3 on J - I over Z_4 and (LA)_2 on J - I over Z_9: -2 (not
                   claimed there).
selftest_full.log  python3 selftest.py --full: 38 gate lines, ALL PASS (321 s).

An independent referee (separate agent, own code) re-derived Theorem 1(d), Lemma 4(G), the witnesses and
Theorem 5 and stress-tested Lemma 4(G) (about 1.7 M evaluations); its scripts are not part of the repository.
