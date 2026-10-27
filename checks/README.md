# Independent checks for the paper

These scripts recompute the numbers quoted in the paper by X. Song and R. Cao on stabilizer entropy cones
(arXiv:2610.18905 and its revised version of October 2026). They read the data of the repository
https://github.com/RuifengCao/quantum-linear-rank-cone (release R48) but not its code, except for the
transcriptions of the published seven-variable inequalities used in `check_sec6.py`.

Requirements: Python 3 with numpy, scipy and sympy; `lrs` (lrslib) for `check_sec4a.py`; a C compiler (`gcc`) for
`q2_gf3.c`, which `check_sec3.py` compiles and runs.

Run the scripts from the directory `checks/` of the repository, where they find the data automatically, or set
`EPR1_REPO` to the path of a checkout of the repository, e.g.

    EPR1_REPO=/path/to/quantum-linear-rank-cone python3 check_sec3.py

Every script ends with `ALL OK` or `FAILURE`, except `check_sec4b_facet.py` and `check_facets119.py`, which print
counts.

The table gives the parts of the paper that each script checks, in the numbering of the revised version (Results
Secs. II A–F, Discussion Sec. III, Methods Sec. IV, Supplementary Notes 1–11) and, in brackets, the sections of
arXiv:2610.18905v1. Numbers of theorems, propositions, lemmas, remarks and equations refer to the revised version.

| script | paper |
| --- | --- |
| `check_sec3.py`, `q2_gf3.c` | Sec. II B, Example 1 (in Supplementary Note 3.2 of the revised version), Supplementary Note 1 [v1: Sec. III, Example 1, Appendix A] (CSS formula and doubling against state vectors) |
| `check_sec3_lift.py` | Sec. II B, Supplementary Note 1 [v1: Sec. III, Appendix A] on random stabilizer states of Z_d qudits, d = 2, 3, 4, 6, 8, 9, 12, including non-free Lagrangian subgroups: entropy formula, the lift as a subgroup rank function, self-duality, doubling |
| `check_sec4a.py` | Theorem 2, Supplementary Note 2, Sec. II C 1 [v1: Theorem 2, Appendix G, Sec. IV A] on the state (R6) of Linden et al. (as printed, S(ab) = 1 instead of 2 on ray 6 of their Table 1, and its entropy vector is not an extreme ray) and the qubit state of Osthus for that ray |
| `check_sec4b.py`, `check_sec4b_facet.py` | Theorem 3, Sec. II C 3, Supplementary Notes 3.1-3.3 [v1: Theorem 3, Sec. IV C, Remark 2, Appendix B] |
| `check_dfz6.py` | the six-variable inequalities (44)-(63) of Dougherty, Freiling and Zeger (arXiv:0910.0284v3) with the purifier as a variable: Secs. II C 2-II C 3, Supplementary Note 3.4 [v1: Sec. IV B-C, Appendix B 4] |
| `check_facets119.py` | the facet count for the qubit cone of v1 (117 classes; 119 with the 12 rays taken from the literature) [v1: Sec. IV C]; the revised version counts facets for every local dimension (`check_rev_stab5.py`) |
| `check_realizations.py` | the qubit, qutrit and GF(4) graph-state certificates of the known extreme-ray orbits of Stab_5, each recomputed [v1: Sec. IV C, Appendix B 2, Sec. VII]; the counts it prints (629/12, 96 = 16 + 80, 67/13, 73/7, 11, 26 and 122) are those of v1, which the revised version replaces by realizations for every local dimension (`check_rev_stab5.py`) |
| `check_sec5.py` | Secs. II D 2 and II E, Supplementary Notes 4 and 8 [v1: Sec. V, Appendices C and D] (including Lemma 1 on perturbed configurations over GF(2), GF(3), GF(5), Z_4 and Z_9) |
| `check_sec6.py` | Secs. II D 3-II D 5, Supplementary Notes 5, 6 and 9 [v1: Sec. VI, Appendices E, F and H] (including Theorems 7 and 9 of Pena-Macias and Sarria with n = 2) |
| `check_rev_sec4.py` | not in v1: Theorem 3, the facet statement for every local dimension: all six-vertex graph states with 0/1 weights over GF(2) and, through rational cut ranks (minors in {0, ±1, ±2}), over every odd prime field; 760 and 835 distinct vectors of rank 31, the same 140 on the hyperplane of (★), rank 30; Sec. II C 2, Supplementary Note 3.1 |
| `check_rev_stab5.py` | not in v1: Sec. II C 3, Supplementary Note 3.2, Discussion: the 641 known extreme-ray orbits of Stab_5 lie in the cone of every local dimension (the 566 unimodular realizations and the 75 covering certificates of the repository, re-verified without Smith forms: rational rank and one nonzero maximal minor of each cut block by fraction-free elimination, then ranks modulo the prime factors of that minor); explicit qubit and qutrit realizations of every orbit; 119 single-CI classes are facets of Stab^(d)_5 for every d, among them the class of (★) |
| `check_rev_sec5.py` | not in v1: Lemma 3 (values on the configuration over GF(l) and Z_N for several choices of j(.), random and perturbed arrangements with a negative control, the proof step by step with C computed explicitly), the doubled-code witnesses of Theorem 7 with 2p + 2 parties, and the witnesses of all four cases of Theorem 8 with at most 2 max(d, d') + 2 parties (entropies from the codes); Sec. II E, Supplementary Notes 4, 8.1 and 8.5 |
| `check_rev_sec6.py` | not in v1: the expansion of f_12 (Eq. (S7)) and of U_Z, the values on the Steane, non-Fano and seven-point states (also against state vectors over Z_2, Z_3, Z_4, Z_6), the Smith invariants and values of Corollary 4 for d = 2..32, Theorem 6 and Remark 2; tests of f_12 >= 0 and of every step of the proof of Proposition 3 (Supplementary Lemmas 7.1-7.5, (LA), Eq. (S6)) on random and structured qubit states in normal form, and of f^(a) >= 0 on Z_4 and Z_8 states; Sec. II D 6, Supplementary Note 7. Also the code states of Remark 5 (R32 Proposition 26): in the CSS states of the binary simplex codes of dimension 4 and 5 and of the p-ary simplex codes of dimension 3 (p = 3, 5, 7), seven qudits have S = r_F, resp. r_NF, and (LB), resp. (LA), with t = 2, balanced over the seven parties as in the proof of Proposition 4, is -1 there (exact ranks; state vectors of the 15-qubit and 13-qutrit states). The group version of the proof of Theorem 6 is tested by gate G40b of the repository, not here |
| `check_rev_qss.py` | not in v1: Lemmas 5 and 6 and Proposition 6 (secret sharing), Sec. II F, Methods, Supplementary Note 10: the access structures of the cube states and their contractions (the Fano and non-Fano access structures of Jafari and Khazaei), (LB) and (LA) with t = 2 on all seven-element minors, the exact verification of the two stored dual certificates (`data/s38_qss_certificates.json`) with negative controls, the scheme with information ratio 7/5 for Gamma_F of Proposition 6 (iii) (part (h): a (lambda, omega) = (6, 1) decomposition, seven relabeled cube schemes and the zero-sum code, over Z_d for d = 3, 5, 7, 9, 11, 13, 15, 21, 25, 27, 45 on all 128 player sets, its CSS state, with a negative control), the scheme with information ratio 7/5 for Gamma_NF with qubits (part (i): the access structures Gamma^S, the GF(4) code that differs from the cube code in two entries, its port and matroid, the seven sub-schemes with lambda = 6 and omega = 1, the combined code over GF(4) and in binary on all 128 player sets, the CSS states, two negative controls), the finite part of the bound 7/5 for decompositions of Gamma_F into ideal schemes (part (j)), the schemes with information ratio 2 (part (f): the binary scheme for Gamma_NF printed in Supplementary Note 10, for a secret of one qubit; also the decomposition scheme for Gamma_F of R48, now superseded), the distances of the witnesses from the cones (Secs. II D 4 and II D 6, Methods), the upper bound 1/12 for the Steane state (Methods; part (g)), and the classification of ideal schemes of Proposition 6 (iv) (Supplementary Note 10.5; part (k): the 2 x 2 x 255 certificates of the forcing lemma in `data/s39_ideal_certificates.json`, rebuilt and verified exactly with negative controls, the normal form over Z_2, Z_4, Z_8 and its kernel, and the cube codes over Z_d, d = 9, ..., 45, as ideal schemes for Gamma_NF), and the counts of Sec. II F and the Introduction (part (l): the correctable erasures of Steane's code, exactly the sets that contain no line, 35 of them beyond the distance, and the average-share bound sum_P S(P) >= (28/3) S(R) from Proposition 6 (i)) |

`check_rev_sec4.py`, `check_rev_sec5.py` and `check_rev_sec6.py` need no repository data and run in under a minute;
`check_rev_stab5.py` reads `data/stab5_extreme_reps.npy`, `data/s32_stab5_unimodular.npz`, `data/s32_stab5_cover.npz`
and `data/s17_ineq_classes.npy` and runs in about a minute; `check_rev_qss.py` reads `data/s38_qss_certificates.json`,
`data/s39_ideal_certificates.json` and, for the distances, `data/s31_six_party.npz`.
