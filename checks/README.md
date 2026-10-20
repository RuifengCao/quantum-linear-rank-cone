# Independent checks for the paper

These scripts recompute the numbers quoted in the paper by X. Song and R. Cao on stabilizer entropy cones
(arXiv:2610.18905 and its revised version of October 2026). They read the data of the repository
https://github.com/RuifengCao/quantum-linear-rank-cone (release R47) but not its code, except for the
transcriptions of the published seven-variable inequalities used in `check_sec6.py` and `check_rev_qss.py`.

Requirements: Python 3 with numpy, scipy and sympy; `lrs` (lrslib) for `check_sec4a.py`; a C compiler (`gcc`) for
`q2_gf3.c`, which `check_sec3.py` compiles and runs.

Run the scripts from the directory `checks/` of the repository, where they find the data automatically, or set
`EPR1_REPO` to the path of a checkout of the repository, e.g.

    EPR1_REPO=/path/to/quantum-linear-rank-cone python3 check_sec3.py

Every script ends with `ALL OK` or `FAILURE`, except `check_sec4b_facet.py` and `check_facets119.py`, which print
counts.

The table gives the parts of the paper that each script checks, in the numbering of the revised version (Results
Secs. II A–F, Discussion Sec. III, Methods Sec. IV, Supplementary Notes 1–10) and, in brackets, the sections of
arXiv:2610.18905v1. Numbers of theorems, lemmas, remarks and equations refer to the revised version.

| script | paper |
| --- | --- |
| `check_sec3.py`, `q2_gf3.c` | Sec. II B, Example 1, Supplementary Note 1 [v1: Sec. III, Example 1, Appendix A] (CSS formula and doubling against state vectors) |
| `check_sec3_lift.py` | Sec. II B, Supplementary Note 1 [v1: Sec. III, Appendix A] on random stabilizer states of Z_d qudits, d = 2, 3, 4, 6, 8, 9, 12, including non-free Lagrangian subgroups: entropy formula, the lift as a subgroup rank function, self-duality, doubling |
| `check_sec4a.py` | Theorem 2, Supplementary Note 8, Sec. II C 1 [v1: Theorem 2, Appendix G, Sec. IV A] on the state (R6) of Linden et al. (as printed, S(ab) = 1 instead of 2 on ray 6 of their Table 1, and its entropy vector is not an extreme ray) and the qubit state of Osthus for that ray |
| `check_sec4b.py`, `check_sec4b_facet.py` | Theorem 3, Sec. II C 3, Supplementary Notes 2.1-2.3 [v1: Theorem 3, Sec. IV C, Remark 2, Appendix B] |
| `check_dfz6.py` | the six-variable inequalities (44)-(63) of Dougherty, Freiling and Zeger (arXiv:0910.0284v3) with the purifier as a variable: Sec. II C 2-3, Supplementary Note 2.4 [v1: Sec. IV B-C, Appendix B 4] |
| `check_facets119.py` | facet count of Sec. II C 3 and Supplementary Note 2.2 (117 classes; 119 with the 12 rays taken from the literature) [v1: Sec. IV C] |
| `check_realizations.py` | qubit and qutrit realizations of the 641 known extreme-ray orbits of Stab_5: Sec. II C 3, Supplementary Note 2.2, Discussion [v1: Sec. IV C, Appendix B 2, Sec. VII] (every stored graph-state certificate recomputed over GF(2), GF(3) or GF(4); the counts 629/12, 96 = 16 + 80, 67/13, 73/7, 11, 26 and 122; the record of the R30 qutrit search) |
| `check_sec5.py` | Sec. II D, Supplementary Notes 3 and 4 [v1: Sec. V, Appendices C and D] (including Lemma 1 on perturbed configurations over GF(2), GF(3), GF(5), Z_4 and Z_9) |
| `check_sec6.py` | Sec. II E 1-4, Supplementary Notes 5, 6 and 9 [v1: Sec. VI, Appendices E, F and H] (including Theorems 7 and 9 of Pena-Macias and Sarria with n = 2) |
| `check_rev_sec5.py` | not in v1: Lemma 2 (values on the configuration over GF(l) and Z_N for several choices of j(.), random and perturbed arrangements with a negative control, the proof step by step with C computed explicitly), the doubled-code witnesses of Theorem 4 with 2p + 2 parties, and the witnesses of all four cases of Theorem 5 with at most 2 max(d, d') + 2 parties (entropies from the codes); Sec. II D, Supplementary Note 3 |
| `check_rev_sec6.py` | not in v1: the expansion of f_12 (Eq. (S8)) and of U_Z, the values on the Steane, non-Fano and seven-point states (also against state vectors over Z_2, Z_3, Z_4, Z_6), the Smith invariants and values of Corollary 4 for d = 2..32, Theorem 8 and Remark 4; tests of f_12 >= 0 and of every step of the proof of Theorem 7 (Lemmas 7.1-7.5, (LA), Eq. (S7)) on random and structured qubit states in normal form, and of f^(a) >= 0 on Z_4 and Z_8 states; Sec. II E 5, Supplementary Note 7. The group version of the proof of Theorem 8 is tested by gate G40b of the repository, not here |
| `check_rev_qss.py` | not in v1: Lemma 5 and Proposition 5 (secret sharing), Sec. II F, Methods, Supplementary Note 10: the access structures of the cube states and their contractions (the Fano and non-Fano access structures of Jafari and Khazaei), (LB) and (LA) with t = 2 on all seven-element minors, the exact verification of the two stored dual certificates (`data/s38_qss_certificates.json`) with negative controls, and the distances of the witnesses from the cones (Sec. II E 3, Sec. II E 5, Methods) |

`check_rev_sec5.py` and `check_rev_sec6.py` need no repository data and run in under a minute; `check_rev_qss.py` reads
`data/s38_qss_certificates.json` and, for the distances, `data/dfz_char_dependent.csv` and `data/s31_six_party.npz`.
