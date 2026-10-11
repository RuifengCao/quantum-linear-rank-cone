# Independent checks for the paper

These scripts recompute the numbers quoted in the paper "Stabilizer entropy cones as self-dual slices of linear rank cones:
a five-party gap and dependence on the local dimension". They read the data of the repository
https://github.com/RuifengCao/quantum-linear-rank-cone (release R37) but not its code, except for the
transcriptions of the published seven-variable inequalities used in `check_sec6.py`.

Requirements: Python 3 with numpy, scipy and sympy; `lrs` (lrslib) for `check_sec4a.py`; a C compiler (`gcc`) for
`q2_gf3.c`, which `check_sec3.py` compiles and runs.

Run the scripts from the directory `checks/` of the repository, where they find the data automatically, or set
`EPR1_REPO` to the path of a checkout of the repository, e.g.

    EPR1_REPO=/path/to/quantum-linear-rank-cone python3 check_sec3.py

Every script ends with `ALL OK` or `FAILURE`, except `check_sec4b_facet.py` and `check_facets119.py`, which print
counts.

| script | paper |
| --- | --- |
| `check_sec3.py`, `q2_gf3.c` | Sec. III, Example 1, Appendix A (CSS formula and doubling against state vectors) |
| `check_sec3_lift.py` | Sec. III, Appendix A on random stabilizer states of Z_d qudits, d = 2, 3, 4, 6, 8, 9, 12, including non-free Lagrangian subgroups: entropy formula, the lift as a subgroup rank function, self-duality, doubling |
| `check_sec4a.py` | Theorem 2, Appendix G; Sec. IV A on the state (R6) of Linden et al. (as printed, S(ab) = 1 instead of 2 on ray 6 of their Table 1, and its entropy vector is not an extreme ray) and the qubit state of Osthus for that ray |
| `check_sec4b.py`, `check_sec4b_facet.py` | Theorem 3, Sec. IV C, Remark 2, Appendix B |
| `check_dfz6.py` | the six-variable inequalities (44)-(63) of Dougherty, Freiling and Zeger (arXiv:0910.0284v3) with the purifier as a variable: Sec. IV B-C, Appendix B 4 |
| `check_facets119.py` | facet count of Sec. IV C (117 classes; 119 with the 12 rays taken from the literature) |
| `check_realizations.py` | qubit and qutrit realizations of the 641 known extreme-ray orbits of Stab_5: Sec. IV C, Appendix B 2, Sec. VII (every stored graph-state certificate recomputed over GF(2), GF(3) or GF(4); the counts 629/12, 96 = 16 + 80, 67/13, 73/7, 11, 26 and 122; the record of the R30 qutrit search) |
| `check_sec5.py` | Sec. V, Appendices C and D (including Lemma 1 on perturbed configurations over GF(2), GF(3), GF(5), Z_4 and Z_9) |
| `check_sec6.py` | Sec. VI, Appendices E, F and H (including Theorems 7 and 9 of Pena-Macias and Sarria with n = 2) |
