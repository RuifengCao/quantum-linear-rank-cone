R30 sandbox runs on the undecided orbits (2-core sandbox, 2026-10-01 BST; input: results/2026-10-02_a2-status).

small10.npy           the 10 undecided orbits with N <= 20 qubits at lambda = 1, in the order of rows
                      273, 396, 434, 484, 517, 807, 970, 973, 981, 995 of
                      results/2026-10-02_a2-status/undecided_reps.npy (row 981 = R27 frontier ray 294)
s23_small10.json/.log scripts/s23_ci_dfz.py on those 10 rows: CI + Ingleton + 128,040 DFZ instances
                      involving Z, every pair: all feasible (floating-point verdicts)
s23_sample16.json/.log the same on 16 rows drawn at random (numpy default_rng(20261001)) from the
                      orbits with N > 20: all feasible
s10_*.npz.json        realisation searches on small10.npy (keys "row:lambda" of best residuals):
                      q2_l2 = qubits lambda 2, 150 s; q3_l2 = qutrits lambda 2, 120 s;
                      q7_l1 = p 7 lambda 1, 60 s; q5_l2 = p 5 lambda 2, 60 s; none realised
search_small10.log    console output of the searches (the first p = 3 run was stopped and restarted
                      with two threads and 120 s)
