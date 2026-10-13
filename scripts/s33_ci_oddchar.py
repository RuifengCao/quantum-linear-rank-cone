#!/usr/bin/env python3
"""s33: common information + seven-variable inequalities valid in odd characteristic (R38; exploratory screen).

For a five-party ray r and a pair (X, Y) of sets of the six systems A..F, every GF(p) arrangement realising
h_norm(r) (Theorem 1) extends by Z = U_X cap U_Y.  When p is odd, the seven subspaces A..F, Z satisfy, in every
relabelling, the seven-variable linear rank inequalities that hold in odd characteristic:
    DFZ (65) (arXiv:1311.4601, Thm 8.6),  Lemma 4 (LB) with t = 2,  BKL (6.17) as corrected (pitfall 21),
as transcribed in data/dfz_char_dependent.csv, scripts/s28_generalized_cubes.py and data/bkl_char_dependent.csv
(the forms of Proposition 8 that hold in odd characteristic; the validity of the corrected (6.17) rests on BKL's proof
with the dropped terms restored, see the paper, Appendix H).  If the LP

    g on the 127 subsets of A..F, Z;  g = h_norm(r) on A..F;  g(XZ) = g(X), g(YZ) = g(Y), g(Z) = I(X; Y);
    Shannon elemental inequalities on seven elements;  every relabelled instance of the three forms

is infeasible, r lies outside Stab^(p)_5 for every odd p.  The instances (8,400 distinct rows) enter as cutting
planes: solve, add the instances the solution violates (at most 300 per round), repeat.  Verdicts are floating-point
LP verdicts; an infeasible pair would still have to be certified exactly.  R38 ran this on the 115 known extreme
rays of Stab_5 without a unimodular certificate (pairs as in s17 --pairs disjoint): every pair was feasible.

  python3 scripts/s33_ci_oddchar.py --rows 456,488 [--pairs disjoint|all] [--out result.json]
"""
import argparse, itertools, json, os, sys, time
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import csr_matrix, lil_matrix, vstack

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..'))
import s17_ci_test as T
import s30_dimension_profile as D30

NV = 127
ODD_FORMS = ('DFZ (65)', '(LB)_2', 'BKL (6.17) corrected')


def instances():
    """distinct rows of the odd-characteristic forms under all 5,040 bijections of their seven variables onto A..F, Z
    (bits 0..6)."""
    forms = dict(D30.six_forms())
    rows, seen = [], set()
    for name in ODD_FORMS:
        nz = [(m, int(c)) for m, c in enumerate(forms[name]) if c and m]
        for pi in itertools.permutations(range(7)):
            row = {}
            for m, c in nz:
                mm = sum(1 << pi[v] for v in range(7) if m >> v & 1)
                row[mm] = row.get(mm, 0) + c
            key = tuple(sorted((k, v) for k, v in row.items() if v))
            if key not in seen:
                seen.add(key); rows.append(dict(key))
    A = lil_matrix((len(rows), NV))
    for k, d in enumerate(rows):
        for m, c in d.items():
            A[k, m - 1] = c
    return csr_matrix(A)


def test_pair(h, X, Y, INST, cap=300, maxit=60):
    """'feasible' / 'infeasible' / 'undetermined' / 'error<status>' and the number of instances added."""
    Z = 1 << 6
    eq = [({m: 1}, h[m]) for m in range(1, 64)] + [({X | Z: 1}, h[X]), ({Y | Z: 1}, h[Y]), ({Z: 1}, h[X] + h[Y] - h[X | Y])]
    Aeq = lil_matrix((len(eq), NV))
    for k, (d, _b) in enumerate(eq):
        for m, c in d.items():
            Aeq[k, m - 1] = c
    Aeq, beq = csr_matrix(Aeq), np.array([b for _d, b in eq], dtype=float)
    active = []
    for _ in range(maxit):
        Aub = T.A_ELEM if not active else vstack([T.A_ELEM, -INST[active]])
        res = linprog(np.zeros(NV), A_ub=Aub, b_ub=np.zeros(Aub.shape[0]), A_eq=Aeq, b_eq=beq,
                      bounds=[(None, None)] * NV, method='highs', options=T.HIGHS_OPTS)
        if res.status == 2:
            return 'infeasible', len(active)
        if res.status != 0:
            return f'error{res.status}', len(active)
        viol = INST @ res.x
        have = set(active)
        new = [int(i) for i in np.argsort(viol)[:cap] if viol[i] < -1e-7 and int(i) not in have]
        if not new:
            return 'feasible', len(active)
        active += new
    return 'undetermined', len(active)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--rows', required=True, help='comma-separated rows of data/stab5_extreme_reps.npy')
    ap.add_argument('--pairs', default='disjoint', choices=('disjoint', 'all'))
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    reps = np.load(os.path.join(HERE, '..', 'data', 'stab5_extreme_reps.npy')).astype(np.int64)
    t0 = time.time()
    INST = instances()
    print(f'{INST.shape[0]} distinct odd-characteristic instances on A..F, Z ({time.time() - t0:.0f} s)', flush=True)
    out = {}
    for k in [int(x) for x in a.rows.split(',')]:
        h = T.h_norm(reps[k])
        pairs = T.pairs_for(h, a.pairs)
        t1, verdicts, infeasible = time.time(), {}, []
        for X, Y in pairs:
            v, _na = test_pair(h, X, Y, INST)
            verdicts[v] = verdicts.get(v, 0) + 1
            if v == 'infeasible':
                infeasible.append([X, Y])
        out[k] = dict(pairs=len(pairs), verdicts=verdicts, infeasible=infeasible)
        print(f'row {k}: {len(pairs)} {a.pairs} pairs, verdicts {verdicts} ({time.time() - t1:.0f} s)', flush=True)
    tot = {}
    for v in out.values():
        for key, c in v['verdicts'].items():
            tot[key] = tot.get(key, 0) + c
    print(f'summary: {len(out)} rays, verdicts over all pairs {tot}, rays with an infeasible pair: '
          f'{sum(1 for v in out.values() if v["infeasible"])}')
    if a.out:
        json.dump(dict(pairs=a.pairs, forms=list(ODD_FORMS), instances=int(INST.shape[0]), rays=out), open(a.out, 'w'))


if __name__ == '__main__':
    main()
