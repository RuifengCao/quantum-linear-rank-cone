#!/usr/bin/env python3
"""s25: two simultaneous common informations (exploratory, R31).

For a ray r and two pairs of disjoint party sets (X1, Y1), (X2, Y2), every stabilizer state in normal form
extends its arrangement by Z1 = W_X1 cap W_Y1 and Z2 = W_X2 cap W_Y2 at the same time: an 8-element
linear polymatroid on A..F, Z1, Z2 with g(X_i Z_i) = h(X_i), g(Y_i Z_i) = h(Y_i), g(Z_i) = I(X_i; Y_i).
The LP "extend h_norm(r) that way and satisfy every Shannon inequality on the 8 elements" (only the
1,560 elemental rows that touch Z1 or Z2 are needed) is solved for every unordered pair of pairs with
|X_i|, |Y_i| <= --maxsize and I(X_i; Y_i) > 0.  Dougherty-Freiling-Zeger (arXiv:0910.0284 Section 6)
report hundreds of six-variable linear rank inequalities that need two common informations; inequality
(61) there uses the pairs (A, B) and (E, DF), which this test covers.

  python3 scripts/s25_two_ci.py rays.npy --index 273 981 [--maxsize 2] [--out out.json]

Floating-point LP verdicts.  R31 ran the 10 smallest undecided orbits and 6 larger ones: all feasible
(results/2026-10-02_r31-sandbox/).  No certificate path exists yet for an infeasible pair of pairs.
"""
import argparse, itertools, json, os, sys, time, warnings
import numpy as np
from scipy.optimize import linprog, OptimizeWarning
from scipy.sparse import lil_matrix, csr_matrix

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import s17_ci_test as T

warnings.filterwarnings('ignore', message='Unrecognized options detected|Unknown solver options', category=OptimizeWarning)
NV8, Z1, Z2 = 255, 1 << 6, 1 << 7


def _ub_matrix():
    rows = [r for r in T.elemental_rows(8) if any(m & (Z1 | Z2) for m in r)]   # rows on A..F alone hold already
    A = lil_matrix((len(rows), NV8))
    for k, d in enumerate(rows):
        for m, c in d.items():
            A[k, m - 1] = -c
    return csr_matrix(A)


A_UB = _ub_matrix()


def small_pairs(h, maxsize=2):
    subs = [m for m in range(1, 64) if bin(m).count('1') <= maxsize]
    return [(X, Y) for X in subs for Y in subs if Y > X and not X & Y and h[X] + h[Y] - h[X | Y] > 0]


def status(h, P1, P2):
    rows, b = [], []
    for m in range(1, 64):
        rows.append({m: 1}); b.append(h[m])
    for (X, Y), Z in ((P1, Z1), (P2, Z2)):
        rows += [{X | Z: 1}, {Y | Z: 1}, {Z: 1}]
        b += [h[X], h[Y], h[X] + h[Y] - h[X | Y]]
    Aeq = lil_matrix((len(rows), NV8))
    for k, d in enumerate(rows):
        for m, c in d.items():
            Aeq[k, m - 1] = c
    res = linprog(np.zeros(NV8), A_ub=A_UB, b_ub=np.zeros(A_UB.shape[0]), A_eq=csr_matrix(Aeq),
                  b_eq=np.array(b, dtype=float), bounds=[(None, None)] * NV8, method='highs', options=T.HIGHS_OPTS)
    return int(res.status)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('rays')
    ap.add_argument('--index', type=int, nargs='*', default=None)
    ap.add_argument('--maxsize', type=int, default=2)
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    R = np.load(a.rays).astype(np.int64)
    idx = list(range(R.shape[0])) if a.index is None else a.index
    out = []
    for i in idx:
        t0 = time.time(); h = T.h_norm(R[i]); P = small_pairs(h, a.maxsize)
        n, verdict, bad, trouble = 0, 'feasible', None, 0
        for P1, P2 in itertools.combinations(P, 2):
            st = status(h, P1, P2); n += 1
            if st == 2:
                verdict, bad = 'INFEASIBLE', [list(P1), list(P2)]
                break
            if st != 0:
                trouble += 1
        out.append({'index': int(i), 'verdict': verdict, 'pairs': len(P), 'lps': n, 'bad': bad, 'solver_troubles': trouble})
        print(f'ray {i}: {verdict}{" at " + str(bad) if bad else ""}; {len(P)} pairs, {n} pairs of pairs, '
              f'{trouble} solver troubles, {time.time() - t0:.0f}s', flush=True)
    if a.out:
        json.dump({'rays': a.rays, 'maxsize': a.maxsize, 'results': out}, open(a.out, 'w'), indent=1)


if __name__ == '__main__':
    main()
