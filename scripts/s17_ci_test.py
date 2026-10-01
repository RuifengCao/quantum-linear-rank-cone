#!/usr/bin/env python3
"""s17: common-information (Ingleton-type) test of five-party rays against the
linear representability that every stabilizer state must have.

Background (proved in the R27 report, cited: Gross-Walter arXiv:1302.6990):
for a stabilizer state in normal form (party p owns exactly S_p qubits, always
achievable without changing the entropy vector) the function
    h(X) = S(X) + sum_{p in X} S_p ,   X subset of the six parties A..E,F,
is the rank function of a GF(2) subspace arrangement.  Hence for every pair of
party sets (X, Y) the arrangement extends by Z = A_X cap A_Y: a 7-element
polymatroid with h(XZ) = h(X), h(YZ) = h(Y), h(Z) = I(X;Y).  If the LP
"extend h to 7 elements satisfying Shannon + these CI equalities" is infeasible,
its Farkas dual is a linear inequality valid for all stabilizer entropy vectors
and violated by the ray: the ray is NOT in the stabilizer cone.

  python3 scripts/s17_ci_test.py rays.npy [--out result.json]
Prints, per ray, the first infeasible pair (X,Y) or 'feasible for all pairs'.
Infeasibility found here is a floating-point LP verdict; the exact certificate
is produced and checked by s17_certify.py.
"""
import argparse, itertools, json, sys
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import lil_matrix, csr_matrix

N7 = 7
FULL6 = 63


def S6(r):
    """entropy of all 63 nonempty subsets of the six parties (purifier = bit 5) from a 31-vector."""
    S = np.zeros(64, dtype=np.int64)
    for m in range(1, 64):
        if m & 32:
            c = FULL6 ^ m
            S[m] = 0 if c == 0 else r[c - 1]
        else:
            S[m] = r[m - 1]
    return S


def h_norm(r):
    S = S6(r)
    sing = [S[1 << p] for p in range(6)]
    h = np.zeros(64, dtype=np.int64)
    for m in range(1, 64):
        h[m] = S[m] + sum(sing[p] for p in range(6) if m >> p & 1)
    return h


def elemental_rows(n=N7):
    """Shannon elemental inequalities as rows over masks 1..2^n-1 (row . x >= 0)."""
    full = (1 << n) - 1
    rows = []
    for i in range(n):
        rows.append({full: 1, full ^ (1 << i): -1} if (full ^ (1 << i)) else {full: 1})
    for i in range(n):
        for j in range(i + 1, n):
            rest = full ^ (1 << i) ^ (1 << j)
            K = rest
            while True:
                d = {}
                for m, c in ((K | 1 << i, 1), (K | 1 << j, 1), (K | 1 << i | 1 << j, -1), (K, -1)):
                    if m:
                        d[m] = d.get(m, 0) + c
                rows.append(d)
                if K == 0:
                    break
                K = (K - 1) & rest
    return rows


ELEM = elemental_rows()
NV = (1 << N7) - 1
A_ELEM = lil_matrix((len(ELEM), NV))
for k, d in enumerate(ELEM):
    for m, c in d.items():
        A_ELEM[k, m - 1] = -c          # -row.x <= 0
A_ELEM = csr_matrix(A_ELEM)
Z = 1 << 6


def ci_feasible(h, X, Y):
    eq_rows, eq_b = [], []
    for m in range(1, 64):
        eq_rows.append({m: 1}); eq_b.append(h[m])
    eq_rows.append({X | Z: 1}); eq_b.append(h[X])
    eq_rows.append({Y | Z: 1}); eq_b.append(h[Y])
    eq_rows.append({Z: 1}); eq_b.append(h[X] + h[Y] - h[X | Y])
    Aeq = lil_matrix((len(eq_rows), NV))
    for k, d in enumerate(eq_rows):
        for m, c in d.items():
            Aeq[k, m - 1] = c
    res = linprog(np.zeros(NV), A_ub=A_ELEM, b_ub=np.zeros(A_ELEM.shape[0]),
                  A_eq=csr_matrix(Aeq), b_eq=np.array(eq_b, dtype=float),
                  bounds=[(None, None)] * NV, method='highs')
    return res.status == 0


def pairs_for(h, mode='all'):
    """candidate CI pairs (X, Y) with I(X;Y) > 0 and neither set containing the other.
    mode 'disjoint': only X & Y == 0, ordered by the pair types that were infeasible
    most often on the R27 frontier ((1,2), (1,3), (2,2), (2,3), (1,4), (1,1)); every one
    of the 139 R27 exclusions has a disjoint infeasible pair.  'all' is the complete test."""
    out = []
    subs = range(1, 64)
    for X in subs:
        for Y in subs:
            if Y <= X or (X & Y) == X or (X & Y) == Y:
                continue
            if mode == 'disjoint' and (X & Y):
                continue
            if h[X] + h[Y] - h[X | Y] > 0:
                out.append((X, Y))
    if mode == 'disjoint':
        rank = {(1, 2): 0, (1, 3): 1, (2, 2): 2, (2, 3): 3, (1, 4): 4, (1, 1): 5}
        pc = lambda m: bin(m).count('1')
        out.sort(key=lambda p: rank.get(tuple(sorted((pc(p[0]), pc(p[1])))), 9))
    return out


_MODE, _FIRST = 'all', False


def _test_one(item):
    i, r = item
    h = h_norm(np.asarray(r, dtype=np.int64))
    bad = []
    for X, Y in pairs_for(h, _MODE):
        if not ci_feasible(h, X, Y):
            bad.append((X, Y))
            if _FIRST:
                break
    return i, bad


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('rays')
    ap.add_argument('--out', default=None)
    ap.add_argument('--limit', type=int, default=None)
    ap.add_argument('--first-only', action='store_true', help='stop at the first infeasible pair per ray')
    ap.add_argument('--pairs', default='all', choices=['all', 'disjoint'],
                    help="'all' = complete single-CI test; 'disjoint' = fast screen (~4x fewer LPs)")
    ap.add_argument('--workers', type=int, default=1)
    a = ap.parse_args()
    R = np.load(a.rays).astype(np.int64)
    if a.limit:
        R = R[:a.limit]
    global _MODE, _FIRST
    _MODE, _FIRST = a.pairs, a.first_only
    results = []
    items = list(enumerate(R.tolist()))
    if a.workers > 1:
        from multiprocessing import get_context
        with get_context('fork').Pool(a.workers) as pool:
            it = pool.imap(_test_one, items, chunksize=4)
            for i, bad in it:
                results.append({'index': i, 'infeasible_pairs': bad})
                print(f'ray {i}: ' + (f'INFEASIBLE for {len(bad)} pair(s), first {bad[0]}' if bad else f'feasible for all {a.pairs} pairs'), flush=True)
    else:
        for item in items:
            i, bad = _test_one(item)
            results.append({'index': i, 'infeasible_pairs': bad})
            print(f'ray {i}: ' + (f'INFEASIBLE for {len(bad)} pair(s), first {bad[0]}' if bad else f'feasible for all {a.pairs} pairs'), flush=True)
    if a.out:
        json.dump(results, open(a.out, 'w'))


if __name__ == '__main__':
    main()
