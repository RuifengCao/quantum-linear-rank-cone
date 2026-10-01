#!/usr/bin/env python3
"""s23: common information + Ingleton + DFZ five-variable linear rank inequalities on the extension (R30).

The R29 test (s17_ci_test.py --extension ingleton) extends h_norm(r) by Z = A_X cap A_Y and asks
for a 7-element polymatroid satisfying Shannon and every conditional Ingleton instance that
involves Z.  This script adds the 24 five-variable linear rank inequalities of Dougherty, Freiling
and Zeger (arXiv:0910.0284, eqs. (1)-(24)) and the four Ingleton forms (36)-(39), as transcribed in
data/dfz_ref28.csv, instantiated on the 7-element ground set A..F, Z: the five variables go to five
distinct elements, and the instance is contracted by any subset K of the two remaining elements,

    sum_m c_m * ( g(phi(m) | K) - g(K) ) >= 0 .

Every instance is valid for every linear polymatroid over any field (DFZ; contraction by K keeps
the arrangement linear), hence for the extension of a stabilizer state of any prime dimension.
149,520 distinct instances; the 128,040 that involve Z enter the LP as cutting planes: solve with
Shannon + Ingleton, add the DFZ rows the solution violates (at most 400 per round), repeat.

  python3 scripts/s23_ci_dfz.py rays.npy [--index 0 3 7] [--out result.json]

Verdicts are floating-point LP verdicts.  A 'feasible for all pairs' verdict means the final LP
solution of every pair satisfies all 128,040 rows (to 1e-6).  No exact certificate path exists yet
for an infeasible pair (none has been found: R30 ran the 10 smallest undecided orbits).
"""
import argparse, csv, itertools, json, os, sys, time, warnings
import numpy as np
from scipy.optimize import linprog, OptimizeWarning
from scipy.sparse import lil_matrix, csr_matrix, vstack

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import s17_ci_test as T

warnings.filterwarnings('ignore', message='Unrecognized options detected|Unknown solver options', category=OptimizeWarning)
DFZ_CSV = os.path.join(HERE, '..', 'data', 'dfz_ref28.csv')


def load_dfz(path=DFZ_CSV):
    """{ineq_id: {mask over ABCDE: coefficient}}, row . h >= 0 for every linear rank function h."""
    let = {c: i for i, c in enumerate('ABCDE')}
    ineqs = {}
    lines = (ln for ln in open(path) if not ln.lstrip().startswith('#'))
    for row in csv.DictReader(lines):
        m = 0
        for ch in row['subset'].strip():
            m |= 1 << let[ch.upper()]
        d = ineqs.setdefault(row['ineq_id'], {})
        d[m] = d.get(m, 0) + int(row['coeff'])
    return {k: {m: v for m, v in d.items() if v} for k, d in ineqs.items()}


def instance(c, phi, K):
    """the instance of the five-variable form c with variable j -> element phi[j], contracted by K."""
    row = {}
    for m, v in c.items():
        M = K
        for j in range(5):
            if m >> j & 1:
                M |= 1 << phi[j]
        row[M] = row.get(M, 0) + v
        if K:
            row[K] = row.get(K, 0) - v
    return {M: v for M, v in row.items() if v}


def dfz_instances(n=7, ineqs=None):
    """distinct instances on n elements: list of (row dict over masks, descriptor (ineq_id, phi, K))."""
    ineqs = load_dfz() if ineqs is None else ineqs
    seen = {}
    for iid, c in ineqs.items():
        for phi in itertools.permutations(range(n), 5):
            rest = [e for e in range(n) if e not in phi]
            for kk in range(1 << len(rest)):
                K = sum(1 << rest[i] for i in range(len(rest)) if kk >> i & 1)
                row = instance(c, phi, K)
                key = tuple(sorted(row.items()))
                if key not in seen:
                    seen[key] = (iid, phi, K)
    return [(dict(k), d) for k, d in seen.items()]


def z_rows_matrix(z=6):
    rows = [r for r, _ in dfz_instances() if any(M >> z & 1 for M in r)]
    D = lil_matrix((len(rows), T.NV))
    for k, r in enumerate(rows):
        for M, v in r.items():
            D[k, M - 1] = v
    return csr_matrix(D)


def eq_system(h, X, Y):
    rows, b = [], []
    for m in range(1, 64):
        rows.append({m: 1}); b.append(h[m])
    rows += [{X | T.Z: 1}, {Y | T.Z: 1}, {T.Z: 1}]
    b += [h[X], h[Y], h[X] + h[Y] - h[X | Y]]
    A = lil_matrix((len(rows), T.NV))
    for k, r in enumerate(rows):
        for M, v in r.items():
            A[k, M - 1] = v
    return csr_matrix(A), np.array(b, dtype=float)


def test_pair(h, X, Y, D, active, cap=400, maxit=60, tol=1e-6):
    """('feasible' | 'infeasible' | 'solver<status>' | 'maxit', rounds); active: set of D rows, grown in place."""
    Aeq, beq = eq_system(h, X, Y)
    for it in range(maxit):
        Aub = vstack([T.A_UB, -D[sorted(active)]]) if active else T.A_UB
        res = linprog(np.zeros(T.NV), A_ub=Aub, b_ub=np.zeros(Aub.shape[0]), A_eq=Aeq, b_eq=beq,
                      bounds=[(None, None)] * T.NV, method='highs', options=T.HIGHS_OPTS)
        if res.status == 2:
            return 'infeasible', it
        if res.status != 0:
            return f'solver{res.status}', it
        val = D @ res.x
        viol = np.flatnonzero(val < -tol)
        if viol.size == 0:
            return 'feasible', it
        active.update(viol[np.argsort(val[viol])[:cap]].tolist())
    return 'maxit', maxit


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('rays')
    ap.add_argument('--index', type=int, nargs='*', default=None)
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    T.set_extension('ingleton')
    t0 = time.time()
    D = z_rows_matrix()
    print(f'{D.shape[0]} DFZ instances involving Z ({time.time() - t0:.0f}s)', flush=True)
    R = np.load(a.rays).astype(np.int64)
    idx = list(range(R.shape[0])) if a.index is None else a.index
    out = []
    for i in idx:
        t1 = time.time(); h = T.h_norm(R[i]); active = set(); rounds = 0; trouble = 0
        verdict, bad = 'feasible for all pairs', None
        pairs = T.pairs_for(h, 'all')
        for X, Y in pairs:
            st, it = test_pair(h, X, Y, D, active)
            rounds += it
            if st == 'infeasible':
                verdict, bad = 'INFEASIBLE', (X, Y)
                break
            if st != 'feasible':
                trouble += 1
        if trouble and verdict != 'INFEASIBLE':
            verdict = 'UNDETERMINED'
        out.append({'index': int(i), 'verdict': verdict, 'pair': bad, 'pairs': len(pairs), 'cut_rounds': rounds,
                    'dfz_rows_activated': len(active), 'solver_troubles': trouble})
        print(f'ray {i}: {verdict}{" at " + str(bad) if bad else ""}; {len(pairs)} pairs, {rounds} cutting rounds, '
              f'{len(active)} DFZ rows activated, {trouble} solver troubles, {time.time() - t1:.0f}s', flush=True)
    if a.out:
        json.dump({'rays': a.rays, 'dfz_rows_with_Z': int(D.shape[0]), 'results': out}, open(a.out, 'w'), indent=1)


if __name__ == '__main__':
    main()
