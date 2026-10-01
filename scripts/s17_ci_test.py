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
is produced and checked by s17_certify.py.  An LP that ends with any HiGHS
status other than 0 (feasible) or 2 (infeasible) is reported as a solver error
and never counted as a verdict.

Threads (R28.1): HiGHS, scipy's LP solver, starts a pool of about nproc/2
threads in every process that solves an LP.  With one LP worker per core that
is ~nproc^2/2 threads on a large server (about 21,000 on 208 cores), so every
LP here runs with options HIGHS_OPTS = {'threads': 1}.  scipy >= 1.11 forwards
the option to HiGHS; older scipy drops it, and main() then caps the number of
workers.  All LPs in one process must use the same setting: once HiGHS has
started its pool, a call asking for a different thread count fails (scipy
1.17: status 4).  s17_certify.py and epr1kit/core.py therefore use the same
options.
"""
import argparse, itertools, json, os, sys, warnings
import numpy as np
import scipy
from scipy.optimize import linprog, OptimizeWarning
from scipy.sparse import lil_matrix, csr_matrix

HIGHS_OPTS = {'threads': 1}
# scipy >= 1.11 calls the option 'unrecognized' and passes it to HiGHS verbatim; older
# scipy drops it ('Unknown solver options').  Either way the warning would be printed
# once per worker process.
warnings.filterwarnings('ignore', message='Unrecognized options detected|Unknown solver options', category=OptimizeWarning)


def _scipy_pins_threads():
    """True when this scipy forwards HIGHS_OPTS to HiGHS (measured: 1.9.3, 1.10.1 do not;
    1.11.4, 1.14.1, 1.17.1 do)."""
    try:
        return tuple(int(x) for x in scipy.__version__.split('.')[:2]) >= (1, 11)
    except ValueError:
        return True


SCIPY_PINS_THREADS = _scipy_pins_threads()

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


def ingleton_rows(n=N7, z=6, with_desc=False):
    """(R29) conditional Ingleton instances on the n-element extension that involve element z:
    A, B, C, D distinct singletons, K any subset of the remaining elements, z in {A, B, C, D} or in K;
    rows over masks (row . x >= 0):
        g(ABK) + g(ACK) + g(ADK) + g(BCK) + g(BDK) - g(AK) - g(BK) - g(CDK) - g(ABCK) - g(ABDK) >= 0.
    Valid for every linear polymatroid over any field (Ingleton 1971; contracting by K keeps the
    arrangement linear), hence for the extension A..F, Z = A_X cap A_Y of a stabilizer state."""
    rows, desc = [], []
    for quad in itertools.combinations(range(n), 4):
        rest = [e for e in range(n) if e not in quad]
        w, x, y, v = quad
        for (a, b), (c, d) in (((w, x), (y, v)), ((y, v), (w, x)), ((w, y), (x, v)), ((x, v), (w, y)),
                               ((w, v), (x, y)), ((x, y), (w, v))):
            for kk in range(1 << len(rest)):
                K = sum(1 << rest[i] for i in range(len(rest)) if kk >> i & 1)
                if z not in quad and not (K >> z & 1):
                    continue
                A, B, C, D = 1 << a, 1 << b, 1 << c, 1 << d
                row = {}
                for m, coef in ((A | B, 1), (A | C, 1), (A | D, 1), (B | C, 1), (B | D, 1),
                                (A, -1), (B, -1), (C | D, -1), (A | B | C, -1), (A | B | D, -1)):
                    row[m | K] = row.get(m | K, 0) + coef
                rows.append(row); desc.append(('I', a, b, c, d, K))
    return (rows, desc) if with_desc else rows


def _rows_matrix(rows):
    A = lil_matrix((len(rows), NV))
    for k, d in enumerate(rows):
        for m, c in d.items():
            A[k, m - 1] = -c          # -row.x <= 0
    return csr_matrix(A)


ELEM = elemental_rows()
NV = (1 << N7) - 1
A_ELEM = _rows_matrix(ELEM)
A_UB = A_ELEM                    # inequality rows used by ci_status (see set_extension)
EXTENSION = 'shannon'
Z = 1 << 6


def set_extension(kind):
    """'shannon' (default, the R27 test): Shannon inequalities on the 7-element extension;
    'ingleton' (R29): Shannon + every conditional Ingleton instance involving Z.  Call before forking workers."""
    global A_UB, EXTENSION
    if kind == 'shannon':
        A_UB = A_ELEM
    elif kind == 'ingleton':
        A_UB = _rows_matrix(ELEM + ingleton_rows())
    else:
        raise ValueError(kind)
    EXTENSION = kind


def ci_status(h, X, Y):
    """HiGHS status of the CI-extension LP for the pair (X, Y): 0 = feasible,
    2 = infeasible; any other value is solver trouble, not a verdict."""
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
    res = linprog(np.zeros(NV), A_ub=A_UB, b_ub=np.zeros(A_UB.shape[0]),
                  A_eq=csr_matrix(Aeq), b_eq=np.array(eq_b, dtype=float),
                  bounds=[(None, None)] * NV, method='highs', options=HIGHS_OPTS)
    return int(res.status)


def ci_feasible(h, X, Y):
    return ci_status(h, X, Y) == 0


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
MAX_ERR_PER_RAY = 20


def _test_one(item):
    """(index, infeasible pairs, solver errors as (X, Y, status)) for one ray."""
    i, r = item
    h = h_norm(np.asarray(r, dtype=np.int64))
    bad, err = [], []
    for X, Y in pairs_for(h, _MODE):
        st = ci_status(h, X, Y)
        if st == 2:
            bad.append((X, Y))
            if _FIRST:
                break
        elif st != 0:
            err.append((X, Y, st))
            if len(err) >= MAX_ERR_PER_RAY:
                break
    return i, bad, err


def _test_chunk(chunk):
    return [_test_one(item) for item in chunk]


def _line(i, bad, err, mode):
    if bad:
        s = f'INFEASIBLE for {len(bad)} pair(s), first {bad[0]}'
        if err:
            s += f' (also {len(err)} LP solver error(s))'
    elif err:
        s = f'UNDETERMINED: LP solver error for {len(err)} pair(s), first (X, Y, status) = {err[0]}'
    else:
        s = f'feasible for all {mode} pairs'
    return f'ray {i}: ' + s


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('rays')
    ap.add_argument('--out', default=None)
    ap.add_argument('--limit', type=int, default=None)
    ap.add_argument('--first-only', action='store_true', help='stop at the first infeasible pair per ray')
    ap.add_argument('--pairs', default='all', choices=['all', 'disjoint'],
                    help="'all' = complete single-CI test; 'disjoint' = fast screen (~4x fewer LPs)")
    ap.add_argument('--workers', type=int, default=1)
    ap.add_argument('--extension', default='shannon', choices=['shannon', 'ingleton'],
                    help="inequalities imposed on the 7-element extension: 'shannon' (R27) or "
                         "'ingleton' (R29: + conditional Ingleton instances involving Z)")
    ap.add_argument('--result-timeout', type=float, default=None,
                    help='abort if no result arrives for this many seconds (a dead worker would otherwise make '
                         'the pool wait forever); default 900 (shannon) / 1800 (ingleton, ~4x slower LPs)')
    ap.add_argument('--chunk', type=int, default=None, help='rays per pool task; default 4 (shannon) / 1 (ingleton)')
    a = ap.parse_args()
    R = np.load(a.rays).astype(np.int64)
    if a.limit:
        R = R[:a.limit]
    global _MODE, _FIRST
    _MODE, _FIRST = a.pairs, a.first_only
    set_extension(a.extension)
    if a.result_timeout is None:
        a.result_timeout = 900.0 if a.extension == 'shannon' else 1800.0
    if a.chunk is None:
        a.chunk = 4 if a.extension == 'shannon' else 1
    workers = a.workers
    if workers > 1 and not SCIPY_PINS_THREADS:
        per = max(1, ((os.cpu_count() or 1) + 1) // 2)
        cap = max(1, 8192 // (per + 2))
        if workers > cap:
            print(f'NOTE: scipy {scipy.__version__} cannot pin HiGHS to one thread (needs scipy >= 1.11), '
                  f'so every worker starts ~{per} threads; using {cap} workers instead of {workers}. '
                  f'To use all cores: python3 -m pip install -U "scipy>=1.11"', flush=True)
            workers = cap
    results = []
    count = {'err': 0}

    def record(i, bad, err):
        results.append({'index': i, 'infeasible_pairs': bad, 'lp_errors': err})
        print(_line(i, bad, err, a.pairs), flush=True)
        count['err'] += bool(err)
        if len(results) == 16 and count['err'] >= 8:
            raise SystemExit('ABORT: LP solver errors on most of the first 16 rays -- a systemic HiGHS '
                             'problem, not a property of the rays; send this log')

    items = list(enumerate(R.tolist()))
    if workers > 1:
        from multiprocessing import get_context, TimeoutError as PoolTimeout
        chunks = [items[k:k + a.chunk] for k in range(0, len(items), a.chunk)]
        with get_context('fork').Pool(workers) as pool:
            # chunks of a.chunk rays, one chunk per task: imap with chunksize=1 returns an iterator
            # whose next() takes a timeout (with chunksize > 1 it is a plain generator)
            it = pool.imap(_test_chunk, chunks)
            for _ in range(len(chunks)):
                try:
                    out = it.next(timeout=a.result_timeout)
                except PoolTimeout:
                    raise SystemExit(f'ABORT: no result for {a.result_timeout:.0f} s -- a worker process '
                                     'probably died (see messages above); not waiting forever')
                for i, bad, err in out:
                    record(i, bad, err)
    else:
        for item in items:
            record(*_test_one(item))
    n_inf = sum(1 for x in results if x['infeasible_pairs'])
    n_und = sum(1 for x in results if x['lp_errors'] and not x['infeasible_pairs'])
    print(f'summary: {len(results)} rays: {n_inf} infeasible, {len(results) - n_inf - n_und} feasible for all '
          f'{a.pairs} pairs, {n_und} undetermined (LP solver errors); extension {EXTENSION} '
          f'({A_UB.shape[0]} inequality rows); scipy {scipy.__version__}, '
          f'HiGHS threads per worker: {"1" if SCIPY_PINS_THREADS else "not pinnable"}', flush=True)
    if a.out:
        json.dump(results, open(a.out, 'w'))


if __name__ == '__main__':
    main()
