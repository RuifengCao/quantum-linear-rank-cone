#!/usr/bin/env python3
"""s42: two common informations plus conditional Ingleton on the 8-element extension, with nesting (Gate 1, E2).

For a ray r of QLR_5, h = h_norm(r) (R32, Theorem 1; s17_ci_test.h_norm) is the rank function of a subspace
arrangement W_A..W_F over a prime field whenever r is a stabilizer entropy vector.  For two pairs (X1, Y1) and
(X2, Y2) of disjoint subsets, the arrangement extends by Z1 = W_X1 cap W_Y1 and then by Z2 = W_X2 cap W_Y2, where
X2 or Y2 may contain Z1 (nesting; the extended arrangement is again linear).  The 8-element arrangement
A..F, Z1, Z2 then satisfies every Shannon inequality, every conditional Ingleton instance, and the
common-information equalities

    g(Xi Zi) = g(Xi),  g(Yi Zi) = g(Yi),  g(Zi) = g(Xi) + g(Yi) - g(Xi Yi)      (i = 1, 2),

with g = h on the subsets of A..F.  The Ingleton instances are set-valued: A, B, C, D are disjoint nonempty unions
of elements and K a union of the remaining ones,

    g(AK) + g(BK) + g(CDK) + g(ABCK) + g(ABDK) <= g(ABK) + g(ACK) + g(ADK) + g(BCK) + g(BDK).

Ingleton holds for any four subspaces (Ingleton 1971), and contracting by W_K keeps an arrangement linear, so every
instance holds for every linear arrangement over every field.  R29 (s17 --extension ingleton) used the instances
with A, B, C, D single elements; these are among them.  Only the instances that involve Z1 or Z2 vary (119,616 of
them).  If the linear program "find such a g" is infeasible, its Farkas dual is a linear inequality that every
stabilizer state of every prime dimension satisfies and r violates.  R30 is the case of one common information
with single-element Ingleton; s25 the case of two without Ingleton.

Verdicts are floating point (HiGHS).  A row reported infeasible is a candidate exclusion only: it counts after an
exact certificate, which this script does not produce.  A row stops at its first infeasible pair of pairs.

The Ingleton rows enter by cutting planes: the LP starts with the Shannon rows and adds the Ingleton rows its
solution violates (at most 400 per round) until none is violated (feasible) or the LP becomes infeasible
(infeasible with all rows, since it is infeasible with a subset); the rows added for one pair of pairs are kept for
the next ones of the same row (at most 2,000, then the set starts afresh).  --self-test checks every row on random
subspace arrangements over GF(2) and GF(3), known answers (rows excluded in R30, known extreme rays), and
cutting planes against the LP with all rows at once.

  python3 scripts/s42_ci2_ingleton.py RAYS.npy [--rows 3,7] [--maxsize 2] [--nested] [--workers W] [--hours H]
                                        --out OUT.json [--resume]
  python3 scripts/s42_ci2_ingleton.py --self-test

Pairs of pairs: the first pair (X1, Y1) runs over disjoint nonempty subsets of A..F with at most --maxsize elements
each and I(X1; Y1) > 0; the second pair runs over the same list (unordered pairs of distinct pairs, as in s25)
and, with --nested, over every pair of disjoint subsets of A..F, Z1 with at most --maxsize elements each that
involves Z1.  Rows are processed in parallel, one row per task, largest N first, and the output is rewritten after
every row; --resume skips the rows already in it.  With --hours H no row is started and no LP is begun after H
hours; rows cut off then are left out of the output (unfinished) and a later --resume does them again.
"""
import argparse, itertools, json, multiprocessing as mp, os, random, sys, time, warnings

for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):   # one thread per worker process
    os.environ.setdefault(_v, '1')
import numpy as np
from scipy.optimize import linprog, OptimizeWarning
from scipy.sparse import coo_matrix, csr_matrix, lil_matrix, vstack

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s17_ci_test as T

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
warnings.filterwarnings('ignore', message='Unrecognized options detected|Unknown solver options', category=OptimizeWarning)
N8, NV8, Z1, Z2 = 8, 255, 1 << 6, 1 << 7
ZZ = Z1 | Z2


def ingleton_matrix8():
    """the set-valued conditional Ingleton instances on 8 elements that involve Z1 or Z2, as rows of A with
    A g <= 0 (duplicates removed)."""
    data, ri, ci, seen, k = [], [], [], set(), 0
    for assign in itertools.product(range(6), repeat=N8):    # element e goes to A, B, C, D, K or nowhere
        m = [0] * 6
        for e, t in enumerate(assign):
            m[t] |= 1 << e
        A, B, C, D, K = m[:5]
        if not (A and B and C and D) or A > B or C > D or not ((A | B | C | D | K) & ZZ):
            continue
        row = {}
        for mm, coef in ((A | B, 1), (A | C, 1), (A | D, 1), (B | C, 1), (B | D, 1),
                         (A, -1), (B, -1), (C | D, -1), (A | B | C, -1), (A | B | D, -1)):
            row[mm | K] = row.get(mm | K, 0) + coef
        key = tuple(sorted((x, c) for x, c in row.items() if c))
        if not key or key in seen:
            continue
        seen.add(key)
        for x, c in key:
            ri.append(k); ci.append(x - 1); data.append(-c)
        k += 1
    return csr_matrix(coo_matrix((data, (ri, ci)), shape=(k, NV8)))


def _matrix(rows):
    A = lil_matrix((len(rows), NV8))
    for k, d in enumerate(rows):
        for m, c in d.items():
            A[k, m - 1] = -c                    # A x <= 0  <=>  row . x >= 0
    return csr_matrix(A)


A_SH = A_ING = A_ALL = None


def setup():
    """build the inequality matrices once (before forking workers)."""
    global A_SH, A_ING, A_ALL
    if A_SH is None:
        A_SH = _matrix([r for r in T.elemental_rows(N8) if any(m & ZZ for m in r)])
        A_ING = ingleton_matrix8()
        A_ALL = vstack([A_SH, A_ING]).tocsr()


def _equalities(h, P1, P2):
    rows, b = [], []
    for m in range(1, 64):
        rows.append({m: 1}); b.append(float(h[m]))
    for (X, Y), Z in ((P1, Z1), (P2, Z2)):
        rows.append({X | Z: 1, X: -1}); b.append(0.0)
        rows.append({Y | Z: 1, Y: -1}); b.append(0.0)
        d = {Z: 1}
        for m, c in ((X, -1), (Y, -1), (X | Y, 1)):
            d[m] = d.get(m, 0) + c
        rows.append(d); b.append(0.0)
    Aeq = lil_matrix((len(rows), NV8))
    for k, d in enumerate(rows):
        for m, c in d.items():
            if c:
                Aeq[k, m - 1] = c
    return csr_matrix(Aeq), np.array(b)


def _lp(A, Aeq, b):
    return linprog(np.zeros(NV8), A_ub=A, b_ub=np.zeros(A.shape[0]), A_eq=Aeq, b_eq=b,
                   bounds=[(None, None)] * NV8, method='highs', options=T.HIGHS_OPTS)


def status(h, P1, P2, cut=True, active=None):
    """HiGHS status of the LP: 0 feasible, 2 infeasible, anything else (9: cutting planes did not settle in 100
    rounds) is solver trouble and no verdict.  active: boolean mask of Ingleton rows to start from, updated in place
    (run_row carries it from one pair of pairs to the next; any set of valid rows gives the same verdict)."""
    Aeq, b = _equalities(h, P1, P2)
    if not cut:
        return int(_lp(A_ALL, Aeq, b).status)
    if active is None:
        active = np.zeros(A_ING.shape[0], dtype=bool)
    for _ in range(100):
        A = vstack([A_SH, A_ING[np.where(active)[0]]]).tocsr() if active.any() else A_SH
        res = _lp(A, Aeq, b)
        if res.status != 0:
            return int(res.status)
        viol = A_ING @ res.x                    # a violated row has A_ING x > 0
        new = np.where((viol > 1e-7) & ~active)[0]
        if len(new) == 0:
            return 0
        active[new[np.argsort(-viol[new])][:400]] = True
    return 9


def first_pairs(h, maxsize):
    subs = [m for m in range(1, 64) if bin(m).count('1') <= maxsize]
    return [(X, Y) for X in subs for Y in subs if Y > X and not X & Y and h[X] + h[Y] - h[X | Y] > 0]


def nested_pairs(maxsize):
    """pairs of disjoint subsets of A..F, Z1 with at most maxsize elements each that involve Z1."""
    subs = [m for m in range(1, 128) if bin(m).count('1') <= maxsize]
    return [(X, Y) for X in subs for Y in subs if Y > X and not X & Y and (X | Y) & Z1]


def combos(h, maxsize, nested):
    P = first_pairs(h, maxsize)
    out = list(itertools.combinations(P, 2))
    if nested:
        Q = nested_pairs(maxsize)
        out += [(p, q) for p in P for q in Q]
    return out


def run_row(task):
    k, r, maxsize, nested, deadline = task
    h = T.h_norm(np.asarray(r, dtype=np.int64))
    C = combos(h, maxsize, nested)
    t0 = time.time(); trouble = 0; done = 0; bad = None
    active = np.zeros(A_ING.shape[0], dtype=bool)
    for P1, P2 in C:
        if deadline and time.time() > deadline:
            return dict(row=k, unfinished=True, done=done, combos=len(C))
        if active.sum() > 2000:                 # keep the LPs small
            active[:] = False
        st = status(h, P1, P2, active=active); done += 1
        if st == 2:
            bad = [list(P1), list(P2)]
            break
        if st != 0:
            trouble += 1
    return dict(row=k, N=int(h[63]), combos=len(C), done=done, infeasible=bad, trouble=trouble,
                seconds=round(time.time() - t0, 1))


# ------------------------------------------------------------------------------------------------- self-test
def _rank_mod(M, p):
    M = M.copy() % p
    r = 0
    for c in range(M.shape[1]):
        piv = next((i for i in range(r, M.shape[0]) if M[i, c]), None)
        if piv is None:
            continue
        M[[r, piv]] = M[[piv, r]]
        M[r] = (M[r] * pow(int(M[r, c]), p - 2, p)) % p
        for i in range(M.shape[0]):
            if i != r and M[i, c]:
                M[i] = (M[i] - M[i, c] * M[r]) % p
        r += 1
        if r == M.shape[0]:
            break
    return r


def self_test(log=print):
    """(1) every Shannon and Ingleton row holds on 20 random subspace arrangements of 8 elements over GF(2) and
    GF(3); (2) three rows excluded in R30 by CI + Ingleton at (X, Y) stay infeasible with (X, Y) first and a second
    pair added (plain and nested), and two of the 641 extreme rays of data/ are feasible on 20 sampled pairs of pairs
    each (maxsize 2, nested included); (3) on four of these LPs, cutting planes and the LP with all rows at once
    agree."""
    setup()
    rng, nrng = random.Random(42), np.random.default_rng(42)
    worst = 0.0
    for t in range(20):
        p, dim = (2, 3)[t % 2], int(nrng.integers(4, 9))
        gens = [nrng.integers(0, p, size=(int(nrng.integers(1, 4)), dim)) for _ in range(N8)]
        g = np.array([_rank_mod(np.vstack([gens[e] for e in range(N8) if m >> e & 1]), p) for m in range(1, 256)],
                     dtype=float)
        worst = max(worst, float((A_ALL @ g).max()))
    ok = worst <= 0
    d = json.load(open(os.path.join(ROOT, 'data', 's22_ingleton_exclusions_a2i.json')))
    known, n = [], 0
    for c in rng.sample(d['certificates'], 3):
        h = T.h_norm(np.array(c['ray'], dtype=np.int64))
        for P2 in (rng.choice([p for p in first_pairs(h, 2) if p != (c['X'], c['Y'])]), rng.choice(nested_pairs(2))):
            known.append((h, (c['X'], c['Y']), P2, 2))
    K = np.load(os.path.join(ROOT, 'data', 'stab5_extreme_reps.npy')).astype(np.int64)
    for i in rng.sample(range(len(K)), 2):
        h = T.h_norm(K[i])
        known += [(h, P1, P2, 0) for P1, P2 in rng.sample(combos(h, 2, True), 20)]
    for h, P1, P2, want in known:
        ok &= status(h, P1, P2) == want
    for h, P1, P2, want in known[:2] + known[-2:]:
        ok &= status(h, P1, P2, cut=False) == want
    log(f's42 self-test: rows hold on 20 random arrangements over GF(2) and GF(3) (max of A g = {worst:g}); '
        f'{len(known)} LPs with known answers (6 infeasible, 40 feasible), 4 of them also with all rows at once: {ok}')
    return ok


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('rays', nargs='?')
    ap.add_argument('--rows', default=None, help='comma-separated rows of the rays file (default: all)')
    ap.add_argument('--maxsize', type=int, default=2)
    ap.add_argument('--nested', action='store_true')
    ap.add_argument('--workers', type=int, default=None, help='default: all cores')
    ap.add_argument('--hours', type=float, default=None, help='start no row and no LP after this many hours')
    ap.add_argument('--out', default=None)
    ap.add_argument('--resume', action='store_true')
    ap.add_argument('--self-test', action='store_true')
    a = ap.parse_args()
    if a.self_test:
        sys.exit(0 if self_test() else 1)
    if not a.rays or not a.out:
        ap.error('rays file and --out are required')
    t_start = time.time()
    deadline = t_start + 3600 * a.hours if a.hours else None
    if not self_test():
        sys.exit('self-test failed; nothing run')
    R = np.load(a.rays).astype(np.int64)
    rows = [int(x) for x in a.rows.split(',')] if a.rows else list(range(len(R)))
    res = {}
    if a.resume and os.path.exists(a.out):
        res = json.load(open(a.out)).get('results', {})
        print(f'resume: {len(res)} rows already in {a.out}', flush=True)
    todo = [k for k in rows if str(k) not in res]
    todo.sort(key=lambda k: -int(T.h_norm(R[k])[63]))          # largest N first
    nw = a.workers or os.cpu_count() or 1
    print(f'{len(todo)} rows to test ({len(rows)} requested), maxsize {a.maxsize}, nested {a.nested}, '
          f'{nw} workers, {A_SH.shape[0]} Shannon + {A_ING.shape[0]} set-valued Ingleton rows'
          + (f', deadline {a.hours} h' if a.hours else ''), flush=True)
    meta = dict(rays=a.rays, rows=len(rows), maxsize=a.maxsize, nested=a.nested, ingleton='set-valued',
                cutting_planes=True, note='rows refer to the input rays file; verdicts are floating point; '
                'an infeasible row is a candidate exclusion until certified exactly')

    def save(final=False):
        inf = sorted(int(k) for k, v in res.items() if v['infeasible'])
        summary = dict(tested=len(res), infeasible_rows=inf, trouble_lps=sum(v['trouble'] for v in res.values()),
                       lps=sum(v['done'] for v in res.values()), complete=final and len(res) == len(rows))
        json.dump(dict(meta, summary=summary, results=res), open(a.out, 'w'), indent=1)
        return summary

    unfinished = 0
    with mp.get_context('fork').Pool(nw) as pool:
        it = pool.imap_unordered(run_row, [(k, R[k].tolist(), a.maxsize, a.nested, deadline) for k in todo])
        for out in it:
            if out.get('unfinished'):
                unfinished += 1
                continue
            res[str(out['row'])] = {x: out[x] for x in ('N', 'combos', 'done', 'infeasible', 'trouble', 'seconds')}
            v = 'INFEASIBLE at ' + str(out['infeasible']) if out['infeasible'] else 'feasible'
            print(f"row {out['row']} (N = {out['N']}): {v}; {out['done']}/{out['combos']} LPs, {out['trouble']} solver "
                  f"troubles, {out['seconds']} s; {len(res)} of {len(rows)} rows done ({time.time() - t_start:.0f} s)",
                  flush=True)
            save()
    s = save(final=True)
    print(f"tested {s['tested']} of {len(rows)} rows ({unfinished} cut off by the deadline), {s['lps']} LPs, "
          f"{s['trouble_lps']} with solver trouble; infeasible (candidate exclusions): {len(s['infeasible_rows'])} "
          f"{s['infeasible_rows']} ({time.time() - t_start:.0f} s)", flush=True)
    print('COMPLETE' if s['complete'] else 'INCOMPLETE: run again with --resume', flush=True)


if __name__ == '__main__':
    setup()
    main()
