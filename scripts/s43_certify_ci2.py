#!/usr/bin/env python3
"""s43: exact certificates for the rows that s42_ci2_ingleton.py reports infeasible (two common informations).

For a ray r and the pair of pairs P1 = (X1, Y1), P2 = (X2, Y2) recorded by s42, find a Farkas certificate of the
infeasibility in exact rational arithmetic: nonnegative multipliers on Shannon elemental rows and set-valued
conditional Ingleton rows of the 8-element extension A..F, Z1, Z2, and free multipliers on the six
common-information equalities, whose combination cancels on every mask that meets Z1 or Z2 and leaves a form F on
A..F with F(h_norm(r)) < 0.  F(h_norm(S)) >= 0 then holds for every stabilizer state of every prime dimension
(the extension of a subspace arrangement by the two intersections is a subspace arrangement), so r is not in the
stabilizer cone.  Every certificate is checked, before it is written and again by --check, with
epr1kit.stabcert.verify_exclusion_ci2, which rebuilds each row from its descriptor without the row-building code of
s42 and of this script.

Method: the cutting-plane LP of s42 (Shannon rows, then the violated Ingleton rows) is run until it is infeasible;
the Farkas dual is solved over the rows it used (all rows if that fails), and its support is solved exactly.

  python3 scripts/s43_certify_ci2.py RAYS.npy --ci2 ci2_ingleton.json --out certs.json [--workers W]
  python3 scripts/s43_certify_ci2.py RAYS.npy --check certs.json
  python3 scripts/s43_certify_ci2.py --self-test

Descriptors: ('H', i) g(all) - g(all - i); ('M', i, j, K) I(i; j | K); ('I', A, B, C, D, K) the Ingleton form
I(A;B|CK) + I(A;B|DK) + I(C;D|K) - I(A;B|K) of four disjoint nonempty unions A, B, C, D given K; ('C', k) the k-th
equality g(X1 Z1) - g(X1), g(Y1 Z1) - g(Y1), g(Z1) - g(X1) - g(Y1) + g(X1 Y1), then the same for P2 and Z2.
"""
import argparse, itertools, json, math, multiprocessing as mp, os, random, sys, time, warnings

for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
from fractions import Fraction
import numpy as np
from scipy.optimize import linprog, OptimizeWarning
from scipy.sparse import coo_matrix, csr_matrix, vstack

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, HERE); sys.path.insert(0, ROOT)
from s17_certify import exact_solve
from epr1kit import stabcert

warnings.filterwarnings('ignore', message='Unrecognized options detected|Unknown solver options', category=OptimizeWarning)
N8, NV8, Z1, Z2 = 8, 255, 64, 128
ZZ = Z1 | Z2
FULL = 255
OPTS = {'threads': 1}


def row_of(d):
    t = d[0]
    if t == 'H':
        return {FULL: 1, FULL ^ (1 << d[1]): -1}
    if t == 'M':
        i, j, K = d[1:]
        row = {}
        for m, c in (((1 << i) | K, 1), ((1 << j) | K, 1), ((1 << i) | (1 << j) | K, -1), (K, -1)):
            if m:
                row[m] = row.get(m, 0) + c
        return row
    A, B, C, D, K = d[1:]
    row = {}
    for m, c in ((A | B, 1), (A | C, 1), (A | D, 1), (B | C, 1), (B | D, 1),
                 (A, -1), (B, -1), (C | D, -1), (A | B | C, -1), (A | B | D, -1)):
        row[m | K] = row.get(m | K, 0) + c
    return {m: c for m, c in row.items() if c}


def ci_rows(P1, P2):
    out = []
    for (X, Y), Z in ((P1, Z1), (P2, Z2)):
        out.append({X | Z: 1, X: -1})
        out.append({Y | Z: 1, Y: -1})
        d = {Z: 1}
        for m, c in ((X, -1), (Y, -1), (X | Y, 1)):
            d[m] = d.get(m, 0) + c
        out.append({m: c for m, c in d.items() if c})
    return out


D_SH = D_ING = R_SH = R_ING = A_SH = A_ING = None


def setup():
    global D_SH, D_ING, R_SH, R_ING, A_SH, A_ING
    if D_SH is not None:
        return
    D_SH = [('H', i) for i in range(N8)]
    for i, j in itertools.combinations(range(N8), 2):
        rest = [e for e in range(N8) if e not in (i, j)]
        for kk in range(1 << len(rest)):
            D_SH.append(('M', i, j, sum(1 << rest[t] for t in range(len(rest)) if kk >> t & 1)))
    D_SH = [d for d in D_SH if any(m & ZZ for m in row_of(d))]
    D_ING, seen = [], set()
    for assign in itertools.product(range(6), repeat=N8):
        m = [0] * 6
        for e, t in enumerate(assign):
            m[t] |= 1 << e
        A, B, C, D, K = m[:5]
        if not (A and B and C and D) or A > B or C > D or not ((A | B | C | D | K) & ZZ):
            continue
        key = tuple(sorted(row_of(('I', A, B, C, D, K)).items()))
        if key and key not in seen:
            seen.add(key)
            D_ING.append(('I', A, B, C, D, K))
    R_SH = [row_of(d) for d in D_SH]
    R_ING = [row_of(d) for d in D_ING]

    def mat(rows):
        ri, ci, vv = [], [], []
        for k, row in enumerate(rows):
            for m, c in row.items():
                ri.append(k); ci.append(m - 1); vv.append(-c)
        return csr_matrix(coo_matrix((vv, (ri, ci)), shape=(len(rows), NV8)))
    A_SH, A_ING = mat(R_SH), mat(R_ING)


def _primal_active(h, P1, P2):
    """cutting-plane primal; returns the boolean mask of Ingleton rows in use when it becomes infeasible, or None."""
    rows, b = [], []
    for m in range(1, 64):
        rows.append({m: 1}); b.append(float(h[m]))
    rows += ci_rows(P1, P2); b += [0.0] * 6
    ri, ci, vv = [], [], []
    for k, row in enumerate(rows):
        for m, c in row.items():
            ri.append(k); ci.append(m - 1); vv.append(c)
    Aeq = csr_matrix(coo_matrix((vv, (ri, ci)), shape=(len(rows), NV8)))
    active = np.zeros(A_ING.shape[0], dtype=bool)
    for _ in range(100):
        A = vstack([A_SH, A_ING[np.where(active)[0]]]).tocsr() if active.any() else A_SH
        res = linprog(np.zeros(NV8), A_ub=A, b_ub=np.zeros(A.shape[0]), A_eq=Aeq, b_eq=np.array(b),
                      bounds=[(None, None)] * NV8, method='highs', options=OPTS)
        if res.status == 2:
            return active
        if res.status != 0:
            return None
        viol = A_ING @ res.x
        new = np.where((viol > 1e-7) & ~active)[0]
        if len(new) == 0:
            return None
        active[new[np.argsort(-viol[new])][:400]] = True
    return None


def _farkas(h, P1, P2, idx_ing):
    """dual LP over the Shannon rows, the Ingleton rows idx_ing and the six equalities; exact solution or None."""
    descs = D_SH + [D_ING[k] for k in idx_ing]
    rows = R_SH + [R_ING[k] for k in idx_ing]
    K = len(rows)
    cols = rows + ci_rows(P1, P2)
    zmasks = [m for m in range(1, 256) if m & ZZ]
    zi = {m: t for t, m in enumerate(zmasks)}
    ri, ci, vv = [], [], []
    phi = np.zeros(len(cols))
    for k, row in enumerate(cols):
        for m, c in row.items():
            if m & ZZ:
                ri.append(zi[m]); ci.append(k); vv.append(c)
            else:
                phi[k] += c * int(h[m])
    for k in range(K):
        ri.append(len(zmasks)); ci.append(k); vv.append(1)
    Aeq = csr_matrix(coo_matrix((vv, (ri, ci)), shape=(len(zmasks) + 1, len(cols))))
    beq = np.zeros(len(zmasks) + 1); beq[-1] = 1
    res = linprog(phi, A_eq=Aeq, b_eq=beq, bounds=[(0, None)] * K + [(None, None)] * 6, method='highs', options=OPTS)
    if res.status != 0 or res.fun > -1e-9:
        return None
    x = res.x
    supp = [k for k in range(K) if x[k] > 1e-10] + list(range(K, K + 6))
    A_s = [[cols[k].get(m, 0) for k in supp] for m in zmasks] + [[phi[k] for k in supp]]
    b_s = [0] * len(zmasks) + [-1]
    sol = exact_solve(A_s, b_s, [x[k] / (-res.fun) for k in supp])
    if sol is None or any(sol[t] < 0 for t, k in enumerate(supp) if k < K):
        return None
    return [(descs[k] if k < K else ('C', k - K), sol[t]) for t, k in enumerate(supp) if sol[t] != 0]


def s_form(terms, P1, P2):
    """integer coefficients (gcd 1) of F'(S) = F(h_norm(S)) over the 31 visible coordinates."""
    cis = ci_rows(P1, P2)
    tot = {}
    for d, v in terms:
        for m, c in (cis[d[1]] if d[0] == 'C' else row_of(d)).items():
            tot[m] = tot.get(m, 0) + v * c
    F = {m: c for m, c in tot.items() if c != 0 and not m & ZZ}
    coef = []
    for i in range(31):
        e = [0] * 31; e[i] = 1
        he = stabcert.h_norm(e)
        coef.append(sum(c * he[m] for m, c in F.items()))
    den = 1
    for c in coef:
        den = den * Fraction(c).denominator // math.gcd(den, Fraction(c).denominator)
    ints = [int(c * den) for c in coef]
    g = 0
    for c in ints:
        g = math.gcd(g, abs(c))
    return [c // g for c in ints] if g else ints


def certify(r, P1, P2):
    """(terms, F_S) with exact Fractions, checked by epr1kit.stabcert.verify_exclusion_ci2, or None."""
    setup()
    P1, P2 = tuple(int(t) for t in P1), tuple(int(t) for t in P2)
    h = stabcert.h_norm(r)
    active = _primal_active(h, P1, P2)
    tries = []
    if active is not None:
        tries.append(list(np.where(active)[0]))
    tries.append(list(range(len(D_ING))))
    for idx in tries:
        terms = _farkas(h, P1, P2, idx)
        if terms is None:
            continue
        F_S = s_form(terms, P1, P2)
        ok, val = stabcert.verify_exclusion_ci2(r, P1, P2, terms, F_S=F_S)
        if ok:
            return terms, F_S
    return None


def _work(task):
    k, r, P1, P2 = task
    t0 = time.time()
    out = certify(np.array(r, dtype=np.int64), P1, P2)
    if out is None:
        return dict(index=k, ok=False, seconds=round(time.time() - t0, 1))
    terms, F_S = out
    return dict(index=k, ok=True, ray=[int(x) for x in r], P1=list(P1), P2=list(P2), F_S=[int(c) for c in F_S],
                F_S_at_ray=int(sum(c * int(x) for c, x in zip(F_S, r))),
                terms=[[list(d), f'{v.numerator}/{v.denominator}'] for d, v in terms], seconds=round(time.time() - t0, 1))


def check_file(rays, path, log=print):
    R = np.load(rays).astype(np.int64)
    d = json.load(open(path))
    bad = 0
    for c in d['certificates']:
        r = R[c['index']]
        terms = [(tuple(t[0]), Fraction(t[1])) for t in c['terms']]
        ok, val = stabcert.verify_exclusion_ci2(r, c['P1'], c['P2'], terms, F_S=c['F_S'])
        ok &= [int(x) for x in r] == c['ray'] and sum(a * int(b) for a, b in zip(c['F_S'], r)) == c['F_S_at_ray'] < 0
        bad += not ok
    log(f"{len(d['certificates'])} certificates re-checked with epr1kit.stabcert.verify_exclusion_ci2: "
        f"{'ALL OK' if not bad else str(bad) + ' FAILED'}")
    return bad == 0


def self_test(log=print):
    """rows excluded in R30 (one common information) stay excluded with a second one added, plain and nested: each
    must get a certificate that the checker (epr1kit.stabcert.verify_exclusion_ci2) accepts and rejects when tampered
    with in four ways (an inequality multiplier made negative, a common-information term dropped, an inequality
    multiplier doubled, the certificate applied to a known extreme ray); a known extreme ray must get none."""
    setup()
    sys.path.insert(0, HERE)
    import s42_ci2_ingleton as E
    rng = random.Random(43)
    d = json.load(open(os.path.join(ROOT, 'data', 's22_ingleton_exclusions_a2i.json')))
    K = np.load(os.path.join(ROOT, 'data', 'stab5_extreme_reps.npy')).astype(np.int64)
    ok, n, rej = True, 0, 0
    for c in rng.sample(d['certificates'], 2):
        r = np.array(c['ray'], dtype=np.int64)
        hh = [int(x) for x in stabcert.h_norm(r)]
        P1 = (c['X'], c['Y'])
        for P2 in (rng.choice([p for p in E.first_pairs(hh, 2) if p != P1]), rng.choice(E.nested_pairs(2))):
            out = certify(r, P1, P2); n += 1
            if out is None:
                ok = False
                continue
            terms, F_S = out
            k = max((i for i, (t, v) in enumerate(terms) if t[0] != 'C'), key=lambda i: terms[i][1])
            j = max((i for i, (t, v) in enumerate(terms) if t[0] == 'C'), key=lambda i: abs(terms[i][1]))
            neg = list(terms); neg[k] = (neg[k][0], -neg[k][1])
            dbl = list(terms); dbl[k] = (dbl[k][0], 2 * dbl[k][1])
            tampered = [(r, neg), (r, terms[:j] + terms[j + 1:]), (r, dbl), (K[rng.randrange(len(K))], terms)]
            for rr, tt in tampered:
                rejected = not stabcert.verify_exclusion_ci2(rr, P1, P2, tt, F_S=F_S)[0]
                ok &= rejected; rej += rejected
    i = rng.randrange(len(K))
    P1, P2 = rng.choice(E.combos([int(x) for x in stabcert.h_norm(K[i])], 2, True))
    ok &= certify(K[i], P1, P2) is None
    log(f's43 self-test: {n} certificates for rows excluded in R30 with a second common information, accepted by the '
        f'checker; {rej} of {4 * n} tampered versions rejected; none for a known extreme ray: {ok}')
    return ok


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('rays', nargs='?')
    ap.add_argument('--ci2', default=None, help='output of s42_ci2_ingleton.py')
    ap.add_argument('--out', default=None)
    ap.add_argument('--workers', type=int, default=None)
    ap.add_argument('--check', default=None)
    ap.add_argument('--self-test', action='store_true')
    a = ap.parse_args()
    if a.self_test:
        sys.exit(0 if self_test() else 1)
    if a.check:
        sys.exit(0 if check_file(a.rays, a.check) else 1)
    if not (a.rays and a.ci2 and a.out):
        ap.error('rays, --ci2 and --out are required')
    setup()
    R = np.load(a.rays).astype(np.int64)
    src = json.load(open(a.ci2))
    tasks = [(int(k), R[int(k)].tolist(), v['infeasible'][0], v['infeasible'][1])
             for k, v in src['results'].items() if v['infeasible']]
    tasks.sort()
    print(f'{len(tasks)} rows reported infeasible in {a.ci2}', flush=True)
    certs, failed, t0 = [], [], time.time()
    with mp.get_context('fork').Pool(a.workers or os.cpu_count() or 1) as pool:
        for out in pool.imap_unordered(_work, tasks):
            if out['ok']:
                certs.append(out)
                print(f"row {out['index']}: certified ({len(out['terms'])} terms), F'(r) = {out['F_S_at_ray']} "
                      f"({out['seconds']} s)", flush=True)
            else:
                failed.append(out['index'])
                print(f"row {out['index']}: NO exact certificate ({out['seconds']} s)", flush=True)
    certs.sort(key=lambda c: c['index'])
    json.dump(dict(rays=a.rays, source=a.ci2, certificates=certs, failed=sorted(failed),
                   note='rows refer to the rays file; check with --check'), open(a.out, 'w'), indent=1)
    print(f'{len(certs)} certified, {len(failed)} failed ({time.time() - t0:.0f} s) -> {a.out}', flush=True)


if __name__ == '__main__':
    main()
