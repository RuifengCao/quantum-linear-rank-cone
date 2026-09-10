"""EPR-1 R3 core: verified recipes from R1/R2, ported without absolute-path deps.

Conventions (frozen, do not change):
  * visible parties A..E = bits 0..4;  purifier F = bit 5;  FULL6 = 63
  * 31 coordinates: index = mask - 1 over nonempty subsets of {A..E}
  * purity projection: ground mask m over 6 subsystems -> rep = m if F not in m
    else complement(m); coordinate = rep - 1; full set / empty -> dropped
  * canonical sha of a row family: sha256 of np.unique(int64 rows, axis=0)
    cast to int8 bytes (all shipped families have |coef| <= 9)
All builders reproduce the R2 shipped arrays byte-identically (see selftest).
"""
import itertools, hashlib, os
import numpy as np

NVIS, DIM, F_BIT, FULL6 = 5, 31, 5, 63
DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data')

# ---------------------------------------------------------------- data access
def load(name):
    return np.load(os.path.join(DATA, name + '.npy')).astype(np.int64)

def sha_rows(a):
    u = np.unique(np.asarray(a, dtype=np.int64), axis=0)
    if u.size and int(np.abs(u).max()) > 127:
        raise ValueError('sha_rows: |coef| > 127 breaks the int8 canonical encoding')
    return hashlib.sha256(np.ascontiguousarray(u.astype(np.int8)).tobytes()).hexdigest()

# ---------------------------------------------------------------- purity & perms
def purity_index():
    P = np.full(64, -1, dtype=np.int64)
    for m in range(1, FULL6):
        rep = m if not (m >> F_BIT) & 1 else (FULL6 ^ m)
        P[m] = rep - 1
    return P

def perms31_s5():
    """120 coordinate permutations of the 31 coords (classical side, no purity)."""
    out = []
    for p in itertools.permutations(range(5)):
        arr = np.zeros(DIM, dtype=np.int64)
        for m in range(1, 32):
            pm = 0
            for i in range(5):
                if (m >> i) & 1:
                    pm |= 1 << p[i]
            arr[m - 1] = pm - 1
        out.append(arr)
    return np.stack(out)

def perms31_s6():
    """720 coordinate permutations induced by S6 on the purified coords."""
    out = []
    for p in itertools.permutations(range(6)):
        arr = np.zeros(DIM, dtype=np.int64)
        for X in range(1, 32):
            gm = 0
            for i in range(5):
                if (X >> i) & 1:
                    gm |= 1 << p[i]
            rep = gm if not (gm >> F_BIT) & 1 else (FULL6 ^ gm)
            arr[X - 1] = rep - 1
        out.append(arr)
    return np.stack(out)

def class_reps(rows, perms, offset=64):
    """S-orbit canonical forms. Returns (list of canonical bytes per row, dict canon->first row idx)."""
    R64 = np.asarray(rows, dtype=np.int64)
    if R64.size and not (0 <= int(R64.min()) + offset and int(R64.max()) + offset <= 255):
        raise ValueError(f'class_reps: values+offset outside uint8 (min {R64.min()}, '
                         f'max {R64.max()}, offset {offset})')
    R = (R64 + offset).astype(np.uint8)
    best = [None] * R.shape[0]
    for k in range(perms.shape[0]):
        blk = R[:, perms[k]]
        for i in range(R.shape[0]):
            b = blk[i].tobytes()
            if best[i] is None or b < best[i]:
                best[i] = b
    reps = {}
    for i, b in enumerate(best):
        reps.setdefault(b, i)
    return best, reps

# ---------------------------------------------------------------- elementals
def _addrow(dst, pairs):
    v = [0] * 32
    for m, c in pairs:
        v[m] += c
    if any(v[1:]):
        dst.add(tuple(v[1:]))

def elementals5():
    """classical side: submodularity(80), grouped weak monotonicity(195), monotonicity(5)."""
    submod, wm, mono = set(), set(), set()
    for i, j in itertools.combinations(range(5), 2):
        others = [k for k in range(5) if k not in (i, j)]
        for r in range(4):
            for Ks in itertools.combinations(others, r):
                K = sum(1 << k for k in Ks)
                _addrow(submod, [((1 << i) | K, 1), ((1 << j) | K, 1),
                                 ((1 << i) | (1 << j) | K, -1)] + ([(K, -1)] if K else []))
    for assign in itertools.product(range(4), repeat=5):
        Am = sum(1 << i for i in range(5) if assign[i] == 0)
        Bm = sum(1 << i for i in range(5) if assign[i] == 1)
        Cm = sum(1 << i for i in range(5) if assign[i] == 2)
        if not Am or not Bm or not Cm:
            continue
        _addrow(wm, [(Am | Cm, 1), (Bm | Cm, 1), (Am, -1), (Bm, -1)])
    for i in range(5):
        _addrow(mono, [(31, 1), (31 & ~(1 << i), -1)])
    arr = lambda s: np.array(sorted(s), dtype=np.int64)
    return arr(submod), arr(wm), arr(mono)

def ssa6_pure(P=None):
    """elemental SSA on the 6 purified subsystems, projected to 31 coords (<=240 rows)."""
    if P is None:
        P = purity_index()
    rows = set()
    for i, j in itertools.combinations(range(6), 2):
        rest = [k for k in range(6) if k not in (i, j)]
        for r in range(5):
            for Ks in itertools.combinations(rest, r):
                K = sum(1 << k for k in Ks)
                acc = [0] * DIM
                for m, c in (((1 << i) | K, 1), ((1 << j) | K, 1),
                             ((1 << i) | (1 << j) | K, -1), (K, -1)):
                    idx = P[m] if m else -1
                    if idx >= 0:
                        acc[idx] += c
                if any(acc):
                    rows.add(tuple(acc))
    return np.array(sorted(rows), dtype=np.int64)

# ---------------------------------------------------------------- substitution builders
def template_terms(T):
    return [[(m + 1, int(t[m])) for m in range(DIM) if t[m]] for t in T]

def subst_assignment_range(terms, P, lo, hi, nonempty):
    """workhorse: substitute templates over raw assignments a in [lo,hi) of {0..5}^6.
    a encodes base-6 digits; digit e = group of element e (5 = unused).
    Returns a set of int8 row bytes."""
    cap = max((sum(abs(cf) for _, cf in tl) for tl in terms), default=0)
    if cap > 127:
        raise ValueError(f'template coefficient mass {cap} > 127: int8 rows would wrap')
    rows = set()
    for code in range(lo, hi):
        a, x = [0] * 6, code
        for e in range(6):
            a[e] = x % 6
            x //= 6
        G = [0] * 5
        for e in range(6):
            if a[e] < 5:
                G[a[e]] |= 1 << e
        if nonempty and any(g == 0 for g in G):
            continue
        tbl = [0] * 32
        for fm in range(1, 32):
            g = 0
            for i in range(5):
                if (fm >> i) & 1:
                    g |= G[i]
            tbl[fm] = g
        for tl in terms:
            acc = [0] * DIM
            for fm, c in tl:
                idx = P[tbl[fm]]
                if idx >= 0:
                    acc[idx] += c
            if any(acc):
                rows.add(np.array(acc, dtype=np.int8).tobytes())
    return rows

def rows_from_bytes(byteset):
    if not byteset:
        return np.zeros((0, DIM), dtype=np.int64)
    return np.unique(np.frombuffer(b''.join(sorted(byteset)), dtype=np.int8)
                     .reshape(-1, DIM).astype(np.int64), axis=0)

def build_qlr(variant, workers=1):
    """variant in {'v3','pure','pure2','pure28'}; parallel over the assignment space.
    pure28 (templates27 + Ingleton(39)) is the CURRENT cone; pure is historical."""
    if variant == 'v3':
        SUB, WM, MONO = elementals5()
        return np.unique(np.vstack([load('M_LR'), SUB, WM]), axis=0)
    T = load('templates28') if variant.endswith('28') else load('templates27')
    terms = template_terms(T)
    P = purity_index()
    nonempty = variant in ('pure', 'pure28')
    total = 6 ** 6
    if workers <= 1:
        rows = subst_assignment_range(terms, P, 0, total, nonempty)
    else:
        from .parallel import pmap_byteset
        step = max(1, total // (workers * 8))
        chunks = [(lo, min(lo + step, total)) for lo in range(0, total, step)]
        rows = pmap_byteset(_subst_worker, [(terms, lo, hi, nonempty) for lo, hi in chunks],
                            workers=workers)
    for r in ssa6_pure(P):
        rows.add(np.array(r, dtype=np.int8).tobytes())
    return rows_from_bytes(rows)

def _subst_worker(args):
    terms, lo, hi, nonempty = args
    return subst_assignment_range(terms, purity_index(), lo, hi, nonempty)

def build_clr():
    SUB, WM, MONO = elementals5()
    return np.unique(np.vstack([load('M_LR'), SUB, MONO]), axis=0)

# ---------------------------------------------------------------- graph states
def graphstate_vecs6():
    """all simple GF(2) graphs on 6 vertices -> cut-rank entropy vectors, purified (760)."""
    edges = list(itertools.combinations(range(6), 2))
    P = purity_index()
    reps = [m for m in range(1, FULL6) if not (m >> F_BIT) & 1]   # 31 rep masks
    out = set()
    for code in range(1 << len(edges)):
        A = np.zeros((6, 6), dtype=np.uint8)
        for b, (i, j) in enumerate(edges):
            if (code >> b) & 1:
                A[i, j] = A[j, i] = 1
        v = np.zeros(DIM, dtype=np.int64)
        for m in reps:
            S = [i for i in range(6) if (m >> i) & 1]
            Sc = [i for i in range(6) if not (m >> i) & 1]
            v[P[m]] = gf2_rank(A[np.ix_(S, Sc)])
        out.add(v.tobytes())
    return np.unique(np.frombuffer(b''.join(sorted(out)), dtype=np.int64)
                     .reshape(-1, DIM), axis=0)

def gf2_rank(A):
    rows = [int(''.join(map(str, r)), 2) for r in A] if A.size else []
    r = 0
    for c in range(A.shape[1] if A.size else 0):
        bit = 1 << (A.shape[1] - 1 - c)
        piv = next((i for i in range(r, len(rows)) if rows[i] & bit), None)
        if piv is None:
            continue
        rows[r], rows[piv] = rows[piv], rows[r]
        for i in range(len(rows)):
            if i != r and rows[i] & bit:
                rows[i] ^= rows[r]
        r += 1
    return r

# ---------------------------------------------------------------- judge & ranks
def judge(H, V, chunk=4096):
    """returns (min value per H-row over V, count of violating pairs per H-row)."""
    H = np.asarray(H, dtype=np.int64)
    V = np.asarray(V, dtype=np.int64)
    mins = np.full(H.shape[0], np.iinfo(np.int64).max, dtype=np.int64)
    cnts = np.zeros(H.shape[0], dtype=np.int64)
    for s in range(0, V.shape[0], chunk):
        S = H @ V[s:s + chunk].T
        mins = np.minimum(mins, S.min(axis=1))
        cnts += (S < 0).sum(axis=1)
    return mins, cnts

P61 = (1 << 31) - 1     # Mersenne prime, int64-safe products

def rank_mod_p(A, p=P61):
    A = (np.asarray(A, dtype=np.int64) % p).copy()
    r, rows, cols = 0, A.shape[0], A.shape[1]
    for c in range(cols):
        nz = np.nonzero(A[r:, c])[0]
        if nz.size == 0:
            continue
        piv = r + nz[0]
        A[[r, piv]] = A[[piv, r]]
        A[r] = (A[r] * pow(int(A[r, c]), p - 2, p)) % p
        col = A[:, c].copy(); col[r] = 0
        A -= np.outer(col, A[r]); A %= p
        r += 1
        if r == rows:
            break
    return r

def rank_exact(A):
    """fraction-free Bareiss with python ints (small matrices: tight-row sets x 31)."""
    M = [[int(x) for x in row] for row in np.asarray(A)]
    n, m = len(M), len(M[0]) if len(M) else 0
    r, prev = 0, 1
    for c in range(m):
        piv = next((i for i in range(r, n) if M[i][c]), None)
        if piv is None:
            continue
        M[r], M[piv] = M[piv], M[r]
        for i in range(r + 1, n):
            for j in range(c + 1, m):
                num = M[i][j] * M[r][c] - M[i][c] * M[r][j]
                q, rem = divmod(num, prev)
                if rem:            # Bareiss guarantees exact division; a remainder
                    raise ArithmeticError('Bareiss exactness violated -- report this')
                M[i][j] = q
            M[i][c] = 0
        prev = M[r][c]
        r += 1
        if r == n:
            break
    return r

def tight_rank_report(H, rays, exact_confirm=True):
    """for each ray: (#tight rows, rank of tight rows, extreme?). Extreme iff rank == DIM-1."""
    H = np.asarray(H, dtype=np.int64)
    out = []
    for ri, r in enumerate(np.asarray(rays, dtype=np.int64)):
        tight = H[(H @ r) == 0]
        rk = rank_mod_p(tight)
        if exact_confirm and rk != DIM - 1:
            rk = rank_exact(tight)
        out.append((ri, int(tight.shape[0]), int(rk), rk == DIM - 1))
    return out

# ---------------------------------------------------------------- cone I/O & certificates
def write_normaliz(H, path):
    H = np.asarray(H, dtype=np.int64)
    with open(path, 'w') as f:
        f.write(f'amb_space {DIM}\ninequalities {H.shape[0]}\n')
        for row in H:
            f.write(' '.join(map(str, row)) + '\n')
        f.write('ExtremeRays\n')

def write_lrs(H, path):
    H = np.asarray(H, dtype=np.int64)
    with open(path, 'w') as f:
        f.write('cone\nH-representation\nbegin\n')
        f.write(f'{H.shape[0]} {DIM + 1} integer\n')
        for row in H:
            f.write('0 ' + ' '.join(map(str, row)) + '\n')
        f.write('end\n')

def decompose_certificate(r, H, pool, exclude_parallel=True):
    """LP: r = sum lambda_i v_i, lambda>=0, v_i in pool restricted to rows tight wherever r is
    tight (facet-restriction argument). Feasible -> non-extremality certificate."""
    from scipy.optimize import linprog
    H = np.asarray(H, dtype=np.int64); r = np.asarray(r, dtype=np.int64)
    tight = H[(H @ r) == 0]
    i0 = int(np.nonzero(r)[0][0])
    V = []
    for v in np.asarray(pool, dtype=np.int64):
        if np.any(tight @ v != 0):
            continue
        if exclude_parallel:
            # v parallel to r  <=>  v * r[i0] == r * v[i0] elementwise (cross-multiplication);
            # also drops v == 0.  Keeping multiples of r would make the LP trivially feasible.
            if np.array_equal(v * int(r[i0]), r * int(v[i0])):
                continue
        V.append(v)
    if len(V) < 2:
        return None
    V = np.array(V, dtype=np.int64)
    res = linprog(np.zeros(V.shape[0]), A_eq=V.T.astype(float), b_eq=r.astype(float),
                  bounds=[(0, None)] * V.shape[0], method='highs')
    return (V, res.x) if res.status == 0 else None


def rank_exact_gram(rows):
    """Exact rank via the Gram matrix: rank(T) = n_cols - dim ker(T^T T) over Q.
    Fraction elimination on the small square Gram (milliseconds even for tens of
    thousands of rows).  int64-safe for our families (|entries| small).
    Added R11; R7 report claimed this earlier -- that patch silently missed
    (errata C4).  Bareiss rank_exact remains the cross-check."""
    from fractions import Fraction
    T = np.asarray(rows, dtype=np.int64)
    G = T.T @ T
    M = [[Fraction(int(x)) for x in row] for row in G]
    n = G.shape[0]
    r = 0
    for cc in range(n):
        piv = next((i for i in range(r, n) if M[i][cc] != 0), None)
        if piv is None:
            continue
        M[r], M[piv] = M[piv], M[r]
        pv = M[r][cc]
        M[r] = [x / pv for x in M[r]]
        for i in range(n):
            if i != r and M[i][cc] != 0:
                f = M[i][cc]
                M[i] = [a - f * b for a, b in zip(M[i], M[r])]
        r += 1
    return r


def sha_rows_wide(rows):
    """Canonical sha256 for integer families with |coef| up to 32767 (int16
    little-endian encoding; rows sorted).  Added R15 for the QLR_5 partial
    catalogue whose coordinates reach 162 (> int8 guard of sha_rows)."""
    import hashlib
    A = np.asarray(rows, dtype=np.int64)
    if np.abs(A).max() > 32767:
        raise ValueError('sha_rows_wide: |coef| > 32767')
    A = A.astype('<i2')
    A = A[np.lexsort(A.T[::-1])]
    return hashlib.sha256(A.tobytes()).hexdigest()
