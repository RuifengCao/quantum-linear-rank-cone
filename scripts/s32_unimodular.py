#!/usr/bin/env python3
"""s32: unimodular realisations of the known extreme rays of Stab_5 (R38; roadmap item S5(a)).

A realisation of kappa * r is *unimodular* if it is an integer symmetric matrix W with zero diagonal on
N = kappa * (r(A) + ... + r(F)) vertices, party q owning kappa * r(q) consecutive vertices (purifier F last), such
that every visible cut block W[Q_X, Q_{N \\ X}] (X a nonempty set of parties among A..E) has Smith invariants in
{0, 1} and rational rank kappa * r(X).  Such a matrix has the same cut ranks over Q and over every finite field, and
its image in Z_d^{|Q_X|} is free of rank kappa * r(X) for every d.  By the cut-rank formula of the paper
(Appendix A 1: S(X) = log |W[Q_X, Q_{N \\ X}] Z_d^{N \\ X}| for a graph state of Z_d qudits), W mod d is a graph state
of Z_d qudits with entropy vector kappa * r in units of log d, for EVERY d >= 2 at once.  This is the five-party
version of the four-party certificates behind Theorem 2 of the paper (Theorem 7 of the R32 note,
data/s30_stab4_unimodular.npz).

Data: data/s32_stab5_unimodular.npz holds one certificate per solved row of data/stab5_extreme_reps.npy (the 641
known extreme-ray orbits of Stab_5): row, rep (= the catalogue row), kappa, sizes (party sizes, = kappa * singles),
N, W (int8, padded to 64 x 64), method (the strategy that found it) and detail (strategy, seed, time budget).

Search (`--search`): cc/s32_unimodular.c anneals the inter-party entries over {-1, 0, 1} so that the cut ranks are
right over GF(2), GF(3) and modulo 2^31 - 1 (a stand-in for Q), adding any prime below 1000 (and 2^61 - 1) that fails
once these are satisfied.  The search proves nothing by itself; every reported matrix is certified below.  Strategies,
run in the given order on the rows still unsolved:
  plain             start from the stored qubit certificate of the row (data/s16_realisations.npz,
                    results/2026-10-02_a2-ingleton/a2i_q2.npz, data/s24_qubit_l2_certs.npz; first all signs +1, then
                    random signs) and change one entry at a time; random starts for rows without one;
  neutral           the same with neutral moves: party-local unimodular basis changes and integer local
                    complementations, which keep every cut block in its Smith class and so never change the residual;
  ring-eisenstein, ring-gauss, ring-golden, ring-sqrt2, ring-sqrtm2
                    restriction of scalars from Z[theta] (theta^2 = s theta + t: Z[omega], Z[i], Z[(1 + sqrt 5)/2],
                    Z[sqrt 2], Z[sqrt -2]); every entry x + y theta becomes the symmetric 2 x 2 block G M(x + y theta),
                    G = [[0, 1], [1, s]] the Gram matrix of the unimodular trace form; kappa = 2;
  qutrit-start      start from a stored qutrit certificate (data/s24_qutrit_certs.npz, data/s21_qudit_realisations.npz):
                    the Eisenstein ring from a qutrit graph state at multiplicity 1 (Z[omega] / (sqrt -3) = GF(3)), or the
                    neutral search at kappa = 2 from one at multiplicity 2;
  gf3-search-start  first a GF(3)-only search at multiplicity 1 from random starts, then the Eisenstein ring from it.

Certification (`--check`, used by gate G36 of selftest.py) is exact and does not use the search code:
  M1  a diagonal form over Z by integer row and column operations: exactly kappa * r(X) nonzero entries, all +-1;
  M2  (with --full) the rank over Q by Fraction elimination, and t x t minors with gcd 1 (a minor D at the rational
      pivots and, for every prime q dividing D, a minor that is nonzero mod q), so the t-th determinantal divisor is 1;
  plus (with --full) the ranks over GF(p) for p <= 19, and the state vector of W mod d for d = 4 (certificates with at
  most 8 vertices) and d = 6 (at most 7 vertices): on every cut the Schmidt coefficients are flat and their number is
  d^(kappa r(X)).

  python3 scripts/s32_unimodular.py --check [--full]
  python3 scripts/s32_unimodular.py --search --seconds 10 --out FILE [--resume] [--rows 0,1,2] [--strategies ...]
  python3 scripts/s32_unimodular.py --self-test       # neutral moves keep Smith forms; ring expansion; Z_d formula
  python3 scripts/s32_unimodular.py --cover --out FILE --cover-out FILE [--rows ...] [--seconds S] [--rounds R]
                                    [--rings ring-eisenstein,...] [--kappa2]   # covering certificates (cover_search)
  python3 scripts/s32_unimodular.py --check-cover     # bad primes recomputed, GF(p) pieces, coverage (gate G37)
"""
import argparse, itertools, json, os, random, shutil, subprocess, sys, tempfile, time
from fractions import Fraction
from math import gcd
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, '..')
DATA = os.path.join(ROOT, 'data')
OUT_DEFAULT = os.path.join(DATA, 's32_stab5_unimodular.npz')
COVER_DEFAULT = os.path.join(DATA, 's32_stab5_cover.npz')
BIG31 = 2147483647
KERNEL = os.path.join(ROOT, 'cc', 's32_unimodular.c')
MAXN = 64
CROSS_PRIMES = (2, 3, 5, 7, 11, 13, 17, 19)

# ring strategies: theta^2 = s theta + t
RINGS = {'ring-eisenstein': (-1, -1), 'ring-gauss': (0, -1), 'ring-golden': (1, 1), 'ring-sqrt2': (0, 2),
         'ring-sqrtm2': (0, -2)}


# ----------------------------------------------------------------------------------------------- catalogue and starts
def singles(r):
    return [int(r[(1 << q) - 1]) for q in range(5)] + [int(r[30])]


_IDX = None


def canon(v):
    """lexicographically smallest image of a 31-vector (mask order) under S6 acting on A..F with purity."""
    global _IDX
    if _IDX is None:
        _IDX = np.zeros((720, 31), dtype=np.int64)
        for a, p in enumerate(itertools.permutations(range(6))):
            inv = [0] * 6
            for i, j in enumerate(p):
                inv[j] = i
            for M in range(1, 32):
                N = sum(1 << inv[i] for i in range(6) if M >> i & 1)
                _IDX[a, M - 1] = (63 ^ N) if N >> 5 & 1 else N
    ext = np.zeros(64, dtype=np.int64)
    ext[1:32] = v
    imgs = ext[_IDX]
    return tuple(int(x) for x in imgs[np.lexsort(imgs.T[::-1])[0]])


def catalogue():
    return np.load(os.path.join(DATA, 'stab5_extreme_reps.npy')).astype(np.int64)


def qubit_certificates(reps):
    """{row: [(lam, sizes, adjacency 0/1 list of lists, source)]} from the stored qubit certificates."""
    row_of = {canon(r): k for k, r in enumerate(reps)}
    out = {}

    def add(rep, lam, sizes, adj, tag):
        k = row_of[canon(np.asarray(rep, dtype=np.int64))]
        assert [int(x) for x in rep] == [int(x) for x in reps[k]], 'certificate is not on the catalogue representative'
        out.setdefault(k, []).append((int(lam), [int(s) for s in sizes], adj, tag))

    for path, tag in ((os.path.join(DATA, 's16_realisations.npz'), 's16'),
                      (os.path.join(ROOT, 'results', '2026-10-02_a2-ingleton', 'a2i_q2.npz'), 'a2i_q2')):
        z = np.load(path)
        for k in range(len(z['rep'])):
            sizes = [int(s) for s in z['sizes'][k]]; nv = sum(sizes)
            rows = [int(x) for x in z['rows'][k][:nv]]
            add(z['rep'][k], z['lam'][k], sizes, [[(rows[a] >> b) & 1 for b in range(nv)] for a in range(nv)], f'{tag}:{k}')
    z = np.load(os.path.join(DATA, 's24_qubit_l2_certs.npz'))
    for k in range(len(z['rep'])):
        sizes = [int(s) for s in z['sizes'][k]]; nv = sum(sizes)
        add(z['rep'][k], z['lam'][k], sizes, z['G'][k][:nv, :nv].astype(int).tolist(), f's24:{k}')
    return out


# ----------------------------------------------------------------------------------------------- exact certification
def cut_blocks(W, sizes):
    """(mask, block) for the 31 visible cuts; W a list of lists, parties in consecutive blocks, purifier last."""
    owner = [q for q, c in enumerate(sizes) for _ in range(c)]
    N = len(owner)
    for m in range(1, 32):
        X = [v for v in range(N) if owner[v] < 5 and m >> owner[v] & 1]
        Xs = set(X)
        Y = [v for v in range(N) if v not in Xs]
        yield m, [[W[a][b] for b in Y] for a in X]


def diagonal_form(B):
    """absolute values of the nonzero entries of a diagonal form U B V (U, V unimodular over Z)."""
    A = [list(map(int, row)) for row in B]
    m = len(A); n = len(A[0]) if m else 0
    out, r = [], 0
    while r < min(m, n):
        best = None
        for i in range(r, m):
            for j in range(r, n):
                v = A[i][j]
                if v and (best is None or abs(v) < best[0]):
                    best = (abs(v), i, j)
        if best is None:
            break
        _, i, j = best
        A[r], A[i] = A[i], A[r]
        for row in A:
            row[r], row[j] = row[j], row[r]
        while True:
            p, moved = A[r][r], False
            for i in range(r + 1, m):
                if A[i][r]:
                    q = A[i][r] // p
                    if q:
                        A[i] = [a - q * b for a, b in zip(A[i], A[r])]
                    if A[i][r]:
                        A[r], A[i] = A[i], A[r]; moved = True; break
            if moved:
                continue
            for j in range(r + 1, n):
                if A[r][j]:
                    q = A[r][j] // p
                    if q:
                        for row in A:
                            row[j] -= q * row[r]
                    if A[r][j]:
                        for row in A:
                            row[r], row[j] = row[j], row[r]
                        moved = True; break
            if moved:
                continue
            break
        out.append(abs(A[r][r])); r += 1
    return out


def _pivots(B, q=None):
    """pivot rows and columns of a row-echelon reduction over Q (q None) or GF(q)."""
    if q is None:
        M = [[Fraction(x) for x in row] for row in B]
    else:
        M = [[x % q for x in row] for row in B]
    m = len(M); n = len(M[0]) if m else 0
    rows, prow, pcol, r = list(range(m)), [], [], 0
    for c in range(n):
        piv = next((i for i in range(r, m) if M[i][c] != 0), None)
        if piv is None:
            continue
        M[r], M[piv] = M[piv], M[r]; rows[r], rows[piv] = rows[piv], rows[r]
        for i in range(r + 1, m):
            if M[i][c] != 0:
                if q is None:
                    f = M[i][c] / M[r][c]
                    M[i] = [a - f * b for a, b in zip(M[i], M[r])]
                else:
                    f = (M[i][c] * pow(M[r][c], -1, q)) % q
                    M[i] = [(a - f * b) % q for a, b in zip(M[i], M[r])]
        prow.append(rows[r]); pcol.append(c); r += 1
        if r == m:
            break
    return prow, pcol


def _det(A):
    """exact determinant of a square integer matrix (Bareiss)."""
    A = [list(map(int, row)) for row in A]
    n = len(A)
    if n == 0:
        return 1
    sign, prev = 1, 1
    for k in range(n - 1):
        if A[k][k] == 0:
            sw = next((i for i in range(k + 1, n) if A[i][k]), None)
            if sw is None:
                return 0
            A[k], A[sw] = A[sw], A[k]; sign = -sign
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                A[i][j] = (A[i][j] * A[k][k] - A[i][k] * A[k][j]) // prev
        prev = A[k][k]
    return sign * A[n - 1][n - 1]


def _prime_factors(n):
    n = abs(n); out, p = [], 2
    while p * p <= n:
        if n % p == 0:
            out.append(p)
            while n % p == 0:
                n //= p
        p += 1 if p == 2 else 2
    if n > 1:
        out.append(n)
    return out


def minors_certificate(B, t):
    """None unless the block has rational rank t and t x t minors with gcd 1; else the list of (rows, cols, det)."""
    if t == 0:
        return [] if not any(any(row) for row in B) else None
    prow, pcol = _pivots(B)
    if len(prow) != t:
        return None
    sub = lambda R, C: [[B[i][j] for j in C] for i in R]
    D = _det(sub(prow, pcol))
    out = [(prow, pcol, D)]
    for q in _prime_factors(D):
        r2, c2 = _pivots(B, q)
        if len(r2) < t:
            return None
        D2 = _det(sub(r2[:t], c2[:t]))
        if D2 % q == 0:
            return None
        out.append((r2[:t], c2[:t], D2))
    g = 0
    for _r, _c, d in out:
        g = gcd(g, d)
    return out if g == 1 else None


def structure_ok(r, kappa, sizes, W):
    N = sum(sizes)
    if list(sizes) != [kappa * s for s in singles(r)] or N > MAXN or len(W) != N:
        return False
    owner = [q for q, c in enumerate(sizes) for _ in range(c)]
    return all(len(row) == N for row in W) and all(
        W[a][b] == W[b][a] and (W[a][b] == 0 if owner[a] == owner[b] else True) for a in range(N) for b in range(N))


def certify(r, kappa, sizes, W, full=False):
    """dict of exact verdicts for one certificate (all True = unimodular realisation of kappa * r)."""
    res = {'structure': structure_ok(r, kappa, sizes, W), 'M1': True}
    if full:
        res.update(M2=True, gfp=True)
    if not res['structure']:
        return res
    for m, B in cut_blocks(W, sizes):
        t = kappa * int(r[m - 1])
        d = diagonal_form(B) if B and B[0] else []
        res['M1'] &= len(d) == t and all(x == 1 for x in d)
        if full:
            res['M2'] &= minors_certificate(B, t) is not None
            res['gfp'] &= all((len(_pivots(B, p)[0]) if B and B[0] else 0) == t for p in CROSS_PRIMES)
    return res


def zd_state_check(r, kappa, sizes, W, d, tol=1e-8):
    """W mod d as a graph state of Z_d qudits: for every cut, flat Schmidt spectrum with d^(kappa r(X)) terms."""
    N = sum(sizes)
    owner = [q for q, c in enumerate(sizes) for _ in range(c)]
    x = np.array(list(itertools.product(range(d), repeat=N)), dtype=np.int64)        # basis states, qudit 0 first
    Wm = np.array(W, dtype=np.int64) % d
    phase = np.einsum('ki,ij,kj->k', x, np.triu(Wm, 1), x) % d
    psi = np.exp(2j * np.pi * phase / d) / d ** (N / 2)
    psi = psi.reshape([d] * N)
    for m in range(1, 32):
        X = [v for v in range(N) if owner[v] < 5 and m >> owner[v] & 1]
        Y = [v for v in range(N) if v not in X]
        M = np.transpose(psi, X + Y).reshape(d ** len(X), d ** len(Y))
        s = np.linalg.svd(M, compute_uv=False)
        p = s[s > tol] ** 2
        want = d ** (kappa * int(r[m - 1]))
        if len(p) != want or np.abs(p - 1.0 / want).max() > 1e-6:
            return False
    return True


def content_sha(path=OUT_DEFAULT):
    """sha256 of the stored certificates (every array except the provenance texts method and detail, in file order),
    independent of the zip container."""
    import hashlib
    Z = np.load(path, allow_pickle=False)
    h = hashlib.sha256()
    for key in sorted(Z.files):
        if key in ('method', 'detail'):
            continue
        a = Z[key]
        h.update(key.encode())
        if a.dtype.kind == 'U':
            h.update('\x1f'.join(str(x) for x in a).encode())
        else:
            h.update(np.ascontiguousarray(a, dtype=a.dtype.newbyteorder('<')).tobytes())
    return h.hexdigest()


def load(path=OUT_DEFAULT):
    Z = np.load(path, allow_pickle=False)
    out = []
    for i in range(len(Z['row'])):
        N = int(Z['N'][i])
        out.append(dict(row=int(Z['row'][i]), rep=Z['rep'][i].astype(np.int64), kappa=int(Z['kappa'][i]),
                        sizes=[int(s) for s in Z['sizes'][i]], W=Z['W'][i][:N, :N].astype(int).tolist(),
                        method=str(Z['method'][i]), detail=str(Z['detail'][i]) if 'detail' in Z.files else ''))
    return out


def save(certs, path):
    certs = sorted(certs, key=lambda c: c['row'])
    K = len(certs)
    W = np.zeros((K, MAXN, MAXN), dtype=np.int8)
    for i, c in enumerate(certs):
        N = len(c['W'])
        W[i, :N, :N] = np.array(c['W'], dtype=np.int8)
    np.savez_compressed(path, row=np.array([c['row'] for c in certs], dtype=np.int16),
                        rep=np.array([c['rep'] for c in certs], dtype=np.int16),
                        kappa=np.array([c['kappa'] for c in certs], dtype=np.int8),
                        sizes=np.array([c['sizes'] for c in certs], dtype=np.int16),
                        N=np.array([len(c['W']) for c in certs], dtype=np.int16), W=W,
                        method=np.array([c['method'] for c in certs]),
                        detail=np.array([c.get('detail', c['method']) for c in certs]))


def bad_primes(r, kappa, sizes, W):
    """None unless every cut block has rational rank kappa * r(X); else the sorted list of the primes p for which
    W mod p is not a realisation of kappa * r: the primes dividing a Smith invariant of some cut block (a block of
    rational rank t has rank t over GF(p) exactly when p does not divide its t-th determinantal divisor, the product
    of the nonzero entries of any diagonal form)."""
    if not structure_ok(r, kappa, sizes, W):
        return None
    S = set()
    for m, B in cut_blocks(W, sizes):
        t = kappa * int(r[m - 1])
        d = diagonal_form(B) if B and B[0] else []
        if len(d) != t:
            return None
        for x in d:
            if x != 1:
                S |= set(_prime_factors(x))
    return sorted(S)


def gf_realisation_ok(r, kappa, sizes, W, p):
    """W mod p realises kappa * r over GF(p) (structure and every cut rank)."""
    return structure_ok(r, kappa, sizes, W) and all(
        (len(_pivots(B, p)[0]) if B and B[0] else 0) == kappa * int(r[m - 1]) for m, B in cut_blocks(W, sizes))


def covered(pieces):
    """pieces: list of ('Q', bad prime list) or ('GF', p).  True when every prime is covered: the intersection of the
    bad-prime sets of the Q pieces lies in the set of primes of the GF pieces (at least one Q piece is needed)."""
    qs = [set(b) for kind, b in pieces if kind == 'Q']
    if not qs:
        return False
    common = set.intersection(*qs)
    return common <= {p for kind, p in pieces if kind == 'GF'}


ZD_MAX_VERTICES = {4: 8, 6: 7}      # state vectors of d^N amplitudes: 65,536 and 279,936


def check(path=OUT_DEFAULT, full=False, zd=(4, 6), log=print):
    """exact re-verification of the stored certificates.  Returns a summary dict."""
    reps = catalogue()
    certs = load(path)
    rows = [c['row'] for c in certs]
    summ = dict(stored=len(certs), distinct_rows=len(set(rows)), rep_ok=True, ok=0, kappa={}, zd_checked=0, zd_ok=0,
                failures=[])
    t0 = time.time()
    for c in certs:
        k = c['row']
        rep_ok = 0 <= k < len(reps) and np.array_equal(c['rep'], reps[k])
        summ['rep_ok'] &= rep_ok
        res = certify(reps[k], c['kappa'], c['sizes'], c['W'], full=full)
        good = rep_ok and all(res.values())
        ds = [d for d in zd if len(c['W']) <= ZD_MAX_VERTICES.get(d, 0)] if full and good else []
        if ds:
            summ['zd_checked'] += 1
            zok = all(zd_state_check(reps[k], c['kappa'], c['sizes'], c['W'], d) for d in ds)
            summ['zd_ok'] += zok
            good &= zok
        if good:
            summ['ok'] += 1
            summ['kappa'][c['kappa']] = summ['kappa'].get(c['kappa'], 0) + 1
        else:
            summ['failures'].append((k, res))
    summ['seconds'] = round(time.time() - t0, 1)
    return summ


def save_cover(pieces, path):
    """pieces: dicts with row, kind ('Q' or 'GF'), prime (GF: p; Q: 0), kappa, sizes, W, bad (Q: list), detail."""
    pieces = sorted(pieces, key=lambda c: (c['row'], c['kind'] != 'Q', c['prime']))
    K = len(pieces)
    W = np.zeros((K, MAXN, MAXN), dtype=np.int8)
    for i, c in enumerate(pieces):
        N = len(c['W'])
        W[i, :N, :N] = np.array(c['W'], dtype=np.int8)
    np.savez_compressed(path, row=np.array([c['row'] for c in pieces], dtype=np.int16),
                        kind=np.array([c['kind'] for c in pieces]),
                        prime=np.array([c['prime'] for c in pieces], dtype=np.int64),
                        kappa=np.array([c['kappa'] for c in pieces], dtype=np.int8),
                        sizes=np.array([c['sizes'] for c in pieces], dtype=np.int16),
                        N=np.array([len(c['W']) for c in pieces], dtype=np.int16), W=W,
                        bad=np.array([','.join(map(str, c.get('bad', []))) for c in pieces]),
                        detail=np.array([c.get('detail', '') for c in pieces]))


def load_cover(path=COVER_DEFAULT):
    Z = np.load(path, allow_pickle=False)
    out = []
    for i in range(len(Z['row'])):
        N = int(Z['N'][i])
        bad = [int(x) for x in str(Z['bad'][i]).split(',') if x]
        out.append(dict(row=int(Z['row'][i]), kind=str(Z['kind'][i]), prime=int(Z['prime'][i]), kappa=int(Z['kappa'][i]),
                        sizes=[int(x) for x in Z['sizes'][i]], W=Z['W'][i][:N, :N].astype(int).tolist(), bad=bad,
                        detail=str(Z['detail'][i])))
    return out


def check_cover(path=COVER_DEFAULT, unimodular_path=OUT_DEFAULT):
    """exact re-verification of the covering certificates: every Q piece has the right rational cut ranks and exactly
    the stored bad primes, every GF piece is a realisation over GF(p), and the pieces of each row cover every prime.
    Returns a summary dict."""
    reps = catalogue()
    pieces = load_cover(path)
    uni = {c['row'] for c in load(unimodular_path)} if unimodular_path and os.path.exists(unimodular_path) else set()
    rows = {}
    summ = dict(pieces=len(pieces), pieces_ok=0, rows=0, rows_covered=0, overlap_unimodular=0, failures=[],
                open_primes={})
    t0 = time.time()
    for c in pieces:
        r = reps[c['row']]
        sizes_ok = c['sizes'] == [c['kappa'] * x for x in singles(r)]
        if c['kind'] == 'Q':
            b = bad_primes(r, c['kappa'], c['sizes'], c['W'])
            ok = sizes_ok and b is not None and b == c['bad'] and c['prime'] == 0
            entry = ('Q', b or [])
        else:
            ok = sizes_ok and c['prime'] in (2, 3, 5, 7, 11, 13) and gf_realisation_ok(r, c['kappa'], c['sizes'], c['W'], c['prime'])
            entry = ('GF', c['prime'])
        summ['pieces_ok'] += ok
        if not ok:
            summ['failures'].append((c['row'], c['kind'], c['prime']))
        rows.setdefault(c['row'], []).append(entry if ok else ('bad', None))
    for k, ps in rows.items():
        good = [p for p in ps if p[0] != 'bad']
        cov = len(good) == len(ps) and covered(good)
        summ['rows'] += 1
        summ['rows_covered'] += cov
        summ['overlap_unimodular'] += k in uni
        if not cov:
            summ['failures'].append((k, 'not covered'))
    summ['seconds'] = round(time.time() - t0, 1)
    return summ


# ----------------------------------------------------------------------------------------------- search
def build_kernel(tmp):
    exe = os.path.join(tmp, 's32_unimodular')
    flags = ['-O3', '-march=native', '-fopenmp']
    p = subprocess.run(['gcc'] + flags + ['-o', exe, KERNEL, '-lm'], capture_output=True, text=True)
    if p.returncode != 0:          # no OpenMP: one thread
        p = subprocess.run(['gcc', '-O3', '-march=native', '-o', exe, KERNEL, '-lm'], capture_output=True, text=True)
        if p.returncode != 0:
            raise SystemExit('gcc failed:\n' + p.stderr)
    return exe


def job_line(k, kappa, sizes, r, start=None):
    t = [kappa * int(r[m - 1]) for m in range(1, 32)]
    s = f'{k} {kappa} ' + ' '.join(map(str, sizes)) + ' ' + ' '.join(map(str, t))
    if start is None:
        return s + ' 0'
    return s + ' 1 ' + ' '.join(str(int(x)) for row in start for x in row)


def parse_found(text):
    out = []
    for line in text.splitlines():
        tok = line.split()
        if len(tok) > 8 and tok[2] == 'FOUND':
            N = int(tok[3])
            out.append((int(tok[0]), int(tok[1]), [[int(x) for x in row.split(',')] for row in tok[8:8 + N]], float(tok[4])))
    return out


def qutrit_certificates(reps):
    """{row: [(lam, sizes, W over GF(3) lifted to {-1, 0, 1})]} from data/s24_qutrit_certs.npz and
    data/s21_qudit_realisations.npz (only certificates on the catalogue representative)."""
    row_of = {canon(r): k for k, r in enumerate(reps)}
    out = {}
    for path, key in ((os.path.join(DATA, 's24_qutrit_certs.npz'), 'G'), (os.path.join(DATA, 's21_qudit_realisations.npz'), 'W')):
        z = np.load(path)
        for i in range(len(z['rep'])):
            if int(z['p'][i]) != 3:
                continue
            k = row_of.get(canon(z['rep'][i].astype(np.int64)))
            if k is None or [int(x) for x in z['rep'][i]] != [int(x) for x in reps[k]]:
                continue
            sizes = [int(s) for s in z['sizes'][i]]; nv = sum(sizes)
            W = [[x if x <= 1 else x - 3 for x in row] for row in z[key][i][:nv, :nv].astype(int).tolist()]
            out.setdefault(k, []).append((int(z['lam'][i]), sizes, W))
    return out


# strategy -> (kernel ring argument, sweeps, p_neutral, active primes, extra checks)
STRATEGIES = {
    'plain': ('0', '50', '0', '2,3,2147483647', '1'),
    'neutral': ('0', '200', '0.5', '2,3,2147483647', '1'),
    'qutrit-start': (None, '200', '0.5', '2,3,2147483647', '1'),
    'gf3-search-start': ('-1,-1', '200', '0.5', '2,3,2147483647', '1'),
}
for _name, (_s, _t) in RINGS.items():
    STRATEGIES[_name] = ('%d,%d' % (_s, _t), '200', '0.5', '2,3,2147483647', '1')


def kernel(exe, lines, seconds, seed, ring, sweeps, p_neu, primes, extra, dump=0):
    args = [exe, str(seconds), str(seed), '1', '1.0', '0.05', sweeps, '0.6', str(dump), primes, extra, ring, p_neu]
    return subprocess.run(args, input='\n'.join(lines) + '\n', capture_output=True, text=True, check=True).stdout


def run_strategy(exe, name, rows, reps, qcerts, seconds, seed, log, tcerts=None):
    """run one strategy on the given rows; returns (found list, raw kernel output)."""
    ring, sweeps, p_neu, primes, extra = STRATEGIES[name]
    t0 = time.time()
    if name == 'qutrit-start':
        # lambda_3 = 1: Eisenstein ring from the qutrit certificate (Z[omega] / (sqrt -3) = GF(3)); lambda_3 = 2: plain
        ring1, plain2 = [], []
        for k in rows:
            if k in tcerts:
                lam, sizes, W = min(tcerts[k], key=lambda c: c[0])
                (ring1 if lam == 1 else plain2).append(job_line(k, lam, sizes, reps[k], W))
        text = (kernel(exe, ring1, seconds, seed, '-1,-1', sweeps, p_neu, primes, extra) if ring1 else '') + \
               (kernel(exe, plain2, seconds, seed, '0', sweeps, p_neu, primes, extra) if plain2 else '')
        njobs = len(ring1) + len(plain2)
    elif name == 'gf3-search-start':
        # a GF(3) realisation at lambda = 1 from random starts, then the Eisenstein ring from it
        hunt = kernel(exe, [job_line(k, 1, singles(reps[k]), reps[k]) for k in rows], seconds, seed, '0', sweeps, p_neu,
                      '3', '0')
        starts = [job_line(k, kap, singles(reps[k]), reps[k], W) for k, kap, W, _s in parse_found(hunt)]
        text = kernel(exe, starts, seconds, seed + 1, ring, sweeps, p_neu, primes, extra) if starts else ''
        njobs = len(starts)
    else:
        lines = []
        for k in rows:
            r = reps[k]
            if k in qcerts:
                lam, sizes, adj, _src = min(qcerts[k], key=lambda c: c[0])
                if ring == '0':
                    lines.append(job_line(k, lam, sizes, r, adj))
                elif lam == 1:
                    lines.append(job_line(k, 1, sizes, r, adj))
                else:                       # ring jobs describe multiplicity 1 over R: random start
                    lines.append(job_line(k, 1, singles(r), r))
            else:
                lines.append(job_line(k, 1, singles(r), r))
        text = kernel(exe, lines, seconds, seed, ring, sweeps, p_neu, primes, extra)
        njobs = len(lines)
    found = parse_found(text)
    log(f'  {name} seed {seed}: {njobs} jobs, {len(found)} reported, {time.time() - t0:.0f} s')
    return found, text


def search(rows, strategies, seconds, seed, out, resume, log=print, raw=None):
    reps = catalogue()
    qcerts = qubit_certificates(reps)
    tcerts = qutrit_certificates(reps)
    certs = {c['row']: c for c in load(out)} if (resume and os.path.exists(out)) else {}
    tmp = tempfile.mkdtemp(prefix='s32_')
    try:
        exe = build_kernel(tmp)
        for i, name in enumerate(strategies):
            todo = [k for k in rows if k not in certs]
            if not todo:
                break
            found, text = run_strategy(exe, name, todo, reps, qcerts, seconds, seed + 2 * i, log, tcerts)
            if raw:
                with open(raw, 'a') as f:
                    f.write(f'# strategy {name} seed {seed + 2 * i} seconds {seconds}\n' + text)
            for k, kappa, W, secs in found:
                r = reps[k]
                sizes = [kappa * s for s in singles(r)]
                if all(certify(r, kappa, sizes, W).values()):
                    old = certs.get(k)
                    if old is None or old['kappa'] > kappa:
                        certs[k] = dict(row=k, rep=r, kappa=kappa, sizes=sizes, W=W, method=name,
                                        detail=f'{name}:seed{seed + 2 * i}:seconds{seconds}')
                else:
                    log(f'  row {k}: reported matrix FAILED exact certification (kept out)')
            save(list(certs.values()), out)
            log(f'  after {name}: {len(certs)} rows certified')
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    return certs


def doubled(W, s):
    """two copies of a realisation of r at multiplicity 1 (party q's 2 s[q] vertices: copy 1, then copy 2): a start
    at multiplicity 2."""
    off1 = [sum(s[:q]) for q in range(6)]
    owner = [q for q in range(6) for _ in range(s[q])]
    idx = [[2 * off1[owner[v]] + c * s[owner[v]] + (v - off1[owner[v]]) for v in range(len(W))] for c in range(2)]
    W2 = [[0] * (2 * len(W)) for _ in range(2 * len(W))]
    for c in range(2):
        for a in range(len(W)):
            for b in range(len(W)):
                W2[idx[c][a]][idx[c][b]] = W[a][b]
    return W2


def minimal_cover(qpieces, gf_primes):
    """a smallest subset of the rational pieces whose bad-prime sets intersect inside gf_primes (None if none)."""
    for size in range(1, len(qpieces) + 1):
        best = None
        for sub in itertools.combinations(qpieces, size):
            common = set.intersection(*[set(c['bad']) for c in sub])
            if common <= set(gf_primes):
                key = (sum(c['kappa'] for c in sub), sum(len(c['bad']) for c in sub))
                if best is None or key < best[0]:
                    best = (key, list(sub))
        if best:
            return best[1]
    return None


def cover_search(rows, seconds, seed, cover_out, unimodular_path=OUT_DEFAULT, rounds=3, log=print, raw=None,
                 rings=(), kappa2=False):
    """covering certificates for rows without a unimodular certificate (R38; search order of R39).  Each round, for
    every row not yet covered:
      - without a rational realisation: a rational search (residual over 2^31 - 1 only), plain, from the stored qubit
        certificate at multiplicity 1 (random start without one), and the same in ring mode (kappa 2) for every ring
        strategy in `rings`; at multiplicity 2 from the stored qubit certificate at multiplicity 2 if the row has no
        qubit certificate at multiplicity 1, and (with kappa2) from two copies of the one at multiplicity 1;
      - with one: a rational search that must also be right mod every uncovered prime (at most four), started from
        the row's rational realisation with the fewest uncovered and bad primes.  (R38 started these searches from the
        qubit certificate, one prime at a time; starting from a rational realisation is what closed most rows in R39.)
    Every reported matrix has its bad primes recomputed exactly.  One without bad primes is a unimodular realisation;
    it is certified as in --check and added to the unimodular file.  A row is covered when the bad-prime sets of its
    rational realisations intersect inside the primes of its stored qubit (p = 2) and qutrit (p = 3) certificates;
    a smallest such set of rational realisations is stored, with the GF(p) certificates it uses."""
    reps = catalogue()
    qcerts, tcerts = qubit_certificates(reps), qutrit_certificates(reps)
    uni = {c['row']: c for c in load(unimodular_path)} if os.path.exists(unimodular_path) else {}
    have = load_cover(cover_out) if os.path.exists(cover_out) else []
    have_rows = {c['row'] for c in have}
    todo = [k for k in rows if k not in uni and k not in have_rows]
    newuni = {}
    gf = {}
    for k in todo:
        r = reps[k]
        if k in qcerts:
            lam, sizes, adj, src = min(qcerts[k], key=lambda c: c[0])
            gf.setdefault(k, []).append(dict(row=k, kind='GF', prime=2, kappa=lam, sizes=sizes, W=adj, bad=[],
                                             detail=f'stored qubit certificate {src}'))
        if k in tcerts:
            lam, sizes, W3 = min(tcerts[k], key=lambda c: c[0])
            gf.setdefault(k, []).append(dict(row=k, kind='GF', prime=3, kappa=lam, sizes=sizes, W=W3, bad=[],
                                             detail='stored qutrit certificate'))
    qpieces = {k: [] for k in todo}
    tmp = tempfile.mkdtemp(prefix='s32c_')

    def need(k):
        if k in newuni:
            return []
        if not qpieces[k]:
            return None
        common = set.intersection(*[set(c['bad']) for c in qpieces[k]])
        return sorted(common - {c['prime'] for c in gf.get(k, [])})
    try:
        exe = build_kernel(tmp)
        for rnd in range(rounds + 1):
            # the kernel takes one ring and one active set of primes per call
            calls = {}
            for k in todo:
                n = need(k)
                if n == []:
                    continue
                r = reps[k]
                if n is None:
                    lam = min(c[0] for c in qcerts[k]) if k in qcerts else None
                    start = min(qcerts[k], key=lambda c: c[0])[2] if lam == 1 else None
                    for ring in ['0'] + ['%d,%d' % RINGS[name] for name in rings]:
                        calls.setdefault((ring, str(BIG31)), []).append(job_line(k, 1, singles(r), r, start))
                    start2 = (next(c[2] for c in qcerts[k] if c[0] == 2) if lam == 2 else
                              doubled(start, singles(r)) if kappa2 and lam == 1 else None)
                    if start2 is not None:
                        calls[('0', str(BIG31))].append(job_line(k, 2, [2 * x for x in singles(r)], r, start2))
                else:
                    c = min(qpieces[k], key=lambda c: (len(set(c['bad']) & set(n)), len(c['bad']), c['kappa']))
                    ps = ','.join(map(str, [BIG31] + n[:4]))
                    calls.setdefault(('0', ps), []).append(job_line(k, c['kappa'], c['sizes'], r, c['W']))
            if not calls:
                break
            for (ring, ps), lines in calls.items():
                # sweeps per anneal: 400 for the plain rational searches, 200 otherwise (as in the R39 campaign)
                sweeps = '400' if ring == '0' and ps == str(BIG31) else '200'
                text = kernel(exe, lines, seconds, seed + 7 * rnd, ring, sweeps, '0.5', ps, '0')
                if raw:
                    with open(raw, 'a') as f:
                        f.write(f'# cover round {rnd} ring {ring} primes {ps} seed {seed + 7 * rnd} seconds {seconds}\n'
                                + text)
                for k, kap, W, _secs in parse_found(text):
                    r = reps[k]; sizes = [kap * x for x in singles(r)]
                    b = bad_primes(r, kap, sizes, W)
                    if b is None:
                        continue
                    detail = f'cover round {rnd}, ring {ring}, active primes {ps}, seed {seed + 7 * rnd}'
                    if not b:
                        if all(certify(r, kap, sizes, W).values()) and (k not in newuni or newuni[k]['kappa'] > kap):
                            newuni[k] = dict(row=k, rep=r, kappa=kap, sizes=sizes, W=W, method='cover', detail=detail)
                        continue
                    old = need(k)
                    if old is not None and set(b) >= set(old):
                        continue                   # no progress
                    qpieces[k].append(dict(row=k, kind='Q', prime=0, kappa=kap, sizes=sizes, W=W, bad=b, detail=detail))
            log(f'  cover round {rnd}: {sum(len(v) for v in calls.values())} jobs in {len(calls)} calls; covered so '
                f'far {sum(1 for k in todo if need(k) == [])} of {len(todo)} ({len(newuni)} by unimodular realisations)')
        new = []
        for k in todo:
            if k not in newuni and need(k) == []:
                gfp = [c['prime'] for c in gf.get(k, [])]
                keep = minimal_cover(qpieces[k], gfp)
                common = set.intersection(*[set(c['bad']) for c in keep])
                new += keep + [c for c in gf.get(k, []) if c['prime'] in common]
        save_cover(have + new, cover_out)
        if newuni:
            save(list(uni.values()) + list(newuni.values()), unimodular_path)
        log(f'  covering certificates stored for {len({c["row"] for c in have + new})} rows; '
            f'{len(newuni)} new unimodular certificates')
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    return have + new


# ----------------------------------------------------------------------------------------------- self-test
def self_test(trials=300, seed=32, log=print):
    """(1) party-local integer basis changes and integer local complementations keep the Smith form of every cut
    block (the neutral moves of the kernel); (2) the ring expansion G M(v) is symmetric and G unimodular; (3) the
    Z_d entropy formula against state vectors on random small integer graph states."""
    rng = random.Random(seed)

    def snf(B):
        d = sorted(diagonal_form(B)) if B and B[0] else []
        changed = True
        while changed:
            changed = False
            for i in range(len(d)):
                for j in range(i + 1, len(d)):
                    g = gcd(d[i], d[j]); l = d[i] * d[j] // g
                    if (g, l) != (d[i], d[j]):
                        d[i], d[j] = g, l; changed = True
            d.sort()
        return tuple(d)
    sizes = [2, 2, 3, 1, 2, 3]; N = sum(sizes)
    owner = [q for q, c in enumerate(sizes) for _ in range(c)]
    bad = 0
    for _ in range(trials):
        W = [[0] * N for _ in range(N)]
        for i in range(N):
            for j in range(i + 1, N):
                if owner[i] != owner[j]:
                    W[i][j] = W[j][i] = rng.choice([-1, 0, 0, 1, 2])
        prof = [snf(B) for _m, B in cut_blocks(W, sizes)]
        c = rng.randrange(N); a = rng.choice([-1, 1, 2])
        W2 = [[0 if owner[u] == owner[x] else W[u][x] - a * W[u][c] * W[x][c] for x in range(N)] for u in range(N)]
        q = rng.choice([q for q in range(6) if sizes[q] >= 2])
        u1, u2 = rng.sample([i for i in range(N) if owner[i] == q], 2); e = rng.choice([-1, 1, 3])
        W3 = [row[:] for row in W]
        for x in range(N):
            if owner[x] != q:
                W3[u2][x] += e * W[u1][x]; W3[x][u2] = W3[u2][x]
        bad += [snf(B) for _m, B in cut_blocks(W2, sizes)] != prof or [snf(B) for _m, B in cut_blocks(W3, sizes)] != prof
    ok1 = bad == 0
    log(f'neutral moves keep every cut Smith form: {trials - bad}/{trials} random trials')
    ok2 = True
    for s, t in RINGS.values():
        G = np.array([[0, 1], [1, s]])
        for x, y in itertools.product(range(-2, 3), repeat=2):
            M = np.array([[x, t * y], [y, x + s * y]])         # multiplication by x + y theta in the basis {1, theta}
            ok2 &= np.array_equal(G @ M, (G @ M).T)
        ok2 &= round(abs(np.linalg.det(G))) == 1
    log(f'ring expansion: G M(v) symmetric and G unimodular for {len(RINGS)} rings: {ok2}')
    ok3, n3 = True, 0
    for d in (4, 6):
        for _ in range(4):
            sz = [1, 1, 1, 1, 1, 1] if d == 4 else [1, 1, 1, 1, 1, 0]
            Nn = sum(sz); own = [q for q, c in enumerate(sz) for _ in range(c)]
            W = [[0] * Nn for _ in range(Nn)]
            for i in range(Nn):
                for j in range(i + 1, Nn):
                    W[i][j] = W[j][i] = rng.randrange(d)
            # entropies from the formula: log_d |image of the block over Z_d|
            psi_ok = True
            x = np.array(list(itertools.product(range(d), repeat=Nn)), dtype=np.int64)
            phase = np.einsum('ki,ij,kj->k', x, np.triu(np.array(W) % d, 1), x) % d
            psi = (np.exp(2j * np.pi * phase / d) / d ** (Nn / 2)).reshape([d] * Nn)
            for m in range(1, 32):
                X = [v for v in range(Nn) if own[v] < 5 and m >> own[v] & 1]
                Y = [v for v in range(Nn) if v not in X]
                if not X or not Y:
                    continue
                B = [[W[a][b] for b in Y] for a in X]
                img = {tuple(sum(B[i][j] * y[j] for j in range(len(Y))) % d for i in range(len(X)))
                       for y in itertools.product(range(d), repeat=len(Y))}
                s = np.linalg.svd(np.transpose(psi, X + Y).reshape(d ** len(X), d ** len(Y)), compute_uv=False)
                p = s[s > 1e-8] ** 2
                psi_ok &= len(p) == len(img) and np.abs(p - 1.0 / len(img)).max() < 1e-8
            ok3 &= psi_ok; n3 += 1
    log(f'Z_d cut formula |W_X Z_d^Y| vs state vectors (d = 4, 6; {n3} random graph states, every cut): {ok3}')
    return ok1 and ok2 and ok3


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--check', action='store_true')
    ap.add_argument('--full', action='store_true', help='with --check: M2, GF(p) ranks, Z_d state vectors')
    ap.add_argument('--self-test', action='store_true')
    ap.add_argument('--search', action='store_true')
    ap.add_argument('--rows', default=None, help='comma-separated catalogue rows (default: all 641)')
    ap.add_argument('--strategies', default='plain,neutral,ring-eisenstein,qutrit-start,gf3-search-start,neutral',
                    help='comma-separated, run in order on the rows still unsolved: ' + ', '.join(STRATEGIES))
    ap.add_argument('--seconds', type=float, default=10.0)
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--out', default=OUT_DEFAULT)
    ap.add_argument('--resume', action='store_true')
    ap.add_argument('--raw', default=None, help='append the raw kernel output here')
    ap.add_argument('--cover', action='store_true', help='covering certificates for the rows without a unimodular one')
    ap.add_argument('--cover-out', default=COVER_DEFAULT)
    ap.add_argument('--rounds', type=int, default=3, help='with --cover: rounds after the first')
    ap.add_argument('--rings', default='', help='with --cover: ring strategies (comma-separated, e.g. '
                    'ring-eisenstein,ring-gauss) also tried on rows without a rational realisation')
    ap.add_argument('--kappa2', action='store_true', help='with --cover: rows without a rational realisation also get '
                    'a search at multiplicity 2 from two copies of their qubit certificate')
    ap.add_argument('--check-cover', action='store_true')
    a = ap.parse_args()
    if a.self_test:
        ok = self_test()
        print('ALL OK' if ok else 'FAILURE')
    if a.search:
        rows = [int(x) for x in a.rows.split(',')] if a.rows else list(range(len(catalogue())))
        certs = search(rows, a.strategies.split(','), a.seconds, a.seed, a.out, a.resume, raw=a.raw)
        print(f'certified rows: {len(certs)} of {len(rows)}')
    if a.cover:
        rows = [int(x) for x in a.rows.split(',')] if a.rows else list(range(len(catalogue())))
        cover_search(rows, a.seconds, a.seed, a.cover_out, unimodular_path=a.out, rounds=a.rounds, raw=a.raw,
                     rings=[x for x in a.rings.split(',') if x], kappa2=a.kappa2)
    if a.check_cover:
        s = check_cover(a.cover_out, a.out)
        print(json.dumps({k: v for k, v in s.items() if k != 'failures'}))
        for f in s['failures'][:20]:
            print('FAILED', f)
        print('ALL OK' if s['pieces_ok'] == s['pieces'] and s['rows_covered'] == s['rows'] and not s['overlap_unimodular']
              else 'FAILURE')
    if a.check:
        s = check(a.out, full=a.full)
        print(json.dumps({k: v for k, v in s.items() if k != 'failures'}))
        for k, res in s['failures'][:20]:
            print('FAILED row', k, res)
        print('ALL OK' if s['ok'] == s['stored'] and s['rep_ok'] and s['distinct_rows'] == s['stored'] else 'FAILURE')


if __name__ == '__main__':
    main()
