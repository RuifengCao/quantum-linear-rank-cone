"""Certificate checkers for A2 (R27): stabilizer realisations and stabilizer-cone exclusions.

Realisation certificate  (lam, sizes, adjacency A)  for a ray r:
    A is an N x N symmetric GF(2) matrix with zero diagonal, N = sum(sizes), party p
    owning sizes[p] consecutive qubits (parties A..E, purifier F last); the graph state
    |A> has entropy vector lam * r, i.e. rank_GF2 A[X, X^c] = lam * r_X for all 31 X.
    (Hein-Eisert-Briegel, quant-ph/0307130, Prop. 3 / Eq. (43).)

Exclusion certificate  (X, Y, y)  for a ray r:
    y is a nonnegative rational combination of Shannon elemental inequalities on the
    seven elements A..F,Z whose Z-terms are eliminated by the common-information
    equalities g(XZ)=h(X), g(YZ)=h(Y), g(Z)=h(X)+h(Y)-h(XY); the resulting inequality
    F(h) >= 0 holds for every linear polymatroid on A..F, hence F'(S) = F(h_norm(S)) >= 0
    for every stabilizer state (see README, A2); the check is F'(r) < 0 in exact arithmetic.
"""
from fractions import Fraction
import numpy as np

SING = [0, 1, 3, 7, 15]
FULL6 = 63
ZB = 1 << 6


# ---------------------------------------------------------------- realisations
def gf2_rank(M):
    M = (np.array(M, dtype=np.uint8) & 1).copy()
    rows, cols = M.shape
    r = 0
    for c in range(cols):
        piv = next((i for i in range(r, rows) if M[i, c]), None)
        if piv is None:
            continue
        if piv != r:
            M[[r, piv]] = M[[piv, r]]
        for i in range(rows):
            if i != r and M[i, c]:
                M[i] ^= M[r]
        r += 1
        if r == rows:
            break
    return r


def rows_to_adjacency(rows64, N):
    A = np.zeros((N, N), dtype=np.uint8)
    for i in range(N):
        v = int(rows64[i])
        for j in range(N):
            A[i, j] = (v >> j) & 1
    return A


def cut_vector(A, sizes):
    party = np.repeat(np.arange(6), sizes)
    out = np.zeros(31, dtype=np.int64)
    for m in range(1, 32):
        inX = np.array([p < 5 and ((m >> p) & 1) == 1 for p in party])
        out[m - 1] = gf2_rank(A[np.ix_(inX, ~inX)]) if inX.any() and (~inX).any() else 0
    return out


def verify_realisation(r, lam, sizes, rows64):
    sizes = [int(s) for s in sizes]
    N = sum(sizes)
    A = rows_to_adjacency(rows64, N)
    if not np.array_equal(A, A.T) or A.diagonal().any():
        return False
    return bool(np.array_equal(cut_vector(A, sizes), int(lam) * np.asarray(r, dtype=np.int64)))


def gfp_rank(M, p):
    """(R29) rank over GF(p), p prime: plain-Python elimination on lists of ints (independent of
    scripts/s21_qudit_search.py, which uses numpy row operations)."""
    rows = [[int(x) % p for x in row] for row in M]
    rank, ncols = 0, (len(rows[0]) if rows else 0)
    for c in range(ncols):
        piv = next((i for i in range(rank, len(rows)) if rows[i][c]), None)
        if piv is None:
            continue
        rows[rank], rows[piv] = rows[piv], rows[rank]
        inv = pow(rows[rank][c], p - 2, p)
        rows[rank] = [(x * inv) % p for x in rows[rank]]
        for i in range(len(rows)):
            if i != rank and rows[i][c]:
                f = rows[i][c]
                rows[i] = [(x - f * y) % p for x, y in zip(rows[i], rows[rank])]
        rank += 1
    return rank


def verify_realisation_gfp(r, lam, sizes, W, p):
    """(R29) qudit certificate: W is an N x N symmetric matrix over GF(p) with zero diagonal, party q
    owning sizes[q] consecutive qudits; the weighted graph state has entropy (base p)
    S(X) = rank_GF(p) W[X, X^c], which must equal lam * r_X for all 31 X."""
    sizes = [int(s) for s in sizes]
    N = sum(sizes)
    W = [[int(W[i][j]) for j in range(N)] for i in range(N)]
    if any(W[i][j] != W[j][i] or not 0 <= W[i][j] < p for i in range(N) for j in range(N)) or any(W[i][i] for i in range(N)):
        return False
    party = [q for q in range(6) for _ in range(sizes[q])]
    for m in range(1, 32):
        X = [i for i in range(N) if party[i] < 5 and (m >> party[i]) & 1]
        Y = [i for i in range(N) if not (party[i] < 5 and (m >> party[i]) & 1)]
        rk = gfp_rank([[W[i][j] for j in Y] for i in X], p) if X and Y else 0
        if rk != int(lam) * int(r[m - 1]):
            return False
    return True


def verify_realisation_n(r, n, lam, sizes, W, p):
    """(R32) n-party version of verify_realisation_gfp: parties 0..n-1 are visible, party n is the purifier
    (sizes has n + 1 entries, the purifier's qudits last); W is an N x N symmetric matrix over GF(p) with zero
    diagonal; the weighted graph state must have rank_GF(p) W[X, X^c] == lam * r[m - 1] for every nonempty
    mask m over the n visible parties (r has 2^n - 1 entries).  p = 2 is allowed."""
    sizes = [int(s) for s in sizes]
    if len(sizes) != n + 1 or len(r) != (1 << n) - 1:
        return False
    N = sum(sizes)
    W = [[int(W[i][j]) for j in range(N)] for i in range(N)]
    if any(W[i][j] != W[j][i] or not 0 <= W[i][j] < p for i in range(N) for j in range(N)) or any(W[i][i] for i in range(N)):
        return False
    party = [q for q in range(n + 1) for _ in range(sizes[q])]
    for m in range(1, 1 << n):
        X = [i for i in range(N) if party[i] < n and (m >> party[i]) & 1]
        Y = [i for i in range(N) if not (party[i] < n and (m >> party[i]) & 1)]
        rk = gfp_rank([[W[i][j] for j in Y] for i in X], p) if X and Y else 0
        if rk != int(lam) * int(r[m - 1]):
            return False
    return True


# ---------------------------------------------------------------- exclusions
def entropy6(r):
    S = [0] * 64
    for m in range(1, 64):
        if m & 32:
            c = FULL6 ^ m
            S[m] = 0 if c == 0 else int(r[c - 1])
        else:
            S[m] = int(r[m - 1])
    return S


def h_norm(r):
    S = entropy6(r)
    sing = [S[1 << p] for p in range(6)]
    return [0] + [S[m] + sum(sing[p] for p in range(6) if m >> p & 1) for m in range(1, 64)]


def elemental(n=7):
    full = (1 << n) - 1
    rows = []
    for i in range(n):
        d = {full: 1}
        if full ^ (1 << i):
            d[full ^ (1 << i)] = -1
        rows.append(d)
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


_ELEM = None


def verify_exclusion(r, X, Y, y, F_S=None):
    """y: dict {elemental index: Fraction-like}.  Returns (ok, value of F'(r))."""
    global _ELEM
    if _ELEM is None:
        _ELEM = elemental()
    y = {int(k): Fraction(v) for k, v in y.items()}
    if any(v < 0 for v in y.values()):
        return False, None
    c = {}
    for k, v in y.items():
        for m, e in _ELEM[k].items():
            c[m] = c.get(m, 0) + v * e
    F = {}
    for m, v in c.items():
        if v == 0:
            continue
        if not (m & ZB):
            F[m] = F.get(m, 0) + v
        elif m == (X | ZB):
            F[X] = F.get(X, 0) + v
        elif m == (Y | ZB):
            F[Y] = F.get(Y, 0) + v
        elif m == ZB:
            F[X] = F.get(X, 0) + v; F[Y] = F.get(Y, 0) + v; F[X | Y] = F.get(X | Y, 0) - v
        else:
            return False, None           # a Z-term survived: not a valid derivation
    h = h_norm(r)
    val = sum(v * h[m] for m, v in F.items())
    if F_S is not None:
        # the stored S-form must be a positive multiple of F(h_norm(.)) as a linear form
        coef = []
        for i in range(31):
            e = [0] * 31; e[i] = 1
            he = h_norm(e)
            coef.append(sum(v * he[m] for m, v in F.items()))
        ratios = {Fraction(a) / b for a, b in zip(coef, F_S) if b != 0}
        if any((a != 0) != (b != 0) for a, b in zip(coef, F_S)) or len(ratios) != 1 or next(iter(ratios)) <= 0:
            return False, val
    return val < 0, val


# ---------------------------------------------------------------- R29: Ingleton on the extension
def ingleton_row(a, b, c, d, K):
    """Conditional Ingleton instance on singletons a, b, c, d (distinct, 0..6) given the set K
    (bitmask disjoint from them):  r(ab|K) + r(ac|K) + r(ad|K) + r(bc|K) + r(bd|K)
    >= r(a|K) + r(b|K) + r(cd|K) + r(abc|K) + r(abd|K), written with r(S|K) -> g(S u K)
    (the g(K) terms cancel).  Valid for every subspace arrangement over any field
    (Ingleton 1971; contraction by K keeps the arrangement linear).  Rebuilt here from the
    formula, independently of scripts/s17_ci_test.py."""
    A, B, C, D = 1 << a, 1 << b, 1 << c, 1 << d
    row = {}
    for s, sign in ((A | B, 1), (A | C, 1), (A | D, 1), (B | C, 1), (B | D, 1),
                    (A, -1), (B, -1), (C | D, -1), (A | B | C, -1), (A | B | D, -1)):
        row[s | K] = row.get(s | K, 0) + sign
    return row


def verify_exclusion_ext(r, X, Y, terms, F_S=None):
    """Exclusion certificate with typed terms (R29):  terms = [(descriptor, multiplier), ...],
    descriptor ('E', k) = k-th Shannon elemental inequality on A..F,Z (same order as elemental()),
    ('I', a, b, c, d, K) = ingleton_row(a, b, c, d, K).  Multipliers must be >= 0.
    Returns (ok, value of F'(r)) exactly as verify_exclusion."""
    global _ELEM
    if _ELEM is None:
        _ELEM = elemental()
    c = {}
    for desc, v in terms:
        v = Fraction(v)
        if v < 0:
            return False, None
        if desc[0] == 'E':
            row = _ELEM[int(desc[1])]
        elif desc[0] == 'I':
            a, b, cc, d, K = (int(t) for t in desc[1:])
            if len({a, b, cc, d}) != 4 or not all(0 <= t < 7 for t in (a, b, cc, d)) \
                    or K < 0 or K >= 128 or K & ((1 << a) | (1 << b) | (1 << cc) | (1 << d)):
                return False, None
            row = ingleton_row(a, b, cc, d, K)
        else:
            return False, None
        for m, e in row.items():
            c[m] = c.get(m, 0) + v * e
    F = {}
    for m, v in c.items():
        if v == 0:
            continue
        if not (m & ZB):
            F[m] = F.get(m, 0) + v
        elif m == (X | ZB):
            F[X] = F.get(X, 0) + v
        elif m == (Y | ZB):
            F[Y] = F.get(Y, 0) + v
        elif m == ZB:
            F[X] = F.get(X, 0) + v; F[Y] = F.get(Y, 0) + v; F[X | Y] = F.get(X | Y, 0) - v
        else:
            return False, None
    h = h_norm(r)
    val = sum(v * h[m] for m, v in F.items())
    if F_S is not None:
        coef = []
        for i in range(31):
            e = [0] * 31; e[i] = 1
            he = h_norm(e)
            coef.append(sum(v * he[m] for m, v in F.items()))
        ratios = {Fraction(a) / b for a, b in zip(coef, F_S) if b != 0}
        if any((a != 0) != (b != 0) for a, b in zip(coef, F_S)) or len(ratios) != 1 or next(iter(ratios)) <= 0:
            return False, val
    return val < 0, val


def verify_single_ci_feasible(r, pairs, den, num):
    """(R29) Exact proof that the ray r passes EVERY single-common-information test: for every
    unordered pair {X, Y} of distinct nonempty party sets there is an extension g of h_norm(r) to
    A..F,Z satisfying g(XZ) = h(X), g(YZ) = h(Y), g(Z) = h(X) + h(Y) - h(X u Y) and every Shannon
    elemental inequality on the seven elements.  Stored witnesses (pairs[k] = (X, Y), values
    g(m u Z) = num[k][m] / den[k] for m = 0..63) cover the pairs with I(X;Y) > 0 and neither set
    containing the other; for the remaining pairs the witness is explicit (Z a copy of the smaller
    set, or Z = 0) and is checked here too.  Integer arithmetic throughout.
    Returns (ok, number of pairs checked, number of stored witnesses used)."""
    h = h_norm(r)
    E = elemental()
    M = np.zeros((len(E), 127), dtype=np.int64)
    for k, row in enumerate(E):
        for m, c in row.items():
            M[k, m - 1] = c
    stored = {(int(X), int(Y)): (int(d), [int(t) for t in nm]) for (X, Y), d, nm in zip(pairs, den, num)}
    cols, used = [], 0
    for X in range(1, 64):
        for Y in range(X + 1, 64):
            I = h[X] + h[Y] - h[X | Y]
            if (X & Y) == X or (X & Y) == Y:          # comparable: Z = copy of the smaller set
                Sm = X if (X & Y) == X else Y
                d, gz = 1, [h[m | Sm] if (m | Sm) else 0 for m in range(64)]
            elif I == 0:                               # independent: Z = 0
                d, gz = 1, [h[m] if m else 0 for m in range(64)]
            else:
                if (X, Y) not in stored:
                    return False, None, used
                d, gz = stored[(X, Y)]; used += 1
            if d <= 0 or gz[X] != d * h[X] or gz[Y] != d * h[Y] or gz[0] != d * I:
                return False, None, used
            cols.append([d * h[m] for m in range(1, 64)] + list(gz))
    G = np.array(cols, dtype=np.int64).T                # 127 x pairs, scaled by each pair's denominator
    ok = bool((M @ G >= 0).all())
    return ok, G.shape[1], used


# ---------------------------------------------------------------- compact storage (R28)
def pack_exclusions(certs, maxk=None):
    """certs: list of dicts with ray, X, Y, F_S, y ({k: 'a/b'}).  Returns dict of arrays for np.savez."""
    from fractions import Fraction
    n = len(certs)
    K = maxk or max(len(c['y']) for c in certs)
    yk = np.full((n, K), -1, dtype=np.int16)
    yn = np.zeros((n, K), dtype=np.int64)
    yd = np.ones((n, K), dtype=np.int64)
    for i, c in enumerate(certs):
        for j, (k, v) in enumerate(sorted(c['y'].items(), key=lambda t: int(t[0]))):
            f = Fraction(v)
            yk[i, j] = int(k); yn[i, j] = f.numerator; yd[i, j] = f.denominator
    return {'ray': np.array([c['ray'] for c in certs], dtype=np.int16),
            'X': np.array([c['X'] for c in certs], dtype=np.int8),
            'Y': np.array([c['Y'] for c in certs], dtype=np.int8),
            'F_S': np.array([c['F_S'] for c in certs], dtype=np.int8),
            'y_index': yk, 'y_num': yn, 'y_den': yd}


def unpack_exclusion(Z, i):
    from fractions import Fraction
    y = {int(k): Fraction(int(a), int(b)) for k, a, b in zip(Z['y_index'][i], Z['y_num'][i], Z['y_den'][i]) if k >= 0}
    return (Z['ray'][i].astype(np.int64), int(Z['X'][i]), int(Z['Y'][i]), y, Z['F_S'][i].astype(np.int64).tolist())
