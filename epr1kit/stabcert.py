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
