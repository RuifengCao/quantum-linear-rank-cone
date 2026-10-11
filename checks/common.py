"""Independent helper functions for the paper's numerical claims; they do not use the repository code.

Conventions: a vector on the nonempty subsets of k labelled elements is a dict {frozenset: value} or a list in
mask order (mask m = sum of 2^i over elements i, i = 0..k-1; entry m-1 holds the value of mask m).
"""
import itertools, os
from fractions import Fraction

# path of a checkout of https://github.com/RuifengCao/quantum-linear-rank-cone (release R37); by default the parent
# directory of this file, i.e. the repository itself when these scripts are run from its directory checks/
REPO = os.environ.get('EPR1_REPO', os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def masks(k):
    return range(1, 1 << k)


def bits(m):
    return [i for i in range(m.bit_length()) if m >> i & 1]


def rank_mod_p(rows, p):
    """rank over GF(p) of a list of integer row vectors."""
    M = [[x % p for x in r] for r in rows if any(x % p for x in r)]
    if not M:
        return 0
    ncol = len(M[0])
    r = 0
    for c in range(ncol):
        piv = None
        for i in range(r, len(M)):
            if M[i][c] % p:
                piv = i
                break
        if piv is None:
            continue
        M[r], M[piv] = M[piv], M[r]
        inv = pow(M[r][c], p - 2, p)
        M[r] = [(x * inv) % p for x in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c]:
                f = M[i][c]
                M[i] = [(a - f * b) % p for a, b in zip(M[i], M[r])]
        r += 1
        if r == len(M):
            break
    return r


def rank_Q(rows):
    """rank over Q (exact, Fractions)."""
    M = [[Fraction(x) for x in r] for r in rows if any(r)]
    if not M:
        return 0
    ncol = len(M[0])
    r = 0
    for c in range(ncol):
        piv = None
        for i in range(r, len(M)):
            if M[i][c] != 0:
                piv = i
                break
        if piv is None:
            continue
        M[r], M[piv] = M[piv], M[r]
        inv = 1 / M[r][c]
        M[r] = [x * inv for x in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c] != 0:
                f = M[i][c]
                M[i] = [a - f * b for a, b in zip(M[i], M[r])]
        r += 1
        if r == len(M):
            break
    return r


def graph_state_entropy(adj, owner, nparties, p=2):
    """Entropy vector (units of log p) of the (weighted) graph state with symmetric adjacency matrix adj over
    GF(p); owner[v] = party of vertex v in 0..nparties-1. Returns the list over masks 1..2^nparties - 1:
    S(X) = rank_GF(p) adj[X-vertices, complement-vertices]."""
    nv = len(adj)
    out = []
    for m in masks(nparties):
        inX = [v for v in range(nv) if m >> owner[v] & 1]
        out_ = [v for v in range(nv) if not (m >> owner[v] & 1)]
        rows = [[adj[a][b] for b in out_] for a in inX]
        out.append(rank_mod_p(rows, p) if inX and out_ else 0)
    return out


def column_rank_function(cols, owner, k, p):
    """rank function (list over masks of k elements) of the column groups of a matrix over GF(p);
    cols = list of column vectors, owner[j] = element owning column j."""
    out = []
    for m in masks(k):
        sel = [cols[j] for j in range(len(cols)) if m >> owner[j] & 1]
        out.append(rank_mod_p(sel, p) if sel else 0)
    return out


def set_name(m, letters):
    return ''.join(letters[i] for i in bits(m))
