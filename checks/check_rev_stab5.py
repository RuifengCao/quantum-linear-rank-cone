"""Revised paper, Sec. II E and Supplementary Note 3.2: the 641 known orbits of extreme rays of Stab_5 lie in the
stabilizer cone of every local dimension d, and 119 of the 747 single-CI inequality classes are facets of
Stab^(d)_5 for every d.

Data read (not the repository code): data/stab5_extreme_reps.npy (the 641 orbit representatives),
data/s32_stab5_unimodular.npz (566 unimodular realizations), data/s32_stab5_cover.npz (199 pieces of covering
certificates for the other 75 rows), data/s17_ineq_classes.npy (the 747 inequality classes).

A realization of k*r is an integer symmetric matrix W on N = k*(r(A) + ... + r(F)) vertices, system q owning k*r(q)
of them (r(F) = r(ABCDE) by purity), with zero entries inside each system. For a prime p, the graph state of Z_p
qudits with weights W mod p has S(X) = rank_GF(p) W[Q_X, Q_rest] (units of log p). A prime p is *bad* for W if, on
some cut X, the rank of the cut block modulo p is smaller than its rational rank; W is *unimodular* if no prime is
bad. Bad primes are computed here without Smith forms: fraction-free elimination gives the rational rank t of a cut
block and a nonzero t x t minor m; the t-th determinantal divisor divides m, so a bad prime divides m, and each prime
factor q of m is tested by a rank computation modulo q.

Checks:
  (1) structure of every stored matrix (symmetric, zero inside each system, sizes k*r(q));
  (2) unimodular file: rational rank k*r(X) on all 31 cuts and no bad prime (566 rows; k = 1, 2, 3 for 515, 50, 1);
  (3) cover file: every rational piece has rational ranks k*r(X) and exactly the stored bad primes; every GF(p)
      piece has ranks k*r(X) modulo p; for each of the 75 rows the bad sets of its (one or more) rational pieces meet
      in a set of primes covered by its GF pieces (so the row lies in Stab^(p)_5 for every prime p, hence, by Corollary 2, in
      Stab^(d)_5 for every d); the 566 + 75 rows are disjoint and cover all 641;
  (4) every row has an explicit qubit realization and an explicit qutrit realization among these matrices;
  (5) facet count: the 760 entropy vectors of six-vertex 0/1 graph states over GF(2), which are also cut-rank vectors
      over every odd prime field (check_rev_sec4.py), and the S_6 images of the 641 rows all lie in Stab^(d)_5 for
      every d; no class is violated by them, and exactly 119 of the 747 classes have 30 linearly independent of them
      on their hyperplane (rank modulo 1,000,003, a lower bound on the rational rank); the class of (star) is one
      of them.
Runs in a few minutes.
"""
import itertools, os, sys
from collections import Counter
import numpy as np
import sympy

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import REPO, rank_mod_p, rank_Q

failures = []


def ok(cond, msg):
    print(('OK   ' if cond else 'FAIL ') + msg)
    if not cond:
        failures.append(msg)


FIVE_MASKS = list(range(1, 32))          # subsets of {A,...,E}; A is bit 0


def bareiss_rank_minor(B):
    """rational rank t of the integer matrix B (list of rows) and |m| for a nonzero t x t minor m (1 if t = 0),
    by fraction-free elimination with full pivoting."""
    M = [list(map(int, r)) for r in B]
    nr = len(M)
    nc = len(M[0]) if nr else 0
    prev, t = 1, 0
    for k in range(min(nr, nc)):
        piv = None
        for i in range(k, nr):
            for j in range(k, nc):
                if M[i][j]:
                    piv = (i, j)
                    break
            if piv:
                break
        if piv is None:
            break
        i, j = piv
        M[k], M[i] = M[i], M[k]
        for row in M:
            row[k], row[j] = row[j], row[k]
        pk = M[k][k]
        for i in range(k + 1, nr):
            for j in range(k + 1, nc):
                M[i][j] = (M[i][j] * pk - M[i][k] * M[k][j]) // prev
            M[i][k] = 0
        prev = pk
        t += 1
    return t, abs(prev) if t else 1


def blocks(W, sizes):
    """the 31 cut blocks of W, in mask order of the subsets X of {A,...,E}"""
    owner = [q for q, s in enumerate(sizes) for _ in range(s)]
    out = []
    for X in FIVE_MASKS:
        I = [v for v in range(len(owner)) if X >> owner[v] & 1]
        J = [v for v in range(len(owner)) if not X >> owner[v] & 1]
        out.append([[int(W[a][b]) for b in J] for a in I] if I and J else [])
    return out


def structure_ok(W, sizes):
    n = len(W)
    owner = [q for q, s in enumerate(sizes) for _ in range(s)]
    if len(owner) != n:
        return False
    for a in range(n):
        for b in range(n):
            if W[a][b] != W[b][a]:
                return False
            if owner[a] == owner[b] and W[a][b]:
                return False
    return True


def bad_primes_and_ranks(W, sizes):
    """rational ranks on the 31 cuts and the set of bad primes"""
    ranks, bad = [], set()
    for B in blocks(W, sizes):
        if not B:
            ranks.append(0)
            continue
        t, m = bareiss_rank_minor(B)
        ranks.append(t)
        for q in sympy.factorint(m):
            if rank_mod_p(B, q) < t:
                bad.add(int(q))
    return ranks, bad


def ranks_mod(W, sizes, p):
    return [rank_mod_p(B, p) if B else 0 for B in blocks(W, sizes)]


R = np.load(os.path.join(REPO, 'data', 'stab5_extreme_reps.npy')).astype(np.int64)
ok(R.shape == (641, 31), 'stab5_extreme_reps.npy: %s rows' % (R.shape,))


def expected_sizes(r, k):
    r = [int(x) for x in r]
    return [k * r[(1 << q) - 1] for q in range(5)] + [k * r[30]]


# ---------------------------------------------------------------------------------------------- (1), (2) unimodular
U = np.load(os.path.join(REPO, 'data', 's32_stab5_unimodular.npz'), allow_pickle=True)
uni_rows, uni_kappa = [], Counter()
bad_struct = bad_rank = bad_uni = 0
maxN = 0
for e in range(len(U['row'])):
    row, k, n = int(U['row'][e]), int(U['kappa'][e]), int(U['N'][e])
    sizes = [int(s) for s in U['sizes'][e]]
    W = U['W'][e][:n, :n].astype(int).tolist()
    r = R[row]
    if not (list(U['rep'][e]) == list(r) and sizes == expected_sizes(r, k) and structure_ok(W, sizes)
            and not U['W'][e][n:, :].any() and not U['W'][e][:, n:].any()):
        bad_struct += 1
        continue
    ranks, bad = bad_primes_and_ranks(W, sizes)
    if ranks != [k * int(x) for x in r]:
        bad_rank += 1
    if bad:
        bad_uni += 1
    uni_rows.append(row)
    uni_kappa[k] += 1
    maxN = max(maxN, n)
ok(bad_struct == 0 and bad_rank == 0 and bad_uni == 0,
   'unimodular file: %d realizations; structure failures %d, rank failures %d, with a bad prime %d'
   % (len(U['row']), bad_struct, bad_rank, bad_uni))
ok(len(set(uni_rows)) == len(uni_rows) == 566 and dict(uni_kappa) == {1: 515, 2: 50, 3: 1} and maxN <= 40,
   'unimodular file: %d distinct rows, multiplicities %s, at most %d vertices (expect 566, {1: 515, 2: 50, 3: 1}, <= 40)'
   % (len(set(uni_rows)), dict(sorted(uni_kappa.items())), maxN))

# ---------------------------------------------------------------------------------------------- (3) covering
C = np.load(os.path.join(REPO, 'data', 's32_stab5_cover.npz'), allow_pickle=True)
pieces = {}
nQ = nGF = 0
fails = []
for e in range(len(C['row'])):
    row, k, n, kind = int(C['row'][e]), int(C['kappa'][e]), int(C['N'][e]), str(C['kind'][e])
    sizes = [int(s) for s in C['sizes'][e]]
    W = C['W'][e][:n, :n].astype(int).tolist()
    r = R[row]
    target = [k * int(x) for x in r]
    if not (sizes == expected_sizes(r, k) and structure_ok(W, sizes)
            and not C['W'][e][n:, :].any() and not C['W'][e][:, n:].any()):
        fails.append((e, 'structure'))
        continue
    if kind == 'Q':
        nQ += 1
        ranks, bad = bad_primes_and_ranks(W, sizes)
        stored = set(int(x) for x in str(C['bad'][e]).split(',') if x)
        if ranks != target:
            fails.append((e, 'rational ranks'))
        if bad != stored:
            fails.append((e, 'bad primes %s, stored %s' % (sorted(bad), sorted(stored))))
        pieces.setdefault(row, []).append(('Q', bad))
    elif kind == 'GF':
        nGF += 1
        p = int(C['prime'][e])
        if not sympy.isprime(p) or ranks_mod(W, sizes, p) != target:
            fails.append((e, 'GF(%d) ranks' % p))
        pieces.setdefault(row, []).append(('GF', p))
    else:
        fails.append((e, 'kind %r' % kind))
ok(not fails and nQ == 125 and nGF == 74,
   'cover file: %d rational pieces, %d GF pieces (expect 125, 74); failures: %s' % (nQ, nGF, fails[:5]))

uncovered = []
for row, ps in pieces.items():
    badsets = [b for kind, b in ps if kind == 'Q']
    gf = {p for kind, p in ps if kind == 'GF'}
    if not badsets:
        uncovered.append(row)
        continue
    common_bad = set.intersection(*badsets)
    if not common_bad <= gf:
        uncovered.append(row)
cover_rows = set(pieces)
ok(not uncovered and len(cover_rows) == 75,
   'covering: %d rows, every one covered (uncovered: %s)' % (len(cover_rows), uncovered))
ok(not (cover_rows & set(uni_rows)) and cover_rows | set(uni_rows) == set(range(641)),
   'the 566 unimodular rows and the 75 covered rows are disjoint and together are all 641 rows')

# ---------------------------------------------------------------------------------------------- (4) qubits, qutrits


def has_realization(row, p):
    if row in set(uni_rows):
        return True
    return any((kind == 'Q' and p not in x) or (kind == 'GF' and x == p) for kind, x in pieces[row])


uni_set = set(uni_rows)
q2 = sum(1 for row in range(641) if row in uni_set or has_realization(row, 2))
q3 = sum(1 for row in range(641) if row in uni_set or has_realization(row, 3))
ok(q2 == 641 and q3 == 641, 'explicit realizations: qubits %d rows, qutrits %d rows (expect 641, 641)' % (q2, q3))
TWELVE = [6, 7, 15, 34, 39, 40, 42, 43, 44, 46, 48, 50]
k12 = {int(U['row'][e]): int(U['kappa'][e]) for e in range(len(U['row'])) if int(U['row'][e]) in TWELVE}
ok(sorted(k12) == TWELVE and set(k12.values()) == {2},
   'the 12 rows of the BCHS list without a stored qubit realization before R38 have unimodular realizations, '
   'all at multiplicity 2: %s' % k12)

# ---------------------------------------------------------------------------------------------- (5) facets
EDGES = list(itertools.combinations(range(6), 2))
CUTS = []
for X in FIVE_MASKS:
    I = [i for i in range(5) if X >> i & 1]
    CUTS.append((I, [j for j in range(6) if j not in I]))
cache = {}


def cutvec(adj, p):
    out = []
    for I, J in CUTS:
        blk = tuple(tuple(adj[i][j] for j in J) for i in I)
        key = (blk, p)
        if key not in cache:
            cache[key] = rank_Q([list(b) for b in blk]) if p == 0 else rank_mod_p([list(b) for b in blk], p)
        out.append(cache[key])
    return tuple(out)


vec2, vecQ = set(), set()
for mask in range(1 << 15):
    adj = [[0] * 6 for _ in range(6)]
    for b, (i, j) in enumerate(EDGES):
        if mask >> b & 1:
            adj[i][j] = adj[j][i] = 1
    vec2.add(cutvec(adj, 2))
    vecQ.add(cutvec(adj, 0))
ok(len(vec2) == 760 and vec2 <= vecQ,
   'six-vertex 0/1 graph states: %d qubit vectors, all among the %d rational cut-rank vectors' % (len(vec2), len(vecQ)))


def full(S, M):
    """value on any subset M of the six systems (bit 5 = F) of a vector on the subsets of {A,...,E}"""
    if M in (0, 63):
        return 0
    if M >> 5 & 1:
        M ^= 63
    return S[M - 1]


PERMS = list(itertools.permutations(range(6)))


def images(r):
    out = set()
    for p in PERMS:
        out.add(tuple(full(r, sum(1 << p[i] for i in range(6) if M >> i & 1)) for M in FIVE_MASKS))
    return out


pool = set(vec2)
for row in range(641):
    pool |= images([int(x) for x in R[row]])
P = np.array(sorted(pool), dtype=np.int64)
CL = np.load(os.path.join(REPO, 'data', 's17_ineq_classes.npy')).astype(np.int64)
vmin = min(int((P @ CL[c]).min()) for c in range(len(CL)))
ok(CL.shape == (747, 31) and vmin >= 0,
   'pool of %d vectors in Stab^(d)_5 for every d; min over the pool of the 747 classes: %d' % (len(P), vmin))

PR = 1000003


def rank_mod_big(T):
    M = np.array(T, dtype=np.int64) % PR
    r = 0
    for c in range(M.shape[1]):
        nz = np.nonzero(M[r:, c])[0]
        if len(nz) == 0:
            continue
        i = r + nz[0]
        M[[r, i]] = M[[i, r]]
        inv = pow(int(M[r, c]), PR - 2, PR)
        M[r] = (M[r] * inv) % PR
        f = M[:, c].copy()
        f[r] = 0
        M = (M - np.outer(f, M[r])) % PR
        r += 1
        if r == M.shape[0]:
            break
    return r


facets = []
for c in range(len(CL)):
    T = P[(P @ CL[c]) == 0]
    if len(T) >= 30 and rank_mod_big(T) == 30:
        facets.append(c)
ok(len(facets) == 119, '%d of the 747 classes have 30 linearly independent pool vectors on their hyperplane, so they are '
                       'facets of Stab^(d)_5 for every d (expect 119)' % len(facets))

# the class of the inequality (star) of Theorem 6 (transcribed from the paper) is one of them
STAR_L = {'BCD': 1, 'ABCD': 1, 'ABE': 1, 'BCE': 2, 'ADE': 2, 'CDE': 1}
STAR_R = {'BC': 2, 'ABD': 1, 'ACD': 1, 'AE': 1, 'DE': 1, 'BCDE': 1, 'ABCDE': 1}
star = [0] * 31
for nm, c in STAR_L.items():
    star[sum(1 << 'ABCDE'.index(x) for x in nm) - 1] += c
for nm, c in STAR_R.items():
    star[sum(1 << 'ABCDE'.index(x) for x in nm) - 1] -= c


def form_image(g, p):
    h = [0] * 31
    for M in range(1, 32):
        if g[M - 1]:
            Nm = sum(1 << p[i] for i in range(5) if M >> i & 1)
            if Nm >> 5 & 1:
                Nm ^= 63
            h[Nm - 1] += g[M - 1]
    return h


def primitive(v):
    v = np.array(v, dtype=np.int64)
    return tuple(v // np.gcd.reduce(np.abs(v[v != 0])))


star_imgs = {primitive(form_image(star, p)) for p in PERMS}
star_cls = [c for c in range(len(CL)) if primitive(CL[c]) in star_imgs]
ok(len(star_cls) == 1 and star_cls[0] in facets,
   'the inequality (star) is class %s of the 747 and is among the 119 facet classes' % star_cls)

print('ALL OK' if not failures else 'FAILURE (%d)' % len(failures))
sys.exit(1 if failures else 0)
