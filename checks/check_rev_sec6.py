"""Checks for the revised paper (2026-10), Methods Sec. IV E and Supplementary Note 7: the form f_12 (Eq. (10)) and its
expansion (Eq. (S9)), Proposition 3, Theorem 3 (with the seven-point states), the form g_12 (Supplementary Note 7.2) and
Theorem 7.
Independent of the repository code.

Labels: A, B, C, W, X, Y, Z = 100, 010, 001, 110, 101, 011, 111; Z is the purifier.  A pure-symmetric vector S on the
seven systems is a dict {frozenset: value} on all subsets, with S(empty) = S(all) = 0 and S(T) = S(complement).

What is proved here exactly: the expansion of f_12 and of U_Z; the values on the Steane, non-Fano and seven-point
states (from the codes, and against state vectors numerically); the Smith invariants used in Theorem 3; the values
of Theorem 7 and of g_12 (Supplementary Note 7.2).  What is only tested: f_12 >= 0 on random qubit stabilizer states
and every step of the proof of Proposition 3 on them (Supplementary Lemmas 7.1-7.5, (LA), the assembly), and
f^(a) >= 0 on random Z_4 and Z_8 states.
"""
import itertools, math, os, random, sys
from fractions import Fraction

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import rank_mod_p

ok = True


def report(name, good, detail=''):
    global ok
    ok = ok and good
    print(('OK  ' if good else 'FAIL') + ' ' + name + ('  ' + detail if detail else ''))


L7 = 'ABCWXYZ'
ALL = frozenset(L7)
POINT = {'A': (1, 0, 0), 'B': (0, 1, 0), 'C': (0, 0, 1), 'W': (1, 1, 0), 'X': (1, 0, 1), 'Y': (0, 1, 1),
         'Z': (1, 1, 1)}
LINES = ['ABW', 'ACX', 'AYZ', 'BCY', 'BXZ', 'CWZ']
CORNERS = 'ABCZ'
fs = frozenset
SUBSETS = [fs(c) for r in range(0, 8) for c in itertools.combinations(L7, r)]
PARTY_SUBSETS = [fs(c) for r in range(1, 7) for c in itertools.combinations('ABCWXY', r)]


def lines_through(P):
    return [l for l in LINES if P in l]


# ------------------------------------------------------------------------------------- the form f_12, Eq. (10)
def f12_parts(S):
    k = lambda T: sum(S[fs(P)] for P in T) - S[fs(T)]
    d = lambda P, l: k(l) - k(l.replace(P, ''))

    def e(P, l):
        o = [m for m in lines_through(P) if m != l]
        assert len(o) == 2
        return 2 * S[fs(P)] - d(P, o[0]) - d(P, o[1])
    U = {
        'Z': 4 * S[fs('Z')] - d('Z', 'AYZ') - d('Z', 'BXZ') - d('Z', 'CWZ'),
        'W': 2 * S[fs('W')] - d('W', 'ABW') + e('A', 'ABW') + e('B', 'ABW'),
        'X': 2 * S[fs('X')] - d('X', 'ACX') + e('A', 'ACX') + e('C', 'ACX'),
        'Y': 2 * S[fs('Y')] - d('Y', 'BCY') + e('B', 'BCY') + e('C', 'BCY'),
    }
    cp = [a + b for a, b in itertools.combinations(CORNERS, 2)]
    f = (5 * U['Z'] + 4 * U['W'] + U['X'] + U['Y'] + k('WXY') + 6 * k('ABC') + sum(k(x) for x in cp)
         - sum(k(l) for l in LINES) - 3 * k('ABW') - 3 * k('CWZ'))
    return f, U, k, d, e


def f12(S):
    return f12_parts(S)[0]


def basis_vector(T):
    """pure-symmetric vector with S(T) = S(complement of T) = 1, zero elsewhere (T a set of parties)."""
    S = {X: 0 for X in SUBSETS}
    S[T] = 1
    S[ALL - T] = 1
    return S


# Eq. (S9) of the paper (Supplementary Note 7), transcribed term by term
EXPANSION = {'A': 3, 'B': 3, 'C': 3, 'W': -3, 'AB': -5, 'AC': -2, 'AW': -1, 'AX': -1, 'AY': -5, 'BC': -2, 'BW': -1,
             'BX': -5, 'BY': -1, 'CW': -5, 'CX': -4, 'CY': -4, 'ABC': -6, 'ABW': 10, 'ACX': 7, 'BCY': 7, 'WXY': -1,
             'ABXY': 11, 'ACWY': 11, 'BCWX': 11, 'ABCWX': -5, 'ABCWY': -5, 'ABCXY': -2, 'ABWXY': -1, 'ACWXY': -1,
             'BCWXY': -1, 'ABCWXY': 2}
coef = {T: f12(basis_vector(T)) for T in PARTY_SUBSETS}
nz = {''.join(c for c in L7 if c in T): v for T, v in coef.items() if v != 0}
report('expansion of f_12 (Eq. (10)) on the 63 party subsets equals Eq. (S9) term by term (%d terms)' % len(nz),
       nz == EXPANSION and len(nz) == 31)


def f12_expanded(S):
    return sum(v * S[fs(T)] for T, v in EXPANSION.items())


def U_Z(S):
    return S[fs('Z')] + sum(S[fs(l)] - S[fs(l.replace('Z', ''))] for l in LINES if 'Z' in l)


rng = random.Random(161016)
okU = True
for _ in range(200):
    S = {X: 0 for X in SUBSETS}
    for T in PARTY_SUBSETS:
        S[T] = S[ALL - T] = rng.randint(-5, 9)
    okU = okU and f12_parts(S)[1]['Z'] == U_Z(S) and f12(S) == f12_expanded(S)
report('U_Z = S(Z) + sum over the lines through Z of S(l) - S(l minus Z) (200 random pure-symmetric vectors)', okU)


# ------------------------------------------------------------------------------ codes, orders and CSS entropies
def egcd(a, b):
    if b == 0:
        return (abs(a), 1 if a >= 0 else -1, 0)
    g, x, y = egcd(b, a % b)
    return (g, y, x - (a // b) * y)


def lattice_index(gens, D, N):
    basis = [None] * D
    rows = [list(v) for v in gens] + [[N if i == j else 0 for j in range(D)] for i in range(D)]
    for r in rows:
        r = list(r)
        for i in range(D):
            if r[i] == 0:
                continue
            if basis[i] is None:
                if r[i] < 0:
                    r = [-x for x in r]
                basis[i] = r
                break
            b = basis[i]
            g, x, y = egcd(b[i], r[i])
            nb = [x * bb + y * rr for bb, rr in zip(b, r)]
            nr = [(r[i] // g) * bb - (b[i] // g) * rr for bb, rr in zip(b, r)]
            basis[i] = nb[:i + 1] + [v % N for v in nb[i + 1:]]
            r = nr[:i + 1] + [v % N for v in nr[i + 1:]]
    idx = 1
    for i in range(D):
        idx *= abs(basis[i][i])
    return idx


def group_order(gens, D, N):
    return N ** D // lattice_index(gens, D, N) if gens else 1


def css_orders(cols, N):
    """for a code over Z_N given by its generator matrix columns (one column per system), return the function
    X -> 2^S(X) (exact Fraction), S the entropy vector of the CSS state (Lemma A.5): S = r(X) + r(X^c) - r(all),
    r(X) = log2 of the order of the projection, which equals the order of the column space of the columns in X."""
    D = len(next(iter(cols.values())))
    cache = {}

    def r(X):
        if X not in cache:
            cache[X] = group_order([cols[P] for P in X], D, N)
        return cache[X]
    tot = r(ALL)
    return {X: Fraction(r(X) * r(ALL - X), tot) for X in SUBSETS}


def log2_exact(F):
    """log2 of a power of two (Fraction), as an integer."""
    if F.denominator == 1:
        n = F.numerator
        assert n & (n - 1) == 0
        return n.bit_length() - 1
    return -log2_exact(1 / F)


def entropy_from_ranks(rank):
    """S = lambda_r for a rank function given on all subsets (integers)."""
    return {X: rank[X] + rank[ALL - X] - rank[ALL] for X in SUBSETS}


def rank_function(cols, p):
    return {X: (rank_mod_p([cols[P] for P in X], p) if X else 0) for X in SUBSETS}


COLS = {P: list(v) for P, v in POINT.items()}
S_F = entropy_from_ranks(rank_function(COLS, 2))
S_NF = {p: entropy_from_ranks(rank_function(COLS, p)) for p in (3, 5, 7, 11)}
same_NF = all(S_NF[p] == S_NF[3] for p in S_NF)
SNF = S_NF[3]
vals = (f12(SNF), f12(S_F), f12_parts(SNF)[1], f12_parts(S_F)[1])
diff_cut = [X for X in SUBSETS if S_F[X] != SNF[X]]
report('values: f_12(S_NF) = -1 (S_NF the same over GF(3), GF(5), GF(7), GF(11)), f_12(S_F) = 0, U_Z = U_W = U_X = '
       'U_Y = 1 at both; S_F and S_NF differ only on WXY | ABCZ',
       vals[0] == -1 and vals[1] == 0 and same_NF and all(v == 1 for v in vals[2].values()) and
       all(v == 1 for v in vals[3].values()) and sorted(''.join(c for c in L7 if c in X) for X in diff_cut) ==
       ['ABCZ', 'WXY'])
kNF = f12_parts(SNF)[2]
report('at S_NF: S = 1 on points, k = 0 on pairs and on the triples that are not lines, k = 1 on the six lines',
       all(SNF[fs(P)] == 1 for P in L7) and all(kNF(''.join(c)) == 0 for c in itertools.combinations(L7, 2)) and
       all(kNF(''.join(c)) == (1 if ''.join(c) in LINES else 0) for c in itertools.combinations(L7, 3)))


# ------------------------------------------------------------------- state vectors (numerical, independent)
def css_state_vector(cols, d):
    """|C> = |C|^(-1/2) sum_{c in C} |c> on 7 qudits of dimension d (qudit order A, B, C, W, X, Y, Z)."""
    G = np.array([cols[P] for P in L7]).T          # 3 x 7
    words = set()
    for y in itertools.product(range(d), repeat=G.shape[0]):
        words.add(tuple((np.array(y) @ G) % d))
    psi = np.zeros(d ** 7)
    for w in words:
        idx = 0
        for x in w:
            idx = idx * d + int(x)
        psi[idx] = 1.0
    return psi / np.sqrt(len(words))


def entropy_bits(psi, d, X):
    """von Neumann entropy (bits) of the reduced state on the systems in X."""
    pos = [L7.index(P) for P in L7 if P in X]
    rest = [i for i in range(7) if i not in pos]
    T = psi.reshape([d] * 7).transpose(pos + rest).reshape(d ** len(pos), d ** len(rest))
    s = np.linalg.svd(T, compute_uv=False) ** 2
    s = s[s > 1e-14]
    return float(-(s * np.log2(s)).sum())


worst = 0.0
for d in (2, 3, 4, 6):
    psi = css_state_vector(COLS, d)
    Sd = css_orders(COLS, d)
    for X in PARTY_SUBSETS:
        if len(X) > 3:
            continue
        exact = math.log2(Sd[X].numerator) - math.log2(Sd[X].denominator)
        worst = max(worst, abs(entropy_bits(psi, d, X) - exact))
report('seven-point states over Z_2, Z_3, Z_4, Z_6 from state vectors: entropies agree with the CSS formula '
       '(largest deviation %.1e bits)' % worst, worst < 1e-9)


# ----------------------------------------------------------------- Theorem 3 (seven-point states)
def det(m):
    m = [list(r) for r in m]
    n = len(m)
    if n == 1:
        return m[0][0]
    return sum((-1) ** j * m[0][j] * det([r[:j] + r[j + 1:] for r in m[1:]]) for j in range(n))


def invariant_factors(colset):
    """invariant factors over Z of the 3 x |X| matrix of the columns in X (determinantal divisors)."""
    M = [[COLS[P][i] for P in colset] for i in range(3)]
    out, prev = [], 1
    for kk in range(1, min(3, len(colset)) + 1):
        g = 0
        for rows in itertools.combinations(range(3), kk):
            for cs in itertools.combinations(range(len(colset)), kk):
                g = math.gcd(g, det([[M[r][c] for c in cs] for r in rows]))
        if g == 0:
            break
        out.append(g // prev)
        prev = g
    return out


inv = {''.join(c for c in L7 if c in X): invariant_factors(sorted(X, key=L7.index)) for X in SUBSETS if X}
special = {k: v for k, v in inv.items() if any(x != 1 for x in v)}
two_by_two = {det([[COLS[P][r] for P in pr] for r in rows]) for pr in itertools.combinations(L7, 2)
              for rows in itertools.combinations(range(3), 2)}
three = {''.join(t): det([[COLS[P][r] for P in t] for r in range(3)]) for t in itertools.combinations(L7, 3)}
big3 = {k: v for k, v in three.items() if abs(v) > 1}
report('Smith invariants of the 127 column sets: all 1 except WXY with 1, 1, 2; 2x2 minors in {0, +-1}; the only 3x3 '
       'minor of absolute value > 1 is det(W, X, Y) = %s' % big3.get('WXY'),
       special == {'WXY': [1, 1, 2]} and two_by_two <= {0, 1, -1} and set(big3) == {'WXY'})

bad = []
for dd in range(2, 33):
    Sd = css_orders(COLS, dd)
    F = Fraction(1)
    for T, v in EXPANSION.items():
        F *= Sd[fs(T)] ** v
    if F != Fraction(math.gcd(dd, 2), dd):
        bad.append(dd)
report('Theorem 3: 2^f_12(S^(d)) = gcd(d,2)/d exactly, i.e. f_12(S^(d)) = -log2(d/gcd(d,2)), for d = 2..32',
       not bad, str(bad))

S4 = {X: log2_exact(v) for X, v in css_orders(COLS, 4).items()}
report('S^(4) = S_F + S_NF on all 128 subsets (bits; S_NF normalized to S(P) = 1)',
       all(S4[X] == S_F[X] + SNF[X] for X in SUBSETS))

# ------------------------------------------------------------------------------------------- Theorem 7 values
bad = []
for a in range(1, 7):
    if f12(SNF) + (1 - Fraction(1, a)) * U_Z(SNF) != -Fraction(1, a):
        bad.append(('NF', a))
    for j in range(1, 7):
        Sj = {X: log2_exact(v) for X, v in css_orders(COLS, 2 ** j).items()}
        val = f12(Sj) + (1 - Fraction(1, a)) * U_Z(Sj)
        if val != 1 - Fraction(j, a) or U_Z(Sj) != j or f12(Sj) != 1 - j:
            bad.append((a, j))
        if not all(Sj[X] == S_F[X] + (j - 1) * SNF[X] for X in SUBSETS):
            bad.append(('decomposition', j))
report('Theorem 7: f^(a)(S_NF) = -1/a, U_Z(S^(2^j)) = j, f_12(S^(2^j)) = 1 - j, f^(a)(S^(2^j)) = 1 - j/a, and '
       'S^(2^j) = S_F + (j-1) S_NF (a, j = 1..6)', not bad, str(bad[:4]))


# --------------------------------------------------------------------------------------- g_12, Supplementary Note 7.2
def lam(hf):
    return {X: hf[X] + hf[ALL - X] - hf[ALL] for X in SUBSETS}


def dual(hf):
    return {X: sum(hf[fs(i)] for i in X) + hf[ALL - X] - hf[ALL] for X in SUBSETS}


def g12(hf):
    return f12(lam(hf))


rNF = rank_function(COLS, 3)
rF = rank_function(COLS, 2)
hNF = {X: rNF[X] + dual(rNF)[X] for X in SUBSETS}
okD = True
for _ in range(200):
    hf = {X: (rng.randint(0, 9) if X else 0) for X in SUBSETS}
    okD = okD and g12(dual(hf)) == g12(hf)
report('Supplementary Note 7.2: g_12(r_NF) = g_12(r_NF*) = -1, g_12(r_F) = g_12(r_F*) = 0, g_12(h_NF) = -2, and '
       'g_12 o D = g_12 (200 random set functions)',
       (g12(rNF), g12(dual(rNF)), g12(rF), g12(dual(rF)), g12(hNF)) == (-1, -1, 0, 0, -2) and okD)


# ------------------------------------------- random stabilizer states (tests, not proofs): qubits, Z_4 and Z_8
def graph_state_S(adj, owner, d):
    """entropy vector (log2, exact as a Fraction power of 2 for d a power of 2, else via orders) of the graph
    state with symmetric weight matrix adj over Z_d (App. A): S(X) = log2 |W[Q_X, Q_rest] Z_d^...|."""
    n = len(adj)
    out = {}
    for X in SUBSETS:
        if not X or X == ALL:
            out[X] = Fraction(1)
            continue
        inX = [v for v in range(n) if owner[v] in X]
        rest = [v for v in range(n) if owner[v] not in X]
        cols = [[adj[a][b] for a in inX] for b in rest]
        out[X] = Fraction(group_order(cols, len(inX), d))
    return out


def random_graph(n, d, rng, density=0.5):
    adj = [[0] * n for _ in range(n)]
    for a in range(n):
        for b in range(a + 1, n):
            if rng.random() < density:
                adj[a][b] = adj[b][a] = rng.randrange(1, d)
    return adj


def to_bits(Sord):
    return {X: log2_exact(v) for X, v in Sord.items()}


minf, cnt = 10 ** 9, 0
for m, nsamp in ((1, 3000), (2, 600), (3, 120)):
    for _ in range(nsamp):
        owner = [P for P in L7 for _ in range(m)]
        S = to_bits(graph_state_S(random_graph(7 * m, 2, rng, rng.choice([0.3, 0.5, 0.7])), owner, 2))
        minf = min(minf, f12(S))
        cnt += 1
report('f_12 >= 0 on %d random qubit graph states with 1, 2 or 3 qubits per system (minimum %d)' % (cnt, minf),
       minf >= 0)


def random_code_cols(d, rng, rows=3, nonfree=False):
    cols = {}
    for P in L7:
        v = [rng.randrange(d) for _ in range(rows)]
        if nonfree and rng.random() < 0.5:
            v = [(2 * x) % d for x in v]
        cols[P] = v
    return cols


def perturbed_seven_point(d, rng):
    cols = {P: list(v) for P, v in POINT.items()}
    for P in L7:
        if rng.random() < 0.3:
            cols[P] = [(x + (d // 2) * rng.randrange(2)) % d for x in cols[P]]
    return cols


res = {}
for a in (2, 3):
    d = 2 ** a
    mn, mn12, n = None, None, 0
    for trial in range(900):
        kind = trial % 3
        if kind == 0:
            owner = list(L7)
            S = to_bits(graph_state_S(random_graph(7, d, rng), owner, d))
        elif kind == 1:
            S = {X: log2_exact(v) for X, v in css_orders(random_code_cols(d, rng, rng.choice([3, 4]), True), d).items()}
        else:
            S = {X: log2_exact(v) for X, v in css_orders(perturbed_seven_point(d, rng), d).items()}
        val = f12(S) + (1 - Fraction(1, a)) * U_Z(S)
        mn = val if mn is None else min(mn, val)
        mn12 = f12(S) if mn12 is None else min(mn12, f12(S))
        n += 1
    res[d] = (n, mn, mn12)
report('f^(a) >= 0 on random Z_4 (a = 2) and Z_8 (a = 3) states: graph states, CSS states of random and non-free '
       'codes, perturbed seven-point codes (count, min f^(a), min f_12): %s' % {k: (v[0], str(v[1]), v[2]) for k, v in
                                                                            res.items()},
       all(v[1] >= 0 for v in res.values()) and all(v[2] < 0 for v in res.values()))


# ------------------------------------------------------- the proof of Proposition 3 step by step (qubits; tests)
XM = int('01' * 64, 2)       # x bits at even positions (bit 2q), z bits at odd positions (bit 2q + 1)


def omega(u, v):
    xu, zu = u & XM, (u >> 1) & XM
    xv, zv = v & XM, (v >> 1) & XM
    return (bin(xu & zv).count('1') + bin(zu & xv).count('1')) & 1


def gf2_reduce(basis, v):
    for b in basis:
        if v ^ b < v:
            v ^= b
    return v


def gf2_basis(vecs):
    basis = []
    for v in vecs:
        v = gf2_reduce(basis, v)
        if v:
            basis.append(v)
            basis.sort(reverse=True)
    return basis


def gf2_rank(vecs):
    return len(gf2_basis(vecs))


def supported_in(basis, mask):
    """basis of the elements of span(basis) supported inside mask."""
    pairs = []     # (outside part, full vector), echelon on the outside part
    kernel = []
    for v in basis:
        key, val = v & ~mask, v
        for k2, v2 in pairs:
            if key ^ k2 < key:
                key ^= k2
                val ^= v2
        if key:
            pairs.append((key, val))
            pairs.sort(reverse=True)
        else:
            kernel.append(val)
    return gf2_basis(kernel)


def form_rank(vecs):
    vecs = gf2_basis(vecs)
    rows = []
    for u in vecs:
        r = 0
        for j, v in enumerate(vecs):
            if omega(u, v):
                r |= 1 << j
        rows.append(r)
    return gf2_rank(rows)


def graph_L(adj):
    """generators X_q Z^(row q) of the graph state, as bit vectors (x at bit 2q, z at bit 2q + 1)."""
    n = len(adj)
    gens = []
    for q in range(n):
        v = 1 << (2 * q)
        for r in range(n):
            if adj[q][r]:
                v |= 1 << (2 * r + 1)
        gens.append(v)
    return gf2_basis(gens)


def chain_check(adj, owner):
    return chain_check_L(graph_L(adj), owner)


def chain_check_L(Lb, owner):
    """runs every step of the proof of Proposition 3 on the stabilizer state with Lagrangian span(Lb) (qubit q owned by
    owner[q]); entropies from S(X) = m_X - dim L_X (Lemma A.1). Returns the list of failed steps, or None if the
    state is not in normal form."""
    n = len(owner)
    fails = []
    sysmask = {P: sum(3 << (2 * q) for q in range(n) if owner[q] == P) for P in L7}
    mask = lambda T: sum(sysmask[P] for P in T)
    assert len(Lb) == n and all(omega(u, v) == 0 for u in Lb for v in Lb)
    S = {X: sum(owner.count(P) for P in X) - len(supported_in(Lb, mask(X))) for X in SUBSETS}
    if any(S[X] != S[ALL - X] for X in SUBSETS):
        fails.append('purity')
    if any(S[fs(P)] != owner.count(P) for P in L7):
        return None                                  # not in normal form; skipped
    f, U, k, d, e = f12_parts(S)
    # dim L_T = k(T)
    for T in ['A', 'AB', 'ABW', 'WXY', 'ABCZ'] + LINES:
        if len(supported_in(Lb, mask(T))) != k(T):
            fails.append('dimLT')
    Ll = {l: supported_in(Lb, mask(l)) for l in LINES}
    res_ = lambda vecs, P: gf2_basis([v & sysmask[P] for v in vecs])
    for l in LINES:
        for P in l:
            if len(res_(Ll[l], P)) != d(P, l):
                fails.append('d(P,l)')
    for l1, l2 in itertools.combinations(LINES, 2):
        P = (set(l1) & set(l2)).pop()
        if any(omega(u, v) for u in res_(Ll[l1], P) for v in res_(Ll[l2], P)):
            fails.append('orthogonality')
    for P in CORNERS:                                        # Lemma 7.1
        ls = lines_through(P)
        V = [res_(Ll[l], P) for l in ls]
        dd = [len(x) for x in V]
        s = S[fs(P)]
        for i in range(3):
            j, kk = [x for x in range(3) if x != i]
            if form_rank(V[i]) > 2 * s - dd[j] - dd[kk]:
                fails.append('7.1 rank')
        if gf2_rank(V[0] + V[1] + V[2]) > 4 * s - sum(dd):
            fails.append('7.1 sum')
    for l in LINES:                                          # Lemma 7.2, every choice of Q
        for Q in l:
            P1, P2 = [x for x in l if x != Q]
            if form_rank(res_(Ll[l], Q)) > form_rank(res_(Ll[l], P1)) + form_rank(res_(Ll[l], P2)):
                fails.append('7.2')
    Z0 = gf2_basis([v for l in LINES for v in Ll[l]])
    D = {P: res_(Z0, P) for P in L7}
    delta = {P: len(D[P]) for P in L7}
    for P in 'ZWXY':                                         # Lemma 7.3
        if delta[P] > U[P]:
            fails.append('7.3')
    z = {}
    for T in SUBSETS:
        z[T] = len(supported_in(Z0, mask(T))) if T else 0
    if any(z[T] > k(''.join(T)) for T in SUBSETS if T):      # Lemma 7.4
        fails.append('7.4(i)')
    if any(z[fs(l)] != k(l) for l in LINES):
        fails.append('7.4(ii)')
    if z[fs('ABCZ')] < k('ABW') + k('CWZ') - delta['W']:
        fails.append('7.4(iii)')
    r = {}
    for T in SUBSETS:                                        # Lemma 7.5, built explicitly
        GT = [v for P in T for v in D[P]]
        r[T] = gf2_rank(GT + Z0) - len(Z0)
        if r[T] != sum(delta[P] for P in T) - z[T]:
            fails.append('7.5')
    R = lambda s_: r[fs(s_)]
    alpha = R('A') + R('B') + R('C') - R('ABC')
    gamma = R('ABCZ') - R('ABC')
    zeta = (R('BCY') - R('BC')) + (R('ACX') - R('AC')) + (R('ABW') - R('AB'))
    eta = (R('AYZ') - R('AZ')) + (R('BXZ') - R('BZ')) + (R('CWZ') - R('CZ'))
    slackLA = 2 * R('Z') + 3 * alpha + 3 * gamma + zeta + eta - R('WXY')
    zz = lambda s_: z[fs(s_)]
    cp = [a_ + b_ for a_, b_ in itertools.combinations(CORNERS, 2)]
    G1 = (5 * delta['Z'] + delta['W'] + delta['X'] + delta['Y'] + zz('WXY') + 6 * zz('ABC') + sum(zz(x) for x in cp)
          - 3 * zz('ABCZ') - sum(k(l) for l in LINES))
    if slackLA < 0:
        fails.append('(LA)')
    if slackLA != G1:
        fails.append('assembly identity (S8)')
    G1b = (5 * delta['Z'] + 4 * delta['W'] + delta['X'] + delta['Y'] + k('WXY') + 6 * k('ABC') + sum(k(x) for x in cp)
           - 3 * k('ABW') - 3 * k('CWZ') - sum(k(l) for l in LINES))
    if not (G1 <= G1b <= f):
        fails.append('bounds after (S8)')
    if f < 0:
        fails.append('f12 < 0')
    return fails


nst, failed = 0, {}
for m, nsamp in ((1, 300), (2, 200), (3, 40)):
    tries = 0
    while tries < nsamp:
        owner = [P for P in L7 for _ in range(m)]
        adj = random_graph(7 * m, 2, rng, rng.choice([0.3, 0.5, 0.7]))
        fl_list = chain_check(adj, owner)
        if fl_list is None:
            continue
        tries += 1
        nst += 1
        for fl in fl_list:
            failed[fl] = failed.get(fl, 0) + 1
report('proof of Proposition 3 step by step on %d random qubit graph states in normal form (1, 2 or 3 qubits per system): '
       'dim L_T = k(T), k_P(l), orthogonality, Lemmas 7.1-7.5 (7.5 built explicitly), (LA), the identity (S8), the '
       'bounds after it, f_12 >= 0' % nst, not failed, str(failed))


# structured states near the boundary: Steane copies after random local symplectic maps, and CSS states of
# multi-copy point codes with random perturbations (pushed towards the Fano and non-Fano dependencies)
def nullspace_gf2(rows, n):
    """basis (bit vectors over n bits) of {x : <x, r> = 0 for all rows r}."""
    piv, basis = {}, []
    for r in gf2_basis(rows):
        piv[r.bit_length() - 1] = r
    red = gf2_basis(rows)
    # reduced echelon form
    red.sort(reverse=True)
    for i in range(len(red)):
        for j in range(len(red)):
            if i != j and red[j] >> (red[i].bit_length() - 1) & 1:
                red[j] ^= red[i]
    pivots = {r.bit_length() - 1: r for r in red}
    free = [b for b in range(n) if b not in pivots]
    out = []
    for fb in free:
        x = 1 << fb
        for pb, r in pivots.items():
            if r >> fb & 1:
                x |= 1 << pb
        out.append(x)
    return out


def css_L(code_rows, n):
    """Lagrangian {(x, z): x in C, z in C^perp} of the CSS state of the binary code spanned by code_rows (n bits)."""
    C = gf2_basis(code_rows)
    Cp = nullspace_gf2(C, n)
    spread = lambda v, off: sum(1 << (2 * q + off) for q in range(n) if v >> q & 1)
    return gf2_basis([spread(c, 0) for c in C] + [spread(b, 1) for b in Cp])


def random_symplectic_on(Lb, qubits, rng, steps=12):
    """apply a random product of transvections on the 2|qubits| coordinates of the given qubits to every vector."""
    sub = sum(3 << (2 * q) for q in qubits)
    for _ in range(steps):
        v = 0
        for q in qubits:
            v |= rng.randrange(4) << (2 * q)
        if v == 0:
            continue
        Lb = [u ^ v if omega(u & sub, v) else u for u in Lb]
    return gf2_basis(Lb)


PT = {P: sum(b << i for i, b in enumerate(v)) for P, v in POINT.items()}     # point as a 3-bit mask
nst2, failed2, fvals = 0, {}, []
for trial in range(240):
    kind = trial % 3
    if kind == 0:                      # m copies of the Steane state, local symplectic maps
        m = rng.choice([1, 2])
        owner = [P for P in L7 for _ in range(m)]
        rows = []
        for c in range(m):
            for i in range(3):
                rows.append(sum(1 << (L7.index(P) * m + c) for P in L7 if PT[P] >> i & 1))
        Lb = css_L(rows, 7 * m)
        for P in L7:
            Lb = random_symplectic_on(Lb, [q for q in range(7 * m) if owner[q] == P], rng)
    else:                              # perturbed multi-copy point codes over GF(2)
        m = rng.choice([1, 2])
        owner = [P for P in L7 for _ in range(m)]
        rows = []
        for c in range(m):
            for i in range(3):
                rows.append(sum(1 << (L7.index(P) * m + c) for P in L7 if PT[P] >> i & 1))
        for _ in range(rng.randint(0, 2)):
            rows[rng.randrange(len(rows))] ^= 1 << rng.randrange(7 * m)
        Lb = css_L(rows, 7 * m)
        if kind == 2:
            for P in L7:
                Lb = random_symplectic_on(Lb, [q for q in range(7 * m) if owner[q] == P], rng)
    fl_list = chain_check_L(Lb, owner)
    if fl_list is None:
        continue
    nst2 += 1
    sysmask_ = {P: sum(3 << (2 * q) for q in range(len(owner)) if owner[q] == P) for P in L7}
    S_ = {X: sum(owner.count(P) for P in X) - len(supported_in(Lb, sum(sysmask_[P] for P in X))) for X in SUBSETS}
    fvals.append(f12(S_))
    for fl in fl_list:
        failed2[fl] = failed2.get(fl, 0) + 1
report('proof of Proposition 3 step by step on %d structured qubit states in normal form (Steane copies after random '
       'local symplectic maps; perturbed one- and two-copy point codes): every step holds; f_12 = 0 on %d of them, '
       'minimum %d' % (nst2, fvals.count(0), min(fvals)), not failed2 and nst2 > 50 and min(fvals) >= 0, str(failed2))

# --------- familiar code states as witnesses (R48 addendum 2; R32 Proposition 26; Supplementary Note 8.1)
def lemma1_forms_sets(perm=('A', 'B', 'C', 'Y', 'X', 'W', 'Z')):
    """(LB) and (LA) with t = 2 (the paper's Lemma 1) as dicts {frozenset: coefficient} on the subsets of the seven
    elements, with (A_1, A_2, A_3, B_1, B_2, B_3, C) = perm; for the labels above, A_i = <e_i>, B_i = <c - e_i>, C = <c>."""
    M = 3
    A, B, C = [fs(x) for x in perm[:3]], [fs(x) for x in perm[3:6]], fs(perm[6])
    allA = fs().union(*A)

    def H(X, c=1):
        return {X: c} if X else {}

    def add(*forms):
        out = {}
        for f in forms:
            for k, v in f.items():
                out[k] = out.get(k, 0) + v
        return {k: v for k, v in out.items() if v}

    def sc(f, c):
        return {k: c * v for k, v in f.items()}

    def Hc(x, y):
        return add(H(x | y), H(y, -1))

    def I(x, y):
        return add(H(x), H(y), H(x | y, -1))
    alpha = add(*[H(a) for a in A], H(allA, -1))
    gamma = Hc(C, allA)
    delta = add(*[I(C, allA - A[i]) for i in range(M)])
    eps = add(*[Hc(C, A[i] | B[i]) for i in range(M)])
    zeta = add(*[Hc(B[i], allA - A[i]) for i in range(M)])
    eta = add(*[Hc(B[i], A[i] | C) for i in range(M)])
    HB = H(fs().union(*B))
    LB = add(HB, sc(add(gamma, sc(alpha, M), delta, eps, zeta), M), H(C, -M))
    LA = add(H(C, M - 1), sc(alpha, M), sc(gamma, M), zeta, eta, sc(HB, -1))
    return LB, LA


def evalf(f, h):
    return sum(c * h[X] for X, c in f.items())


def balanced_value(f, h):
    """f_bal(h) = f(h) - sum_i f(1_i) (h(all) - h(all - i)), balanced over the seven elements (the parties), as in the
    proof of Proposition 4 (Supplementary Note 8.1); balancing over all eight systems, as in the proof of Corollary 3,
    would give +20 (LB) and +7 (LA) on the code states below."""
    return evalf(f, h) - sum(sum(c for X, c in f.items() if P in X) * (h[ALL] - h[ALL - fs(P)]) for P in L7)


LB2, LA2 = lemma1_forms_sets()
rF = rank_function(COLS, 2)
rNF = {p: rank_function(COLS, p) for p in (3, 5, 7)}
report('(LB) with t = 2, balanced over the seven elements, is -1 on the Fano rank function r_F; (LA) with t = 2, '
       'balanced, is -1 on the non-Fano rank function over GF(3), GF(5), GF(7); each is >= 0 on the other one',
       balanced_value(LB2, rF) == -1 and all(balanced_value(LA2, rNF[p]) == -1 for p in rNF)
       and all(balanced_value(LB2, rNF[p]) >= 0 for p in rNF) and balanced_value(LA2, rF) >= 0)
report('contrast (Proposition 2): on the Steane vector S_F itself the balanced (LB) with t = 2 is >= 0 (value %d), '
       'and the balanced (LA) is >= 0 on S_NF (value %d)' % (balanced_value(LB2, S_F), balanced_value(LA2, SNF)),
       balanced_value(LB2, S_F) >= 0 and balanced_value(LA2, SNF) >= 0)
rng_fam = random.Random(26)
okv = True
for p, form in ((3, LB2), (5, LB2), (2, LA2)):
    for _ in range(150):
        cols_r = {P: [rng_fam.randrange(p) for _ in range(3)] for P in L7}
        okv &= balanced_value(form, rank_function(cols_r, p)) >= 0
report('tests (not proofs): the balanced (LB) is >= 0 on 300 random rank functions over GF(3) and GF(5), the balanced '
       '(LA) on 150 over GF(2)', okv)


def projective_points(p, m):
    """one representative (first nonzero coordinate 1) of each point of PG(m-1, p)"""
    pts = []
    for v in itertools.product(range(p), repeat=m):
        nz = [x for x in v if x]
        if nz and nz[0] == 1:
            pts.append(tuple(v))
    return pts


def subconfiguration_entropies(p, m):
    """the CSS state of the p-ary simplex code of dimension m (one qudit per point of PG(m-1, p)); the seven parties are
    the points (POINT[P], 0, ..., 0), all other qudits form the purifier.  Returns S on the subsets of the seven parties
    (in units of log p), from S(X) = r(X) + r(rest) - r(all) with r the column rank, and the rank of the purifier."""
    pts = projective_points(p, m)
    party = {P: tuple(POINT[P]) + (0,) * (m - 3) for P in L7}
    assert all(v in pts for v in party.values())
    rest = [v for v in pts if v not in party.values()]
    r_all = rank_mod_p(pts, p)
    S = {}
    for X in SUBSETS:
        inX = [party[P] for P in X]
        out = rest + [party[P] for P in L7 if P not in X]
        S[X] = (rank_mod_p(inX, p) if inX else 0) + rank_mod_p(out, p) - r_all
    return S, rank_mod_p(rest, p), r_all, len(pts)


fam_ok, fam_txt = True, []
for p, m, target in ((2, 4, rF), (2, 5, rF), (3, 3, rNF[3]), (5, 3, rNF[5]), (7, 3, rNF[7])):
    S_fam, r_rest, r_all, npts = subconfiguration_entropies(p, m)
    good = r_rest == r_all and all(S_fam[X] == target[X] for X in SUBSETS)
    form = LB2 if p == 2 else LA2
    good &= balanced_value(form, S_fam) == -1
    fam_ok &= good
    fam_txt.append('p=%d m=%d (%d qudits)' % (p, m, npts))
report('simplex codes: in the CSS state of the binary simplex code of dimension m = 4, 5 the qubits outside a Fano '
       'subplane determine the codeword and the seven qubits of the subplane have S = r_F; for the p-ary simplex code '
       'of dimension 3, p = 3, 5, 7, the seven qudits {0,1}^3 minus 0 have S = r_NF; the balanced (LB), resp. (LA), with '
       't = 2 is -1 there', fam_ok, ', '.join(fam_txt))


def simplex_state_entropies(p, m, maxsize):
    """state vector of the CSS state of the p-ary simplex code of dimension m; entropies (in units of log p) of the
    sets X of the seven parties with |X| <= maxsize, from singular values"""
    pts = projective_points(p, m)
    party_idx = {P: pts.index(tuple(POINT[P]) + (0,) * (m - 3)) for P in L7}
    n = len(pts)
    G = np.array(pts).T
    psi = np.zeros(p ** n)
    for y in itertools.product(range(p), repeat=m):
        w = (np.array(y) @ G) % p
        idx = 0
        for x in w:
            idx = idx * p + int(x)
        psi[idx] = 1.0
    psi /= np.linalg.norm(psi)
    T0 = psi.reshape([p] * n)
    out = {}
    for X in SUBSETS:
        if not X or len(X) > maxsize:
            continue
        pos = [party_idx[P] for P in L7 if P in X]
        rest = [i for i in range(n) if i not in pos]
        T = T0.transpose(pos + rest).reshape(p ** len(pos), p ** len(rest))
        s = np.linalg.svd(T, compute_uv=False) ** 2
        s = s[s > 1e-14]
        out[X] = float(-(s * np.log(s)).sum() / math.log(p))
    return out


dev = 0.0
for p, m, target, maxsize in ((2, 4, rF, 7), (3, 3, rNF[3], 4)):
    num = simplex_state_entropies(p, m, maxsize)
    dev = max(dev, max(abs(num[X] - target[X]) for X in num))
report('state vectors: the 15-qubit binary and the 13-qutrit ternary simplex states give the same entropies on the seven '
       'parties (all subsets, resp. those of size <= 4; largest deviation %.1e)' % dev, dev < 1e-9)

print('ALL OK' if ok else 'FAILURE')
sys.exit(0 if ok else 1)
