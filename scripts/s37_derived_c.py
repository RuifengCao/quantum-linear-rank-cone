#!/usr/bin/env python3
"""s37: fewer parties for two primes, and the T8 / T8' vectors at seven parties (R46;
docs/R32_reduction_and_local_dimension.md, sections 4.13, 4.14 and 4.5).

Lemma 4* (section 4.13).  For subspaces A_1..A_M, B_1..B_M (M = t + 1) the derived subspace
C = (A_1 + B_1) ∩ ... ∩ (A_M + B_M) turns Lemma 4 into two inequalities on the 2M variables alone: (LB*) when the
characteristic does not divide t, (LA*) when it does.  On T(t) = [I_M | J_M - I_M] (A_i = <e_i>, B_i = <c - e_i>) with
M >= 4 they take the value -1 in the other case, and 0 otherwise.  Through doubled codes (Proposition 3) this gives
Theorem 17: for primes p < q, Stab^(p)_n and Stab^(q)_n are incomparable for every n >= 2p + 2 (Theorem 2' had 2p + 3).

Proposition 18 (section 4.14).  2 S_T8' is the entropy vector of a qutrit and of a qubit stabilizer state (two qudits
per system), 2 S_T8 that of a Z_p stabilizer state for every prime p >= 5 (an integer two-copy family, bad primes 2
and 3), and 6 S_T8 that of a qubit state: the eight-element analogue of the Fano / non-Fano pair witnesses nothing at
seven parties.

DFZ arXiv:1401.2507v1 (section 4.5).  Theorem 4.1 (characteristic 3) transcribed twice and checked; Theorem 3.1
(characteristic != 3), misprinted, re-expanded from the bounds of its own proof (pp. 13-19, 25).

Variables: A_1..A_M = bits 0..M-1, B_1..B_M = bits M..2M-1.  For M = 4 these are A, B, C, D, W, X, Y, Z (DFZ's names;
v_1..v_4, w_1..w_4 of section 4.10).  A form f is an integer vector over the subsets with f . h >= 0.

  python3 scripts/s37_derived_c.py --check [--full]    # gates G42 / G42b
  python3 scripts/s37_derived_c.py --forms             # (LB*) and (LA*) for M = 4, term by term
"""
import argparse, hashlib, itertools, math, os, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, '..')
DATA = os.path.join(ROOT, 'data', 's37_derived_forms.npz')
NAMES8 = 'ABCDWXYZ'


# ------------------------------------------------------------------------------------------------ the forms
def _e(N, X, k=1):
    v = np.zeros(N, dtype=np.int64)
    if X:
        v[X] += k
    return v


def _H_cond(N, x, y):                    # H(x | y)
    return _e(N, x | y) - _e(N, y)


def _I(N, x, y, z=0):                    # I(x ; y | z)
    return _e(N, x | z) + _e(N, y | z) - _e(N, x | y | z) - _e(N, z)


def lemma4_star(M, jmap=None):
    """(LB*, LA*) for M = t + 1; jmap[i] = j(i) != i (default i + 1 mod M).  Section 4.13."""
    n = 2 * M
    N = 1 << n
    if jmap is None:
        jmap = [(i + 1) % M for i in range(M)]
    assert all(jmap[i] != i for i in range(M))
    a = lambda *idx: sum(1 << i for i in idx)
    b = lambda *idx: sum(1 << (M + i) for i in idx)
    P = lambda *idx: a(*idx) | b(*idx)
    p = lambda *idx: _e(N, P(*idx))
    allA, allB = (1 << M) - 1, ((1 << M) - 1) << M
    Aneq = lambda i: allA & ~(1 << i)
    L = p(0) + p(1) - p(0, 1)
    for k in range(2, M):
        L = L + p(k) - p(0, k) - p(1, k) + p(0, 1, k)
    kappa = _e(N, N - 1) - _e(N, allA)
    alpha = sum(_e(N, a(k)) for k in range(M)) - _e(N, allA)
    zeta = [_H_cond(N, b(i), Aneq(i)) for i in range(M)]
    delta = [_I(N, P(jmap[i]), b(i)) + _I(N, a(i), Aneq(i), b(i)) for i in range(M)]
    theta = [_I(N, a(i), P(jmap[i])) for i in range(M)]
    HB = _e(N, allB)
    LB = HB - M * L + M * (kappa + M * alpha + sum(delta[i] + zeta[i] for i in range(M)))
    LA = (sum(p(i) - _e(N, a(i)) for i in range(M)) - L + sum(theta) + M * alpha + M * kappa + sum(zeta) - HB)
    return LB, LA


def balanced(f, n):
    """Chan's balancing: f_bal = f - sum_i c_i (e_E - e_{E - i}), c_i = sum of the coefficients of the sets containing i."""
    N = 1 << n
    c = np.array([int(f[[X for X in range(N) if X >> i & 1]].sum()) for i in range(n)], dtype=np.int64)
    fb = f.copy()
    for i in range(n):
        fb[N - 1] -= c[i]
        fb[(N - 1) ^ (1 << i)] += c[i]
    return c, fb


def terms(f, n, names=None):
    names = names or ([f'A{i + 1}' for i in range(n // 2)] + [f'B{i + 1}' for i in range(n // 2)])
    out = [(''.join(names[i] for i in range(n) if X >> i & 1), int(f[X])) for X in range(1, 1 << n) if f[X]]
    out.sort(key=lambda t: (len(t[0]), t[0]))
    return out


def form_sha(f):
    return hashlib.sha256(np.ascontiguousarray(f, dtype='<i8').tobytes()).hexdigest()


# ------------------------------------------------------------------------------------------------ fields and ranks
class GF:
    """GF(p^k) = GF(p)[x]/(x^k - sum_i low[i] x^i); elements are tuples of k coefficients."""
    def __init__(self, p, low=()):
        self.p, self.k = p, max(1, len(low))
        self.low = [c % p for c in low] if low else [0]
        self.els = list(itertools.product(range(p), repeat=self.k))
        self.zero, self.one = tuple([0] * self.k), tuple([1] + [0] * (self.k - 1))
        self.gen = tuple([0, 1] + [0] * (self.k - 2)) if self.k > 1 else (0,)
        nz = [u for u in self.els if u != self.zero]
        self._inv = {}
        for u in nz:
            for v in nz:
                if self.mul(u, v) == self.one:
                    self._inv[u] = v
                    break
        assert len(self._inv) == len(nz), 'not a field'

    def mul(self, u, v):
        p, k = self.p, self.k
        prod = [0] * (2 * k - 1)
        for i, x in enumerate(u):
            if x:
                for j, y in enumerate(v):
                    prod[i + j] = (prod[i + j] + x * y) % p
        for d in range(2 * k - 2, k - 1, -1):
            c = prod[d]
            if c:
                prod[d] = 0
                for i, l in enumerate(self.low):
                    prod[d - k + i] = (prod[d - k + i] + c * l) % p
        return tuple(prod[:k])

    def add(self, u, v):
        return tuple((x + y) % self.p for x, y in zip(u, v))

    def neg(self, u):
        return tuple((-x) % self.p for x in u)

    def inv(self, u):
        return self._inv[u]

    def scalar(self, c):
        return tuple([c % self.p] + [0] * (self.k - 1))

    def restrict(self, vec):
        """restriction of scalars: the GF(p)-span of {x^i vec : i < k}, as k rows over GF(p)."""
        rows, cur = [], list(vec)
        for _ in range(self.k):
            rows.append([c for e in cur for c in e])
            cur = [self.mul(self.gen, e) for e in cur]
        return rows


FIELD_POLYS = {2: (2, ()), 3: (3, ()), 4: (2, (1, 1)), 5: (5, ()), 7: (7, ()), 8: (2, (1, 1, 0)), 9: (3, (-1, 0)),
               11: (11, ()), 13: (13, ())}        # x^2 = x + 1 (GF(4)), x^3 = x + 1 (GF(8)), x^2 = -1 (GF(9))


def field(q):
    p, low = FIELD_POLYS[q]
    return GF(p, low)


def rank_mod(rows, p):
    A = [[x % p for x in r] for r in rows]
    rk, ncol = 0, (len(A[0]) if A else 0)
    for c in range(ncol):
        piv = next((i for i in range(rk, len(A)) if A[i][c]), None)
        if piv is None:
            continue
        A[rk], A[piv] = A[piv], A[rk]
        inv = pow(A[rk][c], p - 2, p)
        A[rk] = [(x * inv) % p for x in A[rk]]
        for i in range(len(A)):
            if i != rk and A[i][c]:
                t = A[i][c]
                A[i] = [(x - t * y) % p for x, y in zip(A[i], A[rk])]
        rk += 1
    return rk


def rank_vector(subs, p):
    """subs: n lists of row vectors over GF(p) (one list per element); returns the rank function on the 2^n subsets."""
    n = len(subs)
    r = np.zeros(1 << n, dtype=np.int64)
    for X in range(1, 1 << n):
        rows = [v for i in range(n) if X >> i & 1 for v in subs[i]]
        r[X] = rank_mod(rows, p) if rows else 0
    return r


def rank_vector_field(vectors, q):
    """rank function (over GF(q)) of one vector per element, through restriction of scalars to the prime field."""
    F = field(q)
    vecs = [[F.scalar(x) if isinstance(x, int) else x for x in v] for v in vectors]
    r = rank_vector([F.restrict(v) for v in vecs], F.p)
    assert (r % F.k == 0).all()
    return r // F.k


def T_config(M):
    """T(t), t = M - 1: A_i = e_i, B_i = c - e_i (integer vectors)."""
    return ([[1 if k == i else 0 for k in range(M)] for i in range(M)]
            + [[0 if k == i else 1 for k in range(M)] for i in range(M)])


# ------------------------------------------------------------------------------------------------ groups Z_N
def _egcd(a, b):
    if b == 0:
        return (abs(a), 1 if a >= 0 else -1, 0)
    g, x, y = _egcd(b, a % b)
    return (g, y, x - (a // b) * y)


def subgroup_order(vectors, mods):
    """|<vectors>| in Z_{m_1} x ... x Z_{m_k} = prod(m) / [Z^k : <vectors> + diag(m) Z^k] (triangular lattice basis)."""
    k = len(mods)
    basis = [None] * k
    cols = [list(v) for v in vectors] + [[m if j == i else 0 for j in range(k)] for i, m in enumerate(mods)]
    for v in cols:
        v = list(v)
        for i in range(k):
            if v[i] == 0:
                continue
            if basis[i] is None:
                basis[i] = v
                break
            x0, y0 = basis[i][i], v[i]
            g, s, u = _egcd(x0, y0)
            nb = [s * c1 + u * c2 for c1, c2 in zip(basis[i], v)]
            v = [(x0 // g) * c2 - (y0 // g) * c1 for c1, c2 in zip(basis[i], v)]
            basis[i] = nb if nb[i] > 0 else [-z for z in nb]
            v = [0] * (i + 1) + [z % m for z, m in zip(v[i + 1:], mods[i + 1:])]
    index = 1
    for i in range(k):
        index *= basis[i][i]
    return math.prod(mods) // index


def zn_values(M, N, base):
    """(LB*, LA*) on T(M - 1) over Z_N, in units of log base (exact as long as the orders are powers of base)."""
    LB, LA = lemma4_star(M)
    cfg = T_config(M)
    n = 2 * M
    r = np.zeros(1 << n, dtype=np.float64)
    for X in range(1, 1 << n):
        r[X] = math.log(subgroup_order([cfg[i] for i in range(n) if X >> i & 1], [N] * M), base)
    return round(float(LB @ r), 9), round(float(LA @ r), 9)


# ------------------------------------------------------------------------------------------------ Theorem 17 witnesses
def doubled_code_entropies(vectors, p):
    """The CSS state of [G | G] over GF(p) (party i: column i; purifier: all copies).  Entropies S(X) of every set X
    of the n parties and the purifier (bit n), from S(X) = r(X) + r(X^c) - r(E) for the code's column matroid."""
    n = len(vectors)
    cols = [[v] for v in vectors]
    def r(X):                                   # X over parties 0..n-1 and the purifier n
        rows = [v for i in range(n) if X >> i & 1 for v in cols[i]]
        if X >> n & 1:
            rows += [v for i in range(n) for v in cols[i]]
        return rank_mod(rows, p) if rows else 0
    full = (1 << (n + 1)) - 1
    rE = r(full)
    return np.array([r(X) + r(full ^ X) - rE for X in range(1 << (n + 1))], dtype=np.int64)


def theorem17_witness(t, p_state, form):
    """balanced form value on the doubled-code state of T(t) over GF(p_state): h_S(X) = S(X) + sum_i S(i)."""
    M = t + 1
    n = 2 * M
    S = doubled_code_entropies(T_config(M), p_state)
    h = np.array([S[X] + sum(S[1 << i] for i in range(n) if X >> i & 1) for X in range(1 << n)], dtype=np.int64)
    r = rank_vector([[v] for v in T_config(M)], p_state)
    modular_ok = np.array_equal(h - r, np.array([sum(S[1 << i] for i in range(n) if X >> i & 1) for X in range(1 << n)]))
    LB, LA = lemma4_star(M)
    f = LB if form == 'LB' else LA
    _c, fb = balanced(f, n)
    return int(fb @ h), int(f @ r), modular_ok


# ------------------------------------------------------------------------------------------------ Proposition 18
def _m8(st):
    return sum(1 << NAMES8.index(ch) for ch in st)


H_SETS = [_m8('BCDW'), _m8('ACDX'), _m8('ABDY'), _m8('ABCZ')]
A_SETS = [255 ^ h for h in H_SETS]
TYPE2 = [_m8(s) for s in ('ABWX', 'ACWY', 'ADWZ', 'BCXY', 'BDXZ', 'CDYZ')]


def pi(r):
    N = len(r)
    return np.array([r[X] + r[(N - 1) ^ X] - r[N - 1] for X in range(N)], dtype=np.int64)


def dual(r):
    N = len(r)
    card = np.array([bin(X).count('1') for X in range(N)])
    return card + r[(N - 1) ^ np.arange(N)] - r[N - 1]


def rho_config(F, rho):
    """v_j = e_j; w_m = (1, .., 1) with t_m = 1 + 1/rho_m in place m (section 4.14)."""
    t = [F.add(F.one, F.inv(x)) for x in rho]
    e = [[F.one if i == j else F.zero for i in range(4)] for j in range(4)]
    w = [[t[j] if i == j else F.one for i in range(4)] for j in range(4)]
    return e + w


def n1_config(F, tau):
    """v_1, v_2, v_3 = e_1, e_2, e_3, v_4 = -(v_1 + v_2 + v_3), w_m = s_m v_m + e_4 with s_m = 1/tau_m."""
    e = [[F.one if i == j else F.zero for i in range(4)] for j in range(3)]
    v = e + [[F.neg(F.one)] * 3 + [F.zero]]
    return v + [[F.mul(F.inv(tau[m]), c) for c in v[m][:3]] + [F.one] for m in range(4)]


def circuit_hyperplanes(r, a):
    """the 4-sets of deficiency a, if r is a times a rank-4 sparse paving matroid on 8 elements; else None."""
    for X in range(1, 256):
        c = bin(X).count('1')
        if c != 4 and r[X] != min(c, 4) * a:
            return None
    dep = {X for X in range(256) if bin(X).count('1') == 4 and r[X] != 4 * a}
    return dep if all(r[X] == 3 * a for X in dep) else None


def permute8(r, perm):
    out = np.zeros(256, dtype=np.int64)
    for X in range(256):
        out[sum(1 << perm[i] for i in range(8) if X >> i & 1)] = r[X]
    return out


def proposition18():
    o = {}
    rT8 = rank_vector_field(T_config(4), 3)
    rT8p = rank_vector_field(T_config(4), 5)
    S8, S8p = pi(rT8), pi(rT8p)
    o['one_cut'] = [(''.join(NAMES8[i] for i in range(8) if X >> i & 1), int(S8[X]), int(S8p[X]))
                    for X in range(256) if S8[X] != S8p[X]]
    # (i) GF(9) and GF(3)
    F9 = field(9)
    m1 = F9.neg(F9.one)
    sub9 = [F9.restrict(v) for v in rho_config(F9, [m1, m1, m1, F9.gen])]
    g9 = rank_vector(sub9, 3)                                  # twice the GF(9) rank function
    cfg3 = rho_config(F9, [m1, m1, m1, F9.one])
    assert all(all(c[1] == 0 for c in v) for v in cfg3)
    r3 = rank_vector([[[c[0] for c in v]] for v in cfg3], 3)
    o['gf9_N_H123'] = circuit_hyperplanes(g9, 2) == set(TYPE2 + H_SETS[:3])
    o['gf3_N_H123_A123'] = circuit_hyperplanes(r3, 1) == set(TYPE2 + H_SETS[:3] + A_SETS[:3])
    g = np.zeros(256, dtype=np.int64)
    for miss in range(4):                     # the index transposition (miss, 4) on v and w simultaneously
        sig = [0, 1, 2, 3]
        sig[3], sig[miss] = sig[miss], sig[3]
        perm = sig + [4 + s for s in sig]
        g += permute8(g9, perm) + permute8(r3, perm)
    o['qutrit_12'] = np.array_equal(pi(g), 12 * S8p) and all(g[1 << i] == 12 for i in range(8))
    # the same as an explicit GF(3) family (columns stacked), to see the state: 48 x 96 generator
    # (ii) GF(4)
    F4 = field(4)
    g4 = rank_vector([F4.restrict(v) for v in rho_config(F4, [F4.one, F4.one, F4.one, F4.gen])], 2)
    o['gf4_N_H123_A4'] = circuit_hyperplanes(g4, 2) == set(TYPE2 + H_SETS[:3] + [A_SETS[3]])
    o['qubit_2'] = np.array_equal(pi(g4), 2 * S8p)
    # (iii) GF(8) and GF(2)
    F8 = field(8)
    x = F8.gen
    x2 = F8.mul(x, x)
    tau = [F8.one, x, x2, F8.add(F8.add(F8.one, x), x2)]
    g8 = rank_vector([F8.restrict(v) for v in n1_config(F8, tau)], 2)     # 3 r_N1
    o['gf8_N1'] = circuit_hyperplanes(g8, 3) == set(TYPE2 + [_m8('ABCD'), _m8('WXYZ')])
    rAG = rank_vector_field(T_config(4), 2)
    o['gf2_AG'] = circuit_hyperplanes(rAG, 1) == set(TYPE2 + H_SETS + A_SETS)
    o['T8_plus_dual'] = np.array_equal(rT8 + dual(rT8), g8 // 3 + rAG)
    o['qubit_6'] = np.array_equal(pi(g8 + 3 * rAG), 6 * S8)
    hT8, hT8p = rT8 + dual(rT8), rT8p + dual(rT8p)
    # (i') h_T8' two-fold over GF(3): v_j = block j, w_m = I in every block but tau_m = I + rho_m^{-1} in block m
    g3 = rank_vector(two_fold_t8p_gf3(), 3)
    o['qutrit_2'] = np.array_equal(g3, hT8p) and np.array_equal(pi(g3), 2 * S8p)
    # (iv) h_T8 two-fold over Q from an integer family; reductions mod every prime p >= 5
    els = two_fold_t8_integer()
    o['rational_two_fold_T8'] = all(rank_Q([v for i in range(8) if X >> i & 1 for v in els[i]]) == hT8[X]
                                    for X in range(1, 256))
    o['bad_primes_T8'] = sorted(bad_primes(els))
    o['two_fold_T8_mod_p'] = {p: bool(np.array_equal(rank_vector(els, p), hT8)) for p in (2, 3, 5, 7, 11, 13, 17)}
    o['qudit_2_T8'] = np.array_equal(pi(rank_vector(els, 5)), 2 * S8)
    return o


TAU_T8P = [[[0, 1], [1, 0]], [[0, 1], [2, 1]], [[0, 1], [2, 1]], [[0, 2], [2, 0]]]   # rho_1..rho_4 over GF(3)


def two_fold_t8p_gf3():
    """eight 2-dimensional subspaces of GF(3)^8 with rank function h_T8' (rows = spanning vectors)."""
    p = 3
    out = []
    I = [[1, 0], [0, 1]]
    Z = [[0, 0], [0, 0]]
    def inv2(M):
        a, b, c, d = M[0][0], M[0][1], M[1][0], M[1][1]
        di = pow((a * d - b * c) % p, p - 2, p)
        return [[d * di % p, -b * di % p], [-c * di % p, a * di % p]]
    tau = [[[(I[i][j] + inv2(r)[i][j]) % p for j in range(2)] for i in range(2)] for r in TAU_T8P]
    def element(Ms):             # 8 x 2 block column -> two spanning vectors of length 8
        return [[Ms[k][i][j] for k in range(4) for i in range(2)] for j in range(2)]
    for j in range(4):
        out.append(element([I if k == j else Z for k in range(4)]))
    for m in range(4):
        out.append(element([tau[m] if k == m else I for k in range(4)]))
    return out


def two_fold_t8_integer():
    """eight 2-dimensional subspaces of Q^8 (integer spanning vectors) with rank function h_T8 over Q and over GF(p)
    for every prime p >= 5.  Coordinates: blocks A, B, C, W; U_X = [I; I; 0; I], U_Y = [I; 0; I; I],
    U_Z = [I; J; BJ; R4], U_D = [S1; I; B; S4] with J = [[1,1],[0,1]], B = [[0,-1],[1,1]], R4 = [[0,0],[4,3]],
    S1 = [[2,-1],[-4,2]], S4 = [[1,0],[0,0]]  (section 4.14)."""
    I = [[1, 0], [0, 1]]
    Z = [[0, 0], [0, 0]]
    J = [[1, 1], [0, 1]]
    B = [[0, -1], [1, 1]]
    BJ = [[0, -1], [1, 2]]
    R4 = [[0, 0], [4, 3]]
    S1 = [[2, -1], [-4, 2]]
    S4 = [[1, 0], [0, 0]]
    def element(Ms):
        return [[Ms[k][i][j] for k in range(4) for i in range(2)] for j in range(2)]
    U = {'A': element([I, Z, Z, Z]), 'B': element([Z, I, Z, Z]), 'C': element([Z, Z, I, Z]), 'W': element([Z, Z, Z, I]),
         'X': element([I, I, Z, I]), 'Y': element([I, Z, I, I]), 'Z': element([I, J, BJ, R4]),
         'D': element([S1, I, B, S4])}
    return [U[c] for c in NAMES8]


def rank_Q(rows):
    """rank over Q of integer row vectors (exact, with fractions)."""
    from fractions import Fraction
    A = [[Fraction(x) for x in r] for r in rows]
    rk = 0
    ncol = len(A[0]) if A else 0
    for c in range(ncol):
        piv = next((i for i in range(rk, len(A)) if A[i][c] != 0), None)
        if piv is None:
            continue
        A[rk], A[piv] = A[piv], A[rk]
        for i in range(len(A)):
            if i != rk and A[i][c] != 0:
                t = A[i][c] / A[rk][c]
                A[i] = [x - t * y for x, y in zip(A[i], A[rk])]
        rk += 1
    return rk


def bad_primes(els):
    """primes dividing a nonzero invariant factor of some set of the integer family: away from them the reduction mod p
    keeps every rank (the rank mod p is the number of invariant factors not divisible by p)."""
    import sympy
    from sympy import Matrix, ZZ
    from sympy.matrices.normalforms import smith_normal_form
    bad = set()
    for X in range(1, 256):
        M = Matrix([v for i in range(8) if X >> i & 1 for v in els[i]])
        S = smith_normal_form(M, domain=ZZ)
        for k in range(min(S.shape)):
            if S[k, k] != 0:
                bad |= {int(q) for q in sympy.factorint(abs(int(S[k, k])))}
    return bad


# ------------------------------------------------------------------------------------------------ DFZ arXiv:1401.2507
def _f8_single(st, k=1):
    return _e(256, _m8(st), k)


def _f8_cond(x, y, k=1):
    v = np.zeros(256, dtype=np.int64)
    v[_m8(x) | _m8(y)] += k
    if y:
        v[_m8(y)] -= k
    return v


_Q = _f8_single('A') + _f8_single('B') + _f8_single('C') + _f8_single('D') - _f8_single('ABCD')


def dfz41():
    """Theorem 4.1 (p. 31): H(A) <= 9H(Z) + 8H(Y) + 5H(X) + 6H(W) - 4H(D) - 12H(C) - 11H(B) - H(A) + 19H(Z|A,B,C)
    + 17H(Y|A,B,D) + 13H(X|A,C,D) + 11H(W|B,C,D) + H(A|W,X,Y,Z) + H(A|B,W,X) + 7H(B|D,X,Z) + 4H(B|C,X,Y)
    + 7H(C|D,Y,Z) + 5H(C|A,W,Y) + 4H(D|A,W,Z) + 29(H(A) + H(B) + H(C) + H(D) - H(A,B,C,D)); characteristic 3."""
    f = sum(_f8_single(s, k) for s, k in (('Z', 9), ('Y', 8), ('X', 5), ('W', 6), ('D', -4), ('C', -12), ('B', -11),
                                            ('A', -1)))
    f = f + sum(_f8_cond(x, y, k) for x, y, k in (('Z', 'ABC', 19), ('Y', 'ABD', 17), ('X', 'ACD', 13), ('W', 'BCD', 11),
                                                  ('A', 'WXYZ', 1), ('A', 'BWX', 1), ('B', 'DXZ', 7), ('B', 'CXY', 4),
                                                  ('C', 'DYZ', 7), ('C', 'AWY', 5), ('D', 'AWZ', 4)))
    return f + 29 * _Q - _f8_single('A')


def dfz41_functional():
    """second transcription, typed separately as a functional and evaluated on the unit vectors."""
    def F(h):
        Hs = lambda *s: h(frozenset(''.join(s)))
        Hc = lambda x, *c: Hs(x, *c) - Hs(*c)
        rhs = (9 * Hs('Z') + 8 * Hs('Y') + 5 * Hs('X') + 6 * Hs('W') - 4 * Hs('D') - 12 * Hs('C') - 11 * Hs('B') - Hs('A')
               + 19 * Hc('Z', 'A', 'B', 'C') + 17 * Hc('Y', 'A', 'B', 'D') + 13 * Hc('X', 'A', 'C', 'D')
               + 11 * Hc('W', 'B', 'C', 'D') + Hc('A', 'W', 'X', 'Y', 'Z') + Hc('A', 'B', 'W', 'X')
               + 7 * Hc('B', 'D', 'X', 'Z') + 4 * Hc('B', 'C', 'X', 'Y') + 7 * Hc('C', 'D', 'Y', 'Z')
               + 5 * Hc('C', 'A', 'W', 'Y') + 4 * Hc('D', 'A', 'W', 'Z')
               + 29 * (Hs('A') + Hs('B') + Hs('C') + Hs('D') - Hs('A', 'B', 'C', 'D')))
        return rhs - Hs('A')
    f = np.zeros(256, dtype=np.int64)
    for X in range(1, 256):
        S = frozenset(NAMES8[i] for i in range(8) if X >> i & 1)
        f[X] = F(lambda T, S=S: 1 if T == S else 0)
    return f


def dfz31_printed():
    """Theorem 3.1 as printed (p. 11); false over every field (section 4.5)."""
    f = sum(_f8_single(s, k) for s, k in (('Z', 8), ('Y', 29), ('X', 3), ('W', 8), ('D', -6), ('C', -17), ('B', -8),
                                            ('A', -17)))
    f = f + sum(_f8_cond(x, y, k) for x, y, k in (('Z', 'ABC', 55), ('Y', 'WXZ', 35), ('X', 'ACD', 50), ('W', 'BCD', 49),
                                                  ('A', 'BDY', 18), ('B', 'DXZ', 7), ('B', 'AWX', 1), ('C', 'DYZ', 7),
                                                  ('C', 'BXY', 7), ('C', 'AWY', 3), ('D', 'AWZ', 6)))
    return f + 49 * _Q - _f8_single('A')


def dfz31_reexpanded():
    """Theorem 3.1 re-expanded from the proof: the bounds Delta_A, Delta_B, Delta_B^, Delta_C, Delta_C^, Delta_C~,
    Delta_D (pp. 13-16), Delta_A* (eq. (63), p. 19) and the first display of p. 25."""
    H, C = _f8_single, _f8_cond
    dA = C('Z', 'ABC') + C('W', 'BCD') + C('X', 'ACD') + C('Y', 'WXZ') + C('A', 'BDY') + _Q
    dB = C('Z', 'ABC') + C('X', 'ACD') + C('B', 'DXZ') + _Q
    dBh = C('W', 'BCD') + C('X', 'ACD') + C('B', 'AWX') + _Q
    dC = 2 * C('Z', 'ABC') + C('W', 'BCD') + C('X', 'ACD') + C('Y', 'WXZ') + C('C', 'DYZ') + _Q
    dCh = C('Z', 'ABC') + C('W', 'BCD') + 2 * C('X', 'ACD') + C('Y', 'WXZ') + C('C', 'BXY') + _Q
    dCt = C('Z', 'ABC') + 2 * C('W', 'BCD') + C('X', 'ACD') + C('Y', 'WXZ') + C('C', 'AWY') + _Q
    dD = C('Z', 'ABC') + C('W', 'BCD') + C('D', 'AWZ') + _Q
    dAs = H('W') + 4 * H('Y') + H('Z') - 2 * H('A') - H('B') - 2 * H('C') - H('D') + 3 * dA + dB + dC + dCt + dD
    rhs = (dAs + H('Z') - H('C') + H('Y') - H('A') + dC + dAs
           + H('W') - H('C') + H('Y') - H('A') + dCt + dAs
           + H('X') - H('C') + H('Y') - H('A') + dCh + dAs
           + H('Z') - H('B') + H('X') - H('C') + H('Y') - H('A') + dAs + dB + dCh
           + H('W') - H('B') + H('X') - H('C') + H('Y') - H('A') + dAs + dBh + dCh)
    return rhs - H('A')


def exhaustive_2dim(f, q, chunk=1_000_000):
    """all arrangements of 8 subspaces of GF(q)^2 (0, the q + 1 lines, the plane); returns (count, #negative, min)."""
    L = q + 1
    Sb = np.array([0] + [1 << i for i in range(L)] + [(1 << L) - 1], dtype=np.int64)
    K = len(Sb)
    total, nneg, mn = 0, 0, None
    support = [X for X in range(1, 256) if f[X]]
    for start in range(0, K ** 8, chunk):
        idx = np.arange(start, min(K ** 8, start + chunk), dtype=np.int64)
        bits = Sb[np.stack([(idx // K ** i) % K for i in range(8)], axis=1)]
        vals = np.zeros(len(idx), dtype=np.int64)
        for X in support:
            u = np.zeros(len(idx), dtype=np.int64)
            for i in range(8):
                if X >> i & 1:
                    u |= bits[:, i]
            vals += f[X] * np.where(u == 0, 0, np.where((u & (u - 1)) == 0, 1, 2))
        total += len(idx)
        nneg += int((vals < 0).sum())
        mn = int(vals.min()) if mn is None else min(mn, int(vals.min()))
    return total, nneg, mn


def relabel_min(f, h):
    """min over the 8! relabellings sigma of sum_X f[X] h[sigma(X)]."""
    perms = np.array(list(itertools.permutations(range(8))), dtype=np.int64)
    P = np.zeros((len(perms), 256), dtype=np.int64)
    for i in range(8):
        P += np.outer(1 << perms[:, i], (np.arange(256) >> i) & 1)
    v = h[P] @ f
    return int(v.min()), int((v < 0).sum())


# ------------------------------------------------------------------------------------------------ subspace tables
def _rref(rows, q, ncol):
    A = [list(r) for r in rows]
    rk = 0
    for c in range(ncol):
        piv = next((i for i in range(rk, len(A)) if A[i][c] % q), None)
        if piv is None:
            continue
        A[rk], A[piv] = A[piv], A[rk]
        inv = pow(A[rk][c], q - 2, q)
        A[rk] = [(x * inv) % q for x in A[rk]]
        for i in range(len(A)):
            if i != rk and A[i][c] % q:
                t = A[i][c]
                A[i] = [(x - t * y) % q for x, y in zip(A[i], A[rk])]
        rk += 1
    return tuple(tuple(r) for r in A[:rk])


class Subspaces4:
    """all subspaces of GF(q)^4 (q prime), with join tables; ranks of 8-tuples by dynamic programming over subsets."""
    def __init__(self, q):
        self.q = q
        vecs = list(itertools.product(range(q), repeat=4))
        index, forms, frontier, jv = {(): 0}, [()], [()], {}
        while frontier:
            new = []
            for F in frontier:
                for vi, v in enumerate(vecs):
                    G = _rref(list(F) + [v], q, 4)
                    if G not in index:
                        index[G] = len(forms)
                        forms.append(G)
                        new.append(G)
                    jv[(index[F], vi)] = index[G]
            frontier = new
        n = len(forms)
        jvec = np.zeros((n, len(vecs)), dtype=np.int32)
        for (a, vi), b in jv.items():
            jvec[a, vi] = b
        vindex = {v: i for i, v in enumerate(vecs)}
        self.join = np.zeros((n, n), dtype=np.int32)
        for b, G in enumerate(forms):
            col = np.arange(n, dtype=np.int32)
            for row in G:
                col = jvec[col, vindex[row]]
            self.join[:, b] = col
        self.dim = np.array([len(F) for F in forms], dtype=np.int64)
        self.index = index
        self.by_dim = {d: np.array([i for i, F in enumerate(forms) if len(F) == d], dtype=np.int32) for d in range(5)}
        self.coord = [index[_rref([tuple(1 if k == i else 0 for k in range(4))], q, 4)] for i in range(4)]

    def ranks(self, idx):
        S = np.zeros((idx.shape[0], 256), dtype=np.int32)
        for X in range(1, 256):
            low = (X & -X).bit_length() - 1
            S[:, X] = self.join[S[:, X & (X - 1)], idx[:, low]]
        return self.dim[S]


def exhaustive_coordinate_lines(f, q, T=None, chunk=200_000):
    """A..D = the coordinate lines of GF(q)^4, W..Z every line: (count, #negative, min)."""
    T = T or Subspaces4(q)
    lines = T.by_dim[1]
    L = len(lines)
    total = L ** 4
    worst, neg = None, 0
    for start in range(0, total, chunk):
        k = np.arange(start, min(total, start + chunk))
        idx = np.zeros((len(k), 8), dtype=np.int32)
        idx[:, :4] = T.coord
        for j in range(4):
            idx[:, 4 + j] = lines[(k // L ** j) % L]
        v = T.ranks(idx) @ f
        worst = int(v.min()) if worst is None else min(worst, int(v.min()))
        neg += int((v < 0).sum())
    return total, neg, worst


def random_arrangements(f, q, count, seed, T=None, batch=50_000):
    """random 8-tuples of subspaces of GF(q)^4 (dimension 1 with probability 0.6, 2: 0.2, 0: 0.1, 3: 0.07, 4: 0.03),
    and coordinate lines A..D with random lines W..Z (half of each batch)."""
    T = T or Subspaces4(q)
    rng = np.random.default_rng(seed)
    worst = None
    done = 0
    while done < count:
        bsz = min(batch, count - done)
        dims = rng.choice(5, size=(bsz, 8), p=[0.1, 0.6, 0.2, 0.07, 0.03])
        idx = np.zeros((bsz, 8), dtype=np.int32)
        for d in range(5):
            pool, mask = T.by_dim[d], dims == d
            idx[mask] = pool[rng.integers(0, len(pool), int(mask.sum()))]
        half = bsz // 2
        idx[:half, :4] = T.coord
        lines = T.by_dim[1]
        idx[:half, 4:] = lines[rng.integers(0, len(lines), (half, 4))]
        v = T.ranks(idx) @ f
        worst = int(v.min()) if worst is None else min(worst, int(v.min()))
        done += bsz
    return worst


# ------------------------------------------------------------------------------------------------ data
def build_data():
    LB4, LA4 = lemma4_star(4)
    LB6, LA6 = lemma4_star(6)
    return {'LB4': LB4, 'LA4': LA4, 'LB6': LB6, 'LA6': LA6, 'DFZ41': dfz41(), 'DFZ31_printed': dfz31_printed(),
            'DFZ31_reexpanded': dfz31_reexpanded()}


def content_sha(path=DATA):
    Z = np.load(path)
    hs = hashlib.sha256()
    for k in sorted(Z.files):
        hs.update(k.encode())
        hs.update(np.ascontiguousarray(Z[k], dtype='<i8').tobytes())
    return hs.hexdigest()


# ------------------------------------------------------------------------------------------------ checks
def check(full_check=False, log=print):
    o = {}
    t0 = time.time()
    # Lemma 4*: values on T(t) over fields, M = 4, 5, 6
    vals = {}
    for M in (4, 5, 6):
        LB, LA = lemma4_star(M)
        t = M - 1
        for q in (2, 3, 4, 5, 7, 8, 9, 11):
            r = rank_vector_field(T_config(M), q)
            ch = FIELD_POLYS[q][0]
            vals[(M, q)] = (int(LB @ r), int(LA @ r), int(r[((1 << M) - 1) << M]))
            exp = (-1 if t % ch == 0 else 0, 0 if t % ch == 0 else -1, M - 1 if t % ch == 0 else M)
            if vals[(M, q)] != exp:
                log(f'unexpected value M={M} q={q}: {vals[(M, q)]} != {exp}')
    o['values_ok'] = all(vals[(M, q)] == (-1 if (M - 1) % FIELD_POLYS[q][0] == 0 else 0,
                                           0 if (M - 1) % FIELD_POLYS[q][0] == 0 else -1,
                                           M - 1 if (M - 1) % FIELD_POLYS[q][0] == 0 else M) for (M, q) in vals)
    o['values'] = {f'M={M},GF({q})': v[:2] for (M, q), v in vals.items()}
    # every choice of j(i) gives the same values on T(3)
    r3, r5 = rank_vector_field(T_config(4), 3), rank_vector_field(T_config(4), 5)
    jv = set()
    for jm in itertools.product(range(4), repeat=4):
        if all(jm[i] != i for i in range(4)):
            LBj, LAj = lemma4_star(4, list(jm))
            jv.add((int(LBj @ r3), int(LBj @ r5), int(LAj @ r3), int(LAj @ r5)))
    o['jmap_values'] = sorted(jv)
    # M = 3: no gain (L = 0 on T(2))
    LB3, LA3 = lemma4_star(3)
    o['M3_values'] = [(int(LB3 @ rank_vector_field(T_config(3), q)), int(LA3 @ rank_vector_field(T_config(3), q)))
                      for q in (2, 3)]
    # balancing coefficients
    o['balancing'] = {M: (balanced(lemma4_star(M)[0], 2 * M)[0].tolist(), balanced(lemma4_star(M)[1], 2 * M)[0].tolist())
                      for M in (4, 6)}
    o['balancing_nonneg'] = all(min(a) >= 0 and min(b) >= 0 for a, b in o['balancing'].values())
    # Theorem 17: the doubled-code witnesses for t = 3 (8 parties) and t = 5 (12 parties, LB only over GF(5))
    w = {}
    for t, cases in ((3, [('LB', 3), ('LB', 2), ('LB', 5), ('LA', 2), ('LA', 5), ('LA', 7), ('LA', 3)]),
                     (5, [('LB', 5), ('LA', 3)])):
        for form, q in cases:
            w[(t, form, q)] = theorem17_witness(t, q, form)
    o['witness'] = {f't={t},{form},GF({q})': v[0] for (t, form, q), v in w.items()}
    o['witness_ok'] = (all(v[2] for v in w.values())
                       and all(v[0] == (-1 if ((form == 'LB') == (t % q == 0)) else 0) for (t, form, q), v in w.items())
                       and all(v[0] == v[1] for v in w.values()))
    # groups: values on T(t) over Z_N (units of log p): (LA*) = -log(N / gcd(t, N)), (LB*) = -log gcd(t, N)
    zn = {}
    for (M, N, base) in ((4, 9, 3), (4, 4, 2), (5, 8, 2), (4, 27, 3), (6, 25, 5)):
        zn[(M, N)] = zn_values(M, N, base)
    o['zn'] = {f'T({M - 1}) over Z_{N}': v for (M, N), v in zn.items()}
    o['zn_ok'] = all(abs(zn[(M, N)][1] + math.log(N // math.gcd(M - 1, N), base)) < 1e-9
                     and abs(zn[(M, N)][0] + math.log(math.gcd(M - 1, N), base)) < 1e-9
                     for (M, N, base) in ((4, 9, 3), (4, 4, 2), (5, 8, 2), (4, 27, 3), (6, 25, 5)))
    # Proposition 18
    o['prop18'] = proposition18()
    P18 = o['prop18']
    o['prop18_ok'] = (P18['one_cut'] == [('ABCD', 3, 4), ('WXYZ', 3, 4)]
                      and all(bool(P18[k]) for k in ('gf9_N_H123', 'gf3_N_H123_A123', 'qutrit_12', 'gf4_N_H123_A4',
                                                       'qubit_2', 'gf8_N1', 'gf2_AG', 'T8_plus_dual', 'qubit_6',
                                                       'qutrit_2', 'rational_two_fold_T8', 'qudit_2_T8'))
                      and P18['bad_primes_T8'] == [2, 3]
                      and P18['two_fold_T8_mod_p'] == {2: False, 3: False, 5: True, 7: True, 11: True, 13: True,
                                                       17: True})
    # DFZ
    f41, f41b = dfz41(), dfz41_functional()
    o['dfz41_transcriptions_agree'] = bool(np.array_equal(f41, f41b))
    o['dfz41_thm42'] = {q: int(f41 @ rank_vector_field(T_config(4), q)) for q in (2, 3, 4, 5, 7, 8, 9, 11, 13)}
    c41, fb41 = balanced(f41, 8)
    o['dfz41_balancing'] = c41.tolist()
    hT8p = rank_vector_field(T_config(4), 5)
    hT8p = hT8p + dual(hT8p)
    hT8 = rank_vector_field(T_config(4), 3)
    hT8 = hT8 + dual(hT8)
    hAG = rank_vector_field(T_config(4), 2)
    hAG = hAG + dual(hAG)
    o['dfz41_relabel_min'] = {'h_T8\'': relabel_min(fb41, hT8p), 'h_T8': relabel_min(fb41, hT8), 'h_AG': relabel_min(fb41, hAG)}
    p31, e31 = dfz31_printed(), dfz31_reexpanded()
    o['dfz31_difference'] = terms(e31 - p31, 8, list(NAMES8))
    h_cx = np.array([1 if (X & (_m8('C') | _m8('X'))) else 0 for X in range(256)])
    o['dfz31_cx'] = (int(p31 @ h_cx), int(e31 @ h_cx))
    o['dfz31_T8'] = {q: int(e31 @ rank_vector_field(T_config(4), q)) for q in (2, 3, 5, 7)}
    # stored data
    D = build_data()
    if os.path.exists(DATA):
        Z = np.load(DATA)
        o['data_ok'] = sorted(Z.files) == sorted(D) and all(np.array_equal(Z[k], D[k]) for k in D)
    else:
        o['data_ok'] = False
    o['sha'] = {k: form_sha(v) for k, v in D.items()}
    o['time_fast'] = round(time.time() - t0, 1)
    if full_check:
        t1 = time.time()
        LB4, LA4 = D['LB4'], D['LA4']
        ex = {}
        for q in (2, 3, 5):
            ex[f'LB*,GF({q})^2'] = exhaustive_2dim(LB4, q)
            ex[f'LA*,GF({q})^2'] = exhaustive_2dim(LA4, q)
        ex['DFZ41,GF(3)^2'] = exhaustive_2dim(f41, 3)
        ex['DFZ31 re-expanded,GF(2)^2'] = exhaustive_2dim(e31, 2)
        ex['DFZ31 re-expanded,GF(5)^2'] = exhaustive_2dim(e31, 5)
        ex['DFZ31 printed,GF(2)^2'] = exhaustive_2dim(p31, 2)
        o['exhaustive_2dim'] = ex
        T2, T3, T5 = Subspaces4(2), Subspaces4(3), Subspaces4(5)
        o['coordinate_lines'] = {'LB*,GF(2)': exhaustive_coordinate_lines(LB4, 2, T2),
                                 'LA*,GF(3)': exhaustive_coordinate_lines(LA4, 3, T3),
                                 'LA*,GF(2) (not claimed)': exhaustive_coordinate_lines(LA4, 2, T2),
                                 'LB*,GF(3) (not claimed)': exhaustive_coordinate_lines(LB4, 3, T3)}
        o['random'] = {'LB*,GF(2)': random_arrangements(LB4, 2, 400_000, 4201, T2),
                       'LB*,GF(5)': random_arrangements(LB4, 5, 400_000, 4202, T5),
                       'LA*,GF(3)': random_arrangements(LA4, 3, 400_000, 4203, T3),
                       'DFZ41,GF(3)': random_arrangements(f41, 3, 200_000, 4204, T3),
                       'DFZ31 re-expanded,GF(5)': random_arrangements(e31, 5, 200_000, 4205, T5)}
        # M = 8 (t = 7, 16 variables): values over GF(7) and GF(3)
        LB8, LA8 = lemma4_star(8)
        o['M8'] = {q: (int(LB8 @ rank_vector_field(T_config(8), q)), int(LA8 @ rank_vector_field(T_config(8), q)))
                   for q in (7, 3)}
        o['time_full'] = round(time.time() - t1, 1)
    return o


EXPECT = {
    'jmap_values': [(-1, 0, 0, -1)],
    'M3_values': [(2, 1), (3, 0)],
    'dfz41_thm42': {2: -1, 3: 0, 4: -1, 5: -1, 7: -1, 8: -1, 9: 0, 11: -1, 13: -1},
    'dfz41_relabel_min': {'h_T8\'': (42, 0), 'h_T8': (42, 0), 'h_AG': (-2, 1344)},
    'dfz31_difference': [('ACD', 4), ('AWY', -4), ('BXY', 4), ('ACDX', -4), ('ACWY', 4), ('BCXY', -4)],
    'dfz31_cx': (-4, 0),
    'dfz31_T8': {2: 34, 3: -1, 5: 34, 7: 34},
}


def check_ok(o, full_check=False):
    ok = (o['values_ok'] and o['witness_ok'] and o['zn_ok'] and o['prop18_ok'] and o['balancing_nonneg']
          and o['dfz41_transcriptions_agree'] and o['data_ok']
          and all(o[k] == v for k, v in EXPECT.items()))
    if full_check:
        ex = o['exhaustive_2dim']
        ok = ok and all(v[1] == 0 for k, v in ex.items() if 'printed' not in k) and ex['DFZ31 printed,GF(2)^2'][1] == 125
        cl = o['coordinate_lines']
        ok = ok and cl['LB*,GF(2)'][1] == 0 and cl['LA*,GF(3)'][1] == 0
        ok = ok and cl['LA*,GF(2) (not claimed)'][1] == 1 and cl['LB*,GF(3) (not claimed)'][1] == 8
        ok = ok and all(v == 0 for v in o['random'].values())
        ok = ok and o['M8'] == {7: (-1, 0), 3: (0, -1)}
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--check', action='store_true')
    ap.add_argument('--full', action='store_true')
    ap.add_argument('--forms', action='store_true')
    ap.add_argument('--write-data', action='store_true')
    a = ap.parse_args()
    if a.write_data:
        np.savez(DATA, **build_data())
        print('wrote', DATA, 'content sha256', content_sha())
    if a.forms:
        LB, LA = lemma4_star(4)
        for name, f in (('(LB*) t = 3, char != 3', LB), ('(LA*) t = 3, char 3', LA)):
            print(name + ':  0 <= ' + ' '.join('%+d h(%s)' % (c, s) for s, c in terms(f, 8, list(NAMES8))))
    if a.check:
        o = check(a.full)
        for k, v in o.items():
            print(f'{k}: {v}')
        print('OK' if check_ok(o, a.full) else 'FAILURE')


if __name__ == '__main__':
    main()
