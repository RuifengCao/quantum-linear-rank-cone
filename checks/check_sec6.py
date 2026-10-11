"""Checks for Sec. VI / Appendices E, F and H (six parties), independent of the repository code except where noted
(the transcriptions of DFZ (65)/(91), BKL (6.8)/(6.17) and Pena-Sarria Example 6 are taken from the repository;
Theorems 7 and 9 of Pena-Macias and Sarria (2021) with n = 2 are transcribed in this script).

Elements A, B, C, W, X, Y, Z = bits 0..6 = the points 100, 010, 001, 110, 101, 011, 111 of GF(2)^3.
"""
import csv, itertools, math, sys
from fractions import Fraction
import numpy as np
import sympy
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from common import REPO, rank_mod_p, rank_Q
sys.path.insert(0, REPO + '/scripts')

L = 'ABCWXYZ'
VEC = [(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0), (1, 0, 1), (0, 1, 1), (1, 1, 1)]
FULL = 127
LINES = ['ABW', 'ACX', 'BCY', 'AYZ', 'BXZ', 'CWZ', 'WXY']


def m(s):
    return sum(1 << L.index(c) for c in s)


def rank_fn(p):
    r = np.zeros(128, dtype=np.int64)
    for X in range(1, 128):
        r[X] = rank_mod_p([list(VEC[i]) for i in range(7) if X >> i & 1], p)
    return r


def dual(h):
    hs = np.zeros(128, dtype=np.int64)
    for X in range(1, 128):
        hs[X] = sum(h[1 << i] for i in range(7) if X >> i & 1) + h[FULL ^ X] - h[FULL]
    return hs


ok = True
rF = rank_fn(2)
rNF = rank_fn(3)
rFs, rNFs = dual(rF), dual(rNF)
hF, hNF = rF + rFs, rNF + rNFs
# values of h_F by orbit type
lines = {m(l) for l in LINES}
def kind(X):
    k = bin(X).count('1')
    if k == 3:
        return 'line' if X in lines else 'triple'
    if k == 4:
        return 'linecompl' if (FULL ^ X) in lines else '4set'
    return str(k)
vals = {}
for X in range(1, 128):
    vals.setdefault(kind(X), set()).add(int(hF[X]))
print('h_F by type:', {k: sorted(v) for k, v in vals.items()})

# Steane |0_L> (Nielsen-Chuang eq. 10.78 codewords) = CSS state of the simplex code, up to qubit order
nc = ['0000000', '1010101', '0110011', '1100110', '0001111', '1011010', '0111100', '1101001']
ncset = {tuple(int(c) for c in w) for w in nc}
# NC codewords are the row space of H with columns 1..7 in binary (bit strings 001..111, MSB first)
H = [[(j >> (2 - k)) & 1 for j in range(1, 8)] for k in range(3)]
rows = set()
for y in itertools.product(range(2), repeat=3):
    rows.add(tuple(sum(y[k] * H[k][j] for k in range(3)) % 2 for j in range(7)))
print('Nielsen-Chuang |0_L> codewords = row space of [binary 1..7]:', rows == ncset)
ok &= rows == ncset
# simplex code in our labelling: columns = VEC; entropy of its CSS state from state vector
code = set()
for y in itertools.product(range(2), repeat=3):
    code.add(tuple(sum(y[k] * VEC[j][k] for k in range(3)) % 2 for j in range(7)))
psi = np.zeros(128)
for c in code:
    psi[int(''.join(map(str, c)), 2)] = 1
psi /= np.linalg.norm(psi)
def ent(X):
    q = [j for j in range(7) if X >> j & 1]
    rest = [j for j in range(7) if j not in q]
    T = psi.reshape([2] * 7).transpose(q + rest).reshape(2 ** len(q), 2 ** len(rest))
    s = np.linalg.svd(T, compute_uv=False) ** 2
    s = s[s > 1e-12]
    return float(-(s * np.log2(s)).sum())
S = np.zeros(128)
for X in range(1, 127):
    S[X] = ent(X)
hS = np.array([S[X] + sum(S[1 << i] for i in range(7) if X >> i & 1) for X in range(128)])
print('state vector: h_S == r_F + r_F*:', np.allclose(hS[1:128], hF[1:128]), ' (h_S of all seven systems:', round(float(hS[127]), 9), ')')
ok &= np.allclose(hS[1:128], hF[1:128])

# ------------------------------------------------------------------ the eight forms
def load_csv(path):
    forms = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith('#') or line.startswith('ineq_id') or not line.strip():
                continue
            name, subset, coeff = line.strip().split(',')
            v = forms.setdefault(name, np.zeros(128, dtype=np.int64))
            v[m(subset)] += int(coeff)
    return forms


DFZ = load_csv(REPO + '/data/dfz_char_dependent.csv')
BKL = load_csv(REPO + '/data/bkl_char_dependent.csv')
print('DFZ forms:', list(DFZ), ' BKL forms:', list(BKL))


def lemma_forms_7(t=2, perm=('A', 'B', 'C', 'Y', 'X', 'W', 'Z')):
    """(LB) and (LA) with t = 2 as vectors over 128 masks; (A1,A2,A3,B1,B2,B3,C) = perm."""
    Mm = t + 1
    Am = [m(perm[i]) for i in range(Mm)]
    Bm = [m(perm[Mm + i]) for i in range(Mm)]
    Cm = m(perm[2 * Mm])
    def vec():
        return np.zeros(128, dtype=np.int64)
    def Hv(X):
        v = vec()
        if X:
            v[X] += 1
        return v
    def Hc(x, y):
        return Hv(x | y) - Hv(y)
    def I(x, y, z=0):
        return Hv(x | z) + Hv(y | z) - Hv(x | y | z) - Hv(z)
    allA = sum(Am)
    alpha = sum(Hv(a) for a in Am) - Hv(allA)
    gamma = Hc(Cm, allA)
    delta = sum(I(Cm, allA ^ Am[i]) for i in range(Mm))
    eps = sum(Hc(Cm, Am[i] | Bm[i]) for i in range(Mm))
    zeta = sum(Hc(Bm[i], allA ^ Am[i]) for i in range(Mm))
    eta = sum(Hc(Bm[i], Am[i] | Cm) for i in range(Mm))
    HB = Hv(sum(Bm))
    LB = HB + Mm * (gamma + Mm * alpha + delta + eps + zeta) - Mm * Hv(Cm)
    LA = (Mm - 1) * Hv(Cm) + Mm * alpha + Mm * gamma + zeta + eta - HB
    return LB, LA


LB2, LA2 = lemma_forms_7()
printed_LB2 = {'A': 9, 'B': 9, 'C': 9, 'ABC': -12, 'Z': 6, 'ABCZ': 3, 'ABZ': -3, 'ACZ': -3, 'BCZ': -3, 'AYZ': 3, 'BXZ': 3,
               'CWZ': 3, 'AY': -3, 'BX': -3, 'CW': -3, 'BCY': 3, 'ACX': 3, 'ABW': 3, 'WXY': 1}
pv = np.zeros(128, dtype=np.int64)
for s_, c in printed_LB2.items():
    pv[m(s_)] += c
print('(LB)_2 equals the printed 19-term form:', np.array_equal(pv, LB2))
ok &= np.array_equal(pv, LB2)

import s28_generalized_cubes as g28                        # Pena-Sarria transcription (repository)
a_ps, b_ps, nps = g28.ps_forms(2, 3)                       # bits: A1..A3 = 0..2, B1..B3 = 3..5, C = 6
# map to our letters with (A1,A2,A3,B1,B2,B3,C) = (A,B,C,Y,X,W,Z)
perm_bits = [m('A'), m('B'), m('C'), m('Y'), m('X'), m('W'), m('Z')]
def remap(v):
    out = np.zeros(128, dtype=np.int64)
    for X in range(128):
        if v[X]:
            out[sum(perm_bits[i] for i in range(7) if X >> i & 1)] += v[X]
    return out
forms = [('DFZ (65)', DFZ['dfz65_odd'], 'odd'), ('DFZ (91)', DFZ['dfz91_even'], 'even'),
         ('(LA)_2', LA2, 'even'), ('(LB)_2', LB2, 'odd'),
         ('Pena-Sarria (a), t=2', remap(np.asarray(a_ps)), 'even'), ('Pena-Sarria (b), t=2', remap(np.asarray(b_ps)), 'odd'),
         ('BKL (6.8)', BKL['bkl68_even'], 'even'), ('BKL (6.17) corrected', BKL['bkl617_odd_corrected'], 'odd')]


def elemental7():
    rows = []
    for i in range(7):
        v = np.zeros(128, dtype=np.int64); v[FULL] += 1; v[FULL ^ (1 << i)] -= 1; rows.append(v)
    for i, j in itertools.combinations(range(7), 2):
        rest = [k for k in range(7) if k not in (i, j)]
        for s in range(1 << 5):
            K = sum(1 << rest[b] for b in range(5) if s >> b & 1)
            v = np.zeros(128, dtype=np.int64)
            v[K | 1 << i] += 1; v[K | 1 << j] += 1; v[K | 1 << i | 1 << j] -= 1
            if K:
                v[K] -= 1
            rows.append(v)
    return np.array(rows)


E = elemental7()


def symmetrise(f):
    g = f.copy()
    for X in range(1, 128):
        c = int(f[X])
        if c:
            for P in range(7):
                if X >> P & 1:
                    g[1 << P] += c
            if FULL ^ X:
                g[FULL ^ X] += c
            g[FULL] -= c
    return g


from scipy.optimize import linprog
perms = list(itertools.permutations(range(7)))
def apply_perm(v, p):
    out = np.zeros(128, dtype=np.int64)
    for X in range(128):
        if v[X]:
            out[sum(1 << p[i] for i in range(7) if X >> i & 1)] += v[X]
    return out
def min_over_relabelings(f, h):
    best = None
    for p in perms:
        val = int(apply_perm(f, p) @ h)
        best = val if best is None else min(best, val)
    return best

for name, f, char in forms:
    g = symmetrise(f)
    res = linprog(np.ones(E.shape[0]), A_eq=E[:, 1:].T.astype(float), b_eq=g[1:].astype(float), bounds=(0, None), method='highs')
    assert res.status == 0, name
    supp = np.flatnonzero(res.x > 1e-9)
    sol, params = sympy.Matrix(E[supp][:, 1:].T.tolist()).gauss_jordan_solve(sympy.Matrix([int(x) for x in g[1:]]))
    sol = sol.subs({s_: 0 for s_ in params})
    fr = [Fraction(int(sympy.fraction(x)[0]), int(sympy.fraction(x)[1])) for x in sol]
    exact = all(x >= 0 for x in fr)
    chk = np.zeros(127, dtype=object)
    for k, x in zip(supp, fr):
        chk = chk + x * E[k][1:].astype(object)
    exact &= all(chk[i] == g[1 + i] for i in range(127))
    den = 1
    for x in fr:
        den = den * x.denominator // math.gcd(den, x.denominator)
    vF, vNF = int(f @ rF), int(f @ rNF)
    print('%-24s valid in %-4s char: f(r_F) = %3d, f(r_NF) = %3d; f+f.D Shannon certificate exact: %s, support %d, '
          'denominator %d' % (name, char, vF, vNF, exact, len([x for x in fr if x != 0]), den))
    ok &= exact
    ok &= (vF == -1 and vNF >= 0) if char == 'odd' else (vNF == -1 and vF >= 0)
# Pena-Macias and Sarria, J. Inf. Optim. Sci. 42 (2021) 345 (arXiv:1903.11587v2), Theorems 7 and 9 with n = 2,
# transcribed here from the arXiv text, under (A1,A2,A3,B1,B2,B3,C) = (A,B,C,Y,X,W,Z)
def _Hv(X):
    v = np.zeros(128, dtype=np.int64)
    if X:
        v[X] += 1
    return v
def _Hc(x, y):
    return _Hv(x | y) - _Hv(y)
def _I(x, y, z=0):
    return _Hv(x | z) + _Hv(y | z) - _Hv(x | y | z) - _Hv(z)
_A = [m('A'), m('B'), m('C')]; _B = [m('Y'), m('X'), m('W')]; _C = m('Z')
_allA, _allB = _A[0] | _A[1] | _A[2], _B[0] | _B[1] | _B[2]
ps9 = (_Hv(_allB) + 3 * (_Hc(_C, _allA) + sum(_I(_allA ^ _A[i], _C) for i in range(3)) + sum(_Hc(_C, _A[i] | _B[i]) for i in range(3))
       + 2 * _I(_A[0], _A[1]) + 3 * _I(_A[0] | _A[1], _A[2]) + sum(_Hc(_B[i], _allA ^ _A[i]) for i in range(3))) - 3 * _Hv(_C))
ps7 = (2 * _I(_allA, _C) + sum(_Hc(_B[i], _allA ^ _A[i]) for i in range(3)) + sum(_Hc(_B[i], _A[i] | _C) for i in range(3))
       + 2 * _I(_A[0], _A[1]) + 3 * (_I(_A[0] | _A[1], _A[2]) + _Hc(_C, _allA) + sum(_I(_allA ^ _A[i], _C) for i in range(3)))
       - _Hv(_allB))
print('Pena-Macias-Sarria Thm 9 (x3): Fano %d, non-Fano %d; Thm 7: Fano %d, non-Fano %d'
      % (int(ps9 @ rF), int(ps9 @ rNF), int(ps7 @ rF), int(ps7 @ rNF)))
ok &= int(ps9 @ rF) == -1 and int(ps9 @ rNF) >= 0 and int(ps7 @ rNF) == -1 and int(ps7 @ rF) >= 0
print('3 x Thm 9 == (LB)_2 - 3 I(A;B):', np.array_equal(ps9, LB2 - 3 * _I(m('A'), m('B'))))
ok &= np.array_equal(ps9, LB2 - 3 * _I(m('A'), m('B')))
for name, f in [('PMS Theorem 9 (x3), n=2', ps9), ('PMS Theorem 7, n=2', ps7)]:
    g = symmetrise(f)
    res = linprog(np.ones(E.shape[0]), A_eq=E[:, 1:].T.astype(float), b_eq=g[1:].astype(float), bounds=(0, None), method='highs')
    assert res.status == 0, name
    supp = np.flatnonzero(res.x > 1e-9)
    sol, params = sympy.Matrix(E[supp][:, 1:].T.tolist()).gauss_jordan_solve(sympy.Matrix([int(x) for x in g[1:]]))
    sol = sol.subs({s_: 0 for s_ in params})
    fr = [Fraction(int(sympy.fraction(x)[0]), int(sympy.fraction(x)[1])) for x in sol]
    chk = np.zeros(127, dtype=object)
    for k, x in zip(supp, fr):
        chk = chk + x * E[k][1:].astype(object)
    exact = all(x >= 0 for x in fr) and all(chk[i] == g[1 + i] for i in range(127))
    print('%-24s f+f.D Shannon certificate exact: %s, support %d' % (name, exact, len([x for x in fr if x != 0])))
    ok &= exact
# Theorem 7 of Pena-Macias-Sarria coincides with (6.8); 3 x Theorem 9 = (b) + 3 I(B;C); the identities of App. E
fd = dict((name, f) for name, f, _c in forms)
ids = [('PMS Theorem 7 (n=2) == (6.8)', np.array_equal(ps7, fd['BKL (6.8)'])),
       ('3 x PMS Theorem 9 == (b) + 3 I(B;C)', np.array_equal(ps9, fd['Pena-Sarria (b), t=2'] + 3 * _I(m('B'), m('C')))),
       ('(6.8) == (a) + I(B;C)', np.array_equal(fd['BKL (6.8)'], fd['Pena-Sarria (a), t=2'] + _I(m('B'), m('C')))),
       ('(b) == (LB)_2 - 3 I(A;B) - 3 I(B;C)', np.array_equal(fd['Pena-Sarria (b), t=2'],
                                                          LB2 - 3 * _I(m('A'), m('B')) - 3 * _I(m('B'), m('C'))))]
for name, val in ids:
    print('App. E identity:', name, val)
    ok &= val
# no two of the eight forms are related by a relabeling (compared after dividing by the gcd of the coefficients)
maps = np.array([[sum(1 << p_[i] for i in range(7) if X >> i & 1) for X in range(128)] for p_ in perms])
def normalized(v):
    g_ = int(np.gcd.reduce(np.abs(v[v != 0])))
    return v // g_
related = []
names8 = [nm for nm, _f, _c in forms]
for a_, b_ in itertools.combinations(range(8), 2):
    fa, fb = normalized(forms[a_][1]), normalized(forms[b_][1])
    for mp in maps:
        img = np.zeros(128, dtype=np.int64); img[mp] = fa
        if np.array_equal(img, fb):
            related.append((names8[a_], names8[b_])); break
print('pairs of the eight forms related by a relabeling (up to scaling):', related)
ok &= not related
# minimum over relabelings on h_F / h_NF for (LB)_2 and DFZ (65), (LA)_2 and DFZ (91)
mins = (min_over_relabelings(LB2, hF), min_over_relabelings(DFZ['dfz65_odd'], hF),
        min_over_relabelings(LA2, hNF), min_over_relabelings(DFZ['dfz91_even'], hNF))
print('min over 5040 relabelings: (LB)_2 on h_F:', mins[0], ' DFZ(65) on h_F:', mins[1],
      ' (LA)_2 on h_NF:', mins[2], ' DFZ(91) on h_NF:', mins[3])
ok &= mins == (17, 19, 6, 16)

# ------------------------------------------------------------------ the transfer and f9
P = lambda s: m(s)
lines_through = {i: [l for l in LINES if L[i] in l] for i in range(7)}
order = 'ABCWXYZ'


def form(d):
    v = np.array([Fraction(0)] * 128, dtype=object)
    for X, c in d.items():
        if X:
            v[X] += Fraction(c)
    return v


def H(X):
    return form({X: 1})


def I(x, y, z=0):
    return H(x | z) + H(y | z) - H(x | y | z) - (H(z) if z else form({}))


def mform(i, line):
    rest = m(line) ^ (1 << i)
    return I(1 << i, rest)


def u(i):
    Ls = lines_through[i]
    tot = form({})
    for a, b in itertools.combinations(range(3), 2):
        Q = (m(Ls[a]) | m(Ls[b])) ^ (1 << i)
        tot = tot + 2 * I(1 << i, Q)
    for a in range(3):
        tot = tot - mform(i, Ls[a])
    return tot / 3


def lbar(i):
    tot = form({})
    for l in lines_through[i]:
        tot = tot + mform(i, l)
    return tot / 3


def upper(X):
    k = bin(X).count('1')
    pts = [i for i in range(7) if X >> i & 1]
    if k == 1:
        return u(pts[0])
    name = ''.join(L[i] for i in pts)
    if any(set(name) == set(l) for l in LINES):
        line = [l for l in LINES if set(l) == set(name)][0]
        first = min(pts, key=lambda i: order.index(L[i]))
        return u(pts[0]) + u(pts[1]) + u(pts[2]) - mform(first, line)
    if X == m('ABCZ'):
        return u(L.index('A')) + u(L.index('B')) + u(L.index('C')) + (u(L.index('Y')) - mform(L.index('Y'), 'BCY')) \
            + (u(L.index('Z')) - mform(L.index('Z'), 'AYZ'))
    raise ValueError(name)


def lower(X):
    pts = [i for i in range(7) if X >> i & 1]
    tot = H(X)
    for i in pts:
        tot = tot - (H(1 << i) - lbar(i))
    return tot


f9 = form({})
for X in range(1, 128):
    c = int(LB2[X])
    if c > 0:
        f9 = f9 + c * upper(X)
    elif c < 0:
        f9 = f9 + c * lower(X)
f9 = 3 * f9
assert all(x.denominator == 1 for x in f9)
f9i = np.array([int(x) for x in f9], dtype=np.int64)
z = np.load(REPO + '/data/s31_six_party.npz')
print('f9 recomputed from the transfer equals the stored f9:', np.array_equal(f9i, z['f9'].astype(np.int64)),
      ' terms:', int((f9i != 0).sum()))
ok &= np.array_equal(f9i, z['f9'].astype(np.int64))
print('f9(h_F) =', int(f9i @ hF), ' f9(r_F) =', int(f9i @ rF), ' f9(r_F*) =', int(f9i @ rFs),
      ' f9(h_NF) =', int(f9i @ hNF), ' f9(r_NF) =', int(f9i @ rNF))
ok &= (int(f9i @ hF), int(f9i @ rF), int(f9i @ rFs), int(f9i @ hNF), int(f9i @ rNF)) == (-3, -3, 0, 81, 33)
ok &= int((f9i != 0).sum()) == 67
# print f9 for the appendix
terms = [(c, ''.join(L[i] for i in range(7) if X >> i & 1)) for X, c in enumerate(f9i) if c]
terms.sort(key=lambda t: (len(t[1]), [L.index(ch) for ch in t[1]]))
print('f9 =', ' '.join('%+d h(%s)' % (c, s_) for c, s_ in terms))

# invariant form on Fano-invariant self-dual vectors: parameters v1, v2, v3L, v3N
def inv_vec(v1, v2, v3L, v3N):
    v7 = Fraction(7, 2) * v1
    v4L = v3L - 3 * v1 + v7            # line complements: v3L = 3 v1 + v4L - v7
    v4N = v3N - 3 * v1 + v7
    v5 = v2 - 2 * v1 + v7
    v6 = v7
    out = np.array([Fraction(0)] * 128, dtype=object)
    for X in range(1, 128):
        k = kind(X)
        out[X] = {'1': v1, '2': v2, 'line': v3L, 'triple': v3N, 'linecompl': v4L, '4set': v4N, '5': v5, '6': v6, '7': v7}[k]
    return out
v1, v2, v3L, v3N = sympy.symbols('v1 v2 v3L v3N')
val = sum(int(f9i[X]) * inv_vec(v1, v2, v3L, v3N)[X] for X in range(1, 128))
val = sympy.expand(val / 3)
alpha, mu, nu = 2 * v1 - v2, v1 + v2 - v3L, v1 + v2 - v3N
beta = v1 - 2 * mu
print('f9/3 on the invariant slice:', val)
print('  = 210 v1 + 30 alpha + 21 nu - 421 mu:', sympy.simplify(val - (210 * v1 + 30 * alpha + 21 * nu - 421 * mu)) == 0)
print('  = 210 beta + 30 alpha + 21 nu - mu:', sympy.simplify(val - (210 * beta + 30 * alpha + 21 * nu - mu)) == 0)
ok &= sympy.simplify(val - (210 * beta + 30 * alpha + 21 * nu - mu)) == 0
# check that h_F is on the slice with alpha = beta = nu = 0, mu = 1
print('h_F parameters: v1, v2, v3L, v3N =', hF[1], hF[3], hF[m("ABW")], hF[m("ABC")])

# ------------------------------------------------------------------ faces (Proposition 4)
def tight_rank(h):
    vals_ = E[:, 1:] @ h[1:]
    assert (vals_ >= 0).all()
    T = E[vals_ == 0][:, 1:]
    return len(T), rank_Q(T.tolist())
expected_faces = {'h_F': (385, 125), 'r_F': (490, 126), 'r_F*': (490, 126), 'h_NF': (367, 125), 'r_NF': (487, 126),
                  'r_NF*': (487, 126)}
for name, h in [('h_F', hF), ('r_F', rF), ('r_F*', rFs), ('h_NF', hNF), ('r_NF', rNF), ('r_NF*', rNFs)]:
    n_t, rk = tight_rank(h)
    print('tight Shannon elementals at %-5s: %3d, rank %d' % (name, n_t, rk))
    ok &= (n_t, rk) == expected_faces[name]

# ------------------------------------------------------------------ Appendix H: the two misprints
def rank_vectors(vecs, p):
    """rank function over GF(p) of the 7 subspaces spanned by vecs[i] (lists of vectors), as a 128-vector."""
    r = np.zeros(128, dtype=np.int64)
    for X in range(1, 128):
        rows = [v for i in range(7) if X >> i & 1 for v in vecs[i]]
        r[X] = rank_mod_p(rows, p) if rows else 0
    return r
e1, e2, e3 = [1, 0, 0], [0, 1, 0], [0, 0, 1]
ex_bkl = {'A': [e1], 'B': [e2], 'C': [e3], 'Z': [[1, 1, 1]], 'W': [], 'X': [], 'Y': []}
vals_printed, vals_corr = [], []
for p_ in (2, 3, 5, 7):
    r_ = rank_vectors([ex_bkl[c] for c in L], p_)
    vals_printed.append(int(BKL['bkl617_odd'] @ r_)); vals_corr.append(int(BKL['bkl617_odd_corrected'] @ r_))
def Ivec(x, y, z=0):
    v = np.zeros(128, dtype=np.int64)
    for X_, c in ((x | z, 1), (y | z, 1), (x | y | z, -1), (z, -1)):
        if X_: v[X_] += c
    return v
ident = np.array_equal(BKL['bkl617_odd_corrected'], LB2 - 6 * Ivec(m('A'), m('B')))
print('App. H, BKL (6.17): printed form on A,B,C = e1,e2,e3, Z = e1+e2+e3 over GF(2,3,5,7):', vals_printed,
      ' corrected:', vals_corr, ' corrected on Fano:', int(BKL['bkl617_odd_corrected'] @ rF),
      ' corrected = (LB)_2 - 6 I(A;B):', ident)
ok &= vals_printed == [-3] * 4 and vals_corr == [6] * 4 and int(BKL['bkl617_odd_corrected'] @ rF) == -1 and ident
# DFZ arXiv:1401.2507v1, Thm 3.1, as transcribed in Appendix H (variables A,B,C,D,W,X,Y,Z); RHS - LHS on C = X = <e1>
def dfz31(h):
    Hs = lambda *a: h(frozenset(''.join(a)))
    Hc = lambda x, *cond: Hs(x, *cond) - Hs(*cond)
    rhs = (8 * Hs('Z') + 29 * Hs('Y') + 3 * Hs('X') + 8 * Hs('W') - 6 * Hs('D') - 17 * Hs('C') - 8 * Hs('B') - 17 * Hs('A')
           + 55 * Hc('Z', 'A', 'B', 'C') + 35 * Hc('Y', 'W', 'X', 'Z') + 50 * Hc('X', 'A', 'C', 'D') + 49 * Hc('W', 'B', 'C', 'D')
           + 18 * Hc('A', 'B', 'D', 'Y') + 7 * Hc('B', 'D', 'X', 'Z') + Hc('B', 'A', 'W', 'X') + 7 * Hc('C', 'D', 'Y', 'Z')
           + 7 * Hc('C', 'B', 'X', 'Y') + 3 * Hc('C', 'A', 'W', 'Y') + 6 * Hc('D', 'A', 'W', 'Z')
           + 49 * (Hs('A') + Hs('B') + Hs('C') + Hs('D') - Hs('A', 'B', 'C', 'D')))
    return rhs - Hs('A')
def h_CX(S):            # C = X = the same one-dimensional subspace, all other subspaces zero
    return 1 if ('C' in S or 'X' in S) else 0
v31 = dfz31(lambda S: h_CX(S) if S else 0)
print('App. H, DFZ Thm 3.1 as transcribed: right-hand side minus left-hand side on C = X, others 0:', v31)
ok &= v31 == -4
print('ALL OK' if ok else 'FAILURE')
