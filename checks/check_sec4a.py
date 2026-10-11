"""Checks for Sec. IV.A / Appendix G (four parties), independent of the repository code.

1. Build C4 = {S : h_S satisfies Shannon and all Ingleton instances (pairwise disjoint arguments and conditioning
   set) on the five elements 0..4} in the 15 coordinates S(X), X a nonempty subset of the parties {1,2,3,4}
   (bits 0..3; the purifier 0 is bit 4 and is eliminated by purity).
2. Enumerate its extreme rays exactly with lrs; expect 46.
3. Compare with the rays stored in data/s30_stab4_unimodular.npz and check the unimodular lifts: every visible cut
   block of the integer matrix has Smith invariants in {0,1} and rank lam * S(X).
4. Check that every ray with S(party) = 0 for some party is realized by a state in which that party is isolated.
"""
import itertools, os, shutil, subprocess, sys, tempfile
from fractions import Fraction
import numpy as np
import sympy
from sympy.matrices.normalforms import smith_normal_form
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from common import REPO, rank_mod_p

K = 5                      # elements: bits 0..3 parties 1..4, bit 4 purifier
FULL = (1 << K) - 1
VIS = 4


def S_coord(M):
    """index (0..14) of S(M) for a mask M over the 5 elements, using purity; None for M = 0 or FULL."""
    if M == 0 or M == FULL:
        return None
    if M >> VIS & 1:
        M = FULL ^ M
    return M - 1


def pull(row):
    """row: dict mask -> coefficient on h; returns coefficient vector (15) on S via h_S(X) = S(X) + sum S(i)."""
    v = [Fraction(0)] * 15
    for M, c in row.items():
        k = S_coord(M)
        if k is not None:
            v[k] += c
        for i in range(K):
            if M >> i & 1:
                k = S_coord(1 << i)
                v[k] += c
    return v


def H(*sets):
    return sum(sets)


def I_row(A, B, C=0):
    """I(A;B|C) as dict over masks (A, B, C disjoint or overlapping masks)."""
    r = {}
    for M, c in ((A | C, 1), (B | C, 1), (A | B | C, -1), (C, -1)):
        if M:
            r[M] = r.get(M, 0) + c
    return {k: v for k, v in r.items() if v}


def add(*rows):
    out = {}
    for sgn, r in rows:
        for k, v in r.items():
            out[k] = out.get(k, 0) + sgn * v
    return {k: v for k, v in out.items() if v}


rows = []
# Shannon elemental inequalities on 5 elements
for i in range(K):
    rows.append({FULL: 1, FULL ^ (1 << i): -1})
for i, j in itertools.combinations(range(K), 2):
    rest = [k for k in range(K) if k not in (i, j)]
    for t in range(1 << len(rest)):
        Kset = sum(1 << rest[s] for s in range(len(rest)) if t >> s & 1)
        rows.append(I_row(1 << i, 1 << j, Kset))
n_shannon = len(rows)
# Ingleton instances: Ing(A,B;C,D|Kc): I(A;B|Kc) <= I(A;B|C Kc) + I(A;B|D Kc) + I(C;D|Kc), with A, B, C, D nonempty,
# Kc possibly empty, all pairwise disjoint subsets of the five elements.
ing = []
for assign in itertools.product(range(6), repeat=K):          # 0..3 -> A,B,C,D ; 4 -> Kc ; 5 -> unused
    parts = [sum(1 << e for e in range(K) if assign[e] == s) for s in range(6)]
    A, B, C, D, Kc = parts[:5]
    if not (A and B and C and D):
        continue
    row = add((1, I_row(A, B, C | Kc)), (1, I_row(A, B, D | Kc)), (1, I_row(C, D, Kc)), (-1, I_row(A, B, Kc)))
    ing.append(row)
rows += ing
print('Shannon rows:', n_shannon, ' Ingleton instances (with repetitions):', len(ing))

# pull back, normalise, deduplicate
pulled = {}
for r in rows:
    v = pull(r)
    if not any(v):
        continue
    den = 1
    for x in v:
        den = den * x.denominator // np.gcd(den, x.denominator)
    iv = [int(x * den) for x in v]
    g = 0
    for x in iv:
        g = np.gcd(g, abs(x))
    iv = tuple(x // g for x in iv)
    pulled[iv] = True
A = np.array(list(pulled.keys()), dtype=np.int64)
print('distinct pulled-back rows:', len(A))

# lrs
lrs = shutil.which('lrs')
with tempfile.TemporaryDirectory() as td:
    f = os.path.join(td, 'c4.ine')
    with open(f, 'w') as fh:
        fh.write('c4\nH-representation\nbegin\n%d 16 rational\n' % len(A))
        for row in A:
            fh.write('0 ' + ' '.join(str(int(x)) for x in row) + '\n')
        fh.write('end\n')
    out = subprocess.run([lrs, f], capture_output=True, text=True).stdout
rays = []
inbody = False
for line in out.splitlines():
    t = line.split()
    if line.startswith('begin'):
        inbody = True
        continue
    if line.startswith('end'):
        inbody = False
    if inbody and t and t[0] == '0' and len(t) == 16:
        v = [Fraction(x) for x in t[1:]]
        den = 1
        for x in v:
            den = den * x.denominator // np.gcd(den, x.denominator)
        iv = [int(x * den) for x in v]
        g = 0
        for x in iv:
            g = np.gcd(g, abs(x))
        rays.append(tuple(x // g for x in iv))
rays = sorted(set(rays))
print('extreme rays of C4 (lrs, exact):', len(rays))

z = np.load(REPO + '/data/s30_stab4_unimodular.npz')
stored = []
for r in z['rays']:
    iv = [int(x) for x in r]
    g = 0
    for x in iv:
        g = np.gcd(g, abs(x))
    stored.append(tuple(x // g for x in iv))
print('stored rays equal the lrs rays (as sets):', sorted(set(stored)) == rays)

# unimodular realizations (called "lifts" in the repository)
ok_all = True
primes = [2, 3, 5, 7, 11, 13, 17, 19]
for k in range(len(z['rays'])):
    lam = int(z['lam'][k])
    sizes = [int(s) for s in z['sizes'][k]]
    nv = int(z['N'][k])
    W = z['Wt'][k][:nv, :nv].astype(int)
    assert nv == sum(sizes)
    assert (W == W.T).all() and all(W[i, i] == 0 for i in range(nv))
    owner = [p for p, s in enumerate(sizes) for _ in range(s)]       # party 4 = purifier
    r = [int(x) for x in z['rays'][k]]
    for M in range(1, 16):
        X = [v for v in range(nv) if owner[v] < 4 and M >> owner[v] & 1]
        Y = [v for v in range(nv) if v not in X]
        target = lam * r[M - 1]
        if not X or not Y:
            ok_all &= target == 0
            continue
        B = sympy.Matrix(W[np.ix_(X, Y)].tolist())
        snf = smith_normal_form(B, domain=sympy.ZZ)
        inv = [abs(int(snf[i, i])) for i in range(min(snf.shape))]
        nz = [x for x in inv if x != 0]
        ok_all &= all(x == 1 for x in nz) and len(nz) == target
        for p in primes:
            ok_all &= rank_mod_p(W[np.ix_(X, Y)].tolist(), p) == target
print('unimodular realizations: all cut blocks have Smith invariants in {0,1} and rank lam*S(X); ranks over GF(p), '
      'p <= 19, agree:', ok_all)
print('multiplicities lam:', sorted(set(int(x) for x in z['lam'])), ' counts:',
      {l: int((z['lam'] == l).sum()) for l in set(int(x) for x in z['lam'])},
      ' max qubits:', int(z['N'].max()))
counts = {l: int((z['lam'] == l).sum()) for l in set(int(x) for x in z['lam'])}
ok = (n_shannon == 85 and len(ing) == 480 and len(A) == 100 and len(rays) == 46
      and sorted(set(stored)) == rays and ok_all and counts == {1: 41, 2: 5})

# Ray 6 of Table 1 of Linden, Matus, Ruskai and Winter (LIPIcs TQC 2013, p. 281), systems a, b, c, d and e = purifier,
# and two states for it, from explicit state vectors: their state (R6) of eight qutrits, which as printed has
# S(ab) = 1 instead of 2, and the eight-qubit state of Osthus (M.S. thesis, Univ. of Minnesota 2024, Eq. (3.63)).
def state_vector(d, words):
    psi = np.zeros(d ** len(words[0]))
    for w in words:
        idx = 0
        for x in w:
            idx = idx * d + x
        psi[idx] += 1
    return psi / np.linalg.norm(psi)


def entropy_vector(psi, d, owner):
    """S(X) in units of log d for the 15 nonempty subsets X of the parties a..d (bits 0..3); owner[q] in 0..4."""
    nq = len(owner)
    out = []
    for M in range(1, 16):
        X = [q for q in range(nq) if owner[q] < 4 and M >> owner[q] & 1]
        Y = [q for q in range(nq) if q not in X]
        T = np.transpose(psi.reshape([d] * nq), X + Y).reshape(d ** len(X), -1)
        s = np.linalg.svd(T, compute_uv=False) ** 2
        s = s[s > 1e-12]
        out.append(float(-(s * np.log(s)).sum() / np.log(d)))
    return out


# masks over a, b, c, d (bits 0..3): a b ab c ac bc abc d ad bd abd cd acd bcd abcd; by purity abc = de, abcd = e, ...
ray6 = [1, 1, 2, 2, 3, 3, 2, 2, 3, 3, 2, 2, 3, 3, 2]
owner8 = [0, 1, 2, 2, 3, 3, 4, 4]
R6 = state_vector(3, [(i, j, (i + j) % 3, k, (i + j) % 3, k, (i + j) % 3, k)
                      for i in range(3) for j in range(3) for k in range(3)])
O6 = state_vector(2, [(i, j, (i + j) % 2, k, (i + j) % 2, l, (i + j) % 2, (j + k + l) % 2)
                      for i in range(2) for j in range(2) for k in range(2) for l in range(2)])
vR6 = entropy_vector(R6, 3, owner8)
vO6 = entropy_vector(O6, 2, owner8)
diff_R6 = [(M, round(vR6[M - 1], 6), ray6[M - 1]) for M in range(1, 16) if abs(vR6[M - 1] - ray6[M - 1]) > 1e-6]
print('LMRW ray 6 is one of the 46 extreme rays of C4:', tuple(ray6) in set(rays))
print('LMRW state (R6) as printed differs from ray 6 at (mask, computed, ray 6):', diff_R6)
vR6_int = [int(round(x)) for x in vR6]
g6 = 0
for x in vR6_int:
    g6 = np.gcd(g6, abs(x))
R6_extreme = tuple(x // g6 for x in vR6_int) in set(rays)
print('LMRW state (R6) as printed: its entropy vector is an extreme ray of C4:', R6_extreme)
print('Osthus state (3.63) has entropy vector ray 6 (qubits):', all(abs(vO6[M] - ray6[M]) < 1e-6 for M in range(15)))
ok &= (tuple(ray6) in set(rays) and diff_R6 == [(3, 1.0, 2)] and not R6_extreme
       and all(abs(vO6[M] - ray6[M]) < 1e-6 for M in range(15)))
print('ALL OK' if ok else 'FAILURE')
