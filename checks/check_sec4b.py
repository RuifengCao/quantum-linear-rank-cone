"""Checks for Sec. IV.B-C / Appendix B (five parties), independent of the repository code.

Coordinates: S on the 31 nonempty subsets of the parties A..E (mask bits 0..4); the purifier F is bit 5; Z is bit 6.
"""
import itertools, json, sys
from fractions import Fraction
import numpy as np
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from common import REPO, rank_Q

L = 'ABCDEFZ'
FULL6 = 63


def m(s):
    return sum(1 << L.index(c) for c in s)


def S_full(S, M):
    """S(M) for a mask M over the six systems A..F (purity)."""
    if M == 0 or M == FULL6:
        return 0
    if M >> 5 & 1:
        M = FULL6 ^ M
    return S[M - 1]


def h_of(S):
    """lift h_S on masks over A..F."""
    h = [0] * 64
    for M in range(1, 64):
        h[M] = S_full(S, M) + sum(S_full(S, 1 << i) for i in range(6) if M >> i & 1)
    return h


def pull_form(F):
    """linear form F on h (dict mask over A..F -> coeff) -> coefficient vector on the 31 S coordinates."""
    v = [Fraction(0)] * 31
    for k in range(31):
        e = [0] * 31
        e[k] = 1
        h = h_of(e)
        v[k] = sum(c * h[M] for M, c in F.items())
    return v


def I(a, b, K=''):
    A, B, C = m(a), m(b), m(K)
    r = {}
    for M, c in ((A | C, 1), (B | C, 1), (A | B | C, -1), (C, -1)):
        if M:
            r[M] = r.get(M, 0) + c
    return r


def Hc(a, K):
    A, C = m(a), m(K)
    return {A | C: 1, C: -1}


def combine(terms):
    out = {}
    for coef, row in terms:
        for M, c in row.items():
            out[M] = out.get(M, 0) + coef * c
    return {M: c for M, c in out.items() if c}


def eliminate_Z(F, X, Y):
    """substitute g(XZ) = h(X), g(YZ) = h(Y), g(Z) = h(X) + h(Y) - h(XY); fail if another Z-term survives."""
    ZB = 1 << 6
    out = {}
    for M, c in F.items():
        if not M & ZB:
            out[M] = out.get(M, 0) + c
        elif M == X | ZB:
            out[X] = out.get(X, 0) + c
        elif M == Y | ZB:
            out[Y] = out.get(Y, 0) + c
        elif M == ZB:
            for MM, cc in ((X, c), (Y, c), (X | Y, -c)):
                out[MM] = out.get(MM, 0) + cc
        else:
            raise ValueError('Z-term survives: ' + str(M))
    return {M: c for M, c in out.items() if c}


ok = True
# ------------------------------------------------------------------ QLR5 facets and the witness ray
Hq = np.load(REPO + '/data/qlr5_H_facets10860.npy')
print('QLR5 facet instances:', Hq.shape)
F300 = np.load(REPO + '/data/a2_frontier300.npy')
r = [int(x) for x in F300[275]]
print('witness r =', r)
vals = Hq.astype(np.int64) @ np.array(r, dtype=np.int64)
tight = Hq[vals == 0]
print('r in QLR5:', bool((vals >= 0).all()), ' tight facets:', len(tight), ' exact rank:', rank_Q(tight.tolist()))
ok &= bool((vals >= 0).all()) and rank_Q(tight.tolist()) == 30

R641 = np.load(REPO + '/data/stab5_extreme_reps.npy')
q2 = [int(x) for x in R641[195]]
v2 = Hq.astype(np.int64) @ np.array(q2, dtype=np.int64)
t2 = Hq[v2 == 0]
print('q2 in QLR5:', bool((v2 >= 0).all()), ' tight rank:', rank_Q(t2.tolist()))
ok &= bool((v2 >= 0).all()) and rank_Q(t2.tolist()) == 30

# ------------------------------------------------------------------ the inequality (star)
lhs = ['BCD', 'ABCD', 'ABE', 'BCE', 'BCE', 'ADE', 'ADE', 'CDE']
rhs = ['BC', 'BC', 'ABD', 'ACD', 'AE', 'DE', 'BCDE', 'ABCDE']
star = [0] * 31
for s in lhs:
    star[m(s) - 1] += 1
for s in rhs:
    star[m(s) - 1] -= 1
print('star on r:', sum(a * b for a, b in zip(star, r)))
ok &= sum(a * b for a, b in zip(star, r)) == -1

# the lift of r satisfies the Shannon inequalities and every Ingleton instance (Eq. (3)) on the six systems
hr = h_of(r)
def Ih(a, b, k):
    return hr[a | k] + hr[b | k] - hr[a | b | k] - hr[k]
sh_min = min([hr[FULL6] - hr[FULL6 ^ (1 << i)] for i in range(6)] +
             [Ih(1 << i, 1 << j, K) for i in range(6) for j in range(i + 1, 6)
              for K in range(64) if not K >> i & 1 and not K >> j & 1])
ing_min, n_ing = None, 0
for assign in itertools.product(range(6), repeat=6):   # 0..3: A, B, C, D of Eq. (3); 4: K; 5: unused
    sets = [0] * 5
    for e, s in enumerate(assign):
        if s < 5:
            sets[s] |= 1 << e
    a, b, c3, d3, k = sets
    if not (a and b and c3 and d3):
        continue
    v = Ih(a, b, c3 | k) + Ih(a, b, d3 | k) + Ih(c3, d3, k) - Ih(a, b, k)
    ing_min = v if ing_min is None else min(ing_min, v)
    n_ing += 1
print('lift of r: min over elemental Shannon', sh_min, '; Ingleton instances (ordered)', n_ing, 'min', ing_min)
ok &= sh_min >= 0 and ing_min >= 0

terms = [(5, Hc('Z', 'ABCDEF')), (2, I('F', 'Z', 'ABCDE'))] + [(1, I(*t)) for t in [
    ('A', 'B', 'Z'), ('A', 'C', 'BFZ'), ('A', 'E', 'BCDZ'), ('A', 'Z', 'BCDEF'), ('A', 'Z', 'BE'), ('B', 'D', 'AEZ'),
    ('B', 'D', 'AZ'), ('B', 'Z', 'ACDEF'), ('B', 'Z', 'AF'), ('C', 'Z', 'BDEF'), ('C', 'Z', 'ADEF'),
    ('C', 'Z', 'ABDE'), ('C', 'Z', 'BE'), ('C', 'Z', 'ABD'), ('D', 'E', 'AFZ'), ('D', 'E', 'BCZ'),
    ('D', 'Z', 'ABCEF'), ('D', 'Z', 'BEF'), ('E', 'F', 'BZ'), ('E', 'Z', 'ABCF'), ('E', 'Z', 'AD'),
    ('F', 'Z', 'AE'), ('F', 'Z', 'AD'), ('F', 'Z', 'BC')]]
print('number of Shannon terms in the certificate:', len(terms))
Fh = eliminate_Z(combine(terms), m('AD'), m('BE'))
Fs = pull_form(Fh)
ratio = {Fraction(a) / b for a, b in zip(Fs, star) if b != 0}
same_support = all((a != 0) == (b != 0) for a, b in zip(Fs, star))
print('certificate pulls back to a positive multiple of (star):', same_support and len(ratio) == 1,
      ' ratio:', ratio)
ok &= same_support and ratio == {Fraction(1)}
print('linear form in h (six elements):', {''.join(L[i] for i in range(6) if M >> i & 1): c for M, c in sorted(Fh.items())})

# ------------------------------------------------------------------ facet check of (star) on qubit stabilizer vectors
perms = list(itertools.permutations(range(6)))


def act(S, p):
    """S6 image: (pS)(X) = S(p^{-1} X)."""
    inv = [0] * 6
    for i, j in enumerate(p):
        inv[j] = i
    out = []
    for M in range(1, 32):
        MM = sum(1 << inv[i] for i in range(6) if M >> i & 1)
        out.append(S_full(S, MM))
    return out


gs = np.load(REPO + '/data/graphstate_vecs.npy')
print('six-qubit graph-state vectors:', gs.shape)
pool = set(tuple(int(x) for x in v) for v in gs)
for v in R641:
    vv = [int(x) for x in v]
    for p in perms:
        pool.add(tuple(act(vv, p)))
pool = [list(v) for v in pool]
st = np.array(star)
vals = np.array(pool) @ st
print('pool size:', len(pool), ' min of (star) on pool:', int(vals.min()))
tightvecs = [v for v, x in zip(pool, vals) if x == 0]
print('tight vectors:', len(tightvecs), ' exact rank:', rank_Q(tightvecs))
ok &= int(vals.min()) >= 0 and rank_Q(tightvecs) == 30

# ------------------------------------------------------------------ one common information is not enough: smallest example
d = json.load(open(REPO + '/data/s22_ingleton_exclusions.json'))
cand = [c for c in d['certificates'] if c['catalogue_index'] == 22884]
print('certificates for catalogue orbit 22884:', len(cand))
c = cand[0]
r2 = c['ray']
print('ray:', r2, ' X, Y =', c['X'], c['Y'])


def elemental7():
    full = 127
    rows = []
    for i in range(7):
        rows.append({full: 1, full ^ (1 << i): -1})
    for i in range(7):
        for j in range(i + 1, 7):
            rest = full ^ (1 << i) ^ (1 << j)
            K = rest
            while True:
                dd = {}
                for M, cc in ((K | 1 << i, 1), (K | 1 << j, 1), (K | 1 << i | 1 << j, -1), (K, -1)):
                    if M:
                        dd[M] = dd.get(M, 0) + cc
                rows.append(dd)
                if K == 0:
                    break
                K = (K - 1) & rest
    return rows


E7 = elemental7()


def ingleton(a, b, cc, dd, K):
    A, B, C, D = 1 << a, 1 << b, 1 << cc, 1 << dd
    row = {}
    for s, sg in ((A | B, 1), (A | C, 1), (A | D, 1), (B | C, 1), (B | D, 1),
                  (A, -1), (B, -1), (C | D, -1), (A | B | C, -1), (A | B | D, -1)):
        row[s | K] = row.get(s | K, 0) + sg
    return row


tt = []
ning = 0
for desc, v in c['terms']:
    v = Fraction(v)
    assert v >= 0
    if desc[0] == 'E':
        tt.append((v, E7[desc[1]]))
    else:
        a, b, cc, dd, K = desc[1:]
        assert len({a, b, cc, dd}) == 4 and not K & ((1 << a) | (1 << b) | (1 << cc) | (1 << dd))
        tt.append((v, ingleton(a, b, cc, dd, K)))
        ning += 1
Fh2 = eliminate_Z(combine(tt), c['X'], c['Y'])
Fs2 = pull_form(Fh2)
claimed_lhs = ['A', 'AB', 'ABC', 'ABC', 'ABC', 'AD', 'CD', 'CD', 'BCD', 'BCD', 'ACE', 'ACE', 'BCE', 'BCE', 'ADE', 'BDE', 'CDE']
claimed_rhs = ['B', 'AC', 'BC', 'BC', 'BC', 'D', 'ABD', 'ABD', 'ACD', 'ACD', 'ACD', 'AE', 'ABCE', 'ABCE', 'DE', 'BCDE']
cl = [0] * 31
for s in claimed_lhs:
    cl[m(s) - 1] += 1
for s in claimed_rhs:
    cl[m(s) - 1] -= 1
ratio2 = {Fraction(a) / b for a, b in zip(Fs2, cl) if b != 0}
same2 = all((a != 0) == (b != 0) for a, b in zip(Fs2, cl))
print('terms:', len(c['terms']), ' Ingleton terms:', ning, ' pulls back to the printed inequality:', same2 and len(ratio2) == 1,
      ratio2, ' value on the ray:', sum(a * b for a, b in zip(cl, r2)))
ok &= same2 and len(ratio2) == 1 and sum(a * b for a, b in zip(cl, r2)) == -1
v3 = Hq.astype(np.int64) @ np.array(r2, dtype=np.int64)
print('ray in QLR5:', bool((v3 >= 0).all()), ' tight rank:', rank_Q(Hq[v3 == 0].tolist()))
ok &= bool((v3 >= 0).all()) and rank_Q(Hq[v3 == 0].tolist()) == 30

# single-CI feasibility for this ray at every pair, exact
W = np.load(REPO + '/data/s22_single_ci_witnesses.npz')
print('witness file keys:', W.files)
EI = np.zeros((len(E7), 127), dtype=np.int64)
for k, row in enumerate(E7):
    for M, cc in row.items():
        EI[k, M - 1] = cc
h = h_of(r2)
ridx = [i for i in range(len(W['rays'])) if [int(x) for x in W['rays'][i]] == r2]
assert len(ridx) == 1
sel = [k for k in range(len(W['den'])) if int(W['ray'][k]) == ridx[0]]
store = {}
for k in sel:
    store[(int(W['pair'][k][0]), int(W['pair'][k][1]))] = (int(W['den'][k]), [int(t) for t in W['num'][k]])
print('stored witnesses for this ray:', len(store))
npairs = 0
bad = 0
for X in range(1, 64):
    for Y in range(X + 1, 64):
        Iv = h[X] + h[Y] - h[X | Y]
        if (X & Y) == X or (X & Y) == Y:
            Sm = X if (X & Y) == X else Y
            den, gz = 1, [h[M | Sm] if (M | Sm) else 0 for M in range(64)]
        elif Iv == 0:
            den, gz = 1, [h[M] if M else 0 for M in range(64)]
        else:
            if (X, Y) not in store:
                bad += 1
                continue
            den, gz = store[(X, Y)]
        g = [den * h[M] for M in range(1, 64)] + list(gz)
        okp = den > 0 and gz[X] == den * h[X] and gz[Y] == den * h[Y] and gz[0] == den * Iv and \
            bool((EI @ np.array(g, dtype=np.int64) >= 0).all())
        bad += 0 if okp else 1
        npairs += 1
print('single-CI feasibility: pairs checked', npairs, ' failures', bad)
ok &= bad == 0

print('ALL OK' if ok else 'FAILURE')
