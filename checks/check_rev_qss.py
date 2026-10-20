#!/usr/bin/env python3
"""check_rev_qss.py -- independent check of the secret-sharing results (Lemma "Stabilizer schemes are linear schemes",
Proposition "Secret sharing", Supplementary Note 10) and of the distances of the witnesses from the cones (Methods).

Independent of the repository code: nothing from epr1kit/ or scripts/ is imported, and the paper source is not read.
The script reads only the data files data/s38_qss_certificates.json (rational dual certificates of the two linear
programs with Lemma 1 (LB), t = 2), data/dfz_char_dependent.csv (transcription of the inequalities (65) and (91) of
Dougherty, Freiling and Zeger, arXiv:1311.4601v1; used for the margins) and data/s31_six_party.npz (key 'f9', the form
f9 on the 128 subsets of A,B,C,W,X,Y,Z in mask order; used for the margins).  The repository is found through the
environment variable EPR1_REPO, or else as the parent directory of the directory of this script (the repository itself
when the script sits in <repo>/checks/).  Requirements: Python 3.8+, numpy.  Running time: well under a minute.

Labels: players A,B,C,W,X,Y,Z = vertices 100,010,001,110,101,011,111 of the cube = bits 0..6; the reference R (the
origin) is bit 7; for mixed schemes the environment E is bit 7 and R is bit 8.  Linear-program variables: h(X) for
nonempty X at index X-1 (mask order), then t at index 2^N - 1; objective c = t.  The forms LB2 and LA2 are Lemma 1 (LB)
and (LA) with t = 2, M = 3, in the variable order A_1, A_2, A_3, B_1, B_2, B_3, C = bits 0..6; in a 'form' descriptor
the form variable b is carried by the system perm[b].  Exact computations use integers and fractions.Fraction; floating
point is used only for the Euclidean norms of (e) and the state-vector cross-check of (a).

What is checked:
 (a) the cube codes (affine functions on {0,1}^3) over GF(2), GF(3), GF(5), GF(7): ranks, the CSS entropy vector
     S = lambda_r, its lift h_S = 2r, self-duality, the access structures Gamma_F (GF(2)) and Gamma_NF (odd p): every set
     authorized or unauthorized, minimal authorized sets, 64/64, complement property, difference only in WXY/ABCZ, the
     same structure for p = 3, 5, 7, ratio 1; the Steane identification; the matroids: AG(3,2) (circuits = the 14 affine
     planes) and R8 = AG(3,2) with exactly the two complementary circuit-hyperplanes {R,W,X,Y} and {A,B,C,Z} relaxed
     (rank over GF(p) = rank over Q), both ports at R; contractions at R are the Fano and non-Fano matroids.
 (b) contractions of the access structures: Gamma_F/q is isomorphic to the Fano port F = 14+25+36+123+156+246+345 of
     Jafari-Khazaei for every q, with minimal triples = the four lines missing q; Gamma_NF/q is isomorphic to
     N = F + 456 for q in {A,B,C,Z}, and every isomorphism maps 456 to the other three of A,B,C,Z; for q in {W,X,Y} it
     is isomorphic to neither; F and N are ports of the Fano and non-Fano matroids; the conditioning lemma (the scheme
     conditioned on a share has entropy function h(Tq) - h(q)) on random linear codes (a test).
 (c) LB2 and LA2 built from Lemma 1; their expansion in the labeling (A_1,A_2,A_3,B_1,B_2,B_3,C) = (A,B,C,Y,X,W,Z) against
     the coefficients printed in the paper (transcribed here); values on the configuration J - I over GF(2) and GF(3);
     minima over all 16 seven-element minors and 5040 labelings of the cube states; tests of LB2 on random and perturbed
     arrangements in groups of odd order (GF(3), GF(5), Z_9, Z_27) and of LA2 on random GF(2) arrangements.
 (d) the two certificates: every row rebuilt from its descriptor and checked to be admissible (signs; access rows
     against the access structure of (a); purity rows as consequences of h* = h; forms LB2 on seven-element minors that
     remove E when E is present), the identity c - sum y a - sum z e >= 0 exactly, bound 43/41, row counts, negative
     controls; the sizes of the linear programs; the controls (LA2 or Shannon only: optimum 1).
 (e) the margins v/||g||_1 and v/||g||_2 for f9 (Steane state), f12 (built from its definition; non-Fano state) and
     DFZ (65)/(91) pulled back through the purifier O of the seven-party cube states.
The validity of Lemma 1, of the results of Jafari-Khazaei and of DFZ is an input; it is only tested here.
"""
import itertools, json, math, os, random, sys
from fractions import Fraction

import numpy as np

REPO = os.environ.get('EPR1_REPO', os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA = os.path.join(REPO, 'data')

ok_all = True


def report(name, good, detail=''):
    global ok_all
    ok_all = ok_all and bool(good)
    print(('OK   ' if good else 'FAIL ') + name + (('  [' + detail + ']') if detail else ''), flush=True)


L7 = 'ABCWXYZ'
VERT = {'A': (1, 0, 0), 'B': (0, 1, 0), 'C': (0, 0, 1), 'W': (1, 1, 0), 'X': (1, 0, 1), 'Y': (0, 1, 1),
        'Z': (1, 1, 1)}
LINES7 = ['ABW', 'ACX', 'AYZ', 'BCY', 'BXZ', 'CWZ', 'WXY']          # lines of the Fano plane


def mask(s, letters=L7):
    m = 0
    for ch in s:
        m |= 1 << letters.index(ch)
    return m


def name(m, letters=L7 + 'R'):
    return ''.join(letters[i] for i in range(len(letters)) if m >> i & 1) or '{}'


def popcount(m):
    return bin(m).count('1')


def rank_mod_p(rows, p):
    """rank over GF(p) of a list of integer vectors (Gaussian elimination)."""
    M = [[x % p for x in r] for r in rows]
    M = [r for r in M if any(r)]
    rk, col = 0, 0
    ncol = len(M[0]) if M else 0
    while rk < len(M) and col < ncol:
        piv = next((i for i in range(rk, len(M)) if M[i][col]), None)
        if piv is None:
            col += 1
            continue
        M[rk], M[piv] = M[piv], M[rk]
        inv = pow(M[rk][col], p - 2, p)
        M[rk] = [(x * inv) % p for x in M[rk]]
        for i in range(len(M)):
            if i != rk and M[i][col]:
                f = M[i][col]
                M[i] = [(a - f * b) % p for a, b in zip(M[i], M[rk])]
        rk += 1
        col += 1
    return rk


def rank_Q(rows):
    """rank over Q (exact)."""
    M = [[Fraction(x) for x in r] for r in rows if any(r)]
    rk, col = 0, 0
    ncol = len(M[0]) if M else 0
    while rk < len(M) and col < ncol:
        piv = next((i for i in range(rk, len(M)) if M[i][col] != 0), None)
        if piv is None:
            col += 1
            continue
        M[rk], M[piv] = M[piv], M[rk]
        for i in range(len(M)):
            if i != rk and M[i][col] != 0:
                f = M[i][col] / M[rk][col]
                M[i] = [a - f * b for a, b in zip(M[i], M[rk])]
        rk += 1
        col += 1
    return rk


def rank_function(cols, p):
    """rank function (list over all masks, mask 0 included) of a list of column vectors over GF(p) (p = 0: over Q)."""
    k = len(cols)
    sel = lambda m: [cols[i] for i in range(k) if m >> i & 1]
    if p == 0:
        return [rank_Q(sel(m)) if m else 0 for m in range(1 << k)]
    return [rank_mod_p(sel(m), p) if m else 0 for m in range(1 << k)]


def zn_log_order(vecs, k, N, base):
    """log_base of the order of the subgroup of Z_N^k generated by vecs (N a power of base), via the index of the
    lattice spanned by vecs and N Z^k in Z^k (integer row echelon form)."""
    rows = [list(v) for v in vecs] + [[N if i == j else 0 for j in range(k)] for i in range(k)]
    det, r = 1, 0
    for c in range(k):
        while True:
            nz = [i for i in range(r, len(rows)) if rows[i][c] != 0]
            piv = min(nz, key=lambda i: abs(rows[i][c]))
            rows[r], rows[piv] = rows[piv], rows[r]
            again = False
            for i in range(r + 1, len(rows)):
                if rows[i][c]:
                    q = rows[i][c] // rows[r][c]
                    rows[i] = [a - q * b for a, b in zip(rows[i], rows[r])]
                    again = again or rows[i][c] != 0
            if not again:
                break
        det *= abs(rows[r][c])
        r += 1
    order, e = N ** k // det, 0
    assert N ** k % det == 0
    while order > 1:
        assert order % base == 0
        order //= base
        e += 1
    return e


def dual(f, n):
    """f*(X) = sum_{i in X} f(i) + f(N minus X) - f(N), f a list over all masks of n elements."""
    full = (1 << n) - 1
    return [sum(f[1 << i] for i in range(n) if X >> i & 1) + f[full ^ X] - f[full] for X in range(1 << n)]


def conn(f, n):
    """connectivity function lambda_f(X) = f(X) + f(N minus X) - f(N)."""
    full = (1 << n) - 1
    return [f[X] + f[full ^ X] - f[full] for X in range(1 << n)]


def lift(S, n):
    """h_S(X) = S(X) + sum_{i in X} S(i)."""
    return [S[X] + sum(S[1 << i] for i in range(n) if X >> i & 1) for X in range(1 << n)]


def is_polymatroid(h, n):
    """all elemental inequalities h(N) >= h(N - i) and I(i;j|K) >= 0, and h(empty) = 0."""
    full = (1 << n) - 1
    if h[0] != 0 or any(h[full] < h[full ^ (1 << i)] for i in range(n)):
        return False
    for i, j in itertools.combinations(range(n), 2):
        rest = [k for k in range(n) if k not in (i, j)]
        for mm in range(1 << len(rest)):
            K = sum(1 << rest[b] for b in range(len(rest)) if mm >> b & 1)
            if h[K | 1 << i] + h[K | 1 << j] - h[K | 1 << i | 1 << j] - h[K] < 0:
                return False
    return True


# ===================================================================================== (a) cube codes and matroids
print('(a) cube codes, access structures, matroids')
CUBE_COLS = [(1,) + VERT[c] for c in L7] + [(1, 0, 0, 0)]       # players bits 0..6, R = origin bit 7
NSYS, FULL8, RB = 8, 255, 1 << 7
cube = {}
for q in (2, 3, 5, 7):
    r = rank_function(CUBE_COLS, q)
    S = conn(r, NSYS)                      # CSS entropy vector lambda_r (fact (b)), units of log q
    h = lift(S, NSYS)
    rstar = [popcount(X) + r[FULL8 ^ X] - r[FULL8] for X in range(256)]
    good = (r[FULL8] == 4 and all(r[1 << i] == 1 for i in range(NSYS)) and rstar == r and h == [2 * x for x in r]
            and dual(h, NSYS) == h and all(S[X] == S[FULL8 ^ X] for X in range(256)) and is_polymatroid(h, NSYS))
    report('GF(%d) cube code: rank 4, r* = r (identically self-dual), S = lambda_r pure-symmetric, h_S = 2r, h* = h, '
           'h a polymatroid' % q, good)
    cube[q] = dict(r=r, S=S, h=h)
rQ = rank_function(CUBE_COLS, 0)
report('rank over GF(p) equals rank over Q on every set of columns for p = 3, 5, 7 (the real affine cube); over GF(2) '
       'it differs', all(cube[q]['r'] == rQ for q in (3, 5, 7)) and cube[2]['r'] != rQ)


def state_vector_entropies(q):
    """entropies (units of log q) of the CSS state of the cube code over GF(q), from the explicit state vector."""
    words = set()
    for a in itertools.product(range(q), repeat=4):
        words.add(tuple(sum(a[k] * col[k] for k in range(4)) % q for col in CUBE_COLS))
    psi = np.zeros(q ** NSYS)
    for w in words:
        psi[sum(w[i] * q ** (NSYS - 1 - i) for i in range(NSYS))] = 1.0
    psi /= np.linalg.norm(psi)
    T = psi.reshape([q] * NSYS)
    out = [0.0] * 256
    for X in range(1, 256):
        A = [i for i in range(NSYS) if X >> i & 1]
        B = [i for i in range(NSYS) if not X >> i & 1]
        M = np.transpose(T, A + B).reshape(q ** len(A), q ** len(B))
        sv = np.linalg.svd(M, compute_uv=False) ** 2
        sv = sv[sv > 1e-12]
        out[X] = float(-(sv * np.log(sv)).sum() / math.log(q))
    return len(words), out


for q in (2, 3):
    nw, Snum = state_vector_entropies(q)
    err = max(abs(Snum[X] - cube[q]['S'][X]) for X in range(256))
    report('GF(%d) cube state: entropies of the explicit state vector (%d codewords) agree with lambda_r' % (q, nw),
           nw == q ** 4 and err < 1e-9, 'max error %.1e' % err)


def access(h, Rbit=7):
    """classify the 128 sets of players: 'auth' if h(RA) = h(A), 'unauth' if h(RA) = h(A) + h(R), else None."""
    Rm = 1 << Rbit
    return {A: ('auth' if h[A | Rm] == h[A] else 'unauth' if h[A | Rm] == h[A] + h[Rm] else None) for A in range(128)}


def minimal(auth_sets):
    return sorted((A for A in auth_sets if not any(B != A and B & A == B for B in auth_sets)),
                  key=lambda m: (popcount(m), m))


ACC = {q: access(cube[q]['h']) for q in (2, 3, 5, 7)}
AUTH = {q: {A for A, k in ACC[q].items() if k == 'auth'} for q in ACC}
MIN_F = sorted(mask(s) for s in LINES7)
MIN_NF = sorted(mask(s) for s in ['ABW', 'ACX', 'AYZ', 'BCY', 'BXZ', 'CWZ', 'ABCZ', 'AWXY', 'BWXY', 'CWXY', 'WXYZ'])
for q in (2, 3, 5, 7):
    acc = ACC[q]
    n_auth = sum(1 for k in acc.values() if k == 'auth')
    n_un = sum(1 for k in acc.values() if k == 'unauth')
    compl = all((acc[A] == 'auth') == (acc[127 ^ A] == 'unauth') for A in range(128))
    mins = sorted(minimal(AUTH[q]))
    report('GF(%d): all 128 player sets authorized or unauthorized; %d/%d (empty set unauthorized); complement property; '
           'minimal authorized sets = %s' % (q, n_auth, n_un, 'the 7 Fano lines' if q == 2 else
                                             '6 lines + ABCZ, AWXY, BWXY, CWXY, WXYZ'),
           all(k is not None for k in acc.values()) and n_auth == 64 and n_un == 64 and acc[0] == 'unauth' and compl
           and mins == (MIN_F if q == 2 else MIN_NF), 'minimal: ' + ' '.join(name(m) for m in minimal(AUTH[q])))
report('Gamma_F = sets of players containing a Fano line',
       AUTH[2] == {A for A in range(128) if any(A & mask(L) == mask(L) for L in LINES7)})
report('Gamma_F and Gamma_NF differ exactly in WXY (F only) and ABCZ (NF only); Gamma_NF the same for p = 3, 5, 7',
       AUTH[2] - AUTH[3] == {mask('WXY')} and AUTH[3] - AUTH[2] == {mask('ABCZ')} and AUTH[3] == AUTH[5] == AUTH[7])
report('ratio 1: S(P) = S(R) = 1 for every player, over GF(2), GF(3), GF(5), GF(7)',
       all(cube[q]['S'][1 << i] == 1 for q in cube for i in range(8)))

dets = set()
for k in range(1, 5):
    for rows in itertools.combinations(range(4), k):
        for cols in itertools.combinations(range(8), k):
            Mx = [[Fraction(CUBE_COLS[c][rr]) for c in cols] for rr in rows]
            d = Fraction(1)
            for i in range(k):
                piv = next((j for j in range(i, k) if Mx[j][i] != 0), None)
                if piv is None:
                    d = Fraction(0)
                    break
                if piv != i:
                    Mx[i], Mx[piv] = Mx[piv], Mx[i]
                    d = -d
                d *= Mx[i][i]
                for j in range(i + 1, k):
                    f = Mx[j][i] / Mx[i][i]
                    Mx[j] = [a - f * b for a, b in zip(Mx[j], Mx[i])]
            dets.add(int(d))
report('all minors of the 4x8 integer generator matrix lie in {0, +-1, +-2}', dets <= {0, 1, -1, 2, -2},
       'values %s' % sorted(dets))

M3 = [[VERT[c][k] for c in L7] for k in range(3)]
simplex = {tuple(sum(a[k] * M3[k][i] for k in range(3)) % 2 for i in range(7)) for a in itertools.product((0, 1), repeat=3)}
code2 = {tuple(sum(a[k] * col[k] for k in range(4)) % 2 for col in CUBE_COLS) for a in itertools.product((0, 1), repeat=4)}
hamming = {v for v in itertools.product((0, 1), repeat=7)
           if all(sum(v[i] * M3[k][i] for i in range(7)) % 2 == 0 for k in range(3))}
even = {v for v in hamming if sum(v) % 2 == 0}
odd = {v for v in hamming if sum(v) % 2 == 1}
ones = lambda S_: {tuple((x + 1) % 2 for x in v) for v in S_}
report('GF(2): extended Hamming code [8,4,4]; codewords with c_R = 0 give the simplex code on the players, the others its '
       'translate by 1111111; these are the even/odd-weight words of the Hamming code (Steane |0_L>, |1_L>)',
       len(code2) == 16 and min(sum(w) for w in code2 if any(w)) == 4 and {w[:7] for w in code2 if w[7] == 0} == simplex
       and {w[:7] for w in code2 if w[7] == 1} == ones(simplex) and even == simplex and odd == ones(simplex))
psi = np.zeros(256)
for w in code2:
    psi[sum(w[i] << i for i in range(8))] = 0.25
tgt = np.zeros(256)
for v in even:
    tgt[sum(v[i] << i for i in range(7))] += 0.25
for v in odd:
    tgt[sum(v[i] << i for i in range(7)) | 128] += 0.25
report('qubit cube state = (|0>_R|0_L> + |1>_R|1_L>)/sqrt(2) for Steane\'s code (numerical)', np.allclose(psi, tgt))


def circuits(r, n):
    deps = {X for X in range(1, 1 << n) if r[X] < popcount(X)}
    return [X for X in sorted(deps) if all((X ^ (1 << i)) not in deps for i in range(n) if X >> i & 1)]


def is_flat(r, X, n):
    return all(r[X | 1 << i] > r[X] for i in range(n) if not X >> i & 1)


rA, rR = cube[2]['r'], cube[3]['r']
circA = circuits(rA, 8)
planes = {X for X in range(256) if popcount(X) == 4 and
          all(sum(CUBE_COLS[i][k] for i in range(8) if X >> i & 1) % 2 == 0 for k in range(1, 4))}
portA = {A for A in range(128) if any(Cc >> 7 & 1 and Cc & ~(A | RB) == 0 for Cc in circA)}
report('GF(2): the matroid is AG(3,2): its circuits are the 14 affine planes (all of size 4); Gamma_F is its port at R',
       set(circA) == planes and len(planes) == 14 and portA == AUTH[2])
diff = [X for X in range(256) if rR[X] != rA[X]]
CH = [mask('WXY') | RB, mask('ABCZ')]
relax_ok = (sorted(diff) == sorted(CH) and all(rR[X] == rA[X] + 1 for X in diff)
            and all(X in circA and is_flat(rA, X, 8) and rA[X] == 3 for X in CH) and CH[0] & CH[1] == 0
            and CH[0] | CH[1] == FULL8)
circR = circuits(rR, 8)
portR = {A for A in range(128) if any(Cc >> 7 & 1 and Cc & ~(A | RB) == 0 for Cc in circR)}
report('GF(3) (= Q): the matroid R8 is AG(3,2) with exactly the two complementary circuit-hyperplanes {R,W,X,Y} and '
       '{A,B,C,Z} relaxed (ranks differ only there, by 1); Gamma_NF is its port at R',
       relax_ok and portR == AUTH[3], 'sets where the ranks differ: %s' % ' '.join(name(X) for X in diff))
rF = rank_function([VERT[c] for c in L7], 2)
rNF = rank_function([VERT[c] for c in L7], 3)
contrA = [rA[U | RB] - rA[RB] for U in range(128)]
contrR = [rR[U | RB] - rR[RB] for U in range(128)]
report('contraction at R: AG(3,2)/R = Fano matroid, R8/R = non-Fano matroid (on A..Z)', contrA == rF and contrR == rNF)

# ============================================================================ (b) contractions and Jafari-Khazaei
print('(b) contractions of the access structures and the structures F, N of Jafari-Khazaei')
JK_F = ['14', '25', '36', '123', '156', '246', '345']
JK_N = JK_F + ['456']


def closure6(minsets):
    ms = [sum(1 << (int(c) - 1) for c in s) for s in minsets]
    return {A for A in range(64) if any(A & m == m for m in ms)}


F6, N6 = closure6(JK_F), closure6(JK_N)


def contraction(auth, q):
    """Gamma/q on the six players other than q, as masks over the list `others`."""
    others = [i for i in range(7) if i != q]
    out = set()
    for A6 in range(64):
        A = sum(1 << others[j] for j in range(6) if A6 >> j & 1)
        if A | (1 << q) in auth:
            out.add(A6)
    return others, out


def isomorphisms(G6, target):
    """bijections sigma (positions 0..5 of `others` -> JK players 0..5) with A in G6 iff sigma(A) in target."""
    out = []
    for sig in itertools.permutations(range(6)):
        img = lambda A: sum(1 << sig[j] for j in range(6) if A >> j & 1)
        if all((A in G6) == (img(A) in target) for A in range(64)):
            out.append(sig)
    return out


okF, lines_ok = True, True
for q in range(7):
    others, G = contraction(AUTH[2], q)
    isoF, isoN = isomorphisms(G, F6), isomorphisms(G, N6)
    okF = okF and len(isoF) > 0 and not isoN
    triples = {sum(1 << others[j] for j in range(6) if A >> j & 1) for A in minimal(G) if popcount(A) == 3}
    lines_ok = lines_ok and triples == {mask(L) for L in LINES7 if not mask(L) >> q & 1}
report('for every player q, Gamma_F/q is isomorphic to F (and not to N); its minimal authorized triples are the four '
       'lines missing q', okF and lines_ok)
corr, okN = {}, True
for q in range(7):
    others, G = contraction(AUTH[3], q)
    isoF, isoN = isomorphisms(G, F6), isomorphisms(G, N6)
    if L7[q] in 'ABCZ':
        imgs = {sum(1 << others[j] for j in range(6) if sig[j] in (3, 4, 5)) for sig in isoN}
        corr[L7[q]] = ' '.join(sorted(name(m) for m in imgs))
        okN = okN and bool(isoN) and not isoF and imgs == {mask('ABCZ') & ~(1 << q)}
    else:
        corr[L7[q]] = 'neither' if not isoF and not isoN else 'ISOMORPHIC'
        okN = okN and not isoF and not isoN
report('Gamma_NF/q is isomorphic to N for q in {A,B,C,Z}, and every isomorphism maps 456 to the other three of A,B,C,Z; '
       'for q in {W,X,Y} it is isomorphic to neither F nor N', okN, ', '.join('%s: %s' % kv for kv in corr.items()))
# F and N are ports of the Fano and non-Fano matroids (port at a point p of the 7-point configuration)
port_ok = True
for p_, rk in ((2, rF), (3, rNF)):
    for q in range(7):
        others = [i for i in range(7) if i != q]
        G = {A6 for A6 in range(64)
             if rk[sum(1 << others[j] for j in range(6) if A6 >> j & 1) | 1 << q]
             == rk[sum(1 << others[j] for j in range(6) if A6 >> j & 1)]}
        isoF, isoN = isomorphisms(G, F6), isomorphisms(G, N6)
        if p_ == 2:
            port_ok = port_ok and bool(isoF)
        else:
            port_ok = port_ok and (bool(isoN) if L7[q] in 'ABCZ' else not isoN and not isoF)
report('F is the port of the Fano matroid at every point; N is the port of the non-Fano matroid at a corner '
       '(A, B, C or Z; at W, X, Y the port is neither)', port_ok)
in_min = all(any(A >> P & 1 for A in minimal(AUTH[q])) for q in (2, 3) for P in range(7))
report('every player lies in a minimal authorized set of Gamma_F and of Gamma_NF; single players are unauthorized',
       in_min and all(ACC[q][1 << P] == 'unauth' for q in (2, 3) for P in range(7)))


def nullspace_left(G, p):
    """basis of {x : x G = 0} over GF(p), G a list of rows (k x n)."""
    k = len(G)
    n = len(G[0]) if G else 0
    M = [[G[i][j] % p for j in range(n)] + [1 if i == l else 0 for l in range(k)] for i in range(k)]
    rk, col = 0, 0
    while rk < k and col < n:
        piv = next((i for i in range(rk, k) if M[i][col]), None)
        if piv is None:
            col += 1
            continue
        M[rk], M[piv] = M[piv], M[rk]
        inv = pow(M[rk][col], p - 2, p)
        M[rk] = [(x * inv) % p for x in M[rk]]
        for i in range(k):
            if i != rk and M[i][col]:
                f = M[i][col]
                M[i] = [(a - f * b) % p for a, b in zip(M[i], M[rk])]
        rk += 1
        col += 1
    return [row[n:] for row in M[rk:]]


rng = random.Random(20261020)
cond_ok = True
for trial in range(24):
    p_ = (2, 3)[trial % 2]
    owner = [rng.randrange(8) for _ in range(10)]               # 10 coordinates distributed among 8 systems
    kdim = rng.randint(3, 7)
    G = [[rng.randrange(p_) for _ in range(10)] for _ in range(kdim)]
    def hfun(Gm):
        return [rank_mod_p([[row[j] for j in range(10) if X >> owner[j] & 1] for row in Gm], p_) if X else 0
                for X in range(256)]
    hh = hfun(G)
    q = rng.randrange(8)
    Gq = [[row[j] for j in range(10) if owner[j] == q] for row in G]
    K = nullspace_left(Gq, p_) if any(owner[j] == q for j in range(10)) else [[1 if i == l else 0 for l in range(kdim)]
                                                                             for i in range(kdim)]
    G0 = [[sum(x[i] * G[i][j] for i in range(kdim)) % p_ for j in range(10)] for x in K]
    h0 = hfun(G0) if G0 else [0] * 256
    cond_ok = cond_ok and all(h0[T] == hh[T | 1 << q] - hh[1 << q] for T in range(256) if not T >> q & 1)
report('test: conditioning a random linear code on a vanishing share has entropy function h(Tq) - h(q) (24 random '
       'codes over GF(2), GF(3))', cond_ok)

# ===================================================================================== (c) the forms LB2 and LA2
print('(c) Lemma 1 with t = 2: the forms LB2 and LA2')


def lemma1_forms(M=3):
    """(LB) and (LA) of Lemma 1 for t = M - 1, variables A_1..A_M = bits 0..M-1, B_1..B_M = bits M..2M-1, C = bit 2M,
    as {mask: coefficient} with g . h >= 0 when valid.  Error terms as in the lemma:
      alpha = sum_k h(A_k) - h(A), gamma = h(C|A), delta_i = I(C; A_{!=i}), eps_i = h(C|A_i B_i),
      zeta_i = h(B_i|A_{!=i}), eta_i = h(B_i|A_i C);
      (LB): h(B_1..B_M) + M[gamma + M alpha + sum_i (delta_i + eps_i + zeta_i)] - M h(C) >= 0,
      (LA): (M-1) h(C) + M alpha + M gamma + sum_i (zeta_i + eta_i) - h(B_1..B_M) >= 0."""
    Ab = [1 << i for i in range(M)]
    Bb = [1 << (M + i) for i in range(M)]
    Cb = 1 << (2 * M)
    Aall = sum(Ab)
    Ball = sum(Bb)

    def add(f, X, c):
        if X:
            f[X] = f.get(X, 0) + c

    def alpha(f, c):
        for a in Ab:
            add(f, a, c)
        add(f, Aall, -c)

    def gamma(f, c):
        add(f, Cb | Aall, c); add(f, Aall, -c)

    def delta(f, i, c):
        Ao = Aall ^ Ab[i]
        add(f, Cb, c); add(f, Ao, c); add(f, Cb | Ao, -c)

    def eps(f, i, c):
        add(f, Cb | Ab[i] | Bb[i], c); add(f, Ab[i] | Bb[i], -c)

    def zeta(f, i, c):
        Ao = Aall ^ Ab[i]
        add(f, Bb[i] | Ao, c); add(f, Ao, -c)

    def eta(f, i, c):
        add(f, Bb[i] | Ab[i] | Cb, c); add(f, Ab[i] | Cb, -c)
    LB, LA = {}, {}
    add(LB, Ball, 1); gamma(LB, M); alpha(LB, M * M); add(LB, Cb, -M)
    for i in range(M):
        delta(LB, i, M); eps(LB, i, M); zeta(LB, i, M)
    add(LA, Cb, M - 1); alpha(LA, M); gamma(LA, M); add(LA, Ball, -1)
    for i in range(M):
        zeta(LA, i, 1); eta(LA, i, 1)
    return {m: c for m, c in LB.items() if c}, {m: c for m, c in LA.items() if c}


LB2, LA2 = lemma1_forms(3)
FORMS = {'LB2': LB2, 'LA2': LA2}
# expansion in the labeling (A_1,A_2,A_3,B_1,B_2,B_3,C) = (A,B,C,Y,X,W,Z), against the coefficients printed in the
# paper (Table of the eight inequalities, columns (LB) and (LA); Eq. (LB2)), transcribed here
TO_AZ = [mask(c) for c in 'ABCYXWZ']


def relabel(form):
    out = {}
    for m, c in form.items():
        out[sum(TO_AZ[b] for b in range(7) if m >> b & 1)] = c
    return out


PRINTED_LB = {'A': 9, 'B': 9, 'C': 9, 'Z': 6, 'AY': -3, 'BX': -3, 'CW': -3, 'ABC': -12, 'ABW': 3, 'ABZ': -3, 'ACX': 3,
              'ACZ': -3, 'AYZ': 3, 'BCY': 3, 'BCZ': -3, 'BXZ': 3, 'CWZ': 3, 'WXY': 1, 'ABCZ': 3}
PRINTED_LA = {'A': 3, 'B': 3, 'C': 3, 'Z': 2, 'AB': -1, 'AC': -1, 'AZ': -1, 'BC': -1, 'BZ': -1, 'CZ': -1, 'ABC': -6,
              'ABW': 1, 'ACX': 1, 'AYZ': 1, 'BCY': 1, 'BXZ': 1, 'CWZ': 1, 'WXY': -1, 'ABCZ': 3}
report('LB2 and LA2 from Lemma 1 (%d and %d terms) equal, in the labeling (A,B,C,Y,X,W,Z), the printed columns (LB) '
       'and (LA) term by term' % (len(LB2), len(LA2)),
       relabel(LB2) == {mask(k): v for k, v in PRINTED_LB.items()} and relabel(LA2) == {mask(k): v for k, v in
                                                                                       PRINTED_LA.items()})


def ev(form, h):
    return sum(c * h[m] for m, c in form.items())


CONF = [(1, 0, 0), (0, 1, 0), (0, 0, 1), (0, 1, 1), (1, 0, 1), (1, 1, 0), (1, 1, 1)]   # e_i, c - e_i, c
conf = {p_: rank_function(CONF, p_) for p_ in (2, 3, 5)}
report('on the configuration J - I (t = 2): LB2 = -1 over GF(2), LB2 >= 0 over GF(3), GF(5); LA2 = -1 over GF(3), '
       'GF(5), LA2 >= 0 over GF(2)',
       ev(LB2, conf[2]) == -1 and ev(LB2, conf[3]) >= 0 and ev(LB2, conf[5]) >= 0 and ev(LA2, conf[3]) == -1
       and ev(LA2, conf[5]) == -1 and ev(LA2, conf[2]) >= 0,
       'LB2: %d, %d, %d; LA2: %d, %d, %d' % (ev(LB2, conf[2]), ev(LB2, conf[3]), ev(LB2, conf[5]), ev(LA2, conf[2]),
                                             ev(LA2, conf[3]), ev(LA2, conf[5])))

PERMS = np.array(list(itertools.permutations(range(7))), dtype=np.int64)
PMAP = np.zeros((len(PERMS), 128), dtype=np.int64)          # mask over form variables -> mask over positions
for m in range(128):
    for b in range(7):
        if m >> b & 1:
            PMAP[:, m] |= 1 << PERMS[:, b]


def form_vector(form):
    v = np.zeros(128, dtype=np.int64)
    for m, c in form.items():
        v[m] = c
    return v


FV = {k: form_vector(v) for k, v in FORMS.items()}


def minors_min(h, fv, n=8):
    """minimum of each form over all seven-element minors (one element contracted or deleted) and all labelings."""
    out = {}
    for nm, g in fv.items():
        best = None
        for x in range(n):
            rest = [y for y in range(n) if y != x]
            for mode in 'cd':
                C = (1 << x) if mode == 'c' else 0
                hm = np.array([h[sum(1 << rest[j] for j in range(7) if mm >> j & 1) | C] - h[C] for mm in range(128)],
                              dtype=np.int64)
                v = int((hm[PMAP] @ g).min())
                best = v if best is None else min(best, v)
        out[nm] = best
    return out


MINV = {q: minors_min(cube[q]['h'], FV) for q in (2, 3, 5, 7)}
report('GF(2) cube state: LB2 has minimum -2 over all 16 minors x 5040 labelings (the programs exclude it); LA2 >= 0 '
       'on all of them', MINV[2]['LB2'] == -2 and MINV[2]['LA2'] >= 0, 'minima %s' % MINV[2])
report('GF(3), GF(5), GF(7) cube states: LB2 >= 0 on all minors and labelings; LA2 has minimum -2',
       all(MINV[q]['LB2'] >= 0 and MINV[q]['LA2'] == -2 for q in (3, 5, 7)),
       '; '.join('GF(%d) %s' % (q, MINV[q]) for q in (3, 5, 7)))


def random_rank_function(field_or_ring, rng, n=7, dim=3):
    """rank function of 7 random subgroups (1-2 generators each, sometimes 0) of R^dim, R = GF(p) or Z_N; returned
    in units of log p (N a power of p)."""
    kind, N, base = field_or_ring
    gens = [[tuple(rng.randrange(N) for _ in range(dim)) for _ in range(rng.choice((0, 1, 1, 2, 2)))] for _ in range(n)]
    out = []
    for m in range(1 << n):
        vecs = [v for i in range(n) if m >> i & 1 for v in gens[i]]
        if not vecs:
            out.append(0)
        elif kind == 'field':
            out.append(rank_mod_p(vecs, N))
        else:
            out.append(zn_log_order(vecs, dim, N, base))
    return out


def perturbed_configuration(field_or_ring, rng, dim=3):
    """the configuration J - I with each subgroup enlarged, with probability 1/2, by a random extra generator."""
    kind, N, base = field_or_ring
    gens = [[v] + ([tuple(rng.randrange(N) for _ in range(dim))] if rng.random() < 0.5 else []) for v in CONF]
    out = []
    for m in range(128):
        vecs = [v for i in range(7) if m >> i & 1 for v in gens[i]]
        out.append(0 if not vecs else rank_mod_p(vecs, N) if kind == 'field' else zn_log_order(vecs, dim, N, base))
    return out


tests = {}
cases_t = [(('field', 3, 3), 'LB2'), (('field', 5, 5), 'LB2'), (('ring', 9, 3), 'LB2'), (('ring', 27, 3), 'LB2'),
           (('field', 2, 2), 'LA2')]
for ring, nm in cases_t:
    mn = None
    for trial in range(30):
        rk = random_rank_function(ring, rng) if trial % 2 == 0 else perturbed_configuration(ring, rng)
        v = int((np.array(rk, dtype=np.int64)[PMAP] @ FV[nm]).min())
        mn = v if mn is None else min(mn, v)
    tests['%s/%s%d' % (nm, 'GF' if ring[0] == 'field' else 'Z', ring[1])] = mn
zconf = {N: [zn_log_order([v for i in range(7) if m >> i & 1 for v in [CONF[i]]], 3, N, b) if m else 0
             for m in range(128)] for N, b in ((9, 3), (4, 2))}
report('test: LB2 >= 0 in all 5040 labelings on 30 random and perturbed arrangements each over GF(3), GF(5), Z_9, Z_27 '
       '(groups of odd order), LA2 >= 0 on 30 over GF(2); on J - I over Z_9 LB2 >= 0, over Z_4 LB2 < 0',
       all(v >= 0 for v in tests.values()) and ev(LB2, zconf[9]) >= 0 and ev(LB2, zconf[4]) < 0,
       'minima %s; J - I: Z_9 %d, Z_4 %d' % (tests, ev(LB2, zconf[9]), ev(LB2, zconf[4])))

# ========================================================================================= (d) the certificates
print('(d) the certificates of the two linear programs')
CASES = {   # witness field, quantum (purity rows allowed), environment E, form, expected bound, expected row counts
    'F_pure_LB2': (2, True, False, 'LB2', '43/41', (72, 109)),
    'F_mixed_LB2': (2, True, True, 'LB2', '43/41', (140, 201)),
}


class Inadmissible(Exception):
    pass


def is_int(x):
    return isinstance(x, int) and not isinstance(x, bool)


class Case:
    def __init__(self, label):
        self.label = label
        self.q, self.quantum, self.env, self.form, self.expected, self.rows = CASES[label]
        self.N = 9 if self.env else 8
        self.FULL = (1 << self.N) - 1
        self.Rbit = self.N - 1
        self.Rm = 1 << self.Rbit
        self.Ebit = 7 if self.env else None
        self.T = self.FULL                         # index of t
        self.acc = ACC[self.q]                     # access structure of the witness (players 0..6)
        h8 = cube[self.q]['h']                     # witness: cube state over GF(q), R = bit 7; E trivial
        self.hw = [h8[(X & 127) | ((X >> 8 & 1) << 7)] for X in range(1 << 9)] if self.env else list(h8)

    def e(self, X, c, d):
        if X:
            d[X - 1] = d.get(X - 1, 0) + c

    def S_row(self, Y, d, sign):
        """sign * S(Y) in the lift coordinates: S(Y) = h(Y) - 1/2 sum_{i in Y} h(i)."""
        self.e(Y, Fraction(sign), d)
        for i in range(self.N):
            if Y >> i & 1:
                self.e(1 << i, Fraction(-sign, 2), d)

    def selfdual_functional(self, Y):
        """E_Y = h*(Y) - h(Y), with h*(Y) = sum_{i in Y} h(i) + h(N minus Y) - h(N)."""
        d = {}
        for i in range(self.N):
            if Y >> i & 1:
                self.e(1 << i, Fraction(1), d)
        self.e(self.FULL ^ Y, Fraction(1), d)
        self.e(self.FULL, Fraction(-1), d)
        self.e(Y, Fraction(-1), d)
        return {k: v for k, v in d.items() if v}

    def ineq(self, desc):
        """row a of an admissible inequality a.x <= 0."""
        N, FULL = self.N, self.FULL
        kind = desc[0]
        d = {}
        if kind == 'mono':
            if len(desc) != 2 or not is_int(desc[1]) or not 0 <= desc[1] < N:
                raise Inadmissible('mono')
            self.e(FULL, -1, d)
            self.e(FULL ^ (1 << desc[1]), 1, d)
        elif kind == 'sub':
            if len(desc) != 4 or not all(is_int(x) for x in desc[1:]):
                raise Inadmissible('sub')
            i, j, K = desc[1:]
            if not (0 <= i < N and 0 <= j < N and i != j and 0 <= K <= FULL and not K & (1 << i | 1 << j)):
                raise Inadmissible('sub')
            for X, c in ((K | 1 << i, -1), (K | 1 << j, -1), (K | 1 << i | 1 << j, 1), (K, 1)):
                self.e(X, c, d)
        elif kind == 'share':
            if len(desc) != 2 or not is_int(desc[1]) or not 0 <= desc[1] < 7:
                raise Inadmissible('share')
            self.e(1 << desc[1], 1, d)
            d[self.T] = -2
        elif kind == 'form':
            if len(desc) != 5:
                raise Inadmissible('form arity')
            nm, removed, modes, perm = desc[1], list(desc[2]), list(desc[3]), list(desc[4])
            if nm != self.form:
                raise Inadmissible('form %s not the form of this case' % nm)
            if (len(removed) != N - 7 or len(set(removed)) != len(removed)
                    or not all(is_int(x) and 0 <= x < N for x in removed)):
                raise Inadmissible('form: removed elements')
            if self.env and self.Ebit not in removed:
                raise Inadmissible('form: minor does not remove E')
            if len(modes) != len(removed) or any(m not in ('c', 'd') for m in modes):
                raise Inadmissible('form: modes')
            rest = sorted(set(range(N)) - set(removed))
            if len(perm) != 7 or not all(is_int(x) for x in perm) or sorted(perm) != rest:
                raise Inadmissible('form: labeling')
            C = sum(1 << x for x, md in zip(removed, modes) if md == 'c')
            for m, c in FORMS[nm].items():            # g(minor) >= 0, minor(T) = h(T u C) - h(C); row = -g(minor)
                Tm = sum(1 << perm[b] for b in range(7) if m >> b & 1)
                self.e(Tm | C, -c, d)
                self.e(C, c, d)
        else:
            raise Inadmissible('unknown inequality kind %r' % (kind,))
        return {k: v for k, v in d.items() if v}

    def eq(self, desc):
        """(row e, right-hand side b) of an admissible equation e.x = b."""
        kind = desc[0]
        d = {}
        if kind == 'sd':
            if not self.quantum:
                raise Inadmissible('purity row in a classical case')
            if len(desc) != 2 or not is_int(desc[1]) or not 0 <= desc[1] <= self.FULL:
                raise Inadmissible('sd')
            X = desc[1]
            self.S_row(self.FULL ^ X, d, 1)        # S(N - X) - S(X) = 0 (purity), written in the lift
            self.S_row(X, d, -1)
            d = {k: v for k, v in d.items() if v}
            EX, EN = self.selfdual_functional(X), self.selfdual_functional(self.FULL)
            comb = dict(EX)
            for k, v in EN.items():
                comb[k] = comb.get(k, 0) - v / 2
            if {k: v for k, v in comb.items() if v} != d:
                raise Inadmissible('sd row is not E_X - E_N/2')
            return d, Fraction(0)
        if kind in ('auth', 'unauth'):
            if len(desc) != 2 or not is_int(desc[1]) or not 1 <= desc[1] <= 127:
                raise Inadmissible(kind)
            A = desc[1]
            if self.acc[A] != kind:
                raise Inadmissible('%s %s does not match the access structure' % (kind, name(A)))
            self.e(A | self.Rm, 1, d)
            self.e(A, -1, d)
            if kind == 'unauth':
                self.e(self.Rm, -1, d)
            return {k: v for k, v in d.items() if v}, Fraction(0)
        if kind == 'norm':
            if len(desc) != 1:
                raise Inadmissible('norm')
            self.e(self.Rm, 1, d)
            return d, Fraction(2)
        raise Inadmissible('unknown equation kind %r' % (kind,))

    def verify(self, cert):
        """(ok, bound, message[, s]): rebuild the rows, check admissibility and signs, check s >= 0."""
        s = {self.T: Fraction(1)}
        try:
            for desc, y in cert['ineq']:
                y = Fraction(y)
                if y > 0:
                    return False, None, 'positive multiplier on %s' % (desc,)
                for k, v in self.ineq(desc).items():
                    s[k] = s.get(k, 0) - y * v
            bound = Fraction(0)
            for desc, z in cert['eq']:
                z = Fraction(z)
                e, b = self.eq(desc)
                for k, v in e.items():
                    s[k] = s.get(k, 0) - z * v
                bound += z * b
        except Inadmissible as ex:
            return False, None, 'inadmissible row: %s' % ex
        neg = [k for k, v in s.items() if v < 0]
        if neg:
            return False, bound, '%d negative entries of c - sum y a - sum z e' % len(neg)
        return True, bound, 'ok', s

    def xw(self, k):
        return Fraction(1) if k == self.T else Fraction(self.hw[k + 1])


def run_verify(case, cert):
    r = case.verify(cert)
    return r[0], r[1], r[2]


cert_all = json.load(open(os.path.join(DATA, 's38_qss_certificates.json')))
report('certificate file has exactly the two cases F_pure_LB2 and F_mixed_LB2', sorted(cert_all) == sorted(CASES),
       ', '.join(sorted(cert_all)))
for label in CASES:
    case = Case(label)
    cert = cert_all[label]
    res = case.verify(cert)
    ok, bound, msg = res[0], res[1], res[2]
    nI, nE = len(cert['ineq']), len(cert['eq'])
    distinct = (len({json.dumps(d) for d, _ in cert['ineq']}) == nI and len({json.dumps(d) for d, _ in cert['eq']}) == nE)
    nonzero = all(Fraction(y) != 0 for _, y in cert['ineq'] + cert['eq'])
    kinds = {}
    for d_, _ in cert['ineq'] + cert['eq']:
        kinds[d_[0]] = kinds.get(d_[0], 0) + 1
    report('%s: every row admissible, y <= 0, identity exact, certified bound %s = expected %s = stored %s; rows %d/%d '
           '(expected %d/%d), distinct and nonzero' % (label, bound, case.expected, cert['bound'], nI, nE, *case.rows),
           ok and bound == Fraction(case.expected) == Fraction(cert['bound']) and (nI, nE) == case.rows and distinct
           and nonzero, msg + '; ' + ', '.join('%s %d' % kv for kv in sorted(kinds.items())))
    s = res[3] if ok else {}
    xw = case.xw
    eq_ok = all(sum(v * xw(k) for k, v in case.eq(d_)[0].items()) == case.eq(d_)[1] for d_, _ in cert['eq'])
    sh_ok = all(sum(v * xw(k) for k, v in case.ineq(d_).items()) <= 0 for d_, _ in cert['ineq'] if d_[0] != 'form')
    nviol = sum(1 for d_, _ in cert['ineq'] if d_[0] == 'form' and sum(v * xw(k) for k, v in case.ineq(d_).items()) > 0)
    nform = sum(1 for d_, _ in cert['ineq'] if d_[0] == 'form')
    report('%s: s(h(R)) = s(t) = 0; the witness (GF(2) cube state, t = 1) satisfies all equations and Shannon/share rows '
           'used and violates %d of the %d form rows used' % (label, nviol, nform),
           ok and s.get(case.Rm - 1, 0) == 0 and s.get(case.T, 0) == 0 and eq_ok and sh_ok and nviol > 0)
    controls = []
    c1 = dict(ineq=cert['ineq'], eq=[[d_, str(Fraction(z) + Fraction(1, 10 ** 9)) if d_[0] == 'norm' else z]
                                     for d_, z in cert['eq']])
    controls.append(('norm multiplier + 1e-9', not run_verify(case, c1)[0]))
    controls.append(('form rows dropped', not run_verify(case, dict(ineq=[r_ for r_ in cert['ineq'] if r_[0][0] != 'form'],
                                                                     eq=cert['eq']))[0]))
    k_acc = next(i for i, (d_, _) in enumerate(cert['eq']) if d_[0] in ('auth', 'unauth'))
    c3 = dict(ineq=cert['ineq'], eq=[[[('unauth' if d_[0] == 'auth' else 'auth'), d_[1]] if i == k_acc else d_, z]
                                     for i, (d_, z) in enumerate(cert['eq'])])
    controls.append(('one access row flipped', not run_verify(case, c3)[0]))
    c4 = dict(ineq=[[[d_[0], 'LA2'] + list(d_[2:]), y] if d_[0] == 'form' else [d_, y] for d_, y in cert['ineq']],
              eq=cert['eq'])
    controls.append(('LA2 in place of LB2', not run_verify(case, c4)[0]))
    k_s = next(i for i, (d_, _) in enumerate(cert['ineq']) if d_[0] == 'sub')
    c5 = dict(ineq=[[d_, str(-Fraction(y))] if i == k_s else [d_, y] for i, (d_, y) in enumerate(cert['ineq'])],
              eq=cert['eq'])
    controls.append(('sign of one multiplier flipped', not run_verify(case, c5)[0]))
    k_f = next(i for i, (d_, _) in enumerate(cert['ineq']) if d_[0] == 'form')
    if case.env:
        d_ = cert['ineq'][k_f][0]
        removed = list(d_[2])
        other_el = next(x for x in range(9) if x not in removed)
        rem2 = sorted([x for x in removed if x != case.Ebit] + [other_el])
        rest = sorted(set(range(9)) - set(rem2))
        c6 = dict(ineq=[[[d_[0], d_[1], rem2, list(d_[3]), rest], y] if i == k_f else [dd, y]
                        for i, (dd, y) in enumerate(cert['ineq'])], eq=cert['eq'])
        controls.append(('a minor that keeps E', not run_verify(case, c6)[0]))
    k_big = max((i for i, (d_, _) in enumerate(cert['ineq']) if d_[0] == 'form'),
                key=lambda i: abs(Fraction(cert['ineq'][i][1])))
    fails = 0
    for sgn in (1, -1):
        c8 = dict(ineq=[[d_, str(Fraction(y) * (1 + sgn * Fraction(1, 10 ** 6)))] if i == k_big else [d_, y]
                        for i, (d_, y) in enumerate(cert['ineq'])], eq=cert['eq'])
        fails += not run_verify(case, c8)[0]
    controls.append(('largest form multiplier perturbed by +-1e-6 (%d of 2 directions fail)' % fails, fails >= 1))
    report('%s: negative controls rejected (%s)' % (label, '; '.join(c for c, _ in controls)), all(g for _, g in controls),
           ', '.join('%s: %s' % (c, 'rejected' if g else 'NOT REJECTED') for c, g in controls if not g))
    s_ac = sum(abs(Fraction(z)) for d_, z in cert['eq'] if d_[0] in ('auth', 'unauth'))
    print('     %s: sum of |multipliers| on the access equations = %.3f' % (label, float(s_ac)))

for env in (False, True):
    N = 9 if env else 8
    shannon = N + math.comb(N, 2) * 2 ** (N - 2)
    removable = [(x,) for x in range(8)] if not env else [tuple(sorted((7, x))) for x in range(9) if x != 7]
    minors = sum(2 ** len(rm) for rm in removable)
    nforms = minors * math.factorial(7)
    pairs = [X for X in range(1 << N) if X <= ((1 << N) - 1) ^ X]
    cs = Case('F_mixed_LB2' if env else 'F_pure_LB2')
    rowsd = [cs.eq(('sd', X))[0] for X in pairs]
    cnt = {}
    for r_ in rowsd:
        for k in r_:
            cnt[k] = cnt.get(k, 0) + 1
    indep = all(any(cnt[k] == 1 for k in r_) for r_ in rowsd)
    exp_i, exp_e, exp_f = (165904, 384, 161280) if env else (82447, 256, 80640)
    report('LP size %s E: %d Shannon + 7 share + %d form rows (%d minors x 5040) = %d inequalities; %d independent purity '
           '+ 127 access + 1 = %d equations' % ('with' if env else 'without', shannon, nforms, minors,
                                                 shannon + 7 + nforms, len(pairs), len(pairs) + 128),
           shannon + 7 + nforms == exp_i and len(pairs) + 128 == exp_e and nforms == exp_f and indep)
report('controls: the program with LA2 in place of LB2 and the program with Shannon inequalities only have optimum '
       'exactly 1 (the GF(2) cube state is feasible with t = 1, LA2 >= 0 on all its minors; t >= 1 on every feasible '
       'point since every player lies in a minimal authorized set)', MINV[2]['LA2'] >= 0 and in_min
       and all(cube[2]['h'][1 << P] == 2 for P in range(8)))

# ======================================================================================================= (e) margins
print('(e) distances of the witnesses from the cones')


def read_dfz(path):
    forms = {}
    for line in open(path):
        line = line.strip()
        if not line or line.startswith('#') or line.startswith('ineq_id'):
            continue
        nm, subset, coeff = [x.strip() for x in line.split(',')]
        forms.setdefault(nm, {})
        forms[nm][mask(subset)] = forms[nm].get(mask(subset), 0) + int(coeff)
    return {k: {m: c for m, c in v.items() if c} for k, v in forms.items()}


DFZ = read_dfz(os.path.join(DATA, 'dfz_char_dependent.csv'))
f9arr = np.load(os.path.join(DATA, 's31_six_party.npz'))['f9']
F9 = {m: int(f9arr[m]) for m in range(128) if f9arr[m]}
hF = [a + b for a, b in zip(rF, dual(rF, 7))]
SF7, SNF7 = conn(rF, 7), conn(rNF, 7)
report('margin inputs: DFZ (65) = -1 on the Fano and (91) = -1 on the non-Fano configuration; f9 has 67 terms and '
       'f9(h_F) = -3 with h_F = r_F + r_F^* computed here; the lift of the Steane state S_F = lambda_{r_F} is h_F',
       ev(DFZ['dfz65_odd'], rF) == -1 and ev(DFZ['dfz91_even'], rNF) == -1 and len(F9) == 67 and ev(F9, hF) == -3
       and lift(SF7, 7) == hF)


def h_form_to_S7(f):
    c = [0] * 128
    for m, v in f.items():
        c[m] += v
        for i in range(7):
            if m >> i & 1:
                c[1 << i] += v
    return c


def fold(c7, pur):
    g = {}
    for X in range(1, 127):
        T = X if not X >> pur & 1 else 127 ^ X
        g[T] = g.get(T, 0) + c7[X]
    return g


def norms(g):
    vals = [v for v in g.values() if v]
    return len(vals), sum(abs(v) for v in vals), math.sqrt(sum(v * v for v in vals)), sum(v * v for v in vals)


margins = {}
c7 = h_form_to_S7(F9)
f9_all = []
for pur in range(7):
    g = fold(c7, pur)
    val = sum(v * SF7[T] for T, v in g.items())
    f9_all.append((-val,) + norms(g)[1:])
    if pur == 6:
        margins['f9, Steane state, purifier Z'] = (-val,) + norms(g)


def f12(S):
    """the form f12 of Eq. (f12) on a pure-symmetric S given on all 128 subsets of the seven systems."""
    k = lambda T: sum(S[1 << L7.index(Q)] for Q in T) - S[mask(T)]
    kP = lambda P, l: k(l) - k(l.replace(P, ''))
    lines6 = ['ABW', 'ACX', 'AYZ', 'BCY', 'BXZ', 'CWZ']

    def eP(P, l):
        others = [m for m in lines6 if P in m and m != l]
        return 2 * S[mask(P)] - kP(P, others[0]) - kP(P, others[1])
    UZ = 4 * S[mask('Z')] - kP('Z', 'AYZ') - kP('Z', 'BXZ') - kP('Z', 'CWZ')
    UW = 2 * S[mask('W')] - kP('W', 'ABW') + eP('A', 'ABW') + eP('B', 'ABW')
    UX = 2 * S[mask('X')] - kP('X', 'ACX') + eP('A', 'ACX') + eP('C', 'ACX')
    UY = 2 * S[mask('Y')] - kP('Y', 'BCY') + eP('B', 'BCY') + eP('C', 'BCY')
    cp = [a + b for a, b in itertools.combinations('ABCZ', 2)]
    return (5 * UZ + 4 * UW + UX + UY + k('WXY') + 6 * k('ABC') + sum(k(x) for x in cp)
            - sum(k(l) for l in lines6) - 3 * k('ABW') - 3 * k('CWZ'))


g12 = {}
for T in range(1, 64):
    S = [0] * 128
    S[T] = S[127 ^ T] = 1
    v = f12(S)
    if v:
        g12[T] = v
val12 = sum(v * SNF7[T] for T, v in g12.items())
report('f12 from its definition: %d terms; f12(S_NF) = %d, f12(S_F) = %d'
       % (len(g12), val12, sum(v * SF7[T] for T, v in g12.items())),
       len(g12) == 31 and val12 == -1 and sum(v * SF7[T] for T, v in g12.items()) == 0)
margins['f12, non-Fano state, purifier Z'] = (-val12,) + norms(g12)
for nm, q in (('dfz65_odd', 2), ('dfz91_even', 3)):
    g = {}
    for T in range(1, 128):
        S = [0] * 256
        S[T] = S[255 ^ T] = 1
        hS = lift(S, 8)
        v = ev(DFZ[nm], [hS[U | 128] - hS[128] for U in range(128)])
        if v:
            g[T] = v
    val = sum(v * cube[q]['S'][T] for T, v in g.items())
    margins['DFZ (%s) pulled back, GF(%d) cube state, purifier O' % (nm[3:5], q)] = (-val,) + norms(g)
print('     %-50s %3s %6s %9s %6s %10s %9s' % ('form, witness', 'v', '||g||1', '||g||2', 'terms', 'v/||g||1', 'v/||g||2'))
for k_, (v, nt, l1, l2, sq) in margins.items():
    print('     %-50s %3d %6d %9.3f %6d %10s %9.5f' % (k_, v, l1, l2, nt, Fraction(v, l1), v / l2))
m9, m12 = margins['f9, Steane state, purifier Z'], margins['f12, non-Fano state, purifier Z']
m65 = margins['DFZ (65) pulled back, GF(2) cube state, purifier O']
m91 = margins['DFZ (91) pulled back, GF(3) cube state, purifier O']
report('f9: v = 3, ||g||_1 = 2856, ||g||_2 = sqrt(348372) = 590.2 for every choice of purifier; 1/952 = 1.05e-3 and '
       'Euclidean %.5f >= 5.0e-3' % (3 / m9[3]),
       m9[0] == 3 and m9[2] == 2856 and m9[4] == 348372 and all(x == (3, 2856, m9[3]) for x in [a[:3] for a in f9_all])
       and Fraction(3, 2856) == Fraction(1, 952) and 5.0e-3 <= 3 / m9[3] < 5.1e-3)
report('f12: v = 1, ||g||_1 = 129, ||g||_2 = sqrt(839) = 28.97; 1/129 and Euclidean %.4f >= 0.034' % (1 / m12[3]),
       m12[0] == 1 and m12[2] == 129 and m12[4] == 839 and 0.034 <= 1 / m12[3] < 0.035)
report('DFZ (65), (91) pulled back through O: v = 2, ||g||_1 = 90 and 105 (1/45, 2/105), ||g||_2 = %.2f and %.2f, '
       'Euclidean %.4f and %.4f' % (m65[3], m91[3], 2 / m65[3], 2 / m91[3]),
       m65[0] == 2 and m91[0] == 2 and m65[2] == 90 and m91[2] == 105 and abs(m65[3] - 20.05) < 0.005
       and abs(m91[3] - 25.90) < 0.005 and 0.099 <= 2 / m65[3] < 0.1 and 0.077 <= 2 / m91[3] < 0.078)

print('ALL OK' if ok_all else 'FAILURE')
sys.exit(0 if ok_all else 1)
