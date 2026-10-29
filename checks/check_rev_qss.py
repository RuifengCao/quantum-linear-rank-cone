#!/usr/bin/env python3
"""check_rev_qss.py -- independent check of the secret-sharing results (Lemma "Stabilizer schemes are linear schemes",
Proposition "Stabilizer schemes", Supplementary Note 10) and of the distances of the witnesses from the cones
(Supplementary Note 7.6).

Independent of the repository code: nothing from epr1kit/ or scripts/ is imported, and the paper source is not read.
The script reads only the data files data/s38_qss_certificates.json (rational dual certificates of the two linear
programs with Lemma 1 (LB), t = 2), data/s39_ideal_certificates.json (rational certificates of the forcing lemma of
part (k)) and data/s31_six_party.npz (key 'f9', the form f9 on the 128 subsets of A,B,C,W,X,Y,Z in mask order; used
for the margins).  The repository is found through the
environment variable EPR1_REPO, or else as the parent directory of the directory of this script (the repository itself
when the script sits in <repo>/checks/).  Requirements: Python 3.8+, numpy.  Running time: about half a minute.

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
 (e) the margins v/||g||_1 and v/||g||_2 for f9 (Steane state) and f12 (built from its definition; non-Fano state).
 (f) the schemes with information ratio 2: the decomposition scheme for Gamma_F over Z_d of R48 (R32 Proposition 23;
     in the paper it is superseded by the scheme of (h)) for d = 3, 5, 7, 9, 15, 21, 25, 27 (all 128 player sets,
     orders of projections computed exactly, and the CSS state), with a negative control; the binary scheme for
     Gamma_NF of Supplementary Note 10 (a secret of one qubit, shares of two; R32 Proposition 23), and its CSS state.
 (g) an upper bound on the distance of the Steane state from the cones of odd local dimension (Supplementary Note 7.6):
     the non-Fano vector relabeled by the 168 collineations of the Fano plane gives S_F plus 1 on a line and its
     complement, and (5/36) times the sum over the 7 lines lies within 1/12 of S_F in every entropy (exact fractions).
 (h) the scheme with information ratio 7/5 for Gamma_F of Proposition 1 (iii) (R32 Proposition 25, R48 addendum 2), a
     (lambda, omega) = (6, 1) decomposition in the sense of van Dijk, Jackson and Martin: seven cube schemes over Z_d,
     relabeled by collineations that map WXY to the seven lines, share the coordinates of a codeword of the zero-sum
     code; for d = 3, 5, 7, 9, 11, 13, 15, 21, 25, 27, 45, all 128 player sets (orders of projections computed exactly)
     and the CSS state; negative control without the uniform coefficient zeta.
 (i) the scheme with information ratio 7/5 for Gamma_NF, qubits and characteristic 2 (Proposition 1 (iii); R32 Lemma 27
     and Proposition 28, R48 addendum 3): the access structures Gamma^S, in which the lines in S and their complements
     exchange their status (Gamma^{} = Gamma_F, Gamma^{WXY} = Gamma_NF); the GF(4) code obtained from the binary cube
     code by changing two entries to alpha, its port Gamma^{WXY,ACX} and its matroid (AG(3,2) with four
     circuit-hyperplanes relaxed); the seven sub-schemes (the cube code and six relabelings of the GF(4) code by
     collineations that fix WXY), with lambda = 6 and omega = 1 on the minimal authorized and maximal unauthorized
     sets of Gamma_NF; the combined code over GF(4) with b = (1,1,1,1,1,alpha,alpha^2) and its binary version
     (restriction of scalars) on all 128 player sets, and the CSS states; two negative controls.
 (j) the finite part of the bound for decompositions of Gamma_F into ideal schemes of odd characteristic (R32
     Proposition 29): for supports of at most six players, a_j - c_j <= 4, 2, 1, 1, 0, 0; the orbit counts used for
     full supports; the bound is attained by the sub-schemes of (h).
 (k) ideal schemes (Proposition 1 (iv); R32 Proposition 30, R48 addendum 5): the 2 x 2 x 255 certificates of the
     forcing lemma (a polymatroid g on R and the seven players with g(R) = g(P) = 1 and the access equalities of
     Gamma_F, resp. Gamma_NF, equals the rank function of AG(3,2), resp. R8): every row rebuilt from its descriptor
     and checked to be admissible, every identity checked in exact arithmetic, the rank function a feasible point,
     negative controls; the orders of the projections of the normal form {(u,v,w,u+v,u+w,v+w,u+v+w)} over Z_{2^a},
     a = 1, 2, 3 (only WXY differs from AG(3,2) when a >= 2, with order 2^{6ma-2m}), also by direct enumeration over
     Z_4^2, and its kernel; the cube codes over Z_d, d = 9, 15, 21, 25, 27, 45, are ideal schemes for Gamma_NF.
 (l) the counts of Supplementary Note 10.6, Sec. II D and the Introduction: the erasure of a set E of qubits of Steane's code [[7,1,3]] is
     correctable (no nontrivial logical operator, i.e. no Hamming codeword outside the simplex code, is supported in
     E) exactly when E contains no line of the Fano plane, exactly when E is unauthorized for Gamma_F of (k); the
     correctable sets by size are 1, 7, 21, 28, 7, so 35 of them (the 28 triples that are not lines and the 7
     complements of lines) lie beyond the distance 3; every point lies on three lines, so summing the line
     inequalities of Proposition 1 (i) gives sum_P S(P) >= (28/3) S(R), and a five-qudit secret needs a share of at
     least ceil(20/3) = 7 qudits, as in the schemes of (h).
 (m) Lemma 3 (Fano lines in groups of odd order; R32 Proposition 31, R48 addendum 7) and the bound 4/3 of
     Proposition 1 (i) for every odd local dimension: the 26-term form built from the statement; on the Fano
     configuration in the labeling of the lemma (lines 014, 025, 036, 234, 135, 126, 456) the left side is 3 and the
     right side 4; the same vectors over GF(3), GF(5), GF(7) satisfy it in all 5040 labelings; tests on random and
     perturbed arrangements over GF(3), GF(5), Z_9, Z_25, Z_27 (>= 0) and Z_4 (< 0); the incidences used in the proof
     for all 28 pairs (line, player off it); the stronger form of Supplementary Note 10.2 and its relation to the
     lemma; the constant 11/2 for approximate schemes; the labels for the Fano port
     F of Jafari-Khazaei (classical abelian schemes) and the transitivity of its automorphisms on its four triples; an
     integer polymatroid on which the Shannon inequalities, all Ingleton instances and the four odd-characteristic
     forms of Supplementary Table S4 hold in all labelings, but the form of Lemma 3 is -37; and one that satisfies
     them and the access equations of the Fano port F with information ratio 23/22.
 (n) the schemes of Proposition 1 (iii) for Z_d qudits with d = 2r, r odd (R32 Proposition 31 (v), R48 addendum 7):
     the code over Z_d whose reduction modulo 2 is five binary cube codes (two vanishing coordinates added per player)
     and modulo r the code of (h) is a scheme for Gamma_F, and the code whose reduction modulo 2 is the binary code of
     (i) and modulo r ten cube codes (four vanishing coordinates added per player) is a scheme for Gamma_NF, both with
     ratio 7/5, for d = 6, 10, 18, 30 (all 128 player sets, orders of projections computed exactly, and the CSS states);
     two negative controls; and, for the ranges of the optimal ratios in Supplementary Note 10.3, the scheme over Z_6 for Gamma_F with a secret of
     dimension 2^7 3^5 (not a whole number of qudits) and information ratio 7 log 6/(7 log 2 + 5 log 3) < 4/3.
 (o) the numbers of the codes for large secrets of Sec. II D and Supplementary Note 10.7 (exact integer arithmetic):
     m = ceil(k log2 d) is the least m with 2^m >= d^k and n = ceil(m log_d 2) the least n with d^n >= 2^m; for d = 3
     and k = 4, 5 this gives (m, n) = (7, 5), (8, 6), hence 7n = 35, 42 qutrits, against ceil(28k/3) = 38, 47 for
     every qutrit stabilizer scheme with such a secret; and n/k -> 1 for d = 3, 5, 6, 7 (n/k < 1.01 at k = 1000).
 (p) the ports of the generalized cube codes (Supplementary Note 10.1): for t = 2, 3, 4 (M = t + 1) and p = 2, 3, 5, 7,
     the port at O of the code with columns (1, v), v in {0, e_i, c - e_i, c}, contracted at the player C, has the
     same authorized sets as the port at c of the configuration (JI) (c in the span), with the labels A_i and B_i
     exchanged; this is the access structure Gamma_p of Corollary 2.9 of Pena-Macias (Kybernetika 2023).
 (q) the planes and good directions used in the proof of Theorem 5 (Supplementary Note 10.7; R48 addendum 13):
     AG(3,2) has 14 planes and R8 12 (the coplanar quadruples of vertices of the real cube, 6 faces and 6 diagonal
     rectangles: all those of AG(3,2) but the tetrahedra RWXY and ABCZ); in both, every set of at most three points
     is independent (A); the good directions (any two pairs of points along them form a plane) are all seven for
     AG(3,2) and 100, 010, 001, 111 for R8, and along the face diagonals the two unions of pairs that are not planes
     of R8 are the tetrahedra; the facts (G1)-(G4) of the proof; the circuits through R are R plus the minimal
     authorized sets (7 and 11), and every element lies in one, so both matroids are connected; the ranks of the 23
     sets RT used in the proof of Supplementary Lemma 10.4, which differ between the two matroids only on RWXY.
 (r) erasures (Supplementary Note 10.6; R48 addendum 13): the uncorrectable erasures of Steane's code and of the
     color code of odd local dimension, i.e. the authorized sets, by size (7, 28, 21, 7, 1 and 6, 29, 21, 7, 1 for
     3..7 sites); P_color - P_Steane = -eps^3 (1-eps)^3 (1-2 eps) for the erasure channel (exact polynomials); and
     an exhaustive search over the qubit stabilizer codes that encode one qubit into seven with distance >= 3, as
     graph states on R and the seven sites: the 1,044 graphs on seven vertices up to isomorphism, generated here by
     marking the orbits of the 5040 permutations and counted with Burnside's lemma, times the 127 nonempty
     neighborhoods of R; the GF(2) rank formula for graph-state entropies, cross-checked on explicit state vectors;
     3,414 graph states give codes of distance >= 3, each with at least 7 uncorrectable erasures of three sites.
The validity of Lemma 1 and of the results of Jafari-Khazaei is an input; it is only tested here. So is the theorem
that every qubit stabilizer state is local-Clifford equivalent to a graph state (Van den Nest, Dehaene and De Moor,
Phys. Rev. A 69, 022316 (2004), Theorem 1), used in (r). Lemma 3 is proved in
the paper (Supplementary Note 10.2); it is only tested here.
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


f9arr = np.load(os.path.join(DATA, 's31_six_party.npz'))['f9']
F9 = {m: int(f9arr[m]) for m in range(128) if f9arr[m]}
hF = [a + b for a, b in zip(rF, dual(rF, 7))]
SF7, SNF7 = conn(rF, 7), conn(rNF, 7)
report('margin inputs: f9 has 67 terms and f9(h_F) = -3 with h_F = r_F + r_F^* computed here; the lift of the '
       'Steane state S_F = lambda_{r_F} is h_F',
       len(F9) == 67 and ev(F9, hF) == -3 and lift(SF7, 7) == hF)


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
print('     %-50s %3s %6s %9s %6s %10s %9s' % ('form, witness', 'v', '||g||1', '||g||2', 'terms', 'v/||g||1', 'v/||g||2'))
for k_, (v, nt, l1, l2, sq) in margins.items():
    print('     %-50s %3d %6d %9.3f %6d %10s %9.5f' % (k_, v, l1, l2, nt, Fraction(v, l1), v / l2))
m9, m12 = margins['f9, Steane state, purifier Z'], margins['f12, non-Fano state, purifier Z']
report('f9: v = 3, ||g||_1 = 2856, ||g||_2 = sqrt(348372) = 590.2 for every choice of purifier; 1/952 = 1.05e-3 and '
       'Euclidean %.5f >= 5.0e-3' % (3 / m9[3]),
       m9[0] == 3 and m9[2] == 2856 and m9[4] == 348372 and all(x == (3, 2856, m9[3]) for x in [a[:3] for a in f9_all])
       and Fraction(3, 2856) == Fraction(1, 952) and 5.0e-3 <= 3 / m9[3] < 5.1e-3)
report('f12: v = 1, ||g||_1 = 129, ||g||_2 = sqrt(839) = 28.97; 1/129 and Euclidean %.4f >= 0.034' % (1 / m12[3]),
       m12[0] == 1 and m12[2] == 129 and m12[4] == 839 and 0.034 <= 1 / m12[3] < 0.035)

# ============================================================================== (f) schemes with information ratio 2
print('(f) schemes with information ratio 2 (Lemma "Linear schemes give stabilizer schemes", Supplementary Note 10)')
ELL = {'A': 'BCY', 'B': 'ACX', 'C': 'ABW', 'Z': 'WXY', 'W': 'AYZ', 'X': 'CWZ', 'Y': 'BXZ'}
INV = {v: k for k, v in VERT.items()}
report('the bijection q -> l_q of the paper: l_q is a line, q not on it, all seven lines used; on A, B, C, Z it is the '
       'polarity q -> {x : q.x = 0}',
       all(ELL[q] in LINES7 and q not in ELL[q] for q in L7) and len(set(ELL.values())) == 7
       and all(set(ELL[q]) == {P for P in L7 if sum(a * b for a, b in zip(VERT[q], VERT[P])) % 2 == 0}
               for q in 'ABCZ'))


def gf2_apply(M, v):
    return tuple(sum(M[i][j] * v[j] for j in range(3)) % 2 for i in range(3))


COLL = []
for e in itertools.product((0, 1), repeat=9):
    M = (e[0:3], e[3:6], e[6:9])
    img = {gf2_apply(M, VERT[P]) for P in L7}
    if len(img) == 7 and (0, 0, 0) not in img:
        COLL.append({P: INV[gf2_apply(M, VERT[P])] for P in L7})
SIG = {}
for q in L7:
    SIG[q] = next(c for c in COLL if {c[P] for P in 'WXY'} == set(ELL[q]))
report('168 collineations of the Fano plane; for each q one with sigma_q(WXY) = l_q',
       len(COLL) == 168 and all(all(''.join(sorted(SIG[q][P] for P in l)) in [''.join(sorted(m)) for m in LINES7]
                                     for l in LINES7) for q in L7))


def decomposition_code(cvec):
    """coordinates of the combined code as integer functionals on the parameters (s in Z^3, a_q in Z^3 for the seven
    q); returns {owner: [functional, ...]} with owners 'R' and the players."""
    cols = {'R': [[1 if j == i else 0 for j in range(24)] for i in range(3)]}
    for P in L7:
        cols[P] = []
        for jq, q in enumerate(L7):
            if q == P:
                continue
            x = VERT[next(Q for Q in L7 if SIG[q][Q] == P)]          # the cube point of P in the scheme of q
            f = [0] * 24
            for i in range(3):
                f[i] = cvec[q][i]
                f[3 + 3 * jq + i] = x[i]
            cols[P].append(f)
    return cols


def factor_small(d):
    out, x, f = {}, d, 2
    while x > 1:
        while x % f == 0:
            out[f] = out.get(f, 0) + 1
            x //= f
        f += 1
    return out


def order_zd(functionals, d):
    """order of the projection: the subgroup of Z_d^k generated by the images of the 24 parameter vectors"""
    if not functionals:
        return 1
    vecs = [[f[i] for f in functionals] for i in range(24)]
    o = 1
    for pr, a in factor_small(d).items():
        o *= pr ** zn_log_order(vecs, len(functionals), pr ** a, pr)
    return o


def scheme_report(cols, d, minsets, label, ratio_num):
    owners = ['R'] + list(L7)
    def o_of(mask8):
        fs = [f for i, w in enumerate(owners) if mask8 >> i & 1 for f in cols[w]]
        return order_zd(fs, d)
    O = {m: o_of(m) for m in range(256)}
    oR = O[1]
    auth = lambda T: any(T & mm == mm for mm in minsets)
    acc_ok = True
    for T in range(128):
        mT = T << 1
        acc_ok &= (O[mT | 1] == O[mT]) if auth(T) else (O[mT | 1] == O[mT] * oR)
    full = 255
    E = {m: Fraction(O[m] * O[full ^ m], O[full]) for m in range(256)}
    css_ok = E[1] == oR
    for T in range(128):
        mT = T << 1
        I = E[1] * E[mT] / E[mT | 1]
        css_ok &= (I == E[1] ** 2) if auth(T) else (I == 1)
    shares_ok = all(O[1 << (i + 1)] == d ** ratio_num for i in range(7)) and oR == d ** 3
    return acc_ok, css_ok, shares_ok


CVEC = {q: VERT[q] for q in L7}
for d in (3, 5, 7, 9, 15, 21, 25, 27):
    cols = decomposition_code(CVEC)
    acc_ok, css_ok, sh_ok = scheme_report(cols, d, MIN_F, 'F', 6)
    report('Z_%d: the decomposition scheme is a scheme for Gamma_F on all 128 player sets (orders of projections from '
           'Smith-type elimination); secret Z_d^3, each player 6 coordinates of full rank; its CSS state has S(R) = 3 log d '
           'and I(R;T) = 2 S(R) on Gamma_F and 0 elsewhere' % d, acc_ok and css_ok and sh_ok)
bad = decomposition_code({q: (1, 0, 0) for q in L7})
acc_bad, _, _ = scheme_report(bad, 3, MIN_F, 'F', 6)
report('negative control: with c_q = (1,0,0) for every q the combined code is not a scheme for Gamma_F', not acc_bad)

NF_G = ['100000001111111', '010000001111001', '001000011011011', '000100000010011', '000010001001110',
        '000001000010110', '000000100101011']        # the generator matrix printed in Supplementary Note 10
NF_cols = {'R': [[int(r[0]) for r in NF_G]]}
for i, P in enumerate(L7):
    NF_cols[P] = [[int(r[1 + 2 * i + j]) for r in NF_G] for j in range(2)]


def gf2_scheme(cols, minsets):
    owners = ['R'] + list(L7)
    rk = {}
    for m in range(256):
        vs = [v for i, w in enumerate(owners) if m >> i & 1 for v in cols[w]]
        rk[m] = rank_mod_p(vs, 2) if vs else 0
    auth = lambda T: any(T & mm == mm for mm in minsets)
    ok = rk[1] == 1 and all(rk[1 << (i + 1)] == 2 for i in range(7))
    for T in range(128):
        mT = T << 1
        ok &= (rk[mT | 1] == rk[mT]) if auth(T) else (rk[mT | 1] == rk[mT] + 1)
    S = {m: rk[m] + rk[255 ^ m] - rk[255] for m in range(256)}
    for T in range(128):
        mT = T << 1
        I = S[1] + S[mT] - S[mT | 1]
        ok &= (I == 2 * S[1]) if auth(T) else (I == 0)
    return ok


report('the binary code of Supplementary Note 10 (secret 1 coordinate, 2 per player) is a scheme for Gamma_NF on all 128 '
       'player sets, and its CSS state is a qubit scheme for Gamma_NF with ratio 2', gf2_scheme(NF_cols, MIN_NF))
report('negative control: the same code is not a scheme for Gamma_F', not gf2_scheme(NF_cols, MIN_F))


print('(g) an upper bound on the distance of the Steane state from the cones of odd local dimension')

LINE_MASKS = {mask(l) for l in LINES7}


def is_collineation(perm):
    """perm maps the i-th point of L7 to the point perm[i]; collineations map lines to lines."""
    return all(mask(''.join(L7[perm[L7.index(c)]] for c in l)) in LINE_MASKS for l in LINES7)


def relabel_vec(S, perm):
    """S' with S'(sigma(X)) = S(X) for the permutation sigma of the seven systems given by perm."""
    out = [0] * 128
    for X in range(128):
        Y = 0
        for i in range(7):
            if X >> i & 1:
                Y |= 1 << perm[i]
        out[Y] = S[X]
    return out


colls = [perm for perm in itertools.permutations(range(7)) if is_collineation(perm)]
wxy = mask('WXY')
relabeled = {}
for perm in colls:
    img = 0
    for i in range(7):
        if wxy >> i & 1:
            img |= 1 << perm[i]
    relabeled.setdefault(img, relabel_vec(SNF7, perm))
delta_ok = True
for lm, vec in relabeled.items():
    for X in range(128):
        extra = 1 if X in (lm, 127 ^ lm) else 0
        delta_ok &= vec[X] == SF7[X] + extra
report('the Fano plane has 168 collineations; for each of the 7 lines l, the non-Fano vector relabeled so that WXY goes '
       'to l equals S_F plus 1 on l and on its complement', len(colls) == 168 and len(relabeled) == 7 and delta_ok)
comb = [Fraction(5, 36) * sum(vec[X] for vec in relabeled.values()) for X in range(128)]
devs = {comb[X] - SF7[X] for X in range(1, 127)}
dev_max = max(abs(x) for x in devs)
report('(5/36) times the sum of the 7 relabeled non-Fano vectors differs from S_F by %s on the nonempty proper subsets, '
       'so it lies within %s of S_F in every entropy (and 1/952 < 1/12)'
       % (', '.join(str(x) for x in sorted(devs)), dev_max),
       devs == {Fraction(-1, 36), Fraction(-1, 18), Fraction(-1, 12), Fraction(1, 12)} and dev_max == Fraction(1, 12)
       and Fraction(1, 952) < Fraction(1, 12))
on_line = [any(X in (lm, 127 ^ lm) for lm in relabeled) for X in range(128)]
best = min(max(abs(c * (SF7[X] + Fraction(int(on_line[X]), 7)) - SF7[X]) for X in range(1, 127))
           for c in [Fraction(k, 360) for k in range(300, 361)])
report('within the family c * (1/7) * sum, the coefficient c = 35/36 (that is, 5/36 for the sum) is optimal on a grid of '
       'step 1/360: best maximal deviation %s' % best, best == Fraction(1, 12))


# ============================================================================== (h) a scheme with information ratio 7/5
print('(h) a scheme with information ratio 7/5 for Gamma_F and every odd d (Proposition 1 (iii))')
SIG_L = {l: next(c for c in COLL if {c[P] for P in 'WXY'} == set(l)) for l in LINES7}
report('for each of the 7 lines l a collineation sigma_l with sigma_l(WXY) = l',
       all(''.join(sorted(SIG_L[l][P] for P in 'WXY')) == ''.join(sorted(l)) for l in LINES7))
B_RAMP = (1, 1, 1, 1, 1, -1, -4)
NPAR_H = 5 + 1 + 21            # the secret u in Z^5, zeta, and a_l in Z^3 for the seven lines


def ramp_code(b, use_zeta=True):
    """functionals on the parameters: the secret u (5 coordinates at R) and, for each player, the value of the
    affine function s_l + a_l . x of the scheme of the line l at the player's cube vertex, where
    s = u_1 (e_2 - e_7) + ... + u_5 (e_6 - e_7) + zeta b."""
    comb = []
    for j in range(7):
        f = [0] * 6
        if 1 <= j <= 5:
            f[j - 1] = 1
        if j == 6:
            f[0:5] = [-1] * 5
        if use_zeta:
            f[5] = b[j]
        comb.append(f)
    cols = {'R': [[1 if k == i else 0 for k in range(NPAR_H)] for i in range(5)]}
    for P in L7:
        cols[P] = []
        for j, l in enumerate(LINES7):
            x = VERT[next(Q for Q in L7 if SIG_L[l][Q] == P)]          # the cube vertex of P in the scheme of l
            f = comb[j][:] + [0] * 21
            for i in range(3):
                f[6 + 3 * j + i] = x[i]
            cols[P].append(f)
    return cols


def order_zd_h(functionals, d):
    if not functionals:
        return 1
    vecs = [[f[i] for f in functionals] for i in range(NPAR_H)]
    o = 1
    for pr, a in factor_small(d).items():
        o *= pr ** zn_log_order(vecs, len(functionals), pr ** a, pr)
    return o


def scheme_report_h(cols, d):
    owners = ['R'] + list(L7)
    O = {}
    for m8 in range(256):
        O[m8] = order_zd_h([f for i, w in enumerate(owners) if m8 >> i & 1 for f in cols[w]], d)
    oR = O[1]
    auth = lambda T: any(T & mm == mm for mm in MIN_F)
    acc = all((O[(T << 1) | 1] == O[T << 1]) if auth(T) else (O[(T << 1) | 1] == O[T << 1] * oR) for T in range(128))
    full = 255
    E = {m8: Fraction(O[m8] * O[full ^ m8], O[full]) for m8 in range(256)}
    css = E[1] == oR and all((E[1] * E[T << 1] / E[(T << 1) | 1] == E[1] ** 2) if auth(T)
                             else (E[1] * E[T << 1] / E[(T << 1) | 1] == 1) for T in range(128))
    sizes = oR == d ** 5 and all(O[1 << (i + 1)] == d ** 7 for i in range(7))
    return acc, css, sizes


zero_sum_basis_ok = sum(B_RAMP) == 0 and B_RAMP[0] == 1
report('b = (1,1,1,1,1,-1,-4): coordinate sum 0, first coordinate 1, so b and e_2 - e_7, ..., e_6 - e_7 form a basis of '
       'the zero-sum code; every coordinate of b is a unit modulo every odd d', zero_sum_basis_ok
       and all(math.gcd(x, d) == 1 for x in B_RAMP for d in (3, 5, 7, 9, 11, 13, 15, 21, 25, 27, 45)))
for d in (3, 5, 7, 9, 11, 13, 15, 21, 25, 27, 45):
    acc, css, sizes = scheme_report_h(ramp_code(B_RAMP), d)
    report('Z_%d: the seven relabeled cube schemes with the zero-sum code are a scheme for Gamma_F on all 128 player '
           'sets; secret Z_d^5, each player 7 coordinates of full rank (ratio 7/5); the CSS state has S(R) = 5 log d and '
           'I(R;T) = 2 S(R) on Gamma_F and 0 elsewhere' % d, acc and css and sizes)
acc_bad, _, _ = scheme_report_h(ramp_code(B_RAMP, use_zeta=False), 3)
report('negative control: without zeta (s = the five secret coordinates only) the code is not a scheme for Gamma_F '
       '(a line complement learns one coordinate of the secret)', not acc_bad)

# ============================================================================ (i) ratio 7/5 for Gamma_NF and qubits
print('(i) a scheme with information ratio 7/5 for Gamma_NF, qubits and characteristic 2 (Proposition 1 (iii))')
# GF(4) = {0, 1, 2, 3} with 2 = alpha, 3 = alpha^2 = alpha + 1; addition is XOR.
G4_EXP = [1, 2, 3]
G4_LOG = {1: 0, 2: 1, 3: 2}


def g4_mul(a, b):
    return 0 if a == 0 or b == 0 else G4_EXP[(G4_LOG[a] + G4_LOG[b]) % 3]


def g4_inv(a):
    return G4_EXP[(-G4_LOG[a]) % 3]


def rank_gf4(rows):
    """rank over GF(4) of a list of vectors with entries in {0, 1, 2, 3}"""
    M = [list(r) for r in rows if any(r)]
    rk = 0
    ncol = len(M[0]) if M else 0
    for col in range(ncol):
        piv = next((i for i in range(rk, len(M)) if M[i][col]), None)
        if piv is None:
            continue
        M[rk], M[piv] = M[piv], M[rk]
        iv = g4_inv(M[rk][col])
        M[rk] = [g4_mul(iv, x) for x in M[rk]]
        for i in range(len(M)):
            if i != rk and M[i][col]:
                f = M[i][col]
                M[i] = [x ^ g4_mul(f, y) for x, y in zip(M[i], M[rk])]
        rk += 1
    return rk


def gamma_swap(S):
    """Gamma^S: a set of players is authorized iff it contains a line not in S, or a line in S and one more player,
    or the complement of a line in S (each line in S and its complement exchange their status)."""
    Sm = {mask(l) for l in S}
    def auth(T):
        if any(T & mask(l) == mask(l) and mask(l) not in Sm for l in LINES7):
            return True
        if any(T & m == m and T != m for m in Sm):
            return True
        return any(T & (127 ^ m) == (127 ^ m) for m in Sm)
    return auth


def minimal_of(auth):
    return sorted(T for T in range(1, 128) if auth(T) and not any(auth(T ^ (1 << i)) for i in range(7) if T >> i & 1))


def expected_minimal(S):
    Sm = [mask(l) for l in S]
    out = [mask(l) for l in LINES7 if mask(l) not in Sm] + [127 ^ m for m in Sm]
    out += [m | (1 << i) for m in Sm for i in range(7) if not m >> i & 1]
    return sorted(out)


ALL_S = [S for k in range(8) for S in itertools.combinations(LINES7, k)]
report('Gamma^{} = Gamma_F and Gamma^{WXY} = Gamma_NF (minimal authorized sets as in (a)); for all 128 sets S of lines, '
       'Gamma^S is self-dual (T authorized iff its complement is not) and its minimal authorized sets are the lines not '
       'in S, the complements of the lines in S and the sets l + x with l in S, x not on l',
       minimal_of(gamma_swap([])) == MIN_F and minimal_of(gamma_swap(['WXY'])) == MIN_NF
       and all(all(gamma_swap(S)(T) != gamma_swap(S)(127 ^ T) for T in range(128)) for S in ALL_S)
       and all(minimal_of(gamma_swap(S)) == expected_minimal(S) for S in ALL_S))
# generator matrices: rows = coordinates, columns R, A, B, C, W, X, Y, Z
E8 = 'RABCWXYZ'
CUBE_SYS = [[1, 0, 0, 0, 1, 1, 1, 0],
            [0, 1, 0, 0, 1, 1, 0, 1],
            [0, 0, 1, 0, 1, 0, 1, 1],
            [0, 0, 0, 1, 0, 1, 1, 1]]
SWAP_G4 = [[1, 0, 0, 0, 1, 1, 1, 2],
           [0, 1, 0, 0, 1, 1, 0, 1],
           [0, 0, 1, 0, 1, 2, 1, 1],
           [0, 0, 0, 1, 0, 1, 1, 1]]
CUBE_AFF = [[1] * 8] + [[0] + [VERT[P][i] for P in 'ABCWXYZ'] for i in range(3)]   # affine functions, R = origin
report('the matrix with columns R, A, B, C = unit vectors and W, X, Y, Z = (1,1,1,0), (1,1,0,1), (1,0,1,1), '
       '(0,1,1,1) generates the binary cube code (same row space over GF(2) as the affine functions); the GF(4) matrix '
       'differs from it in two entries, alpha at (row 1, Z) and (row 3, X)',
       rank_mod_p(CUBE_SYS, 2) == 4 and rank_mod_p(CUBE_AFF, 2) == 4 and rank_mod_p(CUBE_SYS + CUBE_AFF, 2) == 4
       and sum(1 for i in range(4) for j in range(8) if SWAP_G4[i][j] != CUBE_SYS[i][j]) == 2
       and SWAP_G4[0][7] == 2 and SWAP_G4[2][5] == 2)


def columns(G):
    return {E8[j]: tuple(G[i][j] for i in range(4)) for j in range(8)}


def port_of(cols, rankf):
    out = {}
    for T in range(128):
        vs = [cols[P] for i, P in enumerate(L7) if T >> i & 1]
        out[T] = bool(vs) and rankf(vs + [cols['R']]) == rankf(vs)
    return out


def rk8(cols, X, rankf):
    vs = [cols[c] for c in X]
    return rankf(vs) if vs else 0


C_CUBE = columns(CUBE_SYS)
C_SWAP = columns(SWAP_G4)
port_cube = port_of(C_CUBE, rank_gf4)
port_swap = port_of(C_SWAP, rank_gf4)
target_swap = gamma_swap(['WXY', 'ACX'])
report('over GF(4) the cube code is an ideal scheme for Gamma_F (= Gamma^{}), and the GF(4) code is an ideal scheme for '
       'Gamma^{WXY,ACX}: WXY and ACX unauthorized, ABCZ and BWYZ authorized, all other sets as in Gamma_F (128 sets)',
       all(port_cube[T] == gamma_swap([])(T) for T in range(128))
       and all(port_swap[T] == target_swap(T) for T in range(128))
       and all(rank_gf4([C_SWAP[P]]) == 1 for P in E8))
relaxed = [set('RWXY'), set('ABCZ'), set('RACX'), set('BWYZ')]
mat_ok = True
n_planes = 0
for k in range(9):
    for X in itertools.combinations(E8, k):
        rc, rs = rk8(C_CUBE, X, rank_gf4), rk8(C_SWAP, X, rank_gf4)
        n_planes += (k == 4 and rc == 3)
        mat_ok &= rs == (4 if (k == 4 and set(X) in relaxed) else rc)
report('matroids: the cube code gives AG(3,2) (14 four-element circuits); the GF(4) code gives AG(3,2) with the four '
       'circuit-hyperplanes RWXY, ABCZ, RACX, BWYZ relaxed (all 256 ranks)', mat_ok and n_planes == 14)
COLL_WXY = [g for g in COLL if {g[P] for P in 'WXY'} == set('WXY')]
others = [l for l in LINES7 if l != 'WXY']
report('24 collineations fix the line WXY and act transitively on the other six lines',
       len(COLL_WXY) == 24 and {''.join(sorted(g[P] for P in 'ACX')) for g in COLL_WXY} == {''.join(sorted(l))
                                                                                         for l in others})


def relabel_cols(cols, g):
    out = {'R': cols['R']}
    for P in L7:
        out[g[P]] = cols[P]
    return out


PIECES_NF = {'WXY': C_CUBE}
for l in others:
    g = next(g for g in COLL_WXY if {g[P] for P in 'ACX'} == set(l))
    PIECES_NF[l] = relabel_cols(C_SWAP, g)
PORTS_NF = {l: port_of(PIECES_NF[l], rank_gf4) for l in LINES7}
report('the seven sub-schemes: the cube code for Gamma^{} and, for each line l other than WXY, the relabeled GF(4) code '
       'for Gamma^{WXY,l} (128 sets each)',
       all(all(PORTS_NF[l][T] == (gamma_swap([]) if l == 'WXY' else gamma_swap(['WXY', l]))(T) for T in range(128))
           for l in LINES7))
MAXU_NF = sorted(127 ^ m for m in MIN_NF)
lam = min(sum(PORTS_NF[l][M] for l in LINES7) for M in MIN_NF)
omg = max(sum(PORTS_NF[l][U] for l in LINES7) for U in MAXU_NF)
report('(lambda, omega)-decomposition of Gamma_NF: every one of its 11 minimal authorized sets is authorized in at least '
       '6 sub-schemes, every one of its 11 maximal unauthorized sets in at most 1', lam == 6 and omg == 1,
       'lambda = %d, omega = %d' % (lam, omg))
B4 = (1, 1, 1, 1, 1, 2, 3)          # (1,1,1,1,1,alpha,alpha^2)
s4 = 0
for x in B4:
    s4 ^= x
report('b = (1,1,1,1,1,alpha,alpha^2) over GF(4): nonzero coordinates, coordinate sum 0, first coordinate 1',
       all(B4) and s4 == 0 and B4[0] == 1)
NPAR_I = 5 + 1 + 21            # the secret u in GF(4)^5, zeta, and 3 random coordinates for each sub-scheme


def nf_code(b, pieces):
    """GF(4) functionals on the parameters: R holds u; player P holds, for each line l (numbered as in LINES7), the
    coordinate of P in the sub-scheme of l, whose message is (s_l, three random coordinates), with
    s = u_1 (e_2 + e_7) + ... + u_5 (e_6 + e_7) + zeta b (characteristic 2)."""
    comb = []
    for j in range(7):
        f = [0] * 6
        if 1 <= j <= 5:
            f[j - 1] = 1
        if j == 6:
            f[0:5] = [1] * 5
        f[5] = b[j]
        comb.append(f)
    cols = {'R': [[1 if k == i else 0 for k in range(NPAR_I)] for i in range(5)]}
    for P in L7:
        cols[P] = []
        for j, l in enumerate(LINES7):
            v = pieces[l][P]
            f = [g4_mul(v[0], x) for x in comb[j]] + [0] * 21
            for t in range(3):
                f[6 + 3 * j + t] = v[1 + t]
            cols[P].append(f)
    return cols


def to_binary(fs):
    """restriction of scalars: a GF(4)-linear functional on GF(4)^27 gives two GF(2)-linear functionals on GF(2)^54
    (each element written as c0 + c1 alpha)"""
    out = []
    for f in fs:
        for comp in (0, 1):
            row = []
            for a in f:
                for basis in (1, 2):
                    pr = g4_mul(a, basis)
                    row.append((pr >> comp) & 1)
            out.append(row)
    return out


def nf_report(cols, field):
    owners = ['R'] + list(L7)
    if field == 4:
        rk = lambda fs: rank_gf4(fs) if fs else 0
        sec, sh = 5, 7
    else:
        rk = lambda fs: rank_mod_p(to_binary(fs), 2) if fs else 0
        sec, sh = 10, 14
    O = {m8: rk([f for i, w in enumerate(owners) if m8 >> i & 1 for f in cols[w]]) for m8 in range(256)}
    auth = gamma_swap(['WXY'])
    acc = all((O[(T << 1) | 1] == O[T << 1]) if auth(T) else (O[(T << 1) | 1] == O[T << 1] + O[1]) for T in range(128))
    S = {m8: O[m8] + O[255 ^ m8] - O[255] for m8 in range(256)}
    css = S[1] == O[1] and all((S[1] + S[T << 1] - S[(T << 1) | 1] == 2 * S[1]) if auth(T)
                               else (S[1] + S[T << 1] - S[(T << 1) | 1] == 0) for T in range(128))
    sizes = O[1] == sec and all(O[1 << (i + 1)] == sh for i in range(7))
    return acc, css, sizes


acc4, css4, sz4 = nf_report(nf_code(B4, PIECES_NF), 4)
report('GF(4): the seven sub-schemes with the zero-sum code are a scheme for Gamma_NF on all 128 player sets; secret '
       'GF(4)^5, each player 7 coordinates of full rank (ratio 7/5); its CSS state (Galois qudits of dimension 4) has '
       'I(R;T) = 2 S(R) on Gamma_NF and 0 elsewhere', acc4 and css4 and sz4)
acc2, css2, sz2 = nf_report(nf_code(B4, PIECES_NF), 2)
report('binary code by restriction of scalars: a scheme for Gamma_NF with secret of 10 and shares of 14 coordinates; its '
       'CSS state is a qubit scheme with ratio 7/5 (all 128 player sets, ranks over GF(2))', acc2 and css2 and sz2)
bad1, _, _ = nf_report(nf_code((1, 1, 1, 1, 1, 1, 0), PIECES_NF), 4)
bad2, _, _ = nf_report(nf_code(B4, {l: C_CUBE for l in LINES7}), 4)
report('negative controls: with b = (1,1,1,1,1,1,0) (zero sum, b_7 not a unit), and with the cube code as all seven '
       'sub-schemes, the code is not a scheme for Gamma_NF', not bad1 and not bad2)

# ============================================================== (j) decompositions of Gamma_F into ideal schemes
print('(j) decompositions of Gamma_F into ideal schemes give no ratio below 7/5 (Supplementary Note 10)')
LM7 = [mask(l) for l in LINES7]
CM7 = [127 ^ m for m in LM7]
best = {}
viol = 0
for Pm in range(1, 127):
    sz = popcount(Pm)
    tr = [m & Pm for m in LM7]
    cand = [i for i in range(7) if tr[i]]
    for k in range(len(cand) + 1):
        for A in itertools.combinations(cand, k):
            bmin = sum(1 for c in CM7 if any((c & Pm) & tr[i] == tr[i] for i in A))
            best[sz] = max(best.get(sz, -99), k - bmin)
            viol += 5 * sz < 7 * (k - bmin)
report('for every set P_j of at most six players and every set of lines whose traces on P_j are authorized, the '
       'complements forced by monotonicity give a_j - c_j <= 4, 2, 1, 1, 0, 0 for n_j = |P_j| = 6, ..., 1, so '
       '5 n_j >= 7 (a_j - c_j)', viol == 0 and best == {6: 4, 5: 2, 4: 1, 3: 1, 2: 0, 1: 0}, str(dict(sorted(best.items()))))
stab = [g for g in COLL if {g[P] for P in 'WXY'} == set('WXY')]
report('the 168 collineations act transitively on the lines and on their complements, and the stabilizer of a line has '
       'order 24 (so the 168 relabelings of a sub-scheme authorize each line 24 a_j times and each complement 24 c_j '
       'times)', len(stab) == 24 and len({''.join(sorted(g[P] for P in 'WXY')) for g in COLL}) == 7)
pr48 = port_of(columns([[x % 3 for x in r] for r in CUBE_AFF]), lambda vs: rank_mod_p(vs, 3))
a48 = sum(pr48[m] for m in LM7)
b48 = sum(pr48[c] for c in CM7)
report('the bound is attained: the cube code over GF(3) (a sub-scheme of the scheme of (h)) has n = 7, a = 6, c = 1, '
       'so 5 n = 35 = 7 (a - c)', a48 == 6 and b48 == 1)

# ============================================================================ (k) ideal schemes (Proposition 1 (iv))
print('(k) ideal schemes (Proposition 1 (iv)): forcing of the lift and the normal form over Z_{2^a}')
# Data: data/s39_ideal_certificates.json, for 'F' (Gamma_F, matroid AG(3,2)) and 'NF' (Gamma_NF, matroid R8), and for
# every nonempty set X of the eight elements (mask, players bits 0..6 in the order A,B,C,W,X,Y,Z, R = bit 7) a 'lower'
# and an 'upper' certificate {'ineq': [[descriptor, y], ...], 'eq': [[descriptor, z], ...]} with rational y, z.
# Descriptors: ['mono', i]  : g(full) - g(full - i) >= 0;   ['sub', i, j, K] : I(i; j | K) >= 0  (i < j, K disjoint);
#              ['norm', i]  : g(i) = 1;  ['auth', T] : g(RT) - g(T) = 0;  ['unauth', T] : g(RT) - g(T) = 1.
# Lower certificate:  sum y a + sum z e = +1_X and sum z b = r_M(X)   (so g(X) >= r_M(X) for every feasible g);
# upper certificate:  sum y a + sum z e = -1_X and sum z b = -r_M(X)  (so g(X) <= r_M(X)).
CERT39 = json.load(open(os.path.join(DATA, 's39_ideal_certificates.json')))
RM39 = {'F': cube[2]['r'], 'NF': rQ}                     # AG(3,2) over GF(2); R8 = the cube over Q, from (a)
AUTH39 = {'F': {T: ACC[2][T] == 'auth' for T in range(1, 128)}, 'NF': {T: ACC[3][T] == 'auth' for T in range(1, 128)}}
report('the ports at R of AG(3,2) and of R8 are Gamma_F and Gamma_NF: T authorized iff r(RT) = r(T), and r(RT) = r(T) + 1 '
       'for every other nonempty set T of players',
       all((RM39[c][T | 128] == RM39[c][T]) == AUTH39[c][T] and (AUTH39[c][T] or RM39[c][T | 128] == RM39[c][T] + 1)
           for c in ('F', 'NF') for T in range(1, 128)))


def row39_ineq(desc):
    if desc[0] == 'mono' and len(desc) == 2 and 0 <= desc[1] < 8:
        return {255: 1, 255 ^ (1 << desc[1]): -1}
    if desc[0] == 'sub' and len(desc) == 4:
        i, j, K = desc[1], desc[2], desc[3]
        if 0 <= i < j < 8 and 0 <= K < 256 and not K & (1 << i | 1 << j):
            out = {}
            for Y, c in ((K | 1 << i, 1), (K | 1 << j, 1), (K | 1 << i | 1 << j, -1), (K, -1)):
                if Y:
                    out[Y] = out.get(Y, 0) + c
            return out
    return None


def row39_eq(desc, case):
    if desc[0] == 'norm' and len(desc) == 2 and 0 <= desc[1] < 8:
        return {1 << desc[1]: 1}, 1
    if desc[0] in ('auth', 'unauth') and len(desc) == 2 and 0 < desc[1] < 128 \
            and AUTH39[case][desc[1]] == (desc[0] == 'auth'):
        return {desc[1] | 128: 1, desc[1]: -1}, (0 if desc[0] == 'auth' else 1)
    return None


def check39(cert, X, sense, value, case):
    tot, bound = {}, Fraction(0)
    for desc, y in cert['ineq']:
        a, y = row39_ineq(desc), Fraction(y)
        if a is None or y < 0:
            return False
        for Y, c in a.items():
            tot[Y] = tot.get(Y, Fraction(0)) + y * c
    for desc, z in cert['eq']:
        e, z = row39_eq(desc, case), Fraction(z)
        if e is None:
            return False
        for Y, c in e[0].items():
            tot[Y] = tot.get(Y, Fraction(0)) + z * c
        bound += z * e[1]
    return {Y: c for Y, c in tot.items() if c} == {X: Fraction(sense)} and bound == sense * value


stats39 = {}
for case in ('F', 'NF'):
    good = sorted(CERT39[case]) == sorted(str(X) for X in range(1, 256)) and all(
        check39(CERT39[case][str(X)][key], X, s, RM39[case][X], case) for X in range(1, 256)
        for key, s in (('lower', 1), ('upper', -1)))
    feas = is_polymatroid(RM39[case], 8) and all(RM39[case][1 << i] == 1 for i in range(8))
    nrows = [len(CERT39[case][str(X)][k]['ineq']) + len(CERT39[case][str(X)][k]['eq'])
             for X in range(1, 256) for k in ('lower', 'upper')]
    dens = max(Fraction(v).denominator for X in range(1, 256) for k in ('lower', 'upper')
               for part in ('ineq', 'eq') for _, v in CERT39[case][str(X)][k][part])
    stats39[case] = (max(nrows), dens)
    report('%s: the 510 certificates verify exactly (every row rebuilt from its descriptor and admissible for %s), and '
           'the rank function of %s is a feasible point; so every polymatroid g with g(R) = g(P) = 1 and the access '
           'equalities equals r_M (Lemma "Forcing")' % (case, 'Gamma_F' if case == 'F' else 'Gamma_NF',
                                                          'AG(3,2)' if case == 'F' else 'R8'),
           good and feas, 'at most %d rows per certificate, largest denominator %d' % stats39[case])
X184 = mask('WXY') | 128                                   # R W X Y: rank 3 in AG(3,2), 4 in R8
bad1 = check39(CERT39['F'][str(X184)]['lower'], X184, 1, RM39['NF'][X184], 'F')
bad2 = check39(CERT39['NF'][str(X184)]['upper'], X184, -1, RM39['F'][X184], 'NF')
tampered = json.loads(json.dumps(CERT39['F'][str(X184)]['lower']))
tampered['ineq'][0][1] = str(2 * Fraction(tampered['ineq'][0][1]))
bad3 = check39(tampered, X184, 1, RM39['F'][X184], 'F')
bad4 = any(check39(CERT39['F'][str(X)][k], X, s, RM39['F'][X], 'NF') is False for X in range(1, 256)
           for k, s in (('lower', 1), ('upper', -1)))
report('negative controls: r(RWXY) = %d (AG(3,2)) and %d (R8); the AG(3,2) certificate does not prove the R8 value, nor '
       'conversely; doubling one multiplier breaks a certificate; some Gamma_F certificates are inadmissible for '
       'Gamma_NF' % (RM39['F'][X184], RM39['NF'][X184]),
       RM39['F'][X184] == 3 and RM39['NF'][X184] == 4 and not bad1 and not bad2 and not bad3 and bad4)

# the normal form L'' = {(u, v, w, u+v, u+w, v+w, u+v+w)}, u, v, w in V = Z_q^{2m}, q = 2^a: the map is a 0/1 matrix
# tensored with the identity on Z_q^{2m}, so log2 |pr_T L''| = 2m log2 |image of Z_q^3|.
NF_ROWS = {'A': (1, 0, 0), 'B': (0, 1, 0), 'C': (0, 0, 1), 'W': (1, 1, 0), 'X': (1, 0, 1), 'Y': (0, 1, 1),
           'Z': (1, 1, 1)}


def image_log2(rows, q):
    img = {tuple(sum(r[k] * x[k] for k in range(3)) % q for r in rows) for x in itertools.product(range(q), repeat=3)}
    return len(img).bit_length() - 1


simplex_ok = all(tuple(x % 2 for x in NF_ROWS[P]) == VERT[P] for P in L7)
report('for q = 2 the normal form is the binary simplex code (the cube codewords vanishing at R): the coordinate at '
       'each player is its cube vertex applied to (u, v, w)', simplex_ok)
for a in (1, 2, 3):
    q = 2 ** a
    diff = []
    for T in range(1, 128):
        e = 2 * image_log2([NF_ROWS[P] for i, P in enumerate(L7) if T >> i & 1], q)      # m = 1
        if e != 2 * a * (RM39['F'][T | 128] - 1):
            diff.append((name(T, L7), e, 2 * a * (RM39['F'][T | 128] - 1), 2 * a * (RM39['NF'][T | 128] - 1)))
    want = [] if a == 1 else [('WXY', 6 * a - 2, 4 * a, 6 * a)]
    report('Z_%d, m = 1: log2 |pr_T L\'\'| = 2a (r(RT) - 1) for AG(3,2) on every set T of players%s' %
           (q, '' if a == 1 else ' except T = WXY, where it is 6a - 2 = %d, against %d for AG(3,2) and %d for R8'
            % (6 * a - 2, 4 * a, 6 * a)), diff == want, str(diff) if diff != want else '')
# direct check without the tensor shortcut: q = 4, V = Z_4^2 (m = 1), and the kernel of (u,v,w) -> (u+v, u+w, v+w)
V4 = list(itertools.product(range(4), repeat=2))
add4 = lambda x, y: tuple((s + t) % 4 for s, t in zip(x, y))
pts4 = [(u, v, w, add4(u, v), add4(u, w), add4(v, w), add4(add4(u, v), w)) for u in V4 for v in V4 for w in V4]
direct = all(len({tuple(p[i] for i in range(7) if T >> i & 1) for p in pts4}).bit_length() - 1
             == 2 * image_log2([NF_ROWS[P] for i, P in enumerate(L7) if T >> i & 1], 4) for T in range(1, 128))
ker = [p for p in pts4 if p[3] == p[4] == p[5] == (0, 0)]
report('Z_4, V = Z_4^2 by enumeration of all 4096 elements: the orders of all 127 projections agree with the tensor '
       'computation; the kernel of (u,v,w) -> (u+v, u+w, v+w) is {(x,-x,-x) : 2x = 0}, of order 2^{2m} = 4',
       direct and len(ker) == 4 and all(p[1] == tuple(-s % 4 for s in p[0]) == p[2] and all(2 * s % 4 == 0 for s in p[0])
                                        for p in ker))
# the converse direction: the cube state over Z_d, d odd, is an ideal scheme for Gamma_NF
for d in (9, 15, 21, 25, 27, 45):
    O = {}
    for X in range(1, 256):
        cols = [CUBE_COLS[i] for i in range(8) if X >> i & 1]
        vecs = [[c[k] for c in cols] for k in range(4)]
        o = 1
        for pr, e in factor_small(d).items():
            o *= pr ** zn_log_order(vecs, len(cols), pr ** e, pr)
        O[X] = o
    free = all(O[X] == d ** rQ[X] for X in range(1, 256))
    report('Z_%d: the cube code has |pr_X C| = %d^{r_Q(X)} for every set X of the eight vertices, so the cube state is '
           'the R8 scheme: ideal (one qudit per system, all maximally mixed) for Gamma_NF' % (d, d), free)

# ----------------------------------------------------------------------------------------------------------------------
print('(l) erasures of Steane\'s code (Supplementary Note 10.6) and the average share size (Sec. II D, Introduction)')
GEN = [[VERT[P][k] for P in L7] for k in range(3)]                    # 3 x 7: columns are the cube vertices
HAM = [v for v in itertools.product((0, 1), repeat=7)
       if all(sum(r[i] * v[i] for i in range(7)) % 2 == 0 for r in GEN)]
SIMP = {tuple(sum(c[k] * GEN[k][i] for k in range(3)) % 2 for i in range(7)) for c in itertools.product((0, 1), repeat=3)}
LOGICAL = [v for v in HAM if v not in SIMP]                         # nontrivial logical X (and, the code being CSS with
                                                                     # equal X and Z parts, logical Z) operators
line_masks = [mask(l) for l in LINES7]
size_count, beyond, agree, agree_access = {}, [], True, True
for E in range(128):
    corr = not any(all(v[i] == 0 or E >> i & 1 for i in range(7)) for v in LOGICAL)
    noline = not any(lm & E == lm for lm in line_masks)
    agree = agree and corr == noline
    if E:
        agree_access = agree_access and corr == (not AUTH39['F'][E])
    if corr:
        size_count[popcount(E)] = size_count.get(popcount(E), 0) + 1
        if popcount(E) >= 3:
            beyond.append(E)
complements = sorted(127 ^ lm for lm in line_masks)
report('Steane\'s code: %d Hamming codewords, %d in the simplex code, %d nontrivial logical representatives; the erasure '
       'of E is correctable exactly when E contains no line of the Fano plane, exactly when E is unauthorized for '
       'Gamma_F' % (len(HAM), len(SIMP), len(LOGICAL)), len(HAM) == 16 and len(SIMP) == 8 and agree and agree_access)
report('correctable erasures by size 0..7: %s; the %d of size 3 or 4 (beyond the distance 3) are the 28 triples that are '
       'not lines and the 7 complements of lines'
       % ([size_count.get(k, 0) for k in range(8)], len(beyond)),
       [size_count.get(k, 0) for k in range(8)] == [1, 7, 21, 28, 7, 0, 0, 0] and len(beyond) == 35
       and sorted(E for E in beyond if popcount(E) == 4) == complements
       and all(E not in line_masks for E in beyond if popcount(E) == 3))
on_lines = [sum(1 for lm in line_masks if lm >> i & 1) for i in range(7)]
coef = Fraction(7 * 4, on_lines[0])
report('every point lies on %s lines, so the seven line inequalities of Proposition 1 (i), sum_{P in l} S(P) >= 4 S(R), '
       'add up to sum_P S(P) >= %s S(R) (average share >= %s S(R)); a five-qudit secret needs a share of at least %d '
       'qudits' % (on_lines[0], coef, coef / 7, math.ceil(Fraction(20, 3))),
       on_lines == [3] * 7 and coef == Fraction(28, 3) and coef / 7 == Fraction(4, 3) and math.ceil(Fraction(20, 3)) == 7)

# ----------------------------------------------------------------------------------------------------------------------
print('(m) Lemma 3 (Fano lines in groups of odd order) and the bound 4/3 for every odd local dimension (Sec. II D)')
# The form of Lemma 3 on the variables T_0, ..., T_6 = bits 0..6, built from its definition:
#   g = h(1)+h(2)+h(3) - [I(0;123)+I(0;14)+I(0;25)+I(0;36)] + [I(0;4)+I(0;5)+I(0;6)]
#       + [I(0;234)+I(0;135)+I(0;126)+I(0;456)] >= 0,   I(0;X) = h(0)+h(X)-h(0X).
L7F = {}


def _add(form, X, c):
    if X and c:
        form[X] = form.get(X, 0) + c
        if form[X] == 0:
            del form[X]


def _I0(form, X, c):
    _add(form, 1, c); _add(form, X, c); _add(form, X | 1, -c)


_b = lambda *ii: sum(1 << i for i in ii)
for i in (1, 2, 3):
    _add(L7F, _b(i), 1)
_I0(L7F, _b(1, 2, 3), -1)
for i in (1, 2, 3):
    _I0(L7F, _b(i, i + 3), -1)
    _I0(L7F, _b(i + 3), 1)
for SS in ((2, 3, 4), (1, 3, 5), (1, 2, 6), (4, 5, 6)):
    _I0(L7F, _b(*SS), 1)
FV7 = form_vector(L7F)
# the Fano configuration in the labeling of Lemma 3: lines 014, 025, 036, 234, 135, 126, 456
FCONF = [(1, 1, 1), (1, 0, 0), (0, 1, 0), (0, 0, 1), (0, 1, 1), (1, 0, 1), (1, 1, 0)]
lines_lab = [(0, 1, 4), (0, 2, 5), (0, 3, 6), (2, 3, 4), (1, 3, 5), (1, 2, 6), (4, 5, 6)]
rF7 = rank_function(FCONF, 2)
lines_ok = all(rF7[_b(*L)] == 2 for L in lines_lab) and sum(1 for L in itertools.combinations(range(7), 3)
                                                            if rF7[_b(*L)] == 2) == 7
lhs = rF7[_b(1)] + rF7[_b(2)] + rF7[_b(3)]
i0 = lambda h, X: h[1] + h[X] - h[X | 1]
rhs = (i0(rF7, _b(1, 2, 3)) + i0(rF7, _b(1, 4)) + i0(rF7, _b(2, 5)) + i0(rF7, _b(3, 6)) - i0(rF7, _b(4)) - i0(rF7, _b(5))
       - i0(rF7, _b(6)) - i0(rF7, _b(2, 3, 4)) - i0(rF7, _b(1, 3, 5)) - i0(rF7, _b(1, 2, 6)) - i0(rF7, _b(4, 5, 6)))
report('the form of Lemma 3 has %d terms; the configuration (1,1,1), e_1, e_2, e_3, e_2+e_3, e_1+e_3, e_1+e_2 over '
       'GF(2) has exactly the seven lines 014, 025, 036, 234, 135, 126, 456 (the Fano plane); there the left side is %d '
       'and the right side %d, so the form equals %d' % (len(L7F), lhs, rhs, ev(L7F, rF7)),
       lines_ok and lhs == 3 and rhs == 4 and ev(L7F, rF7) == -1)
vals = {}
for q in (2, 3, 5, 7):
    rk = np.array(rank_function(FCONF, q), dtype=np.int64)
    vals[q] = (int(ev(L7F, [int(x) for x in rk])), int((rk[PMAP] @ FV7).min()))
report('the same vectors over GF(3), GF(5), GF(7) (the non-Fano configuration: 4, 5, 6 not collinear): the form is '
       '>= 0 in all 5040 labelings; over GF(2) its minimum over the labelings is -1',
       vals[2] == (-1, -1) and all(vals[q][0] >= 0 and vals[q][1] >= 0 for q in (3, 5, 7)), str(vals))
tests7, ARR7 = {}, []
for ring in (('field', 3, 3), ('field', 5, 5), ('ring', 9, 3), ('ring', 25, 5), ('ring', 27, 3)):
    mn = None
    for trial in range(40):
        rk = random_rank_function(ring, rng) if trial % 2 == 0 else perturbed_configuration(ring, rng)
        ARR7.append(np.array(rk, dtype=np.int64))
        v = int((ARR7[-1][PMAP] @ FV7).min())
        mn = v if mn is None else min(mn, v)
    tests7['%s%d' % ('GF' if ring[0] == 'field' else 'Z', ring[1])] = mn
z4 = [zn_log_order([v for i in range(7) if m >> i & 1 for v in [FCONF[i]]], 3, 4, 2) if m else 0 for m in range(128)]
z9 = [zn_log_order([v for i in range(7) if m >> i & 1 for v in [FCONF[i]]], 3, 9, 3) if m else 0 for m in range(128)]
report('test (the lemma is proved by hand): the form is >= 0 in all 5040 labelings on 40 random and perturbed '
       'arrangements each over GF(3), GF(5), Z_9, Z_25, Z_27 (groups of odd order); on the configuration above it is '
       '>= 0 over Z_9 and negative over Z_4',
       all(v >= 0 for v in tests7.values()) and ev(L7F, z9) >= 0 and ev(L7F, z4) < 0,
       'minima %s; configuration: Z_9 %d, Z_4 %d (units of log 3, log 2)' % (tests7, ev(L7F, z9), ev(L7F, z4)))
# the incidences used in the proof of Proposition 1 (i), for all 28 pairs (line, point off the line)
inc_ok, npairs = True, 0
for lm in [mask(l) for l in LINES7]:
    for qi in range(7):
        if lm >> qi & 1:
            continue
        npairs += 1
        xs = [i for i in range(7) if lm >> i & 1]
        third = {}
        for x in xs:
            th = [l for l in LINES7 if mask(l) >> x & 1 and mask(l) >> qi & 1]
            third[x] = [i for i in range(7) if mask(th[0]) >> i & 1 and i not in (x, qi)][0]
        for perm in itertools.permutations(xs):
            lab = list(perm) + [third[x] for x in perm]          # x_1, x_2, x_3, x_4, x_5, x_6
            setq = lambda ii: sum(1 << lab[i - 1] for i in ii) | (1 << qi)
            pos = [setq(X) for X in ((1, 2, 3), (1, 4), (2, 5), (3, 6))]
            neg = [setq(X) for X in ((4,), (5,), (6,), (2, 3, 4), (1, 3, 5), (1, 2, 6), (4, 5, 6))]
            inc_ok = inc_ok and all(AUTH39['F'][T] for T in pos) and not any(AUTH39['F'][T] for T in neg) \
                and not AUTH39['F'][1 << qi] and sorted(lab[3:]) == sorted(set(range(7)) - set(xs) - {qi})
report('for all %d pairs (line l, player q not on l) and every order of the points x_1, x_2, x_3 of l, with x_{i+3} the '
       'third point of the line through q and x_i: X u {q} is authorized for Gamma_F for X = 123, 14, 25, 36 and '
       'unauthorized for X = 4, 5, 6, 234, 135, 126, 456; q is unauthorized; x_4, x_5, x_6 are the points off l u {q}'
       % npairs, npairs == 28 and inc_ok)
# the stronger form proved in Supplementary Note 10.2:
#   I(0;123) <= h(1|04) + h(2|05) + h(3|06) + I(0;234) + I(0;135) + I(0;126) + I(0;456);
# the form of Lemma 3 equals this one plus I(1;4) + I(2;5) + I(3;6) (an identity of linear forms)
STRONG = {}
_I0(STRONG, _b(1, 2, 3), -1)
for i in (1, 2, 3):
    _add(STRONG, _b(0, i, i + 3), 1); _add(STRONG, _b(0, i + 3), -1)           # + h(i | 0, i+3)
for SS in ((2, 3, 4), (1, 3, 5), (1, 2, 6), (4, 5, 6)):
    _I0(STRONG, _b(*SS), 1)
SUM_I = {}
for i in (1, 2, 3):
    _add(SUM_I, _b(i), 1); _add(SUM_I, _b(i + 3), 1); _add(SUM_I, _b(i, i + 3), -1)
comb = dict(STRONG)
for m, c in SUM_I.items():
    _add(comb, m, c)
report('the stronger form of Supplementary Note 10.2 (%d terms) plus I(1;4) + I(2;5) + I(3;6) is the form of Lemma 3, '
       'term by term; on the Fano configuration its left side is 1 and its right side 0 (value -1); on the odd-'
       'characteristic configuration and on the random and perturbed arrangements above it is >= 0'
       % len(STRONG),
       comb == L7F and ev(STRONG, rF7) == -1 and all(int((np.array(rank_function(FCONF, q), dtype=np.int64)[PMAP]
                                                          @ form_vector(STRONG)).min()) >= 0 for q in (3, 5, 7))
       and all(int((a[PMAP] @ form_vector(STRONG)).min()) >= 0 for a in ARR7))
# the constant 11/2 for approximate schemes: on the contraction by q the right side of Lemma 3 is
# sum over the four positive sets of I(R;X'q) - sum over the seven negative sets + 3 I(R;q)
s_, e_ = Fraction(1), Fraction(1)
worst_rhs = 4 * (2 * s_ - e_) - 7 * e_ + 3 * 0
report('approximate schemes: the right side of Lemma 3 on the contraction by q is the sum of I(R;X\'q) over the four '
       'positive sets minus the sum over the seven negative sets plus 3 I(R;q); with the first >= 2S(R) - eps, the '
       'second <= eps and I(R;q) >= 0 it is >= 8 S(R) - 11 eps, so sum_{P in l} S(P) >= 4 S(R) - (11/2) eps',
       worst_rhs == 8 * s_ - 11 * e_ and (8 * s_ - 11 * e_) / 2 == 4 * s_ - Fraction(11, 2) * e_
       and sum(1 for X in ((1, 2, 3), (1, 4), (2, 5), (3, 6))) == 4 and 7 - 4 == 3)
# classical abelian schemes for the Fano port F of Jafari-Khazaei: the labels of Lemma 3 for the triple 123
F6set = closure6(JK_F)
mk6 = lambda ii: sum(1 << (i - 1) for i in ii)
labels_ok = all(mk6(X) in F6set for X in ((1, 2, 3), (1, 4), (2, 5), (3, 6))) and \
    not any(mk6(X) in F6set for X in ((4,), (5,), (6,), (2, 3, 4), (1, 3, 5), (1, 2, 6), (4, 5, 6)))
autF = [pm for pm in itertools.permutations(range(6))
        if {sum(1 << pm[i] for i in range(6) if A >> i & 1) for A in F6set} == F6set]
triples = [mk6(t) for t in ((1, 2, 3), (1, 5, 6), (2, 4, 6), (3, 4, 5))]
orbit = {sum(1 << pm[i] for i in range(6) if triples[0] >> i & 1) for pm in autF}
cover = [sum(1 for t in triples if t >> i & 1) for i in range(6)]
report('the Fano port F = 14+25+36+123+156+246+345: for the triple 123 the sets 123, 14, 25, 36 are authorized and 4, 5, 6, '
       '234, 135, 126, 456 unauthorized; its %d automorphisms act transitively on the four triples, and every player '
       'lies in two of them, so the four triple inequalities give an average share >= 4/3' % len(autF),
       labels_ok and orbit == set(triples) and cover == [2] * 6)

# Lemma 3 is not a consequence of the known odd-characteristic seven-variable inequalities (Supplementary Note 10.4):
# an integer polymatroid W on seven elements (bit i <-> 'ABCWXYZ'[i], W(N) = 120) satisfies the elemental Shannon
# inequalities, all Ingleton instances I(A;B|K) <= I(A;B|CK) + I(A;B|DK) + I(C;D|K) on disjoint sets and the forms
# (65), (6.17) corrected, (b) (times 3) and LB of Supplementary Table S4 in all 5040 labelings, while the form of
# Lemma 3, with T_0, ..., T_6 = bits 0..6, takes the value -37 on it.  The coefficients below are transcribed from that
# table (sets of the seven points A, B, C, W, X, Y, Z).  W was found by a linear program in floating point and rounded;
# here it is checked exactly.
W120 = [
    0, 39, 38, 77, 38, 77, 76, 115, 38, 77, 76, 115, 76, 115, 114, 115,
    44, 83, 82, 83, 81, 120, 119, 120, 81, 120, 119, 120, 81, 120, 119, 120,
    44, 83, 81, 120, 82, 83, 119, 120, 81, 120, 81, 120, 119, 120, 119, 120,
    82, 120, 119, 120, 119, 120, 119, 120, 119, 120, 119, 120, 119, 120, 119, 120,
    44, 83, 81, 120, 81, 120, 81, 120, 82, 83, 119, 120, 119, 120, 119, 120,
    82, 120, 119, 120, 119, 120, 119, 120, 119, 120, 119, 120, 119, 120, 119, 120,
    82, 120, 119, 120, 119, 120, 119, 120, 119, 120, 119, 120, 119, 120, 119, 120,
    82, 120, 119, 120, 119, 120, 119, 120, 119, 120, 119, 120, 119, 120, 119, 120,
]
ODD4 = {
    '(65)': {'A': 3, 'B': 4, 'C': 3, 'W': 1, 'X': 1, 'Y': 1, 'Z': 1, 'AB': -5, 'AX': -2, 'BC': -5, 'CW': -3, 'WY': -3,
             'XZ': -1, 'YZ': -2, 'ABC': -5, 'ABW': 5, 'ACX': 2, 'AYZ': 2, 'BCY': 5, 'BXZ': 1, 'CWZ': 3, 'WXY': 3},
    '(6.17) corrected': {'A': 3, 'B': 3, 'C': 9, 'Z': 6, 'AB': 6, 'AY': -3, 'BX': -3, 'CW': -3, 'ABC': -12, 'ABW': 3,
                         'ABZ': -3, 'ACX': 3, 'ACZ': -3, 'AYZ': 3, 'BCY': 3, 'BCZ': -3, 'BXZ': 3, 'CWZ': 3, 'WXY': 1,
                         'ABCZ': 3},
    '(b) times 3': {'A': 6, 'B': 3, 'C': 6, 'Z': 6, 'AB': 3, 'AY': -3, 'BC': 3, 'BX': -3, 'CW': -3, 'ABC': -12,
                    'ABW': 3, 'ABZ': -3, 'ACX': 3, 'ACZ': -3, 'AYZ': 3, 'BCY': 3, 'BCZ': -3, 'BXZ': 3, 'CWZ': 3,
                    'WXY': 1, 'ABCZ': 3},
    'LB': {'A': 9, 'B': 9, 'C': 9, 'Z': 6, 'AY': -3, 'BX': -3, 'CW': -3, 'ABC': -12, 'ABW': 3, 'ABZ': -3, 'ACX': 3,
           'ACZ': -3, 'AYZ': 3, 'BCY': 3, 'BCZ': -3, 'BXZ': 3, 'CWZ': 3, 'WXY': 1, 'ABCZ': 3},
}
W7 = np.array(W120, dtype=np.int64)
_LET = 'ABCWXYZ'
odd_vecs = {}
for nm, f in ODD4.items():
    v = np.zeros(128, dtype=np.int64)
    for S_, c_ in f.items():
        v[sum(1 << _LET.index(ch) for ch in S_)] += c_
    odd_vecs[nm] = v
# each transcribed form is -1 on the Fano configuration (A..Z = 100, 010, 001, 110, 101, 011, 111 over GF(2)) and >= 0
# on the same vectors over GF(3)
VEC_AZ = [(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0), (1, 0, 1), (0, 1, 1), (1, 1, 1)]
rF_AZ = np.array(rank_function(VEC_AZ, 2), dtype=np.int64)
rNF_AZ = np.array(rank_function(VEC_AZ, 3), dtype=np.int64)
transc_ok = all(int(rF_AZ @ v) == -1 and int(rNF_AZ @ v) >= 0 for v in odd_vecs.values())
def _I(h_, a_, b_, k_=0):
    return int(h_[a_ | k_] + h_[b_ | k_] - h_[a_ | b_ | k_] - h_[k_])
shannon_ok = W7[0] == 0 and all(W7[127] - W7[127 ^ (1 << i)] >= 0 for i in range(7))
for i_, j_ in itertools.combinations(range(7), 2):
    rest_ = [k for k in range(7) if k not in (i_, j_)]
    for s_ in range(32):
        K_ = sum(1 << rest_[t] for t in range(5) if s_ >> t & 1)
        shannon_ok = shannon_ok and _I(W7, 1 << i_, 1 << j_, K_) >= 0
ingleton_ok, n_ing = True, 0
for assign_ in itertools.product(range(6), repeat=7):        # each element in none, A, B, C, D or K
    SS_ = [0] * 6
    for i_, g_ in enumerate(assign_):
        SS_[g_] |= 1 << i_
    a_, b_, c_, d_, K_ = SS_[1:]
    if a_ and b_ and c_ and d_:
        n_ing += 1
        ingleton_ok = ingleton_ok and (_I(W7, a_, b_, c_ | K_) + _I(W7, a_, b_, d_ | K_) + _I(W7, c_, d_, K_)
                                       - _I(W7, a_, b_, K_)) >= 0
odd_min = {nm: int((W7[PMAP] @ v).min()) for nm, v in odd_vecs.items()}
l7_val = int(ev(L7F, [int(x) for x in W7]))
report('Lemma 3 does not follow from the Shannon inequalities, the Ingleton inequalities (all %d instances with disjoint '
       'nonempty A, B, C, D and any disjoint K) and the forms (65), (6.17) corrected, (b) and LB in all 5040 labelings: '
       'the transcribed forms are -1 on the Fano and >= 0 on the non-Fano configuration; the integer polymatroid W '
       '(W(N) = 120) satisfies all of them (minima %s), and the form of Lemma 3 is %d < 0 on it' % (n_ing, odd_min, l7_val),
       transc_ok and shannon_ok and ingleton_ok and n_ing == 57120 and min(odd_min.values()) >= 0 and l7_val == -37)

# The four forms allow the information ratio 23/22 for the Fano port F of Jafari-Khazaei (Supplementary Note 10.4): an integer polymatroid
# V on the secret (bit 0) and the players 1..6 (bits 1..6, labels of F) with V(secret) = 22 and every share V(i) <= 23
# satisfies the elemental Shannon inequalities, all Ingleton instances, the four forms in all 5040 labelings and the
# access equations of F (V(Xs) = V(X) for qualified X, V(Xs) = V(X) + V(s) otherwise); Lemma 3 rules it out.
V22 = [
    0, 22, 23, 45, 23, 45, 46, 68, 23, 45, 46, 68, 46, 68, 69, 69,
    23, 45, 46, 46, 46, 68, 69, 69, 46, 68, 69, 69, 48, 70, 70, 70,
    23, 45, 46, 68, 46, 46, 69, 69, 46, 68, 48, 70, 69, 69, 70, 70,
    46, 68, 69, 69, 69, 69, 70, 70, 69, 69, 70, 70, 70, 70, 70, 70,
    23, 45, 46, 68, 46, 68, 48, 70, 46, 46, 69, 69, 69, 69, 70, 70,
    46, 68, 69, 69, 69, 69, 70, 70, 69, 69, 70, 70, 70, 70, 70, 70,
    46, 68, 69, 69, 69, 69, 70, 70, 69, 69, 70, 70, 70, 70, 70, 70,
    48, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70,
]
V7 = np.array(V22, dtype=np.int64)
F_JK = [(1, 4), (2, 5), (3, 6), (1, 2, 3), (1, 5, 6), (2, 4, 6), (3, 4, 5)]
FJm = [sum(1 << i for i in t) for t in F_JK]
qual_ = lambda X: any((X & mm) == mm for mm in FJm)
acc_ok = all((V7[X | 1] == V7[X]) if qual_(X) else (V7[X | 1] == V7[X] + V7[1]) for X in range(2, 128, 2))
sh_ok = V7[0] == 0 and all(V7[127] - V7[127 ^ (1 << i)] >= 0 for i in range(7))
for i_, j_ in itertools.combinations(range(7), 2):
    rest_ = [k for k in range(7) if k not in (i_, j_)]
    for s_ in range(32):
        K_ = sum(1 << rest_[t] for t in range(5) if s_ >> t & 1)
        sh_ok = sh_ok and _I(V7, 1 << i_, 1 << j_, K_) >= 0
ing_ok = True
for assign_ in itertools.product(range(6), repeat=7):
    SS_ = [0] * 6
    for i_, g_ in enumerate(assign_):
        SS_[g_] |= 1 << i_
    a_, b_, c_, d_, K_ = SS_[1:]
    if a_ and b_ and c_ and d_:
        ing_ok = ing_ok and (_I(V7, a_, b_, c_ | K_) + _I(V7, a_, b_, d_ | K_) + _I(V7, c_, d_, K_) - _I(V7, a_, b_, K_)) >= 0
odd_min22 = {nm: int((V7[PMAP] @ v).min()) for nm, v in odd_vecs.items()}
shares_ = [int(V7[1 << i]) for i in range(1, 7)]
triples_ok = all(sum(int(V7[1 << i]) for i in t) < 4 * int(V7[1]) for t in F_JK if len(t) == 3)
report('the four forms allow the information ratio 23/22 for the Fano port F: an integer polymatroid with secret 22 and shares %s '
       'satisfies the Shannon inequalities, all Ingleton instances, the four forms in all 5040 labelings (minima %s) and '
       'the access equations of F; every minimal triple has share sum < 4 x secret, as Lemma 3 forbids'
       % (shares_, odd_min22), acc_ok and sh_ok and ing_ok and min(odd_min22.values()) >= 0 and V7[1] == 22
       and max(shares_) == 23 and triples_ok)

# ======================================================= (n) products over prime powers: d = 2r with r odd
print('(n) products over prime powers: ratio 7/5 for Gamma_F and Gamma_NF with Z_d qudits, d = 2r, r odd (Prop. 1 (iii))')


def order_gen(functionals, d):
    """order of the projection: the subgroup of Z_d^k (k = number of functionals) generated by the images of the
    parameter vectors, computed for each prime power p^a of d after reduction modulo p^a"""
    if not functionals:
        return 1
    npar = len(functionals[0])
    o = 1
    for pr, a in factor_small(d).items():
        N = pr ** a
        vecs = [[f[i] % N for f in functionals] for i in range(npar)]
        vecs = [v for v in vecs if any(v)]
        o *= pr ** zn_log_order(vecs, len(functionals), N, pr) if vecs else 1
    return o


def product_code(cols2, colsq, d, q):
    """the code over Z_d (d = 2r, r odd; the argument q is r) whose reduction modulo 2 is the binary code cols2 and whose
    reduction modulo r is the code colsq over Z_r, on disjoint parameters; e2 = r and eq = 1 - r are the idempotents"""
    e2, eq = q % d, (1 - q) % d
    out = {}
    for w in cols2:
        assert len(cols2[w]) == len(colsq[w])
        out[w] = [[e2 * x % d for x in f2] + [eq * x % d for x in fq] for f2, fq in zip(cols2[w], colsq[w])]
    return out


def pad(cols, k):
    """add coordinates that vanish on the code (qudits in the state |0>) until every player holds k"""
    npar = len(cols['R'][0])
    return {w: fs + ([[0] * npar for _ in range(k - len(fs))] if w != 'R' else []) for w, fs in cols.items()}


def cube_copies(c):
    """c copies of the cube code (affine functions s + a.x on {0,1}^3, the origin R holds s); parameters (s, a) per copy"""
    npar = 4 * c
    cols = {'R': [[1 if k == 4 * j else 0 for k in range(npar)] for j in range(c)]}
    for P in L7:
        cols[P] = []
        for j in range(c):
            f = [0] * npar
            f[4 * j] = 1
            for i in range(3):
                f[4 * j + 1 + i] = VERT[P][i]
            cols[P].append(f)
    return cols


def product_report(cols, d, auth, sec, sh):
    owners = ['R'] + list(L7)
    O = {m8: order_gen([f for i, w in enumerate(owners) if m8 >> i & 1 for f in cols[w]], d) for m8 in range(256)}
    oR = O[1]
    acc = all((O[(T << 1) | 1] == O[T << 1]) if auth(T) else (O[(T << 1) | 1] == O[T << 1] * oR) for T in range(128))
    full = 255
    E = {m8: Fraction(O[m8] * O[full ^ m8], O[full]) for m8 in range(256)}
    css = E[1] == oR and all((E[1] * E[T << 1] / E[(T << 1) | 1] == E[1] ** 2) if auth(T)
                             else (E[1] * E[T << 1] / E[(T << 1) | 1] == 1) for T in range(128))
    sizes = oR == d ** sec and all(len(cols[P]) == sh for P in L7)      # shares of sh qudits; padded ones are not full rank
    return acc, css, sizes


AUTH_F = lambda T: any(T & mm == mm for mm in MIN_F)
AUTH_NF = gamma_swap(['WXY'])
NF_BIN = {w: to_binary(fs) for w, fs in nf_code(B4, PIECES_NF).items()}       # the qubit scheme of (i): 10 and 14
for q_ in (3, 5, 9, 15):
    d_ = 2 * q_
    accF, cssF, szF = product_report(product_code(pad(cube_copies(5), 7), ramp_code(B_RAMP), d_, q_), d_, AUTH_F, 5, 7)
    report('Z_%d: five qubit cube states (two |0> added per player) times the Z_%d scheme of (h) is a scheme for Gamma_F '
           'on all 128 player sets; secret Z_d^5 (full rank), each player 7 coordinates (ratio 7/5); its CSS state has '
           'I(R;T) = 2 S(R) on Gamma_F and 0 elsewhere' % (d_, q_), accF and cssF and szF)
    accN, cssN, szN = product_report(product_code(NF_BIN, pad(cube_copies(10), 14), d_, q_), d_, AUTH_NF, 10, 14)
    report('Z_%d: the qubit scheme of (i) times ten Z_%d cube states (four |0> added per player) is a scheme for '
           'Gamma_NF on all 128 player sets; secret Z_d^10 (full rank), each player 14 coordinates (ratio 7/5); its CSS '
           'state has I(R;T) = 2 S(R) on Gamma_NF and 0 elsewhere' % (d_, q_), accN and cssN and szN)
badF, cssF_bad, _ = product_report(product_code(pad(cube_copies(5), 7), pad(cube_copies(5), 7), 6, 3), 6, AUTH_F, 5, 7)
badN, cssN_bad, _ = product_report(product_code(pad(cube_copies(10), 14), pad(cube_copies(10), 14), 6, 3), 6, AUTH_NF,
                                   10, 14)
report('negative controls over Z_6: qubit cube states times qutrit cube states is not a scheme for Gamma_F (the qutrit '
       'cube state is ideal for Gamma_NF, not Gamma_F) and not one for Gamma_NF (the qubit cube state is ideal for '
       'Gamma_F, not Gamma_NF); both as classical schemes and through the entropies of the CSS states',
       not badF and not badN and not cssF_bad and not cssN_bad)


def pad_ref(cols, k):
    """add coordinates at R that vanish on the code (reference qudits in a fixed state) until R holds k"""
    npar = len(cols['R'][0])
    out = dict(cols)
    out['R'] = cols['R'] + [[0] * npar for _ in range(k - len(cols['R']))]
    return out


# For even d the counts in qudits need secrets of whole qudits (Supplementary Note 10.3, the ranges of the optimal ratios): over
# Z_6, seven qubit cube states times the Z_3 scheme of (h), with two qutrits in a fixed state added to R, is a scheme
# for Gamma_F whose secret has dimension 2^7 3^5 and whose shares have 7 qudits each.
mixed6 = product_code(cube_copies(7), pad_ref(ramp_code(B_RAMP), 7), 6, 3)
accM, cssM, _ = product_report(mixed6, 6, AUTH_F, 7, 7)
oRM = order_gen(mixed6['R'], 6)
ratioM = 7 * math.log(6) / math.log(oRM)
report('Z_6, a secret that is not a whole number of qudits: seven qubit cube states times the Z_3 scheme of (h), with '
       'two qutrits in a fixed state added to R, is a scheme for Gamma_F on all 128 player sets (classically and through '
       'the entropies of its CSS state); its secret has order 2^7 3^5 = %d and every share 7 qudits, so its information '
       'ratio is 7 log 6/(7 log 2 + 5 log 3) = %.4f < 4/3 (the example of Supplementary Note 10.3)' % (oRM, ratioM),
       accM and cssM and oRM == 2 ** 7 * 3 ** 5 and all(len(mixed6[P]) == 7 for P in L7)
       and abs(ratioM - 7 * math.log(6) / (7 * math.log(2) + 5 * math.log(3))) < 1e-12 and ratioM < 4 / 3)

# ----------------------------------------------------------------------------------------------------------------------
print('(o) codes for large secrets: the numbers of qudits (Sec. II D, Supplementary Note 10.7)')


def least_m(d, k):
    """least m with 2^m >= d^k, that is, ceil(k log2 d), in exact integer arithmetic"""
    m = 0
    while 2 ** m < d ** k:
        m += 1
    return m


def least_n(d, m):
    """least n with d^n >= 2^m, that is, ceil(m log_d 2), in exact integer arithmetic"""
    n = 0
    while d ** n < 2 ** m:
        n += 1
    return n


rows = []
for k_, m_exp, n_exp, tot_exp, stab_exp in ((4, 7, 5, 35, 38), (5, 8, 6, 42, 47)):
    m_ = least_m(3, k_)
    n_ = least_n(3, m_)
    stab = -((-28 * k_) // 3)          # ceil(28k/3): each of the seven lines holds at least 4k qutrits, 3 lines per player
    rows.append(m_ == m_exp and n_ == n_exp and 7 * n_ == tot_exp and stab == stab_exp and 3 ** n_ >= 2 ** m_ >= 3 ** k_)
report('d = 3, k = 4 and 5: (m, n) = (7, 5) and (8, 6), so the codes use 7n = 35 and 42 qutrits, against ceil(28k/3) = '
       '38 and 47 for every qutrit stabilizer scheme for Gamma_F with a secret of k qutrits; d^n >= 2^m >= d^k',
       all(rows))
lim = []
for d_ in (3, 5, 6, 7):
    k_ = 1000
    n_ = least_n(d_, least_m(d_, k_))
    lim.append(k_ <= n_ < 1.01 * k_)
report('n/k tends to 1: for d = 3, 5, 6, 7 and k = 1000, k <= n < 1.01 k', all(lim))

# ----------------------------------------------------------------------------------------------------------------------
print('(p) ports of the generalized cube codes contracted at C (Supplementary Note 10.1)')


def in_span(vecs, target, p):
    """whether target lies in the GF(p)-span of vecs"""
    return rank_mod_p(vecs + [target], p) == rank_mod_p(vecs, p)


rows = []
for t_ in (2, 3, 4):
    M_ = t_ + 1
    e_ = [[1 if j == i else 0 for j in range(M_)] for i in range(M_)]
    c_ = [1] * M_
    players = [('A', i) for i in range(M_)] + [('B', i) for i in range(M_)]
    for p_ in (2, 3, 5, 7):
        colA = [[1] + e_[i] for i in range(M_)]
        colB = [[1] + [c_[j] - e_[i][j] for j in range(M_)] for i in range(M_)]
        colC = [1] + c_
        colO = [1] + [0] * M_
        same = True
        for r_ in range(len(players) + 1):
            for T in itertools.combinations(players, r_):
                ours = in_span([colA[i] if k == 'A' else colB[i] for k, i in T] + [colC], colO, p_)
                conf = [[c_[j] - e_[i][j] for j in range(M_)] if k == 'A' else e_[i] for k, i in T]
                theirs = in_span(conf, c_, p_) if conf else False
                same = same and ours == theirs
        rows.append(same)
report('for t = 2, 3, 4 and p = 2, 3, 5, 7: the port at O of the generalized cube code, contracted at C, equals the port at '
       'c of the configuration J - I with A_i and B_i exchanged (all 2^(2t+2) sets of the other players)', all(rows))

# ----------------------------------------------------------------------------------------------------------------------
print('(q) planes and good directions of AG(3,2) and R8 (proof of Theorem 5, Supplementary Note 10.7)')
COORD = [VERT[c] for c in L7] + [(0, 0, 0)]                           # bit i -> point of {0,1}^3; R = origin = bit 7
AT = {COORD[i]: i for i in range(8)}
DIRS = [v for v in itertools.product((0, 1), repeat=3) if any(v)]


def vxor(u, v):
    return tuple((a + b) % 2 for a, b in zip(u, v))


def dname(v):
    return ''.join(str(x) for x in v)


def coplanar_real(X):
    pts = [COORD[i] for i in range(8) if X >> i & 1]
    return rank_Q([[p[k] - pts[0][k] for k in range(3)] for p in pts[1:]]) <= 2


def is_face(X):
    pts = [COORD[i] for i in range(8) if X >> i & 1]
    return any(len({p[k] for p in pts}) == 1 for k in range(3))


def dpairs(delta):
    out = []
    for i in range(8):
        j = AT[vxor(COORD[i], delta)]
        if i < j:
            out.append(1 << i | 1 << j)
    return out


def splits_along(X, delta):
    return all(AT[vxor(COORD[i], delta)] in [k for k in range(8) if X >> k & 1] for i in range(8) if X >> i & 1)


TETRA = {mask('WXY') | RB, mask('ABCZ')}
PL = {}
for nm_, r_ in (('AG(3,2)', rA), ('R8', rR)):
    PL[nm_] = {X for X in range(256) if popcount(X) == 4 and r_[X] == 3}
report('planes (sets of four points of rank 3): AG(3,2) has 14, the 4-sets with zero sum; R8 has 12, the 4-sets of '
       'coplanar vertices of the real cube, namely those of AG(3,2) other than the tetrahedra RWXY and ABCZ: 6 faces '
       'and 6 diagonal rectangles',
       PL['AG(3,2)'] == planes and len(planes) == 14
       and PL['R8'] == {X for X in range(256) if popcount(X) == 4 and coplanar_real(X)} and len(PL['R8']) == 12
       and PL['AG(3,2)'] - PL['R8'] == TETRA and sum(1 for X in PL['R8'] if is_face(X)) == 6)
report('(A): in both matroids every set of at most three points is independent',
       all(r_[X] == popcount(X) for r_ in (rA, rR) for X in range(256) if popcount(X) <= 3))
GOOD = {}
for nm_, r_ in (('AG(3,2)', rA), ('R8', rR)):
    GOOD[nm_] = [d for d in DIRS if all(r_[P | Q] == 3 for P, Q in itertools.combinations(dpairs(d), 2))]
bad_unions = {dname(d): {P | Q for P, Q in itertools.combinations(dpairs(d), 2) if rR[P | Q] != 3}
              for d in DIRS if d not in GOOD['R8']}
report('good directions (any two pairs along them form a plane): all seven for AG(3,2); for R8 exactly 100, 010, 001 '
       '(edges) and 111 (body diagonals); along each face diagonal 110, 101, 011, the two unions of two pairs that are '
       'not planes of R8 are the tetrahedra RWXY and ABCZ',
       len(GOOD['AG(3,2)']) == 7 and sorted(GOOD['R8']) == sorted([(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 1)])
       and sorted(bad_unions) == ['011', '101', '110'] and all(v == TETRA for v in bad_unions.values()),
       'R8: ' + ' '.join(dname(d) for d in GOOD['R8']))
g1, g2, g3, g4 = [], [], [], []
for nm_, r_ in (('AG(3,2)', rA), ('R8', rR)):
    gd = GOOD[nm_]
    along = [sum(1 for d in gd if splits_along(X, d)) for X in PL[nm_]]
    g1.append(set(along) == ({3} if nm_ == 'AG(3,2)' else {2}))
    g2.append(all(r_[sum(1 << AT[p] for p in (v, vxor(v, d), vxor(v, c), vxor(vxor(v, d), c)))] == 3
                  for d, c in itertools.permutations(gd, 2) for v in COORD))
    g3.append(rank_mod_p([list(d) for d in gd], 2) == 3)
    g4.append(all(any(not any(P & X == P for P in dpairs(d)) for X in PL[nm_]) for d in gd))
report('(G1) every plane is split into pairs along good directions in exactly %s ways (AG(3,2), R8); (G2) any two good '
       'directions and any point span a plane; (G3) the good directions span GF(2)^3; (G4) for every good direction '
       'some plane contains none of its pairs' % '3 and 2', all(g1) and all(g2) and all(g3) and all(g4),
       'G1-G4 for AG(3,2): %s, for R8: %s' % ([g1[0], g2[0], g3[0], g4[0]], [g1[1], g2[1], g3[1], g4[1]]))
thr = {}
for nm_, r_, mins_ in (('AG(3,2)', rA, MIN_F), ('R8', rR, MIN_NF)):
    circ = circuits(r_, 8)
    throughR = sorted((Cc ^ RB for Cc in circ if Cc & RB), key=lambda m: (popcount(m), m))
    union = 0
    for Cc in circ:
        if Cc & RB:
            union |= Cc
    thr[nm_] = (sorted(throughR) == sorted(mins_), union == FULL8, len(throughR))
report('the circuits through R, without R, are the minimal authorized sets (7 for Gamma_F, 11 for Gamma_NF), and every '
       'element lies in one of them, so both matroids are connected (as used with Lehman\'s theorem)',
       all(a and b for a, b, _ in thr.values()) and thr['AG(3,2)'][2] == 7 and thr['R8'][2] == 11)
LEM = [('RABC', 4, 4), ('RABW', 3, 3), ('RACX', 3, 3), ('RBCY', 3, 3), ('RAB', 3, 3), ('RAC', 3, 3), ('RBC', 3, 3),
       ('RAYZ', 3, 3), ('RBXZ', 3, 3), ('RCWZ', 3, 3), ('RAY', 3, 3), ('RBX', 3, 3), ('RCW', 3, 3), ('RWXY', 3, 4)]
for I_ in 'ABC':
    for P_ in 'WXYZ':
        if (I_, P_) not in (('C', 'W'), ('B', 'X'), ('A', 'Y')):
            LEM.append(('R' + ''.join(sorted(set('ABC') - {I_})) + P_, 4, 4))


def rmask(s):
    return sum(RB if ch == 'R' else 1 << L7.index(ch) for ch in s)


report('Supplementary Lemma 10.4: the ranks of the %d sets RT used in its proof are as stated in both matroids (the nine '
       'sets of step (3) have rank 4), and only RWXY differs: 3 in AG(3,2), 4 in R8' % len(LEM),
       len(LEM) == 23 and all(rA[rmask(s)] == a and rR[rmask(s)] == b for s, a, b in LEM))

# ----------------------------------------------------------------------------------------------------------------------
print('(r) erasures: Steane\'s code, the color code and all qubit [[7,1,3]] stabilizer codes (Supplementary Note 10.6)')
FATAL = {}
for q in (2, 3):
    cnt = [0] * 8
    for E in range(1, 128):
        if E in AUTH[q]:
            cnt[popcount(E)] += 1
    FATAL[q] = cnt
selfdual = all((E in AUTH[q]) != ((127 ^ E) in AUTH[q]) for q in (2, 3) for E in range(1, 127))
report('uncorrectable erasures, i.e. the authorized sets (a set is authorized exactly when its complement is not), by '
       'size 0..7: Steane\'s code %s, the color code of odd local dimension %s' % (FATAL[2], FATAL[3]),
       FATAL[2] == [0, 0, 0, 7, 28, 21, 7, 1] and FATAL[3] == [0, 0, 0, 6, 29, 21, 7, 1] and selfdual)


def pmul(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out


def ppow(a, k):
    out = [1]
    for _ in range(k):
        out = pmul(out, a)
    return out


def padd(a, b, s=1):
    n_ = max(len(a), len(b))
    a, b = a + [0] * (n_ - len(a)), b + [0] * (n_ - len(b))
    out = [x + s * y for x, y in zip(a, b)]
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


EPS, ONE_M = [0, 1], [1, -1]                                          # eps and 1 - eps, coefficients by degree


def p_fail(cnt):
    tot = [0]
    for k in range(8):
        tot = padd(tot, [cnt[k] * c for c in pmul(ppow(EPS, k), ppow(ONE_M, 7 - k))])
    return tot


diffp = padd(p_fail(FATAL[3]), p_fail(FATAL[2]), -1)
target = [-c for c in pmul(pmul(ppow(EPS, 3), ppow(ONE_M, 3)), [1, -2])]
report('failure probabilities on the erasure channel: P_color(eps) - P_Steane(eps) = -eps^3 (1-eps)^3 (1-2eps) as '
       'polynomials (exact integer coefficients)', diffp == padd(target, [0]), 'coefficients %s' % diffp)

EDGES7 = list(itertools.combinations(range(7), 2))                     # the 21 possible edges among the seven sites
EIDX = {e: k for k, e in enumerate(EDGES7)}
PE = np.array([[EIDX[tuple(sorted((s[a], s[b])))] for a, b in EDGES7] for s in itertools.permutations(range(7))],
              dtype=np.int64)                                         # 5040 x 21: image of each edge under each permutation


def n_cycles(perm):
    seen_, c = [False] * len(perm), 0
    for i in range(len(perm)):
        if not seen_[i]:
            c += 1
            j = i
            while not seen_[j]:
                seen_[j] = True
                j = perm[j]
    return c


burnside = sum(2 ** n_cycles(row) for row in PE.tolist())
WTS = np.left_shift(np.int64(1), PE)
seen_g = np.zeros(1 << 21, dtype=bool)
REPS, m_ = [], 0
while m_ < (1 << 21):
    m_ += int(np.argmin(seen_g[m_:]))
    if seen_g[m_]:
        break
    REPS.append(m_)
    seen_g[WTS @ (np.right_shift(np.int64(m_), np.arange(21, dtype=np.int64)) & 1)] = True
    m_ += 1
report('graphs on seven labeled vertices up to isomorphism: %d orbit representatives (every graph lies in the orbit of '
       'exactly one); Burnside\'s lemma gives %d orbits' % (len(REPS), burnside // 5040),
       len(REPS) == 1044 and burnside == 1044 * 5040 and bool(seen_g.all()))


def rank2(rows):
    rows, r = [x for x in rows if x], 0
    while rows:
        p = rows.pop()
        r += 1
        lb = p & -p
        rows = [x ^ p if x & lb else x for x in rows]
        rows = [x for x in rows if x]
    return r


def s_graph(adj, X):
    """entropy (bits) of the subsystem X of the graph state with adjacency masks adj: the GF(2) rank of the
    adjacency submatrix between X and its complement"""
    comp = FULL8 & ~X
    return rank2([adj[i] & comp for i in range(8) if X >> i & 1])


def graph_adj(g, nbr):
    adj = [0] * 8
    for k, (a, b) in enumerate(EDGES7):
        if g >> k & 1:
            adj[a] |= 1 << b
            adj[b] |= 1 << a
    adj[7] = nbr
    for i in range(7):
        if nbr >> i & 1:
            adj[i] |= RB
    return adj


def graph_state_entropies(adj):
    """entropies (bits) of all subsystems of the eight-qubit graph state, from the explicit state vector"""
    psi = np.array([(-1.0) ** (sum(popcount(adj[i] & x) for i in range(8) if x >> i & 1) // 2 % 2) / 16
                    for x in range(256)])
    T = psi.reshape([2] * 8)                                           # axis k <-> qubit 7 - k
    out = [0.0] * 256
    for X in range(1, 256):
        A = [7 - i for i in range(8) if X >> i & 1]
        B = [7 - i for i in range(8) if not X >> i & 1]
        sv = np.linalg.svd(np.transpose(T, A + B).reshape(2 ** len(A), 2 ** len(B)), compute_uv=False) ** 2
        sv = sv[sv > 1e-12]
        out[X] = float(-(sv * np.log2(sv)).sum())
    return out


rng_r = random.Random(20261029)
errs = []
for _ in range(12):
    adj_ = graph_adj(rng_r.randrange(1 << 21), rng_r.randrange(1, 128))
    num = graph_state_entropies(adj_)
    errs.append(max(abs(num[X] - s_graph(adj_, X)) for X in range(1, 256)))
report('graph-state entropies: the GF(2) rank formula agrees with the explicit state vectors of 12 random eight-qubit '
       'graph states on all 255 subsystems', max(errs) < 1e-9, 'max error %.1e' % max(errs))
T12 = [1 << i for i in range(7)] + [1 << i | 1 << j for i, j in itertools.combinations(range(7), 2)]
T3 = [1 << i | 1 << j | 1 << k for i, j, k in itertools.combinations(range(7), 3)]
hist, n_codes = {}, 0
for g in REPS:
    for nbr in range(1, 128):
        adj_ = graph_adj(g, nbr)
        if any(1 + s_graph(adj_, T) - s_graph(adj_, T | RB) for T in T12):   # I(R;T) > 0 for some |T| <= 2
            continue
        n_codes += 1
        f3 = sum(1 for T in T3 if 1 + s_graph(adj_, T) - s_graph(adj_, T | RB) > 0)
        hist[f3] = hist.get(f3, 0) + 1
report('qubit stabilizer codes that encode one qubit into seven with distance >= 3, as graph states on R and the sites '
       '(every stabilizer state is local-Clifford equivalent to a graph state; local Cliffords keep all entropies): '
       '%d of the 1,044 x 127 graph states; uncorrectable erasures of three sites: %s; minimum %d, as for Steane\'s '
       'code, against %d for the color code' % (n_codes, dict(sorted(hist.items())), min(hist), FATAL[3][3]),
       n_codes == 3414 and min(hist) == 7 and FATAL[2][3] == 7 and FATAL[3][3] == 6)

print('ALL OK' if ok_all else 'FAILURE')
sys.exit(0 if ok_all else 1)
