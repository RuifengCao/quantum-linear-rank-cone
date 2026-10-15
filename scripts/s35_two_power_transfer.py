#!/usr/bin/env python3
"""s35: six parties, powers of two -- the symplectic transfer over Z_{2^a} (R44; Theorem 14 of
docs/R32_reduction_and_local_dimension.md, section 4.12).

Theorem 12 (R42) proves f12 >= 0 on the qubit cone Stab^(2)_6 by a "symplectic transfer": Lemmas T1-T5 bound the
pieces of the stabilizer group on the six non-Fano lines, and Lemma 4 (LA)_2 is applied to a dual arrangement.  Over
Z_{2^a} the same lemmas hold with log|.| in place of dim (finite abelian groups with a nondegenerate alternating
form), but the dual arrangement lives in a group of exponent dividing 2^a, where (LA)_2 fails.  Lemma 4(G') replaces it:
if 2^a V = 0, then H(B1, B2, B3) <= (3 - 1/a) H(C) + 3 alpha + 3 gamma + sum_i (zeta_i + eta_i), because the kernel of
the map Phi of the proof contains the diagonal of the 2-torsion of q(C_A), of size >= |q(C_A)|^(1/a).  The transfer
then gives

    f^(a) = f12 + (1 - 1/a) U_Z >= 0   on Stab^(Z_{2^a})_6,

with f^(a)(S_NF) = -1/a and f^(a) = 1 - j/a bits at the seven-point state over Z_{2^j}.  Hence the non-Fano state is
in no Z_{2^a} cone, and the Z_{2^j} seven-point state is in none with a < j.

Notation: the note writes the local dimension as N = 2^a; in this code the exponent is the variable k, and k(T) of the
note is kk in chain().  Everything is an integer: entropies in bits, and a f^(a) = a f12 + (a - 1) U_Z (fk_form).

This script holds checks, not proofs:
  * linear algebra over Z_{2^k}: subgroup sizes (diagonalisation over the local ring), kernels, sums, intersections,
    the rank of an alternating form restricted to a subgroup (the index of its radical);
  * stabilizer states over Z_{2^k} as Lagrangian subgroups L of Z_N^{2m} (N = 2^k): weighted graph states, CSS states
    of random (also non-free) codes, standard non-free Lagrangians after random symplectic maps, families pushed
    towards the seven-point configuration, the seven-point states over Z_{2^j} (j <= k) and direct sums;
    S(X) = log2 |proj_X L| - n_X k (Theorem 1(d), step 1);
  * chain(st): every step of the proof on one state, without compressing it (each quantity of the compressed state
    is computed from preimages: log |L_P| is subtracted once per system): the sizes of the line pieces and of their
    restrictions, orthogonality, Lemmas T1-T4 over Z_{2^k}, the dual arrangement of Lemma T5 built explicitly as the
    images of the D_P in (sum of the D_P) / Z_0, Lemma 4(G') on it, the algebra of the assembly, and a f^(a) >= 0;
  * exact values at S_NF, S_F and the seven-point states over Z_{2^j} and Z_d;
  * state vectors: the entropy formula against Z_4 / Z_8 graph and CSS states (including non-free codes);
  * Lemma 4(G') on random and near-J - I subgroup configurations of Z_{2^k}^m, and the forms it does not claim
    (coefficient 3 - 1/a on exponent 2^(a+1), or 2 on exponent 4) going negative;
  * cc/s38_minform_zn.c: annealing over CSS codes (anneal) and exhaustive families of codes B + 2^s P with P a 0/1
    matrix (lifts), with forms the proof does not give (weakened_forms) as a sensitivity test.

  python3 scripts/s35_two_power_transfer.py --check [--full --trials 60]   # gates G40 / G40b
  python3 scripts/s35_two_power_transfer.py --anneal "2:4:1,1,1,1,1,1,1;3:4:1,1,1,1,1,1,1" --seconds 120 [--pull 0.3]
"""
import argparse, hashlib, itertools, json, math, os, shutil, subprocess, sys, tempfile, time
from fractions import Fraction
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, '..')
sys.path.insert(0, HERE); sys.path.insert(0, ROOT)
import s34_six_party_survey as s34

LINES = {n: s34._mask(t) for n, t in s34.TR_LINES.items()}
PAIRS = [s34._mask(p) for p in s34.TR_PAIRS]
NAMES = s34.TR_NAMES
A_, B_, C_, W_, X_, Y_, Z_ = range(7)
CORNERS = (A_, B_, C_, Z_)
THROUGH = {P: [l for l, t in s34.TR_LINES.items() if P in t] for P in range(7)}


# ----------------------------------------------------------------------------------------------- Z_{2^k} linear algebra
def _lowbit(A, N):
    """lowest set bit of every entry (N for the zero entries); entries in [0, N)."""
    lb = A & (-A)
    return np.where(A == 0, N, lb)


def diagonalise(M, k, track=False):
    """Diagonalise the integer matrix M over Z_{2^k} by unimodular row and column operations, U M V = D.  Returns the
    2-adic valuations of the nonzero diagonal entries (each < k) and, with track, U (row operations only).  Every pivot
    is an entry of least valuation in the remaining rows, so the row reduction never divides by a non-unit."""
    N = 1 << k
    A = np.array(M, dtype=np.int64) % N
    if A.ndim != 2 or A.size == 0:
        r = A.shape[0] if A.ndim == 2 else 0
        return [], (np.eye(r, dtype=np.int64) if track else None)
    r, c = A.shape
    U = np.eye(r, dtype=np.int64) if track else None
    vals, row = [], 0
    while row < r:
        lb = _lowbit(A[row:], N)
        idx = int(np.argmin(lb))
        i, j = divmod(idx, c)
        if lb[i, j] >= N:
            break
        i += row
        v = int(lb[i - row, j]).bit_length() - 1
        if i != row:
            A[[row, i]] = A[[i, row]]
            if track:
                U[[row, i]] = U[[i, row]]
        uinv = pow(int(A[row, j]) >> v, -1, N)
        col = A[:, j].copy()
        col[row] = 0                                   # rows above `row` are zero in column j; rows below have val >= v
        f = ((col >> v) * uinv) % N
        if f.any():
            A = (A - np.outer(f, A[row])) % N
            if track:
                U = (U - np.outer(f, U[row])) % N
        piv = A[row, j]
        A[row] = 0                                     # column operations clear the rest of the pivot row
        A[row, j] = piv
        vals.append(v)
        row += 1
    return vals, U


def log_size(M, k):
    """log2 of the size of the subgroup of Z_{2^k}^c generated by the rows of M."""
    if M is None or not len(M) or not np.asarray(M).shape[1]:
        return 0
    vals, _ = diagonalise(M, k)
    return sum(k - v for v in vals)


def left_kernel(M, k):
    """generators (rows) of {x : x M = 0 mod 2^k}."""
    M = np.asarray(M, dtype=np.int64)
    r = M.shape[0]
    if r == 0:
        return np.zeros((0, 0), dtype=np.int64)
    if M.shape[1] == 0:
        return np.eye(r, dtype=np.int64)
    N = 1 << k
    vals, U = diagonalise(M, k, track=True)
    rows = [((1 << (k - v)) * U[i]) % N for i, v in enumerate(vals) if v > 0]
    rows += [U[i] for i in range(len(vals), r)]
    return np.array(rows, dtype=np.int64).reshape(len(rows), r)


def echelon(G, k):
    """a generating set of at most (number of columns) rows for the row space of G (row echelon form over the local
    ring Z_{2^k}, least-valuation pivot in each column)."""
    G = np.asarray(G, dtype=np.int64)
    N = 1 << k
    if not len(G) or not G.shape[1]:
        return np.zeros((0, G.shape[1] if G.ndim == 2 else 0), dtype=np.int64)
    A = G % N
    A = A[A.any(axis=1)]
    out = []
    for j in range(A.shape[1]):
        if not len(A):
            break
        lb = _lowbit(A[:, j], N)
        i = int(np.argmin(lb))
        if lb[i] >= N:
            continue
        v = int(lb[i]).bit_length() - 1
        uinv = pow(int(A[i, j]) >> v, -1, N)
        piv = A[i].copy()
        rest = np.delete(A, i, axis=0)
        f = ((rest[:, j] >> v) * uinv) % N
        rest = (rest - np.outer(f, piv)) % N
        out.append(piv)
        A = rest[rest.any(axis=1)]
    return np.array(out, dtype=np.int64).reshape(len(out), G.shape[1])


def intersect(G1, G2, k):
    """generators of rowspace(G1) & rowspace(G2)."""
    N = 1 << k
    if not len(G1) or not len(G2):
        return np.zeros((0, max(np.asarray(G1).shape[1], np.asarray(G2).shape[1])), dtype=np.int64)
    K = left_kernel(np.concatenate([G1, (-G2) % N], axis=0), k)
    if not len(K):
        return np.zeros((0, G1.shape[1]), dtype=np.int64)
    return echelon((K[:, :len(G1)] @ G1) % N, k)


def gram(U, V, n, k):
    """the standard alternating form between the rows of U and of V (n qudits, x-part first): U_x V_z^T - U_z V_x^T."""
    N = 1 << k
    return (U[:, :n] @ V[:, n:].T - U[:, n:] @ V[:, :n].T) % N


def form_rank(U, n, k):
    """rank of the alternating form restricted to the subgroup generated by the rows of U: log2 |U / (U & U^perp)|,
    which is log2 of the size of the row space of the Gram matrix."""
    return log_size(gram(U, U, n, k), k) if len(U) else 0


# ----------------------------------------------------------------------------------------------- states
class ZnState:
    """a stabilizer state over Z_{2^k}: generator rows G of its Lagrangian L in Z_N^{2m} (columns x_0..x_{m-1},
    z_0..z_{m-1}) and the system (0..6) of every qudit."""

    def __init__(self, G, owner, k):
        self.k, self.N = k, 1 << k
        self.owner = list(owner)
        self.m = len(owner)
        self.G = echelon(np.asarray(G, dtype=np.int64).reshape(-1, 2 * self.m) % self.N, k)
        self.cols = {P: [v for v in range(self.m) if owner[v] == P] + [self.m + v for v in range(self.m) if owner[v] == P]
                     for P in range(7)}
        self.n = [self.owner.count(P) for P in range(7)]

    def colset(self, T):
        return [c for P in range(7) if T >> P & 1 for c in self.cols[P]]

    def supported(self, T, basis=None):
        """generators of the elements of span(basis) (default: L) supported in the systems of T."""
        Gm = self.G if basis is None else basis
        out = self.colset(127 ^ T)
        if not len(Gm) or not out:
            return Gm
        K = left_kernel(Gm[:, out], self.k)
        return echelon((K @ Gm) % self.N, self.k) if len(K) else np.zeros((0, Gm.shape[1]), dtype=np.int64)

    def proj(self, M, P):
        return M[:, self.cols[P]] if len(M) else np.zeros((0, len(self.cols[P])), dtype=np.int64)

    def is_lagrangian(self):
        n = self.m
        return bool(log_size(self.G, self.k) == n * self.k and not gram(self.G, self.G, n, self.k).any())

    def entropies(self):
        """S(T) = log2 |proj_T L| - n_T k for the 128 subsets T of the seven systems (bits)."""
        S = np.zeros(128, dtype=np.int64)
        for T in range(1, 127):
            S[T] = log_size(self.G[:, self.colset(T)], self.k) - sum(self.n[P] for P in range(7) if T >> P & 1) * self.k
        return S


def _place(rows, cols, width):
    out = np.zeros((len(rows), width), dtype=np.int64)
    if len(rows):
        out[:, cols] = rows
    return out


def graph_state(Gam, owner, k):
    """weighted graph state: L = {(x, x Gam)} for a symmetric Gam over Z_N."""
    m = len(owner)
    Gam = np.asarray(Gam, dtype=np.int64) % (1 << k)
    return ZnState(np.concatenate([np.eye(m, dtype=np.int64), Gam], axis=1), owner, k)


def css_state(M, owner, k):
    """CSS state of the code C spanned by the rows of M over Z_N (columns = qudits): L = C x C^perp."""
    N = 1 << k
    M = np.asarray(M, dtype=np.int64) % N
    m = M.shape[1]
    C = echelon(M, k)
    Cp = left_kernel(C.T, k) if len(C) else np.eye(m, dtype=np.int64)
    Cp = echelon(Cp, k)
    G = np.concatenate([np.concatenate([C, np.zeros_like(C)], axis=1),
                        np.concatenate([np.zeros((len(Cp), m), dtype=np.int64), Cp], axis=1)], axis=0)
    return ZnState(G, owner, k)


def direct_sum(s1, s2):
    """tensor product of two states on the same seven systems."""
    assert s1.k == s2.k
    m1, m2 = s1.m, s2.m
    G1 = np.concatenate([s1.G[:, :m1], np.zeros((len(s1.G), m2), dtype=np.int64),
                         s1.G[:, m1:], np.zeros((len(s1.G), m2), dtype=np.int64)], axis=1)
    G2 = np.concatenate([np.zeros((len(s2.G), m1), dtype=np.int64), s2.G[:, :m2],
                         np.zeros((len(s2.G), m1), dtype=np.int64), s2.G[:, m2:]], axis=1)
    return ZnState(np.concatenate([G1, G2], axis=0), s1.owner + s2.owner, s1.k)


SEVEN = np.array(s34.SEVEN_POINTS, dtype=np.int64).T            # 3 x 7, columns A, B, C, W, X, Y, Z


def seven_point_state(j, k):
    """the seven-point state over Z_{2^j} (j <= k) realised with Z_{2^k} qudits: the CSS state of 2^(k-j) times the
    seven-point code (Z_{2^j} = 2^(k-j) Z_{2^k}), one qudit per system; same entropies as over Z_{2^j}."""
    return css_state((SEVEN << (k - j)) % (1 << k), list(range(7)), k)


def apply_symplectic(G, m, k, rng, moves, qudits=None):
    """random symplectic maps on the given qudits (default: all): shears z += x S and x += z S (S symmetric),
    x -> x A with z -> z A^(-T) for elementary A, unit scalings, and the Fourier map (x, z) -> (-z, x)."""
    N = 1 << k
    G = np.array(G, dtype=np.int64) % N
    q = list(range(m)) if qudits is None else list(qudits)
    for _ in range(moves):
        t = rng.integers(0, 5)
        if t == 0 or t == 1:
            S = np.zeros((m, m), dtype=np.int64)
            for a in q:
                for b in q:
                    if a <= b and rng.random() < 0.5:
                        S[a, b] = S[b, a] = rng.integers(0, N)
            if t == 0:
                G[:, m:] = (G[:, m:] + G[:, :m] @ S) % N
            else:
                G[:, :m] = (G[:, :m] + G[:, m:] @ S) % N
        elif t == 2 and len(q) > 1:
            a, b = rng.choice(q, size=2, replace=False)
            c = int(rng.integers(1, N))
            G[:, b] = (G[:, b] + c * G[:, a]) % N                   # x_b += c x_a
            G[:, m + a] = (G[:, m + a] - c * G[:, m + b]) % N       # z_a -= c z_b
        elif t == 3:
            a = int(rng.choice(q))
            u = int(rng.integers(0, N // 2)) * 2 + 1
            G[:, a] = (G[:, a] * u) % N
            G[:, m + a] = (G[:, m + a] * pow(u, -1, N)) % N
        else:
            a = int(rng.choice(q))
            xa, za = G[:, a].copy(), G[:, m + a].copy()
            G[:, a] = (-za) % N
            G[:, m + a] = xa
    return G


def local_scramble(st, rng, moves=8):
    """random symplectic maps inside every system: the entropies and every quantity of the proof are unchanged."""
    G = st.G.copy()
    for P in range(7):
        qs = [v for v in range(st.m) if st.owner[v] == P]
        if qs:
            G = apply_symplectic(G, st.m, st.k, rng, moves, qs)
    return ZnState(G, st.owner, st.k)


def random_lagrangian(rng, sizes, k, moves=12):
    """a standard Lagrangian of random type (generators 2^a e_i in x and 2^(k-a) e_i in z, a in 0..k, possibly
    non-free) after random global symplectic maps; qudits distributed by sizes."""
    N = 1 << k
    owner = [P for P in range(7) for _ in range(sizes[P])]
    m = len(owner)
    rows = []
    for i in range(m):
        a = int(rng.integers(0, k + 1))
        if a < k:
            v = np.zeros(2 * m, dtype=np.int64); v[i] = 1 << a; rows.append(v)
        if a > 0:
            v = np.zeros(2 * m, dtype=np.int64); v[m + i] = 1 << (k - a); rows.append(v)
    G = apply_symplectic(np.array(rows), m, k, rng, moves)
    perm = rng.permutation(m)                                       # random assignment of qudits to systems
    G = G[:, list(perm) + [m + p for p in perm]]
    return ZnState(G, owner, k)


def random_graph(rng, sizes, k, density=None, diagonal=False):
    N = 1 << k
    owner = [P for P in range(7) for _ in range(sizes[P])]
    m = len(owner)
    dens = rng.uniform(0.2, 0.9) if density is None else density
    Gam = np.where(rng.random((m, m)) < dens, rng.integers(0, N, size=(m, m)), 0)
    if rng.random() < 0.5:                                          # 2-adic weights
        Gam = (Gam << rng.integers(0, k, size=(m, m))) % N
    Gam = np.triu(Gam, 1)
    Gam = Gam + Gam.T
    if diagonal:
        Gam = Gam + np.diag(rng.integers(0, N, size=m))
    return graph_state(Gam, owner, k)


def random_css(rng, sizes, k, rows=None):
    """CSS state of a random code; rows scaled by random powers of 2 (non-free codes)."""
    N = 1 << k
    owner = [P for P in range(7) for _ in range(sizes[P])]
    m = len(owner)
    r = int(rng.integers(1, m + 1)) if rows is None else rows
    M = rng.integers(0, N, size=(r, m))
    M = (M << rng.integers(0, k, size=(r, 1))) % N
    if rng.random() < 0.3:
        M = (M << rng.integers(0, 2, size=(1, m))) % N              # some columns scaled as well
    return css_state(M, owner, k)


def structured_css(rng, k):
    """CSS state of a random family over Z_{2^k} pushed towards the seven-point configuration: A, B, C random cyclic
    or two-generated subgroups of Z_N^r, Z mostly inside A + B + C, and W, X, Y mostly inside (A + B) & (C + Z),
    (A + C) & (B + Z), (B + C) & (A + Z); random 2-adic scalings.  (Over Z_{2^k} the intersections need not be
    cyclic: the seven-point configuration itself, with W, X, Y on 2-torsion, is the tight case.)"""
    N = 1 << k
    while True:
        r = int(rng.integers(3, 6))
        rnd = lambda c: (rng.integers(0, N, size=(c, r)) << rng.integers(0, k, size=(c, 1))) % N
        U = {}
        for P in (A_, B_, C_):
            U[P] = rnd(int(rng.integers(1, 3)))
        ABC = np.concatenate([U[A_], U[B_], U[C_]], axis=0)
        if rng.random() < 0.8:
            U[Z_] = (rng.integers(0, N, size=(int(rng.integers(1, 3)), len(ABC))) @ ABC) % N
        else:
            U[Z_] = rnd(1)
        for D, (a, b, c, d) in ((W_, (A_, B_, C_, Z_)), (X_, (A_, C_, B_, Z_)), (Y_, (B_, C_, A_, Z_))):
            I = intersect(np.concatenate([U[a], U[b]]), np.concatenate([U[c], U[d]]), k)
            cnt = int(rng.integers(1, 3))
            if len(I) and rng.random() < 0.85:
                V = (rng.integers(0, N, size=(cnt, len(I))) @ I) % N
            else:
                AB = np.concatenate([U[a], U[b]])
                V = (rng.integers(0, N, size=(cnt, len(AB))) @ AB) % N
            if rng.random() < 0.15:
                V = np.concatenate([V, rnd(1)], axis=0)
            U[D] = V
        cols, owner = [], []
        for P in range(7):
            for v in U[P]:
                cols.append(v); owner.append(P)
        M = np.array(cols, dtype=np.int64).T
        st = css_state(M, owner, k)
        S = st.entropies()
        if all(S[1 << P] > 0 for P in range(7)):
            return st


def seven_point_perturbed(rng, k):
    """the seven-point code over Z_{2^k} with random extra columns, random column multipliers and an extra row."""
    N = 1 << k
    M = SEVEN.copy()
    mult = np.where(rng.random(7) < 0.3, rng.integers(1, N, size=7), 1)
    M = (M * mult) % N
    if rng.random() < 0.5:
        M = np.concatenate([M, rng.integers(0, N, size=(1, 7)) * (rng.random() < 0.5) * 2 % N], axis=0)
    owner = list(range(7))
    extra = int(rng.integers(0, 3))
    for _ in range(extra):
        P = int(rng.integers(0, 7))
        c = (rng.integers(0, N, size=(M.shape[0], 1)) << int(rng.integers(0, k))) % N
        M = np.concatenate([M, c], axis=1); owner.append(P)
    return css_state(M, owner, k)


# ----------------------------------------------------------------------------------------------- the forms
def uz_form():
    """U_Z = 4 S(Z) - d(Z, AYZ) - d(Z, BXZ) - d(Z, CWZ) as 128 coefficients."""
    F = np.zeros(128, dtype=np.int64)

    def k_(T, c):
        for P in range(7):
            if T >> P & 1:
                F[1 << P] += c
        F[T] -= c
    F[1 << Z_] += 4
    for l in ('AYZ', 'BXZ', 'CWZ'):
        k_(LINES[l], -1)
        k_(LINES[l] & ~(1 << Z_), 1)
    return F


def fk_form(k):
    """a f^(a) = a f12 + (a - 1) U_Z for a = k, as 128 integer coefficients (f^(1) = f12)."""
    return k * s34.transfer_form() + (k - 1) * uz_form()


def components():
    """the pieces of f12 as 128-coefficient vectors: U_Z, U_W, U_X, U_Y, the e-terms inside U_W, U_X, U_Y, the rest
    (k(WXY) + 6 k(ABC) + sum over corner pairs - sum over lines - 3 k(ABW) - 3 k(CWZ)) and S(Z);
    f12 = 5 U_Z + 4 U_W + U_X + U_Y + rest."""
    out = {}

    def new():
        return np.zeros(128, dtype=np.int64)

    def k_(F, T, c):
        for P in range(7):
            if T >> P & 1:
                F[1 << P] += c
        F[T] -= c

    def d_(F, P, l, c):
        k_(F, LINES[l], c)
        k_(F, LINES[l] & ~(1 << P), -c)

    def e_(F, P, l1, l2, c):
        F[1 << P] += 2 * c
        d_(F, P, l1, -c)
        d_(F, P, l2, -c)
    out['UZ'] = uz_form()
    for name, P, l, (q1, a1, b1), (q2, a2, b2) in (('W', W_, 'ABW', (A_, 'ACX', 'AYZ'), (B_, 'BCY', 'BXZ')),
                                                     ('X', X_, 'ACX', (A_, 'ABW', 'AYZ'), (C_, 'BCY', 'CWZ')),
                                                     ('Y', Y_, 'BCY', (B_, 'ABW', 'BXZ'), (C_, 'ACX', 'CWZ'))):
        E = new(); e_(E, q1, a1, b1, 1); e_(E, q2, a2, b2, 1)
        U = new(); U[1 << P] += 2; d_(U, P, l, -1)
        out['e' + name] = E
        out['U' + name] = U + E
    R = new()
    k_(R, s34._mask((W_, X_, Y_)), 1)
    k_(R, s34._mask((A_, B_, C_)), 6)
    for p in PAIRS:
        k_(R, p, 1)
    for T in LINES.values():
        k_(R, T, -1)
    k_(R, LINES['ABW'], -3); k_(R, LINES['CWZ'], -3)
    out['rest'] = R
    SZ = new(); SZ[1 << Z_] = 1
    out['SZ'] = SZ
    return out


def weakened_forms(k):
    """a f^(a) (a = k) and forms that the proof does not give (each is a f^(a) minus a term that is >= 0 on stabilizer
    states), for sensitivity tests: {name: 128 integer coefficients}."""
    c = components()
    base = fk_form(k)
    return {'a f^(a)': base,
            'U_Z coefficient (1 - 1/a) - 1/(2a), times 2a': 2 * base - c['UZ'],
            'f^(a) without the e-terms': base - k * (4 * c['eW'] + c['eX'] + c['eY']),
            'f^(a) with 3 S(Z) in U_Z': base - k * 5 * c['SZ'] - (k - 1) * c['SZ'],
            'f^(a) with U_W instead of 4 U_W': base - k * 3 * c['UW']}


def seven_point_entropy(d):
    """the seven-point state over Z_d (Theorem 1(d), step 3): {T: (a, b)} with S(T) = log2(a / b) bits, where
    a = |proj_T C| |proj_{E-T} C| and b = |C|."""
    sz = s34.seven_point_sizes(d)
    return {T: (sz[T] * sz[127 ^ T], sz[127]) for T in range(128)}


def fk_seven_point(d, k):
    """2^(a f^(a)) (a = k) at the seven-point state over Z_d, exactly (a Fraction); f^(a) = log2 of it / a."""
    F = fk_form(k)
    e = seven_point_entropy(d)
    q = Fraction(1)
    for T in range(1, 127):
        if F[T]:
            num, den = e[T]
            q *= Fraction(num, den) ** int(F[T])
    return q


def lemma_g_prime(H, k):
    """Lemma 4(G') for M = 3, t = 2 on a rank function H (dict or array over the 128 subsets of A1 A2 A3 B1 B2 B3 C =
    bits 0..6 in the order A, B, C, Y, X, W, Z of the transfer): k * RHS - k * LHS, an integer when H is (bits)."""
    a = [1 << A_, 1 << B_, 1 << C_]
    b = [1 << Y_, 1 << X_, 1 << W_]
    c = 1 << Z_
    A = a[0] | a[1] | a[2]
    Hm = lambda T: int(H[T]) if T else 0
    alpha = sum(Hm(x) for x in a) - Hm(A)
    gamma = Hm(A | c) - Hm(A)
    zeta = [Hm(b[i] | (A & ~a[i])) - Hm(A & ~a[i]) for i in range(3)]
    eta = [Hm(b[i] | a[i] | c) - Hm(a[i] | c) for i in range(3)]
    lhs = Hm(b[0] | b[1] | b[2])
    return (3 * k - 1) * Hm(c) + k * (3 * alpha + 3 * gamma + sum(zeta) + sum(eta)) - k * lhs


# ----------------------------------------------------------------------------------------------- the proof, step by step
def chain(st, S=None):
    """Every step of the proof of Theorem 14 on one state over Z_{2^k}, without compressing it.  The compressed state
    (symplectic reduction L_P^perp / L_P on every system) has L'_T = image of L_T, so its quantities are those below
    with log2 |L_P| =: lam_P subtracted once per system: k(T) = log|L_T| - sum_{P in T} lam_P, d(P, l) =
    log|pi_P L_l| - lam_P, 2 S(P) = log|L_P^perp / L_P|, delta_P = log|pi_P Z_0| - lam_P, z(T) = log|(Z_0)_T| -
    sum_{P in T} lam_P; ranks of forms are unchanged (the preimages contain L_P, which lies in their radicals).
    Returns (values, list of failed steps)."""
    k, N = st.k, st.N
    S = st.entropies() if S is None else S
    s = [int(S[1 << P]) for P in range(7)]
    bad = []
    lam = [log_size(st.supported(1 << P), k) for P in range(7)]
    if any(s[P] != st.n[P] * k - lam[P] for P in range(7)):
        bad.append('S(P) != n_P k - log|L_P|')
    kk = lambda T: sum(s[P] for P in range(7) if T >> P & 1) - int(S[T])
    for T in range(1, 128):
        if log_size(st.supported(T), k) - sum(lam[P] for P in range(7) if T >> P & 1) != kk(T):
            bad.append(f'log|L_T| - sum lam != k(T) at {T}')
            break
    d = lambda P, l: kk(LINES[l]) - kk(LINES[l] & ~(1 << P))
    pieces = {l: st.supported(T) for l, T in LINES.items()}
    R = {l: {P: st.proj(pieces[l], P) for P in s34.TR_LINES[l]} for l in LINES}
    for l in LINES:
        for P in s34.TR_LINES[l]:
            if log_size(R[l][P], k) - lam[P] != d(P, l):
                bad.append(f'size of the restriction of {l} at {NAMES[P]}')
    for P in range(7):
        n = st.n[P]
        for i, l1 in enumerate(THROUGH[P]):
            for l2 in THROUGH[P][i + 1:]:
                if len(R[l1][P]) and len(R[l2][P]) and gram(R[l1][P], R[l2][P], n, k).any():
                    bad.append(f'restrictions of {l1} and {l2} not orthogonal at {NAMES[P]}')
    rk = {(l, P): form_rank(R[l][P], st.n[P], k) for l in LINES for P in s34.TR_LINES[l]}
    for P in CORNERS:                                               # Lemma T1
        ls = THROUGH[P]
        for l in ls:
            o = [m_ for m_ in ls if m_ != l]
            if rk[(l, P)] > 2 * s[P] - d(P, o[0]) - d(P, o[1]):
                bad.append(f'Lemma T1 (rank) at {NAMES[P]} for {l}')
        Vs = [R[l][P] for l in ls if len(R[l][P])]
        if Vs and log_size(np.concatenate(Vs), k) - lam[P] > 4 * s[P] - sum(d(P, l) for l in ls):
            bad.append(f'Lemma T1 (sum) at {NAMES[P]}')
    for l, t in s34.TR_LINES.items():                                # Lemma T2
        P1, P2 = [P for P in t if P in CORNERS]
        Qd = next(P for P in t if P not in CORNERS)
        if rk[(l, Qd)] > rk[(l, P1)] + rk[(l, P2)]:
            bad.append(f'Lemma T2 on {l}')
    nz = [x for x in pieces.values() if len(x)]
    Z0 = echelon(np.concatenate(nz), k) if nz else np.zeros((0, 2 * st.m), dtype=np.int64)
    lz0 = log_size(Z0, k)
    Dg = {P: st.proj(Z0, P) for P in range(7)}
    delta = [log_size(Dg[P], k) - lam[P] for P in range(7)]
    e = lambda P, l1, l2: 2 * s[P] - d(P, l1) - d(P, l2)
    U = {Z_: 4 * s[Z_] - d(Z_, 'AYZ') - d(Z_, 'BXZ') - d(Z_, 'CWZ'),
         W_: 2 * s[W_] - d(W_, 'ABW') + e(A_, 'ACX', 'AYZ') + e(B_, 'BCY', 'BXZ'),
         X_: 2 * s[X_] - d(X_, 'ACX') + e(A_, 'ABW', 'AYZ') + e(C_, 'BCY', 'CWZ'),
         Y_: 2 * s[Y_] - d(Y_, 'BCY') + e(B_, 'ABW', 'BXZ') + e(C_, 'ACX', 'CWZ')}
    for P, u in U.items():                                          # Lemma T3
        if delta[P] > u:
            bad.append(f'Lemma T3: delta_{NAMES[P]} = {delta[P]} > {u}')
    zs = {0: 0}
    for T in range(1, 128):
        out = st.colset(127 ^ T)
        zs[T] = lz0 - (log_size(Z0[:, out], k) if out and len(Z0) else 0) - sum(lam[P] for P in range(7) if T >> P & 1)
    if any(zs[T] > kk(T) for T in range(1, 128)):
        bad.append('Lemma T4(i): z(T) > k(T)')
    if any(zs[T] != kk(T) for T in LINES.values()):
        bad.append('Lemma T4(ii): z(l) != k(l)')
    Q = s34._mask((A_, B_, C_, Z_))
    lz = kk(LINES['ABW']) + kk(LINES['CWZ']) - delta[W_]
    if zs[Q] < lz:
        bad.append(f'Lemma T4(iii): z(ABCZ) = {zs[Q]} < {lz}')
    # Lemma T5: the images R_P of the D_P in (sum of the D_P) / Z_0, built explicitly
    width = 2 * st.m
    emb = {P: _place(Dg[P], st.cols[P], width) for P in range(7)}
    H = np.zeros(128, dtype=np.int64)
    for T in range(1, 128):
        gens = [Z0] + [emb[P] for P in range(7) if T >> P & 1 and len(emb[P])]
        gens = [g for g in gens if len(g)]
        H[T] = (log_size(np.concatenate(gens), k) if gens else 0) - lz0
        if H[T] != sum(delta[P] for P in range(7) if T >> P & 1) - zs[T]:
            bad.append(f'Lemma T5 at {T}')
            break
    slack = lemma_g_prime(H, k)                                     # Lemma 4(G') on the explicit arrangement
    if slack < 0:
        bad.append(f"Lemma 4(G') on the dual arrangement: a (RHS - LHS) = {slack}")
    # the algebra of the proof: (LA) in the delta / z / k(T) form, and a f^(a) = that + a (nonnegative brackets)
    zz = lambda *ps: zs[s34._mask(ps)]
    kz = lambda *ps: kk(s34._mask(ps)) - zs[s34._mask(ps)]
    assembled = ((6 * k - 1) * delta[Z_] + k * (delta[W_] + delta[X_] + delta[Y_] + zz(W_, X_, Y_) + 6 * zz(A_, B_, C_)
                                                 + sum(zs[p] for p in PAIRS) - 3 * zs[Q] - sum(kk(T) for T in LINES.values())))
    if assembled != slack:
        bad.append(f'(LA) in the delta/z form: {assembled} != {slack}')
    fk = int(fk_form(k) @ S)
    brackets = ((6 * k - 1) * (U[Z_] - delta[Z_]) + 4 * k * (U[W_] - delta[W_]) + k * (U[X_] - delta[X_])
                + k * (U[Y_] - delta[Y_]) + k * (kz(W_, X_, Y_) + 6 * kz(A_, B_, C_) + sum(kk(p) - zs[p] for p in PAIRS))
                + 3 * k * (zs[Q] - lz))
    if fk != assembled + brackets:
        bad.append(f'a f^(a) = {fk} != assembled {assembled} + brackets {brackets}')
    if fk < 0:
        bad.append(f'a f^(a) = {fk}')
    f12 = int(s34.transfer_form() @ S)
    return dict(kfk=fk, f12=f12, U=U, delta=delta, lam=lam, LA_slack=slack, zQ=(zs[Q], lz)), bad


# ----------------------------------------------------------------------------------------------- state vectors
def state_vector_entropies(st_kind, data, owner, k):
    """entropies (bits, 128 subsets) of a graph state (data = Gam, no diagonal) or a CSS state (data = M) over Z_N from
    its state vector, one qudit per entry of owner."""
    N = 1 << k
    m = len(owner)
    xs = np.array(list(itertools.product(range(N), repeat=m)), dtype=np.int64)
    if st_kind == 'graph':
        Gam = np.asarray(data, dtype=np.int64) % N
        ph = np.zeros(len(xs), dtype=np.int64)
        for i in range(m):
            for j in range(i + 1, m):
                if Gam[i, j]:
                    ph += Gam[i, j] * xs[:, i] * xs[:, j]
        psi = np.exp(2j * np.pi * (ph % N) / N) / np.sqrt(float(N) ** m)
    else:
        M = np.asarray(data, dtype=np.int64) % N
        code = {tuple((np.array(x) @ M) % N) for x in itertools.product(range(N), repeat=M.shape[0])}
        psi = np.zeros(N ** m, dtype=complex)
        for c in code:
            psi[int(np.ravel_multi_index(c, [N] * m))] = 1.0
        psi /= np.sqrt(len(code))
    t = psi.reshape([N] * m)
    S = np.zeros(128)
    for T in range(1, 127):
        a = [v for v in range(m) if T >> owner[v] & 1]
        b = [v for v in range(m) if v not in a]
        if not a or not b:
            continue
        sv = np.linalg.svd(np.transpose(t, a + b).reshape(N ** len(a), N ** len(b)), compute_uv=False) ** 2
        sv = sv[sv > 1e-13]
        S[T] = float(-(sv * np.log2(sv)).sum())
    return S


def state_vector_check(rng, k=2, trials=6):
    """the entropy formula S(T) = log2 |proj_T L| - n_T k against state vectors: random graph states and CSS states of
    random (non-free) codes over Z_{2^k}, five or six qudits on the seven systems.  Returns the largest deviation."""
    N = 1 << k
    dev = 0.0
    for t in range(trials):
        m = 5 if k >= 3 else 6
        owner = sorted(int(x) for x in rng.integers(0, 7, size=m))
        if t % 2 == 0:
            Gam = np.triu(rng.integers(0, N, size=(m, m)), 1)
            Gam = Gam + Gam.T
            st = graph_state(Gam, owner, k)
            Ssv = state_vector_entropies('graph', Gam, owner, k)
        else:
            r = int(rng.integers(1, m))
            M = (rng.integers(0, N, size=(r, m)) << rng.integers(0, k, size=(r, 1))) % N
            st = css_state(M, owner, k)
            Ssv = state_vector_entropies('css', M, owner, k)
        S = st.entropies()
        dev = max(dev, float(np.abs(S[1:127] - Ssv[1:127]).max()))
    return dev


# ----------------------------------------------------------------------------------------------- Lemma 4(G') tests
def _config_H(gens, k):
    """H(T) = log2 |sum of the subgroups gens[i], i in T| for the 128 subsets (seven variables)."""
    H = np.zeros(128, dtype=np.int64)
    for T in range(1, 128):
        rows = [gens[i] for i in range(7) if T >> i & 1 and len(gens[i])]
        H[T] = log_size(np.concatenate(rows), k) if rows else 0
    return H


def lemma_tests(n, rng, k, exponent=None):
    """Lemma 4(G') with coefficient 3 - 1/k on n configurations in Z_{2^e}^r (e = exponent, default k): near-J - I
    families (the seven-point configuration with random multipliers and perturbations) and random subgroups.  Returns
    the minimum of k (RHS - LHS) and the number of configurations."""
    e = k if exponent is None else exponent
    N = 1 << e
    best = None
    order = (A_, B_, C_, Y_, X_, W_, Z_)                            # (A1, A2, A3, B1, B2, B3, C) -> bits
    base = {A_: (1, 0, 0), B_: (0, 1, 0), C_: (0, 0, 1), W_: (1, 1, 0), X_: (1, 0, 1), Y_: (0, 1, 1), Z_: (1, 1, 1)}
    for t in range(n):
        r = int(rng.integers(3, 5))
        gens = [None] * 7
        if t % 2 == 0:
            for P in range(7):
                v = np.zeros(r, dtype=np.int64); v[:3] = base[P]
                rows = [(v * int(rng.integers(1, N)) if rng.random() < 0.3 else v) % N]
                if rng.random() < 0.25:
                    rows.append((rng.integers(0, N, size=r) << int(rng.integers(0, e))) % N)
                if rng.random() < 0.1:
                    rows = []
                gens[P] = np.array(rows, dtype=np.int64).reshape(len(rows), r)
        else:
            for P in range(7):
                c = int(rng.integers(0, 3))
                gens[P] = (rng.integers(0, N, size=(c, r)) << rng.integers(0, e, size=(c, 1))) % N
        H = _config_H(gens, e)
        val = lemma_g_prime(H, k)
        best = val if best is None else min(best, val)
    return best, n


# ----------------------------------------------------------------------------------------------- adversarial search
def folded(F):
    """a form on the 128 subsets folded onto the 63 masks of the visible parties (Z = the purifier, purity)."""
    G = np.zeros(64, dtype=np.int64)
    for i, c in enumerate(F):
        if 0 < i < 127:
            G[i if not i >> 6 & 1 else 127 ^ i] += c
    return G


def anneal(configs, seconds=60.0, seed=4401, mu=0.0, T0=3.0, form=None, singles=False, log=print):
    """Adversarial search (empirical): cc/s38_minform_zn.c minimises a f^(a) (or the given form, a function of a
    returning 128 coefficients) over CSS states of codes over Z_{2^a}; configs = [(a, rows, sizes)], one process each,
    run in parallel.  With mu > 0 the search is pulled towards m S^(2^a) (m = qudits of system A), or, with singles,
    only its single-system entropies towards their largest values n_P a.  Every reported code is rebuilt here, its
    entropies recomputed exactly, and passed through chain (all steps of the proof).  Returns one dict per
    configuration."""
    tmp = tempfile.mkdtemp(prefix='s38_')
    try:
        exe = s34.build('s38_minform_zn', tmp)
        procs = []
        for i, (k, r, sz) in enumerate(configs):
            F = fk_form(k) if form is None else form(k)
            G = folded(F)
            inp = ' '.join(str(int(G[m])) for m in range(1, 64))
            if mu > 0:
                if singles:
                    one = {1 << P: sz[P] * k for P in range(6)}
                    one[63] = sz[6] * k                                 # S(ABCWXY) = S(Z)
                    tg = [one.get(m, -1) for m in range(1, 64)]
                else:
                    e = seven_point_entropy(1 << k)
                    tg = [sz[0] * int(round(math.log2(e[m][0] / e[m][1]))) for m in range(1, 64)]
                inp += f' {mu} ' + ' '.join(str(x) for x in tg)
            args = [exe, str(k), str(r)] + [str(x) for x in sz] + [str(seconds), str(seed + 17 * i), str(T0)]
            pr = subprocess.Popen(args, stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True)
            pr.stdin.write(inp + '\n')
            pr.stdin.close()
            procs.append(pr)
        outs = []
        for pr in procs:
            outs.append(pr.stdout.read())
            pr.wait()
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    res = []
    for (k, r, sz), out in zip(configs, outs):
        t = out.split()
        minF, obj, restarts, steps = int(t[1]), int(t[3]), int(t[5]), int(t[7])
        c = out.splitlines()[1].split()
        rr, m = int(c[1]), int(c[2])
        Mx = np.array([int(x) for x in c[3:3 + rr * m]], dtype=np.int64).reshape(rr, m)
        owner = [P for P in range(7) for _ in range(sz[P])]
        st = css_state(Mx, owner, k)
        S = st.entropies()
        F = fk_form(k) if form is None else form(k)
        f = int(F @ S)
        vals, bad = chain(st, S)
        res.append(dict(k=k, rows=r, sizes=list(sz), min_form=minF, objective=obj, restarts=restarts, steps=steps,
                        form_recomputed=f, kfk=vals['kfk'], f12=vals['f12'], chain_failed=bad,
                        singles=[int(S[1 << P]) for P in range(7)]))
        log(f'  Z_{1 << k}, {r} rows, sizes {list(sz)}: min {minF} (recomputed {f}), objective {obj}, {restarts} restarts, '
            f'{steps} moves; at the minimiser a f^(a) = {vals["kfk"]}, f12 = {vals["f12"]}, S(P) = {res[-1]["singles"]}; '
            f'proof steps there: {"all hold" if not bad else bad}')
    return res


def lifts(k, base, s, forms, log=print):
    """Exhaustive (one qudit per system): every code B + 2^s P over Z_{2^k}, P a 0/1 matrix of the shape of B
    (cc/s38_minform_zn.c, lifts mode).  forms = {name: 128 coefficients}.  Each reported minimiser is rebuilt and
    recomputed here; the minimiser of the first form also passes through chain.  Returns {name: dict}."""
    B = np.asarray(base, dtype=np.int64) % (1 << k)
    r, m = B.shape
    assert m == 7
    names = list(forms)
    inp = ' '.join(' '.join(str(int(x)) for x in folded(forms[n])[1:64]) for n in names)
    inp += ' ' + ' '.join(str(int(x)) for x in B.ravel())
    tmp = tempfile.mkdtemp(prefix='s38_')
    try:
        exe = s34.build('s38_minform_zn', tmp)
        p = subprocess.run([exe, str(k), str(r)] + ['1'] * 7 + ['lifts', str(s), str(len(names))], input=inp + '\n',
                           capture_output=True, text=True)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    lines = p.stdout.splitlines()
    count = int(lines[0].split()[1])
    out = {}
    for n, ln in zip(names, lines[1:]):
        t = ln.split()
        mn, zeros = int(t[3]), int(t[5])
        rr, mm = int(t[7]), int(t[8])
        Mx = np.array([int(x) for x in t[9:9 + rr * mm]], dtype=np.int64).reshape(rr, mm)
        st = css_state(Mx, list(range(7)), k)
        S = st.entropies()
        rec = int(forms[n] @ S)
        res = dict(codes=count, min=mn, recomputed=rec, kfk_at_min=int(fk_form(k) @ S))
        if n == names[0]:
            vals, bad = chain(st, S)
            res.update(zeros=zeros, chain_failed=bad)
        out[n] = res
    log(f'  Z_{1 << k}, base {B.tolist()}, + 2^{s} P: {count} codes; ' +
        '; '.join(f'{n}: min {v["min"]} (recomputed {v["recomputed"]})' for n, v in out.items()) +
        f'; zeros of the first form {out[names[0]]["zeros"]}; proof steps at its minimiser: '
        f'{"all hold" if not out[names[0]]["chain_failed"] else out[names[0]]["chain_failed"]}')
    return out


# ----------------------------------------------------------------------------------------------- the checks
def _random_state(rng, t, k):
    kind = t % 6
    if kind == 0:
        sizes = rng.integers(1, 3, size=7)
        return 'graph', random_graph(rng, sizes, k, diagonal=bool(rng.random() < 0.3))
    if kind == 1:
        return 'structured', structured_css(rng, k)
    if kind == 2:
        sizes = rng.integers(1, 3, size=7)
        return 'css', random_css(rng, sizes, k)
    if kind == 3:
        sizes = rng.integers(0, 3, size=7)
        sizes[rng.integers(0, 7)] += 1
        return 'lagrangian', random_lagrangian(rng, sizes, k)
    if kind == 4:
        return 'seven-point', local_scramble(seven_point_perturbed(rng, k), rng)
    j = int(rng.integers(1, k + 1))
    st = seven_point_state(j, k)
    if rng.random() < 0.5:
        st = direct_sum(st, seven_point_state(int(rng.integers(1, k + 1)), k))
    return 'seven-point sums', local_scramble(st, rng)


def check(full_check=False, trials=60, seed=44, log=print):
    """Exact values (always) and, with full_check, the proof step by step on random states over Z_4 and Z_8, the
    state-vector check of the entropy formula, and the Lemma 4(G') tests."""
    ks = s34.known_states()
    out = {}
    out['nf'] = {k: str(Fraction(int(fk_form(k) @ ks['NF']), k)) for k in (1, 2, 3, 4)}
    out['steane'] = {k: str(Fraction(int(fk_form(k) @ ks['F']), k)) for k in (1, 2, 3, 4)}
    out['f1_is_f12'] = bool(np.array_equal(fk_form(1), s34.transfer_form()))
    sp = {}
    for k in (1, 2, 3, 4):
        for j in range(1, 6):
            q = fk_seven_point(1 << j, k)                           # 2^(a f^(a)) = 2^(a - j) expected
            sp[(k, j)] = q == Fraction(2) ** (k - j)
    out['seven_point_two_powers'] = all(sp.values())
    odd = {}
    for k in (1, 2, 3):
        for d in (3, 5, 6, 7, 9, 10, 12, 15):
            q = fk_seven_point(d, k)                                # 2^(a f^(a)) = 2^(a [2|d]) / d
            odd[(k, d)] = q == Fraction(2 ** (k * (d % 2 == 0)), d)
    out['seven_point_other_d'] = all(odd.values())
    st = {j: seven_point_state(j, 3) for j in (1, 2, 3)}
    out['seven_point_realised'] = all(st[j].is_lagrangian() and
                                      int(fk_form(3) @ st[j].entropies()) == 3 - j for j in (1, 2, 3))
    tight = {}
    for k in (1, 2, 3):                                             # the tight case: every bound attained
        vals, bad = chain(seven_point_state(k, k))
        tight[k] = (not bad and vals['kfk'] == 0 and vals['LA_slack'] == 0
                    and all(vals['delta'][P] == vals['U'][P] == k for P in (Z_, W_, X_, Y_)))
    out['tight'] = tight
    out['uz_sha'] = hashlib.sha256(np.asarray(uz_form(), dtype='<i8').tobytes()).hexdigest()
    if full_check:
        rng = np.random.default_rng(seed)
        res = {}
        for k in (2, 3):
            kinds, failed, mink, tight, minf12 = {}, 0, None, 0, None
            for t in range(trials):
                kind, stt = _random_state(rng, t, k)
                if not stt.is_lagrangian():
                    failed += 1
                    continue
                vals, bad = chain(stt)
                kinds[kind] = kinds.get(kind, 0) + 1
                failed += bool(bad)
                if bad:
                    log(f'  k = {k}, {kind}: {bad}')
                mink = vals['kfk'] if mink is None else min(mink, vals['kfk'])
                minf12 = vals['f12'] if minf12 is None else min(minf12, vals['f12'])
                tight += vals['kfk'] == 0
            res[k] = dict(kinds=kinds, failed=failed, min_kfk=mink, at_zero=tight, min_f12=minf12)
        out['chain'] = res
        out['state_vector_dev'] = max(state_vector_check(rng, 2, 6), state_vector_check(rng, 3, 2))
        lt = {}
        for k, e_ in ((1, 1), (2, 2), (3, 3), (1, 2), (2, 3)):
            lt[f'a={k}, exponent 2^{e_}'] = lemma_tests(150, rng, k, e_)[0]
        out['lemma_tests'] = lt
        quiet = lambda *_x: None
        lf = {}
        for k, B, sh, name in ((2, np.zeros((3, 7), dtype=np.int64), 0, 'Z_4, all 3 x 7 0/1 codes'),
                               (3, SEVEN, 1, 'Z_8, seven-point + 2P')):
            r = lifts(k, B, sh, weakened_forms(k), log=quiet)
            first = r['a f^(a)']
            lf[name] = dict(codes=first['codes'], min=first['min'], zeros=first['zeros'],
                            consistent=all(v['min'] == v['recomputed'] for v in r.values()) and not first['chain_failed'],
                            weakened={n: v['min'] for n, v in r.items() if n != 'a f^(a)'})
        out['lifts'] = lf
    log(f'  s35: {out}')
    return out


def check_ok(o, full_check=False):
    ok = (o['nf'] == {1: '-1', 2: '-1/2', 3: '-1/3', 4: '-1/4'} and o['steane'] == {1: '0', 2: '1/2', 3: '2/3', 4: '3/4'}
          and o['f1_is_f12'] and o['seven_point_two_powers'] and o['seven_point_other_d'] and o['seven_point_realised']
          and all(o['tight'].values()))
    if full_check:
        ok = ok and all(v['failed'] == 0 and v['min_kfk'] >= 0 for v in o['chain'].values()) and o['state_vector_dev'] < 1e-9
        lt = o['lemma_tests']
        ok = ok and lt['a=1, exponent 2^1'] >= 0 and lt['a=2, exponent 2^2'] >= 0 and lt['a=3, exponent 2^3'] >= 0 \
            and lt['a=1, exponent 2^2'] < 0 and lt['a=2, exponent 2^3'] < 0
        lf = o['lifts']
        ok = ok and all(v['min'] == 0 and v['consistent'] for v in lf.values())
        ok = ok and all(m < 0 for m in lf['Z_4, all 3 x 7 0/1 codes']['weakened'].values())
        ok = ok and lf['Z_8, seven-point + 2P']['weakened']['U_Z coefficient (1 - 1/a) - 1/(2a), times 2a'] < 0
    return ok


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--check', action='store_true')
    ap.add_argument('--full', action='store_true')
    ap.add_argument('--trials', type=int, default=60)
    ap.add_argument('--seed', type=int, default=44)
    ap.add_argument('--anneal', default=None,
                    help='semicolon-separated "a:rows:n0,...,n6", e.g. "2:3:1,1,1,1,1,1,1;3:4:1,1,1,1,1,1,1": cc/s38_minform_zn.c')
    ap.add_argument('--seconds', type=float, default=60.0)
    ap.add_argument('--pull', type=float, default=0.0, help='with --anneal: weight of a pull towards m S^(2^a)')
    ap.add_argument('--singles', action='store_true', help='with --pull: pull only S(P) towards n_P a')
    ap.add_argument('--weakened', type=str, default=None,
                    help='with --anneal: minimise a f12 + (a - 1 - c) U_Z instead (c = this value; sensitivity test)')
    a = ap.parse_args()
    if a.anneal:
        cf = []
        for item in a.anneal.split(';'):
            k_, r_, sz_ = item.split(':')
            cf.append((int(k_), int(r_), [int(x) for x in sz_.split(',')]))
        form = None
        if a.weakened is not None:
            c_ = Fraction(a.weakened)
            form = lambda k: (c_.denominator * (k * s34.transfer_form() + (k - 1) * uz_form())
                              - c_.numerator * uz_form())
        r = anneal(cf, seconds=a.seconds, seed=a.seed, mu=a.pull, form=form, singles=a.singles)
        ok = all(x['form_recomputed'] == x['min_form'] and not x['chain_failed'] for x in r)
        if a.weakened is None:
            ok = ok and all(x['min_form'] >= 0 for x in r)
            print('ALL OK' if ok else 'FAILURE')
        else:
            print('weakened form: minima', [x['min_form'] for x in r], '(negative = detected)',
                  'consistent' if ok else 'FAILURE')
    if a.check:
        o = check(full_check=a.full, trials=a.trials, seed=a.seed)
        print('ALL OK' if check_ok(o, a.full) else 'FAILURE')


if __name__ == '__main__':
    main()
