#!/usr/bin/env python3
"""s31: six parties already separate qubits from every odd local dimension (R36).

docs/R32_reduction_and_local_dimension.md, section 4.8 (Theorem 9).

Notation: seven elements A, B, C, W, X, Y, Z = 100, 010, 001, 110, 101, 011, 111 (Fano labelling, bit i <-> NAMES[i]),
Fano lines ABW, ACX, BCY, AYZ, BXZ, CWZ, WXY.  For subspaces U_P with rank function h:

  * F_i^L = U_i cap (U_j + U_k) for each line L = {i, j, k} through i, of dimension m_i^L = I(i; L - i);
  * Fhat_i = sum of the three F_i^L (the derived family), with ranks hhat.

Transfer (Lemmas 0, A-E of the note): every hhat(X) is bounded above and below by explicit linear forms in h.  Lemma 4
(LB)_2, valid in characteristic != 2 and for abelian groups of odd order, applied to the derived family and bounded
term by term (upper bounds where (LB)_2 has positive coefficients, lower bounds where negative) gives an integer form
f9 on seven variables with f9 >= 0 on LR^(q)_7 for every odd q and on Gamma^(d)_7 for every odd d.  The h-vector of
the seven-qubit simplex-code (Fano) state, i.e. the logical |0> of Steane's [[7,1,3]] code, is h_F = r_F + r_F^*, and
f9(h_F) = -3.  With Theorem 1(d) this gives
Stab^(2)_6 not contained in Stab^(Z_d)_6 for every odd d (Theorem 9).

Checks (all exact):
  --derive          recompute f9 from the definitions and compare with data/s31_six_party.npz
  --witness         simplex-code state vector (128 amplitudes): h_S = r_F + r_F^*, f9 values, (f9 + f9 o D)(h_F)
  --invariant       Fano-invariant self-dual form of f9 (and of the DFZ (65) variant)
  --lemmas N        Lemmas 0, A-E and the chain on N random arrangements per field, with the actual subspaces
  --stress N        f9 >= 0 on N random arrangements over GF(3), GF(5), GF(7); --exhaustive adds all
                    {0,1}-point configurations of GF(p)^3 (p = 2: sensitivity, p = 3, 5: validity)
  --face            the faces of the Shannon cone at h_F and h_NF are two-dimensional, cone(r, r^*) (exact ranks)

  python3 scripts/s31_six_party.py --derive --witness --invariant --lemmas 60 --stress 200
  python3 scripts/s31_six_party.py --derive --out data/s31_six_party.npz          # regenerate the stored form
"""
import argparse, itertools, os, sys, time
from fractions import Fraction
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..'))
from epr1kit import core, stabcert
import s27_cube_states as s27
import s28_generalized_cubes as g28

NAMES = 'ABCWXYZ'
LINES = ('ABW', 'ACX', 'BCY', 'AYZ', 'BXZ', 'CWZ', 'WXY')
FULL = 127
LEMMA_TO_FANO = (0, 1, 2, 5, 4, 3, 6)        # Lemma 4 bits A_1, A_2, A_3, B_1, B_2, B_3, C -> A, B, C, Y, X, W, Z
FORM = os.path.join(core.DATA, 's31_six_party.npz')


def mask(s):
    return sum(1 << NAMES.index(c) for c in s)


def label(X):
    return ''.join(NAMES[i] for i in range(7) if X >> i & 1)


LINE_MASKS = tuple(mask(L) for L in LINES)


# ------------------------------------------------------------------ exact linear forms in h (index = mask)
def _zero():
    return [Fraction(0)] * 128


def _add(*vs, c=None):
    out = _zero()
    for k, v in enumerate(vs):
        w = 1 if c is None else c[k]
        for X in range(128):
            if v[X]:
                out[X] += w * v[X]
    return out


def _H(X):
    v = _zero()
    if X:
        v[X] = Fraction(1)
    return v


def _I(a, b):
    """I(a; b) = h(a) + h(b) - h(a u b) for disjoint masks a, b."""
    return _add(_H(a), _H(b), _H(a | b), c=(1, 1, -1))


def lines_through(i):
    return [L for L in LINE_MASKS if L >> i & 1]


def m_form(i, L):
    """m_i^L = I(i; L - i) = dim F_i^L (Lemma 0)."""
    return _I(1 << i, L & ~(1 << i))


def u_form(i):
    """Lemma B: hhat(i) <= u_i = (1/3) [2 sum_{a<b} I(i; Q_ab) - sum_L m_i^L], Q_ab = (L_a u L_b) - i."""
    Ls = lines_through(i)
    parts, coefs = [], []
    for La, Lb in itertools.combinations(Ls, 2):
        parts.append(_I(1 << i, (La | Lb) & ~(1 << i))); coefs.append(Fraction(2, 3))
    for L in Ls:
        parts.append(m_form(i, L)); coefs.append(Fraction(-1, 3))
    return _add(*parts, c=coefs)


def l_form(i):
    """Lemma C: hhat(i) >= l_i = (1/3) sum_L m_i^L."""
    Ls = lines_through(i)
    return _add(*[m_form(i, L) for L in Ls], c=[Fraction(1, 3)] * len(Ls))


def lower_form(X):
    """Lemma C: hhat(X) >= h(X) - sum_{i in X} (h(i) - l_i)."""
    parts, coefs = [_H(X)], [1]
    for i in range(7):
        if X >> i & 1:
            parts += [_H(1 << i), l_form(i)]; coefs += [-1, 1]
    return _add(*parts, c=coefs)


def upper_form(order, credits):
    """Lemmas B, D, E: add the points of `order` one by one; a credit (i, L) subtracts m_i^L when the other two points
    of the line L were added before i (Lemma A: F_i^L lies in Fhat_j + Fhat_k).  hhat(X) <= hhat(set(order)) for X in
    order (monotonicity).  Returns sum_{i in order} u_i - sum_credits m_i^L."""
    pos = {p: k for k, p in enumerate(order)}
    for i, L in credits:
        others = [j for j in range(7) if L >> j & 1 and j != i]
        if not (L >> i & 1 and all(j in pos and pos[j] < pos[i] for j in others)):
            raise ValueError(f'invalid credit ({NAMES[i]}, {label(L)}) for the order {[NAMES[p] for p in order]}')
    parts = [u_form(p) for p in order] + [m_form(i, L) for i, L in credits]
    return _add(*parts, c=[1] * len(order) + [-1] * len(credits))


def recipes_lb2():
    """the upper-bound recipes used for (LB)_2 (positive-coefficient sets: points, lines, ABCZ)."""
    R = {}
    for i in range(7):
        R[1 << i] = ((i,), ())
    for L in LINES:                          # Lemma D: credit to the first point of L (order A, B, C, W, X, Y, Z)
        i, j, k = (NAMES.index(c) for c in L)
        R[mask(L)] = ((j, k, i), ((i, mask(L)),))
    A, B, C, Y, Z = (NAMES.index(c) for c in 'ABCYZ')   # Lemma E: ABCZ <= ABCYZ, credits Y (BCY) and Z (AYZ)
    R[mask('ABCZ')] = ((A, B, C, Y, Z), ((Y, mask('BCY')), (Z, mask('AYZ'))))
    return R


def to_fano(v):
    w = np.zeros(128, dtype=np.int64)
    for X in range(1, 128):
        if v[X]:
            w[sum(1 << LEMMA_TO_FANO[i] for i in range(7) if X >> i & 1)] += int(v[X])
    return w


def lb2_fano():
    """Lemma 4 (LB)_2 in the Fano labelling (A_i = e_i, B_i = c - e_i, C = c)."""
    return to_fano(g28.lemma_forms(2)[1])


def transfer(base, recipes):
    """3 * (sum_{c_X > 0} c_X UB_X + sum_{c_X < 0} c_X LB_X) as an integer vector (asserts integrality)."""
    tot = _zero()
    for X in range(1, 128):
        c = int(base[X])
        if c > 0:
            if X not in recipes:
                raise ValueError(f'no upper-bound recipe for {label(X)}')
            tot = _add(tot, upper_form(*recipes[X]), c=(1, c))
        elif c < 0:
            tot = _add(tot, lower_form(X), c=(1, c))
    out = np.zeros(128, dtype=np.int64)
    for X in range(128):
        y = 3 * tot[X]
        if y.denominator != 1:
            raise ValueError('non-integral coefficient')
        out[X] = int(y)
    return out


def f9():
    return transfer(lb2_fano(), recipes_lb2())


# ------------------------------------------------------------------ generic recipes (for other base forms)
def _credit_order(Xp):
    """an ordering of the points of Xp with the maximal number of credits (DP over subsets)."""
    pts = [i for i in range(7) if Xp >> i & 1]
    best = {0: (0, ())}
    for S in sorted(range(1, 1 << 7), key=lambda s: bin(s).count('1')):
        if S & ~Xp:
            continue
        cand = None
        for i in pts:
            if S >> i & 1:
                prev = S ^ (1 << i)
                if prev not in best:
                    continue
                cr, seq = best[prev]
                Ls = [L for L in lines_through(i) if (L & ~(1 << i)) & ~prev == 0]
                item = (cr + (1 if Ls else 0), seq + ((i, Ls[0] if Ls else None),))
                if cand is None or item[0] > cand[0]:
                    cand = item
        best[S] = cand
    cr, seq = best[Xp]
    return cr, tuple(i for i, _L in seq), tuple((i, L) for i, L in seq if L is not None)


def auto_recipe(X, rF):
    """smallest superset X' of X with an ordering whose bound is tight at h_F: |X'| - credits = r_F(X)."""
    best = None
    for Xp in range(1, 128):
        if Xp & X != X:
            continue
        cr, order, credits = _credit_order(Xp)
        s = bin(Xp).count('1')
        if s - cr == rF[X] and (best is None or s < best[0]):
            best = (s, order, credits)
    if best is None:
        raise ValueError(f'no tight recipe for {label(X)}')
    return best[1], best[2]


def transfer_auto(base):
    rF = s27.dfz_configuration(2)
    return transfer(base, {X: auto_recipe(X, rF) for X in range(1, 128) if base[X] > 0})


# ------------------------------------------------------------------ witness
def h_star(h):
    """Kaced's dual: h*(X) = sum_{P in X} h(P) + h(E - X) - h(E)."""
    h = np.asarray(h)
    return np.array([sum(int(h[1 << P]) for P in range(7) if X >> P & 1) + int(h[FULL ^ X]) - int(h[FULL]) if X else 0
                     for X in range(128)], dtype=np.int64)


def simplex_state():
    """the CSS state of the [7,3] simplex code (generator columns = A..Z), qubit i <-> bit i."""
    cols = [[1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 0], [1, 0, 1], [0, 1, 1], [1, 1, 1]]
    G = np.array(cols).T
    psi = np.zeros(128)
    for m in itertools.product((0, 1), repeat=3):
        c = (np.array(m) @ G) % 2
        psi[sum(int(c[i]) << i for i in range(7))] += 1
    return psi / np.linalg.norm(psi)


def entropies(psi, n=7):
    t = psi.reshape([2] * n)
    S = np.zeros(1 << n)
    for X in range(1, 1 << n):
        A = [i for i in range(n) if X >> i & 1]; B = [i for i in range(n) if not X >> i & 1]
        M = np.transpose(t, [n - 1 - i for i in A] + [n - 1 - i for i in B]).reshape(2 ** len(A), 2 ** len(B))
        s = np.linalg.svd(M, compute_uv=False) ** 2
        s = s[s > 1e-12]
        S[X] = float(-(s * np.log2(s)).sum())
    return S


def witness(f=None):
    f = f9() if f is None else f
    S = entropies(simplex_state())
    hS = np.array([S[X] + sum(S[1 << i] for i in range(7) if X >> i & 1) for X in range(128)])
    rF, rNF = s27.dfz_configuration(2), s27.dfz_configuration(3)
    hF, hNF = rF + h_star(rF), rNF + h_star(rNF)
    from s30_dimension_profile import selfdual_symmetrise
    g = selfdual_symmetrise(f)
    return {'dev': float(np.abs(hS - hF).max()), 'purity': float(max(abs(S[X] - S[FULL ^ X]) for X in range(1, FULL))),
            'f(h_S)': float(f @ hS), 'f(h_F)': int(f @ hF), 'f(r_F)': int(f @ rF), 'f(h_NF)': int(f @ hNF),
            'f(r_NF)': int(f @ rNF), '(f + f o D)(h_F)': int(g @ hF), 'self-dual h_F': bool(np.array_equal(h_star(hF), hF)),
            'terms': int(np.count_nonzero(f[1:])), 'max |coef|': int(np.abs(f).max())}


# ------------------------------------------------------------------ Fano-invariant self-dual form
def orbit(X):
    k = bin(X).count('1')
    if k == 3:
        return 'L' if X in LINE_MASKS else 'N'
    if k == 4:
        return 'M' if (FULL ^ X) in LINE_MASKS else 'K'
    return {1: 'P', 2: 'Q', 5: 'F', 6: 'S', 7: 'E'}[k]


def invariant_form(f):
    """f restricted to Fano-invariant self-dual vectors, in the coordinates (v_P, alpha, nu, mu):
    alpha = 2 v_P - v_Q, mu = v_P + v_Q - v_L, nu = v_P + v_Q - v_N (v_O = value on the orbit O).  Self-duality gives
    v_E = 7 v_P / 2, v_S = v_E, v_F = v_Q - 2 v_P + v_E, v_M = v_L - 3 v_P + v_E, v_K = v_N - 3 v_P + v_E.
    Returns the coefficients of (v_P, alpha, nu, mu) as Fractions; the value at h_F is 2 c_P + c_mu."""
    h = Fraction(7, 2)
    val = {'P': (1, 0, 0, 0), 'Q': (2, -1, 0, 0), 'L': (3, -1, 0, -1), 'N': (3, -1, -1, 0), 'E': (h, 0, 0, 0),
           'S': (h, 0, 0, 0), 'F': (h, -1, 0, 0), 'M': (h, -1, 0, -1), 'K': (h, -1, -1, 0)}
    out = [Fraction(0)] * 4
    for X in range(1, 128):
        if f[X]:
            for k, c in enumerate(val[orbit(X)]):
                out[k] += int(f[X]) * Fraction(c)
    return tuple(out)


# ------------------------------------------------------------------ lemma-by-lemma checks with actual subspaces
def _rref(M, q):
    M = np.array(M, dtype=np.int64) % q
    if M.ndim != 2 or M.shape[0] == 0:
        return M.reshape(0, M.shape[-1] if M.ndim == 2 else 0)
    r = 0
    for c in range(M.shape[1]):
        piv = next((i for i in range(r, M.shape[0]) if M[i, c]), None)
        if piv is None:
            continue
        M[[r, piv]] = M[[piv, r]]
        M[r] = (M[r] * pow(int(M[r, c]), q - 2, q)) % q
        for i in range(M.shape[0]):
            if i != r and M[i, c]:
                M[i] = (M[i] - M[i, c] * M[r]) % q
        r += 1
        if r == M.shape[0]:
            break
    return M[:r]


def _span(vs, q, n):
    vs = [v for v in vs if len(v)]
    return _rref(np.vstack(vs), q) if vs else np.zeros((0, n), dtype=np.int64)


def _intersect(U, W, q, n):
    """U cap W: left kernel of [U; W] via row reduction of [U W | I]."""
    if not len(U) or not len(W):
        return np.zeros((0, n), dtype=np.int64)
    M = np.vstack([U, W]); a = M.shape[0]
    R = np.hstack([M % q, np.eye(a, dtype=np.int64)])
    r = 0
    for c in range(n):
        piv = next((i for i in range(r, a) if R[i, c]), None)
        if piv is None:
            continue
        R[[r, piv]] = R[[piv, r]]
        R[r] = (R[r] * pow(int(R[r, c]), q - 2, q)) % q
        for i in range(a):
            if i != r and R[i, c]:
                R[i] = (R[i] - R[i, c] * R[r]) % q
        r += 1
    K = R[r:, n:]
    return _span([(K[:, :len(U)] @ U) % q], q, n) if len(K) else np.zeros((0, n), dtype=np.int64)


def _random_family(q, n, rng, kind):
    if kind == 0:                                         # dense random, element dimensions 0-3
        return [_rref(rng.integers(0, q, size=(int(rng.integers(0, 4)), n)), q) for _ in range(7)]
    if kind == 1:                                         # sparse random
        return [_rref(rng.integers(1, q, size=(d, n)) * (rng.random((d, n)) < 0.35), q)
                for d in rng.integers(0, 4, size=7)]
    base = rng.integers(0, q, size=(3, n))                # point-like: scaled Fano/non-Fano points plus noise
    pts = [(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0), (1, 0, 1), (0, 1, 1), (1, 1, 1)]
    out = []
    for p in pts:
        v = ((np.array(p) * rng.integers(1, q, size=3)) @ base) % q
        out.append(_rref(np.vstack([v[None, :], rng.integers(0, q, size=(int(rng.integers(0, 2)), n))]), q))
    return out


def _value(v, h):
    return sum(v[X] * h[X] for X in range(1, 128) if v[X])


def lemma_checks(q, trials, rng, f=None):
    """minimum slack of every lemma bound (and of the chain) over `trials` random families over GF(q)."""
    f = f9() if f is None else f
    base = lb2_fano(); R = recipes_lb2()
    U_ = {i: u_form(i) for i in range(7)}; L_ = {i: l_form(i) for i in range(7)}
    UB = {X: upper_form(*R[X]) for X in range(1, 128) if base[X] > 0}
    LBf = {X: lower_form(X) for X in range(1, 128) if base[X] < 0}
    slack = {}

    def rec(key, s):
        slack[key] = min(slack.get(key, s), s)
    for t in range(trials):
        n = int(rng.integers(3, 8))
        U = _random_family(q, n, rng, t % 3)
        h = [0] + [len(_span([U[i] for i in range(7) if X >> i & 1], q, n)) for X in range(1, 128)]
        Fs = {}
        for L in LINE_MASKS:
            for i in range(7):
                if L >> i & 1:
                    j, k = [x for x in range(7) if L >> x & 1 and x != i]
                    Fs[i, L] = _intersect(U[i], _span([U[j], U[k]], q, n), q, n)
                    rec('Lemma 0: dim F_i^L = m_i^L', 0 if len(Fs[i, L]) == _value(m_form(i, L), h) else -1)
        for L in LINE_MASKS:
            for i in range(7):
                if L >> i & 1:
                    j, k = [x for x in range(7) if L >> x & 1 and x != i]
                    S_jk = _span([Fs[j, L], Fs[k, L]], q, n)
                    rec('Lemma A: F_i^L in F_j^L + F_k^L', 0 if len(_span([S_jk, Fs[i, L]], q, n)) == len(S_jk) else -1)
        Fh = [_span([Fs[i, L] for L in lines_through(i)], q, n) for i in range(7)]
        hh = [0] + [len(_span([Fh[i] for i in range(7) if X >> i & 1], q, n)) for X in range(1, 128)]
        for i in range(7):
            rec('Lemma B: hhat(i) <= u_i', _value(U_[i], h) - hh[1 << i])
            rec('Lemma C: hhat(i) >= l_i', hh[1 << i] - _value(L_[i], h))
        for X, v in UB.items():
            key = 'Lemma B (points)' if bin(X).count('1') == 1 else ('Lemma D (lines)' if X in LINE_MASKS else 'Lemma E (ABCZ)')
            rec(key + ': upper bound', _value(v, h) - hh[X])
        for X, v in LBf.items():
            rec('Lemma C: lower bound on negative sets', hh[X] - _value(v, h))
        rec('(LB)_2 on the derived family', int(base @ np.array(hh)))
        rec('f9(h)', int(f @ np.array(h)))
    return slack


LARGE_PRIME = 1_000_003


def face_check():
    """exact ranks of the Shannon elemental inequalities tight at r, r^* and r + r^* for the Fano (GF(2)) and non-Fano
    (GF(3)) matroids.  Upper bound: the tight rows annihilate v (and, at r + r^*, the independent pair r, r^*), so the
    rank is at most 126 (resp. 125).  Lower bound: the rank over GF(1000003) is a lower bound for the rational rank.
    Rank 126 = extreme ray of Gamma_7; rank 125 at r + r^* = a two-dimensional face, which then is cone(r, r^*).
    Returns {name: (number of tight rows, certified rational rank or None, minimum elemental value)}."""
    from s30_dimension_profile import elemental7
    E = elemental7()
    out = {}
    for p, tag in ((2, 'F'), (3, 'NF')):
        r = s27.dfz_configuration(p); rd = h_star(r)
        indep = stabcert.gfp_rank([r[1:].tolist(), rd[1:].tolist()], LARGE_PRIME) == 2
        for name, v, kernel in ((f'r_{tag}', r, [r]), (f'r_{tag}*', rd, [rd]), (f'h_{tag}', r + rd, [r, rd])):
            T = E[(E @ v[1:]) == 0]
            upper = 127 - len(kernel) if all(not (T @ w[1:]).any() for w in kernel) and indep else None
            lower = stabcert.gfp_rank(T.tolist(), LARGE_PRIME)
            out[name] = (int(len(T)), lower if upper == lower else None, int((E @ v[1:]).min()))
    return out


def stress(trials, rng, f=None, fields=(3, 5, 7)):
    f = f9() if f is None else f
    return {p: s27.random_arrangement_min(f, p, trials, rng) for p in fields}


def exhaustive(f=None):
    f = f9() if f is None else f
    return {p: g28.exhaustive01_min(f, 7, 3, p)[0] for p in (2, 3, 5)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--derive', action='store_true')
    ap.add_argument('--out', default=None)
    ap.add_argument('--witness', action='store_true')
    ap.add_argument('--invariant', action='store_true')
    ap.add_argument('--lemmas', type=int, default=0)
    ap.add_argument('--stress', type=int, default=0)
    ap.add_argument('--exhaustive', action='store_true')
    ap.add_argument('--face', action='store_true', help='Shannon faces of h_F and h_NF (exact ranks)')
    ap.add_argument('--seed', type=int, default=36)
    a = ap.parse_args()
    t0 = time.time()
    rng = np.random.default_rng(a.seed)
    f = f9()
    if a.derive:
        print(f'f9: {np.count_nonzero(f[1:])} nonzero coefficients, max |c| = {np.abs(f).max()}', flush=True)
        print('  ' + ', '.join(f'{label(X)}:{int(f[X])}' for X in range(1, 128) if f[X]), flush=True)
        if a.out:
            np.savez_compressed(a.out, f9=f, h_F=s27.dfz_configuration(2) + h_star(s27.dfz_configuration(2)))
            print('wrote', a.out)
        elif os.path.exists(FORM):
            print(f'  matches the stored {os.path.relpath(FORM)}: {np.array_equal(np.load(FORM)["f9"], f)}', flush=True)
    if a.witness:
        for k, v in witness(f).items():
            print(f'witness: {k} = {v}', flush=True)
    if a.invariant:
        names = ('v_P', 'alpha', 'nu', 'mu')
        for lab, g in (('f9 [(LB)_2]', f), ('DFZ (65) variant', transfer_auto(s27.load_dfz_char()['dfz65_odd']))):
            c = invariant_form(g)
            print(f'invariant: {lab}: (1/3) f = ' + ' + '.join(f'{x / 3}*{n_}' for x, n_ in zip(c, names))
                  + f'  (value at h_F: {(2 * c[0] + c[3]) / 3})', flush=True)
    if a.lemmas:
        for q in (2, 3, 5):
            s = lemma_checks(q, a.lemmas, rng, f)
            print(f'lemmas over GF({q}), {a.lemmas} families: ' + '; '.join(f'{k} {v}' for k, v in sorted(s.items()))
                  + f'  ({time.time() - t0:.0f}s)', flush=True)
    if a.stress:
        print(f'stress: min f9 over {a.stress} random arrangements: {stress(a.stress, rng, f)}  ({time.time() - t0:.0f}s)',
              flush=True)
    if a.exhaustive:
        print(f'exhaustive: min f9 over all {{0,1}}-point configurations of GF(p)^3: {exhaustive(f)} '
              f'(expect -3 for p = 2, >= 0 for p = 3, 5)  ({time.time() - t0:.0f}s)', flush=True)
    if a.face:
        for k, (nt, rk, mn) in face_check().items():
            print(f'face: {k}: {nt} tight elemental inequalities, exact rank {rk} (face dimension '
                  f'{127 - rk if rk is not None else "?"}), min elemental {mn}', flush=True)


if __name__ == '__main__':
    main()
