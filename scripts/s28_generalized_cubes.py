#!/usr/bin/env python3
"""s28: stabilizer cones of any two distinct primes are incomparable (R33).

Generalised cube states.  For t >= 2 put M = t + 1, c = (1, ..., 1) in F^M and take the 2M + 2 = 2t + 4 points

    e_1, ..., e_M  (A_1..A_M),   c - e_1, ..., c - e_M  (B_1..B_{t+1}),   c  (C),   0  (the purifier),

one qudit per point, and the code of affine functions on them (generator rows 1, x_1, ..., x_M).  For t = 2 the
points are the whole cube {0,1}^3 and these are the cube states of s27.  Facts checked here, exactly:

  * the code is monomially self-dual over F whenever char F does not divide t - 1, with explicit weights
    d = (-1 on the e_i, +1 on the c - e_i, -(M - 2) on c, M - 2 on 0):  G diag(d) G^T = 0.  Then the matroid is
    identically self-dual, S(X) = 2 r(X) - |X| and h_S = 2r;
  * contracting the purifier (the origin) leaves exactly the configuration A_i = <e_i>, B_i = <c - e_i>, C = <c>
    (the counterexample of Pena-Sarria, arXiv:1905.00003v3, Theorem 3 / Example 6 with M(n, t) = t + 1);
  * Lemma 4 (docs/R32_reduction_and_local_dimension.md, proved there; `lemma_forms`):
      LB, valid over every field whose characteristic does not divide t, and
      LA, valid over every field whose characteristic divides t,
    take the value -1 on that configuration over GF(p) with p | t (LB) resp. over GF(q) with q not dividing t (LA);
  * pulled back through "contract the purifier" and h_S, LA and LB become (2t+3)-party inequalities on S with value
    -2 on the corresponding witness.

Hence, for primes p < q and t = p: the state over GF(p) is outside Stab^(q) (LB is valid for char != p) and the
state over GF(q) is outside Stab^(p) (LA is valid for char p); both codes are self-dual because p and q do not
divide p - 1.  So Stab^(p)_n and Stab^(q)_n are incomparable for every n >= 2p + 3.

Cross-check: Pena-Sarria's own Example 6 inequalities (`ps_forms`, typed from the rendered page 4 of
arXiv:1905.00003v3) give the same -1 / -2 values.  Sanity tests (not proofs): no violation of a valid form on
random subspace arrangements, on all configurations of points spanned by {0,1}-vectors of GF(q)^3 for t = 2
(--exhaustive01, which is sensitive: the wrong-characteristic form is violated there), and on all subspace
arrangements of GF(q)^2 (--exhaustive2; insensitive to characteristic, kept as a transcription check).

  python3 scripts/s28_generalized_cubes.py --t 2 --fields 2 3 5 --exhaustive01      # Lemma 4 on all {0,1}-configurations
  python3 scripts/s28_generalized_cubes.py --t 3 --fields 3 5 7 [--random 300] [--exhaustive2] [--exhaustive01] [--out f.npz]
  python3 scripts/s28_generalized_cubes.py --t 5 --fields 5 7 [11]
  python3 scripts/s28_generalized_cubes.py --dfz1401     # the misprint in DFZ arXiv:1401.2507, Theorem 3.1
"""
import argparse, itertools, os, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..'))
from epr1kit import stabcert
from s27_cube_states import rank_function, dual, contract


# ------------------------------------------------------------------ Pena-Sarria, Example 6
def ps_forms(t, M):
    """Example 6 of arXiv:1905.00003v3 for given t and M = M(n, t) (n = M + t + 2 variables).
    Variables (bits): A_1..A_M = 0..M-1, B_1..B_{t+1} = M..M+t, C = M+t+1.
    Returns (a, b, n): integer vectors over the 2^n masks with  a . H >= 0  [valid if char divides t] and
    b . H >= 0 [valid if char does not divide t; the printed inequality multiplied by M].

    (a)  H(B_[t+1], A_[M]-[t+1]) + (t+2)(M-t-1) H(C)  <=  (M-1) I(A_[M]; C) + (t+2) sum_{i=t+2}^{M} H(A_i)
         + [(t+2)(M-t) - 1] ( H(C | A_[M]) + sum_{i=1}^{M} I(A_[M]-i; C) )
         + sum_{i=1}^{t+1} ( H(B_i | A_[M]-i) + H(B_i | A_i, C) + I(A_[i]; A_[t+1]-[i]) + I(A_[i-1]; A_i) )
    (b)  H(C)  <=  (1/M) H(B_[t+1], A_[M]-[t+1]) + H(C | A_[M]) + sum_{i=1}^{M} I(A_[M]-i; C)
         + sum_{i=2}^{t+1} I(A_[i-1]; A_i) + sum_{i=1}^{t+1} ( H(C | A_i, B_i) + H(B_i | A_[M]-i) + I(A_[i]; A_[M]-[i]) )
    """
    n = M + t + 2
    A = lambda i: 1 << (i - 1)
    B = lambda i: 1 << (M + i - 1)
    Cb = 1 << (M + t + 1)
    Aset = lambda idx: sum(A(i) for i in idx)
    allA = Aset(range(1, M + 1))
    Bt = sum(B(i) for i in range(1, t + 2))
    Arest = Aset(range(t + 2, M + 1))

    def H(c, X, k):
        if X:
            c[X] += k

    def Hc(c, X, Y, k):
        H(c, X | Y, k); H(c, Y, -k)

    def I(c, X, Y, k):
        if X and Y:
            H(c, X, k); H(c, Y, k); H(c, X | Y, -k)
    a = np.zeros(1 << n, dtype=np.int64)
    H(a, Bt | Arest, -1); H(a, Cb, -(t + 2) * (M - t - 1))
    I(a, allA, Cb, M - 1)
    for i in range(t + 2, M + 1):
        H(a, A(i), t + 2)
    K = (t + 2) * (M - t) - 1
    Hc(a, Cb, allA, K)
    for i in range(1, M + 1):
        I(a, allA ^ A(i), Cb, K)
    for i in range(1, t + 2):
        Hc(a, B(i), allA ^ A(i), 1); Hc(a, B(i), A(i) | Cb, 1)
        I(a, Aset(range(1, i + 1)), Aset(range(i + 1, t + 2)), 1); I(a, Aset(range(1, i)), A(i), 1)
    b = np.zeros(1 << n, dtype=np.int64)
    H(b, Cb, -M); H(b, Bt | Arest, 1)
    Hc(b, Cb, allA, M)
    for i in range(1, M + 1):
        I(b, allA ^ A(i), Cb, M)
    for i in range(2, t + 2):
        I(b, Aset(range(1, i)), A(i), M)
    for i in range(1, t + 2):
        Hc(b, Cb, A(i) | B(i), M); Hc(b, B(i), allA ^ A(i), M); I(b, Aset(range(1, i + 1)), Aset(range(i + 1, M + 1)), M)
    return a, b, n


# ------------------------------------------------------------------ Lemma 4 (our own inequalities)
def lemma_forms(t):
    """Lemma 4 of docs/R32_reduction_and_local_dimension.md (proved there) for M = t + 1.
    Variables (bits): A_1..A_M = 0..M-1, B_1..B_M = M..2M-1, C = 2M (the layout of ps_forms with M = t + 1).
    Error terms: alpha = sum_k H(A_k) - H(A_[M]); gamma = H(C | A_[M]); delta_i = I(C; A_{!=i}); eps_i = H(C | A_i, B_i);
    zeta_i = H(B_i | A_{!=i}); eta_i = H(B_i | A_i, C).  Returns (LA, LB, n), integer vectors with . H >= 0:
      LB [char does not divide t]:  H(B_[M]) + M gamma + M^2 alpha + M sum_i (delta_i + eps_i + zeta_i) - M H(C)
      LA [char divides t]:          (M-1) H(C) + M alpha + M gamma + sum_i (zeta_i + eta_i) - H(B_[M])"""
    M = t + 1; n = 2 * M + 1
    A = [1 << i for i in range(M)]; B = [1 << (M + i) for i in range(M)]; C = 1 << (2 * M)
    allA = sum(A)

    def H(v, X, k):
        if X:
            v[X] += k

    def alpha(v, k):
        for a in A:
            H(v, a, k)
        H(v, allA, -k)

    def gamma(v, k):
        H(v, C | allA, k); H(v, allA, -k)

    def delta(v, i, k):
        Ao = allA ^ A[i]; H(v, C, k); H(v, Ao, k); H(v, C | Ao, -k)

    def eps(v, i, k):
        H(v, C | A[i] | B[i], k); H(v, A[i] | B[i], -k)

    def zeta(v, i, k):
        Ao = allA ^ A[i]; H(v, B[i] | Ao, k); H(v, Ao, -k)

    def eta(v, i, k):
        H(v, B[i] | A[i] | C, k); H(v, A[i] | C, -k)
    LB = np.zeros(1 << n, dtype=np.int64)
    H(LB, sum(B), 1); gamma(LB, M); alpha(LB, M * M); H(LB, C, -M)
    for i in range(M):
        delta(LB, i, M); eps(LB, i, M); zeta(LB, i, M)
    LA = np.zeros(1 << n, dtype=np.int64)
    H(LA, C, M - 1); alpha(LA, M); gamma(LA, M); H(LA, sum(B), -1)
    for i in range(M):
        zeta(LA, i, 1); eta(LA, i, 1)
    return LA, LB, n


# ------------------------------------------------------------------ generalised cubes
def configuration(t):
    """M x (2t+3) matrix of the Pena-Sarria counterexample vectors (M = t + 1): A_i = e_i, B_i = c - e_i, C = c."""
    M = t + 1
    I = np.eye(M, dtype=np.int64); c = np.ones(M, dtype=np.int64)
    return np.array([I[:, i] for i in range(M)] + [c - I[:, i] for i in range(M)] + [c]).T


def generator(t):
    """(M+1) x (2t+4) generator of the affine code on the configuration points plus the origin (last column)."""
    V = configuration(t)
    pts = np.hstack([V, np.zeros((V.shape[0], 1), dtype=np.int64)])
    return np.vstack([np.ones((1, pts.shape[1]), dtype=np.int64), pts])


def self_dual_weights(t):
    M = t + 1
    return np.array([-1] * M + [1] * M + [-(M - 2), M - 2], dtype=np.int64)


def weights_ok(t, p):
    """exact check of G diag(d) G^T = 0 over GF(p) with all weights nonzero mod p."""
    G, d = generator(t), self_dual_weights(t)
    return bool(np.all(d % p != 0) and not ((G * d) @ G.T % p).any())


def pullback(f, n):
    """(n)-party inequality on S (2^n - 1 coordinates, purifier = party n) equivalent to f applied to the
    contraction of h_S by the purifier."""
    N = n + 1; FULL = (1 << N) - 1; PUR = n
    coef = np.zeros(1 << N, dtype=np.int64)

    def add_h(X, c):
        coef[X] += c
        for P in range(N):
            if X >> P & 1:
                coef[1 << P] += c
    for T in range(1, 1 << n):
        if f[T]:
            add_h(T | 1 << PUR, int(f[T])); add_h(1 << PUR, -int(f[T]))
    fold = np.zeros(1 << n, dtype=np.int64)
    for X in range(1, FULL):
        Xr = X if not X >> PUR & 1 else FULL ^ X
        fold[Xr] += coef[X]
    return fold[1:]


def witness(t, p):
    """rank function r (all 2^(2t+4) masks), entropy S (visible masks, purity built in) and h_S of the
    generalised cube state over GF(p)."""
    G = generator(t) % p
    N = G.shape[1]; FULL = (1 << N) - 1
    r = rank_function(G, p)
    S_all = np.array([r[X] + r[FULL ^ X] - r[FULL] if 0 < X < FULL else 0 for X in range(1 << N)], dtype=np.int64)
    h = np.array([S_all[X] + sum(S_all[1 << i] for i in range(N) if X >> i & 1) for X in range(1 << N)], dtype=np.int64)
    return r, S_all, h


def analyse(t, fields):
    """witness facts for each field; LA/LB = Lemma 4 (primary), a/b = Pena-Sarria Example 6 (cross-check)."""
    M = t + 1
    a, b, n = ps_forms(t, M)
    LA, LB, n2 = lemma_forms(t)
    assert n2 == n
    FA, FB, Fa, Fb = pullback(LA, n), pullback(LB, n), pullback(a, n), pullback(b, n)
    out = {'n_visible': n, 'FA': FA, 'FB': FB, 'Fa': Fa, 'Fb': Fb, 'fields': {}}
    for p in fields:
        r, S_all, h = witness(t, p)
        N = n + 1
        sd = bool(np.array_equal(dual(r, N), r))
        cfg = rank_function(configuration(t) % p, p)
        con = bool(np.array_equal(contract(r, N, n), cfg))
        S = S_all[1:1 << n]
        out['fields'][p] = dict(self_dual=sd, weights=weights_ok(t, p), contraction=con, h2r=bool(np.array_equal(h, 2 * r)),
                                cfg_LA=int(LA @ cfg), cfg_LB=int(LB @ cfg), cfg_a=int(a @ cfg), cfg_b=int(b @ cfg),
                                FA=int(FA @ S), FB=int(FB @ S), Fa=int(Fa @ S), Fb=int(Fb @ S), S=S)
    return out


def exhaustive01_min(f, n, M, p, limit=2_200_000, rng=None):
    """minimum of f over configurations in which every variable is 0 or the span of one nonzero {0,1}-vector of
    GF(p)^M (all of them when (2^M)^n <= limit, otherwise `limit` random ones).  Returns (min, count, all?)."""
    pts = [np.array(v) for v in itertools.product((0, 1), repeat=M) if any(v)]
    P = len(pts)
    tab = np.zeros(1 << P, dtype=np.int64)
    for S_ in range(1, 1 << P):
        tab[S_] = stabcert.gfp_rank(np.array([pts[i] for i in range(P) if S_ >> i & 1]).T.tolist(), p)
    K = P + 1; total = K ** n
    if total <= limit:
        idx = np.arange(total, dtype=np.int64)
        digits = np.stack([(idx // K ** i) % K for i in range(n)], axis=1); full = True
    else:
        digits = (rng or np.random.default_rng(1)).integers(0, K, size=(limit, n)); full = False
    bits = np.where(digits == 0, 0, np.left_shift(np.int64(1), np.maximum(digits - 1, 0)))
    vals = np.zeros(len(digits), dtype=np.int64)
    for X in range(1, 1 << n):
        if f[X]:
            u = np.zeros(len(digits), dtype=np.int64)
            for i in range(n):
                if X >> i & 1:
                    u |= bits[:, i]
            vals += f[X] * tab[u]
    return int(vals.min()), len(digits), full


def random_min(f, n, p, trials, rng, dmax=5):
    mn = None
    for k in range(trials):
        d = int(rng.integers(1, dmax))
        vecs = [rng.integers(0, p, size=(d, 1 if k % 3 == 0 else int(rng.integers(0, 3)))) for _ in range(n)]
        h = np.array([stabcert.gfp_rank(np.hstack([vecs[i] for i in range(n) if X >> i & 1]).tolist(), p)
                      if X and any(vecs[i].shape[1] for i in range(n) if X >> i & 1) else 0 for X in range(1 << n)])
        v = int(f @ h)
        mn = v if mn is None else min(mn, v)
    return mn


def exhaustive2_min(f, n, q, chunk=400000):
    """minimum of f over all arrangements of subspaces of GF(q)^2 (0, the q+1 lines, the plane)."""
    L = q + 1
    Sb = np.array([0] + [1 << i for i in range(L)] + [(1 << L) - 1], dtype=np.int64)
    K = len(Sb); total = K ** n
    masks = [X for X in range(1, 1 << n) if f[X]]
    mn = None
    for start in range(0, total, chunk):
        idx = np.arange(start, min(total, start + chunk), dtype=np.int64)
        bits = Sb[np.stack([(idx // K ** i) % K for i in range(n)], axis=1)]
        vals = np.zeros(len(idx), dtype=np.int64)
        for X in masks:
            u = np.zeros(len(idx), dtype=np.int64)
            for i in range(n):
                if X >> i & 1:
                    u |= bits[:, i]
            nz = u != 0
            vals += f[X] * np.where(~nz, 0, np.where(nz & ((u & (u - 1)) == 0), 1, 2))
        m = int(vals.min()); mn = m if mn is None else min(mn, m)
    return mn, total


def dfz1401_thm31():
    """Theorem 3.1 of DFZ arXiv:1401.2507v1 exactly as printed (p. 11), variables A, B, C, D, W, X, Y, Z = bits 0..7;
    vector f with f . H >= 0 claimed for characteristic != 3.  It is false as printed (docs, section 4.5)."""
    names = 'ABCDWXYZ'
    m = lambda st: sum(1 << names.index(ch) for ch in st)
    f = np.zeros(256, dtype=np.int64)
    for st, k in (('Z', 8), ('Y', 29), ('X', 3), ('W', 8), ('D', -6), ('C', -17), ('B', -8), ('A', -17)):
        f[m(st)] += k
    for x, y, k in (('Z', 'ABC', 55), ('Y', 'WXZ', 35), ('X', 'ACD', 50), ('W', 'BCD', 49), ('A', 'BDY', 18),
                    ('B', 'DXZ', 7), ('B', 'AWX', 1), ('C', 'DYZ', 7), ('C', 'BXY', 7), ('C', 'AWY', 3), ('D', 'AWZ', 6)):
        f[m(x) | m(y)] += k; f[m(y)] -= k
    for st in 'ABCD':
        f[m(st)] += 49
    f[m('ABCD')] -= 49
    f[m('A')] -= 1                         # left side H(A)
    return f


def dfz1401_report():
    f = dfz1401_thm31()
    # C = X = one line, everything else zero: H(T) = 1 iff T meets {C, X}
    h = np.array([1 if (T & 0b00100100) else 0 for T in range(256)])
    print(f'DFZ arXiv:1401.2507v1 Thm 3.1 as printed, on C = X = a line (all else 0): {int(f @ h)} (a valid inequality gives >= 0)')
    for q in (2, 3):
        L = q + 1
        Sb = np.array([0] + [1 << i for i in range(L)] + [(1 << L) - 1], dtype=np.int64)
        K = len(Sb); idx = np.arange(K ** 8, dtype=np.int64)
        bits = Sb[np.stack([(idx // K ** i) % K for i in range(8)], axis=1)]
        vals = np.zeros(len(idx), dtype=np.int64)
        for X in range(1, 256):
            if f[X]:
                u = np.zeros(len(idx), dtype=np.int64)
                for i in range(8):
                    if X >> i & 1:
                        u |= bits[:, i]
                nz = u != 0
                vals += f[X] * np.where(~nz, 0, np.where(nz & ((u & (u - 1)) == 0), 1, 2))
        print(f'  all {len(idx)} subspace arrangements of GF({q})^2: {int((vals < 0).sum())} violate it (min {int(vals.min())})')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--t', type=int, default=3)
    ap.add_argument('--fields', type=int, nargs='+', default=[3, 5, 7])
    ap.add_argument('--random', type=int, default=0)
    ap.add_argument('--exhaustive2', action='store_true')
    ap.add_argument('--exhaustive01', action='store_true')
    ap.add_argument('--dfz1401', action='store_true', help='show that DFZ arXiv:1401.2507 Thm 3.1 is false as printed')
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    if a.dfz1401:
        dfz1401_report(); return
    t0 = time.time(); t = a.t
    R = analyse(t, a.fields)
    n = R['n_visible']
    print(f't = {t}: {2 * t + 4} qudits, {n} visible parties + purifier; Lemma 4: LA [char | t], LB [char not | t]; '
          f'Pena-Sarria Example 6: (a) [char | t], (b) [char not | t]')
    for p, d in R['fields'].items():
        print(f"  GF({p}): weights {d['weights']}, identically self-dual {d['self_dual']}, h_S == 2r {d['h2r']}, "
              f"contraction == configuration {d['contraction']}; on the configuration LA {d['cfg_LA']}, LB {d['cfg_LB']} "
              f"[(a) {d['cfg_a']}, (b) {d['cfg_b']}]; pulled back on the state LA {d['FA']}, LB {d['FB']} [(a) {d['Fa']}, (b) {d['Fb']}]")
    ps = sorted(R['fields'])
    for i, p in enumerate(ps):
        for q in ps[i + 1:]:
            diff = np.flatnonzero(R['fields'][p]['S'] != R['fields'][q]['S'])
            print(f'  states over GF({p}) and GF({q}) differ in {diff.size} of {2 ** n - 1} coordinates')
    if a.random or a.exhaustive2 or a.exhaustive01:
        fa, fb, n_ = ps_forms(t, t + 1)
        LA, LB, _ = lemma_forms(t)
        rng = np.random.default_rng(28)
        for p in dict.fromkeys(a.fields + [2]):          # the given fields, then GF(2), each once
            div = t % p == 0
            pairs = [('LA', LA), ('(a)', fa)] if div else [('LB', LB), ('(b)', fb)]
            wrong = [('LB', LB), ('(b)', fb)] if div else [('LA', LA), ('(a)', fa)]
            for nm, f in pairs:
                if a.random:
                    print(f'  valid {nm} on {a.random} random GF({p}) arrangements: min {random_min(f, n_, p, a.random, rng)}', flush=True)
                if a.exhaustive2 and (p + 3) ** n_ <= 2.2e7:
                    mn, tot = exhaustive2_min(f, n_, p)
                    print(f'  valid {nm} on all {tot} subspace arrangements of GF({p})^2: min {mn}', flush=True)
                if a.exhaustive01:
                    mn, tot, full = exhaustive01_min(f, n_, t + 1, p, rng=rng)
                    print(f"  valid {nm} on {'all' if full else 'random'} {tot} {{0,1}}-point configurations of GF({p})^{t + 1}: min {mn}", flush=True)
            if a.exhaustive01:
                for nm, f in wrong:
                    mn, tot, full = exhaustive01_min(f, n_, t + 1, p, rng=rng)
                    print(f"  (sensitivity) wrong-characteristic {nm} on the same kind of configurations: min {mn}", flush=True)
    if a.out:
        np.savez_compressed(a.out, t=t, fields=np.array(ps), FA=R['FA'].astype(np.int32), FB=R['FB'].astype(np.int32),
                            Fa=R['Fa'].astype(np.int32), Fb=R['Fb'].astype(np.int32),
                            S=np.stack([R['fields'][p]['S'] for p in ps]).astype(np.int8), generator=generator(t).astype(np.int8))
        print('wrote', a.out)
    print(f'({time.time() - t0:.0f}s)')


if __name__ == '__main__':
    main()
