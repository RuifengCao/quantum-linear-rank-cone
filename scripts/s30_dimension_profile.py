#!/usr/bin/env python3
"""s30: how far down does the dependence on the local dimension go? (R35)

docs/R32_reduction_and_local_dimension.md, section 4.7.  Three exact checks:

  * Theorem 7 (four parties, every d): each of the 46 extreme rays of Stab_4 = h^-1(Shannon + Ingleton) has an
    integer graph state (signs on the edges of the qubit certificate of s26) whose cut blocks have Smith invariants
    0 and 1 only and rational rank lam * r.  Reduced mod p it realises lam * r over GF(p) for EVERY prime p, so
    Stab^(Z_d)_4 = h^-1(Shannon + Ingleton) for every d >= 2 (data/s30_stab4_unimodular.npz, `--n4`).
  * Proposition 6 (one cut): over a field of characteristic l not dividing t - 1, the rank function of the generalised
    cube points (vectors (1, x), x in {e_j, c - e_j, c, 0}) differs from the rational one exactly on the "balanced"
    sets X (one of e_j, c - e_j for every j; a = #{j : c - e_j in X}) with c in X, 0 not in X, l | t - a, a != t, or
    0 in X, c not in X, l | a - 1, a != 1, by exactly 1 (`--one-cut`).
  * Proposition 8 (six parties are blind): for the six known 7-variable characteristic-dependent inequalities f
    (DFZ (65), (91); Lemma 4 (LA)_2, (LB)_2; Pena-Sarria (a), (b) with t = 2), f + f o D (D: h -> h*) is a nonnegative
    rational combination of Shannon elemental inequalities.  The certificates are stored exactly
    (data/s30_selfdual_shannon_certs.npz, `--six`; `--solve` recomputes them with an LP and exact rational solving).

  python3 scripts/s30_dimension_profile.py --n4 --one-cut 2 3 4 5 --six
  python3 scripts/s30_dimension_profile.py --n4-search --out data/s30_stab4_unimodular.npz     # regenerate the lifts
  python3 scripts/s30_dimension_profile.py --solve --out-certs data/s30_selfdual_shannon_certs.npz
  python3 scripts/s30_dimension_profile.py --probe6 --n5-pilot 60                          # sanity check, pilot
"""
import argparse, itertools, os, sys, time
from fractions import Fraction
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..'))
from epr1kit import core, stabcert
import s27_cube_states as s27
import s28_generalized_cubes as g28
import s29_zd_stabilizer as z29

N4_CERTS = os.path.join(core.DATA, 's26_stab4_certs.npz')
N4_LIFTS = os.path.join(core.DATA, 's30_stab4_unimodular.npz')
SIX_CERTS = os.path.join(core.DATA, 's30_selfdual_shannon_certs.npz')
PRIMES = (2, 3, 5, 7, 11, 13)


# ------------------------------------------------------------------ Theorem 7: four parties over every field
def cut_blocks(sizes, n=4):
    party = [q for q in range(n + 1) for _ in range(int(sizes[q]))]
    N = len(party)
    out = []
    for m in range(1, 1 << n):
        X = [i for i in range(N) if party[i] < n and m >> party[i] & 1]
        out.append((m, X, [i for i in range(N) if i not in X]))
    return out


def unimodular_ok(Wt, sizes, r, lam, n=4):
    """every visible cut block of the integer matrix Wt has Smith invariants in {0, 1} and rational rank lam * r."""
    for m, X, Y in cut_blocks(sizes, n):
        want = lam * int(r[m - 1])
        if not X or not Y:
            if want:
                return False
            continue
        d = z29.diag_form([[int(Wt[i][j]) for j in Y] for i in X])
        if len(d) != want or any(abs(x) != 1 for x in d):
            return False
    return True


def lift_search(r, lam, sizes, W, rng, tries=20000):
    """signs on the edges of the 0/1 adjacency W making every cut block unimodular (all +1 first, then exhaustive
    for <= 2^14 sign patterns, else random)."""
    N = W.shape[0]
    E = [(i, j) for i in range(N) for j in range(i + 1, N) if W[i, j]]

    def build(signs):
        Wt = np.zeros((N, N), dtype=np.int64)
        for (i, j), s in zip(E, signs):
            Wt[i, j] = Wt[j, i] = s
        return Wt
    first = build([1] * len(E))
    if unimodular_ok(first, sizes, r, lam):
        return first
    pats = itertools.product((1, -1), repeat=len(E)) if len(E) <= 14 else (rng.choice((1, -1), size=len(E)) for _ in range(tries))
    for signs in pats:
        Wt = build(signs)
        if unimodular_ok(Wt, sizes, r, lam):
            return Wt
    return None


def n4_search(out=None):
    Z = np.load(N4_CERTS)
    rng = np.random.default_rng(4)
    lifts = np.zeros((len(Z['rays']), 8, 8), dtype=np.int8)
    for k in range(len(Z['rays'])):
        N = int(Z['N'][k])
        Wt = lift_search(Z['rays'][k], int(Z['lam'][k]), Z['sizes'][k], np.array(Z['W'][k][:N, :N], dtype=np.int64), rng)
        if Wt is None:
            raise SystemExit(f'no unimodular lift found for ray {k}')
        lifts[k, :N, :N] = Wt
    if out:
        np.savez_compressed(out, rays=Z['rays'], lam=Z['lam'], sizes=Z['sizes'], N=Z['N'], Wt=lifts)
        print('wrote', out)
    return lifts


def n4_check(path=N4_LIFTS):
    """exact checks of the stored lifts: Wt = W (s26) mod 2, unimodular cut blocks, and the repository verifier over
    GF(p) for p = 2, 3, 5, 7, 11, 13.  Returns (number of rays, number passing)."""
    L = np.load(path); C = np.load(N4_CERTS)
    good = 0
    for k in range(len(L['rays'])):
        N = int(L['N'][k]); lam = int(L['lam'][k]); sizes = [int(s) for s in L['sizes'][k]]
        r = [int(x) for x in L['rays'][k]]
        Wt = np.array(L['Wt'][k][:N, :N], dtype=np.int64)
        ok = (np.array_equal(L['rays'][k], C['rays'][k]) and np.array_equal(Wt % 2, C['W'][k][:N, :N])
              and np.array_equal(Wt, Wt.T) and not np.diag(Wt).any() and unimodular_ok(Wt, sizes, r, lam)
              and all(stabcert.verify_realisation_n(r, 4, lam, sizes, (Wt % p).tolist(), p) for p in PRIMES))
        good += ok
    return len(L['rays']), good


# ------------------------------------------------------------------ Proposition 6: one cut
def one_cut_predicted(t, l):
    M = t + 1
    out = set()
    for choice in range(1 << M):          # bit j: c - e_j in X (else e_j)
        base = sum((1 << (M + j)) if choice >> j & 1 else (1 << j) for j in range(M))
        a = bin(choice).count('1')
        if (t - a) % l == 0 and a != t:
            out.add(base | 1 << (2 * M))          # with c, without 0
        if (a - 1) % l == 0 and a != 1:
            out.add(base | 1 << (2 * M + 1))      # with 0, without c
    return out


def one_cut_check(t, primes=PRIMES):
    """for each prime l not dividing t - 1: (sets where rank over GF(l) != rank over Q) == predicted, drops all 1;
    returns {l: (changed sets, differing entropy coordinates, ok)}."""
    G = g28.generator(t); m = G.shape[1]
    ls = [l for l in primes if (t - 1) % l]
    rQ = np.zeros(1 << m, dtype=np.int64); rl = {l: np.zeros(1 << m, dtype=np.int64) for l in ls}
    for X in range(1, 1 << m):
        d = z29.diag_form(G[:, [i for i in range(m) if X >> i & 1]].tolist())
        rQ[X] = sum(1 for x in d if x)
        for l in ls:
            rl[l][X] = sum(1 for x in d if x % l)
    out = {}
    for l in ls:
        diff = np.flatnonzero(rl[l] != rQ)
        pred = one_cut_predicted(t, l)
        ok = set(diff.tolist()) == pred and set((rQ - rl[l])[diff].tolist()) <= {1}
        out[l] = (len(diff), sum(1 for X in pred if not X >> (m - 1) & 1), ok)
    return out


# ------------------------------------------------------------------ Proposition 8: six parties are blind
def elemental7(n=7):
    """the Shannon elemental inequalities on n variables, as rows over the masks 1..2^n - 1 (fixed order)."""
    FULL = (1 << n) - 1
    rows = []
    for i in range(n):
        v = np.zeros(1 << n, dtype=np.int64); v[FULL] += 1; v[FULL ^ (1 << i)] -= 1; rows.append(v)
    for i, j in itertools.combinations(range(n), 2):
        rest = [k for k in range(n) if k not in (i, j)]
        for s in range(1 << len(rest)):
            K = sum(1 << rest[b] for b in range(len(rest)) if s >> b & 1)
            v = np.zeros(1 << n, dtype=np.int64)
            v[K | 1 << i] += 1; v[K | 1 << j] += 1; v[K | 1 << i | 1 << j] -= 1
            if K:
                v[K] -= 1
            rows.append(v)
    return np.array(rows)[:, 1:]


def selfdual_symmetrise(f, n=7):
    """g = f + f o D, where D h = h*, h*(X) = sum_{P in X} h(P) + h(E \\ X) - h(E); g . h = f(h) + f(h*)."""
    FULL = (1 << n) - 1
    g = np.asarray(f, dtype=np.int64).copy()
    for X in range(1, 1 << n):
        c = int(f[X])
        if c:
            for P in range(n):
                if X >> P & 1:
                    g[1 << P] += c
            if FULL ^ X:
                g[FULL ^ X] += c
            g[FULL] -= c
    return g


def six_forms():
    LA, LB, _ = g28.lemma_forms(2)
    a, b, _n = g28.ps_forms(2, 3)
    F = s27.load_dfz_char()
    return [('DFZ (65)', F['dfz65_odd']), ('DFZ (91)', F['dfz91_even']), ('(LA)_2', LA), ('(LB)_2', LB),
            ('Pena-Sarria (a), t=2', a), ('Pena-Sarria (b), t=2', b)]


def six_solve(out=None):
    """LP (HiGHS) for a sparse nonnegative combination, then exact rational solution on its support."""
    from scipy.optimize import linprog
    import sympy as sp
    E = elemental7()
    store = {}
    for k, (name, f) in enumerate(six_forms()):
        g = selfdual_symmetrise(np.asarray(f, dtype=np.int64))
        res = linprog(np.ones(E.shape[0]), A_eq=E.T, b_eq=g[1:].astype(float), bounds=(0, None), method='highs')
        if res.status != 0:
            raise SystemExit(f'{name}: no Shannon certificate')
        supp = np.flatnonzero(res.x > 1e-9)
        sol, params = sp.Matrix(E[supp].T.tolist()).gauss_jordan_solve(sp.Matrix([int(x) for x in g[1:]]))
        sol = sol.subs({s_: 0 for s_ in params})
        fr = [Fraction(int(sp.fraction(x)[0]), int(sp.fraction(x)[1])) for x in sol]
        den = 1
        for x in fr:
            den = den * x.denominator // np.gcd(den, x.denominator)
        num = [int(x * den) for x in fr]
        store[f'idx{k}'] = supp.astype(np.int32); store[f'num{k}'] = np.array(num, dtype=np.int64)
        store[f'den{k}'] = np.array(den, dtype=np.int64)
        print(f'  {name}: support {len(supp)}, denominator {den}, min numerator {min(num)}')
    if out:
        np.savez_compressed(out, **store)
        print('wrote', out)


def six_check(path=SIX_CERTS):
    """exact: den * g == sum num_i * e_i with all num_i >= 0, for each of the six forms.  Returns list of (name, ok)."""
    C = np.load(path); E = elemental7()
    out = []
    for k, (name, f) in enumerate(six_forms()):
        g = selfdual_symmetrise(np.asarray(f, dtype=np.int64))
        idx, num, den = C[f'idx{k}'], C[f'num{k}'], int(C[f'den{k}'])
        lhs = (num[:, None] * E[idx]).sum(axis=0)
        out.append((name, bool((num >= 0).all() and den > 0 and np.array_equal(lhs, den * g[1:]))))
    return out


def probe6(rng=None):
    """sanity check of Proposition 8: f(r + r*) = (f + f o D) . r on all {0,1}-point configurations of GF(p)^3 and on
    2.2 M random ones of GF(p)^4, for the forms over the characteristic where f itself is NOT valid."""
    rng = rng or np.random.default_rng(6)
    forms = dict(six_forms())
    out = []
    for name, p in (('(LB)_2', 2), ('DFZ (65)', 2), ('(LA)_2', 3), ('DFZ (91)', 3), ('(LA)_2', 5), ('DFZ (91)', 5)):
        f = np.asarray(forms[name], dtype=np.int64); g = selfdual_symmetrise(f)
        m0 = g28.exhaustive01_min(f, 7, 3, p)[0]
        m1 = g28.exhaustive01_min(g, 7, 3, p)[0]
        m2 = g28.exhaustive01_min(g, 7, 4, p, rng=rng)[0]
        out.append((name, p, m0, m1, m2))
    return out


def n5_pilot(count, tries=3000, seed=5):
    """pilot (not a result): signs on the qubit certificates (data/s16_realisations.npz) of the first `count` known
    extreme rays of Stab_5 (data/stab5_extreme_reps.npy, catalogue order) making every visible cut block unimodular.
    All +1 first, then a random local search flipping edges inside failing blocks (GF(3), GF(5), GF(7) rank filters),
    then an exact Smith-form check.  Returns (all +1, after search, failures)."""
    Z = np.load(os.path.join(core.DATA, 's16_realisations.npz'))
    want = {tuple(r) for r in np.load(os.path.join(core.DATA, 'stab5_extreme_reps.npy')).tolist()}
    todo, seen = [], set()
    for k in range(len(Z['rep'])):
        key = tuple(Z['rep'][k].tolist())
        if key in want and key not in seen:
            seen.add(key); todo.append(k)
    rng = np.random.default_rng(seed)

    def rank_p(B, p):
        A = np.array(B, dtype=np.int64) % p; r = 0
        for c in range(A.shape[1]):
            piv = next((i for i in range(r, A.shape[0]) if A[i, c]), None)
            if piv is None:
                continue
            A[[r, piv]] = A[[piv, r]]; A[r] = (A[r] * pow(int(A[r, c]), p - 2, p)) % p
            for i in np.nonzero(A[:, c])[0]:
                if i != r:
                    A[i] = (A[i] - A[i, c] * A[r]) % p
            r += 1
            if r == A.shape[0]:
                break
        return r
    plus = searched = 0; fails = []
    for k in todo[:count]:
        sizes = [int(x) for x in Z['sizes'][k]]; N = sum(sizes); lam = int(Z['lam'][k]); r = Z['rep'][k]
        W = np.array([[(int(Z['rows'][k][i]) >> j) & 1 for j in range(N)] for i in range(N)], dtype=np.int64)
        blocks = cut_blocks(sizes, 5); target = [lam * int(r[m - 1]) for m in range(1, 32)]
        E = [(i, j) for i in range(N) for j in range(i + 1, N) if W[i, j]]
        bad_of = lambda Wt: [m for (m, X, Y), tg in zip(blocks, target) if X and Y and
                             any(rank_p(Wt[np.ix_(X, Y)], p) != tg for p in (3, 5, 7))]
        Wt = W.copy(); bad = bad_of(Wt); flips = 0
        while bad and flips < tries:
            flips += 1
            m = bad[rng.integers(len(bad))]; X = set(blocks[m - 1][1]); Y = set(blocks[m - 1][2])
            cand = [(i, j) for (i, j) in E if (i in X and j in Y) or (j in X and i in Y)]
            if not cand:
                break
            i, j = cand[rng.integers(len(cand))]
            Wt[i, j] = Wt[j, i] = -Wt[i, j]
            nb = bad_of(Wt)
            if len(nb) > len(bad) and rng.random() < 0.7:
                Wt[i, j] = Wt[j, i] = -Wt[i, j]
            else:
                bad = nb
        if not bad and unimodular_ok(Wt, sizes, r, lam, n=5):
            if flips:
                searched += 1
            else:
                plus += 1
        else:
            fails.append(int(k))
    return plus, searched, fails


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--n4', action='store_true', help='check the stored unimodular four-party lifts')
    ap.add_argument('--n4-search', action='store_true', help='recompute the lifts')
    ap.add_argument('--out', default=None)
    ap.add_argument('--one-cut', type=int, nargs='*', default=None, help='values of t')
    ap.add_argument('--six', action='store_true', help='check the stored six-party Shannon certificates')
    ap.add_argument('--solve', action='store_true', help='recompute the six-party certificates')
    ap.add_argument('--out-certs', default=None)
    ap.add_argument('--probe6', action='store_true', help='f(r + r*) on point configurations (sanity of Proposition 8)')
    ap.add_argument('--n5-pilot', type=int, default=0, help='five-party lifting pilot on the first K known extreme rays')
    a = ap.parse_args()
    t0 = time.time()
    if a.n4_search:
        n4_search(a.out)
    if a.n4:
        tot, good = n4_check(a.out or N4_LIFTS)
        print(f'four parties: {good}/{tot} extreme rays of h^-1(Shannon + Ingleton) have integer graph states with unimodular '
              f'cut blocks (re-verified over GF(p), p = {", ".join(map(str, PRIMES))})  ({time.time() - t0:.0f}s)', flush=True)
    if a.one_cut is not None:
        for t in a.one_cut:
            res = one_cut_check(t)
            print(f'one cut, t = {t}: ' + '; '.join(f'GF({l}) vs Q: {c} sets change rank, {k} entropy coordinates, as predicted {ok}'
                                                   for l, (c, k, ok) in res.items()) + f'  ({time.time() - t0:.0f}s)', flush=True)
    if a.solve:
        six_solve(a.out_certs)
    if a.six:
        for name, ok in six_check(a.out_certs or SIX_CERTS):
            print(f'six parties: {name}: f + f o D is an exact nonnegative combination of Shannon inequalities: {ok}', flush=True)
    if a.probe6:
        for name, p, m0, m1, m2 in probe6():
            print(f'probe6: {name} over GF({p}): f on all {{0,1}}-point configurations of GF({p})^3: min {m0}; f(r + r*) there: '
                  f'min {m1}; f(r + r*) on 2.2 M random ones of GF({p})^4: min {m2}  ({time.time() - t0:.0f}s)', flush=True)
    if a.n5_pilot:
        plus, searched, fails = n5_pilot(a.n5_pilot)
        print(f'five-party pilot (first {a.n5_pilot} known extreme rays of Stab_5): lifted with all +1: {plus}, after sign '
              f'search: {searched}, not lifted: {len(fails)} (certificate rows {fails})  ({time.time() - t0:.0f}s)', flush=True)


if __name__ == '__main__':
    main()
