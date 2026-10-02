#!/usr/bin/env python3
"""s27: stabilizer entropy cones depend on the local dimension from seven parties on (R32).

Witnesses ("cube states").  Put one qudit on each vertex x of the cube {0,1}^3 and take the code of affine
functions f(x) = a0 + a.x evaluated at the eight vertices (generator rows 1, x1, x2, x3).  Over GF(2) this is
the first-order Reed-Muller code RM(1,3) = the extended Hamming [8,4,4] code (the CSS state is local-Clifford
equivalent to the 3-cube graph state); over GF(p), p odd, it is the same construction restricted to the cube
(not RM(1,3), which would use all p^3 points).  The state is the CSS state sum_{c in C}|c>; the seven nonzero
vertices are the visible parties, in DFZ's order A, B, C, W, X, Y, Z = 100, 010, 001, 110, 101, 011, 111, and the
origin 000 is the purifier.

Facts checked here (exact integer arithmetic unless stated):
  * both matroids are identically self-dual (r* = r), so S(X) = 2 r(X) - |X| and h_S = 2r;
  * contracting the origin gives exactly the seven vectors of DFZ Thms 8.7 / 8.9 -- the Fano configuration over
    GF(2), the non-Fano configuration over GF(p), p odd;
  * all 4 x 4 minors of the generator are in {0, +-1, +-2}, so for every odd p the matroid (and the entropy
    vector) is the same as over GF(3);
  * pulling DFZ (65) (valid in odd characteristic) back through "contract the purifier" and h_S gives a
    seven-party inequality on S valid for every stabilizer state of odd prime dimension; the qubit cube state
    gives -2.  Pulling back DFZ (91) (valid in characteristic 2) gives an inequality valid for every qubit
    stabilizer state; the qudit cube state gives -2.  (Validity rests on Lemma 2 of docs/A2_stabilizer_vs_QLR.md
    and on DFZ arXiv:1311.4601 Thms 8.6 and 8.8, data/dfz_char_dependent.csv.)
  * the two entropy vectors differ in exactly one coordinate: the cut {A,B,C,Z} | {W,X,Y,purifier} (odd- versus
    even-weight vertices), 2 for qubits and 4 for odd p;
  * Proposition 3 (doubled codes): for any configuration G over GF(p), the CSS state of [G | G] (party i holds
    column i, the copies form the purifier) has h_S = r_G + |X| on the visible parties; checked for the Fano,
    non-Fano and T8 configurations, where the balanced (65) / (91) take the value -1.  This gives a second,
    general route from characteristic-dependent n-variable inequalities to n-party separations;
  * optional: entropies against explicit state vectors (2^8 and 3^8 amplitudes); random-state sanity tests
    (not proofs); a six-party annealing probe (exploratory, negative so far).

  python3 scripts/s27_cube_states.py [--state-vector] [--random 1000] [--out data/s27_cube_witnesses.npz]
  python3 scripts/s27_cube_states.py --probe6 300        # six-party probe (exploratory)
"""
import argparse, csv, itertools, os, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..'))
from epr1kit import stabcert

DATA = os.path.join(HERE, '..', 'data')
DFZ_CHAR = os.path.join(DATA, 'dfz_char_dependent.csv')
NAMES = 'ABCWXYZ'                                   # visible parties, bits 0..6
VERTS = [(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0), (1, 0, 1), (0, 1, 1), (1, 1, 1), (0, 0, 0)]   # bit 7 = purifier
N8, FULL8, PUR = 8, 255, 7
ODD_CUT = 0b1000111                                 # {A, B, C, Z}: odd-weight vertices


def cube_generator():
    return np.array([[1] + list(v) for v in VERTS], dtype=np.int64).T          # 4 x 8, rows 1, x1, x2, x3


def rank_function(G, p, n=None):
    n = G.shape[1] if n is None else n
    return np.array([stabcert.gfp_rank(G[:, [j for j in range(n) if X >> j & 1]].tolist(), p) if X else 0
                     for X in range(1 << n)], dtype=np.int64)


def dual(r, n):
    full = (1 << n) - 1
    return np.array([bin(X).count('1') + r[full ^ X] - r[full] for X in range(1 << n)], dtype=np.int64)


def contract(r, n, e):
    rest = [i for i in range(n) if i != e]
    return np.array([r[sum(1 << rest[i] for i in range(n - 1) if X >> i & 1) | 1 << e] - r[1 << e]
                     for X in range(1 << (n - 1))], dtype=np.int64)


def css_S(r, n=N8):
    """entropy (units of log p) of the CSS state of a code with column rank function r, one qudit per party:
    S(X) = r(X) + r(E \\ X) - r(E) on all 2^n masks (purity built in)."""
    full = (1 << n) - 1
    return np.array([r[X] + r[full ^ X] - r[full] if 0 < X < full else 0 for X in range(1 << n)], dtype=np.int64)


def h_of_S(S, n=N8):
    return np.array([S[X] + sum(S[1 << i] for i in range(n) if X >> i & 1) for X in range(1 << n)], dtype=np.int64)


def load_dfz_char(path=DFZ_CHAR):
    """{ineq_id: coefficient vector over the 128 masks of A,B,C,W,X,Y,Z}; vector . H >= 0."""
    out = {}
    lines = (ln for ln in open(path) if not ln.lstrip().startswith('#'))
    for row in csv.DictReader(lines):
        m = sum(1 << NAMES.index(ch) for ch in row['subset'].strip())
        v = out.setdefault(row['ineq_id'], np.zeros(128, dtype=np.int64))
        v[m] += int(row['coeff'])
    return out


def dfz_configuration(p):
    """rank function of the seven vectors of DFZ Thms 8.7/8.9 over GF(p), in the order A, B, C, W, X, Y, Z."""
    V = np.array(VERTS[:7], dtype=np.int64).T                                  # 3 x 7: e1, e2, e3, e1+e2, ...
    return rank_function(V, p)


def pullback(c7):
    """seven-party inequality on S (127 coordinates, mask order over A..Z) equivalent to c7 applied to the
    contraction of h_S by the purifier:  sum_T c7[T] (h_S(T + purifier) - h_S(purifier)) >= 0."""
    coef = np.zeros(256, dtype=np.int64)

    def add_h(X, c):                       # c * h_S(X) = c * (S(X) + sum_{P in X} S(P))
        coef[X] += c
        for P in range(N8):
            if X >> P & 1:
                coef[1 << P] += c
    for T in range(1, 128):
        if c7[T]:
            add_h(T | 1 << PUR, int(c7[T])); add_h(1 << PUR, -int(c7[T]))
    fold = np.zeros(128, dtype=np.int64)
    for X in range(1, FULL8):
        Xr = X if not X >> PUR & 1 else FULL8 ^ X
        fold[Xr] += coef[X]
    return fold[1:]                                                            # coordinate = mask - 1


def minors_set(G):
    return sorted({int(round(np.linalg.det(G[:, list(c)]))) for c in itertools.combinations(range(G.shape[1]), 4)})


def witnesses():
    G = cube_generator()
    out = {}
    for p in (2, 3):
        r = rank_function(G, p)
        S = css_S(r)
        out[p] = dict(r=r, S=S, h=h_of_S(S))
    return G, out


def state_vector_check(G, p):
    """max |S_state-vector - S_formula| over all 254 proper cuts (8 qudits)."""
    import s26_reduction as s26
    psi = s26.css_state(G, p)
    S = css_S(rank_function(G, p))
    return max(abs(s26.sv_entropy(psi, p, N8, [j for j in range(N8) if X >> j & 1]) - S[X]) for X in range(1, FULL8))


def random_arrangement_min(c7, p, trials, rng):
    mn = None
    for t in range(trials):
        d = int(rng.integers(2, 6))
        blocks = [rng.integers(0, p, size=(d, 1 if t % 3 == 0 else int(rng.integers(0, 3)))) for _ in range(7)]
        h = np.array([stabcert.gfp_rank(np.hstack([blocks[i] for i in range(7) if X >> i & 1]).tolist(), p)
                      if X and any(blocks[i].shape[1] for i in range(7) if X >> i & 1) else 0 for X in range(128)])
        v = int(c7 @ h)
        mn = v if mn is None else min(mn, v)
    return mn


def random_state_min(F, p, trials, rng, n=7):
    """min of F . S over random weighted graph states over GF(p) with n visible parties + purifier, 1-3 qudits each."""
    import s26_reduction as s26
    mn = None
    for _ in range(trials):
        sizes = rng.integers(1, 4, size=n + 1)
        owner = [int(q) for q in np.repeat(np.arange(n + 1), sizes)]
        N = len(owner)
        W = (rng.random((N, N)) < rng.uniform(0.15, 0.85)) * rng.integers(1, p, size=(N, N))
        W = np.triu(W, 1); W = (W + W.T) % p
        v = int(F @ s26.graph_entropy(W, owner, n, p))
        mn = v if mn is None else min(mn, v)
    return mn


def balanced_value(f, h, m):
    """value of the balanced version f_bal(h) = f(h) - sum_V f(e_V) H(V | rest) of a form f on m elements."""
    full = (1 << m) - 1
    fe = [int(sum(f[X] for X in range(1 << m) if X >> V & 1)) for V in range(m)]
    return int(f @ h) - sum(fe[V] * int(h[full] - h[full ^ (1 << V)]) for V in range(m))


def doubled_code_h(G, p):
    """Proposition 3 (R32): CSS state of the code generated by [G | G]; party i holds column i, the m copies form
    the purifier.  Returns h_S restricted to the m visible parties (2^m masks), computed from the code."""
    m = G.shape[1]
    r = rank_function(np.hstack([G, G]) % p, p)
    full = (1 << 2 * m) - 1
    S = lambda Y: int(r[Y] + r[full ^ Y] - r[full]) if 0 < Y < full else 0
    S1 = [S(1 << i) for i in range(m)]
    return np.array([S(X) + sum(S1[i] for i in range(m) if X >> i & 1) if X else 0 for X in range(1 << m)], dtype=np.int64)


def doubled_code_checks():
    """h_S of the doubled Fano (GF(2)), non-Fano (GF(3)) and T8 (GF(3)) codes restricted to the visible parties is
    M + |X| (M = rank function of the configuration); balanced (65) / (91) values on it."""
    D = load_dfz_char()
    V = np.array(VERTS[:7], dtype=np.int64).T
    J = np.ones((4, 4), dtype=np.int64); I4 = np.eye(4, dtype=np.int64)
    out = {}
    for name, G, p, key in (('Fano', V, 2, 'dfz65_odd'), ('non-Fano', V, 3, 'dfz91_even'),
                            ('T8', np.hstack([I4, J - I4]), 3, None)):
        m = G.shape[1]
        h = doubled_code_h(G, p)
        M = rank_function(G % p, p)
        ok = bool(np.array_equal(h, M + np.array([bin(X).count('1') for X in range(1 << m)])))
        out[name] = (ok, balanced_value(D[key], h, m) if key else None)
    return out


def symmetrised_minima():
    """six-party check: h_S of the seven-qudit CSS states of the Fano (GF(2)) and non-Fano (GF(3)) codes is r + r*;
    minimum over all 5040 relabellings of (65) resp. (91) on it (negative would separate at n = 6)."""
    D = load_dfz_char()
    out = {}
    for name, p, key in (('Fano', 2, 'dfz65_odd'), ('non-Fano', 3, 'dfz91_even')):
        r = dfz_configuration(p)
        h = r + dual(r, 7)
        c = D[key]
        best = None
        for perm in itertools.permutations(range(7)):
            hp = np.zeros(128, dtype=np.int64)
            for X in range(128):
                hp[sum(1 << perm[i] for i in range(7) if X >> i & 1)] = h[X]
            v = int(c @ hp)
            best = v if best is None else min(best, v)
        out[name] = best
    return out


def probe6(seconds, seed=6):
    """exploratory: anneal 7-party (6 + purifier) GF(3) graph states, 1-2 qutrits per party, minimising
    min over all 5040 relabellings of DFZ (91) evaluated on h_S (seven elements, no contraction)."""
    import s26_reduction as s26
    c91 = load_dfz_char()['dfz91_even']
    idx = np.flatnonzero(c91); coefs = c91[idx]
    P = np.array([[sum(1 << perm[i] for i in range(7) if T >> i & 1) for T in idx]
                  for perm in itertools.permutations(range(7))])
    rng = np.random.default_rng(seed); t0 = time.time(); best = None
    while time.time() - t0 < seconds:
        sizes = rng.integers(1, 3, size=7); owner = [int(q) for q in np.repeat(np.arange(7), sizes)]; N = len(owner)
        W = np.triu(rng.integers(0, 3, size=(N, N)), 1); W = (W + W.T) % 3

        def score(W):
            S6 = s26.graph_entropy(W, owner, 6, 3)                # 63 coordinates, purifier = party 6
            S = np.zeros(128, dtype=np.int64)
            for X in range(1, 127):
                S[X] = S6[(X if not X >> 6 & 1 else 127 ^ X) - 1]
            h = h_of_S(S, 7)
            return int((h[P] @ coefs).min())
        cur = score(W); T = 2.0
        for _ in range(300):
            i, j = rng.choice(N, 2, replace=False)
            if owner[i] == owner[j]:
                continue
            old = W[i, j]; W[i, j] = W[j, i] = (old + rng.integers(1, 3)) % 3
            s = score(W)
            if s <= cur or rng.random() < np.exp((cur - s) / T):
                cur = s
            else:
                W[i, j] = W[j, i] = old
            T *= 0.99
            if time.time() - t0 > seconds:
                break
        best = cur if best is None else min(best, cur)
        if cur < 0:
            break
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--state-vector', action='store_true')
    ap.add_argument('--random', type=int, default=0, help='random-state sanity tests per field')
    ap.add_argument('--probe6', type=float, default=0, help='seconds for the exploratory six-party probe')
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    t0 = time.time()
    G, W = witnesses()
    D = load_dfz_char()
    F65, F91 = pullback(D['dfz65_odd']), pullback(D['dfz91_even'])
    for p, name, c7 in ((2, 'Fano', D['dfz65_odd']), (3, 'non-Fano', D['dfz91_even'])):
        r = W[p]['r']
        print(f"GF({p}) cube code: rank {r[-1]}, identically self-dual {bool(np.array_equal(dual(r, 8), r))}, "
              f"S == 2r - |X|: {all(W[p]['S'][X] == 2 * r[X] - bin(X).count('1') for X in range(1, FULL8))}, "
              f"h_S == 2r: {bool(np.array_equal(W[p]['h'][1:FULL8], 2 * r[1:FULL8]))}, "
              f"contraction at 000 == DFZ {name} configuration: {bool(np.array_equal(contract(r, 8, PUR), dfz_configuration(p)))}, "
              f"#coplanar 4-sets {sum(1 for X in range(256) if bin(X).count('1') == 4 and r[X] == 3)}")
    print('4x4 minors of the generator:', minors_set(G))
    print('DFZ forms on their own configurations: (65) on Fano/GF(2) =', int(D['dfz65_odd'] @ dfz_configuration(2)),
          '; (91) on non-Fano/GF(3) =', int(D['dfz91_even'] @ dfz_configuration(3)))
    S2, S3 = W[2]['S'][1:128], W[3]['S'][1:128]
    print(f'pulled-back (65) [valid for odd p]: qubit cube state {int(F65 @ S2)}, qudit cube state {int(F65 @ S3)}')
    print(f'pulled-back (91) [valid for p = 2]: qubit cube state {int(F91 @ S2)}, qudit cube state {int(F91 @ S3)}')
    diff = np.flatnonzero(S2 != S3)
    print('coordinates where the two cube states differ:',
          [(''.join(NAMES[i] for i in range(7) if (k + 1) >> i & 1), int(S2[k]), int(S3[k])) for k in diff])
    dc = doubled_code_checks()
    print('doubled codes (Proposition 3): h_S on the visible parties == M + |X| and balanced DFZ value: '
          + '; '.join(f'{k}: {v[0]}' + ('' if v[1] is None else f', {v[1]}') for k, v in dc.items()))
    if a.state_vector:
        print(f'state-vector check: GF(2) max dev {state_vector_check(G, 2):.1e}; GF(3) max dev {state_vector_check(G, 3):.1e}')
    if a.random:
        rng = np.random.default_rng(27)
        print(f"random arrangements: (65) min over {a.random} GF(3) = {random_arrangement_min(D['dfz65_odd'], 3, a.random, rng)}, "
              f"GF(5) = {random_arrangement_min(D['dfz65_odd'], 5, a.random, rng)}; "
              f"(91) min over {a.random} GF(2) = {random_arrangement_min(D['dfz91_even'], 2, a.random, rng)}")
        print(f'random graph states (7 parties + purifier, 1-3 qudits each): pulled-back (65) min over {a.random} qutrit '
              f'states = {random_state_min(F65, 3, a.random, rng)}, over {a.random // 2} p=5 states = '
              f'{random_state_min(F65, 5, a.random // 2, rng)}; pulled-back (91) min over {a.random} qubit states = '
              f'{random_state_min(F91, 2, a.random, rng)}')
    if a.probe6:
        sm = symmetrised_minima()
        print(f"six-party check: min over relabellings of (65) on Fano r + r* = {sm['Fano']}, "
              f"of (91) on non-Fano r + r* = {sm['non-Fano']} (negative would separate)")
        print(f'six-party probe ({a.probe6:.0f} s): best min over relabellings of (91) on h_S = {probe6(a.probe6)}')
    if a.out:
        np.savez_compressed(a.out, generator=G.astype(np.int8), S_qubit=S2.astype(np.int8), S_qudit=S3.astype(np.int8),
                            F65_odd=F65.astype(np.int16), F91_even=F91.astype(np.int16),
                            parties=np.array(list(NAMES) + ['000']), vertices=np.array(VERTS, dtype=np.int8))
        print('wrote', a.out)
    print(f'({time.time() - t0:.0f}s)')


if __name__ == '__main__':
    main()
