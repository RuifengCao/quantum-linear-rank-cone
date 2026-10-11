"""Checks for Sec. III / Appendix A on random stabilizer states, independent of the repository code.

For random stabilizer states of Z_d qudits (d = 2, 3, 4, 6, 8, 9, 12; products of one-qudit stabilizer states,
including the non-free ones sum_{k in aZ_d} |k>, followed by random Clifford circuits), with the qudits distributed
among a few systems, the script computes the entropy vector S from the state vector (SVD) and checks
  - Lemma A.1, both forms of the entropy formula, with the Lagrangian subgroup L read off from the state;
  - Lemma A.4: h_S(X) = log|pr_X Lambda|, and h_S is the rank function of explicit character subgroups U_i;
  - that h_S is a polymatroid, Lemma A.2 (h_S* = h_S and lambda(h_S) = 2S), and Corollary 1 (the CSS state of L
    on 2m qudits has entropy vector 2S);
  - that non-free Lagrangian subgroups occur for d = 4, 8, 9, 12 and never for the squarefree d = 2, 3, 6.
Entropies are in bits; all comparisons are to within 1e-8.
"""
import itertools, math, random, sys
import numpy as np

# ----------------------------------------------------------------------------- state vectors
def entropy_bits(psi, d, m, Xq):
    """von Neumann entropy (bits) of the reduced state on the qudit list Xq of the pure state psi."""
    Xq = sorted(Xq)
    if len(Xq) == 0 or len(Xq) == m:
        return 0.0
    Yq = [j for j in range(m) if j not in Xq]
    T = psi.reshape([d] * m)
    M = np.transpose(T, Xq + Yq).reshape(d ** len(Xq), d ** len(Yq))
    s = np.linalg.svd(M, compute_uv=False) ** 2
    s = s[s > 1e-14]
    s = s / s.sum()
    return float(-(s * np.log2(s)).sum())


def entropy_vector(psi, d, m, owner, nsys):
    """dict X(frozenset of systems) -> S(X) in bits, for all subsets of range(nsys)."""
    out = {}
    for r in range(nsys + 1):
        for X in itertools.combinations(range(nsys), r):
            Xq = [j for j in range(m) if owner[j] in X]
            out[frozenset(X)] = entropy_bits(psi, d, m, Xq)
    return out


def apply_F(psi, d, m, a):
    T = psi.reshape([d] * m)
    T = np.fft.ifft(T, axis=a) * math.sqrt(d)  # (F psi)(k) = d^{-1/2} sum_j w^{jk} psi(j)
    return T.reshape(-1)


def apply_P(psi, d, m, a):
    c = d % 2
    ph = np.exp(1j * np.pi * np.array([j * (j + c) for j in range(d)]) / d)
    T = psi.reshape([d] * m)
    shape = [1] * m
    shape[a] = d
    T = T * ph.reshape(shape)
    return T.reshape(-1)


def apply_M(psi, d, m, a, u):
    """|j> -> |u j> on qudit a (u a unit mod d)."""
    T = psi.reshape([d] * m)
    uinv = pow(u, -1, d)
    idx = [(uinv * k) % d for k in range(d)]
    T = np.take(T, idx, axis=a)  # psi'(k) = psi(u^{-1} k)
    return T.reshape(-1)


def apply_SUM(psi, d, m, a, b):
    """|j_a, j_b> -> |j_a, j_b + j_a>."""
    T = psi.reshape([d] * m)
    out = np.empty_like(T)
    for ja in range(d):
        sl = [slice(None)] * m
        sl[a] = ja
        sub = T[tuple(sl)]               # axes: all but a; b shifts index
        bb = b if b < a else b - 1
        out[tuple(sl)] = np.roll(sub, ja, axis=bb)  # psi'(k_b) = psi(k_b - j_a)
    return out.reshape(-1)


def apply_X(psi, d, m, a, s=1):
    T = psi.reshape([d] * m)
    return np.roll(T, s, axis=a).reshape(-1)


def apply_Z(psi, d, m, a, s=1):
    w = np.exp(2j * np.pi / d)
    T = psi.reshape([d] * m)
    shape = [1] * m
    shape[a] = d
    return (T * (w ** (s * np.arange(d))).reshape(shape)).reshape(-1)


def random_clifford(psi, d, m, rng, depth):
    units = [u for u in range(1, d) if math.gcd(u, d) == 1]
    for _ in range(depth):
        g = rng.randrange(6)
        a = rng.randrange(m)
        if g == 0:
            psi = apply_F(psi, d, m, a)
        elif g == 1:
            psi = apply_P(psi, d, m, a)
        elif g == 2:
            psi = apply_M(psi, d, m, a, rng.choice(units))
        elif g == 3 and m > 1:
            b = rng.randrange(m - 1)
            b = b if b < a else b + 1
            psi = apply_SUM(psi, d, m, a, b)
        elif g == 4:
            psi = apply_X(psi, d, m, a, rng.randrange(d))
        else:
            psi = apply_Z(psi, d, m, a, rng.randrange(d))
    return psi


def lagrangian_from_state(psi, d, m, tol=1e-7):
    """L = {v=(x,z) in Z_d^{2m} : |<psi|Z^z X^x|psi>| = 1}, by brute force with FFTs.
    Returns a set of tuples (x_1..x_m, z_1..z_m)."""
    T = psi.reshape([d] * m)
    L = set()
    for x in itertools.product(range(d), repeat=m):
        shifted = np.roll(T, x, axis=tuple(range(m))) if m > 0 else T   # psi(j - x)
        phi = np.conj(T) * shifted
        E = np.fft.ifftn(phi) * (d ** m)  # E(z) = sum_j phi(j) w^{z.j}
        for z in zip(*np.nonzero(np.abs(np.abs(E) - 1) < tol)):
            L.add(tuple(x) + tuple(int(t) for t in z))
    return L


def symp(v, w, m, d):
    x, z = v[:m], v[m:]
    x2, z2 = w[:m], w[m:]
    return (sum(a * b for a, b in zip(x2, z)) - sum(a * b for a, b in zip(x, z2))) % d


# ----------------------------------------------------------------------------- groups
def span(gens, d, r):
    """subgroup of Z_d^r generated by gens (set of tuples), BFS closure."""
    G = {tuple([0] * r)}
    frontier = [tuple([0] * r)]
    gens = [tuple(int(t) % d for t in g) for g in gens]
    while frontier:
        new = []
        for v in frontier:
            for g in gens:
                w = tuple((a + b) % d for a, b in zip(v, g))
                if w not in G:
                    G.add(w)
                    new.append(w)
        frontier = new
    return G


def generators(S, d, r):
    """a small generating set of the subgroup S (a set of tuples) of Z_d^r."""
    gens = []
    cur = {tuple([0] * r)}
    for v in sorted(S):
        if v not in cur:
            gens.append(v)
            cur = span(gens, d, r)
            if len(cur) == len(S):
                break
    return gens


def annihilator_dot(S, d, r):
    """{a in Z_d^r : a.s = 0 for all s in S} (brute force over Z_d^r; r small)."""
    gS = generators(S, d, r) if S else []
    out = []
    for a in itertools.product(range(d), repeat=r):
        if all(sum(x * y for x, y in zip(a, s)) % d == 0 for s in gS):
            out.append(a)
    return out


def proj(S, cols):
    return {tuple(v[c] for c in cols) for v in S}


def isomorphism_type(S, d, r):
    """invariant factors of the finite abelian group S <= Z_d^r (sorted list of orders of cyclic factors),
    computed from element orders counts (sufficient for p-groups via counts of elements killed by k)."""
    # number of elements killed by k for each k | d determines the group up to isomorphism
    divs = [k for k in range(1, d + 1) if d % k == 0]
    return tuple(sum(1 for v in S if all((k * t) % d == 0 for t in v)) for k in divs)


def free_type(d, m):
    """the counts above for Z_d^m"""
    divs = [k for k in range(1, d + 1) if d % k == 0]
    return tuple(math.gcd(k, d) ** m for k in divs)


rng = random.Random(7)

# special one-qudit seeds with non-free stabilizer groups
def seed_state(d):
    """one-qudit stabilizer states: |0>, and uniform superpositions over the subgroup aZ_d for a | d."""
    out = []
    for a in range(1, d + 1):
        if d % a == 0:
            v = np.zeros(d, complex)
            v[[k for k in range(0, d, a)]] = 1
            out.append(v / np.linalg.norm(v))
    return out  # a = d gives |0>, a = 1 gives |+>


def subsets(nsys):
    for r in range(nsys + 1):
        for X in itertools.combinations(range(nsys), r):
            yield frozenset(X)


def run_case(d, m, nsys, depth, ntrials, docss=True):
    worst = {'A1a': 0, 'A1b': 0, 'lift_proj': 0, 'lift_char': 0, 'selfdual': 0, 'lam2S': 0, 'cor1': 0,
             'polymatroid': 0}
    nonfree = 0
    seeds = seed_state(d)
    for trial in range(ntrials):
        # product of random one-qudit seeds, then a random Clifford circuit
        psi = np.array([1.0 + 0j])
        for j in range(m):
            psi = np.kron(psi, rng.choice(seeds))
        psi = random_clifford(psi, d, m, rng, depth)
        owner = [rng.randrange(nsys) for _ in range(m)]
        S = entropy_vector(psi, d, m, owner, nsys)
        L = lagrangian_from_state(psi, d, m)
        assert len(L) == d ** m, (len(L), d ** m)
        Lg = generators(L, d, 2 * m)
        assert all(symp(v, w, m, d) == 0 for v in Lg for w in Lg)
        if isomorphism_type(L, d, 2 * m) != free_type(d, m):
            nonfree += 1
        logd = math.log2(d)
        cols = {i: [j for j in range(m) if owner[j] == i] for i in range(nsys)}
        coords = lambda X: sorted([j for i in X for j in cols[i]] + [m + j for i in X for j in cols[i]])
        # Lemma A.1
        for X in subsets(nsys):
            mX = sum(len(cols[i]) for i in X)
            cX = coords(X)
            LX = [v for v in L if all(v[c] == 0 for c in range(2 * m) if c not in cX)]
            f1 = mX * logd - math.log2(len(LX))
            f2 = math.log2(len(proj(L, cX))) - mX * logd
            worst['A1a'] = max(worst['A1a'], abs(S[X] - f1))
            worst['A1b'] = max(worst['A1b'], abs(S[X] - f2))
        # lift
        hS = {X: S[X] + sum(S[frozenset([i])] for i in X) for X in subsets(nsys)}
        # polymatroid check of h_S
        for X in subsets(nsys):
            for i in range(nsys):
                if i not in X:
                    worst['polymatroid'] = max(worst['polymatroid'], max(0, hS[X] - hS[X | {i}]))
                    for j in range(i + 1, nsys):
                        if j not in X:
                            v = hS[X | {i}] + hS[X | {j}] - hS[X | {i, j}] - hS[X]
                            worst['polymatroid'] = max(worst['polymatroid'], max(0, -v))
        # Lemma A.4 intermediate: h_S(X) = log|pr_X L| - sum_{i in X} log|L_i| = log|pr_X Lambda|
        Li = {}
        for i in range(nsys):
            ci = coords([i])
            Li[i] = {tuple(v[c] for c in ci) for v in L if all(v[c] == 0 for c in range(2 * m) if c not in ci)}
        # Lambda: image of L in prod_i pr_i L / L_i; represent a coset by its sorted tuple of elements (canonical)
        def coset_key(i, w):
            ci = coords([i])
            return min(tuple((a + b) % d for a, b in zip(w, l)) for l in Li[i])
        Lam = set()
        for v in L:
            Lam.add(tuple(coset_key(i, tuple(v[c] for c in coords([i]))) for i in range(nsys)))
        for X in subsets(nsys):
            prL = {tuple(t[i] for i in sorted(X)) for t in Lam}
            worst['lift_proj'] = max(worst['lift_proj'], abs(hS[X] - math.log2(len(prL))))
        # explicit character subgroups U_i (functions on L evaluated on generators of L), Lemmas A.3/A.4
        Ui = {}
        for i in range(nsys):
            ci = coords([i])
            if not ci:
                Ui[i] = []
                continue
            Ai = annihilator_dot(Li[i], d, len(ci))
            Aig = generators(set(Ai), d, len(ci))
            Ui[i] = [tuple(sum(a[k] * g[c] for k, c in enumerate(ci)) % d for g in Lg) for a in Aig]
        for X in subsets(nsys):
            gens = [u for i in X for u in Ui[i]]
            size = len(span(gens, d, len(Lg)))
            worst['lift_char'] = max(worst['lift_char'], abs(hS[X] - math.log2(size)))
        # Lemma A.2: h_S self-dual, lambda(h_S) = 2S
        N = frozenset(range(nsys))
        for X in subsets(nsys):
            hstar = sum(hS[frozenset([i])] for i in X) + hS[N - X] - hS[N]
            worst['selfdual'] = max(worst['selfdual'], abs(hstar - hS[X]))
            lam = hS[X] + hS[N - X] - hS[N]
            worst['lam2S'] = max(worst['lam2S'], abs(lam - 2 * S[X]))
        # Corollary 1: CSS state of L (2m qudits) has entropy vector 2S
        if docss and d ** (2 * m) <= 600000:
            phi = np.zeros(d ** (2 * m), complex)
            for v in L:
                idx = 0
                for t in v:
                    idx = idx * d + t
                phi[idx] = 1
            phi /= np.linalg.norm(phi)
            owner2 = owner + owner  # coordinate j (x part) and m+j (z part) belong to the owner of qudit j
            S2 = entropy_vector(phi, d, 2 * m, owner2, nsys)
            for X in subsets(nsys):
                worst['cor1'] = max(worst['cor1'], abs(S2[X] - 2 * S[X]))
    return worst, nonfree


if __name__ == '__main__':
    # (d, qudits, systems, Clifford depth, trials); d = 6 is squarefree, so its Lagrangians are always free
    cases = [(2, 5, 4, 60, 6), (3, 4, 4, 60, 6), (4, 4, 4, 60, 8), (6, 3, 4, 60, 6), (8, 3, 3, 40, 5),
             (9, 3, 3, 40, 4), (12, 2, 3, 30, 5), (4, 3, 3, 40, 8), (6, 4, 3, 50, 3)]
    allok = True
    nonfree_by_d = {}
    for (d, m, nsys, depth, nt) in cases:
        w, nf = run_case(d, m, nsys, depth, nt)
        ok = all(v < 1e-8 for v in w.values())
        allok &= ok
        nonfree_by_d[d] = nonfree_by_d.get(d, 0) + nf
        print('d=%d m=%d systems=%d trials=%d non-free L: %d ' % (d, m, nsys, nt, nf),
              ' '.join('%s:%.1e' % (k, v) for k, v in w.items()), 'OK' if ok else 'FAIL')
    # non-free Lagrangians occur for d = 4, 8, 9, 12 (divisible by a square) and never for squarefree d
    sq = all(nonfree_by_d[d] > 0 for d in (4, 8, 9, 12)) and all(nonfree_by_d[d] == 0 for d in (2, 3, 6))
    print('non-free Lagrangians by d:', nonfree_by_d, '(expected > 0 for 4, 8, 9, 12 and 0 for 2, 3, 6):', sq)
    allok &= sq
    print('ALL OK' if allok else 'FAILURE')
