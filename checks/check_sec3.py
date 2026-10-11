"""Checks for Section III / Appendix A of the paper (independent of the repository code).

1. q2 (row 195 of data/stab5_extreme_reps.npy): 2*q2 is the entropy vector of the stored 12-qubit graph state.
2. q2 itself is not the cut-rank vector of any 6-vertex qubit graph state (all 2^15 graphs); 2*q2 is the entropy
   vector of the stored 12-qutrit weighted graph state, and q2 is not the cut-rank vector of any 6-vertex weighted
   graph over GF(3) (all 3^15 graphs, C helper q2_gf3.c).
3. CSS entropy formula S = lambda(r_C) and the doubling 2S of Corollary 1, against state vectors, over Z_d for
   d = 2, 3, 4, 6 (composite d included).
"""
import itertools, math, os, random, subprocess, sys, tempfile
import numpy as np
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from common import REPO, masks, rank_mod_p, graph_state_entropy

np.set_printoptions(linewidth=200)
ok = True

# ---------------------------------------------------------------- 1. q2 and its 12-qubit realisation
R = np.load(REPO + '/data/stab5_extreme_reps.npy')
q2 = [int(x) for x in R[195]]
print('q2 =', q2)
z = np.load(REPO + '/data/s16_realisations.npz', allow_pickle=True)
i = [k for k, l in enumerate(z['label']) if 'q2' in str(l)][0]
rows = [int(x) for x in z['rows'][i][:12]]
adj = [[(rows[a] >> b) & 1 for b in range(12)] for a in range(12)]
assert all(adj[a][b] == adj[b][a] for a in range(12) for b in range(12)) and all(adj[a][a] == 0 for a in range(12))
owner = [v // 2 for v in range(12)]          # A..F own two qubits each, F = party 5 = purifier
S = graph_state_entropy(adj, owner, 6, 2)     # over all 63 masks of the six parties
# restrict to masks of A..E (bits 0..4) and compare with 2*q2
S5 = [S[m - 1] for m in masks(5)]
print('12-qubit graph state, entropies on A..E:', S5)
c1 = S5 == [2 * x for x in q2]
# purity check: S(X) = S(complement within six parties)
c1p = all(S[m - 1] == S[(63 ^ m) - 1] for m in masks(6) if m != 63)
print('check 1 (2*q2 realised, pure-symmetric):', c1 and c1p)
ok &= c1 and c1p
print('number of edges:', sum(adj[a][b] for a in range(12) for b in range(a + 1, 12)))

# ---------------------------------------------------------------- 2. q2 at one qubit per party
def gf2_rank_bits(rowbits):
    r = 0
    rows_ = list(rowbits)
    while rows_:
        piv = rows_.pop()
        if piv == 0:
            continue
        r += 1
        low = piv & -piv
        rows_ = [x ^ piv if x & low else x for x in rows_]
    return r

pairs = list(itertools.combinations(range(6), 2))
target = tuple(q2)
hits = 0
cutvecs = set()
for g in range(1 << 15):
    A = [0] * 6
    for k, (a, b) in enumerate(pairs):
        if g >> k & 1:
            A[a] |= 1 << b
            A[b] |= 1 << a
    vec = []
    for m in masks(5):
        X = [v for v in range(6) if m >> v & 1]
        comp = 63 ^ m
        vec.append(gf2_rank_bits([A[v] & comp for v in X]))
    t = tuple(vec)
    cutvecs.add(t)
    if t == target:
        hits += 1
print('check 2: 2^15 six-vertex qubit graph states, distinct entropy vectors:', len(cutvecs), ' hits on q2:', hits)
ok &= hits == 0

# ---------------------------------------------------------------- 2b. q2 with qutrits
# (i) 2*q2 is the entropy vector of the stored 12-qutrit weighted graph state (data/s24_qutrit_certs.npz)
zq = np.load(REPO + '/data/s24_qutrit_certs.npz', allow_pickle=True)
kq = zq['stab5_row'].tolist().index(195)
szq = [int(x) for x in zq['sizes'][kq]]; nvq = sum(szq)
Gq = zq['G'][kq][:nvq, :nvq].astype(int).tolist()
ownq = [p_ for p_, s_ in enumerate(szq) for _ in range(s_)]
Sq = graph_state_entropy(Gq, ownq, 6, 3)
c1q = [Sq[M - 1] for M in masks(5)] == [2 * x for x in q2] and int(zq['p'][kq]) == 3 and nvq == 12
print('check 2b(i): 2*q2 realised by a 12-qutrit weighted graph state (GF(3) cut ranks):', c1q)
ok &= c1q
# (ii) no weighted six-vertex graph over GF(3) has cut-rank vector q2 (exhaustive, C helper q2_gf3.c);
#      sanity: the helper finds the cut-rank vector of a known weighted graph
here = os.path.dirname(os.path.abspath(__file__))
exe = os.path.join(tempfile.mkdtemp(), 'q2_gf3')
subprocess.run(['gcc', '-O2', '-o', exe, os.path.join(here, 'q2_gf3.c')], check=True)
def run_gf3(vec):
    out = subprocess.run([exe], input=' '.join(str(int(x)) for x in vec), capture_output=True, text=True, check=True).stdout.split()
    return int(out[0]), int(out[1])
Wt = [[0, 1, 2, 0, 1, 0], [1, 0, 1, 2, 0, 0], [2, 1, 0, 1, 0, 1], [0, 2, 1, 0, 2, 1], [1, 0, 0, 2, 0, 2], [0, 0, 1, 1, 2, 0]]
St = graph_state_entropy(Wt, list(range(6)), 6, 3)
tot, hit_t = run_gf3([St[M - 1] for M in masks(5)])
tot, hit_q2 = run_gf3(q2)
print('check 2b(ii): %d weighted six-vertex graphs over GF(3); sanity vector found %d times; hits on q2: %d' % (tot, hit_t, hit_q2))
ok &= tot == 3 ** 15 and hit_t >= 1 and hit_q2 == 0

# ---------------------------------------------------------------- 3. CSS formula and doubling, state vectors over Z_d
def entropy_bits(psi, d, m, X):
    """von Neumann entropy (bits) of the reduced state on qudit set X of a pure state psi (shape d^m)."""
    psi = psi.reshape([d] * m)
    Xs = sorted(X)
    Ys = [j for j in range(m) if j not in X]
    M = np.transpose(psi, Xs + Ys).reshape(d ** len(Xs), d ** len(Ys))
    s = np.linalg.svd(M, compute_uv=False) ** 2
    s = s[s > 1e-12]
    return float(-(s * np.log2(s)).sum())


def subgroup(gens, d, m):
    """subgroup of Z_d^m generated by gens (as a set of tuples)."""
    G = {tuple([0] * m)}
    frontier = list(G)
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


def css_state(C, d, m):
    psi = np.zeros(d ** m)
    for c in C:
        idx = 0
        for x in c:
            idx = idx * d + x
        psi[idx] = 1.0
    return psi / np.linalg.norm(psi)


def proj_size(C, cols):
    return len({tuple(c[j] for j in cols) for c in C})


rng = random.Random(3)
worst = 0.0
cases = 0
for d in (2, 3, 4, 6):
    for trial in range(12):
        m = 4 if d <= 4 else 3
        gens = [tuple(rng.randrange(d) for _ in range(m)) for _ in range(rng.randint(1, 3))]
        C = subgroup(gens, d, m)
        psi = css_state(C, d, m)
        owner = [rng.randrange(3) for _ in range(m)]      # three parties 0, 1, 2 (2 = purifier)
        for X in masks(3):
            if X == 7:
                continue
            Xq = [j for j in range(m) if X >> owner[j] & 1]
            Yq = [j for j in range(m) if not (X >> owner[j] & 1)]
            sv = entropy_bits(psi, d, m, Xq) if Xq and Yq else 0.0
            pred = math.log2(proj_size(C, Xq)) + math.log2(proj_size(C, Yq)) - math.log2(len(C)) if Xq and Yq else 0.0
            worst = max(worst, abs(sv - pred))
            cases += 1
print('check 3a: CSS formula S = lambda(r_C) over Z_2, Z_3, Z_4, Z_6: %d cuts, max deviation %.2e' % (cases, worst))
ok &= worst < 1e-9


def graph_state(W, d, m):
    """Z_d graph state sum_x omega^{x W x / 2 ...}: use the stabilizer state with Lagrangian rowspace [I | W];
    psi(x) = omega^{sum_{a<b} W_ab x_a x_b} for symmetric W with zero diagonal."""
    omega = np.exp(2j * np.pi / d)
    psi = np.zeros(d ** m, dtype=complex)
    for idx, x in enumerate(itertools.product(range(d), repeat=m)):
        e = sum(W[a][b] * x[a] * x[b] for a in range(m) for b in range(a + 1, m))
        psi[idx] = omega ** (e % d)
    return psi / np.linalg.norm(psi)


worst2 = 0.0
cases2 = 0
for d in (2, 3, 4):
    for trial in range(6):
        m = 3
        W = [[0] * m for _ in range(m)]
        for a in range(m):
            for b in range(a + 1, m):
                W[a][b] = W[b][a] = rng.randrange(d)
        psi = graph_state(W, d, m)
        owner = [0, 1, 2]
        # Lagrangian L = rowspace over Z_d of [I | W]; as a classical code of length 2m with qudit j owning
        # coordinates j (x part) and m + j (z part)
        gens = [tuple([1 if k == a else 0 for k in range(m)] + [W[a][k] for k in range(m)]) for a in range(m)]
        L = subgroup(gens, d, 2 * m)
        assert len(L) == d ** m
        phi = css_state(L, d, 2 * m)
        for X in masks(3):
            if X == 7:
                continue
            Xq = [j for j in range(m) if X >> owner[j] & 1]
            S = entropy_bits(psi, d, m, Xq)
            Xq2 = Xq + [m + j for j in Xq]
            S2 = entropy_bits(phi, d, 2 * m, Xq2)
            # entropy formula: S(X) = log|pr_X L| - m_X log d
            pred = math.log2(proj_size(L, Xq2)) - len(Xq) * math.log2(d)
            worst2 = max(worst2, abs(S2 - 2 * S), abs(S - pred))
            cases2 += 1
print('check 3b: doubling (CSS state of L has entropy 2S) and entropy formula over Z_2, Z_3, Z_4: '
      '%d cuts, max deviation %.2e' % (cases2, worst2))
ok &= worst2 < 1e-9

print('ALL OK' if ok else 'FAILURE')
