"""Count the single-CI inequality classes (data/s17_ineq_classes.npy) that are facets of the qubit stabilizer
cone, using only qubit vectors recomputed here: all six-vertex graph states and the S6 images of every stored qubit
graph-state certificate (762 in data/s16_realisations.npz, 11 in results/2026-10-02_a2-ingleton/a2i_q2.npz, 16 at
lambda = 2 in data/s24_qubit_l2_certs.npz). Also check that no class is violated by any of these vectors.
With --with-literature, the 12 known extreme rays without stored qubit certificates (rows 6, 7, 15, 34, 39, 40, 42,
43, 44, 46, 48, 50 of data/stab5_extreme_reps.npy) are added to the pool; they lie in the qubit cone through the
hypergraph realizations of arXiv:2006.16292 and the theorem of Walter and Witteveen (J. Math. Phys. 62, 092203)."""
import itertools, sys
import numpy as np
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from common import REPO, graph_state_entropy, masks
def Sfull(S, M):
    if M == 0 or M == 63: return 0
    if M >> 5 & 1: M = 63 ^ M
    return S[M - 1]
def act(S, p):
    inv = [0] * 6
    for i, j in enumerate(p): inv[j] = i
    return tuple(Sfull(S, sum(1 << inv[i] for i in range(6) if M >> i & 1)) for M in range(1, 32))
perms = list(itertools.permutations(range(6)))
pool = set()
pairs = list(itertools.combinations(range(6), 2))
for g in range(1 << 15):
    adj = [[0] * 6 for _ in range(6)]
    for k, (a, b) in enumerate(pairs):
        if g >> k & 1: adj[a][b] = adj[b][a] = 1
    S = graph_state_entropy(adj, list(range(6)), 6, 2)
    pool.add(tuple(S[M - 1] for M in masks(5)))
z = np.load(REPO + '/data/s16_realisations.npz', allow_pickle=True)
for k in range(len(z['rep'])):
    sizes = [int(s) for s in z['sizes'][k]]; nv = sum(sizes)
    rows = [int(x) for x in z['rows'][k][:nv]]
    adj = [[(rows[a] >> b) & 1 for b in range(nv)] for a in range(nv)]
    owner = [p for p, s in enumerate(sizes) for _ in range(s)]
    S = graph_state_entropy(adj, owner, 6, 2)
    lam = int(z['lam'][k]); rep = [int(x) for x in z['rep'][k]]
    assert [S[M - 1] for M in masks(5)] == [lam * x for x in rep]
    for p in perms:
        pool.add(act(rep, p))

# R30 qubit certificates (results/2026-10-02_a2-ingleton/a2i_q2.npz) and the 16 qubit certificates at lambda = 2
# (data/s24_qubit_l2_certs.npz), recomputed here as well
def add_cert(rep, lam, sizes, adj):
    owner = [p for p, s in enumerate(sizes) for _ in range(s)]
    S = graph_state_entropy(adj, owner, 6, 2)
    assert [S[M - 1] for M in masks(5)] == [lam * x for x in rep]
    for p in perms:
        pool.add(act(rep, p))
z = np.load(REPO + '/results/2026-10-02_a2-ingleton/a2i_q2.npz', allow_pickle=True)
for k in range(len(z['rep'])):
    sizes = [int(s) for s in z['sizes'][k]]; nv = sum(sizes)
    rows = [int(x) for x in z['rows'][k][:nv]]
    add_cert([int(x) for x in z['rep'][k]], int(z['lam'][k]), sizes, [[(rows[a] >> b) & 1 for b in range(nv)] for a in range(nv)])
z = np.load(REPO + '/data/s24_qubit_l2_certs.npz', allow_pickle=True)
for k in range(len(z['rep'])):
    sizes = [int(s) for s in z['sizes'][k]]; nv = sum(sizes)
    add_cert([int(x) for x in z['rep'][k]], int(z['lam'][k]), sizes, z['G'][k][:nv, :nv].astype(int).tolist())
# optional: the 12 known extreme rays that lie in the qubit cone via hypergraph realizations (arXiv:2006.16292) and
# Walter-Witteveen; no explicit qubit realization of them is stored
if '--with-literature' in sys.argv:
    E12 = np.load(REPO + '/data/stab5_extreme_reps.npy').astype(np.int64)
    for k in [6, 7, 15, 34, 39, 40, 42, 43, 44, 46, 48, 50]:
        for p in perms:
            pool.add(act([int(x) for x in E12[k]], p))
from collections import Counter
lam_counts = Counter()
for f in ['/data/s16_realisations.npz', '/results/2026-10-02_a2-ingleton/a2i_q2.npz', '/data/s24_qubit_l2_certs.npz']:
    lam_counts.update(int(x) for x in np.load(REPO + f, allow_pickle=True)['lam'])
print('stored qubit realizations by multiplicity lambda:', dict(sorted(lam_counts.items())),
      ' total:', sum(lam_counts.values()))
P = np.array(sorted(pool), dtype=np.int64)
print('pool:', P.shape)
C = np.load(REPO + '/data/s17_ineq_classes.npy').astype(np.int64)
V = P @ C.T
print('minimum over pool and classes:', int(V.min()))
PR = 1000003
def rank_mod(T):
    """rank over GF(PR) of an integer matrix (rows), numpy int64 elimination"""
    M = np.array(T, dtype=np.int64) % PR
    r = 0
    rows, cols = M.shape
    for c in range(cols):
        piv = None
        for i in range(r, rows):
            if M[i, c]:
                piv = i; break
        if piv is None:
            continue
        M[[r, piv]] = M[[piv, r]]
        inv = pow(int(M[r, c]), PR - 2, PR)
        M[r] = (M[r] * inv) % PR
        nz = np.nonzero(M[:, c])[0]
        for i in nz:
            if i != r:
                M[i] = (M[i] - M[i, c] * M[r]) % PR
        r += 1
        if r == rows:
            break
    return r
nfacet = 0
for j in range(C.shape[0]):
    T = P[V[:, j] == 0]
    if len(T) < 30:
        continue
    if np.linalg.matrix_rank(T.astype(float)) < 30:
        continue
    # exact certificate: rank over GF(PR) is a lower bound for the rational rank, which is at most 30
    if rank_mod(T[:3000]) == 30 or rank_mod(T) == 30:
        nfacet += 1
print('classes certified as facets of the qubit cone (exact, rank mod 1000003 = 30):', nfacet)
