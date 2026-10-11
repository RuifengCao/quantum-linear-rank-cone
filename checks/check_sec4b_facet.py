"""Facet check of (star) on the qubit stabilizer cone using only independently recomputed qubit vectors:
(a) all 2^15 six-vertex qubit graph states (one qubit per system), (b) S6 images of qubit graph-state certificates
from data/s16_realisations.npz and data/s24_qubit_l2_certs.npz, each re-verified here by computing cut ranks."""
import itertools, sys
import numpy as np
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from common import REPO, rank_Q, graph_state_entropy, masks
L = 'ABCDEF'
def m(s): return sum(1 << L.index(c) for c in s)
star = [0] * 31
for s in ['BCD', 'ABCD', 'ABE', 'BCE', 'BCE', 'ADE', 'ADE', 'CDE']: star[m(s) - 1] += 1
for s in ['BC', 'BC', 'ABD', 'ACD', 'AE', 'DE', 'BCDE', 'ABCDE']: star[m(s) - 1] -= 1
star = np.array(star)
def Sfull(S, M):
    if M == 0 or M == 63: return 0
    if M >> 5 & 1: M = 63 ^ M
    return S[M - 1]
def act(S, p):
    inv = [0] * 6
    for i, j in enumerate(p): inv[j] = i
    return tuple(Sfull(S, sum(1 << inv[i] for i in range(6) if M >> i & 1)) for M in range(1, 32))
perms = list(itertools.permutations(range(6)))
# (a) six-vertex graph states
pairs = list(itertools.combinations(range(6), 2))
vecs = set()
for g in range(1 << 15):
    adj = [[0] * 6 for _ in range(6)]
    for k, (a, b) in enumerate(pairs):
        if g >> k & 1: adj[a][b] = adj[b][a] = 1
    S = graph_state_entropy(adj, list(range(6)), 6, 2)
    vecs.add(tuple(S[M - 1] for M in masks(5)))
vecs = [list(v) for v in vecs]
vals = np.array(vecs) @ star
tight = [v for v, x in zip(vecs, vals) if x == 0]
print('(a) six-qubit graph states:', len(vecs), ' min(star):', int(vals.min()), ' tight:', len(tight), ' rank:', rank_Q(tight))
print('(a) rank of all six-qubit graph-state vectors (31 = full dimension):', rank_Q(vecs))
# (b) qubit certificates, re-verified
pool = []
z = np.load(REPO + '/data/s16_realisations.npz', allow_pickle=True)
nver = 0
for k in range(len(z['rep'])):
    sizes = [int(s) for s in z['sizes'][k]]; nv = sum(sizes)
    rows = [int(x) for x in z['rows'][k][:nv]]
    adj = [[(rows[a] >> b) & 1 for b in range(nv)] for a in range(nv)]
    owner = [p for p, s in enumerate(sizes) for _ in range(s)]
    S = graph_state_entropy(adj, owner, 6, 2)
    lam = int(z['lam'][k]); rep = [int(x) for x in z['rep'][k]]
    if [S[M - 1] for M in masks(5)] == [lam * x for x in rep]:
        nver += 1; pool.append(rep)
print('(b) s16 certificates re-verified:', nver, 'of', len(z['rep']))
z2 = np.load(REPO + '/data/s24_qubit_l2_certs.npz', allow_pickle=True)
print('   s24 keys:', z2.files)
tightb = []
allb = []
for rep in pool:
    for p in perms:
        v = act(rep, p)
        allb.append(v)
allb = np.array(sorted(set(allb)))
valsb = allb @ star
print('(b) S6 images:', len(allb), ' min(star):', int(valsb.min()))
T = tight + [list(v) for v, x in zip(allb, valsb) if x == 0]
print('tight vectors in (a)+(b):', len(T), ' rank:', rank_Q(T))
