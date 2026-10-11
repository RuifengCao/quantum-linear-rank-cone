"""Six-variable linear rank inequalities of Dougherty, Freiling and Zeger (arXiv:0910.0284v3, Sec. 6, Eqs. (44)-(63))
applied to five-party stabilizer states with the purifier as the sixth variable.

1. Transcription: each inequality is entered term by term from the rendered pages 24-26 of arXiv:0910.0284v3.
2. Sanity of the transcription: every inequality is balanced, and it is never violated on random subspace
   arrangements of six subspaces over GF(2), GF(3) and GF(5).
3. Pull-back: a balanced six-variable form f gives the linear form S -> f(S) on pure-symmetric vectors of the five
   parties A..E and the purifier F, for each of the 720 bijections between the six variables and the six systems.
4. Consistency: no pulled-back form is violated by an entropy vector of a six-qubit graph state, nor by any of the
   641 realised extreme-ray orbits of the catalogue (catalogue status 1).
5. Evaluation on the catalogue of 1,917,706 extreme-ray orbits of QLR_5, on the 300 frontier rays, and on the
   witness r of Theorem 3 (frontier ray 275).
"""
import itertools, os, sys, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import REPO, rank_mod_p, graph_state_entropy, masks

V = 'ABCDEF'
# (coefficient, X, Y, Z) means coefficient * I(X;Y|Z); an inequality is (lhs terms, rhs terms): sum lhs <= sum rhs
DFZ6 = {
    44: ([(1, 'A', 'B', '')], [(1, 'A', 'C', ''), (1, 'B', 'D', 'C'), (1, 'A', 'E', 'D'), (1, 'B', 'F', 'E'), (1, 'A', 'B', 'F')]),
    45: ([(1, 'A', 'B', '')], [(1, 'A', 'C', ''), (1, 'B', 'D', 'C'), (1, 'A', 'E', 'D'), (1, 'A', 'F', 'E'), (1, 'A', 'B', 'F')]),
    46: ([(1, 'A', 'B', '')], [(1, 'A', 'C', ''), (1, 'B', 'D', 'C'), (1, 'E', 'F', 'D'), (1, 'A', 'B', 'E'), (1, 'A', 'B', 'F')]),
    47: ([(1, 'A', 'B', '')], [(1, 'A', 'C', ''), (1, 'D', 'E', 'C'), (1, 'A', 'B', 'D'), (1, 'B', 'F', 'E'), (1, 'A', 'B', 'F')]),
    48: ([(1, 'A', 'B', '')], [(1, 'C', 'D', ''), (1, 'A', 'B', 'C'), (1, 'E', 'F', 'D'), (1, 'A', 'B', 'E'), (1, 'A', 'B', 'F')]),
    49: ([(2, 'A', 'B', '')], [(1, 'A', 'C', ''), (1, 'D', 'EF', 'C'), (1, 'A', 'B', 'D'), (1, 'E', 'F', ''), (1, 'A', 'B', 'E'),
                              (1, 'A', 'B', 'F')]),
    50: ([(2, 'A', 'B', '')], [(1, 'A', 'C', ''), (1, 'B', 'D', 'C'), (1, 'A', 'EF', 'D'), (1, 'E', 'F', ''), (1, 'A', 'B', 'E'),
                              (1, 'A', 'B', 'F')]),
    51: ([(2, 'A', 'B', '')], [(1, 'C', 'D', ''), (1, 'A', 'B', 'C'), (1, 'B', 'EF', 'D'), (1, 'E', 'F', ''), (1, 'A', 'B', 'E'),
                              (1, 'A', 'B', 'F')]),
    52: ([(2, 'A', 'B', '')], [(1, 'CD', 'E', ''), (1, 'C', 'D', ''), (1, 'A', 'F', 'C'), (1, 'A', 'B', 'F'), (1, 'A', 'B', 'D'),
                              (1, 'A', 'B', 'E')]),
    53: ([(3, 'A', 'B', '')], [(1, 'CD', 'EF', ''), (1, 'C', 'D', ''), (1, 'E', 'F', ''), (1, 'A', 'B', 'C'), (1, 'A', 'B', 'D'),
                              (1, 'A', 'B', 'E'), (1, 'A', 'B', 'F')]),
    54: ([(1, 'A', 'BC', '')], [(1, 'D', 'E', ''), (1, 'C', 'F', 'D'), (1, 'A', 'B', 'DF'), (1, 'A', 'B', 'CD'), (1, 'A', 'C', 'BF'),
                               (1, 'A', 'BC', 'E')]),
    55: ([(1, 'AB', 'CD', '')], [(1, 'A', 'CD', ''), (1, 'B', 'E', 'A'), (1, 'B', 'D', 'ACF'), (1, 'D', 'F', 'AE'), (1, 'B', 'C', 'AEF'),
                                (1, 'B', 'C', 'DE'), (1, 'A', 'D', 'BCF'), (1, 'A', 'C', 'BEF'), (1, 'A', 'F', 'BDE')]),
    56: ([(2, 'A', 'B', '')], [(1, 'D', 'F', ''), (1, 'A', 'C', ''), (1, 'B', 'D', 'C'), (1, 'A', 'B', 'F'), (1, 'A', 'E', 'D'),
                              (1, 'A', 'F', 'CD'), (1, 'A', 'B', 'E')]),
    57: ([(1, 'A', 'BC', '')], [(1, 'A', 'C', ''), (1, 'B', 'D', 'C'), (1, 'A', 'F', 'D'), (1, 'A', 'B', 'F'), (1, 'C', 'E', 'BF'),
                               (1, 'A', 'C', 'BE')]),
    58: ([(3, 'AB', 'CDE', '')], [(1, 'A', 'CF', ''), (1, 'AB', 'D', ''), (1, 'AB', 'E', ''), (1, 'C', 'F', 'D'), (1, 'D', 'F', 'E'),
                                 (1, 'A', 'E', 'DF'), (1, 'B', 'C', 'ADF'), (1, 'B', 'D', 'CF'), (1, 'A', 'DE', 'BC'),
                                 (1, 'A', 'D', 'BCE'), (1, 'A', 'C', 'EF'), (1, 'B', 'D', 'AEF'), (1, 'B', 'CD', 'A'),
                                 (1, 'AB', 'E', 'CD'), (1, 'B', 'E', 'ACD'), (1, 'B', 'D', 'CEF'), (1, 'AB', 'C', 'DE')]),
    59: ([(1, 'A', 'B', ''), (1, 'A', 'C', '')], [(1, 'B', 'C', ''), (1, 'A', 'D', ''), (1, 'B', 'E', 'D'), (1, 'C', 'F', 'D'),
                                                (1, 'A', 'B', 'E'), (1, 'A', 'C', 'F')]),
    60: ([(2, 'A', 'BC', ''), (1, 'B', 'CD', '')], [(1, 'A', 'CE', ''), (1, 'A', 'F', ''), (1, 'A', 'C', 'D'), (2, 'A', 'B', 'CF'),
                                                  (1, 'B', 'C', ''), (1, 'E', 'F', 'C'), (2, 'B', 'D', 'CE'), (1, 'C', 'E', 'F'),
                                                  (1, 'A', 'D', 'EF'), (1, 'D', 'E', 'ACF'), (2, 'A', 'F', 'CDE')]),
    61: ([(1, 'A', 'B', '')], [(1, 'A', 'C', ''), (1, 'B', 'D', 'C'), (1, 'E', 'F', 'D'), (1, 'A', 'B', 'E'), (1, 'A', 'C', 'F'),
                              (1, 'B', 'E', 'CF')]),
    62: ([(2, 'AB', 'CDE', '')], [(1, 'AB', 'DE', ''), (1, 'ADF', 'C', ''), (1, 'AF', 'D', 'C'), (1, 'B', 'C', 'DE'), (1, 'A', 'C', 'B'),
                                 (1, 'A', 'D', 'BCE'), (2, 'A', 'C', 'DEF'), (1, 'B', 'C', 'ADE'), (1, 'A', 'E', 'BDF'),
                                 (1, 'B', 'E', 'ACF'), (1, 'B', 'E', 'ADF'), (1, 'B', 'E', 'CD'), (1, 'B', 'D', 'AEF'),
                                 (1, 'A', 'F', 'BDE'), (1, 'A', 'F', 'BCD')]),
    63: ([(2, 'A', 'BC', '')], [(1, 'A', 'B', ''), (1, 'D', 'E', ''), (1, 'A', 'B', 'C'), (1, 'C', 'E', 'B'), (1, 'D', 'F', 'BE'),
                               (1, 'C', 'F', 'D'), (1, 'A', 'B', 'CD'), (1, 'A', 'BC', 'F'), (1, 'A', 'C', 'E')]),
}


def mk(s):
    return sum(1 << V.index(c) for c in s)


def form(ineq):
    """vector f over the 64 subsets of A..F with f(h) = rhs - lhs >= 0"""
    f = np.zeros(64, dtype=np.int64)
    lhs, rhs = ineq
    for sign, terms in ((-1, lhs), (1, rhs)):
        for c, X, Y, Z in terms:
            x, y, z = mk(X), mk(Y), mk(Z)
            for M, s in ((x | z, 1), (y | z, 1), (x | y | z, -1), (z, -1)):
                if M:
                    f[M] += sign * c * s
    return f


forms = {k: form(v) for k, v in DFZ6.items()}
ok = True
bal = all(sum(int(f[M]) for M in range(64) if M >> i & 1) == 0 for f in forms.values() for i in range(6))
print('all 20 transcribed inequalities balanced:', bal)
ok &= bal

# ------------------------------------------------------------------ 2. random subspace arrangements
rng = np.random.default_rng(20261003)


def rank_function(subspaces, p):
    h = np.zeros(64, dtype=np.int64)
    for M in range(1, 64):
        rows = [v for i in range(6) if M >> i & 1 for v in subspaces[i]]
        h[M] = rank_mod_p(rows, p) if rows else 0
    return h


F = np.array([forms[k] for k in sorted(forms)])
worst = {k: 10 ** 9 for k in forms}
ntest = 0
t0 = time.time()
for trial in range(6000):
    p = (2, 3, 5)[trial % 3]
    k = int(rng.integers(2, 6))
    pool = [[int(x) for x in rng.integers(0, p, size=k)] for _ in range(int(rng.integers(k, k + 4)))]
    subspaces = []
    for i in range(6):
        size = int(rng.integers(0, 3))
        subspaces.append([pool[int(j)] for j in rng.choice(len(pool), size=size, replace=False)] if size else [])
    h = rank_function(subspaces, p)
    vals = F @ h
    for idx, kk in enumerate(sorted(forms)):
        worst[kk] = min(worst[kk], int(vals[idx]))
    ntest += 1
print('random arrangements tested: %d (%.0f s); minimum value of each inequality:' % (ntest, time.time() - t0), worst)
ok &= all(v >= 0 for v in worst.values())

# ------------------------------------------------------------------ 3. pull-back to five parties + purifier
def coord(M):
    """index into the 31 coordinates S(X), X a nonempty subset of A..E, of the set M of the six systems (purity)"""
    if M == 0 or M == 63:
        return None
    if M >> 5 & 1:
        M = 63 ^ M
    return M - 1


perms = list(itertools.permutations(range(6)))
rows, labels = [], []
for kk in sorted(forms):
    f = forms[kk]
    seen = set()
    for pi in perms:
        g = np.zeros(31, dtype=np.int64)
        for M in range(1, 64):
            c = int(f[M])
            if c:
                PM = sum(1 << pi[i] for i in range(6) if M >> i & 1)
                j = coord(PM)
                if j is not None:
                    g[j] += c
        t = tuple(g.tolist())
        if t not in seen and any(t):
            seen.add(t)
            rows.append(g)
            labels.append(kk)
G = np.array(rows, dtype=np.int64)
labels = np.array(labels)
print('pulled-back forms (distinct per inequality):', {kk: int((labels == kk).sum()) for kk in sorted(forms)}, ' total', len(G))

# ------------------------------------------------------------------ 4. consistency on realised vectors
pairs = list(itertools.combinations(range(6), 2))
gs = set()
for code in range(1 << 15):
    adj = [[0] * 6 for _ in range(6)]
    for e, (a, b) in enumerate(pairs):
        if code >> e & 1:
            adj[a][b] = adj[b][a] = 1
    S = graph_state_entropy(adj, list(range(6)), 6, 2)
    gs.add(tuple(S[M - 1] for M in masks(5)))
GS = np.array(sorted(gs), dtype=np.int64)
v_gs = int((GS @ G.T).min())
print('six-qubit graph-state vectors: %d, minimum over all pulled-back forms: %d' % (len(GS), v_gs))
ok &= v_gs >= 0

CAT = np.load(REPO + '/results/2026-09-12_a1-campaign/qlr_adj.npz')['reps'].astype(np.float32)
status = np.load(REPO + '/results/2026-10-02_a2-status/catalogue_status.npz')['status']
Gf = G.T.astype(np.float32)
t0 = time.time()
minval = np.empty(CAT.shape[0], dtype=np.float32)
argmin = np.empty(CAT.shape[0], dtype=np.int64)
for s in range(0, CAT.shape[0], 4000):
    V_ = CAT[s:s + 4000] @ Gf
    argmin[s:s + 4000] = V_.argmin(axis=1)
    minval[s:s + 4000] = V_[np.arange(V_.shape[0]), argmin[s:s + 4000]]
print('catalogue evaluated in %.0f s' % (time.time() - t0))
cut = minval < 0
print('catalogue orbits: %d; cut by at least one of (44)-(63): %d' % (CAT.shape[0], int(cut.sum())))
by = {}
for kk in sorted(forms):
    sel = labels == kk
    m_k = np.full(CAT.shape[0], np.inf, dtype=np.float32)
    for s in range(0, CAT.shape[0], 4000):
        m_k[s:s + 4000] = (CAT[s:s + 4000] @ Gf[:, sel]).min(axis=1)
    by[kk] = m_k
cut4448 = np.zeros(CAT.shape[0], dtype=bool)
for kk in range(44, 49):
    cut4448 |= by[kk] < 0
print('cut by (44)-(48): %d' % int(cut4448.sum()))
print('cut per inequality:', {kk: int((by[kk] < 0).sum()) for kk in sorted(forms)})
print('status of cut orbits (0 excluded by CI, 1 realised, 2 undecided, 3 excluded with Ingleton):',
      {c: int(((status == c) & cut).sum()) for c in (0, 1, 2, 3)})
ok &= int(((status == 1) & cut).sum()) == 0
print('undecided orbits newly cut:', int(((status == 2) & cut).sum()))

FR = np.load(REPO + '/data/a2_frontier300.npy').astype(np.int64)
VF = FR @ G.T
fr_cut = np.flatnonzero(VF.min(axis=1) < 0)
print('frontier rays cut:', fr_cut.tolist(), ' by:', {int(i): sorted(set(labels[VF[i] < 0].tolist())) for i in fr_cut})
r = FR[275]
print('witness r (frontier 275):', r.tolist(), ' minimum over all pulled-back forms:', int((G @ r).min()))
print('ALL OK' if ok else 'FAILURE')
