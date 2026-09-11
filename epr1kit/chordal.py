"""Chordality test and simple-tree construction (Hubeny-Rota, arXiv:2412.18018 /
2512.18702 / 2512.24490) for five-party entropy vectors in mask order.

Conventions: parties 0..4 = A..E, purifier = 5.  For W subset of {0..5}:
S_W = v[mask(W)-1] if 5 not in W (S_empty = 0), else S of the complement.
Theorem 1 (2512.24490): a vector obeying SA and SSA is realizable by a
holographic SIMPLE FOREST graph model iff the line graph L_P of its correlation
hypergraph is chordal.  Theorem 2: for IRREDUCIBLE vectors (no vanishing entry,
Q^cap empty for every max-clique) a simple TREE, constructed by Algorithm 1.
Holographic realizability places the ray inside the stabilizer entropy cone
(Hayden-Nezami-Qi-Thomas-Walter-Yang), so a chordal verdict is S7-positive.
We only ever act on positive verdicts; non-chordal says nothing about
stabilizer realizability by non-tree models.
"""
import itertools
import numpy as np

PARTIES = 6
PUR = 5
ALL = (1 << PARTIES) - 1


def S_of(v, W):
    """entropy of the subset W (bitmask over 6 parties, purifier = bit 5)."""
    if W == 0:
        return 0
    if W & (1 << PUR):
        W = ALL ^ W
        if W == 0:
            return 0
    return int(v[W - 1])


def mi(v, Y, Z):
    return S_of(v, Y) + S_of(v, Z) - S_of(v, Y | Z)


def submasks(X):
    """proper nonempty submasks of X."""
    s = (X - 1) & X
    while s:
        yield s
        s = (s - 1) & X


def sa_ssa_ok(v):
    masks = range(1, ALL)          # nonempty proper subsets of the 6 parties
    for Y in masks:
        for Z in masks:
            if Y & Z or Y >= Z:
                continue
            if mi(v, Y, Z) < 0:
                return False
    # SSA: I(X:YZ) >= I(X:Y) for pairwise disjoint nonempty X,Y,Z
    for X in masks:
        for Y in masks:
            if X & Y or X >= Y:
                pass
            if X & Y:
                continue
            for Z in masks:
                if Z & (X | Y):
                    continue
                if mi(v, X, Y | Z) < mi(v, X, Y):
                    return False
    return True


def positive_beta(v, X):
    """beta(X) positive: every bipartition {Y, X\\Y} of X has I > 0."""
    for Y in submasks(X):
        Z = X ^ Y
        if Y > Z:
            continue
        if mi(v, Y, Z) == 0:
            return False
    return True


def hyperedges(v):
    return [X for X in range(1, ALL + 1) if bin(X).count('1') >= 2 and positive_beta(v, X)]


def line_graph(H):
    n = len(H)
    adj = [set() for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if H[i] & H[j]:
                adj[i].add(j); adj[j].add(i)
    return adj


def is_chordal(adj):
    """maximum cardinality search + perfect elimination check."""
    n = len(adj)
    if n == 0:
        return True
    weight = [0] * n; order = []; numbered = [False] * n
    for _ in range(n):
        u = max((i for i in range(n) if not numbered[i]), key=lambda i: weight[i])
        order.append(u); numbered[u] = True
        for w in adj[u]:
            if not numbered[w]:
                weight[w] += 1
    pos = {u: i for i, u in enumerate(order)}
    # reverse of MCS order is a perfect elimination ordering iff chordal
    for u in order:
        later = [w for w in adj[u] if pos[w] > pos[u]]
        if not later:
            continue
        p = min(later, key=lambda w: pos[w])   # earliest later neighbour
        for w in later:
            if w != p and w not in adj[p]:
                return False
    return True


def max_cliques(adj):
    n = len(adj); out = []
    def bk(R, P, X):
        if not P and not X:
            out.append(frozenset(R)); return
        pivot = max(P | X, key=lambda u: len(adj[u] & P))
        for u in list(P - adj[pivot]):
            bk(R | {u}, P & adj[u], X & adj[u]); P = P - {u}; X = X | {u}
    bk(set(), set(range(n)), set())
    return out


def analyse(v, build=True, check_sa=True):
    """returns a dict with the verdicts and (if built & verified) the tree model."""
    v = np.asarray(v, dtype=np.int64)
    res = {'sa_ssa': (sa_ssa_ok(v) if check_sa else True), 'zero_entries': int((v == 0).sum())}
    if not res['sa_ssa']:
        res.update(chordal=None, irreducible=False, verified=False); return res
    H = hyperedges(v); adj = line_graph(H)
    res['n_hyperedges'] = len(H)
    res['chordal'] = is_chordal(adj)
    cl = max_cliques(adj)
    qcap_empty = all(_intersect_all([H[i] for i in Q]) == 0 for Q in cl)
    res['irreducible'] = res['zero_entries'] == 0 and qcap_empty and len(H) > 0
    res['verified'] = False
    if build and res['chordal'] and res['irreducible']:
        model = algorithm1(v, H, adj, cl)
        if model is not None:
            res['verified'] = bool(np.array_equal(mincut_vector(model), v))
            res['model'] = model
    return res


def _intersect_all(ms):
    r = ALL
    for m in ms:
        r &= m
    return r


def algorithm1(v, H, adj, cliques):
    """Hubeny-Rota Algorithm 1: clique tree of L_P -> topological simple tree -> weights."""
    K = len(cliques)
    # weighted clique graph + maximum spanning tree (Kruskal)
    edges = []
    for i in range(K):
        for j in range(i + 1, K):
            w = len(cliques[i] & cliques[j])
            if w > 0:
                edges.append((w, i, j))
    edges.sort(reverse=True)
    parent = list(range(K))
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    tree = []
    for w, i, j in edges:
        ri, rj = find(i), find(j)
        if ri != rj:
            parent[ri] = rj; tree.append((i, j))
    if len(tree) != K - 1:
        return None                     # clique graph disconnected (should not happen when irreducible)
    # leaves for the parties
    T_edges = list(tree)                # bulk vertices 0..K-1
    leaf = {}
    for l in range(PARTIES):
        Ql = frozenset(i for i, X in enumerate(H) if X & (1 << l))
        cont = [q for q, Q in enumerate(cliques) if Ql <= Q]
        if len(cont) != 1:
            return None                 # uniqueness fails -> not irreducible / bug
        leaf[l] = K + l; T_edges.append((cont[0], K + l))
    nV = K + PARTIES
    nbr = [set() for _ in range(nV)]
    for a, b in T_edges:
        nbr[a].add(b); nbr[b].add(a)
    # weights: delete e -> side without purifier -> S_Y
    weights = []
    for a, b in T_edges:
        side = _component(nbr, a, b)
        Y = 0
        for l in range(PARTIES):
            if leaf[l] in side:
                Y |= 1 << l
        if Y & (1 << PUR):
            Y = ALL ^ Y
        weights.append(S_of(v, Y) if Y else 0)
    return {'edges': T_edges, 'weights': weights, 'leaf': leaf, 'n_bulk': K}


def _component(nbr, a, b):
    """vertices reachable from a when edge (a,b) is removed."""
    seen = {a}; st = [a]
    while st:
        u = st.pop()
        for w in nbr[u]:
            if (u == a and w == b) or w in seen:
                continue
            seen.add(w); st.append(w)
    return seen


def mincut_vector(model):
    """min-cut entropies of all 31 subsets of the 5 visible parties on the tree."""
    edges, weights, leaf = model['edges'], model['weights'], model['leaf']
    nV = model['n_bulk'] + PARTIES
    nbr = [dict() for _ in range(nV)]
    for (a, b), w in zip(edges, weights):
        nbr[a][b] = w; nbr[b][a] = w
    root = leaf[PUR]
    out = np.zeros(31, dtype=np.int64)
    INF = 10 ** 12
    for X in range(1, 32):
        side_of_leaf = {leaf[l]: (1 if X & (1 << l) else 0) for l in range(5)}
        side_of_leaf[leaf[PUR]] = 0
        def dp(u, p):
            """returns (cost if u on side 0, cost if u on side 1) for the subtree of u."""
            if u in side_of_leaf:
                s = side_of_leaf[u]
                base = [INF, INF]; base[s] = 0
            else:
                base = [0, 0]
            c0, c1 = base
            for w, wt in nbr[u].items():
                if w == p:
                    continue
                d0, d1 = dp(w, u)
                c0 += min(d0, d1 + wt)
                c1 += min(d1, d0 + wt)
            return c0, c1
        c0, c1 = dp(root, -1)
        out[X - 1] = c0
    return out
