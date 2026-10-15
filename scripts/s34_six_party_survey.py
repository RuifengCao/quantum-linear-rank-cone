#!/usr/bin/env python3
"""s34: six parties, the other direction -- states with at most one qudit per system (R40, R41; roadmap item M3), and
the closed cones (R42, Theorem 12).

Question (the paper, Discussion; R32 note section 5): is Stab^(p)_6 not contained in Stab^(2)_6 for some odd p, i.e. is
there a six-party stabilizer state of odd local dimension whose entropy vector lies outside the qubit cone?  The only
candidate was the non-Fano state h_NF = r_NF + r_NF* (section 4.8); the transfer argument that settles the other
direction for the Steane state does not work for it.  R42 (section 4.11, Theorem 12): yes, for every odd p.  The
inequality f12 >= 0 holds for every qubit stabilizer state of seven systems and f12(S_NF) = -1; its proof applies
Lemma 4 (LA)_2 to a subspace arrangement dual to the sum of the line pieces, whose dimensions are bounded by isotropy
(the "symplectic transfer").  transfer_form builds f12 from its definition, transfer_chain checks every step of the
proof on a given state, transfer_check runs the checks of gates G39 and G39b, transfer_anneal the adversarial search
of cc/s37_minform7.c.

Survey of the states with at most one qudit per system (seven systems: 0..5 visible, 6 the purifier):
  * cc/s35_enum7.c enumerates every weighted graph state over GF(p) on seven vertices, up to scaling of vertices by
    units.  Every stabilizer state over Z_p (p prime) of seven systems whose single-system entropies are at most 1
    (units of log p) is local-Clifford equivalent to such a graph state (a party-local Clifford leaves S(P) qudits per
    system), and entropies are local-unitary invariant; so the distinct cut-rank vectors are exactly the entropy
    vectors of these states.
  * Qutrit-only vectors (in the GF(3) list, not in the GF(2) one) fall into orbits of S_7 (relabelling the seven
    systems; the purifier is any of them).  cc/s34_anneal7.c searches qubit graph states realising 2 S.
  * Qubit-only vectors: the other set difference.

Conventions: a vector is stored on the 63 cuts (nonempty subsets of the visible systems, mask order, bit i = system i);
on all subsets of the seven systems it follows from purity, S(X) = S(complement).  For a seven-variable h (R36),
S(X) = h(X) - sum_{P in X} h(P) / 2; h_F = r_F + r_F* is the Steane state, h_NF = r_NF + r_NF* the non-Fano state
(scripts/s27_cube_states.py, scripts/s31_six_party.py).

Data: data/s34_seven_qudit_vectors.npz
  gf2, gf3          the two lists (uint8, rows of 63 cut ranks, sorted)
  orbit_vector      one representative of each qutrit-only orbit (uint8, 63 cut ranks)
  orbit_weights     a GF(3) graph realising it (int8, the 21 edges (0,1), (0,2), ..., (5,6))
  orbit_size        number of vectors of the orbit in the GF(3) list
  qubit2            a qubit graph state realising 2 S (uint8, 14 x 14, two qubits per system; zero if none)
  has_qubit2        bool
  nonfano           bool: the orbit of h_NF

  python3 scripts/s34_six_party_survey.py --check [--full]     # exact checks (gate G38; --full: G38b)
  python3 scripts/s34_six_party_survey.py --self-test          # conventions, cut ranks vs state vectors, line pieces
  python3 scripts/s34_six_party_survey.py --enumerate --out F  # cc/s35_enum7.c for p = 2 and 3 (p = 3: ~15 min, 2 cores)
  python3 scripts/s34_six_party_survey.py --survey --lists F --out data.npz --seconds 60   # orbits + qubit searches
  python3 scripts/s34_six_party_survey.py --nonfano --lams 2,3 --seconds 600               # h_NF at multiplicity 2, 3
  python3 scripts/s34_six_party_survey.py --transfer [--full --trials 400]                 # Theorem 12 (G39, G39b)
  python3 scripts/s34_six_party_survey.py --transfer-anneal "2,2,2,2,2,2,2;3,3,3,3,3,3,3" --seconds 120 [--pull-nf 1]
"""
import argparse, hashlib, itertools, json, os, random, shutil, subprocess, sys, tempfile, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, '..')
sys.path.insert(0, HERE); sys.path.insert(0, ROOT)
DATA = os.path.join(ROOT, 'data', 's34_seven_qudit_vectors.npz')
EDGES = [(u, v) for u in range(7) for v in range(u + 1, 7)]


# ----------------------------------------------------------------------------------------------- vectors
def entropy_from_h(h):
    h = np.asarray(h, dtype=np.int64)
    S = np.zeros(128, dtype=np.int64)
    for X in range(1, 128):
        two = 2 * h[X] - sum(int(h[1 << P]) for P in range(7) if X >> P & 1)
        assert two % 2 == 0
        S[X] = two // 2
    return S


def known_states():
    """{'F': S of the Steane state, 'NF': S of the non-Fano state} on the 128 subsets."""
    import s27_cube_states as s27
    import s31_six_party as s31
    rF, rNF = s27.dfz_configuration(2), s27.dfz_configuration(3)
    return {'F': entropy_from_h(rF + s31.h_star(rF)), 'NF': entropy_from_h(rNF + s31.h_star(rNF))}


def full(v):
    """S on the 128 subsets of the seven systems from the 63 visible cuts (purity)."""
    S = np.zeros(128, dtype=np.int64)
    for m in range(1, 64):
        S[m] = v[m - 1]
        S[127 ^ m] = v[m - 1]
    return S


_PERMS = None


def _perms():
    global _PERMS
    if _PERMS is None:
        _PERMS = np.zeros((5040, 128), dtype=np.int64)
        for a, p in enumerate(itertools.permutations(range(7))):
            for m in range(128):
                _PERMS[a, m] = sum(1 << p[i] for i in range(7) if m >> i & 1)
    return _PERMS


def orbit63(v):
    """all images of v under S_7, back on the 63 visible cuts (a set of tuples)."""
    S = full(v)
    P = _perms()
    imgs = np.zeros((5040, 128), dtype=np.int64)
    imgs[np.arange(5040)[:, None], P] = S[None, :]
    return {tuple(int(x) for x in row[1:64]) for row in imgs}


def canon(v):
    return min(orbit63(v))


# ----------------------------------------------------------------------------------------------- exact ranks
def rank_p(M, p):
    M = [[int(x) % p for x in row] for row in M]
    r, rows = 0, len(M)
    cols = len(M[0]) if rows else 0
    for c in range(cols):
        piv = next((i for i in range(r, rows) if M[i][c]), None)
        if piv is None:
            continue
        M[r], M[piv] = M[piv], M[r]
        inv = pow(M[r][c], p - 2, p)
        for i in range(rows):
            if i != r and M[i][c]:
                f = M[i][c] * inv % p
                M[i] = [(a - f * b) % p for a, b in zip(M[i], M[r])]
        r += 1
    return r


def cut_vector(A, p, sizes=None):
    """the 63 cut ranks over GF(p) of a graph with adjacency A (systems own consecutive vertices, sizes per system)."""
    sizes = sizes or [1] * 7
    owner = [q for q, c in enumerate(sizes) for _ in range(c)]
    N = len(owner)
    out = []
    for m in range(1, 64):
        X = [v for v in range(N) if owner[v] < 6 and m >> owner[v] & 1]
        Y = [v for v in range(N) if not (owner[v] < 6 and m >> owner[v] & 1)]
        out.append(rank_p([[A[a][b] for b in Y] for a in X], p))
    return tuple(out)


def graph_from_weights(w):
    A = [[0] * 7 for _ in range(7)]
    for (u, v), x in zip(EDGES, w):
        A[u][v] = A[v][u] = int(x)
    return A


def structure_ok(A, sizes):
    owner = [q for q, c in enumerate(sizes) for _ in range(c)]
    N = len(owner)
    return len(A) == N and all(A[i][j] == A[j][i] for i in range(N) for j in range(N)) and \
        all(A[i][j] == 0 for i in range(N) for j in range(N) if owner[i] == owner[j])


# ----------------------------------------------------------------------------------------------- other fields
def field(q):
    """(add, mul, neg, inv) tables of GF(q), q in 2, 3, 4, 5, 7, 8, 9, 11, 13, 16; elements 0..q-1.  q = 2^k: polynomial
    basis modulo x^2 + x + 1, x^3 + x + 1, x^4 + x + 1; q = 9: a + 3b stands for a + b i with i^2 = -1."""
    if q in (2, 3, 5, 7, 11, 13):
        add = [[(a + b) % q for b in range(q)] for a in range(q)]
        mul = [[a * b % q for b in range(q)] for a in range(q)]
    elif q in (4, 8, 16):
        k = q.bit_length() - 1
        poly = {2: 0b111, 3: 0b1011, 4: 0b10011}[k]

        def pm(a, b):
            r = 0
            while b:
                if b & 1:
                    r ^= a
                b >>= 1
                a <<= 1
                if a >> k & 1:
                    a ^= poly
            return r
        add = [[a ^ b for b in range(q)] for a in range(q)]
        mul = [[pm(a, b) for b in range(q)] for a in range(q)]
    elif q == 9:
        add = [[(a % 3 + b % 3) % 3 + 3 * ((a // 3 + b // 3) % 3) for b in range(9)] for a in range(9)]
        mul = [[(a % 3 * (b % 3) - a // 3 * (b // 3)) % 3 + 3 * ((a % 3 * (b // 3) + a // 3 * (b % 3)) % 3)
                for b in range(9)] for a in range(9)]
    else:
        raise ValueError(f'GF({q}) not supported')
    neg = [next(b for b in range(q) if add[a][b] == 0) for a in range(q)]
    inv = [0] + [next(b for b in range(q) if mul[a][b] == 1) for a in range(1, q)]
    assert all(mul[a][mul[b][c]] == mul[mul[a][b]][c] and mul[a][add[b][c]] == add[mul[a][b]][mul[a][c]]
               for a in range(q) for b in range(q) for c in range(q))
    return add, mul, neg, inv


def _reduce(M, ncols, F):
    """row-reduce the list of rows M in place over the field F; returns the pivot columns."""
    add, mul, neg, inv = F
    piv, r = [], 0
    for c in range(ncols):
        p = next((i for i in range(r, len(M)) if M[i][c]), None)
        if p is None:
            continue
        M[r], M[p] = M[p], M[r]
        iv = inv[M[r][c]]
        M[r] = [mul[iv][x] for x in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c]:
                f = neg[M[i][c]]
                M[i] = [add[x][mul[f][y]] for x, y in zip(M[i], M[r])]
        piv.append(c)
        r += 1
    return piv


def cut_vector_q(A, F):
    """the 63 cut ranks over the field F (tables from field) of a graph on seven vertices."""
    out = []
    for m in range(1, 64):
        X = [v for v in range(7) if v < 6 and m >> v & 1]
        Y = [v for v in range(7) if v not in X]
        out.append(len(_reduce([[A[a][b] for b in Y] for a in X], len(Y), F)))
    return tuple(out)


def target_search(t63, q, recheck='all'):
    """every weighted graph over GF(q) on seven vertices in scaling normal form whose cut-rank vector is t63: a second
    implementation of cc/s36_target7.c (same method, see there).  recheck: 'all' re-checks every solution from scratch
    on all 126 cuts, 'first' only the first one.  Returns dict(count, first, nodes)."""
    F = field(q)
    add, mul, neg, inv = F
    T = [0] * 128
    for X in range(1, 127):
        T[X] = t63[(127 ^ X) - 1] if X >> 6 & 1 else t63[X - 1]
    LO = [[0] * 128 for _ in range(7)]
    HI = [[0] * 128 for _ in range(7)]
    for k in range(7):                              # vertices 0..k placed
        placed = (1 << (k + 1)) - 1
        rest = 127 & ~placed
        subs = [S for S in range(128) if S & ~rest == 0]
        for X in range(placed + 1):
            vals = [T[X | S] for S in subs]
            HI[k][X], LO[k][X] = min(vals), max(vals) - (6 - k)
    A = [[0] * 7 for _ in range(7)]
    out = {'count': 0, 'first': None, 'nodes': [0] * 8}

    def dot(e, w):
        s = 0
        for x, y in zip(e, w):
            if x and y:
                s = add[s][mul[x][y]]
        return s

    def dfs(k):                                     # vertices 0..k-1 placed; vertex k joins
        out['nodes'][k] += 1
        if k == 7:
            if recheck == 'all' or out['count'] == 0:
                for X in range(1, 127):
                    rows = [v for v in range(7) if X >> v & 1]
                    cols = [v for v in range(7) if not X >> v & 1]
                    if len(_reduce([[A[i][j] for j in cols] for i in rows], len(cols), F)) != T[X]:
                        raise AssertionError('a solution failed the final check')
            out['count'] += 1
            if out['first'] is None:
                out['first'] = [A[u][v] for (u, v) in EDGES]
            return
        full_ = (1 << k) - 1
        lin, outs = [], []
        for Z in range(full_ + 1):
            Y = full_ & ~Z
            if not Y:
                continue
            rows = [v for v in range(k) if Z >> v & 1]
            cols = [v for v in range(k) if Y >> v & 1]
            M = [[A[i][j] for j in cols] for i in rows]
            piv = _reduce(M, len(cols), F)
            eqs = []                                # w_Y in the row space of A[Z, Y]  <=>  eqs . w = 0
            for j in range(len(cols)):
                if j not in piv:
                    e = [0] * k
                    e[cols[j]] = 1
                    for i, pc in enumerate(piv):
                        e[cols[pc]] = neg[M[i][j]]
                    eqs.append(e)
            r0 = len(piv)
            lo, hi = LO[k][Z | 1 << k], HI[k][Z | 1 << k]
            i_ok, o_ok = lo <= r0 <= hi, lo <= r0 + 1 <= hi
            if i_ok and o_ok:
                continue
            if not (i_ok or o_ok) or (o_ok and not eqs):
                return
            if i_ok:
                lin += eqs
            else:
                outs.append(eqs)
        piv = _reduce(lin, k, F)
        basis = []
        for c in range(k):
            if c not in piv:
                b = [0] * k
                b[c] = 1
                for i, pc in enumerate(piv):
                    b[pc] = neg[lin[i][c]]
                basis.append(b)
        for coeffs in itertools.product(range(q), repeat=len(basis)):
            w = [0] * k
            for c, b in zip(coeffs, basis):
                if c:
                    w = [add[x][mul[c][y]] for x, y in zip(w, b)]
            if next((x for x in w if x), 1) != 1:     # scaling normal form
                continue
            if any(all(dot(e, w) == 0 for e in g) for g in outs):
                continue
            for j in range(k):
                A[k][j] = A[j][k] = w[j]
            dfs(k + 1)
        for j in range(k):
            A[k][j] = A[j][k] = 0

    dfs(1)
    return out


def fields_check(qs=(2, 3, 4), recheck='all', log=print):
    """the Python search for S_F and S_NF over GF(q), q in qs: {q: [solutions for S_F, solutions for S_NF]}."""
    ks = known_states()
    res = {}
    for q in qs:
        res[q] = [target_search([int(x) for x in ks[n][1:64]], q, recheck)['count'] for n in ('F', 'NF')]
        log(f'  GF({q}): S_F {res[q][0]} solutions, S_NF {res[q][1]} (Python search)')
    return res


def enumerate_gf2_numpy():
    """the GF(2) list recomputed in numpy, independently of cc/s35_enum7.c: all 2^21 graphs on seven vertices; the
    rank of the (at most three) rows on the short side of a cut is log2 of the size of their span."""
    n = 1 << 21
    g = np.arange(n, dtype=np.uint32)
    eidx = {e: k for k, e in enumerate(EDGES)}
    R = np.zeros((n, 63), dtype=np.uint8)
    for m in range(1, 64):
        X = [v for v in range(6) if m >> v & 1]
        Y = [v for v in range(7) if v not in X]
        short, long_ = (X, Y) if len(X) <= len(Y) else (Y, X)
        span = [np.zeros(n, dtype=np.uint8)]
        for i in short:
            row = np.zeros(n, dtype=np.uint8)
            for pos, j in enumerate(long_):
                row |= (((g >> eidx[(min(i, j), max(i, j))]) & 1).astype(np.uint8) << pos)
            span = span + [s ^ row for s in span]
        Sp = np.sort(np.stack(span, axis=1), axis=1)
        distinct = 1 + (np.diff(Sp.astype(np.int16), axis=1) != 0).sum(axis=1)
        R[:, m - 1] = np.round(np.log2(distinct)).astype(np.uint8)
    return np.unique(R, axis=0)


# ----------------------------------------------------------------------------------------------- kernels
def build(name, tmp):
    exe = os.path.join(tmp, name)
    src = os.path.join(ROOT, 'cc', name + '.c')
    p = subprocess.run(['gcc', '-O3', '-march=native', '-fopenmp', '-o', exe, src, '-lm'], capture_output=True, text=True)
    if p.returncode != 0:
        p = subprocess.run(['gcc', '-O3', '-march=native', '-o', exe, src, '-lm'], capture_output=True, text=True)
        if p.returncode != 0:
            raise SystemExit('gcc failed:\n' + p.stderr)
    return exe


def enumerate_lists(out, log=print):
    """run cc/s35_enum7.c for p = 2 and p = 3; save the raw lists (vectors, counts, example GF(p) weights)."""
    tmp = tempfile.mkdtemp(prefix='s35_')
    try:
        exe = build('s35_enum7', tmp)
        res = {}
        for p in (2, 3):
            t0 = time.time()
            txt = subprocess.run([exe, str(p)], capture_output=True, text=True, check=True).stdout
            vecs, cnts, wts = [], [], []
            for line in txt.splitlines():
                t = line.split()
                if t[0] == 'V':
                    cnts.append(int(t[1])); vecs.append([int(x) for x in t[2:65]]); wts.append([int(x) for x in t[65:86]])
                elif t[0] == 'TOTAL':
                    total = int(t[1])
            order = np.lexsort(np.array(vecs).T[::-1])
            res[p] = dict(vec=np.array(vecs, dtype=np.uint8)[order], count=np.array(cnts, dtype=np.int64)[order],
                          w=np.array(wts, dtype=np.int8)[order], graphs=total)
            log(f'  GF({p}): {total} graphs in normal form, {len(vecs)} distinct vectors ({time.time() - t0:.0f} s)')
        np.savez_compressed(out, gf2=res[2]['vec'], gf2_count=res[2]['count'], gf2_w=res[2]['w'], gf2_graphs=res[2]['graphs'],
                            gf3=res[3]['vec'], gf3_count=res[3]['count'], gf3_w=res[3]['w'], gf3_graphs=res[3]['graphs'])
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def run_anneal(exe, jobs, seconds, seed, sweeps=100, pneu=0.2):
    """jobs: (id, lam, S on 128 subsets, start or None).  Returns {id: adjacency or None}."""
    lines = []
    for jid, lam, S, start in jobs:
        sizes = [lam * int(S[1 << p]) for p in range(7)]
        s = f'{jid} {lam} ' + ' '.join(map(str, sizes)) + ' ' + ' '.join(str(lam * int(S[m])) for m in range(1, 64))
        s += ' 0' if start is None else ' 1 ' + ' '.join(str(x) for row in start for x in row)
        lines.append(s)
    out = subprocess.run([exe, str(seconds), str(seed), '1.0', '0.05', str(sweeps), str(pneu), '0'],
                         input='\n'.join(lines) + '\n', capture_output=True, text=True, check=True).stdout
    res = {}
    for line in out.splitlines():
        t = line.split()
        if t[2] == 'FOUND':
            N = int(t[3]); rows = [int(x, 16) for x in t[4:4 + N]]
            res[int(t[0])] = [[(rows[i] >> j) & 1 for j in range(N)] for i in range(N)]
        else:
            res[int(t[0])] = None
    return res, out


# ----------------------------------------------------------------------------------------------- survey
def survey(lists, out, seconds, seed=7001, log=print):
    Z = np.load(lists)
    V2 = {tuple(int(x) for x in r) for r in Z['gf2']}
    V3 = {tuple(int(x) for x in r): k for k, r in enumerate(Z['gf3'])}
    cand = sorted(v for v in V3 if v not in V2)
    orbits, seen = {}, set()
    for v in cand:                                  # one orbit computation per orbit
        if v in seen:
            continue
        o = orbit63(v)
        seen |= o
        orbits[min(o)] = sorted(o & set(cand))
    log(f'  qutrit-only vectors {len(cand)} in {len(orbits)} orbits')
    nf = canon(tuple(int(x) for x in known_states()['NF'][1:64]))
    keys = sorted(orbits)
    tmp = tempfile.mkdtemp(prefix='s34_')
    try:
        exe = build('s34_anneal7', tmp)
        res, raw = run_anneal(exe, [(k, 2, full(orbits[key][0]), None) for k, key in enumerate(keys)], seconds, seed)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    K = len(keys)
    vec = np.zeros((K, 63), dtype=np.uint8); wts = np.zeros((K, 21), dtype=np.int8)
    q2 = np.zeros((K, 14, 14), dtype=np.uint8); has = np.zeros(K, dtype=bool); isnf = np.zeros(K, dtype=bool)
    size = np.zeros(K, dtype=np.int64)
    for k, key in enumerate(keys):
        v = orbits[key][0]
        vec[k] = v; wts[k] = Z['gf3_w'][V3[v]]; size[k] = len(orbits[key]); isnf[k] = key == nf
        A = res.get(k)
        S = full(v)
        sizes = [2 * int(S[1 << p]) for p in range(7)]
        if A is not None and len(A) <= 14 and structure_ok(A, sizes) and cut_vector(A, 2, sizes) == tuple(2 * x for x in v):
            q2[k, :len(A), :len(A)] = A; has[k] = True
        log(f'  orbit {k}: size {size[k]}, non-Fano {isnf[k]}, qubit realisation at multiplicity 2: {has[k]}')
    np.savez_compressed(out, gf2=Z['gf2'], gf3=Z['gf3'], orbit_vector=vec, orbit_weights=wts, orbit_size=size,
                        qubit2=q2, has_qubit2=has, nonfano=isnf)
    return raw


def field_runs(qs, controls=2, seed=36, log=print):
    """cc/s36_target7.c over GF(q), q in qs: the number of graphs in normal form with cut-rank vector S_F and S_NF, and
    `controls` random targets (cut-rank vectors of random graphs over GF(q), computed here) that must be found.  Every
    reported example is re-checked here.  Returns {q: dict(F=count, NF=count, controls_found=k)}."""
    ks = known_states()
    rng = random.Random(seed)
    tmp = tempfile.mkdtemp(prefix='s36_')
    res = {}
    try:
        exe = build('s36_target7', tmp)
        for q in qs:
            F = field(q)
            targets = [('F', tuple(int(x) for x in ks['F'][1:64])), ('NF', tuple(int(x) for x in ks['NF'][1:64]))]
            for c in range(controls):
                A = [[0] * 7 for _ in range(7)]
                for (u, v) in EDGES:
                    A[u][v] = A[v][u] = rng.randrange(q)
                targets.append((f'random{c}', cut_vector_q(A, F)))
            r = dict(controls_found=0)
            for name, t in targets:
                t0 = time.time()
                stop = ['1'] if name.startswith('random') else []          # controls: the first solution suffices
                out = subprocess.run([exe, str(q), '1'] + stop, input=' '.join(map(str, t)) + '\n',
                                     capture_output=True, text=True, check=True).stdout
                n = int(next(line for line in out.splitlines() if line.startswith('TOTAL')).split()[2])
                sol = [line.split()[1:] for line in out.splitlines() if line.startswith('SOL')]
                ok = all(cut_vector_q(graph_from_weights(s), F) == t for s in sol) and (n == 0) == (not sol)
                if not ok:
                    raise SystemExit(f'GF({q}) {name}: a reported solution failed the re-check')
                if name.startswith('random'):
                    r['controls_found'] += n > 0
                else:
                    r[name] = n
                log(f'  GF({q}) {name}: ' + ('found, example re-checked' if name.startswith('random') and n else
                                              f'{n} graphs in normal form' + (', example re-checked' if sol else ''))
                    + f' ({time.time() - t0:.1f} s)')
            res[q] = r
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    return res


def nonfano_runs(lams, seconds, seed=6001, log=print):
    S = known_states()['NF']
    tmp = tempfile.mkdtemp(prefix='s34_')
    try:
        exe = build('s34_anneal7', tmp)
        for lam in lams:
            res, raw = run_anneal(exe, [(lam, lam, S, None)], seconds, seed + lam)
            line = [l for l in raw.splitlines() if l.split()[0] == str(lam)][0]
            log(f'  h_NF at multiplicity {lam}: ' + ('FOUND' if res[lam] is not None else 'not found; ' + line))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


# ----------------------------------------------------------------------------------------------- checks
def load(path=DATA):
    return np.load(path, allow_pickle=False)


def content_sha(path=DATA):
    import hashlib
    Z = load(path)
    h = hashlib.sha256()
    for key in sorted(Z.files):
        a = Z[key]
        h.update(key.encode())
        h.update(np.ascontiguousarray(a, dtype=a.dtype.newbyteorder('<')).tobytes())
    return h.hexdigest()


def check(path=DATA, full_check=False, samples=20000, seed=38, log=print):
    """exact checks of the survey data.  Returns a summary dict."""
    Z = load(path)
    t0 = time.time()
    V2 = {tuple(int(x) for x in r) for r in Z['gf2']}
    V3 = {tuple(int(x) for x in r) for r in Z['gf3']}
    s = dict(gf2=len(V2), gf3=len(V3), gf2_rows=int(len(Z['gf2'])), gf3_rows=int(len(Z['gf3'])))
    K = len(Z['orbit_vector'])
    reps = [tuple(int(x) for x in r) for r in Z['orbit_vector']]
    # each representative is realised by its stored GF(3) graph, is in the GF(3) list and not in the GF(2) one
    s['reps_realised_gf3'] = sum(cut_vector(graph_from_weights(Z['orbit_weights'][k]), 3) == reps[k] for k in range(K))
    s['reps_qutrit_only'] = sum(v in V3 and v not in V2 for v in reps)
    # the orbits of the representatives partition the qutrit-only vectors exactly
    orbs = [orbit63(v) for v in reps]
    union = set().union(*orbs)
    s['orbits'] = K
    s['orbits_disjoint'] = sum(len(o) for o in orbs) == len(union)
    s['orbit_sizes_ok'] = all(len(orbs[k] & V3) == int(Z['orbit_size'][k]) and orbs[k] <= V3 for k in range(K))
    s['qutrit_only'] = len(V3 - V2)
    s['orbits_cover_qutrit_only'] = union == (V3 - V2)
    s['orbits_all_singles_one'] = sum(all(full(v)[1 << p] == 1 for p in range(7)) for v in reps)
    # qubit realisations at multiplicity 2
    ok2 = 0
    for k in range(K):
        if Z['has_qubit2'][k]:
            S = full(reps[k]); sizes = [2 * int(S[1 << p]) for p in range(7)]
            N = sum(sizes); A = Z['qubit2'][k][:N, :N].astype(int).tolist()
            ok2 += structure_ok(A, sizes) and cut_vector(A, 2, sizes) == tuple(2 * x for x in reps[k])
    s['qubit2_stored'] = int(Z['has_qubit2'].sum()); s['qubit2_ok'] = ok2
    # the orbit without one is the non-Fano orbit; the qubit-only vectors are the 30 relabellings of the Steane state
    ks = known_states()
    nf = tuple(int(x) for x in ks['NF'][1:64]); st = tuple(int(x) for x in ks['F'][1:64])
    s['without_qubit2'] = [k for k in range(K) if not Z['has_qubit2'][k]]
    s['nonfano_flag_ok'] = all(bool(Z['nonfano'][k]) == (nf in orbs[k]) for k in range(K))
    s['without_qubit2_is_nonfano'] = all(nf in orbs[k] for k in s['without_qubit2'])
    s['qubit_only'] = len(V2 - V3)
    s['qubit_only_is_steane'] = (V2 - V3) == orbit63(st)
    if full_check:
        g2 = enumerate_gf2_numpy()
        s['gf2_reenumerated'] = int(len(g2)); s['gf2_list_equal'] = np.array_equal(g2, Z['gf2'])
        rng = random.Random(seed)
        miss = 0
        for _ in range(samples):
            A = [[0] * 7 for _ in range(7)]
            for (u, v) in EDGES:
                A[u][v] = A[v][u] = rng.randrange(3)
            miss += cut_vector(A, 3) not in V3
        s['gf3_samples'] = samples; s['gf3_samples_missing'] = miss
        s.update(sums_report(Z))
    s['seconds'] = round(time.time() - t0, 1)
    return s


def sums_report(Z, log=print):
    """(1) for how many qutrit-only orbits is 2 v a sum of two vectors of the GF(2) list (a tensor product of two states
    with one qubit per system); (2) the L1 distance from S_NF to the GF(2) list, and from 2 S_NF to the sums of two
    GF(2) vectors (exhaustive over pairs)."""
    G2 = Z['gf2'].astype(np.int16)
    G2set = {r.tobytes() for r in Z['gf2']}
    sums = 0
    for v in Z['orbit_vector'].astype(np.int16):
        D = 2 * v[None, :] - G2
        sums += any(d.astype(np.uint8).tobytes() in G2set for d in D[(D >= 0).all(1)])
    ks = known_states()
    nf = np.array([int(x) for x in ks['NF'][1:64]], dtype=np.int16)
    st = np.array([int(x) for x in ks['F'][1:64]], dtype=np.int16)
    d1 = np.abs(G2 - nf[None, :]).sum(1)
    # 2 S_NF: is some sum a + b of GF(2) vectors at L1 distance 0 or 1, i.e. equal to 2 S_NF or 2 S_NF +- e_c?
    t = 2 * nf
    hits = {0: 0, 1: 0}
    for c in range(-1, 63):
        for sgn in ((0,) if c < 0 else (1, -1)):
            tt = t.copy()
            if c >= 0:
                tt[c] += sgn
            D = tt[None, :] - G2
            n = sum(d.astype(np.uint8).tobytes() in G2set for d in D[(D >= 0).all(1)])
            hits[0 if c < 0 else 1] += n
    d_ff = int(np.abs(2 * st - t).sum())            # S_F + S_F
    best2 = 0 if hits[0] else 1 if hits[1] else d_ff if d_ff == 2 else -1
    return dict(orbits_2v_a_sum=int(sums), nf_l1_to_gf2=int(d1.min()), nf_l1_nearest=int((d1 == d1.min()).sum()),
                nf_nearest_is_steane=bool(np.array_equal(G2[int(d1.argmin())], st)), two_nf_l1_to_sums=best2)


def check_ok(s, full_check=False):
    ok = (s['reps_realised_gf3'] == s['orbits'] == s['reps_qutrit_only'] and s['orbits_disjoint'] and s['orbit_sizes_ok']
          and s['orbits_cover_qutrit_only'] and s['qubit2_ok'] == s['qubit2_stored'] and s['nonfano_flag_ok']
          and s['without_qubit2_is_nonfano'] and s['qubit_only_is_steane'] and s['gf2'] == s['gf2_rows']
          and s['gf3'] == s['gf3_rows'])
    if full_check:
        ok = (ok and s['gf2_list_equal'] and s['gf3_samples_missing'] == 0 and s['orbits_2v_a_sum'] == 0
              and s['nf_l1_to_gf2'] == 1 and s['nf_l1_nearest'] == 1 and s['nf_nearest_is_steane']
              and s['two_nf_l1_to_sums'] == 2)
    return ok


LINES_NF = [(0, 1, 3), (0, 2, 4), (1, 2, 5), (0, 5, 6), (1, 4, 6), (2, 3, 6)]   # ABW, ACX, BCY, AYZ, BXZ, CWZ
LINE_WXY = (3, 4, 5)                                                         # the seventh Fano line


def _nullspace_left(M, p):
    """basis (rows) of {c : c M = 0 mod p} for an integer matrix M (k x n)."""
    k = M.shape[0]
    A = np.concatenate([M % p, np.eye(k, dtype=np.int64)], axis=1)
    n = M.shape[1]
    r = 0
    for c in range(n):
        piv = next((i for i in range(r, k) if A[i, c] % p), None)
        if piv is None:
            continue
        A[[r, piv]] = A[[piv, r]]
        A[r] = A[r] * pow(int(A[r, c]), p - 2, p) % p
        for i in range(k):
            if i != r and A[i, c] % p:
                A[i] = (A[i] - A[i, c] * A[r]) % p
        r += 1
    return A[r:, n:] % p


def _rank_mod(M, p):
    return rank_p(M.tolist(), p) if len(M) else 0


def line_pieces_check(seed=10, transvections=30):
    """Proposition 10, steps 1 and 1', on states that exist: the CSS states of the non-Fano code over GF(3) and of the
    Fano (simplex) code over GF(2), with m = 1 and m = 2 copies, after a random symplectic map on each system (products of
    random transvections).  For every line l the elements L_l of the stabilizer supported on l must form an m-dimensional
    space, and at every point P the restrictions of the L_l (l through P) must be one Lagrangian subspace of GF(p)^{2m}.
    Also checked: the elements supported on WXY (m-dimensional for Fano, 0 for non-Fano).  Returns a dict of flags."""
    rng = np.random.default_rng(seed)
    V = np.array([[1, 0, 0, 1, 1, 0, 1], [0, 1, 0, 1, 0, 1, 1], [0, 0, 1, 0, 1, 1, 1]], dtype=np.int64)  # A B C W X Y Z
    out = {}
    for p, tag, lines in ((3, 'NF', LINES_NF), (2, 'F', LINES_NF + [LINE_WXY])):
        Cperp = _nullspace_left(V.T % p, p)                    # the dual code: d with V d = 0
        for m in (1, 2):
            N = 7 * m
            rows = []
            for j in range(m):                                  # copy j uses qudit P*m + j of system P
                for c in V:
                    v = np.zeros(2 * N, dtype=np.int64); v[[P * m + j for P in range(7)]] = c; rows.append(v)
                for d in Cperp:
                    v = np.zeros(2 * N, dtype=np.int64); v[[N + P * m + j for P in range(7)]] = d; rows.append(v)
            G = np.array(rows) % p
            cols = {P: [P * m + j for j in range(m)] + [N + P * m + j for j in range(m)] for P in range(7)}
            for P in range(7):                                  # random symplectic map on system P
                for _ in range(transvections):
                    w = rng.integers(0, p, size=2 * m)
                    c = int(rng.integers(1, p))
                    loc = G[:, cols[P]]
                    om = (loc[:, :m] @ w[m:] - loc[:, m:] @ w[:m]) % p       # omega_P(row, w)
                    G[:, cols[P]] = (loc + c * om[:, None] * w[None, :]) % p

            def omega(a, b):
                return int((a[:N] @ b[N:] - a[N:] @ b[:N]) % p)
            iso = all(omega(G[i], G[j]) == 0 for i in range(len(G)) for j in range(len(G)))

            def supported(S):
                outside = [c for P in range(7) if P not in S for c in cols[P]]
                K = _nullspace_left(G[:, outside], p)
                return (K @ G) % p if len(K) else np.zeros((0, 2 * N), dtype=np.int64)
            pieces = {l: supported(l) for l in lines}
            dims_ok = all(_rank_mod(B, p) == m for B in pieces.values())
            restr_ok = True
            for P in range(7):
                R = [pieces[l][:, cols[P]] for l in lines if P in l]
                r0 = R[0]
                for Rk in R:
                    if not (_rank_mod(Rk, p) == m == _rank_mod(np.concatenate([r0, Rk]), p)):
                        restr_ok = False
                lag = all(int((a[:m] @ b[m:] - a[m:] @ b[:m]) % p) == 0 for a in r0 for b in r0)
                restr_ok = restr_ok and lag
            wxy = _rank_mod(supported(LINE_WXY), p)
            zplane = np.concatenate([np.zeros((m, m), dtype=np.int64), np.eye(m, dtype=np.int64)], axis=1)
            moved = sum(_rank_mod(np.concatenate([pieces[lines[0]][:, cols[P]], zplane]), p) > m for P in lines[0])
            out[f'{tag}_m{m}'] = dict(isotropic=iso, line_dims=dims_ok, restrictions_equal_lagrangian=restr_ok,
                                      wxy_dim=wxy, wxy_ok=(wxy == (m if tag == 'F' else 0)), moved_points=int(moved))
    out['ok'] = all(v['isotropic'] and v['line_dims'] and v['restrictions_equal_lagrangian'] and v['wxy_ok']
                    for k, v in out.items() if k != 'ok')
    return out


# ----------------------------------------------------------------------------------------------- the closed cones (R42)
TR_LINES = {'ABW': (0, 1, 3), 'ACX': (0, 2, 4), 'BCY': (1, 2, 5), 'AYZ': (0, 5, 6), 'BXZ': (1, 4, 6), 'CWZ': (2, 3, 6)}
TR_PAIRS = [(0, 1), (0, 2), (1, 2), (0, 6), (1, 6), (2, 6)]           # AB, AC, BC, AZ, BZ, CZ: the corner pairs
TR_NAMES = 'ABCWXYZ'


def _mask(ps):
    return sum(1 << p for p in ps)


def transfer_form():
    """The inequality f12(S) >= 0 of Theorem 12 (R32 note, section 4.11) as 128 integer coefficients, index = subset
    mask over A, B, C, W, X, Y, Z = systems 0..6 (f12 = F . S on all 128 subsets).  Built from its definition:
        k(T) = sum_{P in T} S(P) - S(T),   d(P, l) = k(l) - k(l - P)  for a line l through P,
        e_P(l1, l2) = 2 S(P) - d(P, l1) - d(P, l2)  at a corner P, for two of the lines through P,
        U_Z = 4 S(Z) - d(Z, AYZ) - d(Z, BXZ) - d(Z, CWZ),
        U_W = 2 S(W) - d(W, ABW) + e_A(ACX, AYZ) + e_B(BCY, BXZ),
        U_X = 2 S(X) - d(X, ACX) + e_A(ABW, AYZ) + e_C(BCY, CWZ),
        U_Y = 2 S(Y) - d(Y, BCY) + e_B(ABW, BXZ) + e_C(ACX, CWZ),
        f12 = 5 U_Z + 4 U_W + U_X + U_Y + k(WXY) + 6 k(ABC) + sum_{corner pairs} k - sum_{lines} k
              - 3 k(ABW) - 3 k(CWZ).
    The corner pairs are AB, AC, BC, AZ, BZ, CZ, the lines ABW, ACX, BCY, AYZ, BXZ, CWZ."""
    F = np.zeros(128, dtype=np.int64)
    L = {n: _mask(t) for n, t in TR_LINES.items()}

    def S1(P, c):
        F[1 << P] += c

    def k(T, c):
        for P in range(7):
            if T >> P & 1:
                F[1 << P] += c
        F[T] -= c

    def d(P, l, c):
        k(L[l], c)
        k(L[l] & ~(1 << P), -c)

    def e(P, l1, l2, c):
        S1(P, 2 * c)
        d(P, l1, -c)
        d(P, l2, -c)
    A_, B_, C_, W_, X_, Y_, Z_ = range(7)
    S1(Z_, 4 * 5); d(Z_, 'AYZ', -5); d(Z_, 'BXZ', -5); d(Z_, 'CWZ', -5)                      # 5 U_Z
    S1(W_, 2 * 4); d(W_, 'ABW', -4); e(A_, 'ACX', 'AYZ', 4); e(B_, 'BCY', 'BXZ', 4)          # 4 U_W
    S1(X_, 2); d(X_, 'ACX', -1); e(A_, 'ABW', 'AYZ', 1); e(C_, 'BCY', 'CWZ', 1)             # U_X
    S1(Y_, 2); d(Y_, 'BCY', -1); e(B_, 'ABW', 'BXZ', 1); e(C_, 'ACX', 'CWZ', 1)             # U_Y
    k(_mask((W_, X_, Y_)), 1)
    k(_mask((A_, B_, C_)), 6)
    for pr in TR_PAIRS:
        k(_mask(pr), 1)
    for l in L.values():
        k(l, -1)
    k(L['ABW'], -3); k(L['CWZ'], -3)
    return F


def transfer_terms(F=None):
    """f12 folded onto the 63 cuts of the visible parties (Z = the purifier): [(name, coefficient)], sorted."""
    F = transfer_form() if F is None else F
    G = np.zeros(128, dtype=np.int64)
    for i, c in enumerate(F):
        G[i if not i >> 6 & 1 else 127 ^ i] += c
    out = [(''.join(TR_NAMES[p] for p in range(7) if i >> p & 1), int(c)) for i, c in enumerate(G) if c]
    return sorted(out, key=lambda t: (len(t[0]), t[0]))


def _gf2_rank(M):
    return rank_p(M.tolist(), 2) if len(M) and M.shape[1] else 0


class _Lagrangian:
    """a qubit stabilizer state: generator matrix G (N x 2N over GF(2), columns x_0..x_{N-1}, z_0..z_{N-1}) and the
    system of every qubit."""

    def __init__(self, G, owner):
        self.G = np.asarray(G, dtype=np.int64) % 2
        self.owner = list(owner)
        self.N = len(owner)
        self.cols = {P: [v for v in range(self.N) if owner[v] == P] + [self.N + v for v in range(self.N) if owner[v] == P]
                     for P in range(7)}

    def supported(self, T, basis=None):
        """basis rows of the elements of span(basis) (default: L) supported in the systems of T."""
        Gm = self.G if basis is None else basis
        out = [c for P in range(7) if not T >> P & 1 for c in self.cols[P]]
        if not len(Gm):
            return Gm
        if not out:
            return Gm
        K = _nullspace_left(Gm[:, out], 2)
        return (K @ Gm) % 2 if len(K) else np.zeros((0, Gm.shape[1]), dtype=np.int64)

    def proj(self, M, P):
        return M[:, self.cols[P]] if len(M) else np.zeros((0, len(self.cols[P])), dtype=np.int64)

    def entropies(self):
        S = np.zeros(128, dtype=np.int64)
        for T in range(1, 127):
            nT = sum(1 for v in range(self.N) if T >> self.owner[v] & 1)
            S[T] = nT - _gf2_rank(self.supported(T))
        return S


def _form_rank(M, n):
    """rank of the symplectic form restricted to the row space of M (n qubits, x-part first): rank of M J M^T."""
    return _gf2_rank(_omega(M, M, n)) if len(M) else 0


def transfer_chain(st, S=None):
    """Every step of the proof of Theorem 12 on one compressed qubit stabilizer state st (S(P) = number of qubits of P):
    dim L_T = k(T) for all T; the dimensions d(P, l) of the restrictions of the line pieces; their orthogonality at the
    common point of two lines; Lemma T1 (both bounds) at the four corners; Lemma T2 on the six lines; Lemma T3
    (dim D_P <= U_P at Z, W, X, Y); Lemma T4 (i)-(iii); (LA)_2 on the rank function of Lemma T5; and f12(S) >= 0.
    Returns (values, list of failed steps)."""
    S = st.entropies() if S is None else S
    s = [int(S[1 << P]) for P in range(7)]
    L = {n: _mask(t) for n, t in TR_LINES.items()}
    k = lambda T: sum(s[P] for P in range(7) if T >> P & 1) - int(S[T])
    d = lambda P, l: k(L[l]) - k(L[l] & ~(1 << P))
    bad = []
    for T in range(1, 128):
        if _gf2_rank(st.supported(T)) != k(T):
            bad.append(f'dim L_T != k(T) at {T}')
            break
    pieces = {l: st.supported(T) for l, T in L.items()}
    for l, T in L.items():
        for P in TR_LINES[l]:
            if _gf2_rank(st.proj(pieces[l], P)) != d(P, l):
                bad.append(f'dim of the restriction of {l} at {TR_NAMES[P]}')
    through = {P: [l for l, t in TR_LINES.items() if P in t] for P in range(7)}
    for P in range(7):
        n = s[P]
        for i, l1 in enumerate(through[P]):
            for l2 in through[P][i + 1:]:
                V1, V2 = st.proj(pieces[l1], P), st.proj(pieces[l2], P)
                if len(V1) and len(V2) and _omega(V1, V2, n).any():
                    bad.append(f'restrictions of {l1} and {l2} not orthogonal at {TR_NAMES[P]}')
    for P in (0, 1, 2, 6):                                          # Lemma T1 at the corners
        n = s[P]
        ls = through[P]
        for l in ls:
            o = [m for m in ls if m != l]
            if _form_rank(st.proj(pieces[l], P), n) > 2 * n - d(P, o[0]) - d(P, o[1]):
                bad.append(f'Lemma T1 (rank) at {TR_NAMES[P]} for {l}')
        Vs = [st.proj(pieces[l], P) for l in ls if len(pieces[l])]
        if Vs and _gf2_rank(np.concatenate(Vs, axis=0)) > 4 * n - sum(d(P, l) for l in ls):
            bad.append(f'Lemma T1 (sum) at {TR_NAMES[P]}')
    for l, t in TR_LINES.items():                                   # Lemma T2: from the two corners to the diagonal point
        P1, P2 = [P for P in t if P in (0, 1, 2, 6)]
        Qd = next(P for P in t if P in (3, 4, 5))
        if len(pieces[l]) and (_form_rank(st.proj(pieces[l], Qd), s[Qd]) > _form_rank(st.proj(pieces[l], P1), s[P1])
                               + _form_rank(st.proj(pieces[l], P2), s[P2])):
            bad.append(f'Lemma T2 on {l}')
    nz = [x for x in pieces.values() if len(x)]
    Z0 = np.concatenate(nz, axis=0) if nz else np.zeros((0, 2 * st.N), dtype=np.int64)
    D = [_gf2_rank(st.proj(Z0, P)) for P in range(7)]
    z = lambda T: _gf2_rank(st.supported(T, basis=Z0)) if len(Z0) else 0
    A_, B_, C_, W_, X_, Y_, Z_ = range(7)
    e = lambda P, l1, l2: 2 * s[P] - d(P, l1) - d(P, l2)
    U = {Z_: 4 * s[Z_] - d(Z_, 'AYZ') - d(Z_, 'BXZ') - d(Z_, 'CWZ'),
         W_: 2 * s[W_] - d(W_, 'ABW') + e(A_, 'ACX', 'AYZ') + e(B_, 'BCY', 'BXZ'),
         X_: 2 * s[X_] - d(X_, 'ACX') + e(A_, 'ABW', 'AYZ') + e(C_, 'BCY', 'CWZ'),
         Y_: 2 * s[Y_] - d(Y_, 'BCY') + e(B_, 'ABW', 'BXZ') + e(C_, 'ACX', 'CWZ')}
    for P, u in U.items():
        if D[P] > u:
            bad.append(f'Lemma T3: dim D_{TR_NAMES[P]} = {D[P]} > {u}')
    zs = {T: z(T) for T in range(1, 128)}
    if any(zs[T] > k(T) for T in zs):
        bad.append('Lemma T4(i): z(T) > k(T)')
    if any(zs[T] != k(T) for T in L.values()):
        bad.append('Lemma T4(ii): z(l) != k(l)')
    z = lambda T: zs[T] if T else 0
    Q = _mask((A_, B_, C_, Z_))
    lz = k(L['ABW']) + k(L['CWZ']) - D[W_]                         # Lemma T4(iii)
    if z(Q) < lz:
        bad.append(f'Lemma T4(iii): z(ABCZ) = {z(Q)} < {lz}')
    r = lambda T: sum(D[P] for P in range(7) if T >> P & 1) - z(T)
    pair = lambda a, b: _mask((a, b))
    lhs = r(_mask((W_, X_, Y_)))
    rhs = (2 * r(1 << Z_) + 3 * (r(1 << A_) + r(1 << B_) + r(1 << C_) - r(_mask((A_, B_, C_))))
           + 3 * (r(Q) - r(_mask((A_, B_, C_))))
           + r(L['BCY']) - r(pair(B_, C_)) + r(L['ACX']) - r(pair(A_, C_)) + r(L['ABW']) - r(pair(A_, B_))
           + r(L['AYZ']) - r(pair(A_, Z_)) + r(L['BXZ']) - r(pair(B_, Z_)) + r(L['CWZ']) - r(pair(C_, Z_)))
    if lhs > rhs:
        bad.append(f'(LA)_2 on the dual arrangement: {lhs} > {rhs}')
    f = int(transfer_form() @ S)
    if f < 0:
        bad.append(f'f12 = {f}')
    return dict(f12=f, D=D, U=U, zQ=(z(Q), lz), LA=(lhs, rhs)), bad


def _dual_identity(st):
    """Lemma T5 on one state: build the dual arrangement of Z0 explicitly (the column groups of a generator matrix of
    the dot-product complement of Z0 inside the sum of the D_P, each D_P identified with GF(2)^dim through a basis) and
    compare its ranks with dim D_T - z(T) on all 127 sets T."""
    nz = [x for x in (st.supported(_mask(t)) for t in TR_LINES.values()) if len(x)]
    if not nz:
        return True                                                # Z0 = 0: every D_P = 0 and every rank is 0
    Z0 = np.concatenate(nz, axis=0)
    blocks, sizes = [], []
    for P in range(7):
        Pz = st.proj(Z0, P)
        Bs = []                                                    # a basis of D_P = rowspace(Pz), greedily
        for row in Pz:
            if _gf2_rank(np.array(Bs + [row])) > len(Bs):
                Bs.append(row)
        Bm = np.array(Bs, dtype=np.int64).reshape(len(Bs), Pz.shape[1])
        sizes.append(len(Bs))
        if len(Bs):
            coef = np.zeros((len(Pz), len(Bs)), dtype=np.int64)    # Z0 in coordinates of the basis of D_P
            for i, row in enumerate(Pz):                            # solve x Bm = row
                Aug = np.concatenate([Bm.T, row.reshape(-1, 1)], axis=1) % 2
                sol = _solve_gf2(Aug)
                coef[i] = sol
            blocks.append(coef)
        else:
            blocks.append(np.zeros((len(Z0), 0), dtype=np.int64))
    Zc = np.concatenate(blocks, axis=1) % 2
    Gd = _nullspace_left(Zc.T, 2)                                  # rows y with Zc y = 0
    offs = np.cumsum([0] + sizes)
    for T in range(1, 128):
        cols = [c for P in range(7) if T >> P & 1 for c in range(offs[P], offs[P + 1])]
        rd = _gf2_rank(Gd[:, cols]) if cols and len(Gd) else 0
        if rd != sum(sizes[P] for P in range(7) if T >> P & 1) - _gf2_rank(st.supported(T, basis=Z0)):
            return False
    return True


def _solve_gf2(Aug):
    """one solution x of M x = b over GF(2), Aug = [M | b] (consistent systems only)."""
    M = Aug.copy() % 2
    rows, cols = M.shape
    n = cols - 1
    piv, r = [], 0
    for c in range(n):
        p_ = next((i for i in range(r, rows) if M[i, c]), None)
        if p_ is None:
            continue
        M[[r, p_]] = M[[p_, r]]
        for i in range(rows):
            if i != r and M[i, c]:
                M[i] = (M[i] + M[r]) % 2
        piv.append(c)
        r += 1
    assert not any(M[i, n] for i in range(r, rows)), 'inconsistent'
    x = np.zeros(n, dtype=np.int64)
    for i, c in enumerate(piv):
        x[c] = M[i, n]
    return x


def _omega(U, V, n):
    """the standard symplectic form between the rows of U and of V (x-part first, n qubits): U_x V_z^T + U_z V_x^T."""
    return (U[:, :n] @ V[:, n:].T + U[:, n:] @ V[:, :n].T) % 2


def compress(st):
    """symplectic reduction on every system: the elements of L supported on P form an isotropic I_P, and pi_P(L) =
    I_P^perp; L / (sum of the I_P) is a Lagrangian in the sum of the I_P^perp / I_P, with S(P) qubits on P and the same
    entropy vector (section 4.10: the party-local Clifford that factors out n_P - S(P) qubits)."""
    G = st.G % 2
    newcols, owner = [], []
    for P in range(7):
        cols = st.cols[P]
        n = len(cols) // 2
        if n == 0:
            continue
        loc = G[:, cols]
        I = st.proj(st.supported(1 << P), P)
        basis = []
        for row in I:                                               # a basis of I_P
            if _gf2_rank(np.array(basis + [row])) > len(basis):
                basis.append(row)
        comp = []
        for row in loc:                                             # extend to a basis of pi_P(L) = I_P^perp
            if _gf2_rank(np.array(basis + comp + [row])) > len(basis) + len(comp):
                comp.append(row)
        es, fs = [], []
        comp = [np.array(c) % 2 for c in comp]
        while comp:                                                 # symplectic Gram-Schmidt on the complement
            e = comp.pop(0)
            j = next((j for j, c in enumerate(comp) if _omega(e[None], c[None], n)[0, 0]), None)
            assert j is not None, 'degenerate complement'
            f = comp.pop(j)
            comp = [(c + _omega(c[None], f[None], n)[0, 0] * e + _omega(c[None], e[None], n)[0, 0] * f) % 2 for c in comp]
            es.append(e); fs.append(f)
        k = len(es)
        if k:
            E_, F_ = np.array(es), np.array(fs)
            newcols.append(np.concatenate([_omega(loc, F_, n), _omega(loc, E_, n)], axis=1))   # x_j = w(v, f_j), z_j = w(v, e_j)
            owner += [P] * k
    N = len(owner)
    if N == 0:
        return _Lagrangian(np.zeros((0, 0), dtype=np.int64), [])
    Xs, Zs = [], []
    for B in newcols:
        k = B.shape[1] // 2
        Xs.append(B[:, :k]); Zs.append(B[:, k:])
    Gn = np.concatenate(Xs + Zs, axis=1) % 2
    R = []
    for row in Gn:                                                  # a basis of the image
        if _gf2_rank(np.array(R + [row])) > len(R):
            R.append(row)
    return _Lagrangian(np.array(R, dtype=np.int64), owner)


def _css_state(M, owner):
    """the CSS state of the code spanned by the rows of M (columns distributed by owner): L = C x C^perp."""
    M = np.asarray(M, dtype=np.int64) % 2
    N = M.shape[1]
    R = []
    for row in M:
        if _gf2_rank(np.array(R + [row])) > len(R):
            R.append(row)
    C = np.array(R, dtype=np.int64).reshape(len(R), N)
    Cp = _nullspace_left(C.T, 2) if len(C) else np.eye(N, dtype=np.int64)
    G = np.concatenate([np.concatenate([C, np.zeros_like(C)], axis=1),
                        np.concatenate([np.zeros((len(Cp), N), dtype=np.int64), Cp], axis=1)], axis=0) % 2
    order = sorted(range(N), key=lambda v: owner[v])                # qubits grouped by system
    G = G[:, order + [N + v for v in order]]
    return _Lagrangian(G, [owner[v] for v in order])


def _random_structured_css(rng):
    """CSS state of a random GF(2) family pushed towards the non-Fano dependencies: A, B, C random, Z mostly inside
    A + B + C, and W, X, Y mostly inside (A + B) & (C + Z), (A + C) & (B + Z), (B + C) & (A + Z); then compressed.
    Over GF(2) these intersections force the Fano dependency on W, X, Y; the near-misses are the point."""
    while True:
        r = int(rng.integers(3, 7))
        span = lambda Ms: np.concatenate(Ms, axis=0) if Ms else np.zeros((0, r), dtype=np.int64)
        rnd = lambda k: rng.integers(0, 2, size=(k, r))
        U = {}
        for P in (0, 1, 2):
            U[P] = rnd(int(rng.integers(1, 3)))
        ABC = span([U[0], U[1], U[2]])
        U[6] = (rng.integers(0, 2, size=(int(rng.integers(1, 3)), len(ABC))) @ ABC) % 2 if rng.random() < 0.8 else rnd(1)

        def inter(P1, P2):
            Pm, Qm = P1, P2
            K = _nullspace_left(np.concatenate([Pm, Qm], axis=0) % 2, 2)    # c with c [P; Q] = 0
            return (K[:, :len(Pm)] @ Pm) % 2 if len(K) else np.zeros((0, r), dtype=np.int64)
        for D, (a, b, c, d) in ((3, (0, 1, 2, 6)), (4, (0, 2, 1, 6)), (5, (1, 2, 0, 6))):
            I = inter(span([U[a], U[b]]), span([U[c], U[d]]))
            k = int(rng.integers(1, 3))
            if len(I) and rng.random() < 0.85:
                V = (rng.integers(0, 2, size=(k, len(I))) @ I) % 2
            else:
                V = (rng.integers(0, 2, size=(k, len(span([U[a], U[b]])))) @ span([U[a], U[b]])) % 2
            if rng.random() < 0.15:
                V = np.concatenate([V, rnd(1)], axis=0)
            U[D] = V
        cols, owner = [], []
        for P in range(7):
            for v in U[P]:
                cols.append(v); owner.append(P)
        M = np.array(cols, dtype=np.int64).T                        # r x N, columns = the vectors
        st = compress(_css_state(M, owner))
        if st.N and len(set(st.owner)) == 7:
            return st


def _random_compressed_graph_state(rng):
    while True:
        sizes = rng.integers(1, 4, size=7)
        owner = [p for p in range(7) for _ in range(sizes[p])]
        N = len(owner)
        U = (rng.random((N, N)) < rng.uniform(0.15, 0.85)).astype(np.int64)
        Adj = np.triu(U, 1)
        Adj = Adj + Adj.T
        for i in range(N):
            for j in range(N):
                if owner[i] == owner[j]:
                    Adj[i, j] = 0
        st = _Lagrangian(np.concatenate([np.eye(N, dtype=np.int64), Adj], axis=1), owner)
        S = st.entropies()
        if all(S[1 << P] == sizes[P] for P in range(7)):
            return st, S


def _fano_copies(rng, m, transvections=20):
    """m copies of the Fano (simplex) CSS state, after random symplectic maps on every system: the tight case."""
    V = np.array([[1, 0, 0, 1, 1, 0, 1], [0, 1, 0, 1, 0, 1, 1], [0, 0, 1, 0, 1, 1, 1]], dtype=np.int64)
    Cp = _nullspace_left(V.T, 2)
    N = 7 * m
    rows = []
    for j in range(m):
        for c in V:
            v = np.zeros(2 * N, dtype=np.int64); v[[P * m + j for P in range(7)]] = c; rows.append(v)
        for dd in Cp:
            v = np.zeros(2 * N, dtype=np.int64); v[[N + P * m + j for P in range(7)]] = dd; rows.append(v)
    G = np.array(rows) % 2
    owner = [P for P in range(7) for _ in range(m)]
    st = _Lagrangian(G, owner)
    for P in range(7):
        cols = st.cols[P]
        for _ in range(transvections):
            w = rng.integers(0, 2, size=2 * m)
            loc = G[:, cols]
            om = (loc[:, :m] @ w[m:] + loc[:, m:] @ w[:m]) % 2
            G[:, cols] = (loc + om[:, None] * w[None, :]) % 2
    return _Lagrangian(G, owner)


def transfer_folded(F=None):
    """f12 on the 63 masks of the visible parties A, B, C, W, X, Y (index = mask, Z = the purifier folded in by purity)."""
    F = transfer_form() if F is None else F
    G = np.zeros(64, dtype=np.int64)
    for i, c in enumerate(F):
        G[i if not i >> 6 & 1 else 127 ^ i] += c
    return G


def transfer_sha(F=None):
    """sha256 of the 128 coefficients of f12 (int64, little-endian, index = subset mask): the manifest pin."""
    F = transfer_form() if F is None else F
    return hashlib.sha256(np.asarray(F, dtype='<i8').tobytes()).hexdigest()


def transfer_anneal(configs, seconds=60.0, seed=4201, mu=0.0, target=None, T0=3.0, log=print):
    """Adversarial search (empirical): cc/s37_minform7.c minimises F over qubit graph states with the given numbers of
    qubits per system (one process per configuration, run in parallel), optionally pulled towards a target vector
    (63 integers) with weight mu.  Every reported graph is rebuilt here, its entropies recomputed exactly, compressed,
    and passed through transfer_chain (all steps of the proof).  Returns one dict per configuration."""
    tmp = tempfile.mkdtemp(prefix='s37_')
    try:
        exe = build('s37_minform7', tmp)
        G = transfer_folded()
        inp = ' '.join(str(int(G[m])) for m in range(1, 64))
        if target is not None:
            inp += f' {mu} ' + ' '.join(str(int(target[m - 1])) for m in range(1, 64))
        procs = []
        for k, sz in enumerate(configs):
            args = [exe] + [str(x) for x in sz] + [str(seconds), str(seed + 17 * k), str(T0)]
            pr = subprocess.Popen(args, stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True)
            pr.stdin.write(inp + '\n')
            pr.stdin.close()                                        # all processes start searching at once
            procs.append(pr)
        outs = []
        for pr in procs:
            outs.append(pr.stdout.read())
            pr.wait()
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    res = []
    for sz, out in zip(configs, outs):
        t = out.split()
        minF = int(t[1]); obj = int(t[3]); restarts = int(t[5]); steps = int(t[7])
        g = out.splitlines()[1].split()
        N = int(g[1]); rows = [int(x, 16) for x in g[2:2 + N]]
        owner = [p_ for p_ in range(7) for _ in range(sz[p_])]
        Adj = np.array([[(rows[i] >> j) & 1 for j in range(N)] for i in range(N)], dtype=np.int64)
        st = _Lagrangian(np.concatenate([np.eye(N, dtype=np.int64), Adj], axis=1), owner)
        S = st.entropies()
        f = int(transfer_form() @ S)
        c = compress(st)
        vals, bad = transfer_chain(c) if c.N else (dict(f12=f), [])
        res.append(dict(sizes=list(sz), min_f12=minF, objective=obj, restarts=restarts, steps=steps, f12_recomputed=f,
                        chain_failed=bad, singles=[int(S[1 << p_]) for p_ in range(7)]))
        log(f'  sizes {list(sz)}: min f12 {minF} (recomputed {f}), objective {obj}, {restarts} restarts, {steps} moves; '
            f'S(P) at the minimiser {res[-1]["singles"]}; proof steps there: {"all hold" if not bad else bad}')
    return res


def transfer_check(path=DATA, full_check=False, trials=40, seed=42, log=print):
    """Theorem 12: F(S_NF) = -1 and F(S_F) = 0; F >= 0 on every vector of the GF(2) list (all seven-qubit states with
    at most one qubit per system, every relabelling); on the GF(3) list F < 0 only at S_NF.  With full_check: every
    step of the proof (transfer_chain) on random compressed qubit graph states (1-3 qubits per system), on compressed
    CSS states of random GF(2) families pushed towards the non-Fano dependencies, and on one and two copies of the
    Fano state after random local symplectic maps (where f12 = 0); the dual arrangement built explicitly for the first
    eight states (Lemma T5); and compression (symplectic reduction) against the entropies of uncompressed states."""
    F = transfer_form()
    ks = known_states()
    Z = load(path)
    out = dict(f12_nonfano=int(F @ ks['NF']), f12_fano=int(F @ ks['F']), terms=len(transfer_terms(F)))
    g2 = np.array([F @ full(v) for v in Z['gf2'].astype(np.int64)])
    g3 = np.array([F @ full(v) for v in Z['gf3'].astype(np.int64)])
    nf = tuple(int(x) for x in ks['NF'][1:64])
    out.update(gf2_min=int(g2.min()), gf2_zero=int((g2 == 0).sum()), gf3_negative=int((g3 < 0).sum()),
               gf3_negative_is_nonfano=all(tuple(int(x) for x in v) == nf for v, f in zip(Z['gf3'], g3) if f < 0))
    if full_check:
        rng = np.random.default_rng(seed)
        failed, done, minF, tight = 0, 0, None, 0
        dual_ok = 0
        kinds = {'graph': 0, 'structured': 0, 'fano': 0}
        for t in range(trials):
            if t % 4 == 3:
                st = _fano_copies(rng, 1 + (t // 4) % 2)
                S = st.entropies(); kinds['fano'] += 1
            elif t % 4 == 1:
                st = _random_structured_css(rng)
                S = st.entropies(); kinds['structured'] += 1
            else:
                st, S = _random_compressed_graph_state(rng); kinds['graph'] += 1
            vals, bad = transfer_chain(st, S)
            done += 1
            failed += bool(bad)
            minF = vals['f12'] if minF is None else min(minF, vals['f12'])
            tight += vals['f12'] == 0
            if t < 8:
                dual_ok += _dual_identity(st)
        comp_ok = 0
        for t in range(8):                                         # compression keeps the entropies, S(P) qubits each
            sizes = rng.integers(1, 4, size=7)
            owner = [p for p in range(7) for _ in range(sizes[p])]
            N = len(owner)
            U = (rng.random((N, N)) < rng.uniform(0.05, 0.5)).astype(np.int64)
            Adj = np.triu(U, 1)
            st = _Lagrangian(np.concatenate([np.eye(N, dtype=np.int64), Adj + Adj.T], axis=1), owner)
            S = st.entropies()
            c = compress(st)
            comp_ok += bool(c.N and (c.entropies() == S).all() and all(c.owner.count(P) == S[1 << P] for P in range(7)))
        out.update(chain_states=done, chain_kinds=kinds, chain_failed=failed, chain_min_f12=minF, chain_f12_zero=tight,
                   dual_identity=f'{dual_ok}/{min(trials, 8)}', compression=f'{comp_ok}/8')
    log(f'  transfer: {out}')
    return out


def transfer_ok(o, full_check=False):
    ok = (o['f12_nonfano'] == -1 and o['f12_fano'] == 0 and o['gf2_min'] >= 0 and o['gf3_negative'] == 1
          and o['gf3_negative_is_nonfano'])
    if full_check:
        n8 = min(o['chain_states'], 8)
        ok = (ok and o['chain_failed'] == 0 and o['chain_min_f12'] >= 0 and o['dual_identity'] == f'{n8}/{n8}'
              and o['compression'] == '8/8')
    return ok


def self_test(trials=30, seed=34, log=print):
    """(1) the Steane state: entropies from its state vector equal S of h_F; (2) the cut-rank formula against state
    vectors of random qubit graph states with two qubits per system (14 qubits) and of random qutrit graph states with
    one qutrit per system; (3) scaling a vertex by -1 keeps every cut rank over GF(3); (4) line_pieces_check
    (Proposition 10, steps 1 and 1', on the Fano and non-Fano CSS states with one and two qudits per system)."""
    import s31_six_party as s31
    ks = known_states()
    Ssv = s31.entropies(s31.simplex_state())
    ok1 = all(abs(Ssv[X] - ks['F'][X]) < 1e-9 for X in range(1, 128))

    def sv_entropies(A, sizes, d):
        owner = [q for q, c in enumerate(sizes) for _ in range(c)]
        N = len(owner)
        xs = np.array(list(itertools.product(range(d), repeat=N)), dtype=np.int64)
        ph = np.zeros(len(xs), dtype=np.int64)
        for i in range(N):
            for j in range(i + 1, N):
                if A[i][j]:
                    ph += A[i][j] * xs[:, i] * xs[:, j]
        psi = np.exp(2j * np.pi * (ph % d) / d) / np.sqrt(float(d) ** N)
        t = psi.reshape([d] * N)
        out = []
        for m in range(1, 64):
            X = [v for v in range(N) if owner[v] < 6 and m >> owner[v] & 1]
            Y = [v for v in range(N) if v not in X]
            M = np.transpose(t, X + Y).reshape(d ** len(X), d ** len(Y))
            sv = np.linalg.svd(M, compute_uv=False) ** 2
            sv = sv[sv > 1e-12]
            out.append(float(-(sv * np.log(sv)).sum() / np.log(d)))
        return out

    rng = random.Random(seed)
    ok2 = True
    for t in range(trials):
        if t % 2 == 0:
            sizes, d = [2] * 7, 2
        else:
            sizes, d = [1] * 7, 3
        N = sum(sizes)
        owner = [q for q, c in enumerate(sizes) for _ in range(c)]
        A = [[0] * N for _ in range(N)]
        for i in range(N):
            for j in range(i + 1, N):
                if owner[i] != owner[j]:
                    A[i][j] = A[j][i] = rng.randrange(d)
        cr = cut_vector(A, d, sizes)
        ent = sv_entropies(A, sizes, d)
        ok2 &= all(abs(cr[m] - ent[m]) < 1e-8 for m in range(63))
    ok3 = True
    for t in range(200):
        A = [[0] * 7 for _ in range(7)]
        for (u, v) in EDGES:
            A[u][v] = A[v][u] = rng.randrange(3)
        B = [row[:] for row in A]
        x = rng.randrange(7)
        for j in range(7):
            B[x][j] = (-B[x][j]) % 3; B[j][x] = (-B[j][x]) % 3
        ok3 &= cut_vector(A, 3) == cut_vector(B, 3)
    ok4 = line_pieces_check()['ok']
    log(f'  self-test: Steane entropies {ok1}; cut-rank formula vs state vectors ({trials} graph states) {ok2}; '
        f'vertex scaling over GF(3) {ok3}; line pieces (Proposition 10, steps 1 and 1\') {ok4}')
    return ok1 and ok2 and ok3 and ok4


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--check', action='store_true')
    ap.add_argument('--full', action='store_true')
    ap.add_argument('--self-test', action='store_true')
    ap.add_argument('--enumerate', action='store_true')
    ap.add_argument('--survey', action='store_true')
    ap.add_argument('--nonfano', action='store_true')
    ap.add_argument('--lists', default=None)
    ap.add_argument('--out', default=None)
    ap.add_argument('--seconds', type=float, default=60)
    ap.add_argument('--lams', default='2,3')
    ap.add_argument('--fields', default=None, help='comma-separated q: cc/s36_target7.c for S_F, S_NF and controls')
    ap.add_argument('--controls', type=int, default=2)
    ap.add_argument('--fields-python', default=None, help='comma-separated q: the Python search for S_F and S_NF')
    ap.add_argument('--transfer', action='store_true', help='Theorem 12: the inequality f12 and its checks (--full: proof steps)')
    ap.add_argument('--trials', type=int, default=40)
    ap.add_argument('--transfer-anneal', default=None,
                    help='semicolon-separated qubit counts per system, e.g. "2,2,2,2,2,2,2;3,3,3,3,3,3,3": cc/s37_minform7.c')
    ap.add_argument('--pull-nf', type=float, default=0.0, help='with --transfer-anneal: weight of a pull towards m S_NF')
    ap.add_argument('--seed', type=int, default=4201)
    a = ap.parse_args()
    if a.transfer_anneal:
        cf = [[int(x) for x in c.split(',')] for c in a.transfer_anneal.split(';')]
        tg = None
        if a.pull_nf > 0:
            ks = known_states()
            m_ = cf[0][0]
            tg = [m_ * int(ks['NF'][x]) for x in range(1, 64)]
        r = transfer_anneal(cf, seconds=a.seconds, seed=a.seed, mu=a.pull_nf, target=tg)
        print('ALL OK' if all(x['min_f12'] >= 0 and x['f12_recomputed'] == x['min_f12'] and not x['chain_failed'] for x in r)
              else 'FAILURE')
    if a.transfer:
        o = transfer_check(full_check=a.full, trials=a.trials)
        print('terms (Z = purifier):', ' '.join(f'{c:+d}*S({n})' for n, c in transfer_terms()))
        print('ALL OK' if transfer_ok(o, a.full) else 'FAILURE')
    if a.fields:
        print(json.dumps(field_runs([int(x) for x in a.fields.split(',')], a.controls)))
    if a.fields_python:
        print(json.dumps(fields_check([int(x) for x in a.fields_python.split(',')])))
    if a.self_test:
        print('ALL OK' if self_test() else 'FAILURE')
    if a.enumerate:
        enumerate_lists(a.out)
    if a.survey:
        print(survey(a.lists, a.out, a.seconds))
    if a.nonfano:
        nonfano_runs([int(x) for x in a.lams.split(',')], a.seconds)
    if a.check:
        s = check(full_check=a.full)
        print(json.dumps(s))
        print('ALL OK' if check_ok(s, a.full) else 'FAILURE')


if __name__ == '__main__':
    main()
