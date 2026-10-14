#!/usr/bin/env python3
"""s34: six parties, the other direction -- states with at most one qudit per system (R40, R41; roadmap item M3).

Question (the paper, Discussion; R32 note section 5): is Stab^(p)_6 not contained in Stab^(2)_6 for some odd p, i.e. is
there a six-party stabilizer state of odd local dimension whose entropy vector lies outside the qubit cone?  The only
candidate so far is the non-Fano state h_NF = r_NF + r_NF* (section 4.8); the transfer argument that settles the other
direction for the Steane state does not work for it.

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
"""
import argparse, itertools, json, os, random, shutil, subprocess, sys, tempfile, time
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
    a = ap.parse_args()
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
