#!/usr/bin/env python3
"""s36: six parties, one qudit per system over other fields (R45; docs/R32_reduction_and_local_dimension.md, section
4.10, "Other fields", and section 5, two odd primes at n = 6).

Over GF(q), q = p^k, the cut-rank vector of a weighted graph on seven vertices (0..5 visible, 6 the purifier) is the
entropy vector, in units of log q, of a GF(q)-linear stabilizer state with one GF(q)-qudit per system.  Restriction of
scalars (the trace form tr(omega) and a trace-dual basis of GF(q) over GF(p)) turns that state into a stabilizer state
over Z_p with k qudits per system whose entropies in units of log p are k times these; so the vector lies in the
cone Stab^(p)_6.  For prime q these vectors are all entropy vectors of Z_q stabilizer states of seven systems with every
S(P) <= 1 (scripts/s34_six_party_survey.py).

Data: data/s36_one_qudit_fields.npz
  reps      the 136 S_7 orbits met over the 13 fields below (uint8, 63 cut ranks; the lexicographically smallest vector
            of each orbit, sorted)
  sizes     orbit sizes (int64; 36,430 vectors in all)
  fields    2 3 4 5 7 8 9 11 13 16 25 27 49 (int64)
  cert      (13, 136, 21) uint8: for each field and orbit the 21 weights (edge order (0,1), (0,2), ..., (5,6)) of a graph
            over GF(q) whose cut-rank vector is the representative itself; 255 = no graph found
  steane, nonfano   bool (136): the orbit of the Steane vector S_F and of the non-Fano vector S_NF (section 4.8)
How the lists were obtained: GF(2), GF(3) exhaustively (cc/s35_enum7.c, stored in data/s34_seven_qudit_vectors.npz);
GF(4), GF(5) exhaustively with cc/s40_enum7q.c (gate G41b re-runs it); the other fields by sampling random weighted graphs
with cc/s39_sample7q.c (results/2026-10-15_r45-sandbox/), so their lists are complete only empirically.  Field elements
are encoded as in the C programs: base-p digits of the coefficients of 1, x, x^2, ... modulo the polynomials in MODULI.
The C programs are numbered separately from the scripts (cc/s39, cc/s40 belong to this script).

Proposition 15 (the U_{2,n} states need q >= n - 1) is checked through exhaustive target searches (cc/s36_target7.c),
Proposition 16 (the eight-element analogue of the Fano pair: [I_4 | J_4 - I_4], T8, T8', AG(3,2)) through its 70 minors,
a brute-force representability count over GF(2) .. GF(9), and the entropy vectors of the CSS states of T8 and T8'.

  python3 scripts/s36_one_qudit_fields.py --check [--full]    # gates G41 / G41b
  python3 scripts/s36_one_qudit_fields.py --enumerate 5 [--all] [--slice 0 --nslices 1] [--out e5.npz]
        # cc/s40_enum7q.c over GF(q), q in 2 3 4 5 7 8 9: orbits against the table (GF(5): ~5 min on two cores;
        # GF(7): 815,738,564,501 graphs, about 15-30 CPU-hours at the measured rates -- not run)
  python3 scripts/s36_one_qudit_fields.py --sample 7 --samples 1000000 --seed 1 --out s_gf7   # cc/s39_sample7q.c
  python3 scripts/s36_one_qudit_fields.py --build s_gf --out data.npz   # the table from sampler outputs s_gf<q>.bin/.w
  python3 scripts/s36_one_qudit_fields.py --prop16             # Proposition 16: minors, brute force, S_T8 against S_T8'
"""
import argparse, itertools, os, shutil, subprocess, sys, tempfile, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, '..')
sys.path.insert(0, HERE); sys.path.insert(0, ROOT)
import s34_six_party_survey as s34

DATA = os.path.join(ROOT, 'data', 's36_one_qudit_fields.npz')
FIELDS = (2, 3, 4, 5, 7, 8, 9, 11, 13, 16, 25, 27, 49)
EXHAUSTIVE = (2, 3, 4, 5)
MODULI = {4: (2, (1, 1)), 8: (2, (1, 1, 0)), 16: (2, (1, 1, 0, 0)), 9: (3, (1, 0)), 25: (5, (3, 0)),
          27: (3, (2, 2, 0)), 49: (7, (1, 0))}   # (p, (c_0, .., c_{k-1})): x^k + c_{k-1} x^{k-1} + .. + c_0
EDGES = s34.EDGES
SAMPLES = {2: 2 * 10 ** 6, 3: 4 * 10 ** 6, 4: 2 * 10 ** 7, 5: 10 ** 8, 7: 10 ** 8, 8: 10 ** 8, 9: 10 ** 8, 11: 5 * 10 ** 7,
           13: 5 * 10 ** 7, 16: 10 ** 8, 25: 10 ** 8, 27: 10 ** 8, 49: 5 * 10 ** 7}   # graphs sampled for the stored table
SEEDS = {2: 201, 3: 301, 4: 401, 5: 501, 7: 701, 8: 801, 9: 901, 11: 1101, 13: 1301, 16: 1601, 25: 2501, 27: 2701,
         49: 4901}   # their seeds (cc/s39_sample7q.c q SAMPLES[q] SEEDS[q] s_gf<q>, then --build s_gf)


# ----------------------------------------------------------------------------------------------- GF(q) from scratch
def char_deg(q):
    for p in (2, 3, 5, 7, 11, 13):
        k, x = 0, q
        while x % p == 0:
            x //= p
            k += 1
        if x == 1 and k:
            return p, k
    raise ValueError(q)


class GF:
    """GF(q) by tables built from polynomial arithmetic over GF(p) (independent of the C programs); the field axioms
    are checked on construction."""

    def __init__(self, q):
        self.q = q
        self.p, self.k = char_deg(q)
        p, k = self.p, self.k
        c = MODULI[q][1] if k > 1 else ()

        def vec(a):
            return [(a // p ** i) % p for i in range(k)]

        def num(v):
            return sum(int(x) * p ** i for i, x in enumerate(v))

        def mul(a, b):
            if k == 1:
                return a * b % p
            va, vb = vec(a), vec(b)
            prod = [0] * (2 * k - 1)
            for i, x in enumerate(va):
                for j, y in enumerate(vb):
                    prod[i + j] = (prod[i + j] + x * y) % p
            for t in range(2 * k - 2, k - 1, -1):                 # x^k = -(c_0 + c_1 x + ... + c_{k-1} x^{k-1})
                f = prod[t]
                if f:
                    prod[t] = 0
                    for i in range(k):
                        prod[t - k + i] = (prod[t - k + i] - f * c[i]) % p
            return num(prod[:k])

        R = range(q)
        self.add = [[num([(x + y) % p for x, y in zip(vec(a), vec(b))]) for b in R] for a in R]
        self.sub = [[num([(x - y) % p for x, y in zip(vec(a), vec(b))]) for b in R] for a in R]
        self.mul = [[mul(a, b) for b in R] for a in R]
        self.inv = [0] * q
        for a in range(1, q):
            inv = [b for b in R if self.mul[a][b] == 1]
            if len(inv) != 1:
                raise ValueError(f'GF({q}): {a} has no inverse (the modulus is reducible)')
            self.inv[a] = inv[0]
        ad, mu = self.add, self.mul
        assert all(mu[a][mu[b][e]] == mu[mu[a][b]][e] and mu[a][ad[b][e]] == ad[mu[a][b]][mu[a][e]]
                   and ad[a][ad[b][e]] == ad[ad[a][b]][e] for a in R for b in R for e in R)
        assert all(ad[a][b] == ad[b][a] and mu[a][b] == mu[b][a] for a in R for b in R)

    def rank(self, M):
        M = [list(r) for r in M]
        r, cols = 0, (len(M[0]) if M else 0)
        mul, sub, inv = self.mul, self.sub, self.inv
        for c in range(cols):
            piv = next((i for i in range(r, len(M)) if M[i][c]), None)
            if piv is None:
                continue
            M[r], M[piv] = M[piv], M[r]
            iv = inv[M[r][c]]
            M[r] = [mul[iv][x] for x in M[r]]
            for i in range(len(M)):
                if i != r and M[i][c]:
                    f = M[i][c]
                    M[i] = [sub[x][mul[f][y]] for x, y in zip(M[i], M[r])]
            r += 1
        return r


_FIELDS = {}


def field(q):
    if q not in _FIELDS:
        _FIELDS[q] = GF(q)
    return _FIELDS[q]


def cut_vector(F, w):
    """the 63 cut ranks over F of the graph on seven vertices with the 21 weights w (entry m-1: the cut whose visible
    systems are the bits of m)."""
    A = [[0] * 7 for _ in range(7)]
    for (u, v), x in zip(EDGES, w):
        A[u][v] = A[v][u] = int(x)
    out = []
    for m in range(1, 64):
        X = [v for v in range(6) if m >> v & 1]
        Y = [v for v in range(7) if v not in X]
        rows, cols = (X, Y) if len(X) <= len(Y) else (Y, X)
        out.append(F.rank([[A[i][j] for j in cols] for i in rows]))
    return tuple(out)


# ----------------------------------------------------------------------------------------------- orbits
def load(path=DATA):
    return np.load(path, allow_pickle=False)


def content_sha(path=DATA):
    return s34.content_sha(path)


def orbit_map(reps):
    """{vector: orbit index} over the orbits of the representatives, and the orbits themselves."""
    idx, orbs = {}, []
    for i, r in enumerate(reps):
        o = s34.orbit63(np.array(r, dtype=np.int64))
        orbs.append(o)
        for v in o:
            idx[v] = i
    return idx, orbs


def special_vectors():
    """{'symmetric': S = 1 on singles and six-sets, 2 otherwise; 'AME': S = min(|X|, 7 - |X|)} on the 63 cuts."""
    pc = [bin(m).count('1') for m in range(1, 64)]
    return {'symmetric': tuple(1 if c in (1, 6) else 2 for c in pc), 'AME': tuple(min(c, 7 - c) for c in pc)}


def u26_orbit(reps, sizes):
    """the orbit of the vector with one unentangled system and S = min(|X|, 2, 6 - |X|) on the other six (U_{2,6})."""
    for i, r in enumerate(reps):
        S = s34.full(np.array(r, dtype=np.int64))
        zero = [P for P in range(7) if S[1 << P] == 0]
        if len(zero) != 1:
            continue
        rest = [P for P in range(7) if P != zero[0]]
        if all(int(S[sum(1 << P for P in c)]) == min(len(c), 2, 6 - len(c))
               for k in (1, 2, 3) for c in itertools.combinations(rest, k)):
            return i
    return None


def check(path=DATA, full_check=False, log=print):
    """the table (gate G41) and, with full_check, the exhaustive enumerations, the target searches of Proposition 15,
    a sampling run and Proposition 16 (gate G41b).  Returns a summary dict."""
    t0 = time.time()
    Z = load(path)
    reps = [tuple(int(x) for x in r) for r in Z['reps']]
    sizes = [int(x) for x in Z['sizes']]
    fields = [int(x) for x in Z['fields']]
    cert = Z['cert']
    o = {}
    idx, orbs = orbit_map(reps)
    o['orbits'] = len(reps)
    o['vectors'] = len(idx)
    o['reps_canonical'] = all(r == min(orbs[i]) for i, r in enumerate(reps)) and reps == sorted(reps)
    o['sizes_ok'] = [len(x) for x in orbs] == sizes and sum(sizes) == len(idx)
    have = {q: [i for i in range(len(reps)) if cert[fi, i, 0] != 255] for fi, q in enumerate(fields)}
    o['per_field'] = {q: len(v) for q, v in have.items()}
    o['per_field_vectors'] = {q: sum(sizes[i] for i in v) for q, v in have.items()}
    # every certificate: in range, and its graph realises the representative exactly over GF(q)
    bad = 0
    for fi, q in enumerate(fields):
        F = field(q)
        for i in have[q]:
            w = cert[fi, i]
            bad += bool((w >= q).any()) or cut_vector(F, w) != reps[i]
    o['certificates'] = sum(len(v) for v in have.values())
    o['certificates_bad'] = bad
    # GF(2) and GF(3): the certified orbits are exactly the orbits of the exhaustive lists of section 4.10
    D = s34.load()
    for q in (2, 3):
        L = {tuple(int(x) for x in r) for r in D[f'gf{q}']}
        o[f'gf{q}_list_is_union'] = (all(v in idx for v in L) and {idx[v] for v in L} == set(have[q])
                                     and len(L) == o['per_field_vectors'][q])
    # the Steane and non-Fano orbits; the symmetric, AME and U_{2,6} orbits
    ks = s34.known_states()
    st = tuple(int(x) for x in ks['F'][1:64]); nf = tuple(int(x) for x in ks['NF'][1:64])
    iS, iN = idx[st], idx[nf]
    o['steane'], o['nonfano'] = iS, iN
    o['flags_ok'] = (list(np.flatnonzero(Z['steane'])) == [iS] and list(np.flatnonzero(Z['nonfano'])) == [iN])
    sv = special_vectors()
    i_sym, i_ame, i_u26 = idx[sv['symmetric']], idx[sv['AME']], u26_orbit(reps, sizes)
    o['symmetric_fields'] = [q for q in fields if i_sym in have[q]]
    o['ame_fields'] = [q for q in fields if i_ame in have[q]]
    o['u26_fields'] = [q for q in fields if i_u26 is not None and i_u26 in have[q]]
    # only 2 against odd: over the fields with at least seven elements the certified orbits are 134 common ones plus
    # the non-Fano orbit (odd characteristic) or the Steane orbit (characteristic 2)
    big = [q for q in fields if q >= 7]
    common = set.intersection(*(set(have[q]) for q in big))
    o['common'] = len(common)
    o['big_fields_structure'] = all(set(have[q]) == common | ({iS} if q % 2 == 0 else {iN}) for q in big)
    # per characteristic: the union over the tested fields; every orbit missing over a field reappears over a larger
    # field of the same characteristic
    chars = sorted({char_deg(q)[0] for q in fields})
    union = {p: set().union(*(set(have[q]) for q in fields if char_deg(q)[0] == p)) for p in chars}
    o['char_unions'] = {p: len(u) for p, u in union.items()}
    o['odd_unions_equal'] = len({frozenset(union[p]) for p in chars if p > 2}) == 1
    o['char2_vs_odd'] = (sorted(union[2] - union[3]), sorted(union[3] - union[2]))
    reappear = {}
    for q in fields:
        p = char_deg(q)[0]
        miss = union[p] - set(have[q])
        larger = set().union(*(set(have[q2]) for q2 in fields if q2 > q and char_deg(q2)[0] == p)) if miss else set()
        reappear[q] = (len(miss), miss <= larger)
    o['reappear'] = {q: v for q, v in reappear.items() if v[0]}
    # pairs of tested odd characteristics: every orbit realised in one is realised in the other; for p = 3 the GF(3)
    # list is exhaustive and each of its orbits is certified over GF(q), q = 5, 7, 11, 13 themselves
    o['gf3_in_prime_fields'] = {q: set(have[3]) <= set(have[q]) for q in (5, 7, 11, 13)}
    o['in_char3'] = {q: set(have[q]) <= union[3] for q in (5, 7, 11, 13)}
    odd = (3, 5, 7, 11, 13)                 # every orbit met over GF(p) is certified in characteristic q
    o['pairs'] = {f'{p}-{q}': set(have[p]) <= union[q] for p in odd for q in odd if p != q}
    if full_check:
        o.update(enumeration_check(have, idx, reps, log=log))
        o.update(target_check(reps, i_sym, i_u26, log=log))
        o.update(sampling_check(have, idx, log=log))
        o.update(prop16_check(log=log))
    o['seconds'] = round(time.time() - t0, 1)
    log(f'  s36: {o}')
    return o


def check_ok(o, full_check=False):
    ok = (o['orbits'] == 136 and o['vectors'] == 36430 and o['reps_canonical'] and o['sizes_ok']
          and o['certificates_bad'] == 0 and o['certificates'] == 1646 and o['gf2_list_is_union'] and o['gf3_list_is_union']
          and o['per_field'] == {2: 59, 3: 111, 4: 127, 5: 134, 7: 135, 8: 135, 9: 135, 11: 135, 13: 135, 16: 135,
                                 25: 135, 27: 135, 49: 135}
          and o['flags_ok'] and o['common'] == 134 and o['big_fields_structure'] and o['odd_unions_equal']
          and o['char2_vs_odd'] == ([o['steane']], [o['nonfano']])
          and o['reappear'] == {2: (76, True), 3: (24, True), 4: (8, True), 5: (1, True)}
          and all(o['gf3_in_prime_fields'].values()) and all(o['in_char3'].values()) and len(o['pairs']) == 20
          and all(o['pairs'].values())
          and o['symmetric_fields'] == [7, 8, 9, 11, 13, 16, 25, 27, 49]
          and o['u26_fields'] == [5, 7, 8, 9, 11, 13, 16, 25, 27, 49]
          and o['ame_fields'] == [3, 4, 5, 7, 8, 9, 11, 13, 16, 25, 27, 49])
    if full_check:
        ok = ok and all(v['orbits_equal'] and v['outside'] == 0 and v['certificates_bad'] == 0
                        for v in o['enumeration'].values())
        ok = ok and o['enumeration']['GF(2) all']['list_equal'] and o['enumeration']['GF(3)']['list_equal']
        ok = ok and o['targets'] == {'symmetric GF(4)': 0, 'symmetric GF(5)': 0, 'symmetric GF(7)': 1,
                                     'U26 GF(4)': 0, 'U26 GF(5)': 1}
        ok = ok and len(o['sampling']) == 3 and all(v['outside'] == 0 and v['certificates_bad'] == 0 and v['orbits'] >= 120
                                                    for v in o['sampling'].values())
        ok = ok and o['prop16_ok']
    return ok


# ----------------------------------------------------------------------------------------------- the C programs
def run_enum(q, threads=None, all_=False, slice_=0, nslices=1):
    """cc/s40_enum7q.c over GF(q): (vectors, weights, totals)."""
    tmp = tempfile.mkdtemp(prefix='s40_')
    try:
        exe = s34.build('s40_enum7q', tmp)
        args = [exe, str(q), str(threads or os.cpu_count() or 1), 'all' if all_ else 'reduced']
        if nslices > 1:
            args += ['slice', str(slice_), str(nslices)]
        p = subprocess.run(args, capture_output=True, text=True)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    if p.returncode != 0:
        raise RuntimeError(p.stderr)
    V, W, tot = [], [], None
    for ln in p.stdout.splitlines():
        t = ln.split()
        if t[0] == 'V':
            x = [int(a) for a in t[1:]]
            V.append(tuple(x[:63])); W.append(x[63:])
        elif t[0] == 'TOTAL':
            tot = dict(prefixes=int(t[2]), graphs=int(t[3]), distinct=int(t[4]))
    return V, W, tot


def enumeration_check(have, idx, reps, recheck=400, seed=41, log=print):
    """exhaustive lists (cc/s40_enum7q.c): GF(2) without the relabelling reduction must give the stored list itself;
    GF(2), GF(3), GF(4), GF(5) with it, closed under S_7, the orbits certified in the table (and for GF(3) the stored
    list).  Up to `recheck` certificates of the program per field are re-checked here."""
    D = s34.load()
    rng = np.random.default_rng(seed)
    out = {}
    for name, q, all_ in (('GF(2) all', 2, True), ('GF(2)', 2, False), ('GF(3)', 3, False), ('GF(4)', 4, False),
                          ('GF(5)', 5, False)):
        t0 = time.time()
        V, W, tot = run_enum(q, all_=all_)
        F = field(q)
        pick = rng.choice(len(V), size=min(recheck, len(V)), replace=False)
        badc = sum(cut_vector(F, W[i]) != V[i] for i in pick)
        orb = {idx[v] for v in V if v in idx}
        res = dict(tot, outside=sum(v not in idx for v in V), orbits=len(orb), orbits_equal=orb == set(have[q]),
                   certificates_bad=int(badc), certificates_checked=len(pick), seconds=round(time.time() - t0, 1))
        if q in (2, 3):
            L = {tuple(int(x) for x in r) for r in D[f'gf{q}']}
            res['list_equal'] = (set(V) == L) if all_ else (set().union(*(s34.orbit63(np.array(reps[i])) for i in orb)) == L)
        out[name] = res
        log(f'  enumeration {name}: {res}')
    return {'enumeration': out}


def run_target(t63, q, max_print=1, stop_after=0):
    tmp = tempfile.mkdtemp(prefix='s36t_')
    try:
        exe = s34.build('s36_target7', tmp)
        args = [exe, str(q), str(max_print)] + ([str(stop_after)] if stop_after else [])
        p = subprocess.run(args, input=' '.join(str(int(x)) for x in t63) + '\n', capture_output=True, text=True)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    sols = [[int(x) for x in ln.split()[1:]] for ln in p.stdout.splitlines() if ln.startswith('SOL')]
    total = int(next(ln.split()[2] for ln in p.stdout.splitlines() if ln.startswith('TOTAL')))
    return total, sols


def target_check(reps, i_sym, i_u26, log=print):
    """exhaustive target searches (cc/s36_target7.c): the symmetric vector over GF(4), GF(5) (no graph) and GF(7)
    (stop at the first); the U_{2,6} orbit representative over GF(4) (none) and GF(5) (stop at the first).  Every
    solution printed is re-checked here."""
    sv = special_vectors()['symmetric']
    res, ok = {}, True
    for name, t, q, stop in (('symmetric GF(4)', sv, 4, 0), ('symmetric GF(5)', sv, 5, 0), ('symmetric GF(7)', sv, 7, 1),
                             ('U26 GF(4)', reps[i_u26], 4, 0), ('U26 GF(5)', reps[i_u26], 5, 1)):
        total, sols = run_target(t, q, 1, stop)
        ok = ok and all(cut_vector(field(q), s) == tuple(t) for s in sols)
        res[name] = total
    log(f'  target searches: {res}, solutions re-checked {ok}')
    return {'targets': res if ok else {}}


def run_sample(q, samples, seed, prefix):
    tmp = tempfile.mkdtemp(prefix='s39_')
    try:
        exe = s34.build('s39_sample7q', tmp)
        p = subprocess.run([exe, str(q), str(samples), str(seed), prefix], capture_output=True, text=True)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    if p.returncode != 0:
        raise RuntimeError(p.stderr)
    V = np.fromfile(prefix + '.bin', dtype=np.uint8).reshape(-1, 63)
    W = np.fromfile(prefix + '.w', dtype=np.uint8).reshape(-1, 21)
    return V, W


def sampling_check(have, idx, samples=200000, recheck=200, seed=39, log=print):
    """a short run of cc/s39_sample7q.c over GF(7), GF(8), GF(9): every sampled vector lies in an orbit certified over
    that field; `recheck` of the sampler's certificates per field re-checked here."""
    rng = np.random.default_rng(seed)
    out = {}
    tmp = tempfile.mkdtemp(prefix='s39run_')
    try:
        for q in (7, 8, 9):
            V, W = run_sample(q, samples, seed + q, os.path.join(tmp, f's{q}'))
            vs = [tuple(int(x) for x in r) for r in V]
            outside = sum(v not in idx or idx[v] not in have[q] for v in vs)
            pick = rng.choice(len(vs), size=min(recheck, len(vs)), replace=False)
            badc = sum(cut_vector(field(q), W[i]) != vs[i] for i in pick)
            out[f'GF({q})'] = dict(samples=samples, distinct=len(vs), orbits=len({idx[v] for v in vs if v in idx}),
                                   outside=outside, certificates_bad=int(badc))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    log(f'  sampling: {out}')
    return {'sampling': out}


# ----------------------------------------------------------------------------------------------- Proposition 16
T8_TEN = [frozenset({('v', j) for j in range(4) if j != m} | {('w', m)}) for m in range(4)] + \
         [frozenset({('v', k), ('v', l), ('w', k), ('w', l)}) for k, l in itertools.combinations(range(4), 2)]
W_ALL = frozenset(('w', m) for m in range(4))
AG_EXTRA = [frozenset({('v', k)} | {('w', j) for j in range(4) if j != k}) for k in range(4)]
ELEMS = [('v', j) for j in range(4)] + [('w', m) for m in range(4)]


def _det_int(M):
    M = [list(map(int, r)) for r in M]
    n = len(M)
    if n == 0:
        return 1
    return sum((-1) ** c * M[0][c] * _det_int([r[:c] + r[c + 1:] for r in M[1:]]) for c in range(n))


def prop16_minors():
    """the 70 minors of [I_4 | J_4 - I_4] (columns v_1..v_4, w_1..w_4), exactly: {4-set: minor}."""
    cols = {('v', j): [int(i == j) for i in range(4)] for j in range(4)}
    cols.update({('w', m): [int(i != m) for i in range(4)] for m in range(4)})
    return {frozenset(S): _det_int([[cols[e][i] for e in S] for i in range(4)]) for S in itertools.combinations(ELEMS, 4)}


def prop16_brute(q):
    """all configurations over GF(q) with v_1..v_4 the standard basis, w_1 = (0, 1, 1, 1) and w_m (m = 2, 3, 4) with
    entry m zero, the other three nonzero, the first nonzero one 1.  No representation of T8, T8' or AG(3,2) is lost:
    in each of them {v_1..v_4} is a basis, the 4-sets {v_j : j != m} + {w_m} are dependent and no 3-set is, and scaling
    v_2, v_3, v_4 and the w_m gives this form.  A set {v_j : j in A} + {w_m : m in B} is dependent iff the submatrix
    W[rows not in A, columns B] of W = (w_1 .. w_4) has rank < |B|; all minors are computed with the tables of GF(q).
    Returns {matroid: number of configurations representing it} (q - 1 each if the representation is unique up to
    scaling: the normalisation leaves one parameter)."""
    F = field(q)
    MUL, SUB, ADD = (np.array(t, dtype=np.int64) for t in (F.mul, F.sub, F.add))
    cols = []
    for m in (1, 2, 3):
        opts = []
        for vals in itertools.product(range(1, q), repeat=3):
            if vals[0] == 1:
                w = [0] * 4
                for i, x in zip([i for i in range(4) if i != m], vals):
                    w[i] = x
                opts.append(w)
        cols.append(np.array(opts, dtype=np.int64))
    n = len(cols[0])
    g = np.stack(np.meshgrid(np.arange(n), np.arange(n), np.arange(n), indexing='ij'), -1).reshape(-1, 3)
    N = len(g)
    W = np.zeros((N, 4, 4), dtype=np.int64)
    W[:, 1:, 0] = 1
    for m in (1, 2, 3):
        W[:, :, m] = cols[m - 1][g[:, m - 1]]
    cache = {}

    def det(R, C):
        key = (R, C)
        if key not in cache:
            if len(R) == 1:
                cache[key] = W[:, R[0], C[0]]
            else:
                acc = np.zeros(N, dtype=np.int64)
                for j, c in enumerate(C):
                    t = MUL[W[:, R[0], c], det(R[1:], C[:j] + C[j + 1:])]
                    acc = ADD[acc, t] if j % 2 == 0 else SUB[acc, t]
                cache[key] = acc
        return cache[key]

    def dependent(S):
        A = tuple(j for (t, j) in S if t == 'v')
        B = tuple(sorted(m for (t, m) in S if t == 'w'))
        if not B:
            return np.zeros(N, dtype=bool)
        rows = tuple(i for i in range(4) if i not in A)
        dep = np.ones(N, dtype=bool)
        for R in itertools.combinations(rows, len(B)):
            dep &= det(R, B) == 0
        return dep

    ok3 = np.ones(N, dtype=bool)
    for S in itertools.combinations(ELEMS, 3):
        ok3 &= ~dependent(S)
    dep4 = {frozenset(S): dependent(S) for S in itertools.combinations(ELEMS, 4)}
    targets = {'T8': set(T8_TEN) | {W_ALL}, "T8'": set(T8_TEN), 'AG(3,2)': set(T8_TEN) | set(AG_EXTRA)}
    count = {}
    for k, t in targets.items():
        match = ok3.copy()
        for S, d in dep4.items():
            match &= d if S in t else ~d
        count[k] = int(match.sum())
    return count


def css_vector(q):
    """the entropy vector, in units of log q, of the CSS state of the columns of [I_4 | J_4 - I_4] over GF(q), on the
    subsets of the eight systems v_1..v_4 (bits 0..3), w_1..w_4 (bits 4..7): S(X) = r(X) + r(E - X) - r(E) (section 2,
    step 2)."""
    F = field(q)
    cols = [[int(i == j) for i in range(4)] for j in range(4)] + [[int(i != m) for i in range(4)] for m in range(4)]
    r = [F.rank([cols[e] for e in range(8) if X >> e & 1]) if X else 0 for X in range(256)]
    return [r[X] + r[255 ^ X] - r[255] for X in range(256)]


def prop16_check(fields=(2, 3, 4, 5, 7, 8, 9), log=print):
    mn = prop16_minors()
    zero = {S for S, d in mn.items() if d == 0}
    two = {S for S, d in mn.items() if abs(d) == 2}
    three = {S: d for S, d in mn.items() if abs(d) == 3}
    ones = sum(abs(d) == 1 for d in mn.values())
    minors_ok = (zero == set(T8_TEN) and two == set(AG_EXTRA) and three == {W_ALL: -3} and ones == 55
                 and len(mn) == 70)
    brute = {q: prop16_brute(q) for q in fields}
    exp = {q: {'T8': (q - 1) * (q % 3 == 0), "T8'": (q - 1) * (q % 2 != 0 and q % 3 != 0),
               'AG(3,2)': (q - 1) * (q % 2 == 0)} for q in fields}
    s3, s5 = css_vector(3), css_vector(5)                        # S_T8 (GF(3)) and S_T8' (GF(5))
    diff = {X: (s3[X], s5[X]) for X in range(256) if s3[X] != s5[X]}
    one_cut = diff == {0x0F: (3, 4), 0xF0: (3, 4)}
    log(f'  Proposition 16: minors {len(mn)}: zero at the ten sets {zero == set(T8_TEN)}, +-2 at the four AG(3,2) '
        f'sets {two == set(AG_EXTRA)}, -3 at {{w1..w4}} {three == {W_ALL: -3}}, {ones} at +-1; brute force {brute}; '
        f'S_T8 and S_T8\' differ only on the cut {{v}}|{{w}} (3 against 4): {one_cut}')
    return {'prop16_minors_ok': minors_ok, 'prop16_brute': brute, 'prop16_one_cut': one_cut,
            'prop16_ok': minors_ok and brute == exp and one_cut}


# ----------------------------------------------------------------------------------------------- building the table
def build_table(prefix, out, log=print):
    """the table from sampler / enumerator outputs <prefix><q>.bin and <prefix><q>.w (q in FIELDS): orbits of the union,
    representatives, and per field a certificate for each orbit met (a realising graph relabelled to realise the
    representative, re-checked here)."""
    perms = list(itertools.permutations(range(7)))
    P = s34._perms()
    data = {}
    for q in FIELDS:
        V = np.fromfile(f'{prefix}{q}.bin', dtype=np.uint8).reshape(-1, 63).astype(np.int64)
        W = np.fromfile(f'{prefix}{q}.w', dtype=np.uint8).reshape(-1, 21)
        data[q] = (V, W)
    rep_of, reps = {}, []
    for q in FIELDS:
        for v in data[q][0]:
            t = tuple(int(x) for x in v)
            if t in rep_of:
                continue
            o = s34.orbit63(v)
            r = min(o)
            reps.append(r)
            for u in o:
                rep_of[u] = r
    reps.sort()
    index = {r: i for i, r in enumerate(reps)}
    sizes = [len(s34.orbit63(np.array(r))) for r in reps]
    cert = np.full((len(FIELDS), len(reps), 21), 255, dtype=np.uint8)
    for fi, q in enumerate(FIELDS):
        V, W = data[q]
        F = field(q)
        for v, w in zip(V, W):
            r = rep_of[tuple(int(x) for x in v)]
            i = index[r]
            if cert[fi, i, 0] != 255:
                continue
            S = s34.full(v)
            imgs = np.zeros((5040, 128), dtype=np.int64)
            imgs[np.arange(5040)[:, None], P] = S[None, :]
            a = next(a for a in range(5040) if tuple(int(x) for x in imgs[a, 1:64]) == r)
            p = perms[a]
            A = [[0] * 7 for _ in range(7)]
            for (u, x), wt in zip(EDGES, w):
                A[p[u]][p[x]] = A[p[x]][p[u]] = int(wt)
            nw = [A[u][x] for (u, x) in EDGES]
            assert cut_vector(F, nw) == r, (q, i)
            cert[fi, i] = nw
        log(f'  GF({q}): {int((cert[fi, :, 0] != 255).sum())} orbits')
    ks = s34.known_states()
    st = tuple(int(x) for x in ks['F'][1:64]); nf = tuple(int(x) for x in ks['NF'][1:64])
    np.savez_compressed(out, reps=np.array(reps, dtype=np.uint8), sizes=np.array(sizes, dtype=np.int64),
                        fields=np.array(FIELDS, dtype=np.int64), cert=cert,
                        steane=np.array([rep_of[st] == r for r in reps]), nonfano=np.array([rep_of[nf] == r for r in reps]))
    log(f'  {len(reps)} orbits, {sum(sizes)} vectors -> {out}')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--check', action='store_true')
    ap.add_argument('--full', action='store_true')
    ap.add_argument('--enumerate', type=int, default=None, help='q: cc/s40_enum7q.c, orbits against the table')
    ap.add_argument('--all', action='store_true', help='with --enumerate: without the relabelling reduction')
    ap.add_argument('--slice', type=int, default=0, help='with --enumerate: run only slice s of --nslices')
    ap.add_argument('--nslices', type=int, default=1)
    ap.add_argument('--sample', type=int, default=None, help='q: cc/s39_sample7q.c')
    ap.add_argument('--samples', type=int, default=10 ** 6)
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--out', default=None)
    ap.add_argument('--build', default=None, help='prefix of the sampler outputs <prefix><q>.bin/.w')
    ap.add_argument('--prop16', action='store_true')
    a = ap.parse_args()
    if a.enumerate:
        Z = load()
        reps = [tuple(int(x) for x in r) for r in Z['reps']]
        idx, _ = orbit_map(reps)
        fi = [int(x) for x in Z['fields']].index(a.enumerate) if a.enumerate in FIELDS else None
        t0 = time.time()
        V, W, tot = run_enum(a.enumerate, all_=a.all, slice_=a.slice, nslices=a.nslices)
        orb = {idx[v] for v in V if v in idx}
        print(f'GF({a.enumerate}){"" if a.nslices == 1 else f" slice {a.slice} of {a.nslices}"}: {tot}; '
              f'{time.time() - t0:.0f} s; outside the table {sum(v not in idx for v in V)}; orbits {len(orb)}' +
              ('' if fi is None or a.nslices > 1 else
               f'; equal to the certified orbits {orb == set(np.flatnonzero(Z["cert"][fi, :, 0] != 255))}'))
        if a.out:
            np.savez_compressed(a.out, vectors=np.array(V, dtype=np.uint8), weights=np.array(W, dtype=np.uint8))
    if a.sample:
        V, W = run_sample(a.sample, a.samples, a.seed, a.out or f's_gf{a.sample}')
        print(f'GF({a.sample}): {a.samples} graphs, {len(V)} distinct vectors -> {a.out or f"s_gf{a.sample}"}.bin/.w')
    if a.build:
        build_table(a.build, a.out or DATA)
    if a.prop16:
        r = prop16_check()
        print('Proposition 16:', 'ALL OK' if r['prop16_ok'] else 'FAILURE')
    if a.check:
        o = check(full_check=a.full)
        print('ALL OK' if check_ok(o, a.full) else 'FAILURE')


if __name__ == '__main__':
    main()
