#!/usr/bin/env python3
"""s26: stabilizer entropy cones are self-dual slices of linear rank cones (R32).

Notation: n visible parties (bits 0..n-1) and a purifier (bit n); S is pure-symmetric, S(X) = S(E \\ X);
h_S(X) = S(X) + sum_{P in X} S(P) (the h_norm of docs/A2_stabilizer_vs_QLR.md); pi(r)(X) = r(X) + r(E \\ X)
- r(E) (the connectivity function of a polymatroid r).  Theorem 1 of docs/R32_reduction_and_local_dimension.md
states, for every prime p,

    Stab_n^(p) = { S : h_S in LR_{n+1}^(p) } = pi(LR_{n+1}^(p)) = closed cone of CSS-state entropy vectors,

with the proof
  (i)   Lemma 2 of the A2 note: h_S is the rank function of n+1 subspaces over GF(p) for every stabilizer state;
  (ii)  the CSS state sum_{c in rowspace M} |c> (columns of M grouped into parties) has S = pi(r_M);
  (iii) pi(h_S) = 2S for every pure-symmetric S;
  (iv)  so S = pi(h_S)/2 is a CSS vector whenever h_S is representable.

This script checks (ii) and a direct consequence of (i)-(iii) -- the CSS state of the stabilizer (Lagrangian)
code of a graph state has entropy exactly 2S -- against explicit state vectors, and runs the n = 4 consistency
check: on five elements (A..D and the purifier E) the 24 five-variable linear rank inequalities of
Dougherty-Freiling-Zeger (data/dfz_ref28.csv, arXiv:0910.0284 eqs. (1)-(24)) pulled back through h_S are implied
by the pulled-back Shannon and Ingleton (eqs. (36)-(39)) rows, and every extreme ray of the pulled-back
Shannon + Ingleton cone is realised by a qubit graph state.  Hence Stab_4 = h^-1(Shannon + Ingleton) on five
elements, in agreement with Linden-Matus-Ruskai-Winter (arXiv:1302.5453).

  python3 scripts/s26_reduction.py --self-test [--trials 24]
  python3 scripts/s26_reduction.py --n4                       # LP redundancy, certificates in data/ re-verified
  python3 scripts/s26_reduction.py --n4 --lrs                 # + exact extreme-ray enumeration (needs lrs)
  python3 scripts/s26_reduction.py --n4 --lrs --realise --out data/s26_stab4_certs.npz   # regenerate
"""
import argparse, itertools, os, shutil, subprocess, sys, tempfile, time
from fractions import Fraction
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..'))
from epr1kit import stabcert

DATA = os.path.join(HERE, '..', 'data')
CERTS = os.path.join(DATA, 's26_stab4_certs.npz')


# ------------------------------------------------------------------ CSS entropies
def group_rank(G, cols, p):
    return stabcert.gfp_rank([[int(G[i][j]) for j in cols] for i in range(len(G))], p) if cols else 0


def css_entropy(G, owner, n, p):
    """S over the 2^n - 1 masks of the visible parties for the CSS state of rowspace(G); owner[j] = party of
    column j (0..n, n = purifier).  S(X) = r(X) + r(E \\ X) - r(E)."""
    m = len(owner)
    rE = group_rank(G, list(range(m)), p)
    S = np.zeros((1 << n) - 1, dtype=np.int64)
    for X in range(1, 1 << n):
        inX = [j for j in range(m) if owner[j] < n and X >> owner[j] & 1]
        out = [j for j in range(m) if j not in inX]
        S[X - 1] = group_rank(G, inX, p) + group_rank(G, out, p) - rE
    return S


def css_state(G, p):
    G = np.asarray(G, dtype=np.int64) % p
    k, m = G.shape
    psi = np.zeros(p ** m)
    weights = p ** np.arange(m - 1, -1, -1)
    for coeffs in itertools.product(range(p), repeat=k):
        c = (np.array(coeffs, dtype=np.int64) @ G) % p
        psi[int(c @ weights)] += 1.0
    return psi / np.linalg.norm(psi)


def sv_entropy(psi, p, m, cols):
    rest = [j for j in range(m) if j not in cols]
    M = psi.reshape([p] * m).transpose(list(cols) + rest).reshape(p ** len(cols), p ** len(rest))
    s = np.linalg.svd(M, compute_uv=False) ** 2
    s = s[s > 1e-12]
    return float(-(s * np.log(s)).sum() / np.log(p))


def sv_entropy_vector(G, owner, n, p):
    psi = css_state(G, p)
    m = len(owner)
    return np.array([sv_entropy(psi, p, m, [j for j in range(m) if owner[j] < n and X >> owner[j] & 1])
                     for X in range(1, 1 << n)])


def graph_entropy(W, owner, n, p):
    N = len(owner)
    S = np.zeros((1 << n) - 1, dtype=np.int64)
    for X in range(1, 1 << n):
        inX = [i for i in range(N) if owner[i] < n and X >> owner[i] & 1]
        out = [i for i in range(N) if i not in inX]
        S[X - 1] = stabcert.gfp_rank([[int(W[i][j]) for j in out] for i in inX], p) if inX and out else 0
    return S


def self_test(trials=24, seed=26):
    """max deviation of (a) the CSS formula from state-vector entropies on random codes and (b) the
    state-vector entropy of the CSS state of a graph state's Lagrangian code from 2 * (graph-state entropy)."""
    rng = np.random.default_rng(seed)
    dev_a = dev_b = 0.0
    for t in range(trials):
        p = (2, 3)[t % 2]
        n = int(rng.integers(2, 4))
        m = int(rng.integers(3, 8 if p == 2 else 6))
        k = int(rng.integers(1, m))
        G = rng.integers(0, p, size=(k, m))
        owner = [int(x) for x in rng.integers(0, n + 1, size=m)]
        dev_a = max(dev_a, float(np.abs(sv_entropy_vector(G, owner, n, p) - css_entropy(G, owner, n, p)).max()))
        N = int(rng.integers(2, 5 if p == 2 else 4))           # graph state on N qudits
        W = np.triu(rng.integers(0, p, size=(N, N)), 1)
        W = (W + W.T) % p
        gown = [int(x) for x in rng.integers(0, n + 1, size=N)]
        S = graph_entropy(W, gown, n, p)
        L = np.hstack([np.eye(N, dtype=np.int64), W])            # Lagrangian [I | W], x then z coordinates
        lown = gown + gown
        dev_b = max(dev_b, float(np.abs(sv_entropy_vector(L, lown, n, p) - 2 * S).max()))
    return dev_a, dev_b


# ------------------------------------------------------------------ n = 4 consistency check
N4, FULL5 = 4, 31


def pull4(row):
    """row over 5-element masks (h coefficients; bit 4 = purifier) -> coefficients on the 15 S coordinates."""
    v = np.zeros(16)

    def addS(X, c):
        if X == 0 or X == FULL5:
            return
        if X >> N4 & 1:
            X = FULL5 ^ X
        v[X] += c
    for M, c in row.items():
        addS(M, c)
        for P in range(5):
            if M >> P & 1:
                addS(1 << P, c)
    return v[1:]


def shannon5():
    rows = [{FULL5: 1, FULL5 ^ (1 << i): -1} for i in range(5)]
    for i, j in itertools.combinations(range(5), 2):
        rest = [k for k in range(5) if k not in (i, j)]
        for kk in range(1 << 3):
            K = sum(1 << rest[t] for t in range(3) if kk >> t & 1)
            r = {}
            for M, c in ((K | 1 << i, 1), (K | 1 << j, 1), (K | 1 << i | 1 << j, -1), (K, -1)):
                if M:
                    r[M] = r.get(M, 0) + c
            rows.append(r)
    return rows


def distinct_rows(vecs):
    out = {}
    for v in vecs:
        if not np.any(v):
            continue
        a = v.astype(np.int64)
        g = int(np.gcd.reduce(np.abs(a[a != 0])))
        out[tuple(a // g)] = True
    return [np.array(k, dtype=float) for k in out]


def n4_rows():
    import s23_ci_dfz as d
    dfz = d.load_dfz()
    ing = [k for k in dfz if k.startswith('ing')]
    new = [k for k in dfz if not k.startswith('ing')]
    inst = lambda ids: [pull4(d.instance(dfz[i], phi, 0)) for i in ids for phi in itertools.permutations(range(5))]
    SH = distinct_rows([pull4(r) for r in shannon5()])
    ING = distinct_rows(inst(ing))
    NEW = distinct_rows(inst(new))
    return SH, ING, NEW


def lp_implied(A, v):
    from scipy.optimize import linprog
    sing = [1 if bin(X).count('1') == 1 else 0 for X in range(1, 16)]
    res = linprog(v, A_ub=-A, b_ub=np.zeros(len(A)), A_eq=np.array([sing]), b_eq=[1],
                  bounds=[(None, None)] * 15, method='highs')
    return res.status == 0 and res.fun >= -1e-9


def lrs_rays(A):
    """exact extreme rays of {S : A S >= 0} with lrs (None if lrs is not installed)."""
    exe = shutil.which('lrs')
    if exe is None:
        return None
    with tempfile.TemporaryDirectory() as td:
        f = os.path.join(td, 'c4.ine')
        with open(f, 'w') as fh:
            fh.write('c4\nH-representation\nbegin\n%d 16 rational\n' % len(A))
            for row in A:
                fh.write('0 ' + ' '.join(str(int(x)) for x in row) + '\n')
            fh.write('end\n')
        out = subprocess.run([exe, f], capture_output=True, text=True).stdout
    rays = []
    for ln in out.splitlines():
        t = ln.split()
        if len(t) == 16 and t[0] == '0':
            fr = [Fraction(x) for x in t[1:]]
            den = np.lcm.reduce([x.denominator for x in fr])
            iv = np.array([int(x * den) for x in fr], dtype=np.int64)
            g = int(np.gcd.reduce(np.abs(iv[iv != 0])))
            rays.append(iv // g)
    return np.array(sorted(map(tuple, rays)), dtype=np.int64)


def is_extreme(A, r):
    tight = A[np.abs(A @ r) < 1e-9]
    return tight.shape[0] > 0 and np.linalg.matrix_rank(tight) == 14


def purifier_entropy(r):
    return int(r[14])                         # S(E) = S(ABCD), mask 15 -> index 14


def realise(r, rng, lams=(1, 2, 3), restarts=40, iters=2000):
    """random local search over qubit graph states with lam * S(P) qubits per party; returns a certificate."""
    for lam in lams:
        sizes = [lam * int(r[(1 << P) - 1]) for P in range(4)] + [lam * purifier_entropy(r)]
        owner = [P for P in range(5) for _ in range(sizes[P])]
        N = len(owner)
        target = lam * np.asarray(r, dtype=np.int64)
        for _ in range(restarts):
            W = np.triu(rng.integers(0, 2, size=(N, N)), 1); W = W + W.T
            cost = int(np.abs(graph_entropy(W, owner, 4, 2) - target).sum())
            for _it in range(iters):
                if cost == 0:
                    break
                i, j = rng.choice(N, 2, replace=False)
                if owner[i] == owner[j]:
                    continue
                W[i, j] ^= 1; W[j, i] ^= 1
                c2 = int(np.abs(graph_entropy(W, owner, 4, 2) - target).sum())
                if c2 <= cost or rng.random() < 0.02:
                    cost = c2
                else:
                    W[i, j] ^= 1; W[j, i] ^= 1
            if cost == 0:
                return dict(lam=lam, sizes=sizes, W=W.astype(np.uint8))
    return None


def check_certs(path=CERTS, A=None):
    """re-verify data/s26_stab4_certs.npz: every ray valid, extreme, distinct, realised at its lambda."""
    Z = np.load(path)
    rays, lam, sizes, Ws, Ns = Z['rays'].astype(np.int64), Z['lam'], Z['sizes'], Z['W'], Z['N']
    if A is None:
        SH, ING, _ = n4_rows()
        A = np.array(SH + ING)
    ok = len({tuple(x) for x in rays}) == rays.shape[0]
    for k in range(rays.shape[0]):
        r = rays[k]; n_ = int(Ns[k])
        ok = (ok and float((A @ r).min()) >= 0 and is_extreme(A, r)
              and stabcert.verify_realisation_n(r, 4, int(lam[k]), sizes[k], Ws[k][:n_, :n_], 2))
    return ok, rays


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--self-test', action='store_true')
    ap.add_argument('--trials', type=int, default=24)
    ap.add_argument('--n4', action='store_true')
    ap.add_argument('--lrs', action='store_true', help='enumerate the extreme rays exactly with lrs')
    ap.add_argument('--realise', action='store_true', help='search qubit graph states for the lrs rays')
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    if a.self_test:
        t0 = time.time()
        da, db = self_test(a.trials)
        print(f'CSS formula vs state vector: max deviation {da:.1e}; '
              f'CSS state of the Lagrangian code vs 2S: max deviation {db:.1e}  ({a.trials} trials, {time.time()-t0:.0f}s)')
    if a.n4:
        t0 = time.time()
        SH, ING, NEW = n4_rows()
        A = np.array(SH + ING)
        bad = sum(1 for v in NEW if not lp_implied(A, v))
        print(f'n = 4: pulled back to the 15 coordinates: Shannon {len(SH)}, Ingleton forms {len(ING)}, '
              f'DFZ (1)-(24) {len(NEW)} distinct rows; DFZ rows not implied by Shannon + Ingleton: {bad}  ({time.time()-t0:.0f}s)')
        rays = None
        if a.lrs:
            rays = lrs_rays(A)
            if rays is None:
                print('lrs not found: skipping the enumeration')
            else:
                print(f'extreme rays of the pulled-back Shannon + Ingleton cone (lrs, exact): {rays.shape[0]}')
        if a.realise and rays is not None:
            rng = np.random.default_rng(4)
            certs, t1 = [], time.time()
            for r in rays:
                c = realise(r, rng)
                if c is None:
                    print('not realised:', r.tolist()); continue
                certs.append((r, c))
            Nmax = max(int(sum(c['sizes'])) for _, c in certs)
            Wst = np.zeros((len(certs), Nmax, Nmax), dtype=np.uint8)
            for k, (_, c) in enumerate(certs):
                n_ = int(sum(c['sizes'])); Wst[k, :n_, :n_] = c['W']
            out = a.out or CERTS
            np.savez_compressed(out, rays=np.array([r for r, _ in certs], dtype=np.int16),
                                lam=np.array([c['lam'] for _, c in certs], dtype=np.int8),
                                sizes=np.array([c['sizes'] for _, c in certs], dtype=np.int8),
                                N=np.array([sum(c['sizes']) for _, c in certs], dtype=np.int16), W=Wst)
            print(f'realised {len(certs)}/{rays.shape[0]} rays by qubit graph states '
                  f'(lambda counts {np.bincount([c["lam"] for _, c in certs]).tolist()}) -> {out}  ({time.time()-t1:.0f}s)')
        if os.path.exists(a.out or CERTS):
            ok, R = check_certs(a.out or CERTS, A)
            same = rays is None or (R.shape == rays.shape and {tuple(x) for x in R} == {tuple(x) for x in rays})
            print(f'certificates {a.out or CERTS}: {R.shape[0]} rays, all valid, extreme and realised: {ok}'
                  + ('' if rays is None else f'; equal to the lrs ray set: {same}'))


if __name__ == '__main__':
    main()
