#!/usr/bin/env python3
"""s24: graph-state search over GF(q), q = p^k, with field reduction to GF(p) certificates (R31).

A weighted graph state over GF(q) (symmetric N x N matrix W, zero diagonal) is a stabilizer state of N
Galois qudits with S(X) = rank_GF(q) W[X, X^c] in units of log q: its stabilizer is the GF(q)-linear
Lagrangian rowspace [I | W].  Field reduction turns it into an ordinary stabilizer state of k*N qudits of
prime dimension p with k times the entropy:

  * write x-coordinates in the polynomial basis e_j = a^j of GF(q) over GF(p) and z-coordinates in the
    trace-dual basis d_j (Tr(e_i d_j) = delta_ij, Tr = GF(q) -> GF(p) trace); then
    Tr(x z' - z x') is the standard GF(p) symplectic form on the k qudits that replace one Galois qudit;
  * the GF(p)-span of the rows of [I | W] and their multiples by a^1..a^(k-1) is a GF(p) Lagrangian of
    dimension k*N; its projection onto any set of Galois qudits has k times the GF(q) dimension.

So a hit for r at lambda over GF(q) is a GF(p) realisation of k*lambda*r.  The script finds hits with
cc/s24_anneal_gfq.c, re-checks them over GF(q) with its own arithmetic, performs the reduction, brings the
GF(p) state to graph form by local Fourier transforms (x, z) -> (z, -x) on suitable qudits, and accepts it
only if the independent checker of epr1kit.stabcert verifies it (verify_realisation for p = 2,
verify_realisation_gfp for odd p) at k*lambda.

  python3 scripts/s24_galois_search.py rays.npy --q 4 --lams 1 --seconds 60 --out certs.npz
  python3 scripts/s24_galois_search.py data/s21_qudit_realisations.npz --q 4 --seconds 300   # uses the 'rep' field
  python3 scripts/s24_galois_search.py --self-test       # reduction vs. direct ranks on random states

Motivation (R31): the 16 orbits realised only by qutrits in R29 are realised by GF(4) graph states at
lambda = 1, hence by qubit graph states at lambda = 2 (data/s24_qubit_l2_certs.npz).
"""
import argparse, json, os, shutil, subprocess, sys, tempfile
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..'))
from s16_stab_search import party_sizes, omp_threads
from epr1kit import stabcert

MODULUS = {4: (2, 2, [1, 1]), 8: (2, 3, [1, 1, 0]), 9: (3, 2, [2, 0])}   # a^k = sum_i c_i a^i


class GF:
    """GF(q) for q in {2, 3, 4, 5, 7, 8, 9}; element sum_i c_i a^i encoded as sum_i c_i p^i (as in the kernel)."""

    def __init__(self, q):
        if q in (2, 3, 5, 7):
            self.p, self.k, self.mod = q, 1, []
        elif q in MODULUS:
            self.p, self.k, self.mod = MODULUS[q]
        else:
            raise ValueError(f'unsupported q = {q}')
        self.q = q
        self.mul_t = [[self._mul(a, b) for b in range(q)] for a in range(q)]
        self.inv_t = {a: b for a in range(1, q) for b in range(1, q) if self.mul_t[a][b] == 1}
        assert len(self.inv_t) == q - 1, 'not a field'

    def digits(self, a):
        return [(a // self.p ** i) % self.p for i in range(self.k)]

    def enc(self, d):
        return sum((c % self.p) * self.p ** i for i, c in enumerate(d))

    def add(self, a, b):
        return self.enc([x + y for x, y in zip(self.digits(a), self.digits(b))])

    def sub(self, a, b):
        return self.enc([x - y for x, y in zip(self.digits(a), self.digits(b))])

    def _mul(self, a, b):
        da, db = self.digits(a), self.digits(b)
        pr = [0] * (2 * self.k)
        for i in range(self.k):
            for j in range(self.k):
                pr[i + j] = (pr[i + j] + da[i] * db[j]) % self.p
        for deg in range(2 * self.k - 2, self.k - 1, -1):
            c, pr[deg] = pr[deg], 0
            for i in range(self.k):
                pr[deg - self.k + i] = (pr[deg - self.k + i] + c * self.mod[i]) % self.p
        return self.enc(pr[:self.k])

    def mul(self, a, b):
        return self.mul_t[a][b]

    def trace(self, a):
        """Tr(a) = a + a^p + ... + a^(p^(k-1)), an element of the prime field (returned as int mod p)."""
        t, y = 0, a
        for _ in range(self.k):
            t = self.add(t, y)
            z = 1
            for _ in range(self.p):
                z = self.mul(z, y)
            y = z
        assert t < self.p, 'trace not in the prime field'
        return t

    def rank(self, M):
        M = [list(r) for r in M]
        r = 0
        ncol = len(M[0]) if M else 0
        for c in range(ncol):
            piv = next((i for i in range(r, len(M)) if M[i][c]), None)
            if piv is None:
                continue
            M[r], M[piv] = M[piv], M[r]
            iv = self.inv_t[M[r][c]]
            M[r] = [self.mul(x, iv) for x in M[r]]
            for i in range(len(M)):
                if i != r and M[i][c]:
                    f = M[i][c]
                    M[i] = [self.sub(x, self.mul(f, y)) for x, y in zip(M[i], M[r])]
            r += 1
        return r


def rank_mod_p(M, p):
    M = (np.array(M, dtype=np.int64) % p).copy()
    r = 0
    for c in range(M.shape[1] if M.size else 0):
        piv = [i for i in range(r, M.shape[0]) if M[i, c]]
        if not piv:
            continue
        M[[r, piv[0]]] = M[[piv[0], r]]
        M[r] = (M[r] * pow(int(M[r, c]), p - 2, p)) % p
        for i in range(M.shape[0]):
            if i != r and M[i, c]:
                M[i] = (M[i] - M[i, c] * M[r]) % p
        r += 1
    return r


def reduce_to_prime(W, F):
    """GF(q) graph state (N x N symmetric W) -> GF(p) Lagrangian generator matrix (kN x 2kN): x-part, z-part."""
    N, k, p = len(W), F.k, F.p
    e = [p ** j for j in range(k)]                       # polynomial basis a^j
    rows = []
    for r_ in range(N):
        for s in e:                                      # the row times a^s
            x = [s if i == r_ else 0 for i in range(N)]
            z = [F.mul(s, W[r_][i]) for i in range(N)]
            xb = [F.digits(v)[j] for v in x for j in range(k)]             # coordinates in basis e
            zb = [F.trace(F.mul(v, e[j])) for v in z for j in range(k)]    # coordinates in the dual basis
            rows.append(xb + zb)
    return np.array(rows, dtype=np.int64)


def symplectic_ok(L, p):
    n = L.shape[1] // 2
    X, Z = L[:, :n], L[:, n:]
    return not ((X @ Z.T - Z @ X.T) % p).any()


def to_graph(L, p):
    """local Fourier transforms (x, z) -> (z, -x) on the qudits outside the pivot set of the X block, then
    Gamma = X^{-1} Z (mod p); symmetric for a Lagrangian; diagonal cleared (local phases)."""
    n = L.shape[1] // 2
    M = L.copy() % p
    r, piv = 0, []
    for c in range(n):
        rows = [i for i in range(r, M.shape[0]) if M[i, c]]
        if not rows:
            continue
        M[[r, rows[0]]] = M[[rows[0], r]]
        M[r] = (M[r] * pow(int(M[r, c]), p - 2, p)) % p
        for i in range(M.shape[0]):
            if i != r and M[i, c]:
                M[i] = (M[i] - M[i, c] * M[r]) % p
        piv.append(c); r += 1
    for c in (c for c in range(n) if c not in piv):
        x, z = M[:, c].copy(), M[:, n + c].copy()
        M[:, c], M[:, n + c] = z, (-x) % p
    X, Z = M[:, :n], M[:, n:]
    A = np.concatenate([X % p, np.eye(n, dtype=np.int64)], axis=1)
    rr = 0
    for c in range(n):
        rows = [i for i in range(rr, n) if A[i, c]]
        if not rows:
            raise ValueError('X block not invertible')
        A[[rr, rows[0]]] = A[[rows[0], rr]]
        A[rr] = (A[rr] * pow(int(A[rr, c]), p - 2, p)) % p
        for i in range(n):
            if i != rr and A[i, c]:
                A[i] = (A[i] - A[i, c] * A[rr]) % p
        rr += 1
    G = (A[:, n:] @ Z) % p
    if not (G == G.T).all():
        raise ValueError('not symmetric: input was not Lagrangian')
    np.fill_diagonal(G, 0)
    return G


def certify(r, lam, W, q):
    """returns (ok_q, cert) where cert = dict(p, lam, sizes, G) is verified by epr1kit.stabcert, or None."""
    F = GF(q)
    r = np.asarray(r, dtype=np.int64)
    sz = party_sizes(r, lam)
    N = sum(sz)
    party = [t for t in range(6) for _ in range(sz[t])]
    ok_q = len(W) == N and all(W[i][j] == W[j][i] for i in range(N) for j in range(N))
    for m in range(1, 32):
        X = [i for i in range(N) if party[i] < 5 and (m >> party[i]) & 1]
        Y = [i for i in range(N) if i not in X]
        if ok_q and F.rank([[W[i][j] for j in Y] for i in X]) != lam * int(r[m - 1]):
            ok_q = False
    if not ok_q:
        return False, None
    L = reduce_to_prime(W, F)
    if not symplectic_ok(L, F.p) or rank_mod_p(L, F.p) != F.k * N:
        return True, None
    G = to_graph(L, F.p)
    lam_p, sz_p = F.k * lam, party_sizes(r, F.k * lam)
    if F.p == 2:
        rows = np.zeros(64, dtype=np.uint64)
        for a in range(G.shape[0]):
            rows[a] = np.uint64(sum(1 << b for b in range(G.shape[0]) if G[a, b]))
        ok = stabcert.verify_realisation(r, lam_p, sz_p, rows)
    else:
        ok = stabcert.verify_realisation_gfp(r, lam_p, sz_p, G, F.p)
    return True, (dict(p=F.p, lam=lam_p, sizes=sz_p, G=G.astype(np.uint8)) if ok else None)


def self_test(trials=12, seed=24):
    """field reduction vs. direct ranks: for random GF(q) graph states the reduced GF(p) graph state has
    exactly k times the GF(q) cut ranks, for every cut between party blocks."""
    rng = np.random.default_rng(seed)
    worst = 0
    for q in (4, 8, 9):
        F = GF(q)
        # the trace form is non-degenerate: the Gram matrix Tr(e_i e_j) of the polynomial basis is invertible mod p
        e = [F.p ** j for j in range(F.k)]
        if rank_mod_p([[F.trace(F.mul(x, y)) for y in e] for x in e], F.p) != F.k:
            return 10 ** 6
        for _ in range(trials):
            N = int(rng.integers(3, 7))
            W = np.zeros((N, N), dtype=np.int64)
            for i in range(N):
                for j in range(i + 1, N):
                    W[i, j] = W[j, i] = int(rng.integers(0, q))
            L = reduce_to_prime(W.tolist(), F)
            assert symplectic_ok(L, F.p) and rank_mod_p(L, F.p) == F.k * N
            G = to_graph(L, F.p)
            for mask in range(1, 2 ** N - 1):
                X = [i for i in range(N) if mask >> i & 1]
                Y = [i for i in range(N) if not mask >> i & 1]
                Xr = [F.k * i + j for i in X for j in range(F.k)]
                Yr = [F.k * i + j for i in Y for j in range(F.k)]
                d = rank_mod_p(G[np.ix_(Xr, Yr)], F.p) - F.k * F.rank([[int(W[i, j]) for j in Y] for i in X])
                worst = max(worst, abs(d))
    return worst


def build_kernel(tmp):
    src = os.path.join(HERE, '..', 'cc', 's24_anneal_gfq.c')
    exe = os.path.join(tmp, 's24_anneal_gfq')
    flags = ['-O3', '-march=native', '-fopenmp']
    p = subprocess.run(['gcc'] + flags + ['-o', exe, src, '-lm'], capture_output=True, text=True)
    if p.returncode != 0:
        raise SystemExit('kernel build failed: ' + p.stderr[:400])
    return exe, ' '.join(flags)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('rays', nargs='?')
    ap.add_argument('--q', type=int, default=4, choices=[2, 3, 4, 5, 7, 8, 9])
    ap.add_argument('--lams', type=int, nargs='+', default=[1])
    ap.add_argument('--seconds', type=float, default=60.0)
    ap.add_argument('--seed', type=int, default=20261002)
    ap.add_argument('--threads', type=int, default=None)
    ap.add_argument('--index', type=int, nargs='*', default=None)
    ap.add_argument('--out', default=None)
    ap.add_argument('--self-test', action='store_true')
    ap.add_argument('--exhaustive', action='store_true',
                    help='rays with S_p = 1 for all six parties: enumerate all q^15 one-qudit-per-party graph states '
                         '(cc/s24_exhaustive6.c, q prime) instead of annealing; lambda = 1 only')
    a = ap.parse_args()
    if a.self_test:
        w = self_test()
        print(f'field reduction self-test (q = 4, 8, 9; random graph states, every cut): max |rank_p - k*rank_q| = {w}')
        sys.exit(0 if w == 0 else 1)
    R = np.load(a.rays)
    R = (R['rep'] if a.rays.endswith('.npz') else R).astype(np.int64)     # .npz: its 'rep' field
    idx = list(range(R.shape[0])) if a.index is None else a.index
    if a.exhaustive:
        if a.q not in (2, 3, 5, 7):
            raise SystemExit('--exhaustive needs a prime q')
        tmp = tempfile.mkdtemp(prefix='s24x_')
        exe = os.path.join(tmp, 's24_exhaustive6')
        p = subprocess.run(['gcc', '-O3', '-march=native', '-o', exe, os.path.join(HERE, '..', 'cc', 's24_exhaustive6.c')],
                           capture_output=True, text=True)
        if p.returncode != 0:
            raise SystemExit('build failed: ' + p.stderr[:400])
        for i in idx:
            if list(party_sizes(R[i], 1)) != [1] * 6:
                print(f'ray {i}: party sizes {party_sizes(R[i], 1)} at lambda 1 -- not one qudit per party, skipped')
                continue
            out = subprocess.run([exe, str(a.q)], input=' '.join(map(str, R[i].tolist())), capture_output=True, text=True).stdout
            print(f'ray {i}: ' + out.strip().splitlines()[-1], flush=True)
        shutil.rmtree(tmp, ignore_errors=True)
        return
    tmp = tempfile.mkdtemp(prefix='s24_')
    exe, flags = build_kernel(tmp)
    nth = omp_threads(a.threads)
    env = dict(os.environ); env['OMP_NUM_THREADS'] = str(nth)
    print(f'kernel: gcc {flags}; OpenMP threads {nth}', flush=True)
    certs, summary = [], {}
    for lam in a.lams:
        lines = [' '.join(map(str, [i, lam] + list(party_sizes(R[i], lam)) + [lam * int(x) for x in R[i]])) for i in idx]
        out = subprocess.run([exe, str(a.q), str(a.seconds), str(a.seed + lam)], input='\n'.join(lines) + '\n',
                             capture_output=True, text=True, env=env).stdout
        for ln in out.splitlines():
            t = ln.split()
            i = int(t[0])
            if t[2] != 'FOUND':
                summary[f'{i}:{lam}'] = int(t[3])
                continue
            N = int(t[3]); W = [[int(ch) for ch in row] for row in t[4:4 + N]]
            ok_q, cert = certify(R[i], lam, W, a.q)
            if cert is None:
                print(f'WARNING: ray {i} lambda {lam}: hit over GF({a.q}) {"verified" if ok_q else "NOT verified"} '
                      f'but no verified GF(p) certificate; discarded', flush=True)
                continue
            certs.append(dict(index=i, rep=R[i], lam_q=lam, **cert, Wq=np.array(W, dtype=np.uint8)))
            summary[f'{i}:{lam}'] = 0
        n_hit = sum(1 for c in certs if c['lam_q'] == lam)
        print(f'q={a.q} lambda={lam}: {len(idx)} rays tried, {n_hit} certified over GF(p) at {GF(a.q).k * lam} x lambda', flush=True)
    shutil.rmtree(tmp, ignore_errors=True)
    if a.out:
        pad = lambda M, n: np.pad(M, ((0, n - M.shape[0]), (0, n - M.shape[1])))
        np.savez_compressed(a.out, index=np.array([c['index'] for c in certs], dtype=np.int64),
                            rep=np.array([c['rep'] for c in certs], dtype=np.int16).reshape(-1, 31),
                            q=np.full(len(certs), a.q, dtype=np.int8),
                            lam_q=np.array([c['lam_q'] for c in certs], dtype=np.int8),
                            p=np.array([c['p'] for c in certs], dtype=np.int8),
                            lam=np.array([c['lam'] for c in certs], dtype=np.int8),
                            sizes=np.array([c['sizes'] for c in certs], dtype=np.int16).reshape(-1, 6),
                            G=np.array([pad(c['G'], 64) for c in certs], dtype=np.uint8).reshape(-1, 64, 64),
                            Wq=np.array([pad(c['Wq'], 64) for c in certs], dtype=np.uint8).reshape(-1, 64, 64))
        json.dump({'q': a.q, 'lams': a.lams, 'seconds': a.seconds, 'found': len(certs), 'best_residual': summary},
                  open(a.out + '.json', 'w'), indent=1)
        print(f'certificates -> {a.out}')


if __name__ == '__main__':
    main()
