#!/usr/bin/env python3
"""s21: qudit (prime p = 3, 5, 7) graph-state search for five-party rays (R29).

Same normal form as s16 (party q owns lambda * S_q qudits; only inter-party weights
matter), but over GF(p): a weighted qudit graph state with symmetric weight matrix W
(zero diagonal, entries 0..p-1) has S(X) = rank_GF(p) W[X, X^c] in units of log p.
The C kernel cc/s21_anneal_gfp.c anneals over the inter-party weights; every hit is
re-verified here with an independent GF(p) elimination, and `--self-test` checks the
rank formula itself against exact state vectors (random graph states on 6 qutrits,
5 ququints or 4 qudits of dimension 7, every bipartition, singular values of the
reshaped amplitude matrix).

  python3 scripts/s21_qudit_search.py rays.npy --p 3 --lams 1 --seconds 30 --out certs.npz
  python3 scripts/s21_qudit_search.py --self-test

A hit is an exact certificate that lambda * r is the entropy vector (base p) of an
explicit qudit stabilizer state, hence r lies in the stabilizer cone (all primes).
A miss proves nothing.
"""
import argparse, json, os, shutil, subprocess, sys, tempfile, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, '..')
sys.path.insert(0, HERE)
from s16_stab_search import party_sizes, omp_threads

MAXN = 64


def rank_gfp(M, p):
    """independent GF(p) rank (row reduction on an int matrix), p prime."""
    M = (np.array(M, dtype=np.int64) % p).copy()
    rows, cols = M.shape
    r = 0
    for c in range(cols):
        piv = next((i for i in range(r, rows) if M[i, c]), None)
        if piv is None:
            continue
        M[[r, piv]] = M[[piv, r]]
        M[r] = (M[r] * pow(int(M[r, c]), p - 2, p)) % p
        for i in range(rows):
            if i != r and M[i, c]:
                M[i] = (M[i] - M[i, c] * M[r]) % p
        r += 1
        if r == rows:
            break
    return r


def cut_vector_p(W, sizes, p):
    party = np.repeat(np.arange(6), sizes)
    out = np.zeros(31, dtype=np.int64)
    for m in range(1, 32):
        inX = np.array([(q < 5) and ((m >> q) & 1) == 1 for q in party])
        out[m - 1] = rank_gfp(W[np.ix_(inX, ~inX)], p) if inX.any() and (~inX).any() else 0
    return out


def verify_p(r, lam, sizes, W, p):
    W = np.asarray(W, dtype=np.int64)
    if W.shape[0] != sum(sizes) or not np.array_equal(W, W.T) or W.diagonal().any() or ((W < 0) | (W >= p)).any():
        return False
    return bool(np.array_equal(cut_vector_p(W, sizes, p), lam * np.asarray(r, dtype=np.int64)))


def self_test(p, trials=30, seed=1):
    """S(A) of a qudit graph state == rank_GF(p) W[A, B] (log base p), from exact state vectors."""
    n = {3: 6, 5: 5, 7: 4}[p]
    rng = np.random.default_rng(seed)
    omega = np.exp(2j * np.pi / p)
    xs = np.array(np.unravel_index(np.arange(p ** n), (p,) * n)).T          # all basis labels
    worst = 0.0
    for _ in range(trials):
        W = np.triu(rng.integers(0, p, (n, n)), 1); W = W + W.T
        phase = np.einsum('ki,ij,kj->k', xs, np.triu(W, 1), xs) % p
        psi = omega ** phase / p ** (n / 2)
        for m in range(1, 2 ** n - 1):
            A = [i for i in range(n) if m >> i & 1]; B = [i for i in range(n) if not m >> i & 1]
            T = np.transpose(psi.reshape((p,) * n), A + B).reshape(p ** len(A), p ** len(B))
            s = np.linalg.svd(T, compute_uv=False) ** 2
            s = s[s > 1e-12]
            S = float(-(s * np.log(s)).sum() / np.log(p))
            worst = max(worst, abs(S - rank_gfp(W[np.ix_(A, B)], p)))
    return worst, n


def build_kernel(outdir):
    exe = os.path.join(outdir, 's21_anneal_gfp')
    src = os.path.join(ROOT, 'cc', 's21_anneal_gfp.c')
    if shutil.which('gcc') is None:
        raise SystemExit('gcc not found -- install build tools first')
    for flags in (['-O3', '-march=native', '-fopenmp'], ['-O3', '-fopenmp'], ['-O3']):
        if subprocess.run(['gcc', *flags, '-o', exe, src, '-lm'], capture_output=True).returncode == 0:
            return exe, flags
    raise SystemExit('could not compile cc/s21_anneal_gfp.c')


def parse_rows(tokens, N):
    return np.array([[int(ch) for ch in row] for row in tokens[:N]], dtype=np.int64)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('rays', nargs='?')
    ap.add_argument('--p', type=int, default=3, choices=[3, 5, 7])
    ap.add_argument('--lams', type=int, nargs='+', default=[1])
    ap.add_argument('--seconds', type=float, default=30.0)
    ap.add_argument('--seed', type=int, default=20261001)
    ap.add_argument('--threads', type=int, default=None)
    ap.add_argument('--max-n', type=int, default=MAXN)
    ap.add_argument('--limit', type=int, default=None)
    ap.add_argument('--out', default=None)
    ap.add_argument('--self-test', action='store_true')
    a = ap.parse_args()
    if a.self_test:
        bad = 0
        for p in (3, 5, 7):
            w, n = self_test(p)
            print(f'self-test p={p}: max |S_svd - rank_GF({p})| over 30 random {n}-qudit graph states x {2 ** n - 2} cuts = {w:.2e}')
            bad += w > 1e-9
        raise SystemExit(1 if bad else 0)
    if not a.rays:
        ap.error('rays file required (or --self-test)')
    R = np.load(a.rays).astype(np.int64)
    if a.limit:
        R = R[:a.limit]
    tmp = tempfile.mkdtemp(prefix='s21_')
    exe, flags = build_kernel(tmp)
    env = dict(os.environ); nthr = omp_threads(a.threads); env['OMP_NUM_THREADS'] = str(nthr)
    print(f"kernel: gcc {' '.join(flags)}; OpenMP threads {nthr if '-fopenmp' in flags else 1}", flush=True)
    found, best = {}, {}
    t0 = time.time()
    for lam in a.lams:
        todo = [i for i in range(R.shape[0]) if i not in found and lam * (int(R[i, [0, 1, 3, 7, 15]].sum()) + int(R[i, 30])) <= a.max_n]
        if not todo:
            continue
        lines = [' '.join(map(str, [i, lam] + party_sizes(R[i], lam) + [lam * int(x) for x in R[i]])) for i in todo]
        pr = subprocess.run([exe, str(a.p), str(a.seconds), str(a.seed + lam), '1.5', '0.08', '3000'],
                            input='\n'.join(lines) + '\n', capture_output=True, text=True, env=env)
        if pr.returncode != 0:
            raise SystemExit(f'kernel failed: {pr.stderr[:500]}')
        nf = 0
        for ln in pr.stdout.splitlines():
            t = ln.split(); i, l = int(t[0]), int(t[1])
            if t[2] == 'FOUND':
                N = int(t[3]); W = parse_rows(t[4:], N); sz = party_sizes(R[i], l)
                if verify_p(R[i], l, sz, W, a.p):
                    found[i] = (l, sz, W); nf += 1
                else:
                    print(f'WARNING: kernel hit for ray {i} (lambda {l}) failed independent verification -- discarded')
            elif t[2] == 'BEST':
                best[(i, l)] = int(t[3])
        print(f'p={a.p} lambda={lam}: {len(todo)} rays tried, {nf} realised by qudit graph states and independently verified '
              f'(cumulative {len(found)}/{R.shape[0]}; {time.time() - t0:.0f}s)', flush=True)
    shutil.rmtree(tmp, ignore_errors=True)
    if a.out:
        idx = sorted(found)
        Wst = np.zeros((len(idx), MAXN, MAXN), dtype=np.uint8)
        for k, i in enumerate(idx):
            W = found[i][2]; Wst[k, :W.shape[0], :W.shape[0]] = W
        np.savez(a.out, index=np.array(idx, dtype=np.int64),
                 rep=R[idx].astype(np.int16) if idx else np.zeros((0, 31), np.int16),
                 lam=np.array([found[i][0] for i in idx], dtype=np.int8),
                 sizes=np.array([found[i][1] for i in idx], dtype=np.int16).reshape(-1, 6), W=Wst,
                 p=np.full(len(idx), a.p, dtype=np.int8))
        json.dump({'rays': int(R.shape[0]), 'found': len(found), 'p': a.p,
                   'by_lambda': {str(l): sum(1 for i in found if found[i][0] == l) for l in a.lams},
                   'best_residual': {f'{i}:{l}': e for (i, l), e in best.items()}}, open(a.out + '.json', 'w'), indent=1)
        print(f'certificates -> {a.out}')


if __name__ == '__main__':
    main()
