#!/usr/bin/env python3
"""s16: targeted stabilizer search (A2) -- realise extreme rays of QLR5 by qubit graph states.

For each target ray r (31 integers, mask order) and multiplier lambda, the C kernel
cc/s16_anneal.c anneals over graph states on N = lambda * sum_p S_p qubits, where
party p (A..E and the purifier F) owns exactly lambda * S_p qubits.  For a fixed
lambda this search space loses no generality (every stabilizer state is LC-equivalent
to a graph state; each party's marginal is maximally mixed, so surplus qubits factor
out; within-party edges never enter a cut).  Every hit is re-verified here with an
independent GF(2) implementation before it is accepted as a certificate.

  python3 scripts/s16_stab_search.py rays.npy --lams 1 2 --seconds 20 --out certs.npz

A hit is an S7-positive verdict: lambda * r is the entropy vector of an explicit
stabilizer state.  A miss is NOT a negative verdict (annealing is incomplete); the
exhaustive statements for tiny N are made separately (see README).
"""
import argparse, json, os, shutil, subprocess, sys, tempfile, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, '..')
sys.path.insert(0, ROOT)

SING = [0, 1, 3, 7, 15]          # indices of S_A..S_E in mask order


def party_sizes(r, lam):
    return [lam * int(r[i]) for i in SING] + [lam * int(r[30])]


def gf2_rank(M):
    """independent GF(2) rank (row reduction on a uint8 matrix)."""
    M = (np.array(M, dtype=np.uint8) & 1).copy()
    rows, cols = M.shape
    r = 0
    for c in range(cols):
        piv = None
        for i in range(r, rows):
            if M[i, c]:
                piv = i; break
        if piv is None:
            continue
        if piv != r:
            M[[r, piv]] = M[[piv, r]]
        for i in range(rows):
            if i != r and M[i, c]:
                M[i] ^= M[r]
        r += 1
        if r == rows:
            break
    return r


def cut_vector(A, sizes):
    """31 cut ranks of the graph state with adjacency A (N x N, GF(2)) and party sizes."""
    party = np.repeat(np.arange(6), sizes)
    out = np.zeros(31, dtype=np.int64)
    for m in range(1, 32):
        inX = np.array([(p < 5) and ((m >> p) & 1) == 1 for p in party])
        out[m - 1] = gf2_rank(A[np.ix_(inX, ~inX)]) if inX.any() and (~inX).any() else 0
    return out


def rows_to_matrix(hexrows, N):
    A = np.zeros((N, N), dtype=np.uint8)
    for i, h in enumerate(hexrows):
        v = int(h, 16)
        for j in range(N):
            A[i, j] = (v >> j) & 1
    return A


def verify(r, lam, sizes, A):
    """independent certificate check: symmetric, zero diagonal, cut ranks == lam * r."""
    if not np.array_equal(A, A.T) or A.diagonal().any():
        return False
    return bool(np.array_equal(cut_vector(A, sizes), lam * np.asarray(r, dtype=np.int64)))


def build_kernel(outdir):
    exe = os.path.join(outdir, 's16_anneal')
    src = os.path.join(ROOT, 'cc', 's16_anneal.c')
    if shutil.which('gcc') is None:
        raise SystemExit('gcc not found -- install build tools first')
    for flags in (['-O3', '-march=native', '-fopenmp'], ['-O3', '-fopenmp'], ['-O3']):
        if subprocess.run(['gcc', *flags, '-o', exe, src, '-lm'], capture_output=True).returncode == 0:
            return exe
    raise SystemExit('could not compile cc/s16_anneal.c')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('rays')
    ap.add_argument('--lams', type=int, nargs='+', default=[1, 2])
    ap.add_argument('--seconds', type=float, default=20.0, help='annealing budget per (ray, lambda)')
    ap.add_argument('--seed', type=int, default=20260930)
    ap.add_argument('--threads', type=int, default=None)
    ap.add_argument('--max-n', type=int, default=64)
    ap.add_argument('--sweeps', type=float, default=3000.0)
    ap.add_argument('--out', default=None)
    ap.add_argument('--limit', type=int, default=None)
    a = ap.parse_args()

    R = np.load(a.rays).astype(np.int64)
    if a.limit:
        R = R[:a.limit]
    tmp = tempfile.mkdtemp(prefix='s16_')
    exe = build_kernel(tmp)
    env = dict(os.environ)
    if a.threads:
        env['OMP_NUM_THREADS'] = str(a.threads)
    found = {}      # ray index -> (lam, sizes, hexrows)
    best = {}
    t0 = time.time()
    for lam in a.lams:
        todo = [i for i in range(R.shape[0]) if i not in found and lam * (int(R[i, [0, 1, 3, 7, 15]].sum()) + int(R[i, 30])) <= a.max_n]
        if not todo:
            continue
        lines = []
        for i in todo:
            sz = party_sizes(R[i], lam)
            lines.append(' '.join(map(str, [i, lam] + sz + [lam * int(x) for x in R[i]])))
        p = subprocess.run([exe, str(a.seconds), str(a.seed + lam), '1.5', '0.08', str(a.sweeps)],
                           input='\n'.join(lines) + '\n', capture_output=True, text=True, env=env)
        if p.returncode != 0:
            raise SystemExit(f'kernel failed: {p.stderr[:500]}')
        nf = 0
        for ln in p.stdout.splitlines():
            t = ln.split()
            i, l = int(t[0]), int(t[1])
            if t[2] == 'FOUND':
                N = int(t[3]); hexrows = t[4:4 + N]
                sz = party_sizes(R[i], l)
                A = rows_to_matrix(hexrows, N)
                if verify(R[i], l, sz, A):
                    found[i] = (l, sz, hexrows); nf += 1
                else:
                    print(f'WARNING: kernel hit for ray {i} (lambda {l}) failed independent verification -- discarded')
            elif t[2] == 'BEST':
                best[(i, l)] = int(t[3])
        print(f'lambda={lam}: {len(todo)} rays tried, {nf} realised and independently verified '
              f'(cumulative {len(found)}/{R.shape[0]}; {time.time() - t0:.0f}s)', flush=True)
    shutil.rmtree(tmp, ignore_errors=True)
    if a.out:
        idx = sorted(found)
        rows = np.zeros((len(idx), 64), dtype=np.uint64)
        for k, i in enumerate(idx):
            for j, h in enumerate(found[i][2]):
                rows[k, j] = np.uint64(int(h, 16))
        np.savez(a.out, index=np.array(idx, dtype=np.int64), rep=R[idx].astype(np.int16) if idx else np.zeros((0, 31), np.int16),
                 lam=np.array([found[i][0] for i in idx], dtype=np.int8),
                 sizes=np.array([found[i][1] for i in idx], dtype=np.int16).reshape(-1, 6), rows=rows)
        json.dump({'rays': int(R.shape[0]), 'found': len(found),
                   'by_lambda': {str(l): sum(1 for i in found if found[i][0] == l) for l in a.lams},
                   'best_residual': {f'{i}:{l}': e for (i, l), e in best.items()}},
                  open(a.out + '.json', 'w'), indent=1)
        print(f'certificates -> {a.out}')


if __name__ == '__main__':
    main()
