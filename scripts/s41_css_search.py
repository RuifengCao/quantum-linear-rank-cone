#!/usr/bin/env python3
"""s41: CSS search for five-party rays (R52).

For a five-party ray r (the 31 entropies S(X), X a nonempty set of parties among A..E) let S be its pure-symmetric
extension to the six systems A..F (F the purifier) and h = h_S its lift, h(X) = S(X) + sum_{i in X} S(i)
(Theorem 1 of the paper; `s17_ci_test.h_norm`).  If an N x m binary matrix M, N = kappa h(E), with column groups
of widths kappa h(i), has column-group ranks r(X) = kappa h(X) for all 63 nonempty X, then the CSS state of the
row space C of M has

    S_C(X) = r(X) + r(E \\ X) - r(E) = 2 kappa S(X)

(R32 note, proof of Theorem 1, step 2, with the self-duality of h_S; remark "symplectic structure is irrelevant").
So such an M is a qubit realisation of r at lambda = 2 kappa, with party i owning 2 kappa S(i) qubits.  Every qubit
realisation of kappa r gives one ([I | W] for a graph state W), so this search has the size of the graph-state
search at kappa and at least as many solutions, while it produces realisations at 2 kappa.

The kernel cc/s41_css_anneal.c anneals the columns.  A reported M is accepted only after three checks that do not
use the kernel:
  1. the 63 column-group ranks of M, recomputed here with Python integers, equal kappa h;
  2. M is turned into a graph state: row-reduce M to [I | A] on its pivot columns P; the bipartite graph between
     P and the other columns Q with biadjacency A is local-Clifford equivalent to the CSS state (Hadamards on Q),
     so it has the same entropies; its 31 visible cut ranks, recomputed here with Python integers, equal
     2 kappa r;
  3. with at most 64 qubits, epr1kit.stabcert.verify_realisation(r, 2 kappa, sizes, rows) confirms the graph.
Graph certificates with at most 64 qubits are written in the format of s16 (fields index, rep, lam, sizes, rows,
plus cols, m, kappa); all hits, of any size, also go to OUT.plain.npz (index, rep, kappa, cols, m), and OUT.json
records rays, found, by_lambda and best residuals as s16 does.  The output files are rewritten after every
accepted hit, so an interrupted run keeps what it found.

  python3 scripts/s41_css_search.py rays.npy --seconds 60 --out certs.npz [--rows 0,5] [--kappa 1]
  python3 scripts/s41_css_search.py rays.npy --budget 22:3600x4,32:1800,64:300 --out certs.npz
      (rays with N <= 22 get 3,600 s in four copies with different seeds, N <= 32 get 1,800 s, N <= 64 get 300 s;
       N = kappa (S(A) + ... + S(F)) is the qubit count at lambda = kappa)
  python3 scripts/s41_css_search.py --check certs.npz rays.npy      # re-verify a stored file (checks 1-3)
  python3 scripts/s41_css_search.py --self-test                     # CSS formula and graph conversion against
                                                                    # state vectors; a short kernel run

A hit is an exact certificate; a miss proves nothing.
"""
import argparse, itertools, json, os, shutil, subprocess, sys, tempfile, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, '..')
sys.path.insert(0, HERE); sys.path.insert(0, ROOT)
import s17_ci_test as T
from s16_stab_search import omp_threads
from epr1kit import stabcert

KERNEL = os.path.join(ROOT, 'cc', 's41_css_anneal.c')


# ----------------------------------------------------------------------------------------------- exact routines
def rank_bits(vecs):
    """GF(2) rank of Python-int bit vectors."""
    basis = {}
    for v in vecs:
        while v:
            p = v.bit_length() - 1
            if p in basis:
                v ^= basis[p]
            else:
                basis[p] = v
                break
    return len(basis)


def lift(r, kappa=1):
    """kappa * h_S on the masks 0..63 (bit i = system i, F = bit 5)."""
    return [kappa * int(x) for x in T.h_norm(np.asarray(r, dtype=np.int64))]


def singles(r):
    return [int(r[(1 << q) - 1]) for q in range(5)] + [int(r[30])]


def groups(h):
    nb = [h[1 << i] for i in range(6)]
    off = [sum(nb[:i]) for i in range(7)]
    return nb, off


def column_ranks(cols, h):
    nb, off = groups(h)
    return {X: rank_bits([cols[c] for i in range(6) if X >> i & 1 for c in range(off[i], off[i + 1])])
            for X in range(1, 64)}


def to_graph(cols, nrows):
    """bipartite graph state (adjacency as Python-int rows over the m columns) with the entropies of the CSS state
    of the row space of M (columns = Python-int bit vectors over nrows rows)."""
    m = len(cols)
    M = [[(cols[c] >> i) & 1 for c in range(m)] for i in range(nrows)]
    piv, r = [], 0
    for c in range(m):
        p = next((i for i in range(r, nrows) if M[i][c]), None)
        if p is None:
            continue
        M[r], M[p] = M[p], M[r]
        for i in range(nrows):
            if i != r and M[i][c]:
                M[i] = [a ^ b for a, b in zip(M[i], M[r])]
        piv.append(c)
        r += 1
    pivset = set(piv)
    rows = [0] * m
    for a, pc in enumerate(piv):
        for q in range(m):
            if q not in pivset and M[a][q]:
                rows[pc] |= 1 << q
                rows[q] |= 1 << pc
    return rows, r


def cut_ranks(rows, sizes):
    """31 visible cut ranks rank_GF2 Gamma[X, X^c] of a graph given by Python-int adjacency rows (party q owns
    sizes[q] consecutive vertices, F last)."""
    off = [sum(sizes[:q]) for q in range(7)]
    n = off[6]
    out = []
    for X in range(1, 32):
        inside = [v for q in range(5) if X >> q & 1 for v in range(off[q], off[q + 1])]
        cmask = sum(1 << v for v in range(n)) ^ sum(1 << v for v in inside)
        out.append(rank_bits([rows[v] & cmask for v in inside]))
    return out


def verify(r, kappa, cols):
    """checks 1-3 of the docstring; returns (ok, details, graph rows, sizes)."""
    h = lift(r, kappa)
    nb, off = groups(h)
    det = {}
    det['shape'] = len(cols) == off[6] and all(c >> h[63] == 0 for c in cols)
    rk = column_ranks(cols, h) if det['shape'] else {}
    det['column_ranks'] = bool(rk) and all(rk[X] == h[X] for X in range(1, 64))
    rows, rank = to_graph(cols, h[63])
    sizes = [2 * kappa * s for s in singles(r)]
    det['graph_cut_ranks'] = det['column_ranks'] and cut_ranks(rows, sizes) == [2 * kappa * int(x) for x in r]
    if sum(sizes) <= 64:
        det['stabcert'] = bool(stabcert.verify_realisation(r, 2 * kappa, sizes, rows))
    return all(det.values()), det, rows, sizes


# ----------------------------------------------------------------------------------------------- kernel
def build_kernel(outdir):
    exe = os.path.join(outdir, 's41_css_anneal')
    if shutil.which('gcc') is None:
        raise SystemExit('gcc not found -- install build tools first')
    for flags in (['-O3', '-march=native', '-fopenmp'], ['-O3', '-fopenmp'], ['-O3']):
        if subprocess.run(['gcc', *flags, '-o', exe, KERNEL, '-lm'], capture_output=True).returncode == 0:
            return exe, flags
    raise SystemExit('could not compile cc/s41_css_anneal.c')


def parse_budget(spec):
    """'22:3600x4,32:1800,64:300' -> [(22, 3600.0, 4), (32, 1800.0, 1), (64, 300.0, 1)]"""
    out = []
    for part in spec.split(','):
        nmax, rest = part.split(':')
        sec, _, cop = rest.partition('x')
        out.append((int(nmax), float(sec), int(cop) if cop else 1))
    return sorted(out)


def save(out, hits, summary):
    ok64 = [hh for hh in hits if sum(hh['sizes']) <= 64]
    pad = lambda v, n: np.array(list(v) + [0] * (n - len(v)), dtype=np.uint64)
    np.savez(out,
             index=np.array([hh['index'] for hh in ok64], dtype=np.int64),
             rep=np.array([hh['rep'] for hh in ok64], dtype=np.int64).reshape(-1, 31),
             lam=np.array([2 * hh['kappa'] for hh in ok64], dtype=np.int64),
             sizes=np.array([hh['sizes'] for hh in ok64], dtype=np.int64).reshape(-1, 6),
             rows=np.array([pad(hh['rows'], 64) for hh in ok64], dtype=np.uint64).reshape(-1, 64),
             cols=np.array([pad(hh['cols'], 128) for hh in ok64], dtype=np.uint64).reshape(-1, 128),
             m=np.array([len(hh['cols']) for hh in ok64], dtype=np.int64),
             kappa=np.array([hh['kappa'] for hh in ok64], dtype=np.int64))
    plain = out[:-4] + '.plain.npz' if out.endswith('.npz') else out + '.plain.npz'
    np.savez(plain,
             index=np.array([hh['index'] for hh in hits], dtype=np.int64),
             rep=np.array([hh['rep'] for hh in hits], dtype=np.int64).reshape(-1, 31),
             kappa=np.array([hh['kappa'] for hh in hits], dtype=np.int64),
             cols=np.array([pad(hh['cols'], 128) for hh in hits], dtype=np.uint64).reshape(-1, 128),
             m=np.array([len(hh['cols']) for hh in hits], dtype=np.int64))
    with open(out + '.json', 'w') as f:
        json.dump(summary, f, indent=1)


def search(R, rows, kappa, plan, seed, threads, out, sweeps, log=print):
    """plan: list of (row, seconds, copies)."""
    hits, found = [], set()
    best = {}
    summary = dict(rays=len(R), rows=len(rows), kappa=kappa, lam=2 * kappa, seed=seed, found=0, by_lambda={},
                   best_residual={}, found_rows=[], seconds={}, note='rows refer to the input rays file')
    tmp = tempfile.mkdtemp(prefix='s41_')
    try:
        exe, flags = build_kernel(tmp)
        nthr = omp_threads(threads)
        log(f"kernel: gcc {' '.join(flags)}; OpenMP threads {nthr if '-fopenmp' in flags else 1}")
        jobs = sorted([(sec, k) for k, sec, cop in plan for _ in range(cop)], key=lambda t: -t[0])
        lines = []
        for sec, k in jobs:
            h = lift(R[k], kappa)
            if h[63] > 64 or sum(h[1 << i] for i in range(6)) > 128:
                log(f'NOTE row {k}: N = {h[63]} is beyond the kernel (N <= 64); skipped')
                continue
            lines.append(f"{k} {sec} " + ' '.join(str(h[X]) for X in range(1, 64)))
        env = dict(os.environ, OMP_NUM_THREADS=str(nthr))
        t0 = time.time()
        p = subprocess.Popen([exe, str(seed), '1.5', '0.08', str(sweeps), '0.5'], stdin=subprocess.PIPE,
                             stdout=subprocess.PIPE, text=True, env=env)
        p.stdin.write('\n'.join(lines) + '\n'); p.stdin.close()
        done = 0
        for line in p.stdout:
            f = line.split()
            if len(f) < 2:
                continue
            k = int(f[0]); done += 1
            if f[1] == 'FOUND':
                m = int(f[5]); cols = [int(x, 16) for x in f[6:6 + m]]
                if k in found:
                    continue
                ok, det, grows, sizes = verify(R[k], kappa, cols)
                if not ok:
                    log(f'WARNING row {k}: kernel hit failed verification {det}; discarded')
                    continue
                found.add(k)
                hits.append(dict(index=k, rep=[int(x) for x in R[k]], kappa=kappa, cols=cols, rows=grows, sizes=sizes))
                summary['found_rows'].append(k); summary['found'] = len(found)
                summary['seconds'][str(k)] = float(f[2])
                summary['by_lambda'][str(2 * kappa)] = len(found)
                summary['best_residual'].pop(f'{k}:{2 * kappa}', None)
                save(out, hits, summary)
                log(f'row {k}: FOUND after {float(f[2]):.1f} s, verified ({", ".join(det)}); '
                    f'{len(found)} found, {done} of {len(lines)} searches done ({time.time() - t0:.0f} s)')
            elif f[1] == 'BEST':
                if k not in found:
                    e = int(f[2]); key = f'{k}:{2 * kappa}'
                    summary['best_residual'][key] = min(e, summary['best_residual'].get(key, e))
                if done % 50 == 0:
                    log(f'{done} of {len(lines)} searches done, {len(found)} found ({time.time() - t0:.0f} s)')
                    save(out, hits, summary)
            elif f[1] == 'SKIP':
                log(f'NOTE row {k}: kernel skipped ({" ".join(f[2:])})')
        p.wait()
        save(out, hits, summary)
        log(f'lambda={2 * kappa} (CSS, kappa={kappa}): {len(rows)} rays tried, {len(found)} realised and independently '
            f'verified ({time.time() - t0:.0f}s)')
        log(f'certificates -> {out} (graphs with at most 64 qubits), all hits -> '
            f'{out[:-4] + ".plain.npz" if out.endswith(".npz") else out + ".plain.npz"}')
        return hits
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


# ----------------------------------------------------------------------------------------------- check / self-test
def check_file(path, rays_path, log=print):
    R = np.load(rays_path).astype(np.int64)
    plain = path[:-4] + '.plain.npz' if path.endswith('.npz') else path + '.plain.npz'
    Z = np.load(plain)
    bad = 0
    for k, kap, cols, m in zip(Z['index'], Z['kappa'], Z['cols'], Z['m']):
        ok, det, _, _ = verify(R[int(k)], int(kap), [int(c) for c in cols[:int(m)]])
        bad += not ok
        log(f'row {int(k)}: {"ok" if ok else "FAILED"} {det}')
    G = np.load(path)
    for k, rep, lam, sizes, rows in zip(G['index'], G['rep'], G['lam'], G['sizes'], G['rows']):
        ok = bool(stabcert.verify_realisation(rep, int(lam), [int(s) for s in sizes], [int(x) for x in rows]))
        bad += not ok
    log(f'{len(Z["index"])} plain certificates, {len(G["index"])} graph certificates: '
        f'{"ALL OK" if not bad else f"{bad} FAILED"}')
    return bad == 0


def css_state(cols, nrows, m):
    """amplitudes of |C> = sum_{c in rowspace} |c> over m qubits (qubit q = column q, little-endian index)."""
    rows = [sum(((cols[c] >> i) & 1) << c for c in range(m)) for i in range(nrows)]
    words = {0}
    for v in rows:
        words |= {w ^ v for w in words}
    psi = np.zeros(1 << m)
    for w in words:
        psi[w] = 1.0
    return psi / np.linalg.norm(psi)


def graph_state(rows, m):
    psi = np.zeros(1 << m)
    for x in range(1 << m):
        e = sum(bin(rows[a] & x & ~((1 << (a + 1)) - 1)).count('1') for a in range(m) if x >> a & 1)
        psi[x] = (-1) ** e
    return psi / np.linalg.norm(psi)


def entropy(psi, m, part):
    t = psi.reshape([2] * m)          # axis j <-> qubit m-1-j
    axes = [m - 1 - q for q in part]
    rest = [a for a in range(m) if a not in axes]
    s = np.linalg.svd(np.transpose(t, axes + rest).reshape(1 << len(part), -1), compute_uv=False)
    p = s ** 2; p = p[p > 1e-12]
    return float(-(p * np.log2(p)).sum())


def self_test(trials=60, seed=41, log=print):
    rng = np.random.default_rng(seed)
    nbad = 0
    for t in range(trials):
        m = int(rng.integers(4, 11)); k = int(rng.integers(1, m + 1))
        nparties = int(rng.integers(2, 7))
        owner = sorted(int(x) for x in rng.integers(0, nparties, size=m))
        cols = [int(sum(int(b) << i for i, b in enumerate(rng.integers(0, 2, size=k)))) for _ in range(m)]
        psi = css_state(cols, k, m)
        rE = rank_bits(cols)
        for X in range(1, 1 << nparties):
            part = [q for q in range(m) if X >> owner[q] & 1]
            comp = [q for q in range(m) if not X >> owner[q] & 1]
            formula = rank_bits([cols[q] for q in part]) + rank_bits([cols[q] for q in comp]) - rE
            if abs(entropy(psi, m, part) - formula) > 1e-8:
                nbad += 1
        # graph conversion: the bipartite graph state with local Hadamards on the non-pivot qubits is |C>
        grows, rk = to_graph(cols, k)
        g = graph_state(grows, m).reshape([2] * m)
        Hd = np.array([[1, 1], [1, -1]]) / np.sqrt(2)
        # the pivot set, as in to_graph
        M = [[(cols[c] >> i) & 1 for c in range(m)] for i in range(k)]; piv, r = [], 0
        for c in range(m):
            p = next((i for i in range(r, k) if M[i][c]), None)
            if p is None:
                continue
            M[r], M[p] = M[p], M[r]
            for i in range(k):
                if i != r and M[i][c]:
                    M[i] = [a ^ b for a, b in zip(M[i], M[r])]
            piv.append(c); r += 1
        for q in range(m):
            if q not in piv:
                g = np.moveaxis(np.tensordot(Hd, g, axes=([1], [m - 1 - q])), 0, m - 1 - q)
        ov = abs(np.dot(g.reshape(-1), psi))
        if abs(ov - 1) > 1e-8:
            nbad += 1
    log(f'  CSS formula and graph conversion against state vectors: {trials} random codes, '
        f'{"all agree" if not nbad else f"{nbad} disagreements"}')
    # a short kernel run on two known rays of the catalogue (rows 8 and 25 of data/stab5_extreme_reps.npy)
    R = np.load(os.path.join(ROOT, 'data', 'stab5_extreme_reps.npy')).astype(np.int64)
    out = os.path.join(tempfile.mkdtemp(prefix='s41t_'), 'st.npz')
    hits = search(R, [8, 25], 1, [(8, 30.0, 1), (25, 30.0, 1)], 4141, 2, out, 1000, log=lambda s: None)
    ok_k = len(hits) == 2
    log(f'  kernel: rows 8 and 25 of the catalogue (N = 10, 14) realised at lambda = 2 and verified: '
        f'{len(hits)} of 2')
    shutil.rmtree(os.path.dirname(out), ignore_errors=True)
    ok = nbad == 0 and ok_k
    log('self-test ' + ('PASSED' if ok else 'FAILED'))
    return ok


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('rays', nargs='?')
    ap.add_argument('--rows', default=None, help='comma-separated rows of the rays file (default: all)')
    ap.add_argument('--kappa', type=int, default=1, help='search for kappa * h_S; hits realise lambda = 2 kappa')
    ap.add_argument('--seconds', type=float, default=60.0)
    ap.add_argument('--budget', default=None, help='per-N budgets, e.g. 22:3600x4,32:1800,64:300 (overrides --seconds)')
    ap.add_argument('--seed', type=int, default=20261119)
    ap.add_argument('--sweeps', type=float, default=1000.0)
    ap.add_argument('--threads', type=int, default=None)
    ap.add_argument('--out', default=None)
    ap.add_argument('--check', default=None, help='re-verify a stored certificate file (with the rays file)')
    ap.add_argument('--self-test', action='store_true')
    a = ap.parse_args()
    if a.self_test:
        sys.exit(0 if self_test() else 1)
    if a.check:
        sys.exit(0 if check_file(a.check, a.rays) else 1)
    if not a.rays or not a.out:
        ap.error('rays file and --out are required')
    R = np.load(a.rays).astype(np.int64)
    rows = [int(x) for x in a.rows.split(',')] if a.rows else list(range(len(R)))
    if a.budget:
        bands = parse_budget(a.budget)
        plan = []
        for k in rows:
            N = a.kappa * sum(singles(R[k]))
            band = next((b for b in bands if N <= b[0]), None)
            if band is None:
                print(f'NOTE row {k}: N = {N} above every band; skipped', flush=True)
                continue
            plan.append((k, band[1], band[2]))
    else:
        plan = [(k, a.seconds, 1) for k in rows]
    search(R, rows, a.kappa, plan, a.seed, a.threads, a.out, a.sweeps, log=lambda s: print(s, flush=True))


if __name__ == '__main__':
    main()
