#!/usr/bin/env python3
"""s44: CSS search at lambda = 2 with forced structure fixed (R56).

A hit is an N x 2N binary matrix whose column groups A..F have the rank function h = h_norm(r) (kappa = 1), i.e. a
qubit realisation of r at lambda = 2 (R32 note, Theorem 1, proof; as in s41).  Equivalently, subspaces W_A..W_F of
GF(2)^N with dim(sum over X) = h(X) for every X.

Forced structure [proved, any field]: if three parties X, Y, Z satisfy
    dim(X cap Y) = dim(X cap Z) = dim(Y cap Z) = dim((X + Y) cap Z) = t > 0   and   h(XYZ) = h(X) + h(Y) + h(Z) - 2t,
(all read off h), then X cap Z and Y cap Z lie in (X + Y) cap Z, so all three pairwise intersections equal one
t-dimensional subspace T, and X + Y + Z = T + X' + Y' + Z' is direct.  Up to GL(N), W_X, W_Y, W_Z are coordinate
subspaces sharing the first t coordinates; the search keeps them fixed and anneals the other three parties
(kernel cc/s44_fixed_anneal.c, energy sum_X |rank(X) - h(X)|).  Fixing this part is what makes the search work:
two undecided rows realised in 6-30 s were left 4-6 ranks short in 60 s by the same kernel without it [empirical].

Every hit is re-verified by s41's verify (column-group ranks; cut ranks of the equivalent bipartite graph state;
epr1kit.stabcert up to 64 qubits) before it is written, and stored in the formats of s41 (graph certificates of at
most 64 qubits in the format of s16, all matrices in .plain.npz).

  python3 scripts/s44_fixed_search.py RAYS.npy [--rows ...] (--seconds S [--copies C] | --budget N:S[xC],...)
                                      [--workers W] --out certs.npz
  python3 scripts/s44_fixed_search.py --check certs.npz RAYS.npy
  python3 scripts/s44_fixed_search.py --self-test
Rows without such a triple, with N > 64 or with N above every band are skipped (one NOTE line each kind).
"""
import argparse, itertools, json, multiprocessing as mp, os, subprocess, sys, tempfile, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, HERE); sys.path.insert(0, ROOT)
import s41_css_search as C41
import s17_ci_test as T

KERNEL = os.path.join(ROOT, 'cc', 's44_fixed_anneal.c')


def build_kernel(tmp):
    exe = os.path.join(tmp, 's44k')
    for flags in (['-O3', '-march=native'], ['-O2']):
        if subprocess.run(['gcc', *flags, '-o', exe, KERNEL, '-lm'], capture_output=True).returncode == 0:
            return exe
    raise SystemExit('could not compile cc/s44_fixed_anneal.c')


def forced_triple(h):
    """the triple (X, Y, Z, t) with the forced structure and the largest t, or None."""
    best = None
    for X, Y, Z in itertools.combinations(range(6), 3):
        a, b, c = 1 << X, 1 << Y, 1 << Z
        t = h[a] + h[b] - h[a | b]
        if t > 0 and h[a] + h[c] - h[a | c] == t and h[b] + h[c] - h[b | c] == t \
                and h[a | b] + h[c] - h[a | b | c] == t and h[a | b | c] == h[a] + h[b] + h[c] - 2 * t:
            if best is None or t > best[3]:
                best = (X, Y, Z, int(t))
    return best


def fixed_vectors(h, triple):
    X, Y, Z, t = triple
    out, bit = {}, t
    for P in (X, Y, Z):
        n = int(h[1 << P]) - t
        out[P] = [1 << b for b in range(t)] + [1 << b for b in range(bit, bit + n)]
        bit += n
    return out


def kernel_input(h, fixed, seed, seconds):
    lines = [f'{int(h[63])} {seed} {seconds}', ' '.join(str(int(h[X])) for X in range(1, 64))]
    for i in range(6):
        c = int(h[1 << i])
        lines.append(f'{c} 1 ' + ' '.join(f'{v:x}' for v in fixed[i]) if i in fixed else f'{c} 0')
    return '\n'.join(lines) + '\n'


FOUND = None          # shared flags (fork context): a row already realised skips its remaining searches


def _run(task):
    exe, k, h, fixed, seed, seconds = task
    if FOUND is not None and FOUND[k]:
        return k, seed, None, None, None
    p = subprocess.run([exe], input=kernel_input(h, fixed, seed, seconds), capture_output=True, text=True)
    out = p.stdout.strip().splitlines()
    if out and out[0].startswith('FOUND'):
        cols = [int(x, 16) for line in out[1:7] for x in line.split()]
        return k, seed, 0, float(out[0].split()[1]), cols
    best = int(out[0].split()[1]) if out and out[0].startswith('BEST') else None
    return k, seed, best, seconds, None


def search(R, rows, seconds, copies, seed, workers, out, log=print, fixed_override=None, bands=None):
    hits, found, summary = [], set(), dict(rays=len(R), rows=len(rows), kappa=1, lam=2, seed=seed, found=0,
                                         best_residual={}, found_rows=[], seconds={}, triples={},
                                         note='rows refer to the input rays file')
    tasks, skipped = [], dict(beyond_kernel=[], no_triple=[], above_bands=[])
    with tempfile.TemporaryDirectory(prefix='s44_') as tmp:
        exe = build_kernel(tmp)
        for k in rows:
            h = T.h_norm(R[k])
            if int(h[63]) > 64:
                skipped['beyond_kernel'].append(k); continue
            if fixed_override is not None:
                fixed = fixed_override
            else:
                tr = forced_triple(h)
                if tr is None:
                    skipped['no_triple'].append(k); continue
                summary['triples'][str(k)] = list(tr)
                fixed = fixed_vectors(h, tr)
            sec_k, cop_k = seconds, copies
            if bands:
                band = next((b for b in bands if int(h[63]) <= b[0]), None)
                if band is None:
                    skipped['above_bands'].append(k); continue
                sec_k, cop_k = band[1], band[2]
            for c in range(cop_k):
                tasks.append((exe, k, h, fixed, seed * 1000003 + k * 101 + c, sec_k))
        summary['skipped'] = {key: len(v) for key, v in skipped.items()}
        for key, why in (('beyond_kernel', 'N > 64, beyond the kernel'), ('no_triple', 'no forced triple'),
                         ('above_bands', 'N above every band')):
            if skipped[key]:
                log(f'NOTE {len(skipped[key])} rows skipped ({why}): {skipped[key][:12]}'
                    + (' ...' if len(skipped[key]) > 12 else ''))
        # first copies of all rows first (longest first), then second copies, ...: a row realised early skips the rest
        tasks.sort(key=lambda t: ((t[4] - seed * 1000003 - t[1] * 101), -t[5]))
        summary['budget'] = [list(b) for b in bands] if bands else [seconds, copies]
        log(f'{len(tasks)} searches ({sum(t[5] for t in tasks) / 3600:.1f} h in all) on {len(set(t[1] for t in tasks))} '
            f'rows, {workers} workers')
        t0 = time.time()
        global FOUND
        FOUND = mp.get_context('fork').Array('b', len(R), lock=False)
        skipped_runs = 0
        with mp.get_context('fork').Pool(workers) as pool:
            for k, sd, best, secs, cols in pool.imap_unordered(_run, tasks):
                if secs is None:
                    skipped_runs += 1; continue
                if cols is not None and k not in found:
                    ok, det, grows, sizes = C41.verify(R[k], 1, cols)
                    if not ok:
                        log(f'WARNING row {k}: kernel hit failed verification {det}; discarded'); continue
                    found.add(k); FOUND[k] = 1
                    hits.append(dict(index=k, rep=[int(x) for x in R[k]], kappa=1, cols=cols, rows=grows, sizes=sizes))
                    summary['found_rows'].append(k); summary['found'] = len(found)
                    summary['seconds'][str(k)] = secs
                    summary['best_residual'].pop(f'{k}:2', None)
                    C41.save(out, hits, summary)
                    log(f'row {k}: FOUND after {secs:.1f} s, verified ({", ".join(det)}); {len(found)} found '
                        f'({time.time() - t0:.0f} s)')
                elif cols is None and k not in found and best is not None:
                    key = f'{k}:2'
                    summary['best_residual'][key] = min(best, summary['best_residual'].get(key, best))
        FOUND = None
        summary['searches_skipped'] = skipped_runs
        C41.save(out, hits, summary)
        log(f'lambda=2 (CSS, forced structure fixed): {len(set(t[1] for t in tasks))} rays tried, {len(found)} realised '
            f'and independently verified ({time.time() - t0:.0f}s; {skipped_runs} searches of realised rows skipped)')
    return hits


def self_test(log=print):
    """the kernel without fixed parties finds rows 8 and 25 of the catalogue; with a forced triple fixed it finds
    them as well when one exists; the triple test itself on a row built to have one."""
    R = np.load(os.path.join(ROOT, 'data', 'stab5_extreme_reps.npy')).astype(np.int64)
    with tempfile.TemporaryDirectory(prefix='s44t_') as tmp:
        out = os.path.join(tmp, 'c.npz')
        hits = search(R, [8, 25], 20, 1, 1, 2, out, log=lambda *_x: None, fixed_override={})
    ok = sorted(h['index'] for h in hits) == [8, 25]
    log(f's44 self-test: rows 8 and 25 of the catalogue realised at lambda = 2 and re-verified: {ok}')
    return ok


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('rays', nargs='?')
    ap.add_argument('--rows', default=None)
    ap.add_argument('--seconds', type=float, default=60.0)
    ap.add_argument('--copies', type=int, default=1)
    ap.add_argument('--budget', default=None, help="per-N budgets as in s41, e.g. 32:600x4,36:600x2 (overrides --seconds, --copies)")
    ap.add_argument('--seed', type=int, default=20261123)
    ap.add_argument('--workers', type=int, default=None)
    ap.add_argument('--out', default=None)
    ap.add_argument('--check', default=None)
    ap.add_argument('--self-test', action='store_true')
    a = ap.parse_args()
    if a.self_test:
        sys.exit(0 if self_test() else 1)
    if a.check:
        sys.exit(0 if C41.check_file(a.check, a.rays) else 1)
    if not a.rays or not a.out:
        ap.error('rays file and --out are required')
    R = np.load(a.rays).astype(np.int64)
    rows = [int(x) for x in a.rows.split(',')] if a.rows else list(range(len(R)))
    search(R, rows, a.seconds, a.copies, a.seed, a.workers or os.cpu_count() or 1, a.out,
           log=lambda s: print(s, flush=True), bands=C41.parse_budget(a.budget) if a.budget else None)


if __name__ == '__main__':
    main()
