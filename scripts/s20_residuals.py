#!/usr/bin/env python3
"""s20: where does the annealing get stuck?  Residual analysis for unrealised rays (R29).

Runs the s16 annealing kernel (cc/s16_anneal.c, dump_best = 1) on the given rays at a
fixed lambda, several independent runs with different seeds, and recomputes the 31 cut
ranks of the best graph state of every run with the independent Python implementation
of s16_stab_search.  The residual of a run is  Delta_X = rank_X - lambda * r_X  (31
entries, mask order); E = sum |Delta_X| is the kernel's energy.

  python3 scripts/s20_residuals.py rays.npy [--index 0 5 7] --lam 1 --seconds 20 --runs 3 --out res.json

This is a diagnostic: a residual that never reaches 0 is NOT a proof of anything.
"""
import argparse, json, os, shutil, subprocess, sys, tempfile
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from s16_stab_search import build_kernel, omp_threads, party_sizes, rows_to_matrix, cut_vector, verify


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('rays')
    ap.add_argument('--index', type=int, nargs='*', default=None, help='rows of the ray file to use (default: all)')
    ap.add_argument('--lam', type=int, default=1)
    ap.add_argument('--seconds', type=float, default=20.0)
    ap.add_argument('--runs', type=int, default=3)
    ap.add_argument('--seed', type=int, default=20261001)
    ap.add_argument('--threads', type=int, default=None)
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    R = np.load(a.rays).astype(np.int64)
    idx = list(range(R.shape[0])) if a.index is None else a.index
    tmp = tempfile.mkdtemp(prefix='s20_')
    exe, flags = build_kernel(tmp)
    env = dict(os.environ); env['OMP_NUM_THREADS'] = str(omp_threads(a.threads))
    lines = []
    for i in idx:
        sz = party_sizes(R[i], a.lam)
        lines.append(' '.join(map(str, [i, a.lam] + sz + [a.lam * int(x) for x in R[i]])))
    out = {i: [] for i in idx}
    for run in range(a.runs):
        p = subprocess.run([exe, str(a.seconds), str(a.seed + 1000 * run + a.lam), '1.5', '0.08', '3000', '1'],
                           input='\n'.join(lines) + '\n', capture_output=True, text=True, env=env)
        if p.returncode != 0:
            raise SystemExit(f'kernel failed: {p.stderr[:500]}')
        for ln in p.stdout.splitlines():
            t = ln.split()
            i = int(t[0]); sz = party_sizes(R[i], a.lam)
            if t[2] == 'FOUND':
                N = int(t[3]); A = rows_to_matrix(t[4:4 + N], N)
                rec = {'run': run, 'E': 0, 'verified': bool(verify(R[i], a.lam, sz, A)), 'delta': {},
                       'sizes': [int(x) for x in sz], 'rows': t[4:4 + N]}   # a certificate: keep it
            else:
                E = int(t[3]); N = int(t[6]); A = rows_to_matrix(t[7:7 + N], N)
                d = cut_vector(A, sz) - a.lam * R[i]
                assert int(np.abs(d).sum()) == E, f'ray {i}: kernel energy {E} != recomputed {int(np.abs(d).sum())}'
                rec = {'run': run, 'E': E, 'delta': {str(m + 1): int(d[m]) for m in range(31) if d[m]}}
            out[i].append(rec)
        print(f'run {run + 1}/{a.runs} done', flush=True)
    shutil.rmtree(tmp, ignore_errors=True)
    for i in idx:
        print(f'ray {i}: ' + '; '.join(f"E={r['E']} " + ','.join(f'{m}:{v:+d}' for m, v in r['delta'].items())
                                      for r in out[i]))
    if a.out:
        json.dump({'lam': a.lam, 'seconds': a.seconds, 'runs': a.runs, 'kernel_flags': flags,
                   'rays': {str(i): out[i] for i in idx}}, open(a.out, 'w'), indent=1)


if __name__ == '__main__':
    main()
