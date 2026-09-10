#!/usr/bin/env python3
"""s7: expand the realizable-vector pool beyond the 760 six-vertex GF(2) vectors.
  A) full enumeration of F_3-weighted 6-vertex graph states (qutrit; CONDITIONAL layer)
  B) sampling of 12-vertex GF(2) graph states with 2-qubit parties (unconditional)
Compiles and drives cc/fp_cutrank_enum.c and cc/gf2_sampler12.c, merges + dedups,
diffs against the bundled 760, and judges everything against a pure H-rep."""
import argparse, sys, os, subprocess, glob
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import numpy as np
from epr1kit import core

HERE = os.path.dirname(os.path.abspath(__file__))
CC = os.path.join(HERE, '..', 'cc')

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument('--H', required=True, help='pure QLR H-rep npy (validity screen)')
ap.add_argument('--gf2-samples', type=int, default=5_000_000)
ap.add_argument('--seed', type=int, default=20260812)
ap.add_argument('--skip-f3', action='store_true')
ap.add_argument('--f3-limit', type=int, default=0, help='enumerate only the first N configs (smoke test)')
ap.add_argument('--outdir', default='s7_out')
a = ap.parse_args()
os.makedirs(a.outdir, exist_ok=True)
H = np.load(a.H).astype(np.int64)
base = core.load('graphstate_vecs')
base_set = {r.tobytes() for r in base.astype(np.int8)}

def compile_and_run(src, out, args):
    import shutil
    if shutil.which('gcc') is None:
        raise SystemExit('gcc not found -- install build tools (e.g. apt install build-essential '
                         'or conda install -c conda-forge c-compiler) and rerun')
    exe = os.path.join(a.outdir, out)
    ok = False
    for flags in (['-O3', '-fopenmp'], ['-O3']):
        r = subprocess.run(['gcc', *flags, '-o', exe, os.path.join(CC, src)])
        if r.returncode == 0:
            ok = ('-fopenmp' in flags)
            break
    else:
        raise SystemExit(f'gcc failed to compile {src} -- see errors above')
    if not ok:
        print(f'[warn] {src} compiled WITHOUT OpenMP; will run single-threaded', flush=True)
    subprocess.run([exe, *args], check=True)

def collect(pattern):
    parts = [np.fromfile(p, dtype=np.uint8).reshape(-1, 31) for p in sorted(glob.glob(pattern))]
    return np.unique(np.vstack(parts), axis=0).astype(np.int64) if parts else np.zeros((0, 31), np.int64)

# --- B: 12-vertex GF(2), unconditional stabilizer layer ---
compile_and_run('gf2_sampler12.c', 'g12', [str(a.gf2_samples), str(a.seed),
                                           os.path.join(a.outdir, 'g12')])
V12 = collect(os.path.join(a.outdir, 'g12*.bin'))
new12 = np.array([r for r in V12 if r.astype(np.int8).tobytes() not in base_set], dtype=np.int64)
mins, cnts = core.judge(H, V12)
np.save(os.path.join(a.outdir, 'realizable_gf2_12v.npy'), V12)
np.save(os.path.join(a.outdir, 'realizable_gf2_12v_new.npy'), new12)
print(f'[B] 12v GF(2): {V12.shape[0]} unique, {new12.shape[0]} beyond the 760, '
      f'pure violations {int((cnts > 0).sum())} (must be 0)')

# --- A: F_3 full enumeration, conditional layer ---
if not a.skip_f3:
    fargs = [os.path.join(a.outdir, 'f3')]
    if a.f3_limit > 0:
        fargs.append(str(a.f3_limit))
    compile_and_run('fp_cutrank_enum.c', 'fpc', fargs)
    V3 = collect(os.path.join(a.outdir, 'f3*.bin'))
    new3 = np.array([r for r in V3 if r.astype(np.int8).tobytes() not in base_set], dtype=np.int64)
    mins, cnts = core.judge(H, V3)
    np.save(os.path.join(a.outdir, 'realizable_f3_conditional.npy'), V3)
    np.save(os.path.join(a.outdir, 'realizable_f3_conditional_new.npy'), new3)
    print(f'[A] F_3 qutrit (CONDITIONAL layer, see README): {V3.shape[0]} unique, '
          f'{new3.shape[0]} beyond the 760, pure violations {int((cnts > 0).sum())}')
