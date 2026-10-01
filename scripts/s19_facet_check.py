#!/usr/bin/env python3
"""s19: facet check of stabilizer-cone inequalities against known realisable vectors.

An inequality c.S >= 0 valid on the stabilizer cone defines a facet of it as soon as
30 linearly independent realisable vectors attain equality (the cone is
31-dimensional).  This script computes, for every class c, the exact rank of the
known realisable vectors on the hyperplane c.S = 0.  The candidate set is the
S6-closure of the known extreme rays of Stab5 (they span every facet they lie on),
plus the 760 six-qubit graph states and, optionally, the GF(2) 12-qubit pool.

  python3 scripts/s19_facet_check.py data/s17_ineq_classes.npy data/stab5_extreme_reps.npy \
      [--pool results/2026-09-11_s7-pools/realizable_gf2_12v.npy] --out ranks.npy
Rank 30 = certified facet of the stabilizer cone; < 30 = undecided (not a proof of non-facet).
"""
import argparse, os, sys, time
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from epr1kit import core

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument('classes')
ap.add_argument('extreme')
ap.add_argument('--pool', default=None)
ap.add_argument('--out', default=None)
a = ap.parse_args()
p6 = core.perms31_s6()
C = np.load(a.classes).astype(np.int64)
E = np.load(a.extreme).astype(np.int64)
V = np.unique(np.concatenate([E[:, p6[k]] for k in range(720)] + [core.load('graphstate_vecs').astype(np.int64)]), axis=0)
if a.pool:
    V = np.unique(np.vstack([V, np.load(a.pool).astype(np.int64)]), axis=0)
Vf = V.astype(np.float32)
print(f'{C.shape[0]} classes; {V.shape[0]} realisable candidate vectors', flush=True)
t0 = time.time(); ranks = np.zeros(C.shape[0], dtype=np.int8)
for i, c in enumerate(C):
    z = np.flatnonzero(np.abs(Vf @ c.astype(np.float32)) < 0.5)
    if z.size == 0:
        continue
    T = V[z]
    if T.shape[0] > 4000:
        T = T[np.random.default_rng(i).choice(T.shape[0], 4000, replace=False)]
    rk = core.rank_mod_p(T)
    if rk < 30 and z.size > 4000:        # retry with all tight vectors (mod-p rank is exact for 0/1-free? use Gram for safety)
        rk = core.rank_exact_gram(V[z])
    ranks[i] = rk
print(f'done in {time.time() - t0:.0f}s; facets certified: {int((ranks == 30).sum())}/{C.shape[0]}; '
      f'rank distribution: {dict(zip(*np.unique(ranks, return_counts=True)))}', flush=True)
if a.out:
    np.save(a.out, ranks)
