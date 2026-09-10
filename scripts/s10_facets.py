#!/usr/bin/env python3
"""s10: facet classes of a pure QLR H-rep (use QLR_H_pure28.npy).  S6-orbit
reduction, then one redundancy LP per class: h redundant  <=>  h in
cone(all rows minus h's own orbit).  Irredundant classes = facet classes;
reconcile the count against the 31 BCHS inequality classes.
(pure27 measured 83 defining classes in R3; pure28's count is reported live.)"""
import argparse, sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import numpy as np
from scipy.optimize import linprog
from epr1kit import core

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument('--H', required=True, help='pure QLR H-rep npy')
ap.add_argument('--out', default='pure_facet_classes.npz')
ap.add_argument('--limit', type=int, default=0, help='process only the first N classes (smoke test)')
a = ap.parse_args()
H = np.load(a.H).astype(np.int64)
perms = core.perms31_s6()
canon, reps = core.class_reps(H, perms)
print(f'{H.shape[0]} rows -> {len(reps)} S6 classes; running one LP per class')
facet_reps, redundant = [], []
items = sorted(reps.items())
if a.limit > 0:
    items = items[:a.limit]
    print(f'[smoke] limiting to first {len(items)} classes')
for ci, (cb, i0) in enumerate(items):
    keep = np.array([c != cb for c in canon])
    R = H[keep].astype(float)
    t = H[i0].astype(float)
    res = linprog(np.zeros(R.shape[0]), A_eq=R.T, b_eq=t,
                  bounds=[(0, None)] * R.shape[0], method='highs')
    (redundant if res.status == 0 else facet_reps).append(i0)
    if (ci + 1) % 10 == 0:
        print(f'  {ci+1}/{len(items)} done', flush=True)
np.savez(a.out, facet_rows=H[facet_reps], redundant_rows=H[redundant])
print(f'facet classes {len(facet_reps)}  redundant classes {len(redundant)}  '
      f'(BCHS inequality classes: 31)')
