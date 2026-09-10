#!/usr/bin/env python3
"""s2: judge an H-rep against realizable vectors (default: bundled 760 graph-state vectors).
A valid quantum-side H-rep must have zero violations; CLR is expected to violate exactly
the 5 classical monotonicity rows."""
import argparse, sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import numpy as np
from epr1kit import core

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument('H')
ap.add_argument('--vecs', default=None, help='npy of realizable vectors; default bundled 760')
a = ap.parse_args()
H = np.load(a.H).astype(np.int64)
V = np.load(a.vecs).astype(np.int64) if a.vecs else core.load('graphstate_vecs')
mins, cnts = core.judge(H, V)
bad = np.nonzero(cnts)[0]
print(f'H rows {H.shape[0]}  vectors {V.shape[0]}  violated rows {bad.size}')
for i in bad[:20]:
    print(f'  row {i}: min={mins[i]}  #viol={cnts[i]}')
if bad.size > 20:
    print(f'  ... {bad.size-20} more')
