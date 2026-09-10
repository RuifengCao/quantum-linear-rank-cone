#!/usr/bin/env python3
"""s6: sweep the CLR extreme-ray orbit reps (from s5) through the pure QLR cone:
membership + tight-rank extremality.  Reconciles the '59 = 40 + 17 + 2' ledger:
40 CLR-extreme orbits expected to stay QLR-extreme."""
import argparse, sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import numpy as np
from epr1kit import core

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument('orbit_reps', help='npy from s5')
ap.add_argument('--H', required=True, help='pure QLR H-rep npy (from s1 --variant pure28)')
a = ap.parse_args()
R = np.load(a.orbit_reps).astype(np.int64)
H = np.load(a.H).astype(np.int64)
inside = [(H @ r >= 0).all() for r in R]
print(f'orbit reps {R.shape[0]}  QLR members {sum(inside)}')
rep = core.tight_rank_report(H, R[np.array(inside)], exact_confirm=True)
ext = sum(1 for *_x, e in rep if e)
print(f'QLR-extreme among members: {ext}  (ledger expectation: 40)')
for (ri, nt, rk, e), keep in zip(rep, [i for i, k in enumerate(inside) if k]):
    print(f'  rep {keep:3d}: tight {nt:6d} rank {rk} {"EXTREME" if e else "non-extreme"}')
