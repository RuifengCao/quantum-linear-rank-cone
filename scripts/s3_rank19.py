#!/usr/bin/env python3
"""s3: tight-row rank decision for the 19 HEC_5 extreme-ray orbit reps.
Extreme in the cone cut by H  <=>  rank of tight rows == 30.
--extra-rows adds substituted 6-variable inequalities (from s9) before deciding."""
import argparse, sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import numpy as np
from epr1kit import core

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument('H')
ap.add_argument('--rays', default=None, help='npy; default bundled 19 HEC orbit reps')
ap.add_argument('--extra-rows', default=None, help='npy of extra inequality rows to append')
ap.add_argument('--no-exact', action='store_true', help='skip exact Bareiss confirmation')
a = ap.parse_args()
H = np.load(a.H).astype(np.int64)
if a.extra_rows:
    H = np.unique(np.vstack([H, np.load(a.extra_rows).astype(np.int64)]), axis=0)
rays = np.load(a.rays).astype(np.int64) if a.rays else core.load('hec5_rays_maskorder')
rep = core.tight_rank_report(H, rays, exact_confirm=not a.no_exact)
ext = sum(1 for *_x, e in rep if e)
print(f'H rows {H.shape[0]}  ->  extreme {ext}/{len(rep)}')
for ri, nt, rk, e in rep:
    print(f'  ray #{ri+1:2d}: tight {nt:6d}  rank {rk}  {"EXTREME" if e else "non-extreme"}')
