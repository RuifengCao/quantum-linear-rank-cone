#!/usr/bin/env python3
"""s8: facet-restriction decomposition certificate for the non-extreme ray.
Under pure28 the only remaining non-extreme ray is #19 (default target,
data/targets1.npy); #17/#18 are extreme since R3 and cannot decompose.
For ray r: restrict the pool to vectors tight on ALL rows tight at r, drop
multiples of r, ask LP for r = sum lambda_i v_i (lambda >= 0).
Feasible -> explicit non-extremality certificate (saved as npz).
Historical: --rays with data/targets3.npy re-runs the old #17/#18/#19 trio."""
import argparse, sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import numpy as np
from epr1kit import core

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument('--H', required=True)
ap.add_argument('--pool', required=True, nargs='+',
                help='one or more npy files of candidate vectors; multiple files are merged+deduped (e.g. both s7 layers)')
ap.add_argument('--rays', default=None, help='default: bundled targets1 (= ray #19 only)')
ap.add_argument('--out', default='s8_certs.npz')
a = ap.parse_args()
H = np.load(a.H).astype(np.int64)
pool = np.unique(np.vstack([np.load(p).astype(np.int64) for p in a.pool]), axis=0)
print(f'pool: {pool.shape[0]} unique vectors from {len(a.pool)} file(s)')
rays = np.load(a.rays).astype(np.int64) if a.rays else core.load('targets1')
certs = {}
for k, r in enumerate(rays):
    res = core.decompose_certificate(r, H, pool)
    if res is None:
        print(f'ray {k}: no certificate from this pool')
        continue
    V, lam = res
    keep = lam > 1e-9
    certs[f'ray{k}_V'] = V[keep]
    certs[f'ray{k}_lam'] = lam[keep]
    print(f'ray {k}: CERTIFICATE with {int(keep.sum())} pool vectors')
np.savez(a.out, **certs)
