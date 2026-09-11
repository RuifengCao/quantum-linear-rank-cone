#!/usr/bin/env python3
"""s15: chordality filter + simple-tree construction (Hubeny-Rota) for five-party rays.

  python3 scripts/s15_chordal.py <rays.npy> [--out models.npz]

For each ray: SA/SSA check, correlation hypergraph, chordality of its line graph
(Theorem 1 of arXiv:2512.24490: chordal <=> simple-forest realizable),
irreducibility, and for irreducible chordal rays the explicit simple tree from
Algorithm 1 whose min-cut entropies are re-computed and compared coordinate by
coordinate (independent certificate).  Chordal = holographic = inside the
stabilizer cone: an S7-positive verdict.  Non-chordal is NOT a negative verdict.
"""
import argparse, json, os, sys, time
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from epr1kit import chordal

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument('rays')
ap.add_argument('--out', default=None)
ap.add_argument('--limit', type=int, default=None)
ap.add_argument('--skip-sa', action='store_true', help='rays known to lie in QLR5 satisfy SA/SSA; skip the O(63^3) check')
a = ap.parse_args()
R = np.load(a.rays).astype(np.int64)
if a.limit: R = R[:a.limit]
t0 = time.time()
rows = []; models = []
for i, v in enumerate(R):
    r = chordal.analyse(v, check_sa=not a.skip_sa)
    rows.append(r)
    if r.get('verified'):
        m = r['model']; models.append((i, m['edges'], m['weights'], m['n_bulk']))
n = len(rows)
c = sum(1 for r in rows if r.get('chordal'))
irr = sum(1 for r in rows if r.get('irreducible'))
ci = sum(1 for r in rows if r.get('chordal') and r.get('irreducible'))
ver = sum(1 for r in rows if r.get('verified'))
bad = sum(1 for r in rows if r.get('sa_ssa') is False)
print(f'{n} rays in {time.time()-t0:.0f}s: SA/SSA violations {bad}; chordal {c} (irreducible {irr}, chordal&irreducible {ci}); '
      f'simple trees constructed AND verified by min-cut: {ver}/{ci}')
if models and a.out:
    np.savez(a.out, index=np.array([m[0] for m in models]),
             edges=np.array([json.dumps(m[1]) for m in models]),
             weights=np.array([json.dumps(m[2]) for m in models]),
             n_bulk=np.array([m[3] for m in models]))
    print(f'models -> {a.out}')
json.dump([{k: (v if k != 'model' else None) for k, v in r.items()} for r in rows],
          open((a.out or a.rays) + '.chordal.json', 'w'))
