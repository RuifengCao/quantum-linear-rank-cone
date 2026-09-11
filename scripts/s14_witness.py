#!/usr/bin/env python3
"""s14: realisability witness scan (S7 test harness).

For every orbit of the campaign catalogue whose maximum coordinate is <= --max-coord
(default 12, the reach of 12-vertex graph states), expand the S6 orbit and look
for a match (up to a multiple 1..3) in the realisable pools: the GF(2) 12-vertex
pool (unconditional), the F3 layer (CONDITIONAL, reported separately), the 760
six-qubit graph states and the 7,943 CLR rays.  The scan uses 64-bit hashing;
every hit is then re-verified coordinate by coordinate.

  python3 scripts/s14_witness.py --state results/<dir>/qlr_adj.npz \\
      --gf2 results/<dir>/realizable_gf2_12v.npy --f3 results/<dir>/realizable_f3_conditional.npy
Writes <state>.witness.npz (rep, witness, mult, perm, source) and prints a summary.
Interpretation guard: unwitnessed != unrealisable (the pools are random samples).
"""
import argparse, os, sys, time
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from epr1kit import core

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument('--state', required=True)
ap.add_argument('--gf2', default=None)
ap.add_argument('--f3', default=None)
ap.add_argument('--max-coord', type=int, default=12)
ap.add_argument('--out', default=None)
a = ap.parse_args()

Z = np.load(a.state); R = Z['reps']
mx = R.max(axis=1)
sel = np.nonzero(mx <= a.max_coord)[0]
S = R[sel].astype(np.int64)
print(f'catalogue {R.shape[0]:,} orbits; {S.shape[0]:,} with max coordinate <= {a.max_coord}')
pools = {'graph760': core.load('graphstate_vecs').astype(np.int64), 'CLR': core.load('clr5_rays7943').astype(np.int64)}
if a.gf2: pools['GF2'] = np.load(a.gf2).astype(np.int64)
if a.f3: pools['F3cond'] = np.load(a.f3).astype(np.int64)
rng = np.random.default_rng(2026); W = rng.integers(1, 2**62, size=31, dtype=np.int64).astype(np.uint64)
def hsh(A): return (A.astype(np.uint64) * W).sum(axis=1)
ph = {k: np.unique(hsh(v)) for k, v in pools.items()}
p6 = core.perms31_s6()
wit = {k: np.zeros(S.shape[0], dtype=bool) for k in pools}
t0 = time.time()
for c0 in range(0, S.shape[0], 1000):
    B = S[c0:c0 + 1000]; m = B.shape[0]
    orb = np.stack([B[:, p6[k]] for k in range(720)], axis=1).reshape(m * 720, 31)
    for mult in (1, 2, 3):
        hv = hsh(orb * mult).reshape(m, 720)
        for k in pools:
            wit[k][c0:c0 + m] |= np.isin(hv, ph[k]).any(axis=1)
print(f'hash scan {time.time() - t0:.0f}s')
seeds = core.load('qlr_seeds60').astype(np.int64)
cs = set(core.class_reps(seeds, p6, offset=0, width='u16')[0])
ks, _ = core.class_reps(S, p6, offset=0, width='u16'); isseed = np.array([k in cs for k in ks])
anyw = np.zeros(S.shape[0], dtype=bool)
for k, v in wit.items():
    anyw |= v; print(f'  {k}: {int(v.sum())} witnessed ({int((v & ~isseed).sum())} beyond the 60 seeds)')
# exact re-verification of every non-seed hit (unconditional pools first)
psets = {nm: {r.tobytes(): i for i, r in enumerate(v.astype(np.int16))} for nm, v in pools.items()}
rec = []
for i in np.nonzero(anyw & ~isseed)[0]:
    v = S[i]; hit = None
    for nm in [n for n in ('GF2', 'graph760', 'CLR', 'F3cond') if n in pools]:
        for k in range(720):
            u = v[p6[k]]
            for mult in (1, 2, 3):
                jdx = psets[nm].get((u * mult).astype(np.int16).tobytes())
                if jdx is not None: hit = (nm, k, mult, jdx); break
            if hit: break
        if hit: break
    if hit is None:
        print(f'WARNING: hash hit not confirmed exactly for orbit {int(sel[i])}'); continue
    rec.append((int(i), hit))
if rec:
    W35 = np.stack([S[i] for i, _ in rec]); V = np.stack([pools[h[0]][h[3]] for _, h in rec])
    M = np.array([h[2] for _, h in rec]); K = np.array([h[1] for _, h in rec]); SRC = np.array([h[0] for _, h in rec])
    assert all(np.array_equal(W35[t][p6[K[t]]] * M[t], V[t]) for t in range(len(rec)))
    out = a.out or (a.state + '.witness.npz')
    np.savez(out, rep=W35.astype(np.int16), witness=V.astype(np.int16), mult=M.astype(np.int8), perm=K.astype(np.int16), source=SRC)
    unc = int((SRC != 'F3cond').sum())
    print(f'new orbits witnessed and exactly verified: {len(rec)} ({unc} unconditional, {len(rec) - unc} via the conditional F3 layer) -> {out}')
    for t in (3, 4, 5, 6, 8, 10, 12):
        msk = (mx[sel] <= t) & ~isseed
        print(f'  new orbits with max coordinate <= {t}: {int(msk.sum())}, witnessed {int((anyw & msk).sum())}')
else:
    print('no new orbit witnessed')
