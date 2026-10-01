#!/usr/bin/env python3
"""s18: cutting-plane pass of the single-common-information test over the A1' catalogue.

One pass:
  1. take the catalogue orbits that are neither cut by the current inequality classes
     nor realised (realisation certificates supplied with --realised, plus the 59-ray ledger);
  2. run the fast single-CI test (disjoint pairs, stop at the first infeasible pair);
  3. turn every infeasibility into an exact certificate (s17_certify) and re-verify it
     with the independent checker (epr1kit.stabcert);
  4. add the new S-form inequalities to the class list (up to S6) and re-screen.
Repeat until a pass finds no new infeasibility.

  python3 scripts/s18_cutting_plane.py --catalogue results/<dir>/qlr_adj.npz \
      --classes data/s17_ineq_classes.npy --realised data/s16_realisations.npz [more.npz ...] \
      --max-n 64 --workers 24 --passes 3 --outdir s18_out

Outputs (outdir): classes.npy (all classes), exclusions.json (new certificates),
survivors.npy (single-CI feasible, unrealised orbits), pass log.  Survivors are
single-CI feasible over DISJOINT pairs only (the fast mode); run s17_ci_test.py
--pairs all on them for the complete single-CI statement.
"""
import argparse, json, os, sys, time
import numpy as np
from multiprocessing import get_context

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..'))
from epr1kit import core, stabcert
import s17_ci_test as ci
from s17_certify import certify, s_form


def images(classes, p6):
    if len(classes) == 0:
        return np.zeros((0, 31), dtype=np.int64)
    return np.unique(np.concatenate([np.stack([c[np.argsort(p6[k])] for k in range(720)]) for c in classes]), axis=0)


def min_values(V, imgs, chunk=1000):
    M = imgs.T.astype(np.float32); out = np.empty(V.shape[0], dtype=np.float32)
    for i in range(0, V.shape[0], chunk):
        out[i:i + chunk] = (V[i:i + chunk].astype(np.float32) @ M).min(axis=1)
    return out


def test_one(r):
    h = ci.h_norm(np.asarray(r, dtype=np.int64))
    for X, Y in ci.pairs_for(h, 'disjoint'):
        if not ci.ci_feasible(h, X, Y):
            return (X, Y)
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--catalogue', required=True)
    ap.add_argument('--classes', required=True)
    ap.add_argument('--realised', nargs='*', default=[])
    ap.add_argument('--max-n', type=int, default=10**9, help='only test orbits with sum_p S_p <= this')
    ap.add_argument('--limit', type=int, default=None)
    ap.add_argument('--workers', type=int, default=2)
    ap.add_argument('--passes', type=int, default=3)
    ap.add_argument('--outdir', default='s18_out')
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    p6 = core.perms31_s6()
    key = lambda A: core.class_reps(np.asarray(A, dtype=np.int64), p6, offset=0, width='u16')[0]
    R = np.load(a.catalogue)['reps'].astype(np.int64)
    sing = [0, 1, 3, 7, 15]
    N = R[:, sing].sum(axis=1) + R[:, 30]
    classes = [c for c in np.load(a.classes).astype(np.int64)]
    real = set(key(core.load('qlr59_reps')))
    for f in a.realised:
        Z = np.load(f); real |= set(key(Z['rep']))
    alive = np.ones(R.shape[0], dtype=bool)
    m = min_values(R, images(classes, p6)); alive &= m > -0.5
    print(f'catalogue {R.shape[0]}: uncut by {len(classes)} classes: {int(alive.sum())}', flush=True)
    allcerts = []
    for ps in range(1, a.passes + 1):
        cand = np.nonzero(alive & (N <= a.max_n))[0]
        if cand.size:
            ks = key(R[cand]); cand = np.array([c for c, k in zip(cand, ks) if k not in real], dtype=np.int64)
        if a.limit:
            cand = cand[:a.limit]
        t0 = time.time()
        with get_context('fork').Pool(a.workers) as pool:
            res = pool.map(test_one, [R[i].tolist() for i in cand], chunksize=8)
        newc = []
        for i, pr in zip(cand, res):
            if pr is None:
                continue
            out = certify(R[i], pr[0], pr[1])
            if out is None:
                print(f'WARNING: no exact certificate for orbit {int(i)} pair {pr}', flush=True); continue
            Fh, y, _ = out
            cS = np.array([int(x) for x in s_form(Fh)], dtype=np.int64)
            ys = {str(k): str(v) for k, v in y.items()}
            ok, _ = stabcert.verify_exclusion(R[i], pr[0], pr[1], ys, F_S=cS.tolist())
            if not ok:
                print(f'WARNING: certificate for orbit {int(i)} failed independent check -- dropped', flush=True); continue
            allcerts.append({'orbit': int(i), 'ray': R[i].tolist(), 'X': pr[0], 'Y': pr[1], 'F_S': cS.tolist(), 'y': ys})
            newc.append(cS)
        ninf = len(newc)
        if newc:
            C = np.vstack(classes + newc)
            keys, reps = core.class_reps(C, p6, offset=64)
            fresh = [C[j] for j in sorted(reps.values()) if j >= len(classes)]
            classes += fresh
            m = min_values(R[alive], images(fresh, p6)); idx = np.nonzero(alive)[0]; alive[idx[m < -0.5]] = False
        print(f'pass {ps}: tested {cand.size}, infeasible {ninf}, new classes -> {len(classes)}, '
              f'uncut now {int(alive.sum())} ({time.time() - t0:.0f}s)', flush=True)
        if ninf == 0:
            break
    surv = np.nonzero(alive & (N <= a.max_n))[0]
    if surv.size:
        ks = key(R[surv]); surv = np.array([c for c, k in zip(surv, ks) if k not in real], dtype=np.int64)
    np.save(os.path.join(a.outdir, 'classes.npy'), np.stack(classes).astype(np.int8))
    np.save(os.path.join(a.outdir, 'survivors.npy'), R[surv].astype(np.int16))
    np.save(os.path.join(a.outdir, 'uncut_mask_packed.npy'), np.packbits(alive))
    json.dump(allcerts, open(os.path.join(a.outdir, 'exclusions.json'), 'w'))
    print(f'DONE: classes {len(classes)}, new certificates {len(allcerts)}, uncut {int(alive.sum())}, '
          f'unrealised survivors (N <= {a.max_n}) {surv.size}', flush=True)


if __name__ == '__main__':
    main()
