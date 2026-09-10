#!/usr/bin/env python3
"""s4c: symmetry-aware adjacency decomposition -- the engine that actually
delivered the 162 CLR_5 extreme-ray orbits (R7, 2026-08-13).

Method: BFS on the ray-orbit adjacency graph.  For an extreme ray r, its
vertex figure {T x >= 0} (T = rows of CLR_H_fixed tight at r) is enumerated
with lrs; each figure ray d lifts to the adjacent extreme ray via the exact
integer ratio test  r' = (h*.r) d - (h*.d) r.  Every new ray is certified on
the spot (tight-rank 30).  By Balinski connectivity, full closure (every rep
expanded, no new orbit) is a completeness certificate.

R7 status: 162 orbits found (matches DFZ), orbit sizes sum to 7,943 rays
(matches DFZ), 122/162 figures fully expanded.  The ~40 heaviest figures
(tight up to 1,904 -- the hyper-symmetric rays) are as hard as the original
cone; completeness certification therefore goes via (a) rays5 diff against
DFZ's published list, or (b) an mplrs full run, or (c) attacking the giant
figures separately.  See README roadmap.

Usage:
  python3 s4c_adjacency.py --init rays1.npy [rays2.npy ...] --state adj_state.npy
  python3 s4c_adjacency.py --state adj_state.npy --budget 600
State is saved after every expansion; rerun to resume.  Requires lrs on PATH.
"""
import argparse, os, subprocess, sys, time
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from epr1kit import core

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument('--state', required=True)
ap.add_argument('--init', nargs='*', default=None, help='ray .npy files to seed a NEW state')
ap.add_argument('--budget', type=float, default=600, help='seconds for this batch')
ap.add_argument('--rep-timeout', type=float, default=600, help='lrs timeout per vertex figure')
a = ap.parse_args()

Hs = core.load('CLR_H_fixed')
perms = core.perms31_s5()

def canon(v):
    return core.class_reps(v[None, :], perms, offset=0)[0][0]

if a.init:
    A = np.unique(np.vstack([np.load(f).astype(np.int64) for f in a.init]), axis=0)
    g = np.gcd.reduce(np.abs(A), axis=1); g[g == 0] = 1
    A = np.unique(A // g[:, None], axis=0)
    A = A[((Hs @ A.T) >= 0).all(axis=0)]
    cb, reps = core.class_reps(A, perms, offset=0)
    S = {'reps': {bytes(k): A[i].tolist() for k, i in reps.items()}, 'done': []}
    np.save(a.state, S, allow_pickle=True)
    print(f'init: {A.shape[0]} rays -> {len(reps)} orbits -> {a.state}')
    sys.exit(0)

S = np.load(a.state, allow_pickle=True).item()
reps = {bytes(k): np.array(v, dtype=np.int64) for k, v in S['reps'].items()}
done = set(bytes(x) for x in S['done'])

def save():
    np.save(a.state, {'reps': {k: v.tolist() for k, v in reps.items()},
                      'done': list(done)}, allow_pickle=True)

def neighbors(r):
    T = Hs[(Hs @ r) == 0]
    ine = f'fig\nH-representation\nbegin\n{T.shape[0]} 32 rational\n' + \
          '\n'.join('0 ' + ' '.join(map(str, row)) for row in T) + '\nend\n'
    p = subprocess.run(['lrs'], input=ine, capture_output=True, text=True,
                       timeout=a.rep_timeout)
    lin, rows = [], []
    for ln in p.stdout.split('\n'):
        t = ln.strip()
        if t.startswith('linearity'):
            lin = [int(x) for x in t.split()[2:]]; continue
        ps = t.split()
        if len(ps) == 32 and ps[0] == '0' and all(x.lstrip('-').isdigit() for x in ps):
            rows.append([int(x) for x in ps[1:]])
    out = []
    hr = Hs @ r
    for i, d in enumerate(rows, 1):
        if i in lin:
            continue
        d = np.array(d, dtype=np.int64)
        hd = Hs @ d
        if (hd >= 0).all():
            cand = d
        else:
            m = hd < 0
            j = np.argmin(hr[m] / (-hd[m]).astype(float))
            hidx = np.nonzero(m)[0][j]
            cand = int(hr[hidx]) * d - int(hd[hidx]) * r
        g = np.gcd.reduce(np.abs(cand))
        if g == 0:
            continue
        cand //= g
        if cand.sum() < 0:
            cand = -cand
        out.append(cand)
    return out

t0 = time.time()
order = sorted((int((Hs @ v == 0).sum()), k) for k, v in reps.items() if k not in done)
print(f'resume: {len(reps)} orbits, {len(done)} expanded, {len(order)} to go')
for nt, k in order:
    if time.time() - t0 > a.budget:
        break
    try:
        nb = neighbors(reps[k])
    except subprocess.TimeoutExpired:
        print(f'tight={nt}: figure lrs timeout, SKIPPED (closure incomplete)', flush=True)
        continue
    new = 0
    for cand in nb:
        cb = canon(cand)
        if cb not in reps:
            if core.rank_mod_p(Hs[(Hs @ cand) == 0]) != 30:
                print('WARNING: neighbor tight-rank != 30', cand, flush=True); continue
            reps[cb] = cand; new += 1
    done.add(k); save()
    print(f'tight={nt}: {len(nb)} neighbors, +{new} orbits, total {len(reps)}, closed {len(done)}', flush=True)
save()
print(f'STATE: {len(reps)} orbits, {len(done)}/{len(reps)} expanded'
      + ('  == FULL CLOSURE ==' if len(done) == len(reps) else ''))
