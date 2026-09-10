#!/usr/bin/env python3
"""s4c_qlr: symmetry-aware adjacency decomposition for the QLR_5 cone (A1 campaign).

Same method as s4c_adjacency.py (vertex figures via lrs, exact integer ratio-test
lift, on-the-spot rank-30 certification, BFS over ray orbits), generalised to:
  --H     any H-representation (default data/qlr5_H_facets10860.npy = QLR_5 exact,
          irredundant; 10,860 facet instances of the 31 classes)
  --sym   s5 | s6   (QLR_5 has the full S_6 purification symmetry -> fewer orbits)
Usage:
  python3 scripts/s4c_qlr.py --init data/qlr_seeds60.npy --state qlr_adj.npy
  python3 scripts/s4c_qlr.py --state qlr_adj.npy --budget 3600 --rep-timeout 1800
State is saved after every expansion; giants are skipped (reported) on timeout.
A growth curve is appended to <state>.growth.csv (one row per expansion).
"""
import argparse, os, subprocess, sys, time
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from epr1kit import core

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument('--state', required=True)
ap.add_argument('--init', nargs='*', default=None)
ap.add_argument('--H', default=os.path.join(core.DATA, 'qlr5_H_facets10860.npy'))
ap.add_argument('--sym', default='s6', choices=['s5', 's6'])
ap.add_argument('--budget', type=float, default=600)
ap.add_argument('--rep-timeout', type=float, default=600)
ap.add_argument('--max-tight', type=int, default=None, help='skip figures larger than this (triage)')
a = ap.parse_args()

H = np.load(a.H).astype(np.int64)
perms = core.perms31_s6() if a.sym == 's6' else core.perms31_s5()
NP = perms.shape[0]

def canon(v):
    return core.class_reps(v[None, :], perms, offset=0, width='u16')[0][0]

if a.init:
    A = np.unique(np.vstack([np.load(f).astype(np.int64) for f in a.init]), axis=0)
    g = np.gcd.reduce(np.abs(A), axis=1); g[g == 0] = 1
    A = np.unique(A // g[:, None], axis=0)
    A = A[((H @ A.T) >= 0).all(axis=0)]
    cb, reps = core.class_reps(A, perms, offset=0, width='u16')
    S = {'reps': {bytes(k): A[i].tolist() for k, i in reps.items()}, 'done': [], 'skipped': [], 'key_width': 'u16'}
    np.save(a.state, S, allow_pickle=True)
    print(f'init: {A.shape[0]} rays -> {len(reps)} {a.sym.upper()} orbits -> {a.state}')
    sys.exit(0)

S = np.load(a.state, allow_pickle=True).item()
reps = {bytes(k): np.array(v, dtype=np.int64) for k, v in S['reps'].items()}
done = set(bytes(x) for x in S['done'])
skipped = set(bytes(x) for x in S.get('skipped', []))
if S.get('key_width') != 'u16':
    # legacy u8-keyed state (R15): re-key everything with the u16 encoding
    old2new = {k: canon(v) for k, v in reps.items()}
    reps = {old2new[k]: v for k, v in reps.items()}
    done = {old2new[k] for k in done if k in old2new}
    skipped = {old2new[k] for k in skipped if k in old2new}
    print(f're-keyed legacy state to u16 encoding: {len(reps)} orbits', flush=True)

def save():
    np.save(a.state, {'reps': {k: v.tolist() for k, v in reps.items()},
                      'done': list(done), 'skipped': list(skipped), 'key_width': 'u16'}, allow_pickle=True)

def neighbors(r):
    T = H[(H @ r) == 0]
    ine = f'fig\nH-representation\nbegin\n{T.shape[0]} 32 rational\n' + \
          '\n'.join('0 ' + ' '.join(map(str, row)) for row in T) + '\nend\n'
    p = subprocess.run(['lrs'], input=ine, capture_output=True, text=True, timeout=a.rep_timeout)
    lin, rows = [], []
    for ln in p.stdout.split('\n'):
        t = ln.strip()
        if t.startswith('linearity'):
            lin = [int(x) for x in t.split()[2:]]; continue
        ps = t.split()
        if len(ps) == 32 and ps[0] == '0' and all(x.lstrip('-').isdigit() for x in ps):
            rows.append([int(x) for x in ps[1:]])
    out = []
    hr = H @ r
    for i, d in enumerate(rows, 1):
        if i in lin:
            continue
        d = np.array(d, dtype=np.int64)
        hd = H @ d
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
    return out, T.shape[0]

t0 = time.time()
order = sorted((int((H @ v == 0).sum()), k) for k, v in reps.items() if k not in done and k not in skipped)
print(f'resume: {len(reps)} orbits, {len(done)} expanded, {len(skipped)} skipped, {len(order)} queued; tight range {order[0][0] if order else "-"}..{order[-1][0] if order else "-"}', flush=True)
for nt, k in order:
    if time.time() - t0 > a.budget:
        break
    if a.max_tight and nt > a.max_tight:
        continue
    try:
        nb, _ = neighbors(reps[k])
    except subprocess.TimeoutExpired:
        skipped.add(k); save()
        print(f'tight={nt}: figure timeout -> skipped', flush=True); continue
    new = 0; bad = 0
    for cand in nb:
        cb = canon(cand)
        if cb not in reps:
            if core.rank_mod_p(H[(H @ cand) == 0]) != 30:
                bad += 1; continue
            reps[cb] = cand; new += 1
    done.add(k); save()
    with open(a.state + '.growth.csv', 'a') as gl:
        if gl.tell() == 0:
            gl.write('utc,tight,neighbors,new_orbits,total_orbits,expanded,skipped\n')
        gl.write(f'{time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())},{nt},{len(nb)},{new},{len(reps)},{len(done)},{len(skipped)}\n')
    print(f'tight={nt}: {len(nb)} neighbors, +{new} orbits, total {len(reps)}, expanded {len(done)}'
          + (f' (bad {bad})' if bad else ''), flush=True)
save()
print(f'STATE: {len(reps)} orbits, {len(done)} expanded, {len(skipped)} skipped'
      + ('  == FULL CLOSURE ==' if len(done) == len(reps) else ''))
