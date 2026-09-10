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
State (format 2: compact int16/int32 representatives + done/skipped flags; legacy
R15/R17 states are converted on load) is saved after every expansion / chunk;
giants are skipped (reported) on timeout.
A growth curve is appended to <state>.growth.csv (one row per expansion).
--workers N solves N vertex figures in parallel (fork + one lrs each); the budget
is checked between chunks of 4N figures, so wall time may overrun by one chunk.
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
ap.add_argument('--workers', type=int, default=1, help='parallel vertex-figure solves (fork; lrs per worker)')
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
    R0 = np.stack([A[i] for i in reps.values()]).astype(np.int32)
    S = {'format': 2, 'key_width': 'u16', 'reps': R0,
         'done': np.zeros(R0.shape[0], dtype=bool), 'skipped': np.zeros(R0.shape[0], dtype=bool)}
    np.save(a.state, S, allow_pickle=True)
    print(f'init: {A.shape[0]} rays -> {len(reps)} {a.sym.upper()} orbits -> {a.state}')
    sys.exit(0)

S = np.load(a.state, allow_pickle=True).item()
tl = time.time()
if S.get('format') == 2:
    R0 = np.asarray(S['reps'], dtype=np.int64)
    keys0, _ = core.class_reps(R0, perms, offset=0, width='u16')
    reps = {k: R0[i] for i, k in enumerate(keys0)}
    done = {keys0[i] for i in np.nonzero(S['done'])[0]}
    skipped = {keys0[i] for i in np.nonzero(S['skipped'])[0]}
    order_keys = list(keys0)
else:
    # legacy dict-of-lists state (R15/R17): convert, re-keying with the u16 encoding
    legacy = {bytes(k): np.array(v, dtype=np.int64) for k, v in S['reps'].items()}
    old2new = {k: canon(v) for k, v in legacy.items()}
    reps = {old2new[k]: v for k, v in legacy.items()}
    done = {old2new[k] for k in (bytes(x) for x in S['done']) if k in old2new}
    skipped = {old2new[k] for k in (bytes(x) for x in S.get('skipped', [])) if k in old2new}
    order_keys = list(reps.keys())
    print(f'converted legacy state to format 2 (u16 keys): {len(reps)} orbits', flush=True)
print(f'state loaded in {time.time()-tl:.0f}s', flush=True)

def save():
    """compact format 2: int16 (or int32) representatives + boolean done/skipped flags."""
    for k in reps:
        if k not in _pos:
            _pos[k] = len(order_keys); order_keys.append(k)
    R0 = np.stack([reps[k] for k in order_keys])
    R0 = R0.astype(np.int16) if np.abs(R0).max() <= 32767 else R0.astype(np.int32)
    np.save(a.state, {'format': 2, 'key_width': 'u16', 'reps': R0,
                      'done': np.array([k in done for k in order_keys], dtype=bool),
                      'skipped': np.array([k in skipped for k in order_keys], dtype=bool)},
            allow_pickle=True)

_pos = {k: i for i, k in enumerate(order_keys)}

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
print(f'resume: {len(reps)} orbits, {len(done)} expanded, {len(skipped)} skipped, {len(order)} queued; tight range {order[0][0] if order else "-"}..{order[-1][0] if order else "-"}; workers {a.workers}', flush=True)

def growth(nt, nb_count, new):
    with open(a.state + '.growth.csv', 'a') as gl:
        if gl.tell() == 0:
            gl.write('utc,tight,neighbors,new_orbits,total_orbits,expanded,skipped\n')
        gl.write(f'{time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())},{nt},{nb_count},{new},{len(reps)},{len(done)},{len(skipped)}\n')

def merge(k, nt, nb):
    """canonicalise + certify candidates from one expanded representative (main process)."""
    new = bad = 0
    for cand in nb:
        cand = np.asarray(cand, dtype=np.int64)
        cb = canon(cand)
        if cb not in reps:
            if core.rank_mod_p(H[(H @ cand) == 0]) != 30:
                bad += 1; continue
            reps[cb] = cand; new += 1
    done.add(k); growth(nt, len(nb), new)
    print(f'tight={nt}: {len(nb)} neighbors, +{new} orbits, total {len(reps)}, expanded {len(done)}'
          + (f' (bad {bad})' if bad else ''), flush=True)

def _work(item):
    """worker: vertex figure of one representative (globals H, a inherited by fork)."""
    k, nt, v = item
    try:
        nb, _ = neighbors(np.array(v, dtype=np.int64))
        return k, nt, [c.tolist() for c in nb], None
    except subprocess.TimeoutExpired:
        return k, nt, None, 'timeout'

queue = [(k, nt, reps[k].tolist()) for nt, k in order if not (a.max_tight and nt > a.max_tight)]
if a.workers <= 1:
    for k, nt, v in queue:
        if time.time() - t0 > a.budget:
            break
        try:
            nb, _ = neighbors(reps[k])
        except subprocess.TimeoutExpired:
            skipped.add(k); save()
            print(f'tight={nt}: figure timeout -> skipped', flush=True); continue
        merge(k, nt, nb); save()
else:
    from multiprocessing import get_context
    chunk = 4 * a.workers
    i = 0
    with get_context('fork').Pool(a.workers) as pool:   # explicit: Python >= 3.14 defaults to forkserver
        while i < len(queue) and time.time() - t0 <= a.budget:
            items = queue[i:i + chunk]; i += chunk
            for k, nt, nb, err in pool.imap_unordered(_work, items):
                if err:
                    skipped.add(k); print(f'tight={nt}: figure timeout -> skipped', flush=True); continue
                merge(k, nt, nb)
            save()
save()
print(f'STATE: {len(reps)} orbits, {len(done)} expanded, {len(skipped)} skipped'
      + ('  == FULL CLOSURE ==' if len(done) == len(reps) else ''))
