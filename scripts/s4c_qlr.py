#!/usr/bin/env python3
"""s4c_qlr: symmetry-aware adjacency decomposition for the QLR_5 cone (A1 campaign).

Same method as s4c_adjacency.py (vertex figures via lrs, exact integer ratio-test
lift, on-the-spot rank-30 certification, BFS over ray orbits), generalised to:
  --H     any H-representation (default data/qlr5_H_facets10860.npy = QLR_5 exact,
          irredundant; 10,860 facet instances of the 31 classes)
  --sym   s5 | s6   (QLR_5 has the full S_6 purification symmetry -> fewer orbits)
Usage:
  python3 scripts/s4c_qlr.py --init data/qlr_seeds60.npy --state qlr_adj.npz
  python3 scripts/s4c_qlr.py --state qlr_adj.npz --budget 3600 --rep-timeout 1800 --workers 24 --max-orbits 1500000
State (format 2: compact int16/int32 representatives + done/skipped flags; legacy
R15/R17 states are converted on load) is saved after every expansion / chunk;
giants are skipped (reported) on timeout.
A growth curve is appended to <state>.growth.csv (one row per expansion).
--order coord / --queue-max-coord C focus the expansion on the small-coordinate (S7 front line)
representatives; --save-every throttles state writes (large catalogues take ~10 s per save).
--workers N solves N vertex figures in parallel (fork + one lrs each); workers also
canonicalise and pre-certify candidates against their fork-time snapshot, so the
main process only merges (R25: the main process was the bottleneck in focus mode); the budget
is checked between chunks of 4N figures, so wall time may overrun by one chunk.
The queue is rebuilt from the current catalogue whenever it runs dry (R21: the
first server session idled 4.4 h because newly found orbits were never queued).
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
ap.add_argument('--max-orbits', type=int, default=None, help='stop cleanly once the catalogue reaches this many orbits (safety valve for file size)')
ap.add_argument('--order', default='tight', choices=['tight', 'coord'], help='queue order: tight rows ascending (default) or max coordinate ascending then tight')
ap.add_argument('--queue-max-coord', type=int, default=None, help='only expand representatives whose max coordinate is <= this (S7 front line)')
ap.add_argument('--save-every', type=float, default=120, help='seconds between state saves (final save always happens)')
ap.add_argument('--max-rounds', type=int, default=None, help='stop after this many queue rounds (controlled experiments)')
a = ap.parse_args()

H = np.load(a.H).astype(np.int64)
perms = core.perms31_s6() if a.sym == 's6' else core.perms31_s5()
NP = perms.shape[0]

def canon(v):
    return core.class_reps(v[None, :], perms, offset=0, width='u16')[0][0]

def write_state(path, R0, done_arr, skipped_arr):
    """compressed .npz (recommended: ~3x smaller, ~22 MB per 10^6 orbits) or legacy .npy format 2."""
    R0 = R0.astype(np.int16) if np.abs(R0).max() <= 32767 else R0.astype(np.int32)
    if path.endswith('.npz'):
        np.savez_compressed(path, reps=R0, done=done_arr, skipped=skipped_arr,
                            meta=np.array(['format=2;key_width=u16']))
    else:
        np.save(path, {'format': 2, 'key_width': 'u16', 'reps': R0, 'done': done_arr,
                       'skipped': skipped_arr}, allow_pickle=True)

if a.init:
    A = np.unique(np.vstack([np.load(f).astype(np.int64) for f in a.init]), axis=0)
    g = np.gcd.reduce(np.abs(A), axis=1); g[g == 0] = 1
    A = np.unique(A // g[:, None], axis=0)
    A = A[((H @ A.T) >= 0).all(axis=0)]
    cb, reps = core.class_reps(A, perms, offset=0, width='u16')
    R0 = np.stack([A[i] for i in reps.values()]).astype(np.int32)
    write_state(a.state, R0, np.zeros(R0.shape[0], dtype=bool), np.zeros(R0.shape[0], dtype=bool))
    print(f'init: {A.shape[0]} rays -> {len(reps)} {a.sym.upper()} orbits -> {a.state}')
    sys.exit(0)

tl = time.time()
if a.state.endswith('.npz'):
    Z = np.load(a.state)
    S = {'format': 2, 'reps': Z['reps'], 'done': Z['done'], 'skipped': Z['skipped']}
else:
    S = np.load(a.state, allow_pickle=True).item()
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

_last_save = [0.0]

def save(force=False):
    if not force and time.time() - _last_save[0] < a.save_every:
        return
    save_now(); _last_save[0] = time.time()

def save_now():
    """compact format 2: int16 (or int32) representatives + boolean done/skipped flags."""
    for k in reps:
        if k not in _pos:
            _pos[k] = len(order_keys); order_keys.append(k)
    R0 = np.stack([reps[k] for k in order_keys])
    write_state(a.state, R0, np.array([k in done for k in order_keys], dtype=bool),
                np.array([k in skipped for k in order_keys], dtype=bool))

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
    """canonicalise (vectorised over all candidates) + certify new ones (main process)."""
    new = bad = 0
    if nb:
        C = np.asarray(nb, dtype=np.int64)
        keys, _ = core.class_reps(C, perms, offset=0, width='u16')
        seen = set()
        for cand, cb in zip(C, keys):
            if cb in reps or cb in seen:
                continue
            if core.rank_mod_p(H[(H @ cand) == 0]) != 30:
                bad += 1; continue
            reps[cb] = cand; seen.add(cb); new += 1
    done.add(k); growth(nt, len(nb), new)
    print(f'tight={nt}: {len(nb)} neighbors, +{new} orbits, total {len(reps)}, expanded {len(done)}'
          + (f' (bad {bad})' if bad else ''), flush=True)

def _work(item):
    """worker: vertex figure of one representative, then canonicalise + pre-certify the
    candidates against the fork-time snapshot of `reps` (globals H, perms, reps inherited
    by fork).  Returns only candidates unknown to the snapshot, with their keys and
    rank verdicts, so the main process does dictionary work only."""
    k, nt, v = item
    try:
        nb, _ = neighbors(np.array(v, dtype=np.int64))
    except subprocess.TimeoutExpired:
        return k, nt, None, 'timeout'
    if not nb:
        return k, nt, [], None
    C = np.asarray(nb, dtype=np.int64)
    keys, _ = core.class_reps(C, perms, offset=0, width='u16')
    out = []; seen = set()
    for cand, cb in zip(C, keys):
        if cb in reps or cb in seen:
            continue
        seen.add(cb)
        out.append((cb, cand.tolist(), core.rank_mod_p(H[(H @ cand) == 0]) == 30))
    return k, nt, (len(nb), out), None

def build_queue():
    """vectorised: tight counts by chunked BLAS matmul; optional coordinate filter and ordering."""
    keys = [k for k in reps if k not in done and k not in skipped]
    if not keys:
        return []
    R0 = np.stack([reps[k] for k in keys])
    mxc = R0.max(axis=1)
    if a.queue_max_coord:
        sel = mxc <= a.queue_max_coord
        keys = [k for k, s_ in zip(keys, sel) if s_]; R0 = R0[sel]; mxc = mxc[sel]
        if not keys:
            return []
    Hf = H.astype(np.float64); nt = np.empty(R0.shape[0], dtype=np.int64)
    for i in range(0, R0.shape[0], 20000):
        nt[i:i + 20000] = (np.abs(Hf @ R0[i:i + 20000].astype(np.float64).T) < 0.5).sum(axis=0)
    if a.max_tight:
        sel = nt <= a.max_tight
        keys = [k for k, s_ in zip(keys, sel) if s_]; R0 = R0[sel]; nt = nt[sel]; mxc = mxc[sel]
    idx = np.lexsort((nt, mxc)) if a.order == 'coord' else np.lexsort((mxc, nt))
    return [(keys[i], int(nt[i]), R0[i].tolist()) for i in idx]

def budget_left():
    return time.time() - t0 <= a.budget and not (a.max_orbits and len(reps) >= a.max_orbits)

rounds = 0
while budget_left():
    queue = build_queue()
    if not queue:
        print('queue exhausted (every orbit within --max-tight expanded or skipped)', flush=True); break
    if a.max_rounds and rounds >= a.max_rounds:
        print(f'max rounds ({a.max_rounds}) reached', flush=True); break
    rounds += 1
    print(f'queue round {rounds}: {len(queue)} figures (tight {queue[0][1]}..{queue[-1][1]})', flush=True)
    if a.workers <= 1:
        for k, nt, v in queue:
            if not budget_left():
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
            while i < len(queue) and budget_left():
                items = queue[i:i + chunk]; i += chunk
                for k, nt, payload, err in pool.imap_unordered(_work, items):
                    if err:
                        skipped.add(k); print(f'tight={nt}: figure timeout -> skipped', flush=True); continue
                    n_nb, cands = payload if payload else (0, [])
                    new = bad = 0
                    for cb, cand, ok in cands:
                        if cb in reps:
                            continue
                        if not ok:
                            bad += 1; continue
                        reps[cb] = np.asarray(cand, dtype=np.int64); new += 1
                    done.add(k); growth(nt, n_nb, new)
                    print(f'tight={nt}: {n_nb} neighbors, +{new} orbits, total {len(reps)}, expanded {len(done)}'
                          + (f' (bad {bad})' if bad else ''), flush=True)
                save()
save(force=True)
print(f'STATE: {len(reps)} orbits, {len(done)} expanded, {len(skipped)} skipped'
      + ('  == FULL CLOSURE ==' if len(done) == len(reps) else ''))
