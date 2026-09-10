"""EPR-1 R3 parallel layer.

Design: workers return byte-sets (or write shard files); the parent merges.
Resume: shard files are written tmp -> rename; existing complete shards are skipped.
RSS guard: set EPR1_MAX_RSS_GB to abort a worker whose resident set exceeds it.
"""
import os, resource, sys, time
import numpy as np
from multiprocessing import Pool

def _progress(done, total, t0, tag, last=[0.0]):
    now = time.time()
    if done == total or now - last[0] >= 5.0:
        last[0] = now
        rate = done / max(now - t0, 1e-9)
        eta = (f'{(total - done) / rate:.0f}s' if done >= 2 and now - t0 >= 2 else '--')
        print(f'[{tag}] {done}/{total} chunks  {now-t0:.0f}s elapsed  ETA {eta}',
              file=sys.stderr, flush=True)

def _guard():
    lim = os.environ.get('EPR1_MAX_RSS_GB')
    if lim:
        rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1024 ** 2)  # GB (Linux: KB base)
        if rss > float(lim):
            raise MemoryError(f'worker RSS {rss:.1f} GB > EPR1_MAX_RSS_GB={lim}')

_F = None
def _init(fn):
    global _F
    _F = fn

def _call(args):
    _guard()
    return _F(args)

def pmap_byteset(fn, arglist, workers=None):
    """map fn over arglist in a pool; fn returns a set of bytes; union them."""
    workers = workers or os.cpu_count()
    out = set()
    if workers <= 1:
        for a in arglist:
            out |= fn(a)
        return out
    t0, n = time.time(), len(arglist)
    with Pool(workers, initializer=_init, initargs=(fn,), maxtasksperchild=64) as pool:
        for k, s in enumerate(pool.imap_unordered(_call, arglist, chunksize=1), 1):
            out |= s
            _progress(k, n, t0, 'pmap')
    return out

def pmap_shards(fn, arglist, shard_dir, tag, workers=None, dim=31):
    """resumable variant: each task i writes shard_dir/{tag}_{i:05d}.npy (int8 rows).
    Returns merged unique int64 array."""
    os.makedirs(shard_dir, exist_ok=True)
    todo = []
    for i, a in enumerate(arglist):
        p = os.path.join(shard_dir, f'{tag}_{i:05d}.npy')
        if not os.path.exists(p):
            todo.append((i, a, p))
    def work(item):
        i, a, p = item
        rows = fn(a)
        arr = (np.frombuffer(b''.join(sorted(rows)), dtype=np.int8).reshape(-1, dim)
               if rows else np.zeros((0, dim), dtype=np.int8))
        tmp = p + '.tmp'
        with open(tmp, 'wb') as fh:
            np.save(fh, arr)
        os.replace(tmp, p)
        return p
    workers = workers or os.cpu_count()
    t0, n = time.time(), len(todo)
    print(f'[shards:{tag}] {len(arglist)-n}/{len(arglist)} already done, {n} to run',
          file=sys.stderr, flush=True)
    if workers <= 1:
        for k, it in enumerate(todo, 1):
            work(it)
            _progress(k, n, t0, f'shards:{tag}')
    else:
        with Pool(workers, initializer=_init, initargs=(work,), maxtasksperchild=8) as pool:
            for k, _ in enumerate(pool.imap_unordered(_call, todo, chunksize=1), 1):
                _progress(k, n, t0, f'shards:{tag}')
    parts = []
    for i in range(len(arglist)):
        p = os.path.join(shard_dir, f'{tag}_{i:05d}.npy')
        parts.append(np.load(p))
    A = np.vstack(parts).astype(np.int64) if parts else np.zeros((0, dim), dtype=np.int64)
    return np.unique(A, axis=0)
