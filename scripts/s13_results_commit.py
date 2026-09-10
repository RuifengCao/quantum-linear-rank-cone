#!/usr/bin/env python3
"""s13: stage result files into results/<date>_<tag>/ with a manifest.

Usage:
  python3 scripts/s13_results_commit.py --tag a1-batch2 qlr_adj.npy qlr_adj.npy.growth.csv s4c.log
  python3 scripts/s13_results_commit.py --tag a1-batch2 --note "pools in release v0.1" ...

Copies the files, writes manifest.json (size + sha256 per file; for 2-D
integer .npy files also the canonical row-family sha used by data/manifest_shas.json)
and prints the suggested commit message.  Refuses files above --max-mb (default 50).
"""
import argparse, datetime, hashlib, json, os, shutil, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from epr1kit import core

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument('files', nargs='+')
ap.add_argument('--tag', required=True, help='short label, e.g. a1-batch2')
ap.add_argument('--note', default='')
ap.add_argument('--root', default=os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results'))
ap.add_argument('--max-mb', type=float, default=50)
a = ap.parse_args()

day = datetime.date.today().isoformat()
dest = os.path.join(a.root, f'{day}_{a.tag}')
os.makedirs(dest, exist_ok=True)
man = {'date': day, 'tag': a.tag, 'note': a.note, 'files': {}}
for f in a.files:
    if not os.path.isfile(f):
        raise SystemExit(f'not a file: {f}')
    mb = os.path.getsize(f) / 1e6
    src = f
    entry = {}
    if mb > a.max_mb and f.endswith('.npy'):
        # try a lossless downcast of integer arrays (e.g. int64 pools -> int8/int16)
        try:
            arr = np.load(f)
            if np.issubdtype(arr.dtype, np.integer):
                for dt in (np.int8, np.int16, np.int32):
                    if arr.min() >= np.iinfo(dt).min and arr.max() <= np.iinfo(dt).max:
                        small = os.path.join(dest, os.path.basename(f))
                        np.save(small, arr.astype(dt)); src = small
                        entry['downcast_from'] = str(arr.dtype); break
        except Exception as e:
            entry['npy_note'] = f'downcast failed: {e}'
        mb = os.path.getsize(src) / 1e6
    if mb > a.max_mb:
        print(f'SKIPPED {f}: {mb:.1f} MB > {a.max_mb} MB -- attach it to a GitHub Release instead')
        man['files'][os.path.basename(f)] = {'skipped': f'{mb:.1f} MB > limit; attach to a Release'}
        continue
    if src == f:
        shutil.copy2(f, dest)
    entry.update({'bytes': os.path.getsize(src), 'sha256': hashlib.sha256(open(src, 'rb').read()).hexdigest()})
    if f.endswith('.npz'):
        try:
            Z = np.load(src)
            if 'reps' in Z:
                entry.update({'orbits': int(Z['reps'].shape[0]), 'expanded': int(Z['done'].sum()),
                              'skipped': int(Z['skipped'].sum()), 'state_format': 'npz'})
        except Exception as e:
            entry['npz_note'] = f'not summarised: {e}'
    if f.endswith('.npy'):
        try:
            arr = np.load(src, allow_pickle=True)
            if isinstance(arr, np.ndarray) and arr.ndim == 2 and np.issubdtype(arr.dtype, np.integer):
                A = arr.astype(np.int64)
                entry['rows'] = int(A.shape[0])
                entry['sha_rows'] = core.sha_rows(A) if np.abs(A).max() <= 127 else core.sha_rows_wide(A)
            elif isinstance(arr, np.ndarray) and arr.dtype == object:
                S = arr.item()
                if isinstance(S, dict) and 'reps' in S:
                    if S.get('format') == 2:
                        entry['orbits'] = int(len(S['reps'])); entry['expanded'] = int(np.sum(S['done']))
                        entry['skipped'] = int(np.sum(S['skipped'])); entry['state_format'] = 2
                    else:
                        entry['orbits'] = len(S['reps']); entry['expanded'] = len(S.get('done', []))
                        entry['skipped'] = len(S.get('skipped', [])); entry['state_format'] = 1
        except Exception as e:
            entry['npy_note'] = f'not summarised: {e}'
    man['files'][os.path.basename(f)] = entry
json.dump(man, open(os.path.join(dest, 'manifest.json'), 'w'), indent=1)
print(f'staged {len(a.files)} file(s) into {dest}')
print(f'suggested commit message:  results: {a.tag} ({day})')
