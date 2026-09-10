#!/usr/bin/env python3
"""s11: gather everything to bring back for analysis: named outputs + manifest with
row counts and canonical sha256, packed as one tar.gz."""
import argparse, sys, os, json, tarfile, glob
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import numpy as np
from epr1kit import core

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument('paths', nargs='+', help='files and/or globs to collect')
ap.add_argument('--out', default='epr1_server_results.tar.gz')
a = ap.parse_args()
files = sorted({p for pat in a.paths for p in glob.glob(pat)})
man = {}
for p in files:
    e = {'bytes': os.path.getsize(p)}
    if p.endswith('.npy'):
        arr = np.load(p, allow_pickle=False)
        e['shape'] = list(arr.shape)
        if arr.ndim == 2 and arr.shape[1] == 31:
            e['sha256_canonical'] = core.sha_rows(arr)
    man[p] = e
json.dump(man, open('collect_manifest.json', 'w'), indent=1)
with tarfile.open(a.out, 'w:gz') as t:
    for p in files + ['collect_manifest.json']:
        t.add(p)
print(f'{len(files)} files -> {a.out}')
