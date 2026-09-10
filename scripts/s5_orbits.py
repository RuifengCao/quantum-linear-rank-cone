#!/usr/bin/env python3
"""s5: reduce enumerated extreme rays to S5 orbits; reconcile against the
literature count of 162 CLR_5 extreme-ray orbits.  Input auto-detected:
Normaliz .out, lrs .lrs.out (V-rep; '*' comments and the origin vertex are
skipped), or a plain whitespace matrix via --raw.  Rays are gcd-normalized
to primitive integer vectors before orbit reduction."""
import argparse, sys, os, re
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import numpy as np
from epr1kit import core

def parse_lrs_out(path):
    rows = []
    inside = False
    for ln in open(path):
        t = ln.strip()
        if t == 'begin':
            inside = True; continue
        if t == 'end':
            break
        if not inside or not t or t.startswith('*'):
            continue
        parts = t.split()
        if len(parts) != 32 or any(not p.lstrip('-').isdigit() for p in parts):
            continue
        if parts[0] != '0':          # '1 0 ... 0' = origin vertex; rays lead with 0
            continue
        rows.append(list(map(int, parts[1:])))
    if not rows:
        raise SystemExit('no rays parsed from lrs output ' + path)
    return np.array(rows, dtype=np.int64)

def parse_normaliz_out(path):
    txt = open(path).read()
    m = re.search(r'(\d+)\s+extreme rays:\s*\n((?:\s*[-\d ]+\n)+)', txt)
    if not m:
        raise SystemExit('could not find extreme rays block in ' + path)
    n = int(m.group(1))
    rows = [list(map(int, ln.split())) for ln in m.group(2).strip().splitlines()]
    A = np.array(rows, dtype=np.int64)
    assert A.shape == (n, 31), A.shape
    return A

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument('src', help='normaliz .out (or matrix file with --raw)')
ap.add_argument('--raw', action='store_true')
ap.add_argument('--expect', type=int, default=162)
ap.add_argument('--out', default='clr5_extreme_orbit_reps.npy')
a = ap.parse_args()
if a.raw:
    R = np.loadtxt(a.src, dtype=np.int64).reshape(-1, 31)
else:
    head = open(a.src).read(4000)
    R = parse_lrs_out(a.src) if ('begin' in head and 'lrs' in head) else parse_normaliz_out(a.src)
g = np.gcd.reduce(np.abs(R), axis=1); g[g == 0] = 1
R = np.unique(R // g[:, None], axis=0)          # primitive + dedupe
if R.min() < 0 or R.max() > 255:
    raise SystemExit(f'ray values outside [0,255] (min {R.min()}, max {R.max()}); '
                     'inspect the source file -- canonicalization assumes small nonneg ints')
perms = core.perms31_s5()
canon, reps = core.class_reps(R, perms, offset=0)   # rays are nonnegative small ints
orb = sorted(reps.values())
np.save(a.out, R[orb])
print(f'extreme rays {R.shape[0]}  ->  S5 orbits {len(orb)}  (literature: {a.expect})')
print('MATCH' if len(orb) == a.expect else 'MISMATCH -- see README s12 note before trusting CLR_H')
