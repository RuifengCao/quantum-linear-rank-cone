#!/usr/bin/env python3
"""s5b: independent cross-certification of the 162/7,943 enumeration against
DFZ's own published ray list (rays5 from code.ucsd.edu/zeger/linrank/,
manual browser download -- the site blocks robots).

Robust parser: keeps any line containing exactly 31 integers (or 32 with a
leading 0/1 homogenization column, which is dropped); everything else is
skipped.  Rays are gcd-normalized to primitive vectors, deduplicated, then
compared as SETS against the bundled data/clr5_rays7943.

MATCH criterion: identical canonical sha both ways (=> ray sets equal),
plus orbit counts.  A MATCH makes our enumeration and DFZ's mutually
independent certificates of the same 162-orbit answer.
Usage:  python3 scripts/s5b_diff_rays.py <path-to-rays5>
"""
import argparse, os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from epr1kit import core

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument('src', help='downloaded rays5 file (any reasonable text format)')
a = ap.parse_args()

rows = []
skipped = 0
for ln in open(a.src, errors='replace'):
    toks = ln.replace(',', ' ').split()
    ints = []
    ok = True
    for t in toks:
        tt = t.lstrip('+-')
        if tt.isdigit():
            ints.append(int(t))
        else:
            ok = False
            break
    if not ok or len(ints) not in (31, 32):
        skipped += 1
        continue
    if len(ints) == 32:
        if ints[0] not in (0, 1):
            skipped += 1
            continue
        ints = ints[1:]
    rows.append(ints)
if not rows:
    raise SystemExit('no 31-column integer rows found -- send the first 20 lines of the file for a parser fix')
R = np.array(rows, dtype=np.int64)
g = np.gcd.reduce(np.abs(R), axis=1); g[g == 0] = 1
R = np.unique(R // g[:, None], axis=0)
print(f'parsed {len(rows)} rows ({skipped} skipped) -> {R.shape[0]} unique primitive rays')

OURS = core.load('clr5_rays7943').astype(np.int64)
H = core.load('CLR_H_fixed')
val = int(((H @ R.T) >= 0).all(axis=0).sum())
print(f'validity vs CLR_H_fixed: {val}/{R.shape[0]}')

sha_theirs, sha_ours = core.sha_rows(R), core.sha_rows(OURS)
A = {r.tobytes() for r in R.astype(np.int16)}
B = {r.tobytes() for r in OURS.astype(np.int16)}
only_t, only_o = len(A - B), len(B - A)
perms = core.perms31_s5()
orb_t = len(core.class_reps(R, perms, offset=0)[1])
orb_o = len(core.class_reps(OURS, perms, offset=0)[1])
print(f'orbits: theirs {orb_t}, ours {orb_o};  set diff: theirs-only {only_t}, ours-only {only_o}')
if sha_theirs == sha_ours:
    print('MATCH: ray sets identical (canonical sha equal) -- independent cross-certification achieved')
else:
    print('MISMATCH: investigate the differing rays before drawing conclusions')
    sys.exit(1)
