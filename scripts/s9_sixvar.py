#!/usr/bin/env python3
"""s9: substitute HUMAN-TRANSCRIBED 6-variable linear-rank inequalities and re-judge
the stubborn rays (#17/#18/#19).  STRICT RULE: the CSV must be transcribed from the
literature (DFZ 0910.0284 six-variable section / follow-ups).  Do NOT let any tool
invent inequalities; an invalid row silently poisons every rank downstream.

CSV format (data/sixvar_template.csv shows the header):
    ineq_id,subset,coeff        subset over letters ABCDEF, e.g. "ACD"
Each ineq must be balanced per variable (validity screen refuses otherwise unless
--no-balance-check).  Pipeline: 63-dim rows -> 720 bijections onto the 6 purified
subsystems -> purity projection -> unique -> save as --extra-rows for s3."""
import argparse, csv, sys, os, itertools
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import numpy as np
from epr1kit import core

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument('csv')
ap.add_argument('--out', default='sixvar_rows.npy')
ap.add_argument('--no-balance-check', action='store_true')
a = ap.parse_args()
LET = {c: i for i, c in enumerate('ABCDEF')}
ineqs = {}
lines = (ln for ln in open(a.csv) if not ln.lstrip().startswith('#'))
for row in csv.DictReader(lines):
    m = 0
    for ch in row['subset'].strip():
        m |= 1 << LET[ch.upper()]
    ineqs.setdefault(row['ineq_id'], {})
    ineqs[row['ineq_id']][m] = ineqs[row['ineq_id']].get(m, 0) + int(row['coeff'])
P = core.purity_index()
rows = set()
for iid, terms in ineqs.items():
    if not a.no_balance_check:
        for v in range(6):
            s = sum(c for m, c in terms.items() if (m >> v) & 1)
            if s != 0:
                raise SystemExit(f'{iid}: unbalanced in variable {"ABCDEF"[v]} (sum {s}) '
                                 '-- transcription error?')
    for p in itertools.permutations(range(6)):
        acc = [0] * 31
        for m, c in terms.items():
            gm = 0
            for i in range(6):
                if (m >> i) & 1:
                    gm |= 1 << p[i]
            idx = P[gm] if 0 < gm < 63 else -1
            if idx >= 0:
                acc[idx] += c
        if any(acc):
            rows.add(tuple(acc))
A = np.array(sorted(rows), dtype=np.int64)
np.save(a.out, A)
print(f'{len(ineqs)} transcribed inequalities -> {A.shape[0]} substituted rows -> {a.out}')
print('next: scripts/s2_judge.py', a.out, ' (0 violations required), then '
      's3_rank19.py <pureH> --extra-rows', a.out)
