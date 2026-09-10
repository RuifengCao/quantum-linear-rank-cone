#!/usr/bin/env python3
"""s12: re-verify the R3 cone reconciliation against the bundled literature
reference (data/dfz_ref28.csv: DFZ 0910.0284 eqs (1)-(24) + the 4 Ingleton
forms (36)-(39) from the completeness section = 28 non-Shannon classes).

SETTLED in R3 (2026-08-12, reproduced independently on the user's server):
  direction 1: the only reference class NOT implied by cone(M_LR+Shannon)
               is ing39 (the overlapping-slot Ingleton form);
  direction 2: all 27 psitip classes ARE implied by the reference cone.
  VERDICT: PSITIP-INCOMPLETE.  Remedy bundled: CLR_H_fixed (s4 default),
  templates28 / pure28 (headline 18/19, only #19 non-extreme, exact rank 29).
This script recomputes both directions from the CSV; expect the same output.
Historical background (27-class irredundancy, certificates) lives in
data/sep_certs_27.npz and the R3 report."""
import argparse, csv, sys, os, itertools
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import numpy as np
from scipy.optimize import linprog
from epr1kit import core

def all_perm_instances(rows31):
    perms = core.perms31_s5()
    out = set()
    for r in np.asarray(rows31, dtype=np.int64):
        for k in range(perms.shape[0]):
            out.add(tuple(int(x) for x in r[perms[k]]))
    return np.array(sorted(out), dtype=np.int64)

def implied(t, R):
    res = linprog(np.zeros(R.shape[0]), A_eq=R.T.astype(float), b_eq=np.asarray(t, float),
                  bounds=[(0, None)] * R.shape[0], method='highs')
    return res.status == 0

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument('csv', nargs='?', default=os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data', 'dfz_ref28.csv'),
                help='reference CSV over ABCDE; default: bundled dfz_ref28.csv (24 DFZ + 4 Ingleton forms)')
a = ap.parse_args()
LET = {c: i for i, c in enumerate('ABCDE')}
ineqs = {}
lines = (ln for ln in open(a.csv) if not ln.lstrip().startswith('#'))
for row in csv.DictReader(lines):
    m = 0
    for ch in row['subset'].strip():
        m |= 1 << LET[ch.upper()]
    ineqs.setdefault(row['ineq_id'], [0] * 31)
    ineqs[row['ineq_id']][m - 1] += int(row['coeff'])
DFZ = np.array([v for v in ineqs.values()], dtype=np.int64)
print(f'transcribed {DFZ.shape[0]} inequalities (expect 28 = 24 DFZ + 4 Ingleton forms)')
SUB, WM, MONO = core.elementals5()
SH = np.unique(np.vstack([SUB, MONO]), axis=0)
M = core.load('M_LR'); T = core.load('templates27')
DFZ_full = all_perm_instances(DFZ)
missing_from_psitip = [iid for iid, t in zip(ineqs, DFZ)
                       if not implied(t, np.vstack([M, SH]))]
# R3 in-container result (2026-08-12): missing_from_psitip == ['ing39'], missing_from_dfz == []
# => VERDICT PSITIP-INCOMPLETE; remedy applied: CLR_H_fixed = CLR_H + 30 instances of ing39,
#    templates28 = templates27 + ing39, pure28 rebuilt (18/19, only #19 non-extreme, rank 29).
missing_from_dfz = [ti for ti, t in enumerate(T)
                    if not implied(t, np.vstack([DFZ_full, SH]))]
print('DFZ classes NOT implied by psitip cone:', missing_from_psitip or 'none')
print('psitip classes NOT implied by DFZ cone:', missing_from_dfz or 'none')
if not missing_from_psitip and not missing_from_dfz:
    print('VERDICT: cones EQUAL -- CLR_H is CLR_5; the 27 vs 25 gap is a class-'
          'bookkeeping artefact; s4-s6 may proceed unchanged.')
elif missing_from_psitip:
    print('VERDICT: PSITIP-INCOMPLETE -- CLR_H cuts a strictly larger cone than CLR_5.'
          ' Remedy ALREADY BUNDLED in this kit: s4 defaults to CLR_H_fixed'
          ' (= CLR_H + ing39 orbit); QLR side fixed via templates28/pure28.'
          ' Nothing to append manually.')
else:
    print('VERDICT: EXTRACTION-SUSPECT -- psitip rows not implied by the complete list;'
          ' re-audit transcription AND the D1 inclusion-exclusion extraction.')
