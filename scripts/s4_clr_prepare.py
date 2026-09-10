#!/usr/bin/env python3
"""s4: write the CLR_5 H-rep as Normaliz (.in) and lrs (.ine) inputs for
extreme-ray enumeration on the server.  R3 update: s12 verdict is in
(PSITIP-INCOMPLETE, missing Ingleton(39)); default input is CLR_H_fixed."""
import argparse, sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import numpy as np
from epr1kit import core

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument('--H', default=None, help='npy H-rep; default bundled CLR_H_fixed (1,905; R3 corrected)')
ap.add_argument('--prefix', default='clr5')
a = ap.parse_args()
H = np.load(a.H).astype(np.int64) if a.H else core.load('CLR_H_fixed')
core.write_normaliz(H, a.prefix + '.in')
core.write_lrs(H, a.prefix + '.ine')
print(f'wrote {a.prefix}.in / {a.prefix}.ine  ({H.shape[0]} inequalities, dim 31)')
