#!/usr/bin/env python3
"""s1: build QLR H-representations.
variants: pure28 (23,355 rows, 28 templates incl. Ingleton(39)) -- CURRENT cone;
          pure (23,265, 27 templates) and v3 (2,065) -- historical; pure2 (66,577) -- allow-empty."""
import argparse, sys, os, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import numpy as np
from epr1kit import core

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument('--variant', choices=['v3', 'pure', 'pure2', 'pure28'], required=True)
ap.add_argument('--workers', type=int, default=os.cpu_count())
ap.add_argument('--out', default=None)
a = ap.parse_args()
t0 = time.time()
A = core.build_qlr(a.variant, workers=a.workers)
out = a.out or f'QLR_H_{a.variant}.npy'
np.save(out, A)
print(f'{a.variant}: {A.shape[0]} rows  sha={core.sha_rows(A)[:16]}  {time.time()-t0:.0f}s -> {out}')
