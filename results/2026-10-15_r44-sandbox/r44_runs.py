#!/usr/bin/env python3
"""R44 sandbox runs (Theorem 14): exhaustive families of codes over Z_{2^a} (cc/s38_minform_zn.c, lifts mode) and
annealing runs.  Checks, not proofs.

  python3 r44_runs.py lifts3        # 3 x 7 codes: all 0/1 codes over Z_4, Z_8, Z_16; lifts of the seven-point matrix
  python3 r44_runs.py lifts4        # all 4 x 7 0/1 codes over Z_4 (2^28 codes, ~30 min on one core)
  python3 r44_runs.py anneal        # cc/s38_minform_zn.c on CSS codes over Z_4 and Z_8
"""
import os, sys, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'scripts'))
import s35_two_power_transfer as z

SEVEN = z.SEVEN


def lifts3():
    rng = np.random.default_rng(4403)
    fam = [(2, np.zeros((3, 7), dtype=np.int64), 0, 'all 0/1 codes'),
           (2, SEVEN, 1, 'seven-point + 2P'),
           (3, np.zeros((3, 7), dtype=np.int64), 0, 'all 0/1 codes'),
           (3, SEVEN, 1, 'seven-point + 2P'),
           (3, SEVEN, 2, 'seven-point + 4P'),
           (3, 2 * SEVEN, 2, 'Z_4 seven-point in Z_8 + 4P'),
           (4, np.zeros((3, 7), dtype=np.int64), 0, 'all 0/1 codes'),
           (4, SEVEN, 3, 'seven-point + 8P')]
    for t in range(8):
        fam.append((2, rng.integers(0, 2, size=(3, 7)), 1, f'random 0/1 base #{t} + 2P'))
    ok = True
    for a, B, s, name in fam:
        t0 = time.time()
        print(f'# Z_{1 << a}, {name}', flush=True)
        r = z.lifts(a, B, s, z.weakened_forms(a))
        first = next(iter(r.values()))
        ok &= first['min'] >= 0 and first['recomputed'] == first['min'] and not first['chain_failed']
        ok &= all(v['min'] == v['recomputed'] for v in r.values())
        print(f'  ({time.time() - t0:.0f} s)', flush=True)
    print('ALL OK' if ok else 'FAILURE')


def lifts4():
    t0 = time.time()
    r = z.lifts(2, np.zeros((4, 7), dtype=np.int64), 0, z.weakened_forms(2))
    first = next(iter(r.values()))
    print(f'  ({time.time() - t0:.0f} s)')
    print('ALL OK' if first['min'] >= 0 and first['recomputed'] == first['min'] and not first['chain_failed'] else 'FAILURE')


def anneal():
    runs = [dict(cf=[(2, 3, [1] * 7), (3, 3, [1] * 7)], kw={}),
            dict(cf=[(2, 4, [1] * 7), (3, 4, [1] * 7)], kw=dict(mu=1.0, singles=True)),
            dict(cf=[(2, 4, [1] * 7), (3, 4, [1] * 7)], kw=dict(mu=0.3)),
            dict(cf=[(2, 6, [2] * 7), (3, 6, [2] * 7)], kw=dict(mu=0.3)),
            dict(cf=[(2, 5, [2, 1, 2, 1, 2, 1, 2]), (2, 6, [2, 2, 2, 1, 1, 1, 2])], kw=dict(mu=1.0, singles=True))]
    ok = True
    for i, run in enumerate(runs):
        print(f'# run {i}: {run["cf"]} {run["kw"]}', flush=True)
        r = z.anneal(run['cf'], seconds=120, seed=4401 + 100 * i, **run['kw'])
        ok &= all(x['min_form'] >= 0 and x['form_recomputed'] == x['min_form'] and not x['chain_failed'] for x in r)
    print('ALL OK' if ok else 'FAILURE')


if __name__ == '__main__':
    {'lifts3': lifts3, 'lifts4': lifts4, 'anneal': anneal}[sys.argv[1]]()
