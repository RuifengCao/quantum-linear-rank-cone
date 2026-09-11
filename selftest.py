#!/usr/bin/env python3
"""EPR-1 R3 selftest.  Gates reproduce every headline number of R2 from scratch.

fast (~25 s, default): G0 parallel shard/resume smoke, G1 elementals, G2 CLR sha,
                       G2b CLR_H_fixed sha, G4 pure sha, G6 graph-state sha,
                       G7 judge (pure: 0 violations; CLR: exactly 5 rows),
                       G9 reference family (28 classes, sha, ing39 orbit, targets1),
                       G10 pure28 sha + judge 0 + rank 18/19 (only #19, exact rank 29).
full (~2-6 min):       + G3 v3 build/judge/rank (10/19), G5 pure2 sha,
                       G8 pure(27) rank 16/19 with #17/#18/#19 at exact rank 29 (historical).
Run on the server FIRST; a red gate means the environment or a recipe drifted.
"""
import argparse, json, os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from epr1kit import core

def _g0_rows(i):
    return {np.full(31, i, dtype=np.int8).tobytes(),
            np.full(31, i + 10, dtype=np.int8).tobytes()}

MAN = json.load(open(os.path.join(core.DATA, 'manifest_shas.json')))
FAIL = []

def gate(name, ok, detail=''):
    print(f'[{name}] {"PASS" if ok else "FAIL"}  {detail}')
    if not ok:
        FAIL.append(name)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--full', action='store_true')
    ap.add_argument('--workers', type=int, default=1)
    a = ap.parse_args()
    t0 = time.time()

    import tempfile, shutil
    from epr1kit.parallel import pmap_shards
    tmpd = tempfile.mkdtemp(prefix='epr1_g0_')
    try:
        fn = _g0_rows
        A1 = pmap_shards(fn, [0, 1, 2, 3], tmpd, 'g0', workers=2)
        os.remove(os.path.join(tmpd, 'g0_00002.npy'))          # simulate a lost shard
        A2 = pmap_shards(fn, [0, 1, 2, 3], tmpd, 'g0', workers=2)
        ok0 = np.array_equal(A1, A2) and A1.shape == (8, 31)
    finally:
        shutil.rmtree(tmpd, ignore_errors=True)
    gate('G0 parallel', ok0, f'shard+resume merge {A1.shape} (expect (8, 31), stable)')

    SUB, WM, MONO = core.elementals5()
    gate('G1 elementals', (SUB.shape[0], WM.shape[0], MONO.shape[0]) == (80, 195, 5),
         f'{SUB.shape[0]}/{WM.shape[0]}/{MONO.shape[0]} (expect 80/195/5)')

    CLR = core.build_clr()
    gate('G2 CLR', core.sha_rows(CLR) == MAN['CLR_H']['sha256'],
         f'{CLR.shape[0]} rows (expect 1875, sha match)')

    CLRF = np.unique(np.vstack([CLR, core.load('ING39_orbit')]), axis=0)
    gate('G2b CLR_fixed', core.sha_rows(CLRF) == MAN['CLR_H_fixed']['sha256'],
         f'{CLRF.shape[0]} rows (expect 1905, sha match)')

    PURE = core.build_qlr('pure', workers=a.workers)
    gate('G4 pure', core.sha_rows(PURE) == MAN['QLR_H_pure']['sha256'],
         f'{PURE.shape[0]} rows (expect 23265, sha match)')

    GS = core.graphstate_vecs6()
    gate('G6 graphstates', core.sha_rows(GS) == MAN['graphstate_vecs']['sha256'],
         f'{GS.shape[0]} vectors (expect 760, sha match)')

    _, cp = core.judge(PURE, GS)
    _, cc = core.judge(CLR, GS)
    gate('G7 judge', int((cp > 0).sum()) == 0 and int((cc > 0).sum()) == 5,
         f'pure violated rows {(cp>0).sum()} (expect 0); CLR violated rows {(cc>0).sum()} (expect 5)')

    REF = np.load(os.path.join(core.DATA, 'REF28.npy'))
    import itertools as _it
    perms5 = core.perms31_s5()
    canon = core.class_reps(REF, perms5)[0]
    bal = all(sum(int(r[m - 1]) for m in range(1, 32) if (m >> v) & 1) == 0
              for r in REF for v in range(5))
    orb39 = np.unique(np.stack([REF[27][perms5[k]] for k in range(120)]), axis=0)
    t1 = np.array_equal(core.load('targets1'), core.load('hec5_rays_maskorder')[18:19])
    gate('G9 ref28', REF.shape[0] == 28 and len(set(canon)) == 28 and bal
         and core.sha_rows(REF) == MAN['REF28']['sha256']
         and np.array_equal(orb39, np.unique(core.load('ING39_orbit'), axis=0)) and t1,
         f'{REF.shape[0]} ineqs, {len(set(canon))} S5 classes, balanced={bal}, '
         f'sha match, ing39 orbit consistent, targets1 == ray #19')

    P28 = core.build_qlr('pure28', workers=a.workers)
    _, c28 = core.judge(P28, GS)
    gate('G10 pure28', core.sha_rows(P28) == MAN['QLR_H_pure28']['sha256']
         and int((c28 > 0).sum()) == 0,
         f'{P28.shape[0]} rows (expect 23355, sha match); violations {(c28>0).sum()} (0)')
    rep28 = core.tight_rank_report(P28, core.load('hec5_rays_maskorder'))
    e28 = sum(1 for *_x, e in rep28 if e)
    non28 = {ri + 1: rk for ri, _n, rk, e in rep28 if not e}
    gate('G10b rank19/28', e28 == 18 and non28 == {19: 29},
         f'extreme {e28}/19 (expect 18); non-extreme {non28} (expect only #19 at 29)')

    import tempfile as _tf
    canned = ('*lrs\ncone\nV-representation\nbegin\n***** 32 rational\n'
              ' 1  ' + ' 0' * 31 + ' \n 0  ' + ' 2' * 31 + ' \n*cobasis junk\n 0  ' + ' 1' * 31 + ' \nend\n')
    with _tf.NamedTemporaryFile('w', suffix='.lrs.out', delete=False) as fh:
        fh.write(canned); cp = fh.name
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'scripts'))
    import importlib.util as _iu
    spec = _iu.spec_from_file_location('_s5', os.path.join(os.path.dirname(os.path.abspath(__file__)), 'scripts', 's5_orbits.py'))
    ok11 = True
    try:
        import re as _re
        src = open(spec.origin).read()
        ns = {'__file__': spec.origin, '__name__': '_s5'}
        exec(src[:src.index('ap = argparse')], ns)
        A = ns['parse_lrs_out'](cp)
        ok11 = A.shape == (2, 31) and set(A[0]) == {2} and set(A[1]) == {1}
    except Exception as e:
        ok11 = False
    finally:
        os.remove(cp)
    gate('G11 lrs-parse', ok11, f'canned V-rep -> {A.shape if ok11 else "parse error"} (expect (2, 31); vertex+comments skipped)')

    RY = core.load('clr5_rays7943').astype(np.int64)
    RP = core.load('clr5_orbit_reps162')
    okv = bool(((core.load('CLR_H_fixed') @ RY.T) >= 0).all())
    cb12, reps12 = core.class_reps(RY, perms5, offset=0)
    gate('G12 clr5 rays', RY.shape == (7943, 31) and okv and len(reps12) == 162
         and RP.shape == (162, 31),
         f'{RY.shape[0]} rays valid={okv}, {len(reps12)} orbits (expect 7943/True/162)')

    import tempfile as _tf2
    RYt = core.load('clr5_rays7943').astype(np.int64)
    with _tf2.NamedTemporaryFile('w', suffix='.rays5', delete=False) as fh2:
        fh2.write('# DFZ rays (synthetic header)\n7943 31 integer\n')
        for row in RYt[:64]:
            fh2.write(' '.join(map(str, row)) + ' \n')
        fh2.write('end junk trailing line\n'); cp2 = fh2.name
    try:
        rows2, skip2 = [], 0
        for ln in open(cp2):
            toks = ln.replace(',', ' ').split()
            ints, ok2 = [], True
            for t in toks:
                tt = t.lstrip('+-')
                if tt.isdigit(): ints.append(int(t))
                else: ok2 = False; break
            if not ok2 or len(ints) not in (31, 32): skip2 += 1; continue
            if len(ints) == 32:
                if ints[0] not in (0, 1): skip2 += 1; continue
                ints = ints[1:]
            rows2.append(ints)
        ok13 = len(rows2) == 64 and skip2 == 3
    finally:
        os.remove(cp2)
    gate('G13 rays5-parse', ok13, f'{len(rows2)} rows kept, {skip2} junk skipped (expect 64/3)')

    C19 = np.load(os.path.join(core.DATA, 'cert19_exact.npz'))
    q1c, q2c, rayc = C19['q1'].astype(np.int64), C19['q2'].astype(np.int64), C19['ray'].astype(np.int64)
    al = int(C19['alpha'][0]) / int(C19['alpha'][1]); be = int(C19['beta'][0]) / int(C19['beta'][1])
    P28c = core.build_qlr('pure28', workers=a.workers) if 'P28' not in dir() else P28
    okc = (np.array_equal(rayc, core.load('targets1')[0])
           and np.array_equal(q1c, core.load('hec5_rays_maskorder')[3])
           and np.array_equal(1 * q1c + 2 * q2c, rayc)
           and core.rank_exact_gram(P28c[(P28c @ q1c) == 0]) == 30
           and core.rank_exact_gram(P28c[(P28c @ q2c) == 0]) == 30)
    gate('G14 cert19', okc, 'r19 == 1*q1 + 2*q2 exact; q1 == HEC#4; both edges rank 30')

    import json as _js
    T3 = core.load('shc_table3_rays').astype(np.int64)
    M3 = _js.load(open(os.path.join(core.DATA, 'shc_table3_map.json')))['paper_to_ours']
    ct3, _ = core.class_reps(T3, core.perms31_s6())
    co3, _ = core.class_reps(core.load('hec5_rays_maskorder'), core.perms31_s6())
    bij = all(ct3[int(p) - 1] == co3[o - 1] for p, o in M3.items())
    gate('G15 shc-map', bij and len(set(M3.values())) == 19 and M3['18'] == 19,
         f'Table3<->ours S6 bijection {bij}; paper #18 == our #19')

    Q59 = core.load('qlr59_reps').astype(np.int64)
    c59, _ = core.class_reps(Q59, core.perms31_s6())
    A3 = core.load('a3_hyp_ineq44').astype(np.int64).reshape(31)
    p6g = core.perms31_s6()
    A3orb = np.unique(np.stack([A3[p6g[k]] for k in range(720)]), axis=0)
    okv59 = bool(((P28c @ Q59.T) >= 0).all())
    g15min = int((A3orb @ Q59[-1]).min())
    clrmin = int((A3orb @ core.load('clr5_rays7943').astype(np.int64).T).min())
    gate('G16 qlr59+a3', Q59.shape == (59, 31) and len(set(c59)) == 59 and okv59
         and core.sha_rows(Q59) == MAN['qlr59_reps']['sha256']
         and g15min == -1 and clrmin == 0,
         f'59 reps S6-distinct valid; ineq(4.4): G15 min {g15min} (expect -1), CLR rays min {clrmin} (expect 0)')

    HFq = core.load('qlr5_H_facets10860').astype(np.int64)
    RQ = core.load('qlr5_orbits_partial').astype(np.int64)
    rngs = np.random.default_rng(1717); sub = RQ[rngs.choice(RQ.shape[0], 400, replace=False)]
    cRQ, _ = core.class_reps(sub, core.perms31_s6(), offset=0, width='u16')
    okq = bool(((HFq @ RQ.T) >= 0).all()) and len(set(cRQ)) == sub.shape[0] \
        and MAN['qlr5_orbits_partial'].get('s6_distinct') is True
    rng17 = np.random.default_rng(17); samp = rng17.choice(RQ.shape[0], 40, replace=False)
    okr = all(core.rank_mod_p(HFq[(HFq @ RQ[i]) == 0]) == 30 for i in samp)
    gate('G17 qlr-partial', okq and okr and HFq.shape == (10860, 31)
         and core.sha_rows_wide(RQ) == MAN['qlr5_orbits_partial']['sha256_int16'],
         f'{RQ.shape[0]} certified QLR5 orbits (valid; 400-sample S6-distinct; 40-sample rank 30; full checks recorded in manifest); H_facets 10860')

    W35 = np.load(os.path.join(core.DATA, 'qlr5_new_witnessed35.npz'))
    p6w = core.perms31_s6()
    okw = all(np.array_equal(W35['rep'][t].astype(np.int64)[p6w[int(W35['perm'][t])]] * int(W35['mult'][t]),
                             W35['witness'][t].astype(np.int64)) for t in range(W35['rep'].shape[0]))
    okw = okw and bool(((HFq @ W35['rep'].astype(np.int64).T) >= 0).all()) \
        and all(core.rank_mod_p(HFq[(HFq @ W35['rep'][t].astype(np.int64)) == 0]) == 30 for t in range(W35['rep'].shape[0]))
    SMo = core.load('qlr5_small_orbits').astype(np.int64)
    ksm, _ = core.class_reps(SMo, p6w, offset=0, width='u16')
    oks = SMo.shape[0] == 15183 and len(set(ksm)) == SMo.shape[0] and bool(((HFq @ SMo.T) >= 0).all()) \
        and core.sha_rows(SMo) == MAN['qlr5_small_orbits']['sha256']
    gate('G18 witness35+small', okw and oks and int((W35['source'] != 'F3cond').sum()) == 33,
         f"35 witnessed new orbits re-verified exactly ({int((W35['source'] != 'F3cond').sum())} unconditional); small catalogue {SMo.shape[0]} valid, S6-distinct")

    from epr1kit import chordal as _ch
    hecr = core.load('hec5_rays_maskorder').astype(np.int64)
    res19 = [_ch.analyse(v) for v in hecr]
    n_ch = sum(1 for r in res19 if r['chordal']); n_ci = sum(1 for r in res19 if r['chordal'] and r['irreducible'])
    n_ver = sum(1 for r in res19 if r['verified'])
    gate('G19 chordal', n_ch == 9 and n_ci == 6 and n_ver == 6 and all(r['sa_ssa'] for r in res19),
         f'HEC5 rays: chordal {n_ch} (expect 9), chordal&irreducible {n_ci} (6), simple trees verified by min-cut {n_ver}/6')

    if a.full:
        V3 = core.build_qlr('v3')
        _, cv = core.judge(V3, GS)
        rep3 = core.tight_rank_report(V3, core.load('hec5_rays_maskorder'))
        e3 = sum(1 for *_x, e in rep3 if e)
        gate('G3 v3', V3.shape[0] == 2065 and int((cv > 0).sum()) == 0 and e3 == 10,
             f'{V3.shape[0]} rows (expect 2065), violations {(cv>0).sum()} (0), extreme {e3}/19 (10)')

        P2 = core.build_qlr('pure2', workers=a.workers)
        gate('G5 pure2', core.sha_rows(P2) == MAN['QLR_H_pure2']['sha256'],
             f'{P2.shape[0]} rows (expect 66577, sha match)')

        rep = core.tight_rank_report(PURE, core.load('hec5_rays_maskorder'))
        ext = sum(1 for *_x, e in rep if e)
        non = {ri + 1: rk for ri, _n, rk, e in rep if not e}
        gate('G8 rank19', ext == 16 and non == {17: 29, 18: 29, 19: 29},
             f'extreme {ext}/19 (expect 16); non-extreme {non} (expect 17/18/19 at 29)')

    print(f'--- {"ALL PASS" if not FAIL else "FAILED: " + ",".join(FAIL)}  '
          f'({time.time()-t0:.0f}s)')
    sys.exit(1 if FAIL else 0)

if __name__ == '__main__':
    main()
