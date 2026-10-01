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

def _nthreads():
    """OS threads of this process (Linux; '?' elsewhere) -- diagnostic only."""
    try:
        for ln in open('/proc/self/status'):
            if ln.startswith('Threads:'):
                return int(ln.split()[1])
    except OSError:
        pass
    return '?'

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

    W35 = np.load(os.path.join(core.DATA, 'qlr5_new_witnessed.npz'))
    p6w = core.perms31_s6()
    okw = all(np.array_equal(W35['rep'][t].astype(np.int64)[p6w[int(W35['perm'][t])]] * int(W35['mult'][t]),
                             W35['witness'][t].astype(np.int64)) for t in range(W35['rep'].shape[0]))
    okw = okw and bool(((HFq @ W35['rep'].astype(np.int64).T) >= 0).all()) \
        and all(core.rank_mod_p(HFq[(HFq @ W35['rep'][t].astype(np.int64)) == 0]) == 30 for t in range(W35['rep'].shape[0]))
    SMo = core.load('qlr5_small_orbits').astype(np.int64)
    ksm, _ = core.class_reps(SMo, p6w, offset=0, width='u16')
    oks = SMo.shape[0] == MAN['qlr5_small_orbits']['rows'] and len(set(ksm)) == SMo.shape[0] and bool(((HFq @ SMo.T) >= 0).all()) \
        and core.sha_rows(SMo) == MAN['qlr5_small_orbits']['sha256']
    gate('G18 witness+small', okw and oks and int((W35['source'] != 'F3cond').sum()) == MAN['qlr5_new_witnessed']['unconditional'],
         f"{W35['rep'].shape[0]} witnessed new orbits re-verified exactly ({int((W35['source'] != 'F3cond').sum())} unconditional); small catalogue {SMo.shape[0]} valid, S6-distinct")

    from epr1kit import chordal as _ch
    hecr = core.load('hec5_rays_maskorder').astype(np.int64)
    res19 = [_ch.analyse(v) for v in hecr]
    n_ch = sum(1 for r in res19 if r['chordal']); n_ci = sum(1 for r in res19 if r['chordal'] and r['irreducible'])
    n_ver = sum(1 for r in res19 if r['verified'])
    gate('G19 chordal', n_ch == 9 and n_ci == 6 and n_ver == 6 and all(r['sa_ssa'] for r in res19),
         f'HEC5 rays: chordal {n_ch} (expect 9), chordal&irreducible {n_ci} (6), simple trees verified by min-cut {n_ver}/6')

    # G20 (R27): A2 certificates -- realisations re-verified by an independent GF(2) rank,
    # every exclusion re-derived in exact rational arithmetic, inequality classes checked
    # against the 760 six-qubit graph states.
    import json as _js20
    from epr1kit import stabcert as _sc
    RZ = np.load(os.path.join(core.DATA, 's16_realisations.npz'))
    sub20 = list(range(0, RZ['rep'].shape[0], max(1, RZ['rep'].shape[0] // 25)))
    q2i = [k for k in range(RZ['rep'].shape[0]) if str(RZ['label'][k]) == 'q2']
    okr = all(_sc.verify_realisation(RZ['rep'][k].astype(np.int64), RZ['lam'][k], RZ['sizes'][k], RZ['rows'][k])
              for k in sorted(set(sub20 + q2i)))
    EX = _js20.load(open(os.path.join(core.DATA, 's17_exclusions.json')))
    oke = all(_sc.verify_exclusion(np.array(e['ray'], dtype=np.int64), e['X'], e['Y'], e['y'], F_S=e['F_S'])[0]
              for e in EX)
    IQ = core.load('s17_ineq_classes').astype(np.int64)
    G760 = core.load('graphstate_vecs').astype(np.int64)
    p6_20 = core.perms31_s6()
    okg = all(int((G760[:, p6_20[k]] @ IQ.T).min()) >= 0 for k in range(0, 720, 7))
    HF20 = core.load('qlr5_H_facets10860').astype(np.int64)
    okq = all(bool((HF20 @ np.array(e['ray'], dtype=np.int64) >= 0).all()) for e in EX)
    gate('G20 A2 certs', okr and oke and okg and okq and len(q2i) == 1,
         f"realisations re-verified ({len(set(sub20 + q2i))} sampled incl. 2*q2); {len(EX)} exclusions re-derived exactly; "
         f"{IQ.shape[0]} inequality classes vs 760 graph states; excluded rays lie in QLR5")

    # G21 (R28): catalogue-wide cutting plane -- exclusion certificates (sampled; all with --full),
    # every excluded ray inside QLR5, catalogue status accounting.
    X21 = np.load(os.path.join(core.DATA, 's18_exclusions.npz'))
    n21 = X21['ray'].shape[0]
    idx21 = range(n21) if a.full else range(0, n21, 10)
    ok21 = all(_sc.verify_exclusion(*_sc.unpack_exclusion(X21, i)[:4], F_S=_sc.unpack_exclusion(X21, i)[4])[0] for i in idx21)
    ok21 = ok21 and bool(((HF20 @ X21['ray'].astype(np.int64).T) >= 0).all())
    ST21 = _js20.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results', '2026-10-01_s18-cutting-plane', 'summary.json')))
    acc21 = ST21['excluded'] + ST21['realised'] + ST21['undecided'] == ST21['orbits'] == 1917706
    gate('G21 cutting plane', ok21 and acc21 and core.load('s17_ineq_classes').shape[0] == 747,
         f"{len(idx21)}/{n21} exclusion certificates re-derived exactly; catalogue {ST21['excluded']}/{ST21['realised']}/{ST21['undecided']} "
         f"(excluded/realised/undecided); 747 inequality classes")

    # G22 (R28.1): the single-CI LP exactly as the a2-survivors job runs it (HiGHS pinned to one
    # thread).  The headline ray must come out infeasible at (AD, BE) = (9, 18) with status 2, every
    # disjoint pair of the realised q2 feasible with status 0, and no other status anywhere -- a
    # status 4 here means some LP in this process ran with a different HiGHS thread setting.
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'scripts'))
    import s17_ci_test as _ci22
    th0 = _nthreads()
    st275 = _ci22.ci_status(_ci22.h_norm(core.load('a2_frontier300')[275]), 9, 18)
    hq2 = _ci22.h_norm(RZ['rep'][q2i[0]].astype(np.int64))
    pq2 = _ci22.pairs_for(hq2, 'disjoint')
    stq2 = sorted({_ci22.ci_status(hq2, X, Y) for X, Y in pq2})
    th1 = _nthreads()
    gate('G22 CI-LP', st275 == 2 and stq2 == [0],
         f"#275 at (AD,BE): status {st275} (expect 2); q2: {len(pq2)} disjoint pairs, statuses {stq2} (expect [0]); "
         f"scipy {_ci22.scipy.__version__}, HiGHS {'pinned to 1 thread' if _ci22.SCIPY_PINS_THREADS else 'NOT pinnable (scipy < 1.11)'}, "
         f"threads {th0}->{th1}")

    # G23 (R29): third server session (a2-survivors).  Every server certificate re-verified
    # independently; the complete single-CI test recorded as feasible for all 2,566 orbits with no
    # solver error; the 14 realised orbits carry status 1 in the current catalogue status.
    base23 = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results')
    Z23 = np.load(os.path.join(base23, '2026-10-01_a2-survivors', 'a2_surv_certs.npz'))
    U23 = np.load(os.path.join(base23, '2026-10-01_s18-cutting-plane', 'undecided_reps.npy')).astype(np.int64)
    i23 = Z23['index'].astype(np.int64)
    okc23 = i23.size == 14 and all(
        np.array_equal(Z23['rep'][k].astype(np.int64), U23[i23[k]])
        and _sc.verify_realisation(U23[i23[k]], int(Z23['lam'][k]), Z23['sizes'][k], Z23['rows'][k]) for k in range(i23.size))
    CI23 = json.load(open(os.path.join(base23, '2026-10-01_a2-survivors', 'a2_surv_ci.json')))
    oki23 = len(CI23) == U23.shape[0] == 2566 and not any(x['infeasible_pairs'] or x['lp_errors'] for x in CI23)
    so23 = np.load(os.path.join(base23, '2026-10-01_s18-cutting-plane', 'catalogue_status.npz'))['status']
    sn23 = np.load(os.path.join(base23, '2026-10-01_a2-status', 'catalogue_status.npz'))['status']
    oks23 = bool((sn23[np.flatnonzero(so23 == 2)[i23]] == 1).all()) and bool((sn23[so23 != 2] == so23[so23 != 2]).all())
    gate('G23 a2-survivors', okc23 and oki23 and oks23,
         f"{i23.size}/14 server realisations re-verified; single-CI feasible for {len(CI23)}/2566, 0 solver errors; "
         f"all 14 marked realised; earlier verdicts unchanged")

    # G24 (R29): qudit realisations.  The rank formula S(X) = rank_GF(p) W[X, X^c] is checked against
    # exact state vectors (random 5-qutrit graph states, every bipartition); the 16 qutrit certificates
    # are re-verified with the independent plain-Python GF(p) elimination of epr1kit.stabcert.
    import s21_qudit_search as _q24
    w24, n24 = _q24.self_test(3, trials=8)
    Q24 = np.load(os.path.join(core.DATA, 's21_qudit_realisations.npz'))
    okq24 = all(_sc.verify_realisation_gfp(Q24['rep'][k].astype(np.int64), int(Q24['lam'][k]), Q24['sizes'][k],
                                           Q24['W'][k][:int(Q24['sizes'][k].sum()), :int(Q24['sizes'][k].sum())], int(Q24['p'][k]))
                for k in range(Q24['rep'].shape[0]))
    okq24 = okq24 and Q24['rep'].shape[0] == 16 and bool((sn23[Q24['catalogue_index']] == 1).all())
    gate('G24 qudit', w24 < 1e-9 and okq24,
         f"rank formula vs exact state vectors: max error {w24:.1e} ({n24} qutrits, 8 graphs, all cuts); "
         f"{Q24['rep'].shape[0]} qutrit certificates re-verified (p=3, lambda=1), all marked realised")

    # G25 (R29): CI + Ingleton.  18 exclusion certificates (Shannon + conditional Ingleton on the
    # extension A..F,Z) re-derived exactly by the independent checker; exact single-CI feasibility of the
    # same 18 orbits (24,057 stored witnesses + explicit trivial ones), i.e. they refute S7''; every
    # status-3 orbit violates an S6 image of a class in data/s22_ineq_classes.npy while no realised ray
    # does; final accounting and the 630 extreme rays of Stab5.
    X25 = json.load(open(os.path.join(core.DATA, 's22_ingleton_exclusions.json')))['certificates']
    oke25 = len(X25) == 18 and all(
        _sc.verify_exclusion_ext(np.array(c['ray'], dtype=np.int64), c['X'], c['Y'],
                                 [(tuple(d), v) for d, v in c['terms']], F_S=c['F_S'])[0] for c in X25)
    W25 = np.load(os.path.join(core.DATA, 's22_single_ci_witnesses.npz'))
    okw25 = W25['rays'].shape[0] == 18 and all(
        _sc.verify_single_ci_feasible(W25['rays'][k].astype(np.int64), W25['pair'][W25['ray'] == k],
                                      W25['den'][W25['ray'] == k], W25['num'][W25['ray'] == k])[0] for k in range(18))
    okw25 = okw25 and sorted(map(tuple, W25['rays'].tolist())) == sorted(tuple(c['ray']) for c in X25)
    C25 = core.load('s22_ineq_classes').astype(np.int64)
    I25 = np.unique(np.concatenate([C25[:, p6_20[k]] for k in range(720)]), axis=0)
    CAT = np.load(os.path.join(base23, '2026-09-12_a1-campaign', 'qlr_adj.npz'))['reps']
    st3 = np.flatnonzero(sn23 == 3)
    okc25 = bool(((CAT[st3].astype(np.int64) @ I25.T).min(1) < 0).all())
    REAL25 = np.vstack([core.load('stab5_extreme_reps'), np.load(os.path.join(core.DATA, 's16_realisations.npz'))['rep'].astype(np.int64),
                        Q24['rep'].astype(np.int64), G760])
    okr25 = bool(((REAL25 @ I25.T) >= 0).all())
    S25 = json.load(open(os.path.join(base23, '2026-10-01_a2-status', 'summary.json')))
    cnt25 = tuple(int((sn23 == c).sum()) for c in (0, 1, 2, 3))
    acc25 = (cnt25 == (1914541, 630, 2070, 465) and S25['excluded'] == 1915006 and S25['realised'] == 630
             and S25['undecided'] == 2070 and sum(cnt25) == 1917706)
    E25 = core.load('stab5_extreme_reps')[:630]          # R30 appends to this file; the R29 rows come first
    oke25b = E25.shape[0] == 630 and core.sha_rows(E25) == MAN['stab5_extreme_reps_r29']['sha256']
    gate('G25 CI+Ingleton', oke25 and okw25 and okc25 and okr25 and acc25 and oke25b,
         f"{len(X25)} exclusion certificates re-derived; 18 orbits exactly single-CI feasible (S7'' refuted); "
         f"{st3.size} status-3 orbits cut by {C25.shape[0]} classes, 0 realised rays cut; "
         f"catalogue {S25['excluded']}/{S25['realised']}/{S25['undecided']}; {E25.shape[0]} extreme rays of Stab5 (R29)")

    # G26 (R30): fourth server session (a2-ingleton).  Every infeasible verdict of the server's
    # CI+Ingleton LP (904 of 2,070, no solver error) carries an exact certificate re-derived by the
    # independent checker; the 904 inequalities fall into the 628 S6 classes of
    # data/s22_ineq_classes_a2i.npy; no image of a class is violated by a realised vector or cuts a ray
    # the server found feasible for all pairs; the 11 server realisations are re-verified; status update,
    # accounting, and the 641 extreme rays of Stab5 (= the realised catalogue orbits).
    X26 = json.load(open(os.path.join(core.DATA, 's22_ingleton_exclusions_a2i.json')))['certificates']
    CI26 = json.load(open(os.path.join(base23, '2026-10-02_a2-ingleton', 'a2i_ci.json')))
    U26 = np.load(os.path.join(base23, '2026-10-01_a2-status', 'undecided_reps.npy')).astype(np.int64)
    inf26 = sorted(x['index'] for x in CI26 if x['infeasible_pairs'])
    oke26 = (len(CI26) == U26.shape[0] == 2070 and not any(x['lp_errors'] for x in CI26) and len(inf26) == 904
             and [c['index'] for c in X26] == inf26
             and all(np.array_equal(np.array(c['ray'], dtype=np.int64), U26[c['index']])
                     and _sc.verify_exclusion_ext(U26[c['index']], c['X'], c['Y'],
                                                  [(tuple(d), v) for d, v in c['terms']], F_S=c['F_S'])[0] for c in X26))
    C26 = core.load('s22_ineq_classes_a2i').astype(np.int64)

    def _canon26(F):
        F = np.asarray(F, dtype=np.int64)
        return np.unique((F // int(np.gcd.reduce(np.abs(F))))[p6_20], axis=0)[0]
    okk26 = (C26.shape[0] == 628 and core.sha_rows(C26) == MAN['s22_ineq_classes_a2i']['sha256']
             and len({c['class'] for c in X26}) == 628
             and all(np.array_equal(_canon26(c['F_S']), C26[c['class']]) for c in X26))
    I26 = np.unique(C26[:, p6_20].reshape(-1, 31), axis=0)
    If26 = I26.T.astype(np.float32).copy()

    def _min26(V):          # exact in float32: |coef| <= 100, entries <= 30, 31 terms << 2^24
        V = np.asarray(V, dtype=np.int64)
        assert int(np.abs(I26).max()) * max(1, int(np.abs(V).max())) * 31 < 2 ** 24
        return min(float((V[s:s + 128].astype(np.float32) @ If26).min()) for s in range(0, V.shape[0], 128))
    E26 = core.load('stab5_extreme_reps')
    REAL26 = np.vstack([E26, np.load(os.path.join(core.DATA, 's16_realisations.npz'))['rep'].astype(np.int64),
                        Q24['rep'].astype(np.int64), G760])
    feas26 = U26[np.setdiff1d(np.arange(U26.shape[0]), inf26)]
    okr26 = _min26(REAL26) >= 0 and _min26(feas26) >= 0
    Q26 = np.load(os.path.join(base23, '2026-10-02_a2-ingleton', 'a2i_q2.npz'))
    i26 = Q26['index'].astype(np.int64)
    okq26 = i26.size == 11 and not set(i26.tolist()) & set(inf26) and all(
        np.array_equal(Q26['rep'][k].astype(np.int64), U26[i26[k]])
        and _sc.verify_realisation(U26[i26[k]], int(Q26['lam'][k]), Q26['sizes'][k], Q26['rows'][k]) for k in range(i26.size))
    sn26 = np.load(os.path.join(base23, '2026-10-02_a2-status', 'catalogue_status.npz'))['status']
    exp26 = sn23.copy(); st2_26 = np.flatnonzero(sn23 == 2)
    exp26[st2_26[inf26]] = 3; exp26[st2_26[i26]] = 1
    S26 = json.load(open(os.path.join(base23, '2026-10-02_a2-status', 'summary.json')))
    cnt26 = tuple(int((sn26 == c).sum()) for c in (0, 1, 2, 3))
    acc26 = (np.array_equal(sn26, exp26) and cnt26 == (1914541, 641, 1155, 1369)
             and (S26['excluded'], S26['realised'], S26['undecided']) == (1915910, 641, 1155)
             and np.array_equal(np.load(os.path.join(base23, '2026-10-02_a2-status', 'undecided_reps.npy')).astype(np.int64),
                                CAT[sn26 == 2].astype(np.int64)))
    src26 = np.load(os.path.join(core.DATA, 'stab5_extreme_src.npy'))
    oke26b = (E26.shape[0] == 641 and core.sha_rows(E26) == MAN['stab5_extreme_reps']['sha256']
              and sorted(map(tuple, E26[630:].tolist())) == sorted(map(tuple, Q26['rep'].astype(np.int64).tolist()))
              and bool((src26[630:] == 'cert_r30').all())
              and sorted(map(tuple, E26.tolist())) == sorted(map(tuple, CAT[sn26 == 1].astype(np.int64).tolist())))
    gate('G26 a2-ingleton', oke26 and okk26 and okr26 and okq26 and acc26 and oke26b,
         f"{len(X26)}/{len(inf26)} CI+Ingleton certificates re-derived (2070 orbits tested, 0 solver errors), "
         f"{C26.shape[0]} classes, 0 realised vectors and 0 server-feasible rays cut by {I26.shape[0]} images; "
         f"{i26.size}/11 server realisations re-verified; catalogue {S26['excluded']}/{S26['realised']}/{S26['undecided']}; "
         f"{E26.shape[0]} extreme rays of Stab5")

    # G27 (R30): DFZ five-variable linear rank inequalities on the 7-element extension
    # (scripts/s23_ci_dfz.py).  The 28 transcribed forms hold on random GF(2) subspace arrangements;
    # instance counts; every instance on A..F holds for h_norm of every realised vector.
    import s23_ci_dfz as _d27
    D27 = _d27.load_dfz()
    rng27 = np.random.default_rng(27)
    oko27 = True
    for _ in range(60):
        dim = int(rng27.integers(3, 7))
        gens = [rng27.integers(0, 2, size=(int(rng27.integers(1, 4)), dim)) for _j in range(5)]
        hr = {m: _sc.gf2_rank(np.vstack([gens[j] for j in range(5) if m >> j & 1])) for m in range(1, 32)}
        oko27 = oko27 and all(sum(v * hr[m] for m, v in c.items()) >= 0 for c in D27.values())
    inst27 = _d27.dfz_instances(ineqs=D27)
    rows27 = [r for r, _ in inst27 if not any(M >> 6 & 1 for M in r)]
    A27 = np.zeros((len(rows27), 64))
    for k, r in enumerate(rows27):
        for M, v in r.items():
            A27[k, M] = v
    H27 = np.stack([_ci22.h_norm(v) for v in REAL26], axis=1).astype(float)
    okv27 = float((A27 @ H27).min()) >= 0
    gate('G27 DFZ ext', len(D27) == 28 and oko27 and len(inst27) == 149520 and len(rows27) == 21480 and okv27,
         f"{len(D27)} DFZ forms (24 + 4 Ingleton) valid on 60 random GF(2) arrangements; {len(inst27)} instances on "
         f"A..F,Z ({len(inst27) - len(rows27)} with Z); the {len(rows27)} on A..F hold for all {REAL26.shape[0]} realised vectors")

    # G28 (R31): p = 2 versus p = 3.  The 16 orbits that R29 realised only with qutrits are realised by qubits
    # at lambda = 2 (GF(4) graph states at lambda = 1, field-reduced to qubit graph states by
    # scripts/s24_galois_search.py; data/s24_qubit_l2_certs.npz), and the 80 qubit-realised extreme rays of
    # Stab5 with the smallest qubit counts are realised by qutrits (73 at lambda = 1, 7 at lambda = 2;
    # data/s24_qutrit_certs.npz).  The field reduction is checked against direct ranks on random states.
    import s24_galois_search as _g28
    Q28 = np.load(os.path.join(core.DATA, 's24_qubit_l2_certs.npz'))
    ok28a = Q28['rep'].shape[0] == 16 and sorted(Q28['catalogue_index'].tolist()) == sorted(Q24['catalogue_index'].tolist())
    for k in range(Q28['rep'].shape[0]):
        r28 = Q28['rep'][k].astype(np.int64); sz28 = Q28['sizes'][k]; n28 = int(sz28.sum())
        j28 = int(np.flatnonzero(Q24['catalogue_index'] == Q28['catalogue_index'][k])[0])
        rows28 = np.zeros(64, dtype=np.uint64)
        for i28 in range(n28):
            rows28[i28] = np.uint64(sum(1 << b for b in range(n28) if Q28['G'][k][i28, b]))
        ok28a = (ok28a and np.array_equal(r28, Q24['rep'][j28].astype(np.int64)) and int(Q28['lam'][k]) == 2
                 and list(sz28) == list(_g28.party_sizes(r28, 2)) and _sc.verify_realisation(r28, 2, sz28, rows28))
    for k in (0, 7, 15):                     # re-derive from the stored GF(4) matrices with the s24 reduction
        r28 = Q28['rep'][k].astype(np.int64); N28 = int(sum(_g28.party_sizes(r28, 1)))
        okq28, c28 = _g28.certify(r28, 1, Q28['Wq'][k][:N28, :N28].astype(int).tolist(), 4)
        ok28a = ok28a and okq28 and c28 is not None and np.array_equal(c28['G'], Q28['G'][k][:2 * N28, :2 * N28])
    T28 = np.load(os.path.join(core.DATA, 's24_qutrit_certs.npz'))
    E28 = core.load('stab5_extreme_reps').astype(np.int64)
    src28 = np.load(os.path.join(core.DATA, 'stab5_extreme_src.npy'))
    qb28 = [k for k in range(E28.shape[0]) if src28[k] != 'cert_r29q3']
    nq28 = {k: sum(_g28.party_sizes(E28[k], 1)) for k in qb28}
    pick28 = sorted([k for k in qb28 if nq28[k] <= 20], key=lambda k: nq28[k])[:80]
    ok28b = (T28['stab5_row'].tolist() == pick28 and int((T28['lam'] == 1).sum()) == 73
             and int((T28['lam'] == 2).sum()) == 7)
    for k in range(T28['rep'].shape[0]):
        r28 = T28['rep'][k].astype(np.int64); sz28 = T28['sizes'][k]; n28 = int(sz28.sum())
        ok28b = (ok28b and np.array_equal(r28, E28[pick28[k]])
                 and _sc.verify_realisation_gfp(r28, int(T28['lam'][k]), sz28, T28['G'][k][:n28, :n28], 3))
    w28 = _g28.self_test(trials=4)
    gate('G28 p=2 vs p=3', ok28a and ok28b and w28 == 0,
         f"{Q28['rep'].shape[0]}/16 R29 qutrit-only orbits realised by qubits at lambda=2 (GF(4) + field reduction); "
         f"{T28['rep'].shape[0]}/80 smallest qubit-realised extreme rays realised by qutrits "
         f"({int((T28['lam'] == 1).sum())} at lambda=1, {int((T28['lam'] == 2).sum())} at lambda=2); reduction self-test max error {w28}")

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
