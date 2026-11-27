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

    # G29 (R32): the reduction Stab_n = h^-1(LR_{n+1}) = pi(LR_{n+1}) = CSS cone (docs/R32_reduction_and_local_
    # dimension.md, Theorem 1).  (a) The CSS entropy formula S = pi(r) and "the CSS state of a graph state's
    # Lagrangian code has entropy 2S" against explicit state vectors on random small codes / graph states over
    # GF(2) and GF(3); (b) n = 4: the 24 DFZ five-variable inequalities pulled back through h_S are implied by the
    # pulled-back Shannon + Ingleton rows (LP); (c) the 46 extreme rays of that cone (data/s26_stab4_certs.npz)
    # are valid, extreme, pairwise distinct and realised by qubit graph states (stabcert.verify_realisation_n).
    # --full adds the exact lrs enumeration (G29b), which must give exactly these 46 rays.
    import s26_reduction as _r29
    da29, db29 = _r29.self_test(trials=16)
    SH29, ING29, NEW29 = _r29.n4_rows()
    A29 = np.array(SH29 + ING29)
    bad29 = sum(1 for v in NEW29 if not _r29.lp_implied(A29, v))
    okc29, R29 = _r29.check_certs(_r29.CERTS, A29)
    okm29 = R29.shape[0] == MAN['s26_stab4_certs']['rows'] == 46 and core.sha_rows(R29) == MAN['s26_stab4_certs']['sha256']
    gate('G29 reduction', da29 < 1e-9 and db29 < 1e-9 and len(NEW29) == 820 and bad29 == 0 and okc29 and okm29,
         f'state-vector deviations {da29:.0e} (CSS formula) / {db29:.0e} (Lagrangian code = 2S); n=4: {len(NEW29)} '
         f'pulled-back DFZ rows, {bad29} not implied by Shannon + Ingleton; {R29.shape[0]} extreme rays valid, '
         f'extreme and realised by qubit graph states (Stab4 = h^-1(Shannon + Ingleton))')

    # G30 (R32): local-dimension dependence from seven parties on (Theorem 2).  The cube code (affine functions
    # on {0,1}^3) is identically self-dual over GF(2) and GF(3), h_S = 2r, contracting the origin gives exactly the
    # Fano / non-Fano configuration of DFZ arXiv:1311.4601 Thms 8.7 / 8.9, all 4x4 minors lie in {0, +-1, +-2}
    # (same matroid for every odd p), the pulled-back DFZ (65) [odd p] and (91) [p = 2] take the value -2 on the
    # qubit resp. qudit cube state, the entropies agree with explicit state vectors, the two vectors differ only on
    # the cut ABCZ | WXY+purifier, and the stored data reproduce.  Proposition 3: the CSS states of the doubled
    # Fano / non-Fano / T8 codes have h_S = M + |X| on the visible parties (balanced (65)/(91): -1).  Sanity (not a
    # proof): (65) and (91) hold on 40 random GF(3) resp. GF(2) subspace arrangements.
    import s27_cube_states as _c30
    G30, W30 = _c30.witnesses()
    D30 = _c30.load_dfz_char()
    F65, F91 = _c30.pullback(D30['dfz65_odd']), _c30.pullback(D30['dfz91_even'])
    ok30 = all(np.array_equal(_c30.dual(W30[p]['r'], 8), W30[p]['r'])
               and np.array_equal(W30[p]['h'][1:255], 2 * W30[p]['r'][1:255])
               and np.array_equal(_c30.contract(W30[p]['r'], 8, _c30.PUR), _c30.dfz_configuration(p)) for p in (2, 3))
    ok30 = ok30 and set(_c30.minors_set(G30)) <= {-2, -1, 0, 1, 2}
    S2_30, S3_30 = W30[2]['S'][1:128], W30[3]['S'][1:128]
    v65, v91 = int(F65 @ S2_30), int(F91 @ S3_30)
    c65 = int(D30['dfz65_odd'] @ _c30.dfz_configuration(2)); c91 = int(D30['dfz91_even'] @ _c30.dfz_configuration(3))
    sv30 = max(_c30.state_vector_check(G30, 2), _c30.state_vector_check(G30, 3))
    rng30 = np.random.default_rng(30)
    r65 = _c30.random_arrangement_min(D30['dfz65_odd'], 3, 40, rng30)
    r91 = _c30.random_arrangement_min(D30['dfz91_even'], 2, 40, rng30)
    Z30 = np.load(os.path.join(core.DATA, 's27_cube_witnesses.npz'))
    okz30 = (np.array_equal(Z30['S_qubit'], S2_30) and np.array_equal(Z30['S_qudit'], S3_30)
             and np.array_equal(Z30['F65_odd'], F65) and np.array_equal(Z30['F91_even'], F91)
             and core.sha_rows(np.stack([Z30['S_qubit'], Z30['S_qudit'], Z30['F65_odd'], Z30['F91_even']]).astype(np.int64))
             == MAN['s27_cube_witnesses']['sha256']
             and core.sha_rows(np.stack([D30['dfz65_odd'], D30['dfz91_even']])) == MAN['dfz_char_dependent']['sha256'])
    diff30 = np.flatnonzero(S2_30 != S3_30).tolist()
    dc30 = _c30.doubled_code_checks()          # Proposition 3: doubled codes give h_S = M + |X| on the visible parties
    okd30 = dc30['Fano'] == (True, -1) and dc30['non-Fano'] == (True, -1) and dc30['T8'][0]
    gate('G30 cube states', ok30 and v65 == -2 and v91 == -2 and c65 == -1 and c91 == -1 and sv30 < 1e-9
         and r65 >= 0 and r91 >= 0 and okz30 and diff30 == [_c30.ODD_CUT - 1] and okd30,
         f'self-dual, h_S = 2r, contraction = DFZ Fano/non-Fano (forms -> {c65}/{c91}); pulled-back (65) on the '
         f'qubit cube state {v65}, (91) on the qudit cube state {v91} (expect -2/-2); state-vector dev {sv30:.0e}; '
         f'differ only on ABCZ; random-arrangement minima {r65}/{r91} (>= 0); doubled Fano/non-Fano/T8 codes: '
         f'h_S = M + |X| {okd30}, balanced values {dc30["Fano"][1]}/{dc30["non-Fano"][1]}')

    # G31 (R33): any two distinct primes give incomparable stabilizer cones (Theorem 2', docs/R32_reduction_and_local_
    # dimension.md section 4.4).  Generalised cube states for t = 3 (p = 3 against q = 5, 7) and t = 5 (5 against 7):
    # explicit self-duality weights (G diag(d) G^T = 0), identical self-duality, h_S = 2r, contraction of the origin =
    # the configuration A_i = <e_i>, B_i = <c - e_i>, C = <c>; Lemma 4 (proved in the note): LB [valid for char not
    # dividing t] = -2 on the GF(p) state and LA [valid for char dividing t] = -2 on the GF(q) states (pulled back);
    # Pena-Sarria's Example 6 (b)/(a) give the same values (cross-check); the states differ in exactly one coordinate;
    # the stored data reproduce.  G31b (--full): sanity of LA/LB on all {0,1}-point configurations for t = 2.
    import s28_generalized_cubes as _g31
    ok31, vals31, nd31 = True, [], []
    for t31, ps31, path31, key31 in ((3, [3, 5, 7], 's28_gencube_t3.npz', 's28_gencube_t3'), (5, [5, 7], 's28_gencube_t5.npz', 's28_gencube_t5')):
        R31 = _g31.analyse(t31, ps31)
        for q31, d31 in R31['fields'].items():
            ok31 = ok31 and d31['weights'] and d31['self_dual'] and d31['h2r'] and d31['contraction']
            if q31 == t31:
                ok31 = ok31 and d31['FB'] == -2 and d31['FA'] >= 0 and d31['Fb'] == -2 and d31['Fa'] >= 0
                vals31.append(d31['FB'])
            else:
                ok31 = ok31 and d31['FA'] == -2 and d31['FB'] >= 0 and d31['Fa'] == -2 and d31['Fb'] >= 0
                vals31.append(d31['FA'])
        for q31 in ps31[1:]:
            nd31.append(int((R31['fields'][ps31[0]]['S'] != R31['fields'][q31]['S']).sum()))
        Z31 = np.load(os.path.join(core.DATA, path31))
        St31 = np.stack([R31['fields'][q]['S'] for q in ps31])
        rows31 = np.vstack([St31, R31['FA'][None, :], R31['FB'][None, :], R31['Fa'][None, :], R31['Fb'][None, :]]).astype(np.int64)
        ok31 = (ok31 and np.array_equal(Z31['S'], St31) and all(np.array_equal(Z31[k], R31[k]) for k in ('FA', 'FB', 'Fa', 'Fb'))
                and core.sha_rows_wide(rows31) == MAN[key31]['sha256_int16'])
    gate('G31 any two primes', ok31 and nd31 == [1, 1, 1],
         f'generalised cubes t=3 (GF 3/5/7) and t=5 (GF 5/7): weights, self-dual, h_S = 2r, contraction = configuration; '
         f'Lemma 4 violations {vals31} (expect all -2; Pena-Sarria Example 6 agrees); differing coordinates {nd31} (expect 1 each)')

    # G32 (R34): Z_d qudits -- the cone remembers the local dimension (Theorem 5, docs/R32_reduction_and_local_
    # dimension.md section 4.6).  Generalised cube states over Z_N, exactly (units of log p): free code, unit weights
    # with G diag(d) G^T = 0 mod N, identical self-duality, h_S = 2r, contraction of the origin = the J - I
    # configuration over Z_N; Lemma 4(G) values: Z_4, t = 2: (LA)_2 = (LB)_2 = -2 bits (outside every prime cone);
    # Z_8, t = 4: (LA)_4 = -2 bits (Z_8 vs Z_4); Z_9, t = 3: (LA)_3 = -2 trits (Z_9 vs Z_3); partners on the smaller
    # group satisfy the form ((LA)_4 on Z_4, t = 4: 0).  The Z_4 cube state agrees with its state vector and equals
    # qubit cube (bits) + qutrit cube (trits).  Stored data reproduce.  G32b (--full): Z_4, t = 6 ((LA)_6 = -2 bits,
    # the case Z_4 vs Z_6) and sanity tests of Lemma 4(G) on point configurations.
    import s29_zd_stabilizer as _z32
    Z32 = np.load(os.path.join(core.DATA, 's29_zd_witnesses.npz'))
    ok32, vals32 = True, {}
    for (t32, N32), want in (((2, 4), (-2, -2)), ((4, 8), (-2, -4)), ((3, 9), (-2, -2))):
        W32 = _z32.analyse(t32, N32)
        ok32 = (ok32 and W32['free'] and W32['weights'] and W32['self_dual'] and W32['h2r'] and W32['contraction']
                and (W32['FA'], W32['FB']) == want and W32['cfg_LA'] == -1
                and np.array_equal(Z32[f't{t32}_N{N32}'], W32['S'])
                and core.sha_rows_wide(W32['S'][None, :]) == MAN['s29_zd_witnesses'][f'sha256_int16_t{t32}_N{N32}'])
        vals32[f'Z{N32},t={t32}'] = (W32['FA'], W32['FB'])
    P32 = _z32.analyse(4, 4)
    sv32, size32 = _z32.state_vector_dev(2, 4)
    sum32 = _z32.cube_sum_identity()
    gate('G32 Z_d qudits', ok32 and P32['FA'] == 0 and P32['h2r'] and sv32 < 1e-9 and size32 == 256 and sum32,
         f'generalised cubes over Z_4 (t=2), Z_8 (t=4), Z_9 (t=3): free, unit weights, self-dual, h_S = 2r, contraction = '
         f'J - I configuration; (LA, LB) pulled back {vals32} (expect (-2,-2), (-2,-4), (-2,-2)); partner (LA)_4 on Z_4, '
         f't=4: {P32["FA"]} (expect 0); Z_4 cube state-vector dev {sv32:.0e}; Z_4 cube = qubit cube + qutrit cube {sum32}')

    # G33 (R35, R36): how far down the dimension dependence goes (docs/R32_reduction_and_local_dimension.md section
    # 4.7).  Theorem 7: the 46 extreme rays of Stab_4 have integer graph states with unimodular cut blocks (exact Smith
    # forms; also re-verified over GF(p), p <= 13), so Stab_4 is the same for every local dimension.  Proposition 8: for
    # the eight 7-variable characteristic-dependent inequalities considered (DFZ (65), (91); (LA)_2, (LB)_2; Pena-Sarria
    # (a), (b); BKL (6.8), (6.17) corrected), f + f o D is an exact nonnegative combination of Shannon inequalities
    # (stored certificates).  Proposition 6 (one cut): rank changes of the generalised cube points over GF(l) vs Q are
    # exactly the predicted balanced sets, t = 2..5 (t = 6, 7 in G33b).  Stored data match the manifest.
    import hashlib as _h33
    import s30_dimension_profile as _d33
    tot33, good33 = _d33.n4_check()
    six33 = _d33.six_check()
    oc33 = {t: _d33.one_cut_check(t) for t in (2, 3, 4, 5)}
    ok_oc33 = all(ok for res in oc33.values() for (_c, _k, ok) in res.values())
    files33 = all(_h33.sha256(open(os.path.join(core.DATA, f + '.npz'), 'rb').read()).hexdigest() == MAN[f]['file_sha256']
                  for f in ('s30_stab4_unimodular', 's30_selfdual_shannon_certs'))
    gate('G33 dimension profile', tot33 == 46 and good33 == 46 and len(six33) == 8 and all(ok for _n, ok in six33)
         and ok_oc33 and files33,
         f'four parties: {good33}/{tot33} rays lifted to every field (unimodular cut blocks); six parties: '
         f'{sum(ok for _n, ok in six33)}/{len(six33)} forms (expect 8/8) with exact Shannon certificates for f + f o D; '
         f'one cut t=2..5: predicted rank changes {ok_oc33}; stored files match the manifest {files33}')

    # G34 (R36): the seven-variable inequalities of Blasiak-Kleinberg-Lubetzky (arXiv:1108.2489v1, eqs. (6.8) and (6.17),
    # transcribed in data/bkl_char_dependent.csv and checked against the rendered PDF).  (6.8) [characteristic 2] is 0 on
    # the Fano and -1 on the non-Fano configuration; (6.17) [characteristic != 2] is -1 on the Fano and 0 on the
    # non-Fano.  As printed, (6.17) is false over every field: three independent lines A, B, C with Z = <A + B + C>
    # (W = X = Y = 0) give -3 over GF(2), GF(3), GF(5); the corrected form gives 6, (6.8) gives 2.  Transcription sha
    # matches the manifest.  G34b (--full): minima over all {0,1}-point configurations of GF(p)^3.
    K34 = _d33.bkl_check()
    want34 = {'fano2': {'bkl68_even': 0, 'bkl617_odd': -1, 'bkl617_odd_corrected': -1},
              'nonfano3': {'bkl68_even': -1, 'bkl617_odd': 0, 'bkl617_odd_corrected': 0},
              'counterexample': {p: {'bkl68_even': 2, 'bkl617_odd': -3, 'bkl617_odd_corrected': 6} for p in (2, 3, 5)}}
    B34 = _c30.load_dfz_char(_d33.BKL_CHAR)
    sha34 = core.sha_rows(np.stack([B34[k] for k in MAN['bkl_char_dependent']['ids']])) == MAN['bkl_char_dependent']['sha256']
    gate('G34 BKL', K34 == want34 and sha34,
         f'Fano/GF(2): {K34["fano2"]}; non-Fano/GF(3): {K34["nonfano3"]}; three lines + Z = A + B + C over GF(2/3/5): '
         f'{K34["counterexample"][2]} (printed (6.17) false: -3); transcription sha {sha34}')

    # G35 (R36): six parties already separate qubits from every odd local dimension (Theorem 9, docs/R32_reduction_and_
    # local_dimension.md section 4.8).  The transfer form f9 (Lemma 4 (LB)_2 applied to the derived family
    # Fhat_i = sum_L U_i cap (U_j + U_k), each term bounded by Lemmas 0, A-E) is recomputed exactly from its definition
    # and matches the stored data and the manifest.  The seven-qubit simplex-code state has h_S = r_F + r_F^* (state
    # vector, 128 amplitudes) and f9(h_S) = -3; f9(r_F) = -3, f9 = 81 / 33 on the non-Fano h / r; (f9 + f9 o D)(h_F) = -6,
    # so the self-dual symmetrisation is not a Shannon consequence; on Fano-invariant self-dual vectors
    # f9 = 630 v_P + 90 alpha + 63 nu - 1263 mu.  Sanity (not a proof): with the actual subspaces, Lemmas 0, A-E hold on
    # 30 random families over GF(2) and GF(3), and over GF(3) also (LB)_2 on the derived family and f9 are >= 0; f9 >= 0
    # on 40 random arrangements over GF(3) and GF(5).  Exact faces: the Shannon elemental inequalities tight at h_F
    # (h_NF) have rational rank 125 and those tight at r, r^* rank 126, so the face of Gamma_7 at h_F (h_NF) is
    # cone(r, r^*) (upper bound from the kernel, lower bound from the rank over GF(1000003)).
    # G35b (--full): all {0,1}-point configurations of GF(p)^3
    # (min -3 for p = 2, >= 0 for p = 3, 5) and 200 random families per field for the lemma checks (GF(2), GF(3), GF(5)).
    from fractions import Fraction as _Fr
    import s31_six_party as _s35
    f35 = _s35.f9()
    Z35 = np.load(os.path.join(core.DATA, 's31_six_party.npz'))
    okf35 = (np.array_equal(Z35['f9'], f35) and core.sha_rows(f35[None, :]) == MAN['s31_six_party']['sha256']
             and int(np.count_nonzero(f35[1:])) == MAN['s31_six_party']['terms']
             and _h33.sha256(open(os.path.join(core.DATA, 's31_six_party.npz'), 'rb').read()).hexdigest()
             == MAN['s31_six_party']['file_sha256'])
    w35 = _s35.witness(f35)
    okw35 = (w35['dev'] < 1e-9 and w35['purity'] < 1e-9 and abs(w35['f(h_S)'] + 3) < 1e-9 and w35['f(h_F)'] == -3
             and w35['f(r_F)'] == -3 and w35['f(h_NF)'] == 81 and w35['f(r_NF)'] == 33 and w35['(f + f o D)(h_F)'] == -6
             and w35['self-dual h_F'] and np.array_equal(Z35['h_F'], _s35.s27.dfz_configuration(2)
                                                         + _s35.h_star(_s35.s27.dfz_configuration(2))))
    inv35 = _s35.invariant_form(f35)
    rng35 = np.random.default_rng(35)
    lem35 = {q: _s35.lemma_checks(q, 30, rng35, f35) for q in (2, 3)}
    oklem35 = (all(v >= 0 for k, v in lem35[2].items() if k.startswith('Lemma'))
               and all(v >= 0 for v in lem35[3].values()))
    st35 = _s35.stress(40, rng35, f35, fields=(3, 5))
    fc35 = _s35.face_check()
    okfc35 = (all(fc35[k][1] == 126 and fc35[k][2] >= 0 for k in ('r_F', 'r_F*', 'r_NF', 'r_NF*'))
              and all(fc35[k][1] == 125 and fc35[k][2] >= 0 for k in ('h_F', 'h_NF')))
    gate('G35 six parties', okf35 and okw35 and inv35 == (_Fr(630), _Fr(90), _Fr(63), _Fr(-1263)) and oklem35
         and all(v >= 0 for v in st35.values()) and okfc35,
         f'f9 recomputed = stored {okf35} ({w35["terms"]} terms); simplex-code state: h_S = r_F + r_F^* (dev {w35["dev"]:.0e}), '
         f'f9(h_S) = {w35["f(h_S)"]:.0f} (expect -3), f9(r_F) = {w35["f(r_F)"]}, non-Fano h/r {w35["f(h_NF)"]}/{w35["f(r_NF)"]}; '
         f'(f9 + f9 o D)(h_F) = {w35["(f + f o D)(h_F)"]} (not Shannon); invariant form '
         f'{"/".join(str(x) for x in inv35)} on (v_P, alpha, nu, mu); lemma slacks >= 0 {oklem35} (GF(2) chain min '
         f'{lem35[2]["f9(h)"]}); random minima {st35}; Shannon faces of h_F, h_NF = cone(r, r*) (ranks 125/126) {okfc35}')

    # G36 (R38): five parties -- unimodular realisations of the known extreme rays of Stab_5 (roadmap item S5(a);
    # docs/R32_reduction_and_local_dimension.md section 4.9).  Every stored integer matrix W of
    # data/s32_stab5_unimodular.npz is a graph state whose 31 visible cut blocks have Smith invariants 0 and 1 only and
    # rational rank kappa * r(X) (exact diagonal forms over Z, M1), so W mod d realises kappa * r for every d >= 2.  The
    # rows are distinct rows of data/stab5_extreme_reps.npy, stored with the catalogue representative; the file and its
    # content hash match the manifest.  G36b (--full): the second exact method (rational rank and t x t minors with gcd
    # 1, M2), ranks over GF(p) for p <= 19, Z_4 and Z_6 state vectors of the small certificates, and the self-test of
    # scripts/s32_unimodular.py (neutral moves keep Smith forms; ring expansion; the Z_d cut formula).
    import s32_unimodular as _s36
    S36 = _s36.check(full=False)
    M36 = MAN['s32_stab5_unimodular']
    file36 = (_h33.sha256(open(os.path.join(core.DATA, 's32_stab5_unimodular.npz'), 'rb').read()).hexdigest()
              == M36['file_sha256'] and _s36.content_sha() == M36['sha256'])
    gate('G36 five-party unimodular', S36['ok'] == S36['stored'] == S36['distinct_rows'] == M36['rows'] and S36['rep_ok']
         and file36,
         f'{S36["ok"]}/{S36["stored"]} stored certificates exact (unimodular cut blocks), distinct catalogue rows; kappa '
         f'{dict(sorted(S36["kappa"].items()))}; so {S36["ok"]} of the 641 known extreme rays of Stab_5 lie in Stab^(Z_d)_5 '
         f'for every d; file and content match the manifest {file36}')

    # G37 (R38): covering certificates for known extreme rays of Stab_5 without a unimodular realisation
    # (data/s32_stab5_cover.npz; docs/R32_reduction_and_local_dimension.md section 4.9).  For each stored row: rational
    # realisations W_i (integer matrices with the right rational cut ranks) whose bad primes (the primes dividing a
    # Smith invariant of some cut block, recomputed exactly) intersect inside the primes of the stored GF(p)
    # realisations of the row (qubit p = 2, qutrit p = 3).  Then every prime p is covered, so the ray lies in Stab^(p)_5
    # for every prime p and, by Corollary 3 of the paper, in Stab^(Z_d)_5 for every d.  No row is also in the
    # unimodular file; the file and its content hash match the manifest.
    S37 = _s36.check_cover()
    M37 = MAN['s32_stab5_cover']
    file37 = (_h33.sha256(open(os.path.join(core.DATA, 's32_stab5_cover.npz'), 'rb').read()).hexdigest() == M37['file_sha256']
              and _s36.content_sha(_s36.COVER_DEFAULT) == M37['sha256'])
    gate('G37 five-party covering certificates', S37['pieces_ok'] == S37['pieces'] == M37['pieces']
         and S37['rows_covered'] == S37['rows'] == M37['rows'] and S37['overlap_unimodular'] == 0 and file37,
         f'{S37["rows_covered"]}/{S37["rows"]} rows covered for every prime by {S37["pieces_ok"]}/{S37["pieces"]} exact '
         f'pieces (rational realisations with recomputed bad primes, GF(2) / GF(3) realisations); none overlaps the '
         f'unimodular file; together with G36, {S36["ok"] + S37["rows_covered"]} of the 641 known extreme rays lie in '
         f'Stab^(Z_d)_5 for every d; file and content match the manifest {file37}')

    # G38 (R40): six parties, the other direction -- states with at most one qudit per system (roadmap item M3;
    # docs/R32_reduction_and_local_dimension.md section 4.10).  data/s34_seven_qudit_vectors.npz holds every entropy
    # vector of seven-qubit and of seven-qutrit stabilizer states with at most one qudit per system (exhaustive
    # enumeration of weighted graph states, cc/s35_enum7.c), one representative of each S_7 orbit of the qutrit-only
    # vectors with a GF(3) graph realising it, and qubit graph states realising twice the vector.  Checked exactly: the
    # representatives are realised by their graphs over GF(3) and are qutrit-only; their orbits partition the
    # qutrit-only vectors; each stored qubit graph realises 2 S over GF(2); the only orbit without one is the orbit of
    # the non-Fano state h_NF; the qubit-only vectors are the 30 relabellings of the Steane state; file and content match
    # the manifest.  G38b (--full): the GF(2) list recomputed in numpy (all 2^21 graphs, independently of the C code),
    # 20,000 random GF(3) graphs inside the GF(3) list, the tensor-product and distance checks of section 4.10, the
    # Python target search for S_F and S_NF over GF(2), GF(3), GF(4) against the counts of cc/s36_target7.c (manifest;
    # Proposition 10), the line pieces of the Fano and non-Fano CSS states with one and two qudits per system after
    # random local symplectic maps (Proposition 10, steps 1 and 1'; R41), and the self-test of
    # scripts/s34_six_party_survey.py.
    import s34_six_party_survey as _s38
    S38 = _s38.check()
    M38 = MAN['s34_seven_qudit_vectors']
    file38 = (_h33.sha256(open(os.path.join(core.DATA, 's34_seven_qudit_vectors.npz'), 'rb').read()).hexdigest()
              == M38['file_sha256'] and _s38.content_sha() == M38['sha256'])
    gate('G38 six parties, at most one qudit per system', _s38.check_ok(S38) and S38['gf2'] == M38['gf2']
         and S38['gf3'] == M38['gf3'] and S38['orbits'] == M38['orbits'] and file38,
         f'{S38["gf2"]} qubit and {S38["gf3"]} qutrit vectors; {S38["qutrit_only"]} qutrit-only in {S38["orbits"]} orbits '
         f'({S38["orbits_all_singles_one"]} with every S(P) = 1), all realised by their GF(3) graphs; qubit realisations '
         f'of 2 S exact for {S38["qubit2_ok"]}/{S38["qubit2_stored"]} orbits, the remaining orbit is that of h_NF '
         f'{S38["without_qubit2_is_nonfano"]}; qubit-only = the 30 relabellings of the Steane state '
         f'{S38["qubit_only_is_steane"]}; file and content match the manifest {file38}')

    # G39 (R42): six parties, the other direction for the closed cones (Theorem 12; docs/R32_reduction_and_local_dimension.md
    # section 4.11).  The inequality f12 >= 0, valid on every qubit stabilizer state of seven systems (Z the purifier), is
    # built from its definition (transfer_form) and pinned by the manifest (sha256 of the 128 coefficients, and the 31
    # terms on the cuts of the visible parties).  Checked exactly: f12(S_NF) = -1 and f12(S_F) = 0; f12 >= 0 on the GF(2)
    # list of G38 (every seven-qubit state with at most one qubit per system; minimum 0, attained 38 times); on the GF(3)
    # list f12 < 0 only at S_NF.  Corollary 13 (R43): at the seven-point state over Z_d (the CSS state of the seven nonzero
    # vectors of {0,1}^3 over Z_d), 2^f12 = gcd(d, 2)/d exactly for d = 2..16, and the column sets of its generator matrix
    # have Smith invariants 1 except W, X, Y (1, 1, 2).  G39b (--full): every step of the proof (transfer_chain: dim L_T = k(T), the restrictions
    # of the line pieces and their orthogonality, Lemmas T1-T4, (LA)_2 on the arrangement of Lemma T5) on 40 states -- random
    # compressed qubit graph states, compressed CSS states of random GF(2) families pushed towards the non-Fano
    # dependencies, copies of the Steane state after random local symplectic maps -- the dual arrangement of Lemma T5
    # built explicitly, the compression (symplectic reduction) against uncompressed states, and the seven-point states
    # over Z_4, Z_6, Z_8 against their state vectors (Corollary 13).
    T39 = _s38.transfer_check(log=lambda *_x: None)
    M39 = MAN['s34_transfer_form']
    pin39 = _s38.transfer_sha() == M39['sha256'] and [[n, c] for n, c in _s38.transfer_terms()] == M39['terms']
    gate('G39 six parties, the closed cones (Theorem 12)', _s38.transfer_ok(T39) and pin39
         and T39['f12_nonfano'] == M39['nonfano'] and T39['f12_fano'] == M39['fano'] and T39['gf2_zero'] == M39['gf2_zero']
         and T39['gf3_negative'] == M39['gf3_negative'],
         f'f12(S_NF) = {T39["f12_nonfano"]} (expect -1), f12(S_F) = {T39["f12_fano"]} (0); GF(2) list: minimum '
         f'{T39["gf2_min"]} (0), {T39["gf2_zero"]} zeros (38); GF(3) list: {T39["gf3_negative"]} negative (1), the '
         f'non-Fano vector {T39["gf3_negative_is_nonfano"]}; {T39["terms"]} terms, coefficients match the manifest {pin39}; '
         f'seven-point states over Z_d, d = {T39["seven_point_d"]}: 2^f12 = gcd(d, 2)/d {T39["seven_point_closed_form"]} '
         f'(Z_4: {T39["seven_point_z4"]}), Smith invariants {T39["seven_point_smith"]}')

    # G40 (R44): six parties, powers of two (Theorem 14; docs/R32_reduction_and_local_dimension.md section 4.12).  The
    # inequality f^(a) = f12 + (1 - 1/a) U_Z >= 0 holds on Stab^(Z_{2^a})_6 (Z the purifier); U_Z is pinned by the
    # manifest (sha256 of its 128 coefficients), f12 by G39.  Checked exactly: f^(a)(S_NF) = -1/a and f^(a)(S_F) = 1 - 1/a
    # for a = 1..4; f^(1) = f12; at the seven-point state over Z_{2^j}, 2^(a f^(a)) = 2^(a - j) for a = 1..4, j = 1..5,
    # and over Z_d with an odd factor 2^(a f^(a)) = 2^(a [d even]) / d; the seven-point states over Z_2, Z_4, Z_8 realised
    # as Z_8 stabilizer states; and the tight case: on the seven-point state over Z_{2^a} (a = 1, 2, 3) every step of the
    # proof holds (scripts/s35_two_power_transfer.py, chain), f^(a) = 0, Lemma 4(G') is attained and delta_P = U_P = a.
    # G40b (--full): every step of the proof (chain: the sizes of the line pieces and of their restrictions,
    # orthogonality, Lemmas T1-T4 over Z_{2^a}, the dual arrangement of Lemma T5 built explicitly, Lemma 4(G') on it,
    # the algebra of the assembly) on 60 random states over Z_4 and over Z_8 (graph, CSS of non-free codes, non-free
    # Lagrangians after random symplectic maps, families near the seven-point configuration, sums of seven-point
    # states); the entropy formula against state vectors; Lemma 4(G') on 150 configurations per case (valid forms >= 0,
    # the coefficient for a smaller exponent < 0); and two exhaustive families (cc/s38_minform_zn.c): every 3 x 7 0/1
    # code over Z_4 and every lift seven-point + 2P over Z_8 (2^21 codes each), where a f^(a) >= 0 with minimum 0 while
    # forms the proof does not give go negative.
    import s35_two_power_transfer as _s40
    T40 = _s40.check(log=lambda *_x: None)
    M40 = MAN['s35_two_power_form']
    gate('G40 six parties, powers of two (Theorem 14)', _s40.check_ok(T40) and T40['uz_sha'] == M40['uz_sha256']
         and T40['nf'] == {int(x): v for x, v in M40['nonfano'].items()},
         f'f^(a)(S_NF) = {T40["nf"]} (expect -1/a), f^(a)(S_F) = {T40["steane"]} (1 - 1/a); f^(1) = f12 {T40["f1_is_f12"]}; '
         f'seven-point states over Z_(2^j): f^(a) = 1 - j/a {T40["seven_point_two_powers"]}, over Z_d with odd factors '
         f'{T40["seven_point_other_d"]}, realised over Z_8 {T40["seven_point_realised"]}; tight case, every step of the proof '
         f'on the seven-point state over Z_(2^a), a = 1, 2, 3: {T40["tight"]}; U_Z matches the manifest '
         f'{T40["uz_sha"] == M40["uz_sha256"]}')

    # G41 (R45): six parties, one qudit per system over other fields (docs/R32_reduction_and_local_dimension.md section
    # 4.10, "Other fields").  data/s36_one_qudit_fields.npz holds the 136 S_7 orbits of the entropy vectors of GF(q)-linear
    # stabilizer states of seven systems with one qudit each, over 13 fields from GF(2) to GF(49), and for every field
    # and orbit met over it a graph realising the representative.  Checked exactly: the representatives and orbit sizes;
    # all 1,646 certificates, with a separate implementation of GF(q); over GF(2) and GF(3) the certified orbits are those
    # of the exhaustive lists of G38; over every field with at least seven elements the orbits are 134 common ones plus
    # S_NF (odd characteristic) or S_F (characteristic 2); the unions per characteristic (3, 5, 7, 11, 13) are equal;
    # every orbit missing over GF(2), GF(3), GF(4), GF(5) reappears over a larger field of the same characteristic; the
    # GF(3) orbits are certified over GF(5), GF(7), GF(11), GF(13); the symmetric, U_{2,6} and AME orbits are met exactly
    # over the fields with at least 7, 5 and 3 elements (Proposition 15); file and content match the manifest.
    # G41b (--full): cc/s40_enum7q.c, exhaustive over GF(2) (with and without the relabelling reduction), GF(3), GF(4)
    # and GF(5) (~5 min on two cores): closed under S_7, the lists of G38 and the certified orbits of the table, the
    # counts as in the manifest; exhaustive target searches for Proposition 15 (cc/s36_target7.c); a short run of the
    # sampler cc/s39_sample7q.c over GF(7), GF(8), GF(9); Proposition 16 (the 70 minors of [I_4 | J_4 - I_4], brute-force
    # representability of T8, T8', AG(3,2) over GF(2) .. GF(9), and S_T8 against S_T8').
    import s36_one_qudit_fields as _s41
    T41 = _s41.check(log=lambda *_x: None)
    M41 = MAN['s36_one_qudit_fields']
    file41 = (_h33.sha256(open(os.path.join(core.DATA, 's36_one_qudit_fields.npz'), 'rb').read()).hexdigest()
              == M41['file_sha256'] and _s41.content_sha() == M41['sha256'])
    gate('G41 six parties, one qudit per system over 13 fields', _s41.check_ok(T41) and file41
         and {int(k): v for k, v in M41['per_field'].items()} == T41['per_field'],
         f'{T41["orbits"]} orbits, {T41["vectors"]} vectors; certified orbits per field {T41["per_field"]}; certificates '
         f'{T41["certificates"]} re-checked, {T41["certificates_bad"]} wrong (0); GF(2)/GF(3) = the lists of G38 '
         f'{T41["gf2_list_is_union"] and T41["gf3_list_is_union"]}; fields >= 7 elements: {T41["common"]} common orbits + '
         f'S_NF / S_F {T41["big_fields_structure"]}; odd characteristics equal {T41["odd_unions_equal"]}; missing orbits '
         f'reappear in the same characteristic {T41["reappear"]}; symmetric vector over {T41["symmetric_fields"]}; file and '
         f'content match the manifest {file41}')

    # G42 (R46): fewer parties for two primes, and the T8 / T8' vectors (docs/R32_reduction_and_local_dimension.md
    # sections 4.13, 4.14, 4.5).  Lemma 4*: the forms (LB*), (LA*) for M = 4, 6 recomputed and pinned by the manifest;
    # on T(t) = [I_M | J_M - I_M] over GF(2), GF(3), GF(4), GF(5), GF(7), GF(8), GF(9), GF(11), M = 4, 5, 6, (LB*) = -1
    # exactly when the characteristic divides t and (LA*) = -1 exactly when it does not (0 otherwise), for every choice
    # of j(i) at M = 4; nonnegative balancing coefficients; over Z_N, (LA*) = -log(N / gcd(t, N)) and
    # (LB*) = -log gcd(t, N); Theorem 17: the balanced forms on the doubled-code states of T(3) and T(5) (entropies from
    # the [G | G] code), -1 on the witnesses.  Proposition 18: the GF(9), GF(3), GF(4), GF(8), GF(2) configurations and
    # pi(g) = 12 S_T8', 2 S_T8', 6 S_T8; the two-copy GF(3) family of h_T8'; the integer two-copy family of h_T8 (rank
    # function over Q exact, Smith-form primes {2, 3}, equal to h_T8 mod 5, 7, 11, 13, 17 and not mod 2, 3).  DFZ arXiv:1401.2507: Theorem 4.1 transcribed twice, -1 on [I_4 | J_4 - I_4] in
    # every characteristic but 3, its balanced form >= 42 on h_T8' and h_T8 in all 8! relabellings and -2 on h_AG;
    # Theorem 3.1 re-expanded from the bounds of its proof differs from the print in H(X|A,C,D), H(C|B,X,Y), H(C|A,W,Y).
    # G42b (--full): exhaustive over all 8-tuples of subspaces of GF(2)^2, GF(3)^2, GF(5)^2 (both forms, DFZ 4.1 and the
    # re-expanded 3.1; the printed 3.1 fails 125 times over GF(2)^2); exhaustive over coordinate lines A..D and all
    # lines W..Z of GF(2)^4 (LB*) and GF(3)^4 (LA*), with the counts of the configurations violating the forms in the
    # other characteristic (1 and 8: AG(3,2) and T8); 1.6 M random arrangements; M = 8 values over GF(7), GF(3).
    import s37_derived_c as _s42
    T42 = _s42.check(log=lambda *_x: None)
    M42 = MAN['s37_derived_forms']
    file42 = (_h33.sha256(open(os.path.join(core.DATA, 's37_derived_forms.npz'), 'rb').read()).hexdigest()
              == M42['file_sha256'] and _s42.content_sha() == M42['sha256'])
    gate('G42 fewer parties for two primes (Lemma 4*, Theorem 17), T8 / T8\' (Proposition 18)',
         _s42.check_ok(T42) and file42 and T42['sha'] == M42['forms'] and T42['witness'] == M42['witness'],
         f'values on T(t) as claimed {T42["values_ok"]} (M = 4, 5, 6; eight fields); every j(i): {T42["jmap_values"]}; '
         f'Z_N {T42["zn_ok"]}; witnesses {T42["witness"]}; Proposition 18 {T42["prop18_ok"]} (one cut '
         f'{T42["prop18"]["one_cut"]}); DFZ 4.1: transcriptions agree {T42["dfz41_transcriptions_agree"]}, Thm 4.2 values '
         f'{T42["dfz41_thm42"]}, balanced min over 8! relabellings {T42["dfz41_relabel_min"]}; DFZ 3.1 re-expanded - '
         f'printed {T42["dfz31_difference"]}; forms, file and content match the manifest {file42 and T42["sha"] == M42["forms"]}')

    # G43 (R47): stabilizer secret sharing (docs/R32_reduction_and_local_dimension.md sections 4.15 and 4.16).  The
    # bound 4/3 of Proposition 21 is a proof from the results of Jafari and Khazaei; the gate checks the combinatorial
    # facts it uses: contracting one player turns Gamma_F into their Fano access structure F (for every player, with
    # the size-3 minimal sets equal to the Fano lines that miss it) and Gamma_NF into their non-Fano access structure N
    # (for A, B, C, Z, with the triple 456 of N equal to the other three), the minimal authorised sets of Gamma_F are
    # the lines, Gamma_NF is the same for p = 3, 5, 7, and the two access structures differ exactly in WXY and ABCZ.
    # The bound 43/41 for every odd local dimension: the rational dual certificates of the two linear programs (Lemma 4
    # (LB), t = 2, on all seven-element minors; pure and with an environment E) rebuilt row by row from their
    # descriptors, checked to be admissible rows and verified in exact arithmetic; (LB) is negative on the qubit cube
    # state on some minor and (LA) on the odd one, each form nonnegative on the other cube state over all minors and
    # labellings; the margins of section 4.16 (1/952, 1/129, 1/45, 2/105).
    # G43b (--full): the programs and two controls (optimum 1) re-solved with HiGHS, with an exact feasible point
    # attaining 43/41; the programs with the Jafari-Khazaei inequalities on contractions: optimum 4/3, exact feasible
    # points (the inequalities themselves imply t >= 4/3).
    import s38_secret_sharing as _s43
    T43 = _s43.check(log=lambda *_x: None)
    M43 = MAN['s38_qss_certificates']
    file43 = (_h33.sha256(open(os.path.join(core.DATA, 's38_qss_certificates.json'), 'rb').read()).hexdigest()
              == M43['file_sha256'] and _s43.content_sha() == M43['sha256'])
    cases43 = {k: v['certified'] for k, v in T43['cases'].items()}
    rows43 = {k: v['rows'] for k, v in T43['cases'].items()}
    ports43 = {q: v.get('456') for q, v in T43['ports']['NF'].items()}
    gate('G43 secret sharing: contractions, certificates, controls, margins',
         T43['ok'] and file43 and cases43 == M43['bounds'] and rows43 == M43['rows'] and T43['controls'] == M43['controls']
         and ports43 == M43['ports_456'] and {k: v['eps_inf'] for k, v in T43['margins'].items()} == M43['margins'],
         f'Gamma_F/q = F (Jafari-Khazaei) for all q and Gamma_NF/q = N for q in ABCZ: {_s43.ports_ok(T43["ports"])} '
         f'(456 of N: {ports43}); access structures differ in {T43["access"]["F_only"]} / {T43["access"]["NF_only"]}; '
         f'certified bounds {cases43} with rows {rows43}; minima over all minors {T43["controls"]}; margins '
         + ', '.join(f'{k} {v["eps_inf"]}' for k, v in T43['margins'].items())
         + f'; file and content match the manifest {file43}')

    # G44 (R49): undecided orbits realised after submission (results/2026-11-17_r49-sandbox).  data/ is unchanged, so
    # the gates above and the scripts in checks/ still see the 641 orbits of the paper.  Each of the 141 stored
    # certificates is a qubit graph state whose 31 cut ranks over GF(2), recomputed with epr1kit.stabcert, equal lam
    # times its row of the R30 list results/2026-10-02_a2-status/undecided_reps.npy; every such row satisfies the 10,860
    # facets of QLR5 with a tight set of rank 30 (an extreme ray of QLR5, hence of Stab5) and lies in an S6 orbit
    # different from those of data/stab5_extreme_reps.npy and of the other 140 rows; undecided_reps.npy of the folder
    # is the rest of the R30 list in the same order, and status.json gives 782 known / 1,014 undecided.  The unimodular
    # realisations of every_d_sandbox.npz are exact (scripts/s32_unimodular.py check, catalogue extended by the 141
    # rows as with --append), so those rows lie in Stab^(Z_d)_5 for every d.  File hashes match the folder manifest.
    base44 = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results')
    D44 = os.path.join(base44, '2026-11-17_r49-sandbox')
    M44 = _js20.load(open(os.path.join(D44, 'manifest.json')))
    file44 = all(_h33.sha256(open(os.path.join(D44, f), 'rb').read()).hexdigest() == e['sha256']
                 for f, e in M44['files'].items()) and len(M44['files']) == 6
    R44 = np.load(os.path.join(D44, 'realisations_141.npz'), allow_pickle=False)
    U44 = np.load(os.path.join(base44, '2026-10-02_a2-status', 'undecided_reps.npy')).astype(np.int64)
    L44 = np.load(os.path.join(D44, 'undecided_reps.npy')).astype(np.int64)
    rows44 = [int(x) for x in R44['row']]
    keep44 = np.ones(U44.shape[0], dtype=bool)
    keep44[rows44] = False
    ok44a = (len(rows44) == len(set(rows44)) == 141 and np.array_equal(R44['rep'].astype(np.int64), U44[rows44])
             and np.array_equal(L44, U44[keep44]))
    nv44 = sum(_sc.verify_realisation(U44[r], int(R44['lam'][k]), R44['sizes'][k], R44['rows'][k])
               for k, r in enumerate(rows44))
    HF44 = core.load('qlr5_H_facets10860').astype(np.int64)
    V44 = HF44 @ U44[rows44].T
    ne44 = sum(1 for k in range(len(rows44)) if V44[:, k].min() >= 0 and core.rank_mod_p(HF44[V44[:, k] == 0]) == 30)
    c44, _ = core.class_reps(np.vstack([core.load('stab5_extreme_reps').astype(np.int64), U44[rows44]]),
                             core.perms31_s6(), offset=4096, width='u16')
    no44 = len(set(c44[641:])) == 141 and not (set(c44[641:]) & set(c44[:641]))
    ST44 = _js20.load(open(os.path.join(D44, 'status.json')))
    acc44 = (ST44['known_extreme_ray_orbits'] == 641 + len(rows44) == 782 and ST44['undecided'] == L44.shape[0] == 1014
             and ST44['realised_rows'] == rows44)
    EV44 = np.load(os.path.join(D44, 'every_d_sandbox.npz'), allow_pickle=False)
    _s36.APPEND = os.path.join(D44, 'realisations_141.npz')
    try:
        S44 = _s36.check(os.path.join(D44, 'every_d_sandbox.npz'), full=False, log=lambda *_x: None)
    finally:
        _s36.APPEND = None
    ok44e = (S44['ok'] == S44['stored'] == S44['distinct_rows'] == len(ST44['every_d']['rows']) and S44['rep_ok']
             and all(641 <= int(k) < 782 for k in EV44['row'])
             and sorted(rows44[int(k) - 641] for k in EV44['row']) == ST44['every_d']['rows'])
    lam44 = {int(l): int((R44['lam'] == l).sum()) for l in sorted(set(R44['lam'].tolist()))}
    gate('G44 undecided orbits realised after submission', ok44a and nv44 == 141 and ne44 == 141 and no44 and acc44
         and ok44e and file44,
         f'{nv44}/141 certificates re-verified (GF(2) cut ranks, lambda {lam44}); {ne44}/141 extreme rays of QLR5 '
         f'(tight rank 30); new S6 orbits {no44}; {ST44["known_extreme_ray_orbits"]} known / {ST44["undecided"]} '
         f'undecided (data/ keeps the 641 of the paper); every d: {S44["ok"]}/{S44["stored"]} unimodular realisations '
         f'exact; files match the folder manifest {file44}')

    # G45 (R50): the second batch, results/2026-11-18_r50-sandbox, and the certificates for every d of the job s32-new.
    # As in G44: each of the 126 certificates (the job a2-long2) re-verified with epr1kit.stabcert against its row of the
    # R30 list (and of the R49 list that a2-long2 searched); every row an extreme ray of QLR5 with tight rank 30; S6
    # orbits new with respect to data/stab5_extreme_reps.npy, the 141 rows of G44 and each other; undecided_reps.npy the
    # rest of the R30 list; status.json 908 known / 888 undecided.  The unimodular realisations of
    # results/2026-11-18_s32-new/s32n_unimodular.npz (rows: the catalogue extended by the 141 rows) are exact, include
    # those of G44, and match status.json; the hashes of both job folders and of this folder match their manifests.
    D45 = os.path.join(base44, '2026-11-18_r50-sandbox')
    file45 = True
    for d45, n45 in ((D45, 5), (os.path.join(base44, '2026-11-18_s32-new'), 4),
                     (os.path.join(base44, '2026-11-18_a2-long2'), 8)):
        M45 = _js20.load(open(os.path.join(d45, 'manifest.json')))
        file45 = file45 and len(M45['files']) == n45 and all(
            _h33.sha256(open(os.path.join(d45, f), 'rb').read()).hexdigest() == e['sha256'] for f, e in M45['files'].items())
    R45 = np.load(os.path.join(D45, 'realisations_126.npz'), allow_pickle=False)
    rows45 = [int(x) for x in R45['row']]
    U49 = np.load(os.path.join(D44, 'undecided_reps.npy')).astype(np.int64)
    L45 = np.load(os.path.join(D45, 'undecided_reps.npy')).astype(np.int64)
    keep45 = keep44.copy()
    keep45[rows45] = False
    ok45a = (len(rows45) == len(set(rows45)) == 126 and not set(rows45) & set(rows44)
             and np.array_equal(R45['rep'].astype(np.int64), U44[rows45])
             and np.array_equal(U49[[int(x) for x in R45['row49']]], U44[rows45])
             and np.array_equal(L45, U44[keep45]))
    nv45 = sum(_sc.verify_realisation(U44[r], int(R45['lam'][k]), R45['sizes'][k], R45['rows'][k])
               for k, r in enumerate(rows45))
    V45 = HF44 @ U44[rows45].T
    ne45 = sum(1 for k in range(len(rows45)) if V45[:, k].min() >= 0 and core.rank_mod_p(HF44[V45[:, k] == 0]) == 30)
    c45, _ = core.class_reps(U44[rows45], core.perms31_s6(), offset=4096, width='u16')
    no45 = len(set(c45)) == 126 and not (set(c45) & set(c44))
    ST45 = _js20.load(open(os.path.join(D45, 'status.json')))
    acc45 = (ST45['known_extreme_ray_orbits'] == 641 + 141 + 126 == 908 and ST45['undecided'] == L45.shape[0] == 888
             and ST45['realised_rows_r50'] == rows45)
    P45 = os.path.join(base44, '2026-11-18_s32-new', 's32n_unimodular.npz')
    EV45 = np.load(P45, allow_pickle=False)
    _s36.APPEND = os.path.join(D44, 'realisations_141.npz')
    try:
        S45 = _s36.check(P45, full=False, log=lambda *_x: None)
    finally:
        _s36.APPEND = None
    rows45e = sorted(rows44[int(k) - 641] for k in EV45['row'])
    ok45e = (S45['ok'] == S45['stored'] == S45['distinct_rows'] == ST45['every_d']['unimodular'] and S45['rep_ok']
             and all(641 <= int(k) < 782 for k in EV45['row']) and rows45e == ST45['every_d']['rows']
             and set(int(k) for k in EV44['row']) <= set(int(k) for k in EV45['row']))
    lam45 = {int(l): int((R45['lam'] == l).sum()) for l in sorted(set(R45['lam'].tolist()))}
    gate('G45 second batch realised after submission', ok45a and nv45 == 126 and ne45 == 126 and no45 and acc45
         and ok45e and file45,
         f'{nv45}/126 certificates re-verified (GF(2) cut ranks, lambda {lam45}); {ne45}/126 extreme rays of QLR5 '
         f'(tight rank 30); new S6 orbits {no45}; {ST45["known_extreme_ray_orbits"]} known / {ST45["undecided"]} '
         f'undecided; every d: {S45["ok"]}/{S45["stored"]} unimodular realisations exact (s32-new; kappa '
         f'{S45["kappa"]}); files match the manifests {file45}')

    # G46 (R51): the third batch, results/2026-11-18_r51-sandbox, and the certificates for every d of the job s32-new2.
    # As in G45: each of the 73 certificates (the job a2-long3) re-verified with epr1kit.stabcert against its row of the
    # R30 list (and of the R50 list that a2-long3 searched); every row an extreme ray of QLR5 with tight rank 30; S6
    # orbits new with respect to data/stab5_extreme_reps.npy, the rows of G44 and G45 and each other; undecided_reps.npy
    # the rest of the R30 list; status.json 981 known / 815 undecided.  The unimodular realisations of
    # results/2026-11-18_s32-new2/s32n2_unimodular.npz (rows: the catalogue extended by the 141 rows and then the 126)
    # are exact, include those of G45, and match status.json; the hashes of the three new folders match their manifests.
    D46 = os.path.join(base44, '2026-11-18_r51-sandbox')
    file46 = True
    for d46, n46 in ((D46, 5), (os.path.join(base44, '2026-11-18_s32-new2'), 4),
                     (os.path.join(base44, '2026-11-18_a2-long3'), 8)):
        M46 = _js20.load(open(os.path.join(d46, 'manifest.json')))
        file46 = file46 and len(M46['files']) == n46 and all(
            _h33.sha256(open(os.path.join(d46, f), 'rb').read()).hexdigest() == e['sha256'] for f, e in M46['files'].items())
    R46 = np.load(os.path.join(D46, 'realisations_73.npz'), allow_pickle=False)
    rows46 = [int(x) for x in R46['row']]
    L46 = np.load(os.path.join(D46, 'undecided_reps.npy')).astype(np.int64)
    keep46 = keep45.copy()
    keep46[rows46] = False
    ok46a = (len(rows46) == len(set(rows46)) == 73 and not set(rows46) & (set(rows44) | set(rows45))
             and np.array_equal(R46['rep'].astype(np.int64), U44[rows46])
             and np.array_equal(L45[[int(x) for x in R46['row50']]], U44[rows46])
             and np.array_equal(L46, U44[keep46]))
    nv46 = sum(_sc.verify_realisation(U44[r], int(R46['lam'][k]), R46['sizes'][k], R46['rows'][k])
               for k, r in enumerate(rows46))
    V46 = HF44 @ U44[rows46].T
    ne46 = sum(1 for k in range(len(rows46)) if V46[:, k].min() >= 0 and core.rank_mod_p(HF44[V46[:, k] == 0]) == 30)
    c46, _ = core.class_reps(U44[rows46], core.perms31_s6(), offset=4096, width='u16')
    no46 = len(set(c46)) == 73 and not (set(c46) & (set(c44) | set(c45)))
    ST46 = _js20.load(open(os.path.join(D46, 'status.json')))
    acc46 = (ST46['known_extreme_ray_orbits'] == 641 + 141 + 126 + 73 == 981 and ST46['undecided'] == L46.shape[0] == 815
             and ST46['realised_rows_r51'] == rows46)
    P46 = os.path.join(base44, '2026-11-18_s32-new2', 's32n2_unimodular.npz')
    EV46 = np.load(P46, allow_pickle=False)
    _s36.APPEND = [os.path.join(D44, 'realisations_141.npz'), os.path.join(D45, 'realisations_126.npz')]
    try:
        S46 = _s36.check(P46, full=False, log=lambda *_x: None)
    finally:
        _s36.APPEND = None
    rows46e = sorted(rows44[int(k) - 641] if int(k) < 782 else rows45[int(k) - 782] for k in EV46['row'])
    ok46e = (S46['ok'] == S46['stored'] == S46['distinct_rows'] == ST46['every_d']['unimodular'] and S46['rep_ok']
             and all(641 <= int(k) < 908 for k in EV46['row']) and rows46e == ST46['every_d']['rows']
             and set(int(k) for k in EV45['row']) <= set(int(k) for k in EV46['row']))
    lam46 = {int(l): int((R46['lam'] == l).sum()) for l in sorted(set(R46['lam'].tolist()))}
    gate('G46 third batch realised after submission', ok46a and nv46 == 73 and ne46 == 73 and no46 and acc46
         and ok46e and file46,
         f'{nv46}/73 certificates re-verified (GF(2) cut ranks, lambda {lam46}); {ne46}/73 extreme rays of QLR5 '
         f'(tight rank 30); new S6 orbits {no46}; {ST46["known_extreme_ray_orbits"]} known / {ST46["undecided"]} '
         f'undecided; every d: {S46["ok"]}/{S46["stored"]} unimodular realisations exact (s32-new2; kappa '
         f'{S46["kappa"]}); files match the manifests {file46}')

    # G47 (R52): the fourth batch, results/2026-11-19_r52-sandbox: five rows of the R30 list realised at lambda = 2 by
    # the CSS search (scripts/s41_css_search.py).  Each certificate is checked twice: the stored bipartite graph with
    # epr1kit.stabcert, and the stored matrix M with s41's verify (its 63 column-group ranks equal h_S, so the CSS state
    # of its row space has entropy 2S; the graph built from M has cut ranks 2r; stabcert on that graph).  As in G46:
    # rows of the R30 list (and of the R51 list) not realised before; extreme rays of QLR5 with tight rank 30; S6 orbits
    # new with respect to all 981 known before and each other; undecided_reps.npy the rest of the R30 list; status.json
    # 986 known / 810 undecided; the frontier's cut matrix is covered by its greedy cover; hashes match the manifest.
    import s41_css_search as _css47
    D47 = os.path.join(base44, '2026-11-19_r52-sandbox')
    M47 = _js20.load(open(os.path.join(D47, 'manifest.json')))
    file47 = len(M47['files']) == 7 and all(
        _h33.sha256(open(os.path.join(D47, f), 'rb').read()).hexdigest() == e['sha256'] for f, e in M47['files'].items())
    R47 = np.load(os.path.join(D47, 'realisations_5.npz'), allow_pickle=False)
    rows47 = [int(x) for x in R47['row']]
    L47 = np.load(os.path.join(D47, 'undecided_reps.npy')).astype(np.int64)
    keep47 = keep46.copy()
    keep47[rows47] = False
    ok47a = (len(rows47) == len(set(rows47)) == 5 and not set(rows47) & (set(rows44) | set(rows45) | set(rows46))
             and np.array_equal(R47['rep'].astype(np.int64), U44[rows47])
             and np.array_equal(L46[[int(x) for x in R47['row51']]], U44[rows47])
             and np.array_equal(L47, U44[keep47]) and (R47['lam'] == 2).all())
    nv47 = sum(_sc.verify_realisation(U44[r], int(R47['lam'][k]), R47['sizes'][k], R47['rows'][k])
               for k, r in enumerate(rows47))
    nc47 = sum(_css47.verify(U44[r], 1, [int(c) for c in R47['cols'][k][:int(R47['m'][k])]])[0] for k, r in enumerate(rows47))
    V47 = HF44 @ U44[rows47].T
    ne47 = sum(1 for k in range(len(rows47)) if V47[:, k].min() >= 0 and core.rank_mod_p(HF44[V47[:, k] == 0]) == 30)
    c47, _ = core.class_reps(U44[rows47], core.perms31_s6(), offset=4096, width='u16')
    no47 = len(set(c47)) == 5 and not (set(c47) & (set(c44) | set(c45) | set(c46)))
    ST47 = _js20.load(open(os.path.join(D47, 'status.json')))
    acc47 = (ST47['known_extreme_ray_orbits'] == 641 + 141 + 126 + 73 + 5 == 986 and ST47['undecided'] == L47.shape[0] == 810
             and ST47['realised_rows_r52'] == rows47 and ST47['every_d']['open'] == 259)
    F47 = np.load(os.path.join(D47, 'frontier_986.npz'), allow_pickle=False)
    fr47 = (F47['cut'].shape == (F47['classes'].shape[0], 810) and F47['cut'].any(0).all()
            and F47['cut'][F47['greedy_cover']].any(0).all())
    gate('G47 fourth batch (CSS search)', ok47a and nv47 == 5 and nc47 == 5 and ne47 == 5 and no47 and acc47 and fr47
         and file47,
         f'{nv47}/5 graph certificates re-verified at lambda 2 (GF(2) cut ranks); {nc47}/5 CSS matrices re-verified '
         f'(column-group ranks = h_S, graph cut ranks, stabcert); {ne47}/5 extreme rays of QLR5 (tight rank 30); new S6 '
         f'orbits {no47}; {ST47["known_extreme_ray_orbits"]} known / {ST47["undecided"]} undecided; frontier '
         f'{F47["classes"].shape[0]} classes, greedy cover {len(F47["greedy_cover"])} covers all rows {fr47}; files match '
         f'the manifest {file47}')

    # G48 (R53): the fifth batch, results/2026-11-20_r53-sandbox: nine rows of the R30 list realised at lambda = 2 by the
    # job css-sweep (results/2026-11-20_css-sweep), and the certificates for every d of the job s32-new3
    # (results/2026-11-20_s32-new3).  As in G47: each certificate checked twice (the stored graph with epr1kit.stabcert,
    # the stored matrix M with s41's verify) and equal to the hit stored by the job; rows of the R30 list (and of the R52
    # list that css-sweep searched) not realised before; extreme rays of QLR5 with tight rank 30; S6 orbits new with
    # respect to all 986 known before and each other; undecided_reps.npy the rest of the R30 list; status.json 995 known /
    # 801 undecided.  The unimodular realisations of s32n3_unimodular.npz (rows: the catalogue extended by the four
    # realisation files of R49-R52) are exact, include those of s32-new2 and match status.json, open rows included; the
    # frontier's cut matrix is covered by its greedy cover; the hashes of the three new folders match their manifests.
    D48 = os.path.join(base44, '2026-11-20_r53-sandbox')
    file48 = True
    for d48, n48 in ((D48, 5), (os.path.join(base44, '2026-11-20_css-sweep'), 5),
                     (os.path.join(base44, '2026-11-20_s32-new3'), 4)):
        M48 = _js20.load(open(os.path.join(d48, 'manifest.json')))
        file48 = file48 and len(M48['files']) == n48 and all(
            _h33.sha256(open(os.path.join(d48, f), 'rb').read()).hexdigest() == e['sha256'] for f, e in M48['files'].items())
    R48 = np.load(os.path.join(D48, 'realisations_9.npz'), allow_pickle=False)
    rows48 = [int(x) for x in R48['row']]
    L48 = np.load(os.path.join(D48, 'undecided_reps.npy')).astype(np.int64)
    keep48 = keep47.copy()
    keep48[rows48] = False
    CS48 = np.load(os.path.join(base44, '2026-11-20_css-sweep', 'css_sweep.npz'), allow_pickle=False)
    j48 = {int(x): i for i, x in enumerate(CS48['index'])}
    same48 = (sorted(j48) == sorted(int(x) for x in R48['row52']) and all(
        np.array_equal(CS48['cols'][j48[int(q)]], R48['cols'][k]) and np.array_equal(CS48['rows'][j48[int(q)]], R48['rows'][k])
        and np.array_equal(CS48['sizes'][j48[int(q)]], R48['sizes'][k]) for k, q in enumerate(R48['row52'])))
    ok48a = (len(rows48) == len(set(rows48)) == 9
             and not set(rows48) & (set(rows44) | set(rows45) | set(rows46) | set(rows47))
             and np.array_equal(R48['rep'].astype(np.int64), U44[rows48])
             and np.array_equal(L47[[int(x) for x in R48['row52']]], U44[rows48])
             and np.array_equal(L48, U44[keep48]) and (R48['lam'] == 2).all() and same48)
    nv48 = sum(_sc.verify_realisation(U44[r], int(R48['lam'][k]), R48['sizes'][k], R48['rows'][k])
               for k, r in enumerate(rows48))
    nc48 = sum(_css47.verify(U44[r], 1, [int(c) for c in R48['cols'][k][:int(R48['m'][k])]])[0] for k, r in enumerate(rows48))
    V48 = HF44 @ U44[rows48].T
    ne48 = sum(1 for k in range(len(rows48)) if V48[:, k].min() >= 0 and core.rank_mod_p(HF44[V48[:, k] == 0]) == 30)
    c48, _ = core.class_reps(U44[rows48], core.perms31_s6(), offset=4096, width='u16')
    no48 = len(set(c48)) == 9 and not (set(c48) & (set(c44) | set(c45) | set(c46) | set(c47)))
    ST48 = _js20.load(open(os.path.join(D48, 'status.json')))
    acc48 = (ST48['known_extreme_ray_orbits'] == 641 + 141 + 126 + 73 + 5 + 9 == 995
             and ST48['undecided'] == L48.shape[0] == 801 and ST48['realised_rows_r53'] == rows48)
    P48 = os.path.join(base44, '2026-11-20_s32-new3', 's32n3_unimodular.npz')
    EV48 = np.load(P48, allow_pickle=False)
    _s36.APPEND = [os.path.join(D44, 'realisations_141.npz'), os.path.join(D45, 'realisations_126.npz'),
                   os.path.join(D46, 'realisations_73.npz'), os.path.join(D47, 'realisations_5.npz')]
    try:
        S48 = _s36.check(P48, full=False, log=lambda *_x: None)
    finally:
        _s36.APPEND = None
    cat48 = rows44 + rows45 + rows46 + rows47
    rows48e = sorted(cat48[int(k) - 641] for k in EV48['row'])
    open48 = sorted((set(cat48) | set(rows48)) - set(rows48e))
    ok48e = (S48['ok'] == S48['stored'] == S48['distinct_rows'] == ST48['every_d']['unimodular'] == 101 and S48['rep_ok']
             and all(641 <= int(k) < 986 for k in EV48['row']) and rows48e == ST48['every_d']['rows']
             and open48 == ST48['every_d']['open_rows'] and ST48['every_d']['open'] == len(open48) == 253
             and set(int(k) for k in EV46['row']) <= set(int(k) for k in EV48['row']))
    F48 = np.load(os.path.join(D48, 'frontier_995.npz'), allow_pickle=False)
    fr48 = (F48['cut'].shape == (F48['classes'].shape[0], 801) and F48['cut'].any(0).all()
            and F48['cut'][F48['greedy_cover']].any(0).all())
    gate('G48 fifth batch (css-sweep) and s32-new3', ok48a and nv48 == 9 and nc48 == 9 and ne48 == 9 and no48 and acc48
         and ok48e and fr48 and file48,
         f'{nv48}/9 graph certificates re-verified at lambda 2 (GF(2) cut ranks); {nc48}/9 CSS matrices re-verified '
         f'(column-group ranks = h_S, graph cut ranks, stabcert), equal to the hits of css-sweep {same48}; {ne48}/9 '
         f'extreme rays of QLR5 (tight rank 30); new S6 orbits {no48}; {ST48["known_extreme_ray_orbits"]} known / '
         f'{ST48["undecided"]} undecided; every d: {S48["ok"]}/{S48["stored"]} unimodular realisations exact '
         f'(s32-new3; kappa {S48["kappa"]}), {len(open48)} open; frontier {F48["classes"].shape[0]} classes, greedy '
         f'cover {len(F48["greedy_cover"])} covers all rows {fr48}; files match the manifests {file48}')

    # G49 (R54): the sixth batch, results/2026-11-23_r54-sandbox: eight rows of the R30 list realised at lambda = 2 by the
    # job css-sweep2 (results/2026-11-23_css-sweep2), and the certificates for every d of the job s32-new4
    # (results/2026-11-23_s32-new4).  The checks of G48, one batch later: each certificate checked twice and equal to the
    # hit stored by the job; rows of the R30 list (and of the R53 list that css-sweep2 searched) not realised before;
    # extreme rays of QLR5 with tight rank 30; S6 orbits new with respect to all 995 known before and each other;
    # undecided_reps.npy the rest of the R30 list; status.json 1003 known / 793 undecided.  The unimodular realisations of
    # s32n4_unimodular.npz (rows: the catalogue extended by the five realisation files of R49-R53) are exact, include
    # those of s32-new3 and match status.json, open rows included; the frontier's greedy cover; the three manifests.
    D49 = os.path.join(base44, '2026-11-23_r54-sandbox')
    file49 = True
    for d49, n49 in ((D49, 5), (os.path.join(base44, '2026-11-23_css-sweep2'), 5),
                     (os.path.join(base44, '2026-11-23_s32-new4'), 4)):
        M49 = _js20.load(open(os.path.join(d49, 'manifest.json')))
        file49 = file49 and len(M49['files']) == n49 and all(
            _h33.sha256(open(os.path.join(d49, f), 'rb').read()).hexdigest() == e['sha256'] for f, e in M49['files'].items())
    R49 = np.load(os.path.join(D49, 'realisations_8.npz'), allow_pickle=False)
    rows49 = [int(x) for x in R49['row']]
    L49 = np.load(os.path.join(D49, 'undecided_reps.npy')).astype(np.int64)
    keep49 = keep48.copy()
    keep49[rows49] = False
    CS49 = np.load(os.path.join(base44, '2026-11-23_css-sweep2', 'css_sweep2.npz'), allow_pickle=False)
    j49 = {int(x): i for i, x in enumerate(CS49['index'])}
    same49 = (sorted(j49) == sorted(int(x) for x in R49['row53']) and all(
        np.array_equal(CS49['cols'][j49[int(q)]], R49['cols'][k]) and np.array_equal(CS49['rows'][j49[int(q)]], R49['rows'][k])
        and np.array_equal(CS49['sizes'][j49[int(q)]], R49['sizes'][k]) for k, q in enumerate(R49['row53'])))
    ok49a = (len(rows49) == len(set(rows49)) == 8
             and not set(rows49) & (set(rows44) | set(rows45) | set(rows46) | set(rows47) | set(rows48))
             and np.array_equal(R49['rep'].astype(np.int64), U44[rows49])
             and np.array_equal(L48[[int(x) for x in R49['row53']]], U44[rows49])
             and np.array_equal(L49, U44[keep49]) and (R49['lam'] == 2).all() and same49)
    nv49 = sum(_sc.verify_realisation(U44[r], int(R49['lam'][k]), R49['sizes'][k], R49['rows'][k])
               for k, r in enumerate(rows49))
    nc49 = sum(_css47.verify(U44[r], 1, [int(c) for c in R49['cols'][k][:int(R49['m'][k])]])[0] for k, r in enumerate(rows49))
    V49 = HF44 @ U44[rows49].T
    ne49 = sum(1 for k in range(len(rows49)) if V49[:, k].min() >= 0 and core.rank_mod_p(HF44[V49[:, k] == 0]) == 30)
    c49, _ = core.class_reps(U44[rows49], core.perms31_s6(), offset=4096, width='u16')
    no49 = len(set(c49)) == 8 and not (set(c49) & (set(c44) | set(c45) | set(c46) | set(c47) | set(c48)))
    ST49 = _js20.load(open(os.path.join(D49, 'status.json')))
    acc49 = (ST49['known_extreme_ray_orbits'] == 641 + 141 + 126 + 73 + 5 + 9 + 8 == 1003
             and ST49['undecided'] == L49.shape[0] == 793 and ST49['realised_rows_r54'] == rows49)
    P49 = os.path.join(base44, '2026-11-23_s32-new4', 's32n4_unimodular.npz')
    EV49 = np.load(P49, allow_pickle=False)
    _s36.APPEND = [os.path.join(D44, 'realisations_141.npz'), os.path.join(D45, 'realisations_126.npz'),
                   os.path.join(D46, 'realisations_73.npz'), os.path.join(D47, 'realisations_5.npz'),
                   os.path.join(D48, 'realisations_9.npz')]
    try:
        S49 = _s36.check(P49, full=False, log=lambda *_x: None)
    finally:
        _s36.APPEND = None
    cat49 = cat48 + rows48
    rows49e = sorted(cat49[int(k) - 641] for k in EV49['row'])
    open49 = sorted((set(cat49) | set(rows49)) - set(rows49e))
    ok49e = (S49['ok'] == S49['stored'] == S49['distinct_rows'] == ST49['every_d']['unimodular'] == 111 and S49['rep_ok']
             and all(641 <= int(k) < 995 for k in EV49['row']) and rows49e == ST49['every_d']['rows']
             and open49 == ST49['every_d']['open_rows'] and ST49['every_d']['open'] == len(open49) == 251
             and set(int(k) for k in EV48['row']) <= set(int(k) for k in EV49['row']))
    F49 = np.load(os.path.join(D49, 'frontier_1003.npz'), allow_pickle=False)
    fr49 = (F49['cut'].shape == (F49['classes'].shape[0], 793) and F49['cut'].any(0).all()
            and F49['cut'][F49['greedy_cover']].any(0).all())
    gate('G49 sixth batch (css-sweep2) and s32-new4', ok49a and nv49 == 8 and nc49 == 8 and ne49 == 8 and no49 and acc49
         and ok49e and fr49 and file49,
         f'{nv49}/8 graph certificates re-verified at lambda 2 (GF(2) cut ranks); {nc49}/8 CSS matrices re-verified '
         f'(column-group ranks = h_S, graph cut ranks, stabcert), equal to the hits of css-sweep2 {same49}; {ne49}/8 '
         f'extreme rays of QLR5 (tight rank 30); new S6 orbits {no49}; {ST49["known_extreme_ray_orbits"]} known / '
         f'{ST49["undecided"]} undecided; every d: {S49["ok"]}/{S49["stored"]} unimodular realisations exact '
         f'(s32-new4; kappa {S49["kappa"]}), {len(open49)} open; frontier {F49["classes"].shape[0]} classes, greedy '
         f'cover {len(F49["greedy_cover"])} covers all rows {fr49}; files match the manifests {file49}')

    # G50 (R56): Gate 1, results/2026-11-23_r56-sandbox.  (a) The rows of the R54 list that the job ci2-ingleton
    # (results/2026-11-23_ci2-ingleton; scripts/s42_ci2_ingleton.py: two common informations and the set-valued
    # conditional Ingleton inequalities on the 8-element extension, floating point) reported infeasible: each has an
    # exact certificate (scripts/s43_certify_ci2.py), re-checked here by epr1kit.stabcert.verify_exclusion_ci2, which
    # rebuilds every Shannon, Ingleton and common-information row from its descriptor, on the ray of the R54 list; the
    # certificates' (X1, Y1), (X2, Y2) are those of the job.  (b) The qubit realisations of the fixed-structure search
    # (scripts/s44_fixed_search.py and its sandbox prototype): the checks of G49 (s41 verify, stabcert on the stored
    # graph, extreme rays of QLR5 with tight rank 30, S6 orbits new against all 1,003 known before and each other).
    # (c) The qutrit realisations (GF(3), lambda = 2): the 63 column-group ranks of the stored matrix over GF(3)
    # (stabcert.gfp_rank) and the stored weighted qutrit graph (stabcert.verify_realisation_gfp).  (d) No row is both
    # realised and excluded; undecided_reps.npy is the rest of the R54 list (the qutrit rows stay in it); status.json
    # counts.  (e) The frontier: its rows are the undecided rows, its greedy cover cuts every one.  (f) The manifests.
    D50 = os.path.join(base44, '2026-11-23_r56-sandbox')
    J50d = os.path.join(base44, '2026-11-23_ci2-ingleton')
    ST50 = _js20.load(open(os.path.join(D50, 'status.json')))
    RN50, FN50 = ST50['realisations_file'], ST50['frontier']['file']
    file50 = True
    for d50, f50 in ((D50, ['README.txt', 'certificates_ci2.json', RN50, 'realisations_gf3.json', 'undecided_reps.npy',
                            'status.json', FN50]),
                     (J50d, ['ci2_ingleton.json', 'ci2-ingleton.log', 'JOB_STATUS.json'])):
        M50 = _js20.load(open(os.path.join(d50, 'manifest.json')))
        file50 = file50 and sorted(M50['files']) == sorted(f50) and all(
            _h33.sha256(open(os.path.join(d50, f), 'rb').read()).hexdigest() == e['sha256'] for f, e in M50['files'].items())
    from fractions import Fraction as _Fr50
    J50 = _js20.load(open(os.path.join(J50d, 'ci2_ingleton.json')))
    C50 = _js20.load(open(os.path.join(D50, 'certificates_ci2.json')))
    cand50 = sorted(int(k) for k, v in J50['results'].items() if v['infeasible'])
    rows50x = sorted(int(c['index']) for c in C50['certificates'])
    okx50 = 0
    for c in C50['certificates']:
        r = L49[int(c['index'])]
        okc, _v = _sc.verify_exclusion_ci2(r, c['P1'], c['P2'], [(tuple(t[0]), _Fr50(t[1])) for t in c['terms']],
                                           F_S=c['F_S'])
        inf = J50['results'][str(c['index'])]['infeasible']
        okx50 += bool(okc and [int(x) for x in r] == c['ray'] and list(c['P1']) == list(inf[0])
                      and list(c['P2']) == list(inf[1])
                      and sum(a * int(b) for a, b in zip(c['F_S'], r)) == c['F_S_at_ray'] < 0)
    ex50 = (J50['rays'] == 'results/2026-11-23_r54-sandbox/undecided_reps.npy' and len(set(rows50x)) == len(rows50x)
            and sorted(rows50x + [int(x) for x in C50['failed']]) == cand50 and okx50 == len(rows50x))
    R50 = np.load(os.path.join(D50, RN50), allow_pickle=False)
    rows50 = [int(x) for x in R50['row']]
    i49 = np.flatnonzero(keep49)
    ok50a = (len(rows50) == len(set(rows50)) > 0 and keep49[rows50].all()
             and [int(i49[int(x)]) for x in R50['row54']] == rows50
             and np.array_equal(R50['rep'].astype(np.int64), U44[rows50]) and (R50['lam'] == 2).all())
    nv50 = sum(_sc.verify_realisation(U44[r], 2, R50['sizes'][k], R50['rows'][k]) for k, r in enumerate(rows50))
    nc50 = sum(_css47.verify(U44[r], 1, [int(c) for c in R50['cols'][k][:int(R50['m'][k])]])[0] for k, r in enumerate(rows50))
    V50 = HF44 @ U44[rows50].T
    ne50 = sum(1 for k in range(len(rows50)) if V50[:, k].min() >= 0 and core.rank_mod_p(HF44[V50[:, k] == 0]) == 30)
    c50, _ = core.class_reps(U44[rows50], core.perms31_s6(), offset=4096, width='u16')
    no50 = len(set(c50)) == len(rows50) and not (set(c50) & (set(c44) | set(c45) | set(c46) | set(c47) | set(c48) | set(c49)))
    Q50 = _js20.load(open(os.path.join(D50, 'realisations_gf3.json')))
    nq50 = 0
    for q in Q50:
        r = U44[q['row']]
        hq = _sc.h_norm(r)
        off = [0]
        for i in range(6):
            off.append(off[-1] + int(hq[1 << i]))
        cols = [tuple(int(x, 16) for x in s.split(':')) for s in q['cols']]
        Mq = [[1 if pp >> i & 1 else 2 if qq >> i & 1 else 0 for (pp, qq) in cols] for i in range(int(hq[63]))]
        okr = len(cols) == off[6] and all(not pp & qq for pp, qq in cols) and all(
            _sc.gfp_rank([[row[c] for j in range(6) if X >> j & 1 for c in range(off[j], off[j + 1])] for row in Mq], 3)
            == hq[X] for X in range(1, 64))
        nq50 += bool(okr and keep49[q['row']] and int(i49[q['row54']]) == q['row'] and q['p'] == 3 and q['lam'] == 2
                     and _sc.verify_realisation_gfp(r, 2, q['sizes'], q['graph_weights'], 3))
    x30_50 = [int(i49[k]) for k in rows50x]
    q30_50 = [int(q['row']) for q in Q50]
    keep50 = keep49.copy()
    keep50[rows50] = False
    keep50[x30_50] = False
    L50 = np.load(os.path.join(D50, 'undecided_reps.npy')).astype(np.int64)
    acc50 = (not (set(x30_50) & (set(rows50) | set(q30_50))) and keep50[q30_50].all()
             and np.array_equal(L50, U44[keep50])
             and ST50['known_extreme_ray_orbits'] == 1003 + len(rows50) and ST50['realised_rows_r56'] == rows50
             and ST50['undecided'] == L50.shape[0] == 793 - len(rows50) - len(rows50x)
             and ST50['excluded'] == 1915910 + len(rows50x) and ST50['excluded_rows_r56'] == x30_50
             and ST50['qutrit_r56']['rows'] == q30_50
             and ST50['every_d']['open'] == len(ST50['every_d']['open_rows']) == 251 + len(rows50))
    F50 = np.load(os.path.join(D50, FN50), allow_pickle=False)
    fr50 = (np.array_equal(F50['rows'], np.flatnonzero(keep50)) and F50['cut'].shape == (F50['classes'].shape[0], L50.shape[0])
            and F50['cut'].any(0).all() and F50['cut'][F50['greedy_cover']].any(0).all())
    gate('G50 Gate 1: two common informations; fixed-structure and qutrit realisations',
         ex50 and ok50a and nv50 == nc50 == ne50 == len(rows50) and no50 and nq50 == len(Q50) and acc50 and fr50 and file50,
         f'{okx50}/{len(rows50x)} exclusion certificates re-checked exactly (stabcert.verify_exclusion_ci2) = the rows '
         f'ci2-ingleton reported infeasible ({len(cand50)}, failed {len(C50["failed"])}) {ex50}; {nv50}/{len(rows50)} '
         f'graph certificates and {nc50}/{len(rows50)} CSS matrices re-verified at lambda 2, {ne50} extreme rays of '
         f'QLR5, new S6 orbits {no50}; {nq50}/{len(Q50)} qutrit realisations (GF(3) column ranks, qutrit graph); '
         f'{ST50["known_extreme_ray_orbits"]} known / {ST50["undecided"]} undecided {acc50}; frontier '
         f'{F50["classes"].shape[0]} classes, greedy cover {len(F50["greedy_cover"])} cuts every row {fr50}; files '
         f'match the manifests {file50}')

    # G51 (R57): E2 at maxsize 2 completed, results/2026-11-24_r57-sandbox.  The job ci2-ingleton2
    # (results/2026-11-24_ci2-ingleton2) resumed ci2-ingleton: all 793 rows of the R54 list tested, the 433 rows of the
    # first run carried over unchanged, the seven rows realised in R56 feasible (as they must be).  The rows it reports
    # infeasible are the three certified in R56 (with the pairs of their certificates) and the rows certified here, each
    # re-checked by epr1kit.stabcert.verify_exclusion_ci2 on the ray of the R54 list; undecided_reps.npy is the R56 list
    # without them; status.json counts; on the remaining rows the R56 frontier still cuts every row and its greedy cover
    # still covers them all; the manifests of both new folders.
    D51 = os.path.join(base44, '2026-11-24_r57-sandbox')
    J51d = os.path.join(base44, '2026-11-24_ci2-ingleton2')
    file51 = True
    for d51, f51 in ((D51, ['README.txt', 'certificates_ci2.json', 'undecided_reps.npy', 'status.json']),
                     (J51d, ['ci2_ingleton2.json', 'ci2-ingleton2.log', 'JOB_STATUS.json'])):
        M51 = _js20.load(open(os.path.join(d51, 'manifest.json')))
        file51 = file51 and sorted(M51['files']) == sorted(f51) and all(
            _h33.sha256(open(os.path.join(d51, f), 'rb').read()).hexdigest() == e['sha256'] for f, e in M51['files'].items())
    J51 = _js20.load(open(os.path.join(J51d, 'ci2_ingleton2.json')))
    ctl51 = [int(k) for k in R50['row54']] + [int(q['row54']) for q in Q50]
    run51 = (J51['rays'] == J50['rays'] and J51['rows'] == 793 and len(J51['results']) == 793 and J51['summary']['complete']
             and all(J51['results'][k] == v for k, v in J50['results'].items())
             and all(J51['results'][str(k)]['infeasible'] is None for k in ctl51))
    C51 = _js20.load(open(os.path.join(D51, 'certificates_ci2.json')))
    cand51 = sorted(int(k) for k, v in J51['results'].items() if v['infeasible'])
    rows51x = sorted(int(c['index']) for c in C51['certificates'])
    okx51 = 0
    for c in C51['certificates']:
        r = L49[int(c['index'])]
        okc, _v = _sc.verify_exclusion_ci2(r, c['P1'], c['P2'], [(tuple(t[0]), _Fr50(t[1])) for t in c['terms']],
                                           F_S=c['F_S'])
        inf = J51['results'][str(c['index'])]['infeasible']
        okx51 += bool(okc and [int(x) for x in r] == c['ray'] and list(c['P1']) == list(inf[0])
                      and list(c['P2']) == list(inf[1])
                      and sum(a * int(b) for a, b in zip(c['F_S'], r)) == c['F_S_at_ray'] < 0)
    old51 = all(J51['results'][str(int(c['index']))]['infeasible'] == [list(c['P1']), list(c['P2'])]
                for c in C50['certificates'])
    ex51 = (okx51 == len(rows51x) > 0 and not C51['failed'] and old51 and not set(rows50x) & set(rows51x)
            and sorted(rows50x + rows51x) == cand51)
    x30_51 = [int(i49[k]) for k in rows51x]
    keep51 = keep50.copy()
    keep51[x30_51] = False
    L51 = np.load(os.path.join(D51, 'undecided_reps.npy')).astype(np.int64)
    ST51 = _js20.load(open(os.path.join(D51, 'status.json')))
    acc51 = (np.array_equal(L51, U44[keep51]) and not (set(x30_51) & (set(rows50) | set(q30_50)))
             and ST51['excluded'] == ST50['excluded'] + len(rows51x) and ST51['excluded_rows_r57'] == x30_51
             and ST51['known_extreme_ray_orbits'] == ST50['known_extreme_ray_orbits']
             and ST51['undecided'] == L51.shape[0] == ST50['undecided'] - len(rows51x))
    k51 = keep51[np.flatnonzero(keep50)]
    fr51 = bool(F50['cut'][:, k51].any(0).all() and F50['cut'][F50['greedy_cover']][:, k51].any(0).all())
    gate('G51 E2 at maxsize 2 completed: ci2-ingleton2',
         run51 and ex51 and acc51 and fr51 and file51,
         f'all 793 rows tested, the first run carried over, the {len(ctl51)} rows realised in R56 feasible {run51}; '
         f'infeasible rows = the {len(rows50x)} of R56 (same pairs {old51}) + {len(rows51x)} certified here, '
         f'{okx51}/{len(rows51x)} re-checked exactly (stabcert.verify_exclusion_ci2) {ex51}; '
         f'{ST51["known_extreme_ray_orbits"]} known / {ST51["undecided"]} undecided {acc51}; the R56 frontier and its '
         f'greedy cover still cut every remaining row {fr51}; files match the manifests {file51}')

    # G52 (R59): laptop pieces of 2026-11-24, results/2026-11-24_r59-sandbox.  The jobs css-scale
    # (results/2026-11-24_css-scale; eight rows, no realisation) and s44-a (results/2026-11-24_s44-a): the rows s44-a
    # realised are those of the realisations file, each re-checked as in G49 (s41 verify, stabcert on the stored graph,
    # extreme ray of QLR5 with tight rank 30, S6 orbit new against all known orbits) and still undecided in R57;
    # undecided_reps.npy is the R57 list without them; status.json counts; the frontier: its rows are the undecided
    # rows, every row cut, the greedy cover cuts every row; the manifests of the three new folders.
    D52 = os.path.join(base44, '2026-11-24_r59-sandbox')
    ST52 = _js20.load(open(os.path.join(D52, 'status.json')))
    RN52, FN52 = ST52['realisations_file'], ST52['frontier']['file']
    file52 = True
    for d52, f52 in ((D52, ['README.txt', RN52, 'undecided_reps.npy', 'status.json', FN52]),
                     (os.path.join(base44, '2026-11-24_css-scale'),
                      ['JOB_STATUS.json', 'css-scale.log', 'css_scale.npz', 'css_scale.npz.json', 'css_scale.plain.npz']),
                     (os.path.join(base44, '2026-11-24_s44-a'),
                      ['JOB_STATUS.json', 's44-a.log', 's44_a.npz', 's44_a.npz.json', 's44_a.plain.npz'])):
        M52 = _js20.load(open(os.path.join(d52, 'manifest.json')))
        file52 = file52 and sorted(M52['files']) == sorted(f52) and all(
            _h33.sha256(open(os.path.join(d52, f), 'rb').read()).hexdigest() == e['sha256'] for f, e in M52['files'].items())
    CS52 = _js20.load(open(os.path.join(base44, '2026-11-24_css-scale', 'css_scale.npz.json')))
    SA52 = _js20.load(open(os.path.join(base44, '2026-11-24_s44-a', 's44_a.npz.json')))
    R52 = np.load(os.path.join(D52, RN52), allow_pickle=False)
    rows52 = [int(x) for x in R52['row']]
    i50 = np.flatnonzero(keep50)
    ok52a = (len(rows52) == len(set(rows52)) > 0 and keep51[rows52].all()
             and [int(i50[int(x)]) for x in R52['row56']] == rows52
             and sorted(int(x) for x in R52['row56']) == sorted(SA52['found_rows'])
             and np.array_equal(R52['rep'].astype(np.int64), U44[rows52]) and (R52['lam'] == 2).all()
             and CS52['found'] == 0 and CS52['rows'] == 8)
    nv52 = sum(_sc.verify_realisation(U44[r], 2, R52['sizes'][k], R52['rows'][k]) for k, r in enumerate(rows52))
    nc52 = sum(_css47.verify(U44[r], 1, [int(c) for c in R52['cols'][k][:int(R52['m'][k])]])[0] for k, r in enumerate(rows52))
    V52 = HF44 @ U44[rows52].T
    ne52 = sum(1 for k in range(len(rows52)) if V52[:, k].min() >= 0 and core.rank_mod_p(HF44[V52[:, k] == 0]) == 30)
    c52, _ = core.class_reps(U44[rows52], core.perms31_s6(), offset=4096, width='u16')
    no52 = len(set(c52)) == len(rows52) and not (
        set(c52) & (set(c44) | set(c45) | set(c46) | set(c47) | set(c48) | set(c49) | set(c50)))
    keep52 = keep51.copy()
    keep52[rows52] = False
    L52 = np.load(os.path.join(D52, 'undecided_reps.npy')).astype(np.int64)
    acc52 = (np.array_equal(L52, U44[keep52])
             and ST52['known_extreme_ray_orbits'] == ST51['known_extreme_ray_orbits'] + len(rows52)
             and ST52['realised_rows_r59'] == rows52 and ST52['excluded'] == ST51['excluded']
             and ST52['undecided'] == L52.shape[0] == ST51['undecided'] - len(rows52))
    F52 = np.load(os.path.join(D52, FN52), allow_pickle=False)
    fr52 = bool(np.array_equal(F52['rows'], np.flatnonzero(keep52))
                and F52['cut'].shape == (F52['classes'].shape[0], L52.shape[0])
                and F52['cut'].any(0).all() and F52['cut'][F52['greedy_cover']].any(0).all())
    gate('G52 laptop pieces of 2026-11-24: css-scale, s44-a',
         ok52a and nv52 == nc52 == ne52 == len(rows52) and no52 and acc52 and fr52 and file52,
         f'css-scale {CS52["rows"]} rows, {CS52["found"]} realised; s44-a realised {len(rows52)}: {nv52} graph '
         f'certificates and {nc52} CSS matrices re-verified at lambda 2, {ne52} extreme rays of QLR5, new S6 orbits {no52}; '
         f'{ST52["known_extreme_ray_orbits"]} known / {ST52["undecided"]} undecided {acc52}; frontier '
         f'{F52["classes"].shape[0]} classes, greedy cover {len(F52["greedy_cover"])} cuts every row {fr52}; files match '
         f'the manifests {file52}')

    # G53 (R60): laptop pieces of the evening of 2026-11-24, results/2026-11-24_r60-sandbox.  The jobs s44-b and s44-c
    # (the fixed-structure search on the R56 list) and css-sweep3a, 3b and 3c (the plain CSS search on the R57 list),
    # results/2026-11-24_<job>: the rows they realised are those of the realisations file, each re-checked by s41's
    # verify (column-group ranks; cut ranks of the graph it builds, any size), by stabcert on the stored graph where
    # there is one (at most 64 qubits), as an extreme ray of QLR5 with tight rank 30 in an S6 orbit new against all
    # known orbits, still undecided in R59; undecided_reps.npy is the R59 list without them; status.json counts; the
    # frontier: its rows are the undecided rows, every row cut, the greedy cover cuts every row; the manifests of the
    # six new folders.
    D53 = os.path.join(base44, '2026-11-24_r60-sandbox')
    ST53 = _js20.load(open(os.path.join(D53, 'status.json')))
    RN53, FN53 = ST53['realisations_file'], ST53['frontier']['file']
    file53 = True
    J53 = {}
    for d53, f53 in [(D53, ['README.txt', RN53, 'undecided_reps.npy', 'status.json', FN53])] + [
            (os.path.join(base44, '2026-11-24_' + j), ['JOB_STATUS.json', j + '.log'] + [
                j.replace('-', '_') + e for e in ('.npz', '.npz.json', '.plain.npz')])
            for j in ('s44-b', 's44-c', 'css-sweep3a', 'css-sweep3b', 'css-sweep3c')]:
        M53 = _js20.load(open(os.path.join(d53, 'manifest.json')))
        file53 = file53 and sorted(M53['files']) == sorted(f53) and all(
            _h33.sha256(open(os.path.join(d53, f), 'rb').read()).hexdigest() == e['sha256'] for f, e in M53['files'].items())
        if d53 != D53:
            J53[f53[1][:-4]] = _js20.load(open(os.path.join(d53, f53[3])))
    R53 = np.load(os.path.join(D53, RN53), allow_pickle=False)
    rows53 = [int(x) for x in R53['row']]
    job53 = [str(x) for x in R53['job']]
    i51 = np.flatnonzero(keep51)
    gr53 = [bool(x) for x in R53['graph']]
    ok53a = (len(rows53) == len(set(rows53)) > 0 and keep52[rows53].all()
             and [int(i50[int(x)]) for x in R53['row56']] == rows53 and [int(i51[int(x)]) for x in R53['row57']] == rows53
             and sorted(int(x) for k, x in enumerate(R53['row56']) if job53[k] == 's44-b') == sorted(J53['s44-b']['found_rows'])
             and sorted(int(x) for k, x in enumerate(R53['row57']) if job53[k] == 'css-sweep3a')
             == sorted(J53['css-sweep3a']['found_rows'])
             and set(job53) <= {'s44-b', 'css-sweep3a'} and sum(J['found'] for J in J53.values()) == len(rows53)
             and np.array_equal(R53['rep'].astype(np.int64), U44[rows53]) and (R53['lam'] == 2).all()
             and gr53 == [int(R53['sizes'][k].sum()) <= 64 for k in range(len(rows53))])
    nv53 = sum(_sc.verify_realisation(U44[r], 2, R53['sizes'][k], R53['rows'][k]) for k, r in enumerate(rows53)
               if gr53[k] and int(R53['sizes'][k].sum()) <= 64)
    nc53 = sum(_css47.verify(U44[r], 1, [int(c) for c in R53['cols'][k][:int(R53['m'][k])]])[0] for k, r in enumerate(rows53))
    V53 = HF44 @ U44[rows53].T
    ne53 = sum(1 for k in range(len(rows53)) if V53[:, k].min() >= 0 and core.rank_mod_p(HF44[V53[:, k] == 0]) == 30)
    c53, _ = core.class_reps(U44[rows53], core.perms31_s6(), offset=4096, width='u16')
    no53 = len(set(c53)) == len(rows53) and not (
        set(c53) & (set(c44) | set(c45) | set(c46) | set(c47) | set(c48) | set(c49) | set(c50) | set(c52)))
    keep53 = keep52.copy()
    keep53[rows53] = False
    L53 = np.load(os.path.join(D53, 'undecided_reps.npy')).astype(np.int64)
    acc53 = (np.array_equal(L53, U44[keep53])
             and ST53['known_extreme_ray_orbits'] == ST52['known_extreme_ray_orbits'] + len(rows53)
             and ST53['realised_rows_r60'] == rows53 and ST53['excluded'] == ST52['excluded']
             and ST53['undecided'] == L53.shape[0] == ST52['undecided'] - len(rows53))
    F53 = np.load(os.path.join(D53, FN53), allow_pickle=False)
    fr53 = bool(np.array_equal(F53['rows'], np.flatnonzero(keep53))
                and F53['cut'].shape == (F53['classes'].shape[0], L53.shape[0])
                and F53['cut'].any(0).all() and F53['cut'][F53['greedy_cover']].any(0).all())
    gate('G53 laptop pieces of 2026-11-24 (evening): s44-b, s44-c, css-sweep3a-c',
         ok53a and nv53 == sum(gr53) and nc53 == ne53 == len(rows53) and no53 and acc53 and fr53 and file53,
         f'{len(rows53)} realised ({job53.count("s44-b")} by s44-b, {job53.count("css-sweep3a")} by css-sweep3a, none by '
         f's44-c, css-sweep3b, css-sweep3c): {nc53} CSS matrices re-verified at lambda 2 (s41: column-group ranks, cut '
         f'ranks of the graph, any size), {nv53}/{sum(gr53)} stored graphs (at most 64 qubits) by stabcert, {ne53} '
         f'extreme rays of QLR5, new S6 orbits {no53}; {ST53["known_extreme_ray_orbits"]} known / {ST53["undecided"]} '
         f'undecided {acc53}; frontier {F53["classes"].shape[0]} classes, greedy cover {len(F53["greedy_cover"])} cuts '
         f'every row {fr53}; files match the manifests {file53}')

    # G54 (R63): laptop jobs of 2026-11-25/26, results/2026-11-26_r63-sandbox.  (a) ci3-8h (results/2026-11-25_ci3), the
    # first piece of E2 at maxsize 3 on the R60 list: every row of ci3.json feasible, done = combos = s42's own count at
    # maxsize 3 without nesting, no solver trouble, the summary equal to the rows, status.json's ci3 block.  (b)
    # css-sweep3d and css-sweep3e: nothing found, 32 and 14 rows tried.  (c) s32-new5a ... s32-new5d: each certificate
    # file exact (s32_unimodular check on the catalogue extended by the six realisation files of R49-R54) and holding
    # those of s32-new4; their union, as R30 rows, is status.json's every_d rows, and the open rows are those of R60
    # without them.  (d) Counts unchanged from R60.  (e) The manifests of the eight folders.
    D54 = os.path.join(base44, '2026-11-26_r63-sandbox')
    ST54 = _js20.load(open(os.path.join(D54, 'status.json')))
    file54 = True
    for d54, f54 in [(D54, ['README.txt', 'status.json']),
                     (os.path.join(base44, '2026-11-25_ci3'), ['JOB_STATUS.json', 'ci3-8h.log', 'ci3.json'])] + [
            (os.path.join(base44, '2026-11-26_css-sweep3' + x), ['JOB_STATUS.json', f'css-sweep3{x}.log', f'css_sweep3{x}.npz',
                                                                f'css_sweep3{x}.npz.json', f'css_sweep3{x}.plain.npz'])
            for x in 'de'] + [
            (os.path.join(base44, '2026-11-26_s32-new5' + x), ['JOB_STATUS.json', f's32-new5{x}.log', f's32n5{x}_raw.txt',
                                                              f's32n5{x}_unimodular.npz'])
            for x in 'abcd']:
        M54 = _js20.load(open(os.path.join(d54, 'manifest.json')))
        file54 = file54 and sorted(M54['files']) == sorted(f54) and all(
            _h33.sha256(open(os.path.join(d54, f), 'rb').read()).hexdigest() == e['sha256'] for f, e in M54['files'].items())
    import s42_ci2_ingleton as _ci54
    L54n = 'results/2026-11-24_r60-sandbox/undecided_reps.npy'
    J54 = _js20.load(open(os.path.join(base44, '2026-11-25_ci3', 'ci3.json')))
    R54c = J54['results']
    ok54a = ((J54['rays'], J54['rows'], J54['maxsize'], J54['nested']) == (L54n, 777, 3, False)
             and all(e['infeasible'] is None and e['trouble'] == 0
                     and e['done'] == e['combos'] == len(_ci54.combos(_ci22.h_norm(L53[int(k)]), 3, False))
                     for k, e in R54c.items())
             and J54['summary']['tested'] == len(R54c) == ST54['ci3']['rows_done'] == 777 - ST54['ci3']['rows_left']
             and J54['summary']['infeasible_rows'] == [] and J54['summary']['trouble_lps'] == 0
             and J54['summary']['lps'] == sum(e['done'] for e in R54c.values()) == ST54['ci3']['pieces'][0]['lps'])
    ok54b = True
    for x54, n54 in (('d', 32), ('e', 14)):
        Jx54 = _js20.load(open(os.path.join(base44, '2026-11-26_css-sweep3' + x54, f'css_sweep3{x54}.npz.json')))
        Zx54 = np.load(os.path.join(base44, '2026-11-26_css-sweep3' + x54, f'css_sweep3{x54}.npz'), allow_pickle=False)
        ok54b = ok54b and Jx54['found'] == 0 and len(Zx54['index']) == 0 and len(Jx54['best_residual']) == n54
    _s36.APPEND = [os.path.join(D44, 'realisations_141.npz'), os.path.join(D45, 'realisations_126.npz'),
                   os.path.join(D46, 'realisations_73.npz'), os.path.join(D47, 'realisations_5.npz'),
                   os.path.join(D48, 'realisations_9.npz'), os.path.join(D49, 'realisations_8.npz')]
    u54, ok54c = set(), True
    try:
        for x54 in 'abcd':
            P54 = os.path.join(base44, '2026-11-26_s32-new5' + x54, f's32n5{x54}_unimodular.npz')
            S54 = _s36.check(P54, full=False, log=lambda *_x: None)
            r54 = set(int(k) for k in np.load(P54, allow_pickle=False)['row'])
            ok54c = (ok54c and S54['ok'] == S54['stored'] == S54['distinct_rows'] == len(r54) and S54['rep_ok']
                     and set(int(k) for k in EV49['row']) <= r54 and all(641 <= k < 1003 for k in r54))
            u54 |= r54
    finally:
        _s36.APPEND = None
    cat54 = cat49 + rows49
    rows54e = sorted(cat54[k - 641] for k in u54)
    ST53e = ST53['every_d']
    open54 = sorted(set(ST53e['open_rows']) - set(rows54e))
    ok54c = (ok54c and rows54e == ST54['every_d']['rows'] and ST54['every_d']['unimodular'] == len(u54)
             and open54 == ST54['every_d']['open_rows'] and ST54['every_d']['open'] == len(open54)
             and len(u54) + len(open54) == ST54['realised_after_submission'])
    acc54 = all(ST54[k] == ST53[k] for k in ('excluded', 'known_extreme_ray_orbits', 'realised_after_submission',
                                             'undecided'))
    gate('G54 laptop jobs of 2026-11-25/26: ci3-8h, css-sweep3d-e, s32-new5a-d',
         ok54a and ok54b and ok54c and acc54 and file54,
         f'ci3 first piece: {len(R54c)} of 777 rows feasible, done = combos = s42 count at maxsize 3, '
         f'{J54["summary"]["lps"]} LPs, no trouble {ok54a}; css-sweep3d/e nothing found {ok54b}; every d: '
         f'{len(u54)} certificates exact over the four s32-new5 files (s32-new4 included), {len(open54)} open {ok54c}; '
         f'counts unchanged ({ST54["known_extreme_ray_orbits"]} known / {ST54["undecided"]} undecided) {acc54}; files match '
         f'the manifests {file54}')

    # G55 (R64): laptop pieces of the afternoon of 2026-11-26, results/2026-11-26_r64-sandbox.  The jobs s44-d and s44-e
    # (the fixed-structure search again, with new seeds, on the R60 list), results/2026-11-26_<job>: the row s44-d
    # realised is that of the realisations file, re-checked by s41's verify (column-group ranks; cut ranks of the graph
    # it builds), by stabcert on the stored graph (64 qubits), as an extreme ray of QLR5 with tight rank 30 in an S6
    # orbit new against all known orbits, still undecided in R63; s44-e found nothing; undecided_reps.npy is the R60
    # list without it; status.json: counts, and every_d with the new row open; the frontier: its rows are the undecided
    # rows, every row cut, the greedy cover cuts every row; the manifests of the three new folders.
    D55 = os.path.join(base44, '2026-11-26_r64-sandbox')
    ST55 = _js20.load(open(os.path.join(D55, 'status.json')))
    RN55, FN55 = ST55['realisations_file'], ST55['frontier']['file']
    file55 = True
    J55 = {}
    for d55, f55 in [(D55, ['README.txt', RN55, 'undecided_reps.npy', 'status.json', FN55])] + [
            (os.path.join(base44, '2026-11-26_' + j), ['JOB_STATUS.json', j + '.log'] + [
                j.replace('-', '_') + e for e in ('.npz', '.npz.json', '.plain.npz')])
            for j in ('s44-d', 's44-e')]:
        M55 = _js20.load(open(os.path.join(d55, 'manifest.json')))
        file55 = file55 and sorted(M55['files']) == sorted(f55) and all(
            _h33.sha256(open(os.path.join(d55, f), 'rb').read()).hexdigest() == e['sha256'] for f, e in M55['files'].items())
        if d55 != D55:
            J55[f55[1][:-4]] = _js20.load(open(os.path.join(d55, f55[3])))
    R55 = np.load(os.path.join(D55, RN55), allow_pickle=False)
    rows55 = [int(x) for x in R55['row']]
    job55 = [str(x) for x in R55['job']]
    i53 = np.flatnonzero(keep53)
    gr55 = [bool(x) for x in R55['graph']]
    Ze55 = np.load(os.path.join(base44, '2026-11-26_s44-e', 's44_e.npz'), allow_pickle=False)
    ok55a = (len(rows55) == len(set(rows55)) > 0 and keep53[rows55].all()
             and [int(i53[int(x)]) for x in R55['row60']] == rows55
             and sorted(int(x) for x in R55['row60']) == sorted(J55['s44-d']['found_rows'])
             and set(job55) == {'s44-d'} and J55['s44-d']['found'] == len(rows55)
             and J55['s44-e']['found'] == 0 and J55['s44-e']['found_rows'] == [] and len(Ze55['index']) == 0
             and np.array_equal(R55['rep'].astype(np.int64), U44[rows55]) and (R55['lam'] == 2).all()
             and gr55 == [int(R55['sizes'][k].sum()) <= 64 for k in range(len(rows55))])
    nv55 = sum(_sc.verify_realisation(U44[r], 2, R55['sizes'][k], R55['rows'][k]) for k, r in enumerate(rows55)
               if gr55[k] and int(R55['sizes'][k].sum()) <= 64)
    nc55 = sum(_css47.verify(U44[r], 1, [int(c) for c in R55['cols'][k][:int(R55['m'][k])]])[0] for k, r in enumerate(rows55))
    V55 = HF44 @ U44[rows55].T
    ne55 = sum(1 for k in range(len(rows55)) if V55[:, k].min() >= 0 and core.rank_mod_p(HF44[V55[:, k] == 0]) == 30)
    c55, _ = core.class_reps(U44[rows55], core.perms31_s6(), offset=4096, width='u16')
    no55 = len(set(c55)) == len(rows55) and not (
        set(c55) & (set(c44) | set(c45) | set(c46) | set(c47) | set(c48) | set(c49) | set(c50) | set(c52) | set(c53)))
    keep55 = keep53.copy()
    keep55[rows55] = False
    L55 = np.load(os.path.join(D55, 'undecided_reps.npy')).astype(np.int64)
    EV55, EV54 = ST55['every_d'], ST54['every_d']
    acc55 = (np.array_equal(L55, U44[keep55])
             and ST55['known_extreme_ray_orbits'] == ST54['known_extreme_ray_orbits'] + len(rows55)
             and ST55['realised_after_submission'] == ST54['realised_after_submission'] + len(rows55)
             and ST55['realised_rows_r64'] == rows55 and ST55['excluded'] == ST54['excluded']
             and ST55['undecided'] == L55.shape[0] == ST54['undecided'] - len(rows55)
             and EV55['rows'] == EV54['rows'] and EV55['unimodular'] == EV54['unimodular']
             and EV55['open_rows'] == sorted(EV54['open_rows'] + rows55) and EV55['open'] == len(EV55['open_rows'])
             and EV55['unimodular'] + EV55['open'] == ST55['realised_after_submission'])
    F55 = np.load(os.path.join(D55, FN55), allow_pickle=False)
    fr55 = bool(np.array_equal(F55['rows'], np.flatnonzero(keep55))
                and F55['cut'].shape == (F55['classes'].shape[0], L55.shape[0])
                and F55['cut'].any(0).all() and F55['cut'][F55['greedy_cover']].any(0).all())
    gate('G55 laptop pieces of 2026-11-26 (afternoon): s44-d, s44-e',
         ok55a and nv55 == sum(gr55) and nc55 == ne55 == len(rows55) and no55 and acc55 and fr55 and file55,
         f'{len(rows55)} realised by s44-d (R30 rows {rows55}), none by s44-e: {nc55} CSS matrices re-verified at lambda 2 '
         f'(s41: column-group ranks, cut ranks of the graph), {nv55}/{sum(gr55)} stored graphs by stabcert, {ne55} '
         f'extreme rays of QLR5, new S6 orbits {no55}; {ST55["known_extreme_ray_orbits"]} known / {ST55["undecided"]} '
         f'undecided, every d {EV55["unimodular"]} certified / {EV55["open"]} open {acc55}; frontier '
         f'{F55["classes"].shape[0]} classes, greedy cover {len(F55["greedy_cover"])} cuts every row {fr55}; files match '
         f'the manifests {file55}')

    # G56 (R66): laptop pieces of the afternoon of 2026-11-26, results/2026-11-26_r66-sandbox.  The jobs css-sweep4a ...
    # css-sweep4e (the CSS search again, with a new seed, on the 136 rows with N <= 30 of the R60 list, split
    # 27/27/28/27/27), results/2026-11-26_<job>: the rows 4a, 4b and 4c realised are those of the realisations file,
    # re-checked by s41's verify (column-group ranks; cut ranks of the graph it builds), by stabcert on the stored graphs,
    # as extreme rays of QLR5 with tight rank 30 in S6 orbits new against all known orbits, still undecided in R64; 4d and
    # 4e found nothing; undecided_reps.npy is the R64 list without the new rows; status.json: counts, and every_d with
    # the new rows open; the frontier: its rows are the undecided rows, every row cut, the greedy cover cuts every row;
    # the manifests of the six new folders.
    D56 = os.path.join(base44, '2026-11-26_r66-sandbox')
    ST56 = _js20.load(open(os.path.join(D56, 'status.json')))
    RN56, FN56 = ST56['realisations_file'], ST56['frontier']['file']
    JB56 = [f'css-sweep4{x}' for x in 'abcde']
    file56 = True
    J56 = {}
    for d56, f56 in [(D56, ['README.txt', RN56, 'undecided_reps.npy', 'status.json', FN56])] + [
            (os.path.join(base44, '2026-11-26_' + j), ['JOB_STATUS.json', j + '.log'] + [
                j.replace('-', '_') + e for e in ('.npz', '.npz.json', '.plain.npz')]) for j in JB56]:
        M56 = _js20.load(open(os.path.join(d56, 'manifest.json')))
        file56 = file56 and sorted(M56['files']) == sorted(f56) and all(
            _h33.sha256(open(os.path.join(d56, f), 'rb').read()).hexdigest() == e['sha256'] for f, e in M56['files'].items())
        if d56 != D56:
            J56[f56[1][:-4]] = _js20.load(open(os.path.join(d56, f56[3])))
    Q56 = np.load(os.path.join(D56, RN56), allow_pickle=False)
    rows56 = [int(x) for x in Q56['row']]
    job56 = [str(x) for x in Q56['job']]
    gr56 = [bool(x) for x in Q56['graph']]
    i55 = np.flatnonzero(keep55)
    le56 = [k for k in range(len(L53)) if sum(_css47.singles(L53[k])) <= 30]
    cut56 = [0, 27, 54, 82, 109, 136]
    tr56 = {j: sorted(set(int(k.split(':')[0]) for k in J56[j]['best_residual']) | set(J56[j]['found_rows'])) for j in JB56}
    emp56 = all(len(np.load(os.path.join(base44, '2026-11-26_' + j, j.replace('-', '_') + e))['index']) == 0
                for j in JB56 if not J56[j]['found_rows'] for e in ('.npz', '.plain.npz'))
    ok56a = (len(rows56) == len(set(rows56)) > 0 and keep55[rows56].all() and emp56
             and [int(i53[int(x)]) for x in Q56['row60']] == rows56 and [int(i55[int(x)]) for x in Q56['row64']] == rows56
             and all(sorted(int(x) for k, x in enumerate(Q56['row60']) if job56[k] == j) == sorted(J56[j]['found_rows'])
                     for j in JB56) and set(job56) <= set(JB56)
             and len(le56) == 136 and all(tr56[j] == le56[cut56[t]:cut56[t + 1]] for t, j in enumerate(JB56))
             and np.array_equal(Q56['rep'].astype(np.int64), U44[rows56]) and (Q56['lam'] == 2).all()
             and gr56 == [int(Q56['sizes'][k].sum()) <= 64 for k in range(len(rows56))])
    nv56 = sum(_sc.verify_realisation(U44[r], 2, Q56['sizes'][k], Q56['rows'][k]) for k, r in enumerate(rows56) if gr56[k])
    nc56 = sum(_css47.verify(U44[r], 1, [int(c) for c in Q56['cols'][k][:int(Q56['m'][k])]])[0] for k, r in enumerate(rows56))
    V56 = HF44 @ U44[rows56].T
    ne56 = sum(1 for k in range(len(rows56)) if V56[:, k].min() >= 0 and core.rank_mod_p(HF44[V56[:, k] == 0]) == 30)
    c56, _ = core.class_reps(U44[rows56], core.perms31_s6(), offset=4096, width='u16')
    no56 = len(set(c56)) == len(rows56) and not (
        set(c56) & (set(c44) | set(c45) | set(c46) | set(c47) | set(c48) | set(c49) | set(c50) | set(c52) | set(c53)
                    | set(c55)))
    keep56 = keep55.copy()
    keep56[rows56] = False
    L56 = np.load(os.path.join(D56, 'undecided_reps.npy')).astype(np.int64)
    EV56 = ST56['every_d']
    acc56 = (np.array_equal(L56, U44[keep56])
             and ST56['known_extreme_ray_orbits'] == ST55['known_extreme_ray_orbits'] + len(rows56)
             and ST56['realised_after_submission'] == ST55['realised_after_submission'] + len(rows56)
             and ST56['realised_rows_r66'] == rows56 and ST56['excluded'] == ST55['excluded']
             and ST56['undecided'] == L56.shape[0] == ST55['undecided'] - len(rows56)
             and EV56['rows'] == EV55['rows'] and EV56['unimodular'] == EV55['unimodular']
             and EV56['open_rows'] == sorted(EV55['open_rows'] + rows56) and EV56['open'] == len(EV56['open_rows'])
             and EV56['unimodular'] + EV56['open'] == ST56['realised_after_submission'])
    F56 = np.load(os.path.join(D56, FN56), allow_pickle=False)
    fr56 = bool(np.array_equal(F56['rows'], np.flatnonzero(keep56))
                and F56['cut'].shape == (F56['classes'].shape[0], L56.shape[0])
                and F56['cut'].any(0).all() and F56['cut'][F56['greedy_cover']].any(0).all())
    gate('G56 laptop pieces of 2026-11-26 (afternoon): css-sweep4a-e',
         ok56a and nv56 == sum(gr56) and nc56 == ne56 == len(rows56) and no56 and acc56 and fr56 and file56,
         f'{len(rows56)} realised by css-sweep4a-c (R30 rows {rows56}), none by 4d and 4e, the pieces cover the 136 rows '
         f'with N <= 30: {nc56} CSS matrices re-verified at lambda 2 (s41: column-group ranks, cut ranks of the graph), '
         f'{nv56}/{sum(gr56)} stored graphs by stabcert, {ne56} extreme rays of QLR5, new S6 orbits {no56}; '
         f'{ST56["known_extreme_ray_orbits"]} known / {ST56["undecided"]} undecided, every d {EV56["unimodular"]} '
         f'certified / {EV56["open"]} open {acc56}; frontier {F56["classes"].shape[0]} classes, greedy cover '
         f'{len(F56["greedy_cover"])} cuts every row {fr56}; files match the manifests {file56}')

    # G57 (R68): the second piece of E2 at maxsize 3 (job ci3-8h, results/2026-11-27_ci3) on the R60 list,
    # results/2026-11-27_r68-sandbox: the rows of the first piece carried over unchanged; every row of ci3.json feasible,
    # done = combos = s42's own count at maxsize 3 without nesting, no solver trouble, the summary equal to the rows;
    # status.json's ci3 block (both pieces, rows and LPs left); the rows of the R60 list realised in R64 and R66 that the
    # run reached are feasible, the others still to do; counts, every_d, frontier and undecided list those of R66; the
    # manifests of the two new folders.
    D57 = os.path.join(base44, '2026-11-27_r68-sandbox')
    ST57 = _js20.load(open(os.path.join(D57, 'status.json')))
    file57 = True
    for d57, f57 in [(D57, ['README.txt', 'status.json']),
                     (os.path.join(base44, '2026-11-27_ci3'), ['JOB_STATUS.json', 'ci3-8h.log', 'ci3.json'])]:
        M57 = _js20.load(open(os.path.join(d57, 'manifest.json')))
        file57 = file57 and sorted(M57['files']) == sorted(f57) and all(
            _h33.sha256(open(os.path.join(d57, f), 'rb').read()).hexdigest() == e['sha256'] for f, e in M57['files'].items())
    J57 = _js20.load(open(os.path.join(base44, '2026-11-27_ci3', 'ci3.json')))
    R57c, C57 = J57['results'], ST57['ci3']
    left57 = [k for k in range(len(L53)) if str(k) not in R57c]
    ok57a = ((J57['rays'], J57['rows'], J57['maxsize'], J57['nested']) == (L54n, 777, 3, False)
             and all(k in R57c and R57c[k] == v for k, v in R54c.items())
             and all(e['infeasible'] is None and e['trouble'] == 0
                     and e['done'] == e['combos'] == len(_ci54.combos(_ci22.h_norm(L53[int(k)]), 3, False))
                     for k, e in R57c.items())
             and J57['summary']['tested'] == len(R57c) == C57['rows_done'] == 777 - C57['rows_left'] == 777 - len(left57)
             and J57['summary']['infeasible_rows'] == [] and J57['summary']['trouble_lps'] == 0
             and J57['summary']['lps'] == sum(e['done'] for e in R57c.values()) == C57['lps_done']
             == sum(p['lps'] for p in C57['pieces'])
             and C57['pieces'][0] == ST56['ci3']['pieces'][0] and C57['pieces'][1]['rows_done'] == len(R57c) - len(R54c)
             and C57['lps_left'] == sum(len(_ci54.combos(_ci22.h_norm(L53[k]), 3, False)) for k in left57))
    real57 = [int(x) for x in R55['row60']] + [int(x) for x in Q56['row60']]
    ok57b = ({int(k): v['ci3'] for k, v in C57['realised_rows60'].items()}
             == {k: 'feasible' if str(k) in R57c else 'still to do' for k in real57}
             and all(R57c[str(k)]['infeasible'] is None for k in real57 if str(k) in R57c))
    acc57 = (ST57['known_extreme_ray_orbits'] == ST56['known_extreme_ray_orbits'] and ST57['undecided'] == ST56['undecided']
             and ST57['excluded'] == ST56['excluded']
             and ST57['realised_after_submission'] == ST56['realised_after_submission']
             and ST57['every_d'] == ST56['every_d'] and ST57['frontier'] == 'unchanged: results/2026-11-26_r66-sandbox/'
             + ST56['frontier']['file']
             and ST57['undecided_reps'].startswith('unchanged: results/2026-11-26_r66-sandbox/undecided_reps.npy'))
    gate('G57 laptop job of 2026-11-26/27: ci3-8h (second piece)',
         ok57a and ok57b and acc57 and file57,
         f'ci3: {len(R57c)} of 777 rows feasible ({len(R57c) - len(R54c)} new in this piece, the {len(R54c)} of the first '
         f'carried over), done = combos = s42 count at maxsize 3, no solver trouble, {C57["rows_left"]} rows and '
         f'{C57["lps_left"]} LPs left {ok57a}; realised rows of the R60 list reached by the run feasible '
         f'({sum(1 for k in real57 if str(k) in R57c)} of {len(real57)}) {ok57b}; counts unchanged '
         f'({ST57["known_extreme_ray_orbits"]} known / {ST57["undecided"]} undecided) {acc57}; files match the manifests '
         f'{file57}')

    if a.full:
        R29l = _r29.lrs_rays(A29)
        if R29l is None:
            gate('G29b n4-lrs', True, 'lrs not installed: exact enumeration skipped (the stored 46 rays are checked by G29)')
        else:
            gate('G29b n4-lrs', R29l.shape[0] == 46 and {tuple(x) for x in R29l} == {tuple(x) for x in R29},
                 f'lrs: {R29l.shape[0]} extreme rays of h^-1(Shannon + Ingleton) on five elements (expect 46 = stored set)')

        LA31b, LB31b, n31b = _g31.lemma_forms(2)
        mA, cA, fullA = _g31.exhaustive01_min(LA31b, n31b, 3, 2)
        mB, cB, fullB = _g31.exhaustive01_min(LB31b, n31b, 3, 3)
        sA, _c, _f = _g31.exhaustive01_min(LB31b, n31b, 3, 2)
        sB, _c, _f = _g31.exhaustive01_min(LA31b, n31b, 3, 3)
        gate('G31b Lemma 4 t=2', fullA and fullB and mA >= 0 and mB >= 0 and sA < 0 and sB < 0,
             f'all {cA} {{0,1}}-point configurations: LA over GF(2) min {mA}, LB over GF(3) min {mB} (>= 0); '
             f'sensitivity: LB over GF(2) min {sA}, LA over GF(3) min {sB} (< 0)')

        W32b = _z32.analyse(6, 4)
        Q32b = [_z32.analyse(6, q)['FA'] for q in (2, 3)]          # (LA)_6 is valid for exponent dividing 6
        okW32b = (W32b['free'] and W32b['weights'] and W32b['self_dual'] and W32b['h2r'] and W32b['contraction']
                  and W32b['FA'] == -2 and np.array_equal(Z32['t6_N4'], W32b['S'])
                  and core.sha_rows_wide(W32b['S'][None, :]) == MAN['s29_zd_witnesses']['sha256_int16_t6_N4'])
        L32b = _z32.lemma_tests(2_200_000, np.random.default_rng(29), log=lambda *_: None)
        okL32b = all((mn >= 0) == valid for _lab, mn, _c, _f, valid in L32b)
        gate('G32b Z_4 vs Z_6 + Lemma 4(G)', okW32b and all(v >= 0 for v in Q32b) and okL32b,
             f'Z_4, t=6: (LA)_6 = {W32b["FA"]} (expect -2), on GF(2)/GF(3), t=6: {Q32b} (>= 0); Lemma 4(G) minima '
             + ', '.join(f'{lab} {mn}' for lab, mn, *_r in L32b) + ' (claimed-valid forms >= 0, the others < 0)')

        oc33b = {t: _d33.one_cut_check(t) for t in (6, 7)}
        ok33b = all(ok for res in oc33b.values() for (_c, _k, ok) in res.values())
        gate('G33b one cut t=6,7', ok33b,
             'rank changes over GF(l) vs Q as predicted: ' + '; '.join(
                 f't={t}: ' + ', '.join(f'GF({l}) {c} sets' for l, (c, _k, _ok) in res.items()) for t, res in oc33b.items()))

        E34b = _d33.bkl_check(exhaustive=True)['exhaustive']
        gate('G34b BKL exhaustive', E34b == {'bkl68_even over GF(2)': 0, 'bkl617_odd_corrected over GF(3)': 0,
                                             'bkl617_odd (printed) over GF(3)': -3},
             f'minima over all {{0,1}}-point configurations of GF(p)^3: {E34b} (expect 0, 0, -3)')

        ex35b = _s35.exhaustive(f35)
        rng35b = np.random.default_rng(351)
        lem35b = {q: _s35.lemma_checks(q, 200, rng35b, f35) for q in (2, 3, 5)}
        ok35b = (ex35b[2] == -3 and ex35b[3] >= 0 and ex35b[5] >= 0
                 and all(v >= 0 for q in (2, 3, 5) for k, v in lem35b[q].items() if k.startswith('Lemma'))
                 and all(v >= 0 for q in (3, 5) for v in lem35b[q].values()))
        gate('G35b six parties exhaustive', ok35b,
             f'f9 minima over all {{0,1}}-point configurations of GF(p)^3: {ex35b} (expect -3, >= 0, >= 0); lemma checks on '
             f'200 families per field: all slacks >= 0 {ok35b}; chain minima '
             + ', '.join(f'GF({q}) {lem35b[q]["f9(h)"]}' for q in (2, 3, 5)))

        S36b = _s36.check(full=True)
        st36b = _s36.self_test(log=lambda *_x: None)
        gate('G36b five-party unimodular, second method', S36b['ok'] == S36b['stored'] and S36b['zd_ok'] == S36b['zd_checked']
             and st36b,
             f'rational rank + minors with gcd 1, and GF(p) ranks for p <= 19: {S36b["ok"]}/{S36b["stored"]}; Z_4 / Z_6 state '
             f'vectors {S36b["zd_ok"]}/{S36b["zd_checked"]} certificates; s32 self-test {st36b}')

        S38b = _s38.check(full_check=True)
        st38b = _s38.self_test(trials=10, log=lambda *_x: None)
        F38b = _s38.fields_check((2, 3, 4), recheck='first', log=lambda *_x: None)
        fc38 = M38['field_counts']
        f38ok = all(F38b[q] == fc38[str(q)] for q in F38b)
        lp38 = _s38.line_pieces_check()['ok']
        gate('G38b six parties, enumeration re-checked', _s38.check_ok(S38b, True) and st38b and f38ok and lp38,
             f'GF(2) list recomputed in numpy: {S38b["gf2_reenumerated"]} vectors, equal to the stored list '
             f'{S38b["gf2_list_equal"]}; {S38b["gf3_samples"]} random GF(3) graphs, {S38b["gf3_samples_missing"]} '
             f'outside the GF(3) list; 2 S a sum of two GF(2) vectors for {S38b["orbits_2v_a_sum"]} orbits; L1 distance '
             f'S_NF to the GF(2) list {S38b["nf_l1_to_gf2"]} (only the Steane vector: {S38b["nf_nearest_is_steane"]}), '
             f'2 S_NF to sums of two {S38b["two_nf_l1_to_sums"]}; Python target search, (S_F, S_NF) solutions over '
             f'GF(2), GF(3), GF(4): {[F38b[q] for q in (2, 3, 4)]}, equal to the counts of cc/s36_target7.c {f38ok}; '
             f'line pieces of the Fano and non-Fano states (Proposition 10, steps 1 and 1\') {lp38}; s34 self-test {st38b}')

        T39b = _s38.transfer_check(full_check=True, trials=40, log=lambda *_x: None)
        gate('G39b Theorem 12, the proof step by step', _s38.transfer_ok(T39b, True),
             f'{T39b["chain_states"]} states {T39b["chain_kinds"]}: steps failing on {T39b["chain_failed"]} (expect 0); '
             f'minimum f12 {T39b["chain_min_f12"]}, {T39b["chain_f12_zero"]} states at 0; dual arrangement built explicitly '
             f'{T39b["dual_identity"]}; compression keeps the entropies {T39b["compression"]}; seven-point states from '
             f'state vectors: largest deviation {T39b["seven_point_sv_dev"]:.1e}, f12 {T39b["seven_point_sv_f12"]} bits')

        T40b = _s40.check(full_check=True, trials=60, log=lambda *_x: None)
        ch40 = T40b['chain']
        lf40 = T40b['lifts']
        gate('G40b Theorem 14, the proof step by step', _s40.check_ok(T40b, True)
             and all(lf40[n]['zeros'] == M40['lifts_zeros'][n] for n in lf40),
             '; '.join(f'Z_{1 << k}: {v["failed"]} of {sum(v["kinds"].values())} states fail a step (expect 0), minimum '
                       f'a f^(a) {v["min_kfk"]}, {v["at_zero"]} at 0, minimum f12 {v["min_f12"]}' for k, v in ch40.items())
             + f'; state vectors: largest deviation {T40b["state_vector_dev"]:.1e}; Lemma 4(G\') minima {T40b["lemma_tests"]} '
             f'(valid >= 0, smaller exponent < 0); exhaustive: '
             + '; '.join(f'{n}: {v["codes"]} codes, minimum {v["min"]} (0), {v["zeros"]} at 0, weakened forms '
                         f'{v["weakened"]}' for n, v in lf40.items()))

        T41b = _s41.check(full_check=True, log=lambda *_x: None)
        en41 = T41b['enumeration']
        ex41 = {k: [v['graphs'], v['distinct']] for k, v in en41.items()}
        gate('G41b one qudit per system: exhaustive GF(2)-GF(5), Propositions 15 and 16', _s41.check_ok(T41b, True)
             and ex41 == M41['exhaustive'],
             '; '.join(f'{k}: {v["graphs"]} graphs, {v["distinct"]} vectors, {v["orbits"]} orbits = table {v["orbits_equal"]}'
                       for k, v in en41.items())
             + f'; target searches {T41b["targets"]}; sampling {T41b["sampling"]}; Proposition 16: minors '
             f'{T41b["prop16_minors_ok"]}, brute force {T41b["prop16_brute"]}, S_T8/S_T8\' one cut {T41b["prop16_one_cut"]}')

        T42b = _s42.check(full_check=True, log=lambda *_x: None)
        gate('G42b Lemma 4*: exhaustive and random arrangements', _s42.check_ok(T42b, True)
             and {k: list(v) for k, v in T42b['coordinate_lines'].items()} == M42['coordinate_lines'],
             f'all 8-tuples of subspaces of GF(q)^2: {T42b["exhaustive_2dim"]}; coordinate lines + all lines '
             f'{T42b["coordinate_lines"]}; random minima {T42b["random"]}; M = 8 over GF(7), GF(3): {T42b["M8"]}')

        T43b = _s43.check(solve=True, log=lambda *_x: None)
        gate('G43b secret sharing: the linear programs re-solved', T43b['ok'],
             'HiGHS optima ' + ', '.join(f'{k} {v["lp"]:.9f} (exact feasible point, t = {v["primal"]})'
                                         for k, v in T43b['cases'].items())
             + '; controls (expect 1) ' + ', '.join(f'{k} {v:.9f}' for k, v in T43b['control_lp'].items())
             + '; with the Jafari-Khazaei inequalities (expect 4/3) '
             + ', '.join(f'{k} {v["lp"]:.9f} (exact feasible point, t = {v["exact"]})' for k, v in T43b['jk_lp'].items()))

        ok47b = _css47.self_test(log=lambda *_x: None)
        gate('G47b CSS search: formula, graph conversion, kernel', ok47b,
             'CSS entropy formula and the bipartite-graph conversion against state vectors (60 random codes); '
             'cc/s41_css_anneal.c realises rows 8 and 25 of the catalogue at lambda = 2, re-verified')

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
