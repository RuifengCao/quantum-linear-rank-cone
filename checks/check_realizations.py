"""Recompute the realization counts of Sec. IV C, Appendix B 2 and Sec. VII for the 641 known orbits of extreme rays
of Stab_5 (data/stab5_extreme_reps.npy).

Every stored graph-state certificate is recomputed by ranks over GF(2), GF(3) or GF(4): the qubit certificates in
data/s16_realisations.npz, results/2026-10-02_a2-ingleton/a2i_q2.npz and data/s24_qubit_l2_certs.npz (together with
the GF(4) graph states they come from), and the qutrit certificates in data/s24_qutrit_certs.npz and
data/s21_qudit_realisations.npz. The script then checks that
- 629 of the 641 orbits have an explicit qubit realization, and the other 12 are among the 58 orbits of Bao et al.
  (arXiv:2006.16292) other than their vector 15; data/qlr59_reps.npy lists their 59 orbits as 40 + 17 + 2, with
  their graph-state vectors 11 and 15 last;
- 96 orbits have an explicit qutrit realization: 16 first realized with qutrits, whose stored qubit realizations are
  all at twice the scale and come from graph states over GF(4), and 80 others, namely all 67 with total
  single-system entropy N <= 11 and 13 of the 47 with N = 12 (N is the sum of the six single-system entries in the
  integer normalization of the catalog); 73 of the 80 at the original scale and 7 at twice the scale; 11 of the 80
  have no explicit qubit realization; all 96 have N <= 16, while the 641 orbits reach N = 30;
- 26 of the 58 orbits of Bao et al. other than their vector 15 are not among the 96, so 122 of the 641 orbits are
  known to lie in the qutrit cone;
- the record of the shorter qutrit search of release R30 (results/2026-10-02_a2-ingleton/a2i_q3.npz.json) covers
  1,969 of the 2,070 orbits then undecided, realizes none of them, and includes the 11 orbits that the qubit search
  of the same run realized.
"""
import itertools, json, sys
import numpy as np
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from common import REPO, graph_state_entropy, masks

# canonical form of a vector on the nonempty subsets of A..E (mask order) under S6 acting on A..F with purity
IDX = np.zeros((720, 31), dtype=np.int64)
for a, p in enumerate(itertools.permutations(range(6))):
    inv = [0] * 6
    for i, j in enumerate(p):
        inv[j] = i
    for M in range(1, 32):
        N = sum(1 << inv[i] for i in range(6) if M >> i & 1)
        IDX[a, M - 1] = (63 ^ N) if N >> 5 & 1 else N   # value of S at N, by purity a subset of A..E (or empty)


def canon(v):
    ext = np.zeros(64, dtype=np.int64)
    ext[1:32] = v
    imgs = ext[IDX]
    return tuple(imgs[np.lexsort(imgs.T[::-1])[0]])


reps = np.load(REPO + '/data/stab5_extreme_reps.npy').astype(np.int64)
row_of = {canon(r): k for k, r in enumerate(reps)}
assert len(row_of) == 641 == len(reps)


def recompute(rep, lam, sizes, adj, p):
    owner = [s for s, c in enumerate(sizes) for _ in range(c)]
    S = graph_state_entropy(adj, owner, 6, p)
    assert [S[M - 1] for M in masks(5)] == [lam * int(x) for x in rep], 'certificate does not reproduce lam * rep'
    k = row_of.get(canon(np.asarray(rep, dtype=np.int64)))
    assert k is not None, 'realized vector is not one of the 641 orbits'
    return k


# GF(4) = {0, 1, a, 1 + a} encoded as 0, 1, 2, 3 (c0 + 2 c1 for c0 + c1 a), with a^2 = a + 1; addition is XOR
MUL4 = [[0, 0, 0, 0], [0, 1, 2, 3], [0, 2, 3, 1], [0, 3, 1, 2]]
INV4 = {1: 1, 2: 3, 3: 2}


def rank_gf4(rows):
    M = [list(r) for r in rows if any(r)]
    r = 0
    ncol = len(M[0]) if M else 0
    for c in range(ncol):
        piv = next((i for i in range(r, len(M)) if M[i][c]), None)
        if piv is None:
            continue
        M[r], M[piv] = M[piv], M[r]
        s = INV4[M[r][c]]
        M[r] = [MUL4[s][x] for x in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c]:
                f = M[i][c]
                M[i] = [x ^ MUL4[f][y] for x, y in zip(M[i], M[r])]
        r += 1
    return r


ok = True
# qubit certificates
qubit = {}   # row -> list of multiplicities lambda
z = np.load(REPO + '/data/s16_realisations.npz')
for k in range(len(z['rep'])):
    sizes = [int(s) for s in z['sizes'][k]]; nv = sum(sizes)
    rows = [int(x) for x in z['rows'][k][:nv]]
    adj = [[(rows[a] >> b) & 1 for b in range(nv)] for a in range(nv)]
    qubit.setdefault(recompute(z['rep'][k], int(z['lam'][k]), sizes, adj, 2), []).append(int(z['lam'][k]))
z = np.load(REPO + '/results/2026-10-02_a2-ingleton/a2i_q2.npz')
r30_qubit_index = [int(i) for i in z['index']]
for k in range(len(z['rep'])):
    sizes = [int(s) for s in z['sizes'][k]]; nv = sum(sizes)
    rows = [int(x) for x in z['rows'][k][:nv]]
    adj = [[(rows[a] >> b) & 1 for b in range(nv)] for a in range(nv)]
    qubit.setdefault(recompute(z['rep'][k], int(z['lam'][k]), sizes, adj, 2), []).append(int(z['lam'][k]))
z = np.load(REPO + '/data/s24_qubit_l2_certs.npz')
gf4_ok = True
for k in range(len(z['rep'])):
    sizes = [int(s) for s in z['sizes'][k]]; nv = sum(sizes)
    adj = z['G'][k][:nv, :nv].astype(int).tolist()
    qubit.setdefault(recompute(z['rep'][k], int(z['lam'][k]), sizes, adj, 2), []).append(int(z['lam'][k]))
    # the GF(4) graph state it comes from: one GF(4) qudit per unit of entropy, cut ranks lam_q * rep
    rep = [int(x) for x in z['rep'][k]]; lq = int(z['lam_q'][k])
    sq = [lq * rep[i] for i in (0, 1, 3, 7, 15)] + [lq * rep[30]]; nq = sum(sq)
    W = z['Wq'][k][:nq, :nq].astype(int).tolist()
    owner = [s for s, c in enumerate(sq) for _ in range(c)]
    gf4_ok &= int(z['q'][k]) == 4 and int(z['lam'][k]) == 2 * lq and sizes == [2 * s for s in sq]
    gf4_ok &= all(W[a][b] == W[b][a] for a in range(nq) for b in range(nq))
    for M in masks(5):
        X = [v for v in range(nq) if M >> owner[v] & 1]; Y = [v for v in range(nq) if not M >> owner[v] & 1]
        gf4_ok &= rank_gf4([[W[a][b] for b in Y] for a in X]) == lq * rep[M - 1]
print('GF(4) graph states behind the 16 qubit certificates at twice the scale reproduce their vectors:', gf4_ok)
ok &= gf4_ok
no_qubit = sorted(set(range(641)) - set(qubit))
print('orbits with an explicit qubit realization (recomputed over GF(2)):', len(qubit), ' without:', len(no_qubit), no_qubit)
ok &= len(qubit) == 629 and len(no_qubit) == 12

# qutrit certificates
rows16, rows80, lam80 = set(), set(), {}
z = np.load(REPO + '/data/s21_qudit_realisations.npz')
for k in range(len(z['rep'])):
    assert int(z['p'][k]) == 3
    sizes = [int(s) for s in z['sizes'][k]]; nv = sum(sizes)
    rows16.add(recompute(z['rep'][k], int(z['lam'][k]), sizes, z['W'][k][:nv, :nv].astype(int).tolist(), 3))
z = np.load(REPO + '/data/s24_qutrit_certs.npz')
for k in range(len(z['rep'])):
    assert int(z['p'][k]) == 3
    sizes = [int(s) for s in z['sizes'][k]]; nv = sum(sizes)
    r = recompute(z['rep'][k], int(z['lam'][k]), sizes, z['G'][k][:nv, :nv].astype(int).tolist(), 3)
    assert r == int(z['stab5_row'][k])
    rows80.add(r)
    lam80[r] = min(lam80.get(r, 99), int(z['lam'][k]))
all96 = rows16 | rows80
print('qutrit realizations (recomputed over GF(3)): first realized with qutrits', len(rows16), ' others', len(rows80),
      ' overlap', len(rows16 & rows80), ' total', len(all96))
ok &= len(rows16) == 16 and len(rows80) == 80 and not rows16 & rows80
min_lam16 = sorted({min(qubit[r]) for r in rows16 if r in qubit})
print('the 16: all have stored qubit realizations:', rows16 <= set(qubit), ' smallest multiplicity among them:', min_lam16)
ok &= rows16 <= set(qubit) and min_lam16 == [2]
Ntot = {k: int(r[0] + r[1] + r[3] + r[7] + r[15] + r[30]) for k, r in enumerate(reps)}   # A, B, C, D, E, F = S(ABCDE)
others = [k for k in range(641) if k not in rows16]
le11 = {k for k in others if Ntot[k] <= 11}
eq12 = {k for k in others if Ntot[k] == 12}
print('others with N <= 11:', len(le11), ' all among the 80:', le11 <= rows80, ';  others with N = 12:', len(eq12),
      ' among the 80:', len(eq12 & rows80), ';  the 80 with N > 12:', len([k for k in rows80 if Ntot[k] > 12]))
ok &= len(le11) == 67 and le11 <= rows80 and len(eq12) == 47 and len(eq12 & rows80) == 13 and rows80 <= le11 | eq12
print('largest N among the 96:', max(Ntot[k] for k in all96), ';  among the 641:', max(Ntot.values()))
ok &= max(Ntot[k] for k in all96) == 16 and max(Ntot.values()) == 30
lams = [lam80[k] for k in rows80]
print('the 80: qutrit multiplicity 1:', lams.count(1), ' 2:', lams.count(2), ';  without explicit qubit realization:',
      len(rows80 - set(qubit)))
ok &= lams.count(1) == 73 and lams.count(2) == 7 and len(rows80 - set(qubit)) == 11

# the 59 orbits of Bao et al.
q59 = np.load(REPO + '/data/qlr59_reps.npy').astype(np.int64)
rows59 = [row_of.get(canon(r)) for r in q59]
assert None not in rows59 and len(set(rows59)) == 59
rows58 = rows59[:58]   # all except their vector 15, the last entry
print('the 12 without explicit qubit realization are among the 58 other than vector 15:', set(no_qubit) <= set(rows58))
ok &= set(no_qubit) <= set(rows58)
not96 = [k for k in rows58 if k not in all96]
print('their graph-state vectors 11 and 15 are among the 96:', rows59[57] in all96 and rows59[58] in all96,
      ';  of the 58 other than vector 15, not among the 96:', len(not96), ';  known in the qutrit cone:', len(all96 | set(rows58)))
ok &= rows59[57] in all96 and rows59[58] in all96 and len(not96) == 26 and len(all96 | set(rows58)) == 122

# record of the shorter qutrit search of release R30
rec = json.load(open(REPO + '/results/2026-10-02_a2-ingleton/a2i_q3.npz.json'))
tried = {int(key.split(':')[0]) for key in rec['best_residual']}
print('R30 qutrit search: p =', rec['p'], ' orbits then undecided:', rec['rays'], ' searched:', len(tried), ' realized:', rec['found'],
      ';  includes the 11 realized by qubits in the same run:', set(r30_qubit_index) <= tried)
ok &= rec['p'] == 3 and rec['rays'] == 2070 and len(tried) == 1969 and rec['found'] == 0 and set(r30_qubit_index) <= tried
print('ALL OK' if ok else 'FAILURE')
