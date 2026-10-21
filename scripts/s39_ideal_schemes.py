#!/usr/bin/env python3
"""s39: ideal stabilizer secret sharing schemes for Gamma_F and Gamma_NF, for every local dimension d (R48 addendum 5).

Statement (docs/R32_reduction_and_local_dimension.md section 4.20, Proposition 30; the paper's Proposition 6 (iv)).
For Z_d qudits, Gamma_F has an ideal stabilizer scheme if and only if d = 2, and Gamma_NF if and only if d is odd.

Setting and notation as in scripts/s38_secret_sharing.py: the players A, B, C, W, X, Y, Z = 100, 010, 001, 110, 101,
011, 111 are bits 0..6 and the reference R (the origin) is bit 7.  A set T of players is authorised iff
h(RT) = h(T) and unauthorised iff h(RT) = h(T) + h(R), for the lift h of the scheme.

Forcing (this script).  Let g be a polymatroid on the eight elements R, A, ..., Z with g(R) = g(P) = 1 for every
player P, g(RT) = g(T) for every authorised set T of players and g(RT) = g(T) + 1 for every unauthorised one.  Then
g = r_M, the rank function of the matroid M whose port at R is the access structure: M = AG(3,2) for Gamma_F and
M = R8 for Gamma_NF.  For every nonempty set X of the eight elements, two linear programs over such g (variables
g(X), unrestricted; elemental Shannon inequalities; the normalisations and the access equalities) give
min g(X) = max g(X) = r_M(X).  Each bound is certified by rational multipliers: y >= 0 on elemental inequalities
a.g >= 0 and z on equalities e.g = b with
    s * 1_X = sum_k y_k a_k + sum_j z_j e_j     (an identity of linear forms on the 255 coordinates),
so that s * g(X) >= sum_j z_j b_j for every feasible g, with s = +1 (lower bound r_M(X)) or s = -1 (upper bound,
sum_j z_j b_j = -r_M(X)).  The certificates are stored in data/s39_ideal_certificates.json; --verify rebuilds every
row from its descriptor and checks every identity in exact arithmetic.

Use in the proof.  The lift of an ideal stabilizer scheme, deleted to R and the players and divided by h(R), is such
a g (every system is maximally mixed with the entropy of the reference).  The rest of the proof is linear algebra over
Z_{2^a} (normal form of the projected stabilizer group; see the note): for M = AG(3,2) it needs a = 1, and for
M = R8 it is impossible.  normal_form() computes the orders of the projections of the normal form for small a, m.

  python3 scripts/s39_ideal_schemes.py --verify                 # exact checks (seconds)
  python3 scripts/s39_ideal_schemes.py --certify --out data/s39_ideal_certificates.json   # regenerate (about 20 s)
"""
import argparse, hashlib, itertools, json, os, sys, time
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', 'data', 's39_ideal_certificates.json')
NAMES = 'ABCWXYZ'
PTS = [(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0), (1, 0, 1), (0, 1, 1), (1, 1, 1), (0, 0, 0)]   # bit 7: R
N, R, FULL = 8, 7, 255
PLAYERS = (1 << 7) - 1


def _rank_mod(rows, p):
    M = [list(r) for r in rows]
    rk, ncol = 0, len(M[0]) if M else 0
    for c in range(ncol):
        piv = next((i for i in range(rk, len(M)) if M[i][c] % p), None)
        if piv is None:
            continue
        M[rk], M[piv] = M[piv], M[rk]
        inv = pow(M[rk][c], -1, p)
        M[rk] = [(x * inv) % p for x in M[rk]]
        for i in range(len(M)):
            if i != rk and M[i][c] % p:
                f = M[i][c]
                M[i] = [(a - f * b) % p for a, b in zip(M[i], M[rk])]
        rk += 1
    return rk


def matroid_rank(case):
    """rank functions (on bitmasks 0..255) of AG(3,2) (case 'F', the vectors (1,x) over GF(2)) and of R8 (case 'NF',
    the same vectors over GF(3), equal to their rank over Q since all minors lie in {0, +-1, +-2})."""
    p = 2 if case == 'F' else 3
    out = [0] * 256
    for X in range(1, 256):
        out[X] = _rank_mod([(1,) + PTS[b] for b in range(8) if X >> b & 1], p)
    return out


def authorised(case):
    r = matroid_rank(case)
    return {T: r[T | 1 << R] == r[T] for T in range(1, 128)}


# ---------------------------------------------------------------------------------------------- rows from descriptors
def ineq_row(desc):
    """elemental Shannon inequality a.g >= 0 as {X: coeff} (X nonzero bitmask), or None if inadmissible."""
    if desc[0] == 'mono':
        i = desc[1]
        if not 0 <= i < N:
            return None
        return {FULL: 1, FULL ^ 1 << i: -1}
    if desc[0] == 'sub':
        i, j, K = desc[1:]
        if not (0 <= i < j < N) or K & (1 << i | 1 << j) or K & ~FULL:
            return None
        d = {}
        for X, c in ((K | 1 << i, 1), (K | 1 << j, 1), (K | 1 << i | 1 << j, -1), (K, -1)):
            if X:
                d[X] = d.get(X, 0) + c
        return {X: c for X, c in d.items() if c}
    return None


def eq_row(desc, auth):
    """equality e.g = b as ({X: coeff}, b), or None if inadmissible for the access structure auth."""
    if desc[0] == 'norm':
        i = desc[1]
        if not 0 <= i < N:
            return None
        return {1 << i: 1}, 1
    if desc[0] in ('auth', 'unauth'):
        T = desc[1]
        if not 0 < T < 128 or auth[T] != (desc[0] == 'auth'):
            return None
        return {T | 1 << R: 1, T: -1}, (0 if desc[0] == 'auth' else 1)
    return None


def all_rows(auth):
    ineq = [('mono', i) for i in range(N)]
    for i, j in itertools.combinations(range(N), 2):
        rest = [k for k in range(N) if k not in (i, j)]
        for m in range(1 << len(rest)):
            ineq.append(('sub', i, j, sum(1 << rest[b] for b in range(len(rest)) if m >> b & 1)))
    eq = [('norm', i) for i in range(N)] + [('auth' if auth[T] else 'unauth', T) for T in range(1, 128)]
    return ineq, eq


# ---------------------------------------------------------------------------------------------- exact verification
def verify_cert(cert, X, sense, value, auth):
    """True iff the certificate proves sense * g(X) >= sense * value on the feasible set (sense = +1 or -1)."""
    s = {}
    for desc, y in cert['ineq']:
        y = Fraction(y)
        a = ineq_row(tuple(desc))
        if a is None or y < 0:
            return False
        for k, v in a.items():
            s[k] = s.get(k, Fraction(0)) + y * v
    bound = Fraction(0)
    for desc, z in cert['eq']:
        z = Fraction(z)
        r = eq_row(tuple(desc), auth)
        if r is None:
            return False
        e, b = r
        for k, v in e.items():
            s[k] = s.get(k, Fraction(0)) + z * v
        bound += z * b
    target = {X: Fraction(sense)}
    if {k: v for k, v in s.items() if v} != target:
        return False
    return bound == sense * value


def feasible_rank_function(case):
    """r_M itself satisfies every row (exactly), so min g(X) <= r_M(X) <= max g(X)."""
    r, auth = matroid_rank(case), authorised(case)
    ineq, eq = all_rows(auth)
    for d in ineq:
        if sum(c * r[X] for X, c in ineq_row(d).items()) < 0:
            return False
    for d in eq:
        e, b = eq_row(d, auth)
        if sum(c * r[X] for X, c in e.items()) != b:
            return False
    return True


# ---------------------------------------------------------------------------------------------- certification
def certify_case(case, log=print):
    import numpy as np
    from scipy.optimize import linprog
    from scipy.sparse import coo_matrix
    from sympy.polys.matrices import DomainMatrix
    from sympy import QQ
    r, auth = matroid_rank(case), authorised(case)
    ineq, eq = all_rows(auth)
    irow = [ineq_row(d) for d in ineq]
    erow = [eq_row(d, auth) for d in eq]

    def mat(rows):
        ri, ci, vi = [], [], []
        for i, row in enumerate(rows):
            for X, c in row.items():
                ri.append(i); ci.append(X - 1); vi.append(float(c))
        return coo_matrix((vi, (ri, ci)), shape=(len(rows), FULL)).tocsr()
    A = mat(irow)
    E = mat([e for e, _ in erow])
    b = np.array([float(bb) for _, bb in erow])
    out = {}
    for X in range(1, 256):
        out[str(X)] = {'set': ''.join(NAMES[k] if k < 7 else 'R' for k in range(8) if X >> k & 1), 'value': r[X]}
        for sense, key in ((1, 'lower'), (-1, 'upper')):
            c = np.zeros(FULL); c[X - 1] = sense
            res = linprog(c, A_ub=-A, b_ub=np.zeros(A.shape[0]), A_eq=E, b_eq=b, bounds=[(None, None)] * FULL,
                          method='highs')
            if res.status != 0 or abs(res.fun - sense * r[X]) > 1e-7:
                raise ValueError(f'{case} X={X} {key}: LP value {res.fun} != {sense * r[X]}')
            y0 = -res.ineqlin.marginals            # >= 0 for the rows a.g >= 0
            z0 = res.eqlin.marginals
            S = [i for i, v in enumerate(y0) if v > 1e-10]
            # exact repair on the support: sum_S y a + sum_j z e = sense * 1_X, sum_j z b = sense * r[X]
            unk = [('y', i) for i in S] + [('z', j) for j in range(len(eq))]
            cols = {u: n for n, u in enumerate(unk)}
            rowsX = sorted({k for i in S for k in irow[i]} | {k for e, _ in erow for k in e} | {X})
            rpos = {k: n for n, k in enumerate(rowsX)}
            M = [[Fraction(0)] * len(unk) for _ in range(len(rowsX) + 1)]
            for i in S:
                for k, v in irow[i].items():
                    M[rpos[k]][cols[('y', i)]] += v
            for j, (e, bb) in enumerate(erow):
                for k, v in e.items():
                    M[rpos[k]][cols[('z', j)]] += v
                M[len(rowsX)][cols[('z', j)]] = Fraction(bb)
            rhs = [Fraction(sense if k == X else 0) for k in rowsX] + [Fraction(sense * r[X])]
            u0 = [Fraction(float(y0[i])).limit_denominator(10 ** 6) for i in S] + \
                 [Fraction(float(z0[j])).limit_denominator(10 ** 6) for j in range(len(eq))]
            resid = [rhs[q] - sum((M[q][n] * u0[n] for n in range(len(unk)) if M[q][n]), Fraction(0))
                     for q in range(len(M))]
            if any(resid):
                aug = DomainMatrix([[QQ(x.numerator, x.denominator) for x in M[q]] +
                                    [QQ(resid[q].numerator, resid[q].denominator)] for q in range(len(M))],
                                   (len(M), len(unk) + 1), QQ)
                Rr, piv = aug.rref()
                if len(unk) in piv:
                    raise ValueError(f'{case} X={X} {key}: inconsistent repair system')
                Rm = Rr.to_Matrix()
                for q, pc in enumerate(piv):
                    v = Rm[q, len(unk)]
                    u0[pc] += Fraction(int(v.p), int(v.q))
            cert = {'ineq': [[list(ineq[i]), str(u0[cols[('y', i)]])] for i in S if u0[cols[('y', i)]] != 0],
                    'eq': [[list(eq[j]), str(u0[cols[('z', j)]])] for j in range(len(eq)) if u0[cols[('z', j)]] != 0]}
            cert = json.loads(json.dumps(cert))
            if not verify_cert(cert, X, sense, r[X], auth):
                raise ValueError(f'{case} X={X} {key}: certificate does not verify')
            out[str(X)][key] = cert
        if X % 32 == 0:
            log(f'  {case}: {X}/255 sets certified')
    return out


# ---------------------------------------------------------------------------------------------- the normal form
def normal_form(a, m=1):
    """orders (as exponents of 2) of the projections pr_T L'' of the normal form
    L'' = {(u, v, w, u+v, u+w, v+w, u+v+w) : u, v, w in Z_{2^a}^{2m}} for the sets T of players, against
    2ma (r_M(RT) - 1) for M = AG(3,2) and M = R8.  Only the sets whose order is not that of AG(3,2) are returned."""
    q = 2 ** a
    blocks = [(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0), (1, 0, 1), (0, 1, 1), (1, 1, 1)]   # coefficients of u, v, w
    rF, rNF = matroid_rank('F'), matroid_rank('NF')
    out = []
    for T in range(1, 128):
        # the projection is the image of (Z_q^{2m})^3 under the 0/1 matrix with rows blocks[P], P in T, tensored with
        # the identity; its order is the product of the orders of the images of Z_q^3 (Smith normal form), to the
        # power 2m
        rows = [blocks[P] for P in range(7) if T >> P & 1]
        order = _image_order(rows, q)
        e = 2 * m * order
        if e != 2 * m * a * (rF[T | 1 << R] - 1):
            out.append((''.join(NAMES[P] for P in range(7) if T >> P & 1), e, 2 * m * a * (rF[T | 1 << R] - 1),
                        2 * m * a * (rNF[T | 1 << R] - 1)))
    return out


def _image_order(rows, q):
    """log2 of the order of the image of Z_q^3 under the linear map x -> (r . x)_{r in rows}, q a power of two."""
    import itertools as it
    img = set()
    for x in it.product(range(q), repeat=3):
        img.add(tuple(sum(r[k] * x[k] for k in range(3)) % q for r in rows))
    n = len(img)
    return n.bit_length() - 1


# ---------------------------------------------------------------------------------------------- driver
def content_sha(path=DATA):
    with open(path, 'rb') as f:
        return hashlib.sha256(f.read()).hexdigest()


def check(path=DATA, log=print):
    with open(path) as f:
        data = json.load(f)
    ok = True
    for case in ('F', 'NF'):
        auth = authorised(case)
        r = matroid_rank(case)
        good = all(verify_cert(data[case][str(X)][key], X, s, r[X], auth)
                   for X in range(1, 256) for key, s in (('lower', 1), ('upper', -1)))
        feas = feasible_rank_function(case)
        log(f'{case}: {"OK" if good and feas else "FAIL"}  all 255 values forced to the rank function of '
            f'{"AG(3,2)" if case == "F" else "R8"} (510 exact certificates; the rank function is feasible)')
        ok &= good and feas
    for a in (1, 2, 3):
        diff = normal_form(a)
        want = [] if a == 1 else [('WXY', 2 * (3 * a - 1), 4 * a, 6 * a)]
        good = diff == want
        log(f'normal form over Z_{2 ** a} (m = 1): {"OK" if good else "FAIL"}  sets whose order differs from AG(3,2): '
            f'{diff if diff else "none"}')
        ok &= good
    return ok


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--verify', action='store_true')
    ap.add_argument('--certify', action='store_true')
    ap.add_argument('--out', default=DATA)
    a = ap.parse_args()
    if a.certify:
        t0 = time.time()
        data = {case: certify_case(case) for case in ('F', 'NF')}
        with open(a.out, 'w') as f:
            json.dump(data, f, separators=(',', ':'))
            f.write('\n')
        print(f'wrote {a.out} ({time.time() - t0:.0f}s), sha256 {content_sha(a.out)}')
    if a.verify or not a.certify:
        ok = check(a.out if a.certify else DATA)
        print('ALL OK' if ok else 'FAILED')
        sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
