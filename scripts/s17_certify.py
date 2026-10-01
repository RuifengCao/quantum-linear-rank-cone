#!/usr/bin/env python3
"""s17_certify: exact certificate for a common-information infeasibility found by s17_ci_test.

Given a ray r (31 integers, mask order) and a pair of party sets (X, Y) for which
the LP "extend h_norm(r) by Z = common information of X and Y" is infeasible,
produce an explicit linear inequality valid for every GF(p)-linear polymatroid on
the six parties A..E,F and violated by h_norm(r); translate it into an inequality
F'(S) >= 0 on 31-dimensional pure-state entropy vectors (valid for every
stabilizer state, see the R27 report) and check everything in exact rational
arithmetic:

  sum_k y_k e_k(g) = F(h)   with y_k >= 0 (Fractions), e_k Shannon elemental
                            inequalities on the 7 elements A..F,Z, all Z-terms
                            eliminated by the three CI equalities
                            g(XZ)=h(X), g(YZ)=h(Y), g(Z)=h(X)+h(Y)-h(XY);
  F(h_norm(r)) < 0.

  python3 scripts/s17_certify.py --ray-file rays.npy --index 14 --X 4 --Y 10 --out cert.json
"""
import argparse, json, sys, os
from fractions import Fraction
import numpy as np
from scipy.optimize import linprog

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from s17_ci_test import elemental_rows, h_norm, S6, HIGHS_OPTS   # same HiGHS thread setting as s17_ci_test (see there)

Zb = 1 << 6


def substitute(c, X, Y):
    """map a linear form over 7-element masks to a form over 6-element masks using the CI equalities;
    returns (F, residual) where residual collects Z-masks that cannot be substituted."""
    F = {}
    resid = {}
    for m, v in c.items():
        if v == 0:
            continue
        if not (m & Zb):
            F[m] = F.get(m, 0) + v
        elif m == (X | Zb):
            F[X] = F.get(X, 0) + v
        elif m == (Y | Zb):
            F[Y] = F.get(Y, 0) + v
        elif m == Zb:
            F[X] = F.get(X, 0) + v; F[Y] = F.get(Y, 0) + v; F[X | Y] = F.get(X | Y, 0) - v
        else:
            resid[m] = v
    return {m: v for m, v in F.items() if v != 0}, resid


def exact_solve(Arows, b, guess):
    """solve A y = b over Q; free variables fixed to rationalised guesses. Returns list of Fractions or None."""
    n = len(guess)
    M = [[Fraction(x) for x in row] + [Fraction(bb)] for row, bb in zip(Arows, b)]
    piv_cols = []; r = 0
    for c in range(n):
        p = next((i for i in range(r, len(M)) if M[i][c] != 0), None)
        if p is None:
            continue
        M[r], M[p] = M[p], M[r]
        pv = M[r][c]
        M[r] = [x / pv for x in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c] != 0:
                f = M[i][c]
                M[i] = [a - f * bb for a, bb in zip(M[i], M[r])]
        piv_cols.append(c); r += 1
    for i in range(r, len(M)):
        if M[i][n] != 0:
            return None           # inconsistent
    y = [None] * n
    free = [c for c in range(n) if c not in piv_cols]
    for c in free:
        y[c] = Fraction(guess[c]).limit_denominator(1000)
    for i, c in enumerate(piv_cols):
        y[c] = M[i][n] - sum(M[i][f] * y[f] for f in free)
    return y


def certify(r, X, Y):
    h = h_norm(np.asarray(r, dtype=np.int64))
    E = elemental_rows()
    K = len(E)
    zmasks = [m for m in range(1, 128) if (m & Zb) and m not in (Zb, X | Zb, Y | Zb)]
    # phi_k = contribution of e_k to F(h) after substitution
    phi = []
    for d in E:
        F, _ = substitute(d, X, Y)
        phi.append(sum(v * int(h[m]) for m, v in F.items()))
    Aeq = [[d.get(m, 0) for d in E] for m in zmasks]
    beq = [0] * len(zmasks)
    Aeq.append([1] * K); beq.append(1)                   # normalisation sum y = 1
    res = linprog(np.array(phi, dtype=float), A_eq=np.array(Aeq, dtype=float), b_eq=np.array(beq, dtype=float),
                  bounds=[(0, None)] * K, method='highs', options=HIGHS_OPTS)
    if res.status != 0 or res.fun > -1e-9:
        return None
    ysol = res.x
    supp = [k for k in range(K) if ysol[k] > 1e-10]
    # exact: equalities on the support + phi.y = res.fun rationalised -> scale so that phi.y = -1
    A_s = [[row[k] for k in supp] for row in Aeq[:-1]]
    A_s.append([phi[k] for k in supp]); b_s = [0] * (len(A_s) - 1) + [-1]
    y_s = exact_solve(A_s, b_s, [ysol[k] / (-res.fun) for k in supp])
    if y_s is None or any(v < 0 for v in y_s):
        return None
    y = {supp[i]: y_s[i] for i in range(len(supp)) if y_s[i] != 0}
    # exact recomposition
    c = {}
    for k, yk in y.items():
        for m, v in E[k].items():
            c[m] = c.get(m, 0) + yk * v
    F, resid = substitute(c, X, Y)
    assert not any(v != 0 for v in resid.values()), 'Z-terms not eliminated'
    val = sum(v * int(h[m]) for m, v in F.items())
    assert val < 0, 'not violated'
    # integer scaling
    import math
    den = 1
    for v in list(F.values()) + list(y.values()):
        den = den * v.denominator // math.gcd(den, v.denominator)
    F_int = {m: int(v * den) for m, v in F.items()}
    g = 0
    for v in F_int.values():
        g = math.gcd(g, abs(v))
    F_int = {m: v // g for m, v in F_int.items()}
    y_int = {k: v * den / g for k, v in y.items()}
    return F_int, y_int, val * den / g


def s_form(F_int):
    """coefficients of F'(S) = F(h_norm(S)) over the 31 visible-party coordinates (mask order)."""
    coef = np.zeros(31, dtype=object)
    for i in range(31):
        e = np.zeros(31, dtype=np.int64); e[i] = 1
        h = h_norm(e)
        coef[i] = sum(v * int(h[m]) for m, v in F_int.items())
    return coef


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--ray-file', required=True)
    ap.add_argument('--index', type=int, required=True)
    ap.add_argument('--X', type=int, required=True)
    ap.add_argument('--Y', type=int, required=True)
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    r = np.load(a.ray_file).astype(np.int64)[a.index]
    out = certify(r, a.X, a.Y)
    if out is None:
        raise SystemExit('no exact certificate found (LP feasible or rationalisation failed)')
    F_int, y_int, val = out
    coef = s_form(F_int)
    print('exact certificate OK; F(h_norm(r)) =', val, '; support of y:', len(y_int), 'elemental inequalities')
    print("F'(S) coefficients (mask order):", [int(x) for x in coef])
    print("F'(r) =", int(sum(int(c) * int(v) for c, v in zip(coef, r))))
    if a.out:
        json.dump({'ray': [int(x) for x in r], 'X': a.X, 'Y': a.Y,
                   'F_h': {str(m): v for m, v in sorted(F_int.items())},
                   'y': {str(k): str(v) for k, v in sorted(y_int.items())},
                   'F_S': [int(x) for x in coef]}, open(a.out, 'w'), indent=1)


if __name__ == '__main__':
    main()
