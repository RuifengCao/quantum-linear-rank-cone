#!/usr/bin/env python3
"""check_rev_blind_pm.py -- the inequality of Corollary 2.7 of Pena-Macias (Kybernetika 59 (2023)) with n = 3 is blind
on the self-dual slice (Supplementary Note 5).

Independent of the repository code: nothing from epr1kit/ or scripts/ is imported, no data file is read, and the paper
source is not read.  Requirements: Python 3.8+, numpy, scipy.  Running time: a few seconds.

The inequality, as transcribed from the journal PDF (the PDF itself is not used here), for subspaces A_1, A_2, A_3,
B_1, B_2, B_3 and C of a vector space over a field whose characteristic does not divide n - 1 = 2, with
h(X|Y) = h(XY) - h(Y) and I(X;Y) = h(X) + h(Y) - h(XY):
  h(C) <= (1/4)(h(A_1) + h(A_2) + h(A_3)) + h(C|A_1 A_2 A_3) + sum_i h(C|A_i B_i) + sum_i I(C; A_j, j != i)
          + sum_i I(C; B_i) + sum_i I(C; A_j, j != i, B_i) + I(C; B_1 B_2 B_3).
Written as f(h) >= 0 with f = right side - left side.

What is checked, exactly (integers and fractions.Fraction):
 (a) f on the configuration A_i = <e_i>, B_i = <c - e_i>, C = <c> of (JI), t = 2: -1/4 over GF(2) (the Fano
     configuration) and 3/4 over GF(3) and GF(5); so f is not a consequence of the Shannon inequalities.
 (b) blindness: with (f o D)(h) = f(h*), h*(X) = sum_{i in X} h(i) + h(N - X) - h(N), the form 4(f + f o D) is a
     nonnegative rational combination of the elemental Shannon inequalities on seven variables; the multipliers are
     found by nonnegative least squares in floating point and then verified exactly.  Hence f(h) >= 0 for every
     self-dual polymatroid h, in every labeling of the variables (D commutes with relabelings), and so for the lift of
     every pure quantum state of six parties and a purifier.
"""
import itertools, sys
from fractions import Fraction

import numpy as np
from scipy.optimize import nnls

ok_all = True


def report(name, good, detail=''):
    global ok_all
    ok_all = ok_all and bool(good)
    print(('OK   ' if good else 'FAIL ') + name + (('  [' + detail + ']') if detail else ''), flush=True)


n = 7
FULL = (1 << n) - 1
A, B, C = [0, 1, 2], [3, 4, 5], 6


def mask(*groups):
    m = 0
    for g in groups:
        for i in g:
            m |= 1 << i
    return m


def form():
    """f as a dict {mask: Fraction}, f(h) = sum_m f[m] h(m), with h(empty) = 0"""
    f = {}

    def add(coef, m):
        if m:
            f[m] = f.get(m, Fraction(0)) + Fraction(coef)

    def cond(x, y):                     # h(x|y)
        add(1, mask(x, y)); add(-1, mask(y))

    def mi(x, y):                       # I(x;y)
        add(1, mask(x)); add(1, mask(y)); add(-1, mask(x, y))
    for i in range(3):
        add(Fraction(1, 4), mask([A[i]]))
    cond([C], A)
    for i in range(3):
        cond([C], [A[i], B[i]])
    for i in range(3):
        mi([C], [A[j] for j in range(3) if j != i])
    for i in range(3):
        mi([C], [B[i]])
    for i in range(3):
        mi([C], [A[j] for j in range(3) if j != i] + [B[i]])
    mi([C], B)
    add(-1, mask([C]))
    return {m: c for m, c in f.items() if c}


F = form()


def evaluate(f, h):
    return sum(c * h[m] for m, c in f.items())


def rank_mod_p(rows, p):
    M = [[x % p for x in r] for r in rows]
    r = 0
    ncol = len(M[0]) if M else 0
    for col in range(ncol):
        piv = next((i for i in range(r, len(M)) if M[i][col]), None)
        if piv is None:
            continue
        M[r], M[piv] = M[piv], M[r]
        inv = pow(M[r][col], p - 2, p)
        M[r] = [(x * inv) % p for x in M[r]]
        for i in range(len(M)):
            if i != r and M[i][col]:
                t = M[i][col]
                M[i] = [(x - t * y) % p for x, y in zip(M[i], M[r])]
        r += 1
    return r


print('(a) values on the configuration (JI) with t = 2')
e = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
c = (1, 1, 1)
vals = {}
for p in (2, 3, 5):
    vec = [e[0], e[1], e[2]] + [tuple((c[k] - e[i][k]) % p for k in range(3)) for i in range(3)] + [c]
    h = {m: rank_mod_p([vec[i] for i in range(n) if m >> i & 1], p) for m in range(1, FULL + 1)}
    vals[p] = evaluate(F, h)
report('f = -1/4 over GF(2) (the Fano configuration), 3/4 over GF(3) and GF(5)',
       vals[2] == Fraction(-1, 4) and vals[3] == Fraction(3, 4) and vals[5] == Fraction(3, 4),
       'values %s' % {p: str(v) for p, v in vals.items()})

print('(b) blindness: 4(f + f o D) is a nonnegative combination of elemental Shannon inequalities')
# f o D as a form: f(h*) = sum_m f[m] h*(m)
G = {}
for m, cf in F.items():
    for i in range(n):
        if m >> i & 1:
            G[1 << i] = G.get(1 << i, Fraction(0)) + cf
    comp = FULL & ~m
    if comp:
        G[comp] = G.get(comp, Fraction(0)) + cf
    G[FULL] = G.get(FULL, Fraction(0)) - cf
g = {m: 4 * (F.get(m, Fraction(0)) + G.get(m, Fraction(0))) for m in range(1, FULL + 1)}
g = {m: v for m, v in g.items() if v}
elem = []                                # elemental inequalities as dicts {mask: int}, theta(h) >= 0
for i in range(n):
    elem.append({FULL: 1, FULL & ~(1 << i): -1})
for i, j in itertools.combinations(range(n), 2):
    rest = [k for k in range(n) if k not in (i, j)]
    for r in range(len(rest) + 1):
        for K in itertools.combinations(rest, r):
            k = mask(K)
            t = {}
            for mm, s in ((k | 1 << i, 1), (k | 1 << j, 1), (k | 1 << i | 1 << j, -1), (k, -1)):
                if mm:
                    t[mm] = t.get(mm, 0) + s
            elem.append({mm: s for mm, s in t.items() if s})
Amat = np.zeros((FULL, len(elem)))
for kk, th in enumerate(elem):
    for mm, s in th.items():
        Amat[mm - 1, kk] = s
gvec = np.array([float(g.get(m, 0)) for m in range(1, FULL + 1)])
lam, resid = nnls(Amat, gvec, maxiter=50000)
lamq = [Fraction(x).limit_denominator(1000) if x > 1e-9 else Fraction(0) for x in lam]
comb = {}
for kk, lk in enumerate(lamq):
    if lk:
        for mm, s in elem[kk].items():
            comb[mm] = comb.get(mm, Fraction(0)) + lk * s
comb = {m: v for m, v in comb.items() if v}
nz = sum(1 for x in lamq if x)
report('4(f + f o D) = sum_k c_k theta_k exactly, with rational c_k >= 0 (theta_k the %d elemental inequalities)' % len(elem),
       comb == g and all(x >= 0 for x in lamq),
       '%d nonzero multipliers; denominators %s; floating-point residual %.1e' %
       (nz, sorted({x.denominator for x in lamq if x}), resid))

print('ALL OK' if ok_all else 'FAILURE')
sys.exit(0 if ok_all else 1)
