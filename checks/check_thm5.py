"""The five-party gap of the paper (Theorem 5 except its facet statement) on its own: standard library only, no data
files, under a second.

Checks, in exact integer and rational arithmetic and without the code of the repository:
 (1) the 31 facet classes of QLR_5 listed below generate, under the group S_6 of permutations of the five parties and
     the purifier, 10,860 distinct inequalities: the irredundant description of QLR_5 used in the paper (Methods Sec.
     IV I). The 31 class representatives are copied from data/qlr5_H_facets10860.npy of the repository, which this
     script does not read; 31 is the number of facet classes of QLR_5 given by Bao, Cheng, Hernandez-Cuenca and Su.
 (2) the witness r of Theorem 5 (Supplementary Table S2) satisfies all of them, and the 105 that are tight at r have
     rank 30, so r spans an extreme ray of QLR_5;
 (3) inequality (star) of the paper takes the value -1 at r: its left side is 24 and its right side 25;
 (4) the 26 Shannon terms of Supplementary Note 3.1, on the six systems A..F and one more variable Z, all with positive
     coefficients, sum to the balanced six-variable inequality (star-h) of Methods Sec. IV H after the substitutions
     that make Z a common information of AD and BE, g(ADZ) = h(AD), g(BEZ) = h(BE) and g(Z) = h(AD) + h(BE) - h(ABDE),
     and (star-h) pulls back to (star) through the lift h(X) = S(X) + sum of S(i), i in X; so (star) holds whenever the
     lift has such a common information, in particular on the lift of every stabilizer state, where (U_A + U_D) and
     (U_B + U_E) have one (Methods Sec. IV H).
Controls: the ray q2 of Supplementary Note 3.2 (entropy 1 on every system and 2 on every pair) lies in QLR_5, is extreme
and satisfies (star); the 26 terms without the first one, or with Z a common information of AB and DE instead, leave
terms with Z; and h(ABCDEF) pulls back to S(A) + ... + S(F).
The facet statement of Theorem 5 is checked by check_rev_sec4.py, which needs no data either.
Coordinates: S on the 31 nonempty subsets of the parties A..E as bit masks (A = bit 0, ..., E = bit 4), in increasing
order of the mask; the purifier F is bit 5 and Z is bit 6. Ends with ALL OK or FAILURE.
"""
from fractions import Fraction
from itertools import permutations

L = 'ABCDEFZ'
FULL6 = 63
CLASSES = [
    (-3, -1, 2, 0, -2, 0, 2, 0, 1, -1, 0, 0, 0, 0, -2, 0, 2, 1, -1, 0, 1, -1, 0, 0, -1, 1, 0, -1, 1, 2, 0),
    (-2, -2, 1, 0, -1, 0, 1, 0, 1, 1, 0, 0, 0, 0, -1, 0, 1, 1, 0, 0, 0, 0, -1, 0, -1, 0, 0, -2, 2, 2, 0),
    (-2, -2, 2, 0, 0, 0, -1, 0, 1, 1, -1, 0, 0, 0, 0, 0, 1, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 1, 0),
    (-2, -1, 1, 0, -2, 0, 1, 0, 1, 0, 0, 0, 1, 0, -1, 0, 1, 0, 0, 0, 1, 0, -1, 0, -1, 1, 0, -2, 1, 2, 0),
    (-2, -1, 2, -1, 2, -1, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 1, 1, -1, 1, -1, 0, 0, 0, 0, -1, 1, -1, 1, 1, 0),
    (-2, -1, 2, 0, -1, 0, 1, 0, 1, -1, 0, 0, 0, 0, -1, 0, 1, 1, -1, 0, 0, 0, 0, 0, -1, 1, 0, -1, 1, 1, 0),
    (-2, -1, 2, 0, 0, -1, 1, 0, 1, -1, 0, 0, 0, 0, -1, 0, 1, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 1, 0),
    (-2, 0, -3, 0, 0, 0, 1, 0, 1, 0, 2, -1, 1, 0, -1, 0, 1, 0, 2, -1, 1, 0, -1, 0, -1, -2, 0, 2, -1, 2, 0),
    (-2, 0, -2, 0, 1, -1, 1, 0, 1, 0, 1, 0, -1, 0, 0, 0, 1, 0, 1, 0, -1, 0, 0, -2, 1, 0, -1, 2, 0, 1, 0),
    (-1, -1, 0, -1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 1, 1, -1, 0, 0, 0, 0, 0, -1, -1, 1, 0, 1, 1, 0),
    (-1, -1, 0, 0, -1, 0, 1, 0, 0, 0, 1, 0, 0, 0, -1, 0, 1, 1, -1, 0, 0, -1, 0, -1, 1, 1, -1, 0, 1, 1, 0),
    (-1, -1, 1, -1, 1, 1, -1, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 1, 0, 1, 0, -1, 0, 0, 0, -1, 1, 0, 0, 1, 0),
    (-1, -1, 1, 0, -1, -1, 1, 0, 0, 0, -1, 0, 0, 0, -1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0),
    (-1, -1, 1, 0, -1, -1, 1, 0, 0, 0, -1, 0, 0, 0, -1, 0, 0, 1, 0, 0, 1, 1, -1, 0, 1, 0, 0, -1, 1, 1, 0),
    (-1, -1, 1, 0, -1, 0, 1, 0, 0, -1, 1, 0, 0, 0, -1, 0, 0, 1, 0, 0, 1, -1, 0, 0, 1, 0, -1, -1, 1, 1, 0),
    (-1, -1, 1, 0, -1, 0, 1, 0, 0, 0, 0, 0, 0, 0, -1, 0, 1, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, -1, 1, 1, 0),
    (-1, -1, 1, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, -1, -1, 1, 1, 0),
    (-1, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 1, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 1, 0),
    (-1, 0, -2, 0, 0, 0, 1, 0, 0, 0, 1, -2, 2, 0, -1, 0, 0, 0, 2, 0, 1, -1, 0, 0, 1, -1, 0, 1, -2, 1, 0),
    (-1, 0, -2, 0, 0, 0, 1, 0, 0, 0, 1, -1, 1, 0, -1, 0, 1, -1, 2, 0, -1, -1, 0, 0, 0, 1, -1, 1, 0, 1, 0),
    (-1, 0, -1, 0, 0, 0, -1, 0, 0, 0, 1, -1, 1, 0, 0, 0, 0, 0, 1, -1, 1, 0, 0, 0, 1, -1, 0, 2, -2, 1, 0),
    (-1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 1, -1, 1, 0, -1, 0, 1, -1, 1, 0, -1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0),
    (-1, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, -1, 0, 1, 0, 0, 0, 0, 1, -1, -1, 1, 0, 0, 0, 0, 1, 0),
    (-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0),
    (0, 0, -2, 0, -2, -1, 2, 0, 0, 0, 1, 0, 1, 0, -1, 0, 0, 0, 1, 0, 1, 0, -1, -1, 0, 1, 0, 1, 0, 0, 0),
    (0, 0, -2, 0, -1, 0, 2, 0, -1, 0, 2, 0, -1, -1, 0, 0, 0, 0, 1, -1, 1, 1, -1, -1, 1, 1, -1, 1, 0, 0, 0),
    (0, 0, -1, 0, -1, 0, 0, 0, -1, 0, 1, 0, 1, -1, 0, 0, 0, -1, 1, -1, 1, 1, -1, 0, 0, 1, 0, 1, 0, 0, 0),
    (0, 0, -1, 0, -1, 0, 1, 0, -1, 0, 1, 0, 1, -1, -1, 0, 0, -1, 1, 0, 0, 1, 0, 0, 1, 1, -1, 0, 0, 0, 0),
    (0, 0, -1, 0, -1, 0, 1, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1, 0, 1, 0, -1, -1, 0, 1, 0, 1, 0, 0, 0),
    (0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1),
]
R = (2, 2, 3, 1, 3, 3, 4, 1, 3, 3, 4, 2, 4, 3, 3, 2, 3, 3, 3, 3, 4, 3, 3, 3, 3, 4, 3, 3, 3, 3, 2)


def m(s):
    return sum(1 << L.index(c) for c in s)


def coord(M):
    """coordinate of a nonempty proper subset M of the six systems (purity: S(M) = S(complement))."""
    return (FULL6 ^ M if M >> 5 & 1 else M) - 1


def rank_Q(rows):
    M = [[Fraction(x) for x in r] for r in rows if any(r)]
    rank = 0
    for c in range(len(M[0]) if M else 0):
        piv = next((i for i in range(rank, len(M)) if M[i][c] != 0), None)
        if piv is None:
            continue
        M[rank], M[piv] = M[piv], M[rank]
        M[rank] = [x / M[rank][c] for x in M[rank]]
        for i in range(len(M)):
            if i != rank and M[i][c] != 0:
                f = M[i][c]
                M[i] = [a - f * b for a, b in zip(M[i], M[rank])]
        rank += 1
    return rank


ok = True


def check(cond, msg):
    global ok
    ok &= bool(cond)
    print(('ok    ' if cond else 'FAIL  ') + msg)


def value(f, S):
    return sum(a * b for a, b in zip(f, S))


# (1) the facets of QLR_5 from the 31 classes
images = [[coord(sum(1 << p[i] for i in range(6) if M >> i & 1)) for M in range(1, 32)] for p in permutations(range(6))]
facets = set()
for c in CLASSES:
    for img in images:
        g = [0] * 31
        for j in range(31):
            g[img[j]] += c[j]
        facets.add(tuple(g))
facets = sorted(facets)
check(len(CLASSES) == 31 and len(facets) == 10860, '31 classes give %d facet inequalities of QLR_5' % len(facets))


def in_cone_and_extreme(S, name):
    vals = [value(f, S) for f in facets]
    tight = [f for f, v in zip(facets, vals) if v == 0]
    rk = rank_Q(tight)
    check(min(vals) >= 0, '%s satisfies all of them (minimum value %d)' % (name, min(vals)))
    check(rk == 30, '%d of them are tight at %s, of rank %d, so %s spans an extreme ray' % (len(tight), name, rk, name))
    return len(tight)


# (2) r lies in QLR_5 and spans an extreme ray
check(in_cone_and_extreme(R, 'r') == 105, '105 facets are tight at r')

# (3) (star) at r
LHS = ['BCD', 'ABCD', 'ABE', 'BCE', 'BCE', 'ADE', 'ADE', 'CDE']
RHS = ['BC', 'BC', 'ABD', 'ACD', 'AE', 'DE', 'BCDE', 'ABCDE']
star = [0] * 31
for s in LHS:
    star[coord(m(s))] += 1
for s in RHS:
    star[coord(m(s))] -= 1
lv, rv = sum(R[coord(m(s))] for s in LHS), sum(R[coord(m(s))] for s in RHS)
check((lv, rv) == (24, 25) and value(star, R) == -1, '(star) at r: left side %d, right side %d' % (lv, rv))


# (4) the certificate of (star)
def I(a, b, K=''):
    A, B, C = m(a), m(b), m(K)
    out = {}
    for M, c in ((A | C, 1), (B | C, 1), (A | B | C, -1), (C, -1)):
        if M:
            out[M] = out.get(M, 0) + c
    return out


def Hc(a, K):
    return {m(a) | m(K): 1, m(K): -1}


def combine(terms):
    out = {}
    for coef, row in terms:
        for M, c in row.items():
            out[M] = out.get(M, 0) + coef * c
    return {M: c for M, c in out.items() if c}


def eliminate_Z(F, X, Y):
    """g(XZ) = h(X), g(YZ) = h(Y), g(Z) = h(X) + h(Y) - h(XY); None if another term with Z survives."""
    ZB, out = 1 << 6, {}
    rules = {X | ZB: {X: 1}, Y | ZB: {Y: 1}, ZB: {X: 1, Y: 1, X | Y: -1}}
    for M, c in F.items():
        if M & ZB and M not in rules:
            return None
        for MM, cc in (rules[M] if M & ZB else {M: 1}).items():
            out[MM] = out.get(MM, 0) + c * cc
    return {M: c for M, c in out.items() if c}


def pull_back(F):
    """a linear form in the lift h(M) = S(M) + sum of S(i), i in M, as a form in the 31 coordinates of S."""
    if F is None:
        return None
    v = [0] * 31
    for M, c in F.items():
        if M != FULL6:                  # S of all six systems is 0 for a pure state
            v[coord(M)] += c
        for i in range(6):
            if M >> i & 1:
                v[coord(1 << i)] += c
    return v


terms = [(5, Hc('Z', 'ABCDEF')), (2, I('F', 'Z', 'ABCDE'))] + [(1, I(*t)) for t in [
    ('A', 'B', 'Z'), ('A', 'C', 'BFZ'), ('A', 'E', 'BCDZ'), ('A', 'Z', 'BCDEF'), ('A', 'Z', 'BE'), ('B', 'D', 'AEZ'),
    ('B', 'D', 'AZ'), ('B', 'Z', 'ACDEF'), ('B', 'Z', 'AF'), ('C', 'Z', 'BDEF'), ('C', 'Z', 'ADEF'),
    ('C', 'Z', 'ABDE'), ('C', 'Z', 'BE'), ('C', 'Z', 'ABD'), ('D', 'E', 'AFZ'), ('D', 'E', 'BCZ'),
    ('D', 'Z', 'ABCEF'), ('D', 'Z', 'BEF'), ('E', 'F', 'BZ'), ('E', 'Z', 'ABCF'), ('E', 'Z', 'AD'),
    ('F', 'Z', 'AE'), ('F', 'Z', 'AD'), ('F', 'Z', 'BC')]]
check(len(terms) == 26 and all(c > 0 for c, _ in terms), '26 Shannon terms, all with positive coefficients')
STAR_H = {m(s): 1 for s in ('ABCD', 'ABE', 'BCE', 'ADE', 'ABF', 'BCF', 'ADF', 'AEF')}
STAR_H.update({m(s): -1 for s in ('BC', 'ABD', 'AE', 'ABCDE', 'AF', 'ABCF', 'BEF', 'ADEF')})
check(eliminate_Z(combine(terms), m('AD'), m('BE')) == STAR_H,
      'with Z a common information of AD and BE they sum to (star-h) of the Methods')
check(pull_back(STAR_H) == star, '(star-h) pulls back to (star)')

# controls
Q2 = tuple(1 if bin(M).count('1') in (1, 5) else 2 if bin(M).count('1') in (2, 4)
           else 2 if M in [m(s) for s in ('ABC', 'ABD', 'ACD', 'ABE', 'BCE')] else 3 for M in range(1, 32))
check(in_cone_and_extreme(Q2, 'q2') == 295, '295 facets are tight at q2')
check(value(star, Q2) >= 0, '(star) at q2: %d, not violated' % value(star, Q2))
check(eliminate_Z(combine(terms[1:]), m('AD'), m('BE')) is None, 'control: without the first term, a term with Z survives')
check(eliminate_Z(combine(terms), m('AB'), m('DE')) is None,
      'control: with Z a common information of AB and DE, terms with Z survive')
check(pull_back({FULL6: 1}) == [1 if bin(M).count('1') == 1 or M == 31 else 0 for M in range(1, 32)],
      'control: h(ABCDEF) pulls back to S(A) + ... + S(F)')

print('ALL OK' if ok else 'FAILURE')
