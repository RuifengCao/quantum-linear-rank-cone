"""Revised paper, Theorem 6 (Sec. II E, Methods Sec. IV H, Supplementary Note 3.1): the five-party inequality (star) defines a facet
of Stab^(d)_5 for every local dimension d, and of Stab_5.

Six-vertex graph states with weights 0 and 1, one qudit for each of the parties A, ..., E and the purifier F. For such
a state S(X) is the rank of the cut block adj[X, complement] over GF(p), in units of log p. A cut block has at most
three rows or at most three columns, and every minor of a 0/1 matrix of size at most 3 lies in {0, +-1, +-2}; so for
odd p its rank over GF(p) is its rank over Q, and the vectors do not depend on the odd prime p.

Checks:
  (1) the minors of all 0/1 matrices of size 1, 2, 3 lie in {0, +-1, +-2};
  (2) over GF(3), GF(5) and GF(7), the cut ranks of all 2^15 graphs equal the rational cut ranks;
  (3) GF(2): 760 distinct vectors of rank 31; (star) >= 0 on all of them; 140 attain equality, of rank 30;
  (4) odd p (rational ranks): 835 distinct vectors of rank 31, containing the 760 qubit vectors; (star) >= 0 on all
      of them; the same 140 vectors attain equality, of rank 30.
With the validity of (star) on Stab_5 (Theorem 6) and Corollary 2 (Stab^(p)_5 contained in Stab^(d)_5 for p | d),
(3) and (4) give the facet statement for every d. Needs no repository data; runs in about 15 seconds.
"""
import itertools, os, sys
from fractions import Fraction
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import rank_Q, rank_mod_p

SYSTEMS = 'ABCDEF'          # F is the purifier
FIVE = [frozenset(c) for k in range(1, 6) for c in itertools.combinations('ABCDE', k)]
assert len(FIVE) == 31

# (star), transcribed from the paper:
#   S(BCD) + S(ABCD) + S(ABE) + 2 S(BCE) + 2 S(ADE) + S(CDE)
#   >= 2 S(BC) + S(ABD) + S(ACD) + S(AE) + S(DE) + S(BCDE) + S(ABCDE)
LHS = {'BCD': 1, 'ABCD': 1, 'ABE': 1, 'BCE': 2, 'ADE': 2, 'CDE': 1}
RHS = {'BC': 2, 'ABD': 1, 'ACD': 1, 'AE': 1, 'DE': 1, 'BCDE': 1, 'ABCDE': 1}
STAR = [0] * 31
for name, c in LHS.items():
    STAR[FIVE.index(frozenset(name))] += c
for name, c in RHS.items():
    STAR[FIVE.index(frozenset(name))] -= c

failures = []


def ok(cond, msg):
    print(('OK   ' if cond else 'FAIL ') + msg)
    if not cond:
        failures.append(msg)


def det(M):
    n = len(M)
    if n == 1:
        return M[0][0]
    return sum((-1) ** j * M[0][j] * det([row[:j] + row[j + 1:] for row in M[1:]]) for j in range(n))


# (1) minors of small 0/1 matrices
dets = set()
for n in (1, 2, 3):
    for entries in itertools.product((0, 1), repeat=n * n):
        dets.add(det([list(entries[i * n:(i + 1) * n]) for i in range(n)]))
ok(dets <= {0, 1, -1, 2, -2}, 'minors of 0/1 matrices of size <= 3: %s' % sorted(dets))

EDGES = list(itertools.combinations(range(6), 2))
CUTS = []
for X in FIVE:
    I = [SYSTEMS.index(x) for x in sorted(X)]
    J = [j for j in range(6) if j not in I]
    CUTS.append((I, J))

rank_cache = {}


def cut_ranks(adj, p):
    out = []
    for I, J in CUTS:
        block = tuple(tuple(adj[i][j] for j in J) for i in I)
        key = (block, p)
        if key not in rank_cache:
            rows = [list(r) for r in block]
            rank_cache[key] = rank_Q(rows) if p == 0 else rank_mod_p(rows, p)
        out.append(rank_cache[key])
    return tuple(out)


vec2, vecQ = set(), set()
mismatch = {3: 0, 5: 0, 7: 0}
for mask in range(1 << len(EDGES)):
    adj = [[0] * 6 for _ in range(6)]
    for b, (i, j) in enumerate(EDGES):
        if mask >> b & 1:
            adj[i][j] = adj[j][i] = 1
    v2 = cut_ranks(adj, 2)
    vQ = cut_ranks(adj, 0)
    vec2.add(v2)
    vecQ.add(vQ)
    for p in mismatch:
        if cut_ranks(adj, p) != vQ:
            mismatch[p] += 1
ok(all(c == 0 for c in mismatch.values()),
   'cut ranks over GF(3), GF(5), GF(7) equal the rational cut ranks for all 2^15 graphs: mismatches %s' % mismatch)


def star(v):
    return sum(c * x for c, x in zip(STAR, v))


def report(vecs, name, n_expect):
    vecs = sorted(vecs)
    vals = [star(v) for v in vecs]
    tight = [v for v, s in zip(vecs, vals) if s == 0]
    r_all, r_tight = rank_Q([list(v) for v in vecs]), rank_Q([list(v) for v in tight])
    ok(len(vecs) == n_expect and r_all == 31 and min(vals) >= 0 and len(tight) == 140 and r_tight == 30,
       '%s: %d distinct vectors (expect %d), rank %d (expect 31); min (star) = %d; %d on the hyperplane (expect 140), '
       'rank %d (expect 30)' % (name, len(vecs), n_expect, r_all, min(vals), len(tight), r_tight))
    return set(tight)


t2 = report(vec2, 'GF(2)', 760)
tQ = report(vecQ, 'odd p (rational ranks)', 835)
ok(vec2 <= vecQ and t2 == tQ,
   'the 760 qubit vectors are among the 835, and the same 140 vectors attain equality: %s, %s' % (vec2 <= vecQ, t2 == tQ))

print('ALL OK' if not failures else 'FAILURE (%d)' % len(failures))
sys.exit(1 if failures else 0)
