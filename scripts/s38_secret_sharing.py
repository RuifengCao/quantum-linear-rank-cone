#!/usr/bin/env python3
"""s38: quantum secret sharing with stabilizer states -- the information ratio depends on the local dimension (R47).

Setting.  Eight systems: the seven nonzero vertices A, B, C, W, X, Y, Z = 100, 010, 001, 110, 101, 011, 111 of the
cube {0,1}^3 (the players, bits 0..6) and the origin R (bit 7), the reference system, maximally entangled with the
secret.  A mixed-state scheme is purified by one more system E; then E is bit 7 and R is bit 8.  For the lift
h(X) = S(X) + sum_{P in X} S(P) of a stabilizer state,
    I(R:A) = h(R) + h(A) - h(RA)            for every set A of players,
so a set A is authorised (I(R:A) = 2 S(R)) iff h(RA) = h(A), and unauthorised (I(R:A) = 0) iff
h(RA) = h(A) + h(R); and S(P)/S(R) = h(P)/h(R).  The information ratio max_P log|P| / log|secret| is at least
max_P S(P)/S(R).

Access structures.  The cube state over GF(2) (the Steane code with a reference qubit) and over GF(p), p odd, are
pure-state schemes with S(P) = S(R) = 1 for every player (ratio 1).  Over GF(2) a set of players is authorised iff it
contains a line of the Fano plane (Gamma_F, the port of AG(3,2)); over GF(p), p odd, the authorised sets are the same
except that WXY is unauthorised and ABCZ authorised (Gamma_NF, the port of R8).

The bound 4/3 (docs/R32_reduction_and_local_dimension.md section 4.15, Proposition 21) is a proof: in normal form the
lift of a stabilizer scheme of prime local dimension p is the entropy function of a GF(p)-linear classical scheme for
the same access structure, contracting a player turns Gamma_F into the Fano access structure F and Gamma_NF (players
A, B, C, Z) into the non-Fano access structure N of Jafari and Khazaei (J. Cryptol. 34 (2021); IACR ePrint 2019/575,
Propositions 4.4 and 4.5), and their bounds give 4/3.  ports() checks the two contraction statements.

Linear programs (the bound 43/41 for every odd local dimension).  Minimise t subject to: h a polymatroid (elemental
inequalities); h self-dual (purity of the state of R, the players and E); the access equalities of the witness;
h(R) = 2; h(P) <= 2t for the players; and Lemma 4 (LB) with t = 2 (the paper's Lemma 1), g . h' >= 0, on the minors
h' of h obtained by removing one element (eight systems) or E and one more element (nine systems), each removed
element contracted or deleted, in every labelling of the remaining seven elements.  (LB) with t = 2 holds for every
subgroup arrangement in an abelian group of odd order, and minors of subgroup arrangements are subgroup arrangements
in groups whose exponent divides the original one, so min t bounds the information ratio of every stabilizer scheme
for Gamma_F of odd local dimension.  The optimum is 43/41, pure and mixed.

Certificates.  Rational multipliers y <= 0 on inequality rows (a.x <= 0) and z on equality rows (e.x = b) with
c - A_ub^T y - A_eq^T z >= 0 componentwise (x >= 0); then t = c.x >= z.b for every feasible x.  They are stored in
data/s38_qss_certificates.json; --verify rebuilds every used row from its descriptor, checks that it is one of the
admissible rows of the case, and checks the identity in exact arithmetic.

  python3 scripts/s38_secret_sharing.py --verify              # exact checks (seconds)
  python3 scripts/s38_secret_sharing.py --solve               # re-solve the linear programs with HiGHS (minutes)
  python3 scripts/s38_secret_sharing.py --certify --out data/s38_qss_certificates.json   # regenerate (needs SymPy)
"""
import argparse, hashlib, itertools, json, os, sys, time
from fractions import Fraction
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..'))
from epr1kit import core
import s27_cube_states as s27
import s28_generalized_cubes as s28

DATA = os.path.join(core.DATA, 's38_qss_certificates.json')
NAMES = 'ABCWXYZ'
# name: (field of the witness cube state, purity constraints, environment system, forms)
CASES = {
    'F_pure_LB2': (2, True, False, ('LB2',)),
    'F_mixed_LB2': (2, True, True, ('LB2',)),
}
# controls (solved with --solve, gate G43b): the form of the other characteristic, or Shannon inequalities only; the
# GF(2) cube state is feasible with t = 1, and the optimum is 1
CONTROLS = {
    'F_pure_LA2': (2, True, False, ('LA2',)),
    'F_pure_shannon': (2, True, False, ()),
}
EXPECTED = {'F_pure_LB2': '43/41', 'F_mixed_LB2': '43/41'}


def forms():
    """(LB) and (LA) of Lemma 4 with t = 2 (variables A_1, A_2, A_3, B_1, B_2, B_3, C = bits 0..6; the programs use
    every labelling), as {mask: coefficient} with g . h >= 0 when valid."""
    LA, LB, n = s28.lemma_forms(2)
    assert n == 7
    return {'LB2': {m: int(LB[m]) for m in range(128) if LB[m]}, 'LA2': {m: int(LA[m]) for m in range(128) if LA[m]}}


def witness(p, env):
    """h = 2 r of the cube state over GF(p): players bits 0..6, R bit 7; with env, E = bit 7 (trivial), R = bit 8."""
    h8 = 2 * s27.rank_function(s27.cube_generator(), p)
    if not env:
        return h8
    h9 = np.zeros(512, dtype=np.int64)
    for X in range(512):
        h9[X] = h8[(X & 127) | ((X >> 8 & 1) << 7)]
    return h9


class Case:
    def __init__(self, name):
        self.name = name
        self.p, self.pure, self.env, self.forms = CASES[name] if name in CASES else CONTROLS[name]
        self.N = 9 if self.env else 8
        self.R = self.N - 1
        self.E = 7 if self.env else None
        self.FULL = (1 << self.N) - 1
        self.T = self.FULL                       # index of t; h(X) at index X - 1
        self.h = witness(self.p, self.env)
        self.g = forms()

    # ---------------------------------------------------------------- rows from descriptors
    def h_row(self, terms, sign=1):
        out = {}
        for X, c in terms.items():
            if X and c:
                out[X - 1] = out.get(X - 1, 0) + sign * c
        return {k: v for k, v in out.items() if v}

    def ub_row(self, desc):
        """an admissible inequality row a with a.x <= 0, or None."""
        kind, N, FULL = desc[0], self.N, self.FULL
        if kind == 'mono':
            i = desc[1]
            return self.h_row({FULL: 1, FULL ^ 1 << i: -1}, -1)
        if kind == 'sub':
            i, j, K = desc[1:]
            if not (0 <= i < j < N) or K & (1 << i | 1 << j) or K & ~FULL:
                return None
            d = {}
            for X, c in ((K | 1 << i, 1), (K | 1 << j, 1), (K | 1 << i | 1 << j, -1), (K, -1)):
                d[X] = d.get(X, 0) + c
            return self.h_row(d, -1)
        if kind == 'share':
            P = desc[1]
            if not 0 <= P < 7:
                return None
            return {(1 << P) - 1: 1, self.T: -2}
        if kind == 'form':
            name, removed, modes, perm = desc[1], tuple(desc[2]), tuple(desc[3]), tuple(desc[4])
            if name not in self.forms or len(removed) != N - 7 or sorted(set(removed)) != list(removed):
                return None
            if self.env and self.E not in removed:
                return None
            rest = [x for x in range(N) if x not in removed]
            if sorted(perm) != rest or any(m not in 'cd' for m in modes) or len(modes) != len(removed):
                return None
            Cm = sum(1 << removed[a] for a in range(len(removed)) if modes[a] == 'c')
            d = {}
            for m, c in self.g[name].items():
                T = 0
                for b in range(7):
                    if m >> b & 1:
                        T |= 1 << perm[b]
                d[T | Cm] = d.get(T | Cm, 0) + c
                if Cm:
                    d[Cm] = d.get(Cm, 0) - c
            return self.h_row(d, -1)
        return None

    def eq_row(self, desc):
        """an admissible equality row (e, b) with e.x = b, or None."""
        kind, R, FULL = desc[0], self.R, self.FULL
        if kind == 'sd':
            if not self.pure:
                return None
            X = desc[1]
            Xc = FULL ^ X
            e = {}
            def add(Y, c):
                if Y:
                    e[Y - 1] = e.get(Y - 1, 0) + c
            add(Xc, Fraction(1)); add(X, Fraction(-1))
            for P in range(self.N):
                add(1 << P, Fraction(-1, 2) if Xc >> P & 1 else Fraction(1, 2))
            return {k: v for k, v in e.items() if v}, Fraction(0)
        if kind in ('auth', 'unauth'):
            A = desc[1]
            if not 0 < A < 128:
                return None
            h = self.h
            if kind == 'auth' and h[A | 1 << R] == h[A]:
                return self.h_row({A | 1 << R: 1, A: -1}), Fraction(0)
            if kind == 'unauth' and h[A | 1 << R] == h[A] + h[1 << R]:
                return self.h_row({A | 1 << R: 1, A: -1, 1 << R: -1}), Fraction(0)
            return None
        if kind == 'norm':
            return {(1 << R) - 1: 1}, Fraction(2)
        return None

    # ---------------------------------------------------------------- exact verification
    def verify(self, cert):
        s = {self.T: Fraction(1)}
        for desc, y in cert['ineq']:
            y = Fraction(y)
            a = self.ub_row(desc)
            if a is None or y > 0:
                return False, None, f'inadmissible inequality row {desc} or positive multiplier'
            for k, v in a.items():
                s[k] = s.get(k, Fraction(0)) - y * v
        bound = Fraction(0)
        for desc, z in cert['eq']:
            z = Fraction(z)
            r = self.eq_row(desc)
            if r is None:
                return False, None, f'inadmissible equality row {desc}'
            e, b = r
            for k, v in e.items():
                s[k] = s.get(k, Fraction(0)) - z * v
            bound += z * b
        neg = [k for k, v in s.items() if v < 0]
        if neg:
            return False, None, f'{len(neg)} negative reduced costs'
        return True, bound, 'ok'

    # ---------------------------------------------------------------- the full linear program
    def all_rows(self):
        N, FULL = self.N, self.FULL
        ub, eq = [], []
        for i in range(N):
            ub.append(('mono', i))
        for i, j in itertools.combinations(range(N), 2):
            rest = [k for k in range(N) if k not in (i, j)]
            for m in range(1 << len(rest)):
                ub.append(('sub', i, j, sum(1 << rest[b] for b in range(len(rest)) if m >> b & 1)))
        for P in range(7):
            ub.append(('share', P))
        removable = [(x,) for x in range(N)] if not self.env else \
            [tuple(sorted((self.E, x))) for x in range(N) if x != self.E]
        for name in self.forms:
            for removed in removable:
                rest = [x for x in range(N) if x not in removed]
                for md in itertools.product('cd', repeat=len(removed)):
                    for perm in itertools.permutations(rest):
                        ub.append(('form', name, removed, md, perm))
        if self.pure:
            for X in range(0, FULL + 1):
                if X <= FULL ^ X:
                    eq.append(('sd', X))
        for A in range(1, 128):
            for kind in ('auth', 'unauth'):
                if self.eq_row((kind, A)) is not None:
                    eq.append((kind, A))
        eq.append(('norm',))
        return ub, eq

    def lp(self):
        from scipy.sparse import coo_matrix
        ubd, eqd = self.all_rows()
        def build(rows):
            r, c, v = [], [], []
            for i, row in enumerate(rows):
                for k, x in row.items():
                    r.append(i); c.append(k); v.append(float(x))
            return coo_matrix((v, (r, c)), shape=(len(rows), self.FULL + 1)).tocsr()
        ubr = [self.ub_row(d) for d in ubd]
        eqr = [self.eq_row(d) for d in eqd]
        return ubd, eqd, ubr, eqr, build(ubr), build([e for e, _b in eqr]), np.array([float(b) for _e, b in eqr])

    def solve(self):
        from scipy.optimize import linprog
        ubd, eqd, ubr, eqr, A_ub, A_eq, b_eq = self.lp()
        c = np.zeros(self.FULL + 1); c[self.T] = 1
        res = linprog(c, A_ub=A_ub, b_ub=np.zeros(A_ub.shape[0]), A_eq=A_eq, b_eq=b_eq,
                      bounds=[(0, None)] * (self.FULL + 1), method='highs-ds')
        return res, (ubd, eqd, ubr, eqr)

    def primal_exact(self, res, rows, den=10 ** 6):
        """Round the HiGHS primal point to fractions with denominators <= den and check it exactly: x >= 0, every
        inequality row a.x <= 0 and every equality row e.x = b of the program.  Returns (feasible, t)."""
        ubd, eqd, ubr, eqr = rows
        x = [Fraction(float(v)).limit_denominator(den) for v in res.x]
        if any(v < 0 for v in x):
            return False, None
        for a in ubr:
            if sum((Fraction(c) * x[k] for k, c in a.items()), Fraction(0)) > 0:
                return False, None
        for e, b in eqr:
            if sum((Fraction(c) * x[k] for k, c in e.items()), Fraction(0)) != b:
                return False, None
        return True, x[self.T]

    def certify(self, tol_y=1e-10, tol_s=1e-9):
        """exact repair of the HiGHS dual (needs SymPy): keep the support, impose zero reduced cost on the tight
        columns and the rational target, and solve exactly for a correction on the pivot unknowns."""
        from sympy.polys.matrices import DomainMatrix
        from sympy import QQ
        res, (ubd, eqd, ubr, eqr) = self.solve()
        tau = Fraction(res.fun).limit_denominator(10 ** 4)
        y0, z0 = res.ineqlin.marginals, res.eqlin.marginals
        S = [i for i, v in enumerate(y0) if v < -tol_y]
        J = list(range(len(eqd)))
        s0 = np.zeros(self.FULL + 1); s0[self.T] = 1.0
        for i in S:
            for k, a in ubr[i].items():
                s0[k] -= y0[i] * float(a)
        for j in J:
            for k, a in eqr[j][0].items():
                s0[k] -= z0[j] * float(a)
        K0 = [k for k in range(self.FULL + 1) if abs(s0[k]) < tol_s]
        kpos = {k: r for r, k in enumerate(K0)}
        unk = [('y', i) for i in S] + [('z', j) for j in J]
        col = {u: n for n, u in enumerate(unk)}
        u0 = [Fraction(float(y0[i])).limit_denominator(10 ** 8) for i in S] + \
             [Fraction(float(z0[j])).limit_denominator(10 ** 8) for j in J]
        nrow = len(K0) + 1
        M = [[Fraction(0)] * len(unk) for _ in range(nrow)]
        for i in S:
            for k, a in ubr[i].items():
                if k in kpos:
                    M[kpos[k]][col[('y', i)]] += Fraction(a)
        for j in J:
            for k, a in eqr[j][0].items():
                if k in kpos:
                    M[kpos[k]][col[('z', j)]] += Fraction(a)
            if eqr[j][1]:
                M[nrow - 1][col[('z', j)]] = eqr[j][1]
        rhs = [Fraction(1) if k == self.T else Fraction(0) for k in K0] + [tau]
        resid = [rhs[r] - sum((M[r][n] * u0[n] for n in range(len(unk)) if M[r][n]), Fraction(0)) for r in range(nrow)]
        aug = DomainMatrix([[QQ(x.numerator, x.denominator) for x in M[r]] + [QQ(resid[r].numerator, resid[r].denominator)]
                            for r in range(nrow)], (nrow, len(unk) + 1), QQ)
        Rr, piv = aug.rref()
        if len(unk) in piv:
            raise ValueError('inconsistent repair system')
        Rm = Rr.to_Matrix()
        u = list(u0)
        for r, pc in enumerate(piv):
            v = Rm[r, len(unk)]
            u[pc] += Fraction(int(v.p), int(v.q))
        cert = dict(ineq=[[list(ubd[i]), str(u[col[('y', i)]])] for i in S if u[col[('y', i)]] != 0],
                    eq=[[list(eqd[j]), str(u[col[('z', j)]])] for j in J if u[col[('z', j)]] != 0])
        cert = json.loads(json.dumps(cert))                      # tuples -> lists, as stored
        ok, bound, msg = self.verify(cert)
        if not ok or bound != tau:
            raise ValueError(f'repair failed: {msg}')
        return tau, cert

    # ---------------------------------------------------------------- facts about the witnesses
    def witness_facts(self):
        """ideal: h(P) = h(R) = 2 for every player; access structure perfect; the own-characteristic form holds on
        every minor and the other one fails (value -2 on the natural minor)."""
        h = self.h
        ideal = all(h[1 << P] == 2 for P in range(7)) and h[1 << self.R] == 2
        kinds = {'auth': 0, 'unauth': 0}
        for A in range(1, 128):
            for k in kinds:
                if self.eq_row((k, A)) is not None:
                    kinds[k] += 1
        return dict(ideal=bool(ideal), auth=kinds['auth'], unauth=kinds['unauth'], perfect=kinds['auth'] + kinds['unauth'] == 127)


def controls():
    """The two cube states against (LB) and (LA) with t = 2 on all seven-element minors (one of the eight systems
    contracted or deleted, all 5040 labellings of the other seven): the minimum of g . h' over these minors h' of
    h = 2 r, for the cube state over GF(2) (Gamma_F) and over GF(3) (Gamma_NF).  Expected: (LB) negative on the GF(2)
    cube state and nonnegative on the GF(3) one, (LA) the other way round.  The negative value shows that the linear
    program excludes the ideal qubit scheme; the nonnegative ones are consistency checks of the validity of the forms
    (they are not used as proofs)."""
    g = forms()
    perms = np.array(list(itertools.permutations(range(7))), dtype=np.int64)
    out = {}
    for p in (2, 3):
        h = witness(p, False)
        for name in ('LB2', 'LA2'):
            mn = None
            for x in range(8):
                rest = np.array([y for y in range(8) if y != x], dtype=np.int64)
                img = rest[perms]                                   # system carrying form variable b, per labelling
                bit = [np.left_shift(1, img[:, b]) for b in range(7)]
                for mode in 'cd':
                    C = (1 << x) if mode == 'c' else 0
                    val = np.zeros(len(perms), dtype=np.int64)
                    for m, c in g[name].items():
                        T = np.zeros(len(perms), dtype=np.int64)
                        for b in range(7):
                            if m >> b & 1:
                                T |= bit[b]
                        val += c * (h[T | C] - h[C])
                    v = int(val.min())
                    mn = v if mn is None else min(mn, v)
            out[f'{name} on GF({p})'] = mn
    return out


def controls_ok(c):
    return c['LB2 on GF(2)'] < 0 and c['LA2 on GF(2)'] >= 0 and c['LA2 on GF(3)'] < 0 and c['LB2 on GF(3)'] >= 0


# The access structures of Jafari and Khazaei (IACR ePrint 2019/575, Section 4): participants 1..6,
# F = 14 + 25 + 36 + 123 + 156 + 246 + 345 (the port of the Fano matroid) and N = F + 456 (a port of the non-Fano matroid).
JK_F = ('14', '25', '36', '123', '156', '246', '345')
JK_N = JK_F + ('456',)
LINES = ('ABW', 'ACX', 'BCY', 'AYZ', 'BXZ', 'CWZ', 'WXY')


def _access(p):
    h = witness(p, False)
    return {A for A in range(128) if h[A | 128] == h[A]}


def _minimal(G):
    return {A for A in G if not any(B != A and B & A == B for B in G)}


def _iso(G, players, minimal_sets):
    """a bijection players -> '123456' that maps the access structure G (masks over the 7 players, here a contraction
    on 6 of them) onto the one generated by minimal_sets, or None."""
    target = set()
    for k in range(7):
        for A in itertools.combinations('123456', k):
            if any(set(m) <= set(A) for m in minimal_sets):
                target.add(frozenset(A))
    for perm in itertools.permutations('123456'):
        f = dict(zip(players, perm))
        if {frozenset(f[i] for i in range(7) if A >> i & 1) for A in G} == target:
            return f
    return None


def ports():
    """Contractions by one player: Gamma_F / q is isomorphic to F for every q; Gamma_NF / q is isomorphic to N for
    q in {A, B, C, Z}, with the triple 456 of N corresponding to the other three of A, B, C, Z (and for q in
    {W, X, Y} it is isomorphic to neither).  Also: the size-3 minimal sets of Gamma_F / q are the Fano lines that miss
    q; the minimal authorised sets of Gamma_F are the seven lines; and Gamma_NF is the same for p = 3, 5, 7."""
    out = dict(F={}, NF={})
    for key, p in (('F', 2), ('NF', 3)):
        G = _access(p)
        for q in range(7):
            players = [i for i in range(7) if i != q]
            Gq = {A ^ (1 << q) for A in G if A >> q & 1}           # the contraction: A u {q} authorised
            fF, fN = _iso(Gq, players, JK_F), _iso(Gq, players, JK_N)
            rec = dict(F=fF is not None, N=fN is not None)
            if fN is not None:
                inv = {v: k for k, v in fN.items()}
                rec['456'] = ''.join(sorted(NAMES[inv[c]] for c in '456'))
            if fF is not None:
                tri = sorted(''.join(NAMES[i] for i in range(7) if A >> i & 1) for A in _minimal(Gq) if bin(A).count('1') == 3)
                rec['triples_are_lines_missing_q'] = tri == sorted(L for L in LINES if NAMES[q] not in L)
            out[key][NAMES[q]] = rec
    lines = sorted(''.join(NAMES[i] for i in range(7) if A >> i & 1) for A in _minimal(_access(2) - {0}))
    out['F_minimal_are_lines'] = lines == sorted(LINES)
    out['NF_same_for_p'] = _access(3) == _access(5) == _access(7)
    return out


def ports_ok(r):
    F, NF = r['F'], r['NF']
    return (all(v['F'] and not v['N'] and v['triples_are_lines_missing_q'] for v in F.values())
            and all(NF[q]['N'] and not NF[q]['F'] and NF[q]['456'] == ''.join(sorted(set('ABCZ') - {q})) for q in 'ABCZ')
            and all(not NF[q]['N'] and not NF[q]['F'] for q in 'WXY')
            and r['F_minimal_are_lines'] and r['NF_same_for_p'])


def jk_lp():
    """The programs of CASES / CONTROLS with the inequalities of Jafari and Khazaei on the
    contractions by one player instead of (LB): for Gamma_F, for every player q and every Fano line L missing q,
    sum_{x in L} (h(xq) - h(q)) >= 4 (h(Rq) - h(q)) (their Proposition 4.4 (II), odd characteristic); for Gamma_NF and
    q in {A, B, C, Z}, the sum over the other three >= 4 (h(Rq) - h(q)) (Proposition 4.5 (III), characteristic 2).
    With or without self-duality the optimum is 4/3: the added inequalities imply t >= 4/3 (some line, or three of
    A, B, C, Z, has h-sum >= 4 h(R) = 8, and h(P) <= 2t), and the rounded HiGHS solution is checked to be exactly
    feasible with t = 4/3.  So purity does not raise their bound in these programs."""
    from scipy.optimize import linprog
    from scipy.sparse import coo_matrix, vstack
    bit = lambda c: 1 << NAMES.index(c)
    out = {}
    for name, p, pure in (('F_pure_JK', 2, True), ('F_classical_JK', 2, False), ('NF_pure_JK', 3, True),
                          ('NF_classical_JK', 3, False)):
        CONTROLS[name] = (p, pure, False, ())
        c = Case(name)
        ubd, eqd, ubr, eqr, A_ub, A_eq, b_eq = c.lp()
        rows = []
        groups = ([(q, L) for q in NAMES for L in LINES if q not in L] if p == 2 else
                  [(q, ''.join(x for x in 'ABCZ' if x != q)) for q in 'ABCZ'])
        for q, L in groups:
            Q = bit(q); d = {}
            for x in L:
                d[bit(x) | Q] = d.get(bit(x) | Q, 0) + 1; d[Q] = d.get(Q, 0) - 1
            d[(1 << c.R) | Q] = d.get((1 << c.R) | Q, 0) - 4; d[Q] = d.get(Q, 0) + 4
            rows.append(c.h_row(d, -1))
        r_, c_, v_ = [], [], []
        for i, row in enumerate(rows):
            for k, x in row.items():
                r_.append(i); c_.append(k); v_.append(float(x))
        A = vstack([A_ub, coo_matrix((v_, (r_, c_)), shape=(len(rows), c.FULL + 1))]).tocsr()
        cv = np.zeros(c.FULL + 1); cv[c.T] = 1
        res = linprog(cv, A_ub=A, b_ub=np.zeros(A.shape[0]), A_eq=A_eq, b_eq=b_eq, bounds=[(0, None)] * (c.FULL + 1),
                      method='highs-ds')
        feas, tx = c.primal_exact(res, (ubd, eqd, list(ubr) + rows, eqr))
        out[name] = dict(lp=float(res.fun), exact=str(tx) if feas else None)
        del CONTROLS[name]
    return out


def access_difference():
    """the authorised sets of Gamma_F and Gamma_NF differ exactly in WXY and ABCZ."""
    hF, hN = witness(2, False), witness(3, False)
    aF = {A for A in range(1, 128) if hF[A | 128] == hF[A]}
    aN = {A for A in range(1, 128) if hN[A | 128] == hN[A]}
    lab = lambda A: ''.join(NAMES[i] for i in range(7) if A >> i & 1)
    lines = {'ABW', 'ACX', 'BCY', 'AYZ', 'BXZ', 'CWZ', 'WXY'}
    contains_line = {A for A in range(1, 128) if any(all(NAMES.index(ch) in [i for i in range(7) if A >> i & 1]
                                                         for ch in L) for L in lines)}
    return dict(F_only=sorted(map(lab, aF - aN)), NF_only=sorted(map(lab, aN - aF)), F_is_lines=aF == contains_line)


def margins():
    """How far the six- and seven-party witnesses are from the cones they leave, on the cut entropies of the pure
    state (purifier folded in).  For a valid form g . S >= 0 and a witness with g . S_w = -v < 0, every S in the cone
    has max_X |S(X) - S_w(X)| >= v / ||g||_1 and Euclidean distance >= v / ||g||_2 (cut coordinates).
    f9: Steane state, every odd local dimension (bits); f12: non-Fano state, qubits (units of log p); DFZ (65) and
    (91) pulled back through the purifier: the cube states (units of log p)."""
    import s34_six_party_survey as s34

    def fold7(F):
        G = np.zeros(128, dtype=np.int64)
        for i, c in enumerate(F):
            if 0 < i < 127:
                G[i if not i >> 6 & 1 else 127 ^ i] += c
        return G

    def h_to_S7(f):
        F = np.zeros(128, dtype=np.int64)
        for X, c in enumerate(f):
            if c and X:
                F[X] += c
                for P in range(7):
                    if X >> P & 1:
                        F[1 << P] += c
        return F
    d = np.load(os.path.join(core.DATA, 's31_six_party.npz'))
    rows = {'f9 (Steane)': (fold7(h_to_S7(d['f9'])), int(d['f9'] @ d['h_F'])),
            'f12 (non-Fano)': (fold7(s34.transfer_form()), -1)}
    for name in ('dfz65_odd', 'dfz91_even'):
        rows[f'{name} (cube)'] = (np.asarray(s27.pullback(s27.load_dfz_char()[name])), -2)
    out = {}
    for k, (g, v) in rows.items():
        l1 = int(np.abs(g).sum())
        out[k] = dict(value=v, l1=l1, eps_inf=str(Fraction(-v, l1)), dist2=float(-v / np.sqrt((g.astype(float) ** 2).sum())))
    return out


def load(path=DATA):
    return json.load(open(path))


def content_sha(path=DATA):
    d = load(path)
    return hashlib.sha256(json.dumps({k: d[k] for k in sorted(d)}, sort_keys=True).encode()).hexdigest()


def check(path=DATA, solve=False, log=print):
    d = load(path)
    out = dict(cases={}, ok=True)
    for name in CASES:
        t0 = time.time()
        c = Case(name)
        ok, bound, msg = c.verify(d[name])
        rec = dict(certified=str(bound) if ok else None, expected=EXPECTED[name], stored=d[name]['bound'],
                   rows=[len(d[name]['ineq']), len(d[name]['eq'])], msg=msg, witness=c.witness_facts())
        rec['ok'] = ok and str(bound) == EXPECTED[name] == d[name]['bound'] and rec['witness']['perfect']
        if solve:
            res, rows = c.solve()
            rec['lp'] = float(res.fun)
            feas, tx = c.primal_exact(res, rows)
            rec['primal'] = str(tx) if feas else None
            rec['ok'] = rec['ok'] and abs(res.fun - float(Fraction(EXPECTED[name]))) < 1e-9 and rec['primal'] == EXPECTED[name]
        out['cases'][name] = rec
        out['ok'] = out['ok'] and rec['ok']
        log(f'  {name}: certified {rec["certified"]} (expected {EXPECTED[name]}), {rec["rows"][0]} inequality and '
            f'{rec["rows"][1]} equality rows' + (f', LP {rec["lp"]:.12f}, exact feasible point with t = {rec["primal"]}'
                                                 if solve else '') + f' [{time.time() - t0:.0f}s]')
    if solve:
        out['control_lp'] = {}
        for name in CONTROLS:
            res, _rows = Case(name).solve()
            out['control_lp'][name] = float(res.fun)
            log(f'  control {name}: LP {res.fun:.12f} (expected 1)')
        out['ok'] = out['ok'] and all(abs(v - 1) < 1e-9 for v in out['control_lp'].values())
    out['access'] = access_difference()
    out['ok'] = out['ok'] and out['access'] == dict(F_only=['WXY'], NF_only=['ABCZ'], F_is_lines=True)
    out['controls'] = controls()
    out['ok'] = out['ok'] and controls_ok(out['controls'])
    out['ports'] = ports()
    out['ok'] = out['ok'] and ports_ok(out['ports'])
    if solve:
        out['jk_lp'] = jk_lp()
        out['ok'] = out['ok'] and all(abs(v['lp'] - 4 / 3) < 1e-9 and v['exact'] == '4/3' for v in out['jk_lp'].values())
    out['margins'] = margins()
    expect = {'f9 (Steane)': '1/952', 'f12 (non-Fano)': '1/129', 'dfz65_odd (cube)': '1/45', 'dfz91_even (cube)': '2/105'}
    out['margins_ok'] = {k: v['eps_inf'] for k, v in out['margins'].items()} == expect
    out['ok'] = out['ok'] and out['margins_ok']
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--verify', action='store_true')
    ap.add_argument('--solve', action='store_true')
    ap.add_argument('--certify', action='store_true')
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    if a.certify:
        out = {}
        for name in CASES:
            t0 = time.time()
            tau, cert = Case(name).certify()
            out[name] = dict(bound=str(tau), **cert)
            print(f'{name}: {tau} ({len(cert["ineq"])} + {len(cert["eq"])} rows) [{time.time() - t0:.0f}s]', flush=True)
        if a.out:
            json.dump(out, open(a.out, 'w'), separators=(',', ':'))
            print('written', a.out)
    if a.verify or a.solve or not a.certify:
        r = check(solve=a.solve)
        print('access structures differ in', r['access'])
        print('minimum over all seven-element minors and labellings:', r['controls'])
        print('contractions by one player: Gamma_F/q = F of Jafari-Khazaei for every q; Gamma_NF/q = N for q in ABCZ:',
              ports_ok(r['ports']), {q: v.get('456') for q, v in r['ports']['NF'].items()})
        if 'jk_lp' in r:
            print('programs with the Jafari-Khazaei inequalities on contractions (optimum, exact feasible point):', r['jk_lp'])
        for k, v in r['margins'].items():
            print(f'margin {k}: value {v["value"]}, ||g||_1 = {v["l1"]}, eps_inf = {v["eps_inf"]}, '
                  f'Euclidean distance >= {v["dist2"]:.5f}')
        print('ALL OK' if r['ok'] else 'FAILED')
        sys.exit(0 if r['ok'] else 1)


if __name__ == '__main__':
    main()
