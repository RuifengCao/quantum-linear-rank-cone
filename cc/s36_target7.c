/* s36: exhaustive search for weighted graphs over GF(q) on seven vertices (systems 0..5 visible, 6 the purifier) whose
 * cut-rank vector equals a given target.  q = 2, 3, 4, 5, 7, 8, 9, 11, 13, 16.  Over GF(q) these cut ranks are the
 * entropies (units of log q) of the GF(q)-linear stabilizer states of seven systems with one qudit each.
 *
 * Target: 63 integers on stdin, entry m-1 = the cut whose visible systems are the bits of m (as in s34 and s35),
 * extended to all 128 subsets by purity.
 *
 * Symmetry: the scaling normal form of s35 (the edge from each vertex to its smallest earlier neighbour has weight 1).
 * Every scaling orbit contains a graph in normal form, so the search is complete.
 *
 * Search: vertices join in the order 0..6.  When vertex k joins with weights w on 0..k-1, a bipartition (X, Y) of
 * 0..k with k in X gets rank A[X, Y] = rank A[X-k, Y] + [w_Y not in the row space of A[X-k, Y]].  Every cut of the
 * final graph that splits 0..k as (X, Y) has rank between this rank and this rank plus the number of vertices still
 * to join, and must equal its target value.  So each bipartition either forces w_Y into the row space (a linear
 * condition on w), forces it out, leaves it free, or ends the branch.  The linear conditions are solved first; only
 * the vectors of their solution space that are in normal form are enumerated, then tested against the other
 * conditions.  After vertex 6 the conditions are exact; every solution is still re-checked from scratch on all
 * 126 proper nonempty subsets.
 *
 * Output: "SOL w01 w02 ... w56" for the first max_print solutions (edge order (0,1), (0,2), ..., (5,6)), one
 * "LEVEL k nodes" line per vertex, and "TOTAL q solutions" (with stop_after > 0 the search ends after that many
 * solutions; used for control targets, which can have very many).
 * Arguments: q [max_print [stop_after]]
 * Build: gcc -O3 -march=native -o s36_target7 s36_target7.c
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define NV 7
static int Q;
static unsigned char ADD[16][16], SUB[16][16], MUL[16][16], INV[16];
static int T[128];
static int LO[NV][128], HI[NV][128];
static unsigned char A[NV][NV];
static long long nodes[NV], nsol;
static int max_print;
static long long stop_after;

static int gf2_mul(int a, int b, int poly, int deg) {
    int r = 0;
    for (int i = 0; i < deg; i++) if (b >> i & 1) r ^= a << i;
    for (int i = 2 * deg - 2; i >= deg; i--) if (r >> i & 1) r ^= poly << (i - deg);
    return r;
}

static int setup_field(int q) {
    Q = q;
    if (q == 2 || q == 3 || q == 5 || q == 7 || q == 11 || q == 13) {
        for (int a = 0; a < q; a++) for (int b = 0; b < q; b++) {
            ADD[a][b] = (a + b) % q; SUB[a][b] = (a - b + q) % q; MUL[a][b] = a * b % q;
        }
    } else if (q == 4 || q == 8 || q == 16) {
        int poly = q == 4 ? 0x7 : q == 8 ? 0xB : 0x13, deg = q == 4 ? 2 : q == 8 ? 3 : 4;
        for (int a = 0; a < q; a++) for (int b = 0; b < q; b++) {
            ADD[a][b] = SUB[a][b] = a ^ b; MUL[a][b] = gf2_mul(a, b, poly, deg);
        }
    } else if (q == 9) {  /* GF(3)[i]/(i^2 + 1); element a0 + 3 a1 stands for a0 + a1 i */
        for (int a = 0; a < 9; a++) for (int b = 0; b < 9; b++) {
            int a0 = a % 3, a1 = a / 3, b0 = b % 3, b1 = b / 3;
            ADD[a][b] = (a0 + b0) % 3 + 3 * ((a1 + b1) % 3);
            SUB[a][b] = (a0 - b0 + 3) % 3 + 3 * ((a1 - b1 + 3) % 3);
            MUL[a][b] = ((a0 * b0 - a1 * b1) % 3 + 3) % 3 + 3 * ((a0 * b1 + a1 * b0) % 3);
        }
    } else return 0;
    /* field axioms: inverses, associativity, distributivity */
    for (int a = 1; a < q; a++) {
        INV[a] = 0;
        for (int b = 1; b < q; b++) if (MUL[a][b] == 1) INV[a] = b;
        if (!INV[a]) return 0;
    }
    for (int a = 0; a < q; a++) for (int b = 0; b < q; b++) for (int c = 0; c < q; c++)
        if (MUL[MUL[a][b]][c] != MUL[a][MUL[b][c]] || MUL[a][ADD[b][c]] != ADD[MUL[a][b]][MUL[a][c]] ||
            ADD[ADD[a][b]][c] != ADD[a][ADD[b][c]] || ADD[SUB[a][b]][b] != a) return 0;
    return 1;
}

/* RREF of the rows `rows` (count nr) and columns `cols` (count nc) of A; R gets the reduced rows (indexed by column
 * position), piv the pivot positions; returns the rank */
static int rref_block(const int *rows, int nr, const int *cols, int nc, unsigned char R[NV][NV], int *piv) {
    for (int i = 0; i < nr; i++) for (int j = 0; j < nc; j++) R[i][j] = A[rows[i]][cols[j]];
    int r = 0;
    for (int c = 0; c < nc && r < nr; c++) {
        int p = -1;
        for (int i = r; i < nr; i++) if (R[i][c]) { p = i; break; }
        if (p < 0) continue;
        if (p != r) for (int j = 0; j < nc; j++) { unsigned char t = R[r][j]; R[r][j] = R[p][j]; R[p][j] = t; }
        unsigned char iv = INV[R[r][c]];
        for (int j = 0; j < nc; j++) R[r][j] = MUL[R[r][j]][iv];
        for (int i = 0; i < nr; i++) if (i != r && R[i][c]) {
            unsigned char f = R[i][c];
            for (int j = 0; j < nc; j++) R[i][j] = SUB[R[i][j]][MUL[f][R[r][j]]];
        }
        piv[r++] = c;
    }
    return r;
}

static int cut_rank_full(int X) {
    int rows[NV], cols[NV], nr = 0, nc = 0, piv[NV];
    unsigned char R[NV][NV];
    for (int v = 0; v < NV; v++) { if (X >> v & 1) rows[nr++] = v; else cols[nc++] = v; }
    return rref_block(rows, nr, cols, nc, R, piv);
}

/* parity checks of the row space of A[Z, Y] in coordinates 0..k-1: rows of `eq` (count returned) */
static int parity_checks(int Z, int Y, int k, unsigned char eq[][NV], int *rank) {
    int rows[NV], cols[NV], nr = 0, nc = 0, piv[NV], ispiv[NV] = {0};
    unsigned char R[NV][NV];
    for (int v = 0; v < k; v++) { if (Z >> v & 1) rows[nr++] = v; else if (Y >> v & 1) cols[nc++] = v; }
    int r = rref_block(rows, nr, cols, nc, R, piv);
    *rank = r;
    for (int i = 0; i < r; i++) ispiv[piv[i]] = 1;
    int ne = 0;
    for (int j = 0; j < nc; j++) if (!ispiv[j]) {
        memset(eq[ne], 0, NV);
        eq[ne][cols[j]] = 1;
        for (int i = 0; i < r; i++) eq[ne][cols[piv[i]]] = SUB[0][R[i][j]];
        ne++;
    }
    return ne;
}

static void dfs(int k);

static void join(int k) {
    /* vertex k joins; collect the conditions on its weights w (coordinates 0..k-1) */
    if (k < 1 || k >= NV) return;
    static unsigned char EQ[NV][1 << NV][NV][NV];   /* per level: OUT groups */
    static int EQn[NV][1 << NV];
    unsigned char lin[64 * NV][NV];
    int nlin = 0, nout = 0, full = (1 << k) - 1;
    for (int Z = 0; Z <= full; Z++) {
        int Y = full & ~Z;
        if (!Y) continue;
        int X = Z | 1 << k, r0;
        unsigned char eq[NV][NV];
        int ne = parity_checks(Z, Y, k, eq, &r0);
        int lo = LO[k][X], hi = HI[k][X];
        int in_ok = r0 >= lo && r0 <= hi, out_ok = r0 + 1 >= lo && r0 + 1 <= hi;
        if (in_ok && out_ok) continue;
        if (!in_ok && !out_ok) return;
        if (in_ok) { for (int e = 0; e < ne; e++) memcpy(lin[nlin++], eq[e], NV); }
        else {
            if (ne == 0) return;  /* w_Y is always in the row space */
            memcpy(EQ[k][nout], eq, sizeof eq); EQn[k][nout] = ne; nout++;
        }
    }
    /* solution space of the linear conditions */
    int piv[NV], np = 0, ispiv[NV] = {0};
    for (int c = 0; c < k && np < nlin; c++) {
        int p = -1;
        for (int i = np; i < nlin; i++) if (lin[i][c]) { p = i; break; }
        if (p < 0) continue;
        if (p != np) { unsigned char t[NV]; memcpy(t, lin[p], NV); memcpy(lin[p], lin[np], NV); memcpy(lin[np], t, NV); }
        unsigned char iv = INV[lin[np][c]];
        for (int j = 0; j < k; j++) lin[np][j] = MUL[lin[np][j]][iv];
        for (int i = 0; i < nlin; i++) if (i != np && lin[i][c]) {
            unsigned char f = lin[i][c];
            for (int j = 0; j < k; j++) lin[i][j] = SUB[lin[i][j]][MUL[f][lin[np][j]]];
        }
        piv[np++] = c; ispiv[c] = 1;
    }
    int fr[NV], d = 0;
    unsigned char B[NV][NV];
    for (int c = 0; c < k; c++) if (!ispiv[c]) fr[d++] = c;
    for (int t = 0; t < d; t++) {
        memset(B[t], 0, NV);
        B[t][fr[t]] = 1;
        for (int i = 0; i < np; i++) B[t][piv[i]] = SUB[0][lin[i][fr[t]]];
    }
    long long total = 1;
    for (int t = 0; t < d; t++) total *= Q;
    for (long long idx = 0; idx < total; idx++) {
        unsigned char w[NV] = {0};
        long long x = idx;
        for (int t = 0; t < d; t++) {
            int c = x % Q; x /= Q;
            if (c) for (int j = 0; j < k; j++) w[j] = ADD[w[j]][MUL[c][B[t][j]]];
        }
        int first = -1;
        for (int j = 0; j < k; j++) if (w[j]) { first = j; break; }
        if (first >= 0 && w[first] != 1) continue;   /* scaling normal form */
        int ok = 1;
        for (int g = 0; g < nout && ok; g++) {
            int any = 0;
            for (int e = 0; e < EQn[k][g] && !any; e++) {
                unsigned char s = 0;
                for (int j = 0; j < k; j++) s = ADD[s][MUL[EQ[k][g][e][j]][w[j]]];
                if (s) any = 1;
            }
            ok = any;
        }
        if (!ok) continue;
        for (int j = 0; j < k; j++) A[k][j] = A[j][k] = w[j];
        dfs(k + 1);
        if (stop_after && nsol >= stop_after) break;
    }
    for (int j = 0; j < k; j++) A[k][j] = A[j][k] = 0;
}

static void dfs(int k) {
    nodes[k - 1]++;
    if (k == NV) {
        for (int X = 1; X < 127; X++) if (cut_rank_full(X) != T[X]) {
            fprintf(stderr, "internal error: a solution failed the final check\n"); exit(2);
        }
        if (nsol < max_print) {
            printf("SOL");
            for (int u = 0; u < NV; u++) for (int v = u + 1; v < NV; v++) printf(" %d", A[u][v]);
            printf("\n");
        }
        nsol++;
        return;
    }
    join(k);
}

int main(int argc, char **argv) {
    if (argc < 2 || !setup_field(atoi(argv[1]))) {
        fprintf(stderr, "usage: s36_target7 q [max_print [stop_after]] < target, q in 2 3 4 5 7 8 9 11 13 16\n");
        return 1;
    }
    max_print = argc > 2 ? atoi(argv[2]) : 3;
    stop_after = argc > 3 ? atoll(argv[3]) : 0;
    int t[64];
    for (int m = 1; m < 64; m++) if (scanf("%d", &t[m]) != 1) { fprintf(stderr, "need 63 integers\n"); return 1; }
    T[0] = T[127] = 0;
    for (int X = 1; X < 127; X++) T[X] = (X >> 6 & 1) ? t[127 ^ X] : t[X];
    /* bounds: vertices 0..k placed; for X within 0..k, min and max of T over the completions, max minus the number
     * of vertices still to come */
    for (int k = 0; k < NV; k++) {
        int placed = (1 << (k + 1)) - 1, rest = 127 & ~placed, rem = NV - 1 - k;
        for (int X = 0; X <= placed; X++) {
            int lo = 1 << 30, hi = -1;
            for (int S = rest;; S = (S - 1) & rest) {
                int v = T[X | S];
                if (v < lo) lo = v;
                if (v > hi) hi = v;
                if (!S) break;
            }
            HI[k][X] = lo; LO[k][X] = hi - rem;
        }
    }
    memset(A, 0, sizeof A);
    dfs(1);
    for (int k = 0; k < NV; k++) printf("LEVEL %d %lld\n", k + 1, nodes[k]);
    printf("TOTAL %d %lld\n", Q, nsol);
    return 0;
}
