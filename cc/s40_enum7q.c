/* s40: every entropy vector of a GF(q)-linear stabilizer state of seven systems with one qudit each (R45;
 * docs/R32_reduction_and_local_dimension.md, section 4.10, "Other fields"): the cut-rank vectors of all weighted graphs
 * over GF(q) on seven vertices (0..5 visible, 6 the purifier), q in 2 3 4 5 7 8 9.  Field elements are encoded as in
 * s39_sample7q.c (base-p digits of the coefficients modulo the same irreducible polynomials).
 *
 * Reductions.  Each keeps at least one graph of every class, so the output, closed under S_7, is the complete list.
 *   Scaling: a vertex times a unit keeps every cut rank; the weights from vertex k to 0..k-1 are zero or have their
 *     first nonzero entry 1 (the normal form of s35_enum7.c and s36_target7.c).
 *   Relabelling (skipped with the argument "all"): vertex 6 has the largest degree of the seven, and the degrees of
 *     0..5 in the graph on 0..5 do not increase.  Every graph can be relabelled so (take a vertex of largest degree
 *     as 6, sort the others by their degree without it), and relabelling permutes the vector by an element of S_7.
 *     So the output is then a set of vectors meeting every S_7 orbit, not closed under S_7: close it before comparing
 *     lists (scripts/s36_one_qudit_fields.py does).
 * Method.  For the graph A6 on 0..5 and a cut whose visible part is X (Y = the rest of 0..5),
 *     S(X) = rank A[X, Y + 6] = rank A6[X, Y] + [a_X not in the column space of A6[X, Y]],
 * a = the weights of vertex 6.  For each graph on 0..5 the column spaces that are not all of GF(q)^X are enumerated
 * once, as bitsets over GF(q)^6 (vector index = sum a_j q^j); every admissible a then costs one bit test per such cut.
 * The cut ranks are at most 3 and are packed two bits per cut.
 * Output (stdout): "V t1 ... t63 w01 ... w56" per distinct vector (entry m-1 = the cut whose visible systems are the
 * bits of m; the 21 weights of one graph realising it, edge order (0,1), (0,2), ..., (5,6)), then
 * "TOTAL q prefixes graphs distinct" (prefixes = graphs on 0..5 kept, graphs = seven-vertex graphs examined).
 * Slices: with "slice s n" only the weight vectors of vertex 5 with index = s mod n are taken (to split a long run
 * over machines; the union of the n slices is the full run; slices differ in size).
 * Usage: s40_enum7q q [threads] [all|reduced] [slice s n]   (default: one thread, reduced, no slicing)
 * Build: gcc -O3 -march=native -fopenmp -o s40_enum7q s40_enum7q.c
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif

static int Q, PCH, DEG;
static uint8_t ADD[16][16], SUB[16][16], MUL[16][16], INV[16];
static int modpoly[4];   /* c_0..c_{k-1} of the monic irreducible x^k + c_{k-1} x^{k-1} + ... + c_0 */
static void digits(int a, int *d) { for (int i = 0; i < DEG; i++) { d[i] = a % PCH; a /= PCH; } }
static int undigits(const int *d) { int a = 0; for (int i = DEG - 1; i >= 0; i--) a = a * PCH + d[i]; return a; }
static int fmul(int a, int b) {
    int da[4], db[4], prod[8] = {0};
    digits(a, da); digits(b, db);
    for (int i = 0; i < DEG; i++) for (int j = 0; j < DEG; j++) prod[i + j] = (prod[i + j] + da[i] * db[j]) % PCH;
    for (int t = 2 * DEG - 2; t >= DEG; t--) {
        int c = prod[t];
        if (!c) continue;
        prod[t] = 0;
        for (int i = 0; i < DEG; i++) prod[t - DEG + i] = ((prod[t - DEG + i] - c * modpoly[i]) % PCH + PCH) % PCH;
    }
    return undigits(prod);
}
static int setup(int q) {
    Q = q;
    if (q == 2 || q == 3 || q == 5 || q == 7) { PCH = q; DEG = 1; }
    else if (q == 4) { PCH = 2; DEG = 2; modpoly[0] = 1; modpoly[1] = 1; }               /* x^2 + x + 1 */
    else if (q == 8) { PCH = 2; DEG = 3; modpoly[0] = 1; modpoly[1] = 1; modpoly[2] = 0; } /* x^3 + x + 1 */
    else if (q == 9) { PCH = 3; DEG = 2; modpoly[0] = 1; modpoly[1] = 0; }               /* x^2 + 1 */
    else return 0;
    for (int a = 0; a < q; a++) for (int b = 0; b < q; b++) {
        int da[4], db[4], s[4], d[4];
        digits(a, da); digits(b, db);
        for (int i = 0; i < DEG; i++) { s[i] = (da[i] + db[i]) % PCH; d[i] = (da[i] - db[i] + PCH) % PCH; }
        ADD[a][b] = undigits(s); SUB[a][b] = undigits(d);
        MUL[a][b] = DEG == 1 ? a * b % PCH : fmul(a, b);
    }
    for (int a = 1; a < q; a++) for (int b = 1; b < q; b++) if (MUL[a][b] == 1) INV[a] = b;
    for (int a = 1; a < q; a++) if (!INV[a]) return 0;
    for (int a = 0; a < q; a++) for (int b = 0; b < q; b++) for (int c = 0; c < q; c++)
        if (MUL[a][MUL[b][c]] != MUL[MUL[a][b]][c] || MUL[a][ADD[b][c]] != ADD[MUL[a][b]][MUL[a][c]]) return 0;
    return 1;
}

static int QP[7], Q6;
static uint8_t *NFW[7];             /* normal-form weight vectors of vertex k to 0..k-1 (k entries, stride 6) */
static int NNF[7];
static int NA;                      /* candidates for vertex 6: normal forms of GF(q)^6, sorted by support size, desc. */
static uint8_t *AV;                 /* their entries (stride 6) */
static int *AZ, *AS;                /* support size, support mask */
static uint32_t *PROJ;              /* PROJ[i * 64 + m] = index of the restriction of candidate i to the bits of m */

static int vindex(const uint8_t *v) { int x = 0; for (int j = 5; j >= 0; j--) x = x * Q + v[j]; return x; }

static void normal_forms(void) {
    QP[0] = 1;
    for (int j = 1; j <= 6; j++) QP[j] = QP[j - 1] * Q;
    Q6 = QP[6];
    for (int k = 1; k <= 6; k++) {
        NFW[k] = calloc((size_t)QP[k] * 6, 1);
        NNF[k] = 0;
        for (int x = 0; x < QP[k]; x++) {
            uint8_t e[6] = {0};
            int y = x, first = -1;
            for (int j = 0; j < k; j++) { e[j] = (uint8_t)(y % Q); y /= Q; if (e[j] && first < 0) first = j; }
            if (first >= 0 && e[first] != 1) continue;
            memcpy(NFW[k] + (size_t)NNF[k] * 6, e, 6);
            NNF[k]++;
        }
    }
    NA = NNF[6];
    AV = malloc((size_t)NA * 6); AZ = malloc(sizeof(int) * NA); AS = malloc(sizeof(int) * NA);
    int n = 0;
    for (int z = 6; z >= 0; z--)
        for (int i = 0; i < NNF[6]; i++) {
            const uint8_t *e = NFW[6] + (size_t)i * 6;
            int s = 0, c = 0;
            for (int j = 0; j < 6; j++) if (e[j]) { s |= 1 << j; c++; }
            if (c != z) continue;
            memcpy(AV + (size_t)n * 6, e, 6); AZ[n] = c; AS[n] = s; n++;
        }
    PROJ = malloc(sizeof(uint32_t) * (size_t)NA * 64);
    for (int i = 0; i < NA; i++)
        for (int m = 0; m < 64; m++) {
            uint8_t e[6];
            for (int j = 0; j < 6; j++) e[j] = (m >> j & 1) ? AV[(size_t)i * 6 + j] : 0;
            PROJ[(size_t)i * 64 + m] = (uint32_t)vindex(e);
        }
}

/* --------------------------------------------------------------------------------------------- hash table of vectors */
typedef struct { uint64_t k0, k1; uint8_t w[21]; uint8_t used; } entry;
typedef struct { entry *t; uint64_t cap, n; } table;
static uint64_t mix(uint64_t a, uint64_t b) {
    uint64_t h = a * 0x9E3779B97F4A7C15ULL ^ (b + 0x632BE59BD9B4E019ULL + (a << 6) + (a >> 2));
    h ^= h >> 29; h *= 0xBF58476D1CE4E5B9ULL; h ^= h >> 32;
    return h;
}
static void tinit(table *T, uint64_t cap) { T->cap = cap; T->n = 0; T->t = calloc(cap, sizeof(entry)); }
static void tadd(table *T, uint64_t k0, uint64_t k1, const uint8_t *w);
static void tgrow(table *T) {
    table U; tinit(&U, T->cap * 2);
    for (uint64_t i = 0; i < T->cap; i++) if (T->t[i].used) tadd(&U, T->t[i].k0, T->t[i].k1, T->t[i].w);
    free(T->t); *T = U;
}
static void tadd(table *T, uint64_t k0, uint64_t k1, const uint8_t *w) {
    uint64_t h = mix(k0, k1) & (T->cap - 1);
    while (T->t[h].used) {
        if (T->t[h].k0 == k0 && T->t[h].k1 == k1) return;
        h = (h + 1) & (T->cap - 1);
    }
    T->t[h].used = 1; T->t[h].k0 = k0; T->t[h].k1 = k1; memcpy(T->t[h].w, w, 21);
    if (++T->n * 2 > T->cap) tgrow(T);
}

/* --------------------------------------------------------------------------------------------- one graph on 0..5 */
#define LCAP 4096
typedef struct {
    uint64_t *bits;                 /* 64 bitsets of Q6 bits (index m) */
    int W;
    int *elems;                     /* the indices set, to clear them afterwards */
    uint64_t lkey[LCAP]; uint32_t lgen[LCAP]; int lfirst[LCAP]; int lslot[LCAP]; int nl; uint32_t gen;
    table T;
    uint64_t prefixes, graphs;
} work;

static void process(work *wk, uint8_t A[6][6], const int *deg, int all) {
    int act[64], from[65], nact = 0, nel = 0;
    uint64_t base0 = 0, base1 = 0;
    for (int m = 1; m < 64; m++) {
        int nx = 0;
        uint8_t R[6][6];            /* the columns A6[X, y], y in Y, as vectors of length 6 (zero off X) */
        int nr = 0;
        for (int j = 0; j < 6; j++) if (m >> j & 1) nx++;
        for (int y = 0; y < 6; y++) {
            if (m >> y & 1) continue;
            for (int j = 0; j < 6; j++) R[nr][j] = (m >> j & 1) ? A[j][y] : 0;
            nr++;
        }
        int r = 0;                  /* row-reduce R: the nonzero rows span the column space of A6[X, Y] */
        for (int c = 0; c < 6 && r < nr; c++) {
            if (!(m >> c & 1)) continue;
            int piv = -1;
            for (int i = r; i < nr; i++) if (R[i][c]) { piv = i; break; }
            if (piv < 0) continue;
            if (piv != r) for (int j = 0; j < 6; j++) { uint8_t t = R[r][j]; R[r][j] = R[piv][j]; R[piv][j] = t; }
            uint8_t iv = INV[R[r][c]];
            for (int j = 0; j < 6; j++) R[r][j] = MUL[iv][R[r][j]];
            for (int i = 0; i < nr; i++) if (i != r && R[i][c]) {
                uint8_t f = R[i][c];
                for (int j = 0; j < 6; j++) R[i][j] = SUB[R[i][j]][MUL[f][R[r][j]]];
            }
            r++;
        }
        if (m - 1 < 32) base0 |= (uint64_t)r << (2 * (m - 1)); else base1 |= (uint64_t)r << (2 * (m - 33));
        if (r == nx) continue;      /* the column space is all of GF(q)^X: a_X is always in it */
        uint64_t *bs = wk->bits + (size_t)m * wk->W;
        from[nact] = nel;           /* enumerate the span of R[0..r-1] */
        uint8_t span[729][6];
        int ns = 1;
        memset(span[0], 0, 6);
        for (int t = 0; t < r; t++) {
            int cur = ns;
            for (int c = 1; c < Q; c++)
                for (int e = 0; e < cur; e++) {
                    for (int j = 0; j < 6; j++) span[ns][j] = ADD[span[e][j]][MUL[c][R[t][j]]];
                    ns++;
                }
        }
        for (int e = 0; e < ns; e++) {
            int x = vindex(span[e]);
            bs[x >> 6] |= 1ULL << (x & 63);
            wk->elems[nel++] = x;
        }
        act[nact++] = m;
    }
    /* the admissible a: deg_H(6) = |supp a| >= deg_H(j) = deg[j] + [a_j != 0] for every j (relabelling reduction) */
    int vm[7];
    for (int z = 0; z <= 6; z++) { vm[z] = 0; for (int j = 0; j < 6; j++) if (deg[j] >= z) vm[z] |= 1 << j; }
    int dmax = 0;
    for (int j = 0; j < 6; j++) if (deg[j] > dmax) dmax = deg[j];
    wk->gen++; wk->nl = 0;
    for (int i = 0; i < NA; i++) {
        if (!all) {
            if (AZ[i] < dmax) break;            /* candidates are sorted by support size */
            if (AS[i] & vm[AZ[i]]) continue;
        }
        const uint32_t *pr = PROJ + (size_t)i * 64;
        uint64_t mask = 0;
        for (int t = 0; t < nact; t++) {
            int m = act[t];
            uint32_t x = pr[m];
            if (!((wk->bits[(size_t)m * wk->W + (x >> 6)] >> (x & 63)) & 1)) mask |= 1ULL << (m - 1);
        }
        wk->graphs++;
        uint64_t h = (mask * 0x9E3779B97F4A7C15ULL) >> 52;     /* 12 bits */
        while (wk->lgen[h] == wk->gen && wk->lkey[h] != mask) h = (h + 1) & (LCAP - 1);
        if (wk->lgen[h] != wk->gen) {
            if (wk->nl >= LCAP / 2) { fprintf(stderr, "local table full\n"); exit(3); }
            wk->lgen[h] = wk->gen; wk->lkey[h] = mask; wk->lfirst[h] = i; wk->lslot[wk->nl++] = (int)h;
        }
    }
    for (int s = 0; s < wk->nl; s++) {
        int h = wk->lslot[s];
        uint64_t mask = wk->lkey[h], k0 = base0, k1 = base1;
        for (int m = 1; m < 64; m++) if (mask >> (m - 1) & 1) {
            if (m - 1 < 32) k0 += 1ULL << (2 * (m - 1)); else k1 += 1ULL << (2 * (m - 33));
        }
        uint8_t w[21];
        const uint8_t *a = AV + (size_t)wk->lfirst[h] * 6;
        int k = 0;
        for (int u = 0; u < 7; u++) for (int v = u + 1; v < 7; v++) w[k++] = v == 6 ? a[u] : A[u][v];
        tadd(&wk->T, k0, k1, w);
    }
    from[nact] = nel;
    for (int t = 0; t < nact; t++) {    /* clear the bitsets */
        uint64_t *bs = wk->bits + (size_t)act[t] * wk->W;
        for (int e = from[t]; e < from[t + 1]; e++) bs[wk->elems[e] >> 6] = 0;
    }
    wk->prefixes++;
}

static int popc(const uint8_t *w, int k) { int c = 0; for (int j = 0; j < k; j++) c += w[j] != 0; return c; }

int main(int argc, char **argv) {
    if (argc < 2 || !setup(atoi(argv[1]))) {
        fprintf(stderr, "usage: s40_enum7q q [threads] [all|reduced] [slice s n]   (q in 2 3 4 5 7 8 9)\n");
        return 1;
    }
    int threads = argc > 2 ? atoi(argv[2]) : 1;
    int all = 0, sl = 0, nsl = 1;
    for (int i = 3; i < argc; i++) {
        if (!strcmp(argv[i], "all")) all = 1;
        else if (!strcmp(argv[i], "reduced")) all = 0;
        else if (!strcmp(argv[i], "slice") && i + 2 < argc) { sl = atoi(argv[i + 1]); nsl = atoi(argv[i + 2]); i += 2; }
        else { fprintf(stderr, "unknown argument %s\n", argv[i]); return 1; }
    }
    if (nsl < 1 || sl < 0 || sl >= nsl) { fprintf(stderr, "bad slice\n"); return 1; }
    normal_forms();
#ifdef _OPENMP
    if (threads > 0) omp_set_num_threads(threads);
#endif
    int nth = 1;
#ifdef _OPENMP
    nth = omp_get_max_threads();
#endif
    work *wks = calloc(nth, sizeof(work));
    for (int t = 0; t < nth; t++) {
        wks[t].W = (Q6 + 63) / 64;
        wks[t].bits = calloc((size_t)64 * wks[t].W, sizeof(uint64_t));
        wks[t].elems = malloc(sizeof(int) * 64 * 729);
        tinit(&wks[t].T, 1u << 16);
    }
    #pragma omp parallel for schedule(dynamic, 1)
    for (int i5 = 0; i5 < NNF[5]; i5++) {
        if (i5 % nsl != sl) continue;
        int tid = 0;
#ifdef _OPENMP
        tid = omp_get_thread_num();
#endif
        work *wk = &wks[tid];
        uint8_t A[6][6];
        memset(A, 0, sizeof A);
        const uint8_t *w5 = NFW[5] + (size_t)i5 * 6;
        for (int j = 0; j < 5; j++) A[5][j] = A[j][5] = w5[j];
        int deg[6];
        deg[5] = popc(w5, 5);
        for (int i4 = 0; i4 < NNF[4]; i4++) {
            const uint8_t *w4 = NFW[4] + (size_t)i4 * 6;
            deg[4] = popc(w4, 4) + (A[5][4] != 0);
            if (!all && deg[4] < deg[5]) continue;
            for (int j = 0; j < 4; j++) A[4][j] = A[j][4] = w4[j];
            for (int i3 = 0; i3 < NNF[3]; i3++) {
                const uint8_t *w3 = NFW[3] + (size_t)i3 * 6;
                deg[3] = popc(w3, 3) + (A[5][3] != 0) + (A[4][3] != 0);
                if (!all && deg[3] < deg[4]) continue;
                for (int j = 0; j < 3; j++) A[3][j] = A[j][3] = w3[j];
                for (int i2 = 0; i2 < NNF[2]; i2++) {
                    const uint8_t *w2 = NFW[2] + (size_t)i2 * 6;
                    deg[2] = popc(w2, 2) + (A[5][2] != 0) + (A[4][2] != 0) + (A[3][2] != 0);
                    if (!all && deg[2] < deg[3]) continue;
                    for (int j = 0; j < 2; j++) A[2][j] = A[j][2] = w2[j];
                    for (int i1 = 0; i1 < NNF[1]; i1++) {
                        const uint8_t *w1 = NFW[1] + (size_t)i1 * 6;
                        deg[1] = popc(w1, 1) + (A[5][1] != 0) + (A[4][1] != 0) + (A[3][1] != 0) + (A[2][1] != 0);
                        if (!all && deg[1] < deg[2]) continue;
                        A[1][0] = A[0][1] = w1[0];
                        deg[0] = (A[1][0] != 0) + (A[2][0] != 0) + (A[3][0] != 0) + (A[4][0] != 0) + (A[5][0] != 0);
                        if (!all && deg[0] < deg[1]) continue;
                        process(wk, A, deg, all);
                    }
                }
            }
        }
    }
    table F; tinit(&F, 1u << 17);
    uint64_t prefixes = 0, graphs = 0;
    for (int t = 0; t < nth; t++) {
        prefixes += wks[t].prefixes; graphs += wks[t].graphs;
        for (uint64_t h = 0; h < wks[t].T.cap; h++) if (wks[t].T.t[h].used)
            tadd(&F, wks[t].T.t[h].k0, wks[t].T.t[h].k1, wks[t].T.t[h].w);
    }
    for (uint64_t h = 0; h < F.cap; h++) if (F.t[h].used) {
        printf("V");
        for (int m = 1; m < 64; m++) {
            uint64_t k = m - 1 < 32 ? F.t[h].k0 >> (2 * (m - 1)) : F.t[h].k1 >> (2 * (m - 33));
            printf(" %d", (int)(k & 3));
        }
        for (int e = 0; e < 21; e++) printf(" %d", F.t[h].w[e]);
        printf("\n");
    }
    printf("TOTAL %d %llu %llu %llu\n", Q, (unsigned long long)prefixes, (unsigned long long)graphs, (unsigned long long)F.n);
    return 0;
}
