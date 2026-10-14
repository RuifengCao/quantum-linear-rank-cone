/* s35: all entropy vectors of weighted graph states over GF(p) (p = 2, 3, 5) on seven vertices, one qudit per system
 * (systems 0..5 visible, 6 the purifier).  Every stabilizer state of seven qudits over Z_p is local-Clifford
 * equivalent to such a graph state, and entropies are local-unitary invariant, so the distinct cut-rank vectors
 * (63 cuts: the nonempty subsets of the visible systems) are exactly the entropy vectors of seven-qudit stabilizer
 * states over Z_p with at most one qudit per system.
 *
 * Symmetry used: scaling a vertex by a unit u (row and column times u) keeps every cut rank.  Normal form: for every
 * vertex v with a neighbour u < v, the edge to its smallest such neighbour has weight 1.  Every scaling orbit has at
 * least one graph in normal form (fix the weights along these edges vertex by vertex), so enumerating the normal forms
 * is complete; a few orbits are visited more than once, which only costs time.
 *
 * Output: one line per distinct vector: "V count t1 ... t63 w01 w02 ... w56" (the weights of the first graph found,
 * edge order (0,1), (0,2), ..., (5,6)), then "TOTAL graphs distinct".
 * Arguments: p [threads]
 * Build: gcc -O3 -march=native -fopenmp -o s35_enum7 s35_enum7.c
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#define NV 7
#define NE 21
static int P;
static int EU[NE], EV[NE];
static int cut_rows[64][7], cut_nr[64], cut_cols[64][7], cut_nc[64];

/* rank over GF(P) of the a x b matrix M (a <= 3 after transposition), entries already reduced mod P */
static inline int rank_small(int M[3][7], int a, int b) {
    int r = 0;
    for (int c = 0; c < b && r < a; c++) {
        int piv = -1;
        for (int i = r; i < a; i++) if (M[i][c]) { piv = i; break; }
        if (piv < 0) continue;
        if (piv != r) for (int k = 0; k < b; k++) { int t = M[r][k]; M[r][k] = M[piv][k]; M[piv][k] = t; }
        int inv = 1;
        for (int x = 1; x < P; x++) if ((M[r][c] * x) % P == 1) { inv = x; break; }
        for (int i = r + 1; i < a; i++) if (M[i][c]) {
            int f = (P - (M[i][c] * inv) % P) % P;
            for (int k = c; k < b; k++) M[i][k] = (M[i][k] + f * M[r][k]) % P;
        }
        r++;
    }
    return r;
}

typedef struct { uint64_t k0, k1; uint64_t count; int8_t w[NE]; int used; } Entry;
typedef struct { Entry *tab; uint64_t cap, n; } Table;

static void tab_init(Table *T, uint64_t cap) { T->cap = cap; T->n = 0; T->tab = calloc(cap, sizeof(Entry)); }
static void tab_add(Table *T, uint64_t k0, uint64_t k1, uint64_t cnt, const int8_t *w);
static void tab_grow(Table *T) {
    Table U; tab_init(&U, T->cap * 2);
    for (uint64_t i = 0; i < T->cap; i++) if (T->tab[i].used) tab_add(&U, T->tab[i].k0, T->tab[i].k1, T->tab[i].count, T->tab[i].w);
    free(T->tab); *T = U;
}
static void tab_add(Table *T, uint64_t k0, uint64_t k1, uint64_t cnt, const int8_t *w) {
    if (2 * (T->n + 1) > T->cap) tab_grow(T);
    uint64_t h = (k0 * 0x9E3779B97F4A7C15ULL) ^ (k1 * 0xC2B2AE3D27D4EB4FULL);
    uint64_t i = h & (T->cap - 1);
    while (T->tab[i].used && !(T->tab[i].k0 == k0 && T->tab[i].k1 == k1)) i = (i + 1) & (T->cap - 1);
    if (!T->tab[i].used) { T->tab[i].used = 1; T->tab[i].k0 = k0; T->tab[i].k1 = k1; T->tab[i].count = 0; memcpy(T->tab[i].w, w, NE); T->n++; }
    T->tab[i].count += cnt;
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: %s p\n", argv[0]); return 2; }
    P = atoi(argv[1]);
    if (P != 2 && P != 3 && P != 5) { fprintf(stderr, "p must be 2, 3 or 5\n"); return 2; }
    int e = 0;
    for (int u = 0; u < NV; u++) for (int v = u + 1; v < NV; v++) { EU[e] = u; EV[e] = v; e++; }
    for (int m = 1; m < 64; m++) {
        cut_nr[m] = cut_nc[m] = 0;
        for (int v = 0; v < NV; v++) {
            if (v < 6 && (m >> v & 1)) cut_rows[m][cut_nr[m]++] = v; else cut_cols[m][cut_nc[m]++] = v;
        }
    }
    /* per vertex v >= 1: the admissible weight vectors (w_0v .. w_{v-1}v) in normal form */
    int nopt[NV]; int (*opt[NV])[NV];
    for (int v = 1; v < NV; v++) {
        int total = 1; for (int k = 0; k < v; k++) total *= P;
        opt[v] = malloc(sizeof(int[NV]) * total); nopt[v] = 0;
        for (int code = 0; code < total; code++) {
            int d[NV], c = code, first = -1;
            for (int k = 0; k < v; k++) { d[k] = c % P; c /= P; }
            for (int k = 0; k < v; k++) if (d[k]) { first = k; break; }
            if (first >= 0 && d[first] != 1) continue;
            memcpy(opt[v][nopt[v]++], d, sizeof(int) * v);
        }
    }
    uint64_t ngraphs = 1;
    for (int v = 1; v < NV; v++) ngraphs *= (uint64_t)nopt[v];
    fprintf(stderr, "p = %d: %llu graphs in normal form\n", P, (unsigned long long)ngraphs);
    int nth = 1;
#ifdef _OPENMP
    nth = omp_get_max_threads();
#endif
    Table *tabs = malloc(sizeof(Table) * nth);
    for (int t = 0; t < nth; t++) tab_init(&tabs[t], 1 << 16);
    /* split over the options of vertices 6 and 5 (outer loops) */
    long outer = (long)nopt[6] * nopt[5];
    #pragma omp parallel for schedule(dynamic, 1)
    for (long o = 0; o < outer; o++) {
        int tid = 0;
#ifdef _OPENMP
        tid = omp_get_thread_num();
#endif
        int A[NV][NV]; memset(A, 0, sizeof A);
        int i6 = (int)(o / nopt[5]), i5 = (int)(o % nopt[5]);
        for (int k = 0; k < 6; k++) A[k][6] = A[6][k] = opt[6][i6][k];
        for (int k = 0; k < 5; k++) A[k][5] = A[5][k] = opt[5][i5][k];
        for (int i4 = 0; i4 < nopt[4]; i4++) {
            for (int k = 0; k < 4; k++) A[k][4] = A[4][k] = opt[4][i4][k];
            for (int i3 = 0; i3 < nopt[3]; i3++) {
                for (int k = 0; k < 3; k++) A[k][3] = A[3][k] = opt[3][i3][k];
                for (int i2 = 0; i2 < nopt[2]; i2++) {
                    for (int k = 0; k < 2; k++) A[k][2] = A[2][k] = opt[2][i2][k];
                    for (int i1 = 0; i1 < nopt[1]; i1++) {
                        A[0][1] = A[1][0] = opt[1][i1][0];
                        uint64_t k0 = 0, k1 = 0;
                        for (int m = 1; m < 64; m++) {
                            int M[3][7], a, b;
                            if (cut_nr[m] <= cut_nc[m]) {
                                a = cut_nr[m]; b = cut_nc[m];
                                for (int i = 0; i < a; i++) for (int j = 0; j < b; j++) M[i][j] = A[cut_rows[m][i]][cut_cols[m][j]];
                            } else {
                                a = cut_nc[m]; b = cut_nr[m];
                                for (int i = 0; i < a; i++) for (int j = 0; j < b; j++) M[i][j] = A[cut_cols[m][i]][cut_rows[m][j]];
                            }
                            uint64_t r = (uint64_t)rank_small(M, a, b);
                            if (m < 32) k0 |= r << (2 * m); else k1 |= r << (2 * (m - 32));
                        }
                        int8_t w[NE];
                        for (int q = 0; q < NE; q++) w[q] = (int8_t)A[EU[q]][EV[q]];
                        tab_add(&tabs[tid], k0, k1, 1, w);
                    }
                }
            }
        }
    }
    Table all; tab_init(&all, 1 << 16);
    for (int t = 0; t < nth; t++) {
        for (uint64_t i = 0; i < tabs[t].cap; i++) if (tabs[t].tab[i].used)
            tab_add(&all, tabs[t].tab[i].k0, tabs[t].tab[i].k1, tabs[t].tab[i].count, tabs[t].tab[i].w);
    }
    uint64_t tot = 0;
    for (uint64_t i = 0; i < all.cap; i++) if (all.tab[i].used) {
        Entry *E = &all.tab[i]; tot += E->count;
        printf("V %llu", (unsigned long long)E->count);
        for (int m = 1; m < 64; m++) {
            int r = m < 32 ? (int)(E->k0 >> (2 * m) & 3) : (int)(E->k1 >> (2 * (m - 32)) & 3);
            printf(" %d", r);
        }
        for (int q = 0; q < NE; q++) printf(" %d", E->w[q]);
        printf("\n");
    }
    printf("TOTAL %llu %llu\n", (unsigned long long)tot, (unsigned long long)all.n);
    return 0;
}
