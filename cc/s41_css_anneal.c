/* s41: CSS search for five-party entropy rays (R52).
 *
 * Target: the rank function h of six column groups A..F (bit i of a mask = group i, F = the purifier last),
 * given on the 63 nonempty masks.  For a pure-symmetric five-party vector S and its lift
 * h_S(X) = S(X) + sum_{i in X} S(i) (Theorem 1 of the paper; R32 note, proof of Theorem 1, step 2 and the remark
 * "symplectic structure is irrelevant"), any N x m binary matrix M with N = h(E) rows and column groups of widths
 * h(i), whose column-group ranks are r(X) = h(X) for all X, gives the CSS state of the row space of M, with entropy
 * vector S_C(X) = r(X) + r(E \ X) - r(E) = 2 S(X): a qubit realisation of S at lambda = 2 (party i owns
 * h(i) = 2 S(i) qubits).  Every qubit realisation at lambda = 1 gives such an M ([I | W] for a graph state W), so
 * the search space has the size of the lambda = 1 graph-state search and at least as many solutions.  For a
 * multiple kappa * h the same construction gives lambda = 2 kappa.
 *
 * Search: simulated annealing on the columns of groups B..F; group A is fixed to the first h(A) unit vectors
 * (row operations do not change r).  Energy E = sum_X |r(X) - h(X)| over the 63 masks.  Moves (within one group i
 * at a time, in bursts): with probability pxor add to one column of group i one column (or the sum of two columns)
 * of other groups, which moves the subspace of group i towards theirs; otherwise flip one bit.  A change of
 * group i only changes the 32 masks containing i; for each subset Y of the other five groups the kernel keeps an
 * echelon basis of their span and the residuals of the columns of group i modulo it, so
 * r({i} u Y) = r(Y) + rank(residuals).  Every E = 0 configuration is printed; the Python wrapper
 * (scripts/s41_css_search.py) re-verifies it with independent code before it is accepted.
 *
 * Input (stdin), one target per line:   id seconds h_1 h_2 ... h_63      (h_X for the masks X = 1..63)
 * Arguments: seed [T0 T1 sweeps_per_restart pxor [dump_best]]
 * Output (stdout), one line per target:
 *   id FOUND seconds restarts moves m col_0 ... col_{m-1}       (columns as hex words; bit k = row k)
 *   id BEST E restarts moves seconds                            (dump_best = 0, the default)
 *   id BEST E restarts moves seconds m col_0 ... col_{m-1}      (dump_best = 1: the best configuration seen)
 *   id SKIP reason
 * Limits: N <= 64 rows, m <= 128 columns, every group of width <= 64.
 * Build: gcc -O3 -march=native -fopenmp -o s41_css_anneal s41_css_anneal.c -lm
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>
#ifdef _OPENMP
#include <omp.h>
#endif

typedef uint64_t u64;
typedef struct { int id; double seconds; int h[64]; } Job;

static inline u64 rng_next(u64 *s) { u64 x = *s; x ^= x >> 12; x ^= x << 25; x ^= x >> 27; *s = x; return x * 0x2545F4914F6CDD1DULL; }
static inline double rng_unit(u64 *s) { return (rng_next(s) >> 11) * (1.0 / 9007199254740992.0); }
static double now_s(void) { struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts); return ts.tv_sec + ts.tv_nsec * 1e-9; }

/* rank of a set of vectors (xor basis keyed by the highest bit) */
static inline int rank_of(const u64 *v, int n) {
    u64 b[64]; u64 have = 0; int r = 0;
    for (int k = 0; k < n; k++) {
        u64 x = v[k];
        while (x) { int p = 63 - __builtin_clzll(x); if (have >> p & 1) x ^= b[p]; else { b[p] = x; have |= 1ULL << p; r++; break; } }
    }
    return r;
}

typedef struct {
    int N, m, nb[6], off[6], h[64], others[6][5];
    u64 col[128];
    int rk[64];
    /* cache for the moving group */
    int cur;
    int br[32]; int bpiv[32][64]; u64 bvec[32][64];
    u64 res[32][64];
    int Xof[32];
} State;

static int full_energy(State *S) {
    int E = 0; u64 tmp[128];
    for (int X = 1; X < 64; X++) {
        int n = 0;
        for (int j = 0; j < 6; j++) if (X >> j & 1) for (int c = 0; c < S->nb[j]; c++) tmp[n++] = S->col[S->off[j] + c];
        S->rk[X] = rank_of(tmp, n); E += abs(S->rk[X] - S->h[X]);
    }
    return E;
}

static inline u64 reduce(const State *S, int y, u64 v) {
    for (int k = 0; k < S->br[y]; k++) if (v >> S->bpiv[y][k] & 1) v ^= S->bvec[y][k];
    return v;
}

static void build_cache(State *S, int i) {
    S->cur = i;
    for (int y = 0; y < 32; y++) {
        u64 b[64]; u64 have = 0; int r = 0;
        int X = 1 << i;
        for (int t = 0; t < 5; t++) if (y >> t & 1) {
            int j = S->others[i][t]; X |= 1 << j;
            for (int c = 0; c < S->nb[j]; c++) {
                u64 x = S->col[S->off[j] + c];
                while (x) { int p = 63 - __builtin_clzll(x); if (have >> p & 1) x ^= b[p]; else { b[p] = x; have |= 1ULL << p; r++; break; } }
            }
        }
        S->Xof[y] = X; S->br[y] = 0;
        for (int p = 63; p >= 0; p--) if (have >> p & 1) { S->bpiv[y][S->br[y]] = p; S->bvec[y][S->br[y]] = b[p]; S->br[y]++; }
        for (int c = 0; c < S->nb[i]; c++) S->res[y][c] = reduce(S, y, S->col[S->off[i] + c]);
    }
}

static void randomise(State *S, u64 *rs) {
    u64 mask = (S->N == 64) ? ~0ULL : ((1ULL << S->N) - 1);
    for (int c = 0; c < S->nb[0]; c++) S->col[c] = 1ULL << c;
    for (int j = 1; j < 6; j++) for (int c = 0; c < S->nb[j]; c++) { u64 v; do v = rng_next(rs) & mask; while (!v); S->col[S->off[j] + c] = v; }
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: %s seed [T0 T1 sweeps pxor [dump_best]] < targets\n", argv[0]); return 2; }
    u64 seed0 = strtoull(argv[1], 0, 10);
    double T0 = argc > 2 ? atof(argv[2]) : 1.5, T1 = argc > 3 ? atof(argv[3]) : 0.08;
    double sweeps = argc > 4 ? atof(argv[4]) : 1000.0, pxor = argc > 5 ? atof(argv[5]) : 0.5;
    int dump_best = argc > 6 ? atoi(argv[6]) : 0;

    int cap = 1024, nj = 0; Job *jobs = malloc(cap * sizeof(Job));
    while (1) {
        Job J; if (scanf("%d %lf", &J.id, &J.seconds) != 2) break;
        J.h[0] = 0;
        for (int X = 1; X < 64; X++) if (scanf("%d", &J.h[X]) != 1) return 3;
        if (nj == cap) { cap *= 2; jobs = realloc(jobs, cap * sizeof(Job)); }
        jobs[nj++] = J;
    }

    #pragma omp parallel for schedule(dynamic, 1)
    for (int jj = 0; jj < nj; jj++) {
        Job *J = &jobs[jj];
        State *S = malloc(sizeof(State));
        char *line = malloc(1 << 13); int len = 0;
        memcpy(S->h, J->h, sizeof S->h);
        S->N = S->h[63]; S->m = 0;
        int bad = 0;
        for (int i = 0; i < 6; i++) { S->nb[i] = S->h[1 << i]; S->off[i] = S->m; S->m += S->nb[i]; if (S->nb[i] > 64 || S->nb[i] < 0) bad = 1; }
        if (S->N > 64 || S->N <= 0 || S->m > 128 || bad || S->nb[0] > S->N) {
            snprintf(line, 1 << 13, "%d SKIP N=%d m=%d\n", J->id, S->N, S->m);
            #pragma omp critical
            { fputs(line, stdout); fflush(stdout); }
            free(S); free(line); continue;
        }
        for (int i = 0; i < 6; i++) { int t = 0; for (int j = 0; j < 6; j++) if (j != i) S->others[i][t++] = j; }
        u64 rs = seed0 ^ (0x9E3779B97F4A7C15ULL * (u64)(J->id + 1)) ^ ((u64)jj << 48);
        for (int w = 0; w < 8; w++) rng_next(&rs);
        int freebits = 0; for (int j = 1; j < 6; j++) freebits += S->nb[j] * S->N;
        long L = (long)(sweeps * (freebits > 0 ? freebits : 1));
        double t0 = now_s(), secs = 0; long restarts = 0, moves = 0; int bestE = 1 << 30, found = 0;
        u64 bestcol[128]; memset(bestcol, 0, sizeof bestcol);
        u64 mask = (S->N == 64) ? ~0ULL : ((1ULL << S->N) - 1);
        while (!found) {
            randomise(S, &rs); int E = full_energy(S); restarts++; S->cur = -1;
            if (E < bestE) { bestE = E; memcpy(bestcol, S->col, sizeof bestcol); }
            if (E == 0) { found = 1; break; }
            double lr = log(T1 / T0); long s = 0; int timeout = 0;
            while (s < L) {
                int i = 1 + (int)(rng_next(&rs) % 5);
                if (S->nb[i] == 0) { s++; continue; }
                if (i != S->cur) build_cache(S, i);
                double T = T0 * exp(lr * (double)s / (double)L);
                for (int k = 0; k < 64 && s < L; k++, s++) {
                    int c = (int)(rng_next(&rs) % S->nb[i]); u64 v = S->col[S->off[i] + c];
                    if (rng_unit(&rs) < pxor) {
                        int j = S->others[i][rng_next(&rs) % 5]; if (!S->nb[j]) continue;
                        u64 a = S->col[S->off[j] + (int)(rng_next(&rs) % S->nb[j])];
                        if (rng_unit(&rs) < 0.5) { int j2 = S->others[i][rng_next(&rs) % 5]; if (S->nb[j2]) a ^= S->col[S->off[j2] + (int)(rng_next(&rs) % S->nb[j2])]; }
                        v ^= a;
                    } else v ^= 1ULL << (rng_next(&rs) % S->N);
                    v &= mask;
                    if (!v) continue;
                    moves++;
                    int d = 0, nr[32]; u64 rv[32];
                    for (int y = 0; y < 32; y++) {
                        rv[y] = reduce(S, y, v);
                        u64 tmp[64]; int n = 0;
                        for (int c2 = 0; c2 < S->nb[i]; c2++) tmp[n++] = (c2 == c) ? rv[y] : S->res[y][c2];
                        nr[y] = S->br[y] + rank_of(tmp, n);
                        int X = S->Xof[y];
                        d += abs(nr[y] - S->h[X]) - abs(S->rk[X] - S->h[X]);
                    }
                    if (d <= 0 || rng_unit(&rs) < exp(-d / T)) {
                        S->col[S->off[i] + c] = v;
                        for (int y = 0; y < 32; y++) { S->res[y][c] = rv[y]; S->rk[S->Xof[y]] = nr[y]; }
                        E += d;
                        if (E < bestE) { bestE = E; memcpy(bestcol, S->col, sizeof bestcol); }
                        if (E == 0) { found = 1; break; }
                    }
                }
                if (found) break;
                secs = now_s() - t0;
                if (secs > J->seconds) { timeout = 1; break; }
            }
            if (found || timeout) break;
        }
        secs = now_s() - t0;
        if (found) len = snprintf(line, 1 << 13, "%d FOUND %.3f %ld %ld %d", J->id, secs, restarts, moves, S->m);
        else len = snprintf(line, 1 << 13, "%d BEST %d %ld %ld %.3f", J->id, bestE, restarts, moves, secs);
        if (found || dump_best) {
            u64 *out = found ? S->col : bestcol;
            if (!found) len += snprintf(line + len, (1 << 13) - len, " %d", S->m);
            for (int q = 0; q < S->m; q++) len += snprintf(line + len, (1 << 13) - len, " %llx", (unsigned long long)out[q]);
        }
        snprintf(line + len, (1 << 13) - len, "\n");
        #pragma omp critical
        { fputs(line, stdout); fflush(stdout); }
        free(S); free(line);
    }
    free(jobs);
    return 0;
}
