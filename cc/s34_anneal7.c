/* s34: qubit graph-state search for six-party entropy vectors (seven systems: parties 0..5 visible, system 6 the
 * purifier).  Target: integers t_m for the 63 masks m of visible parties (bit i = party i), already multiplied by
 * lambda.  System p owns n_p consecutive qubits (n_p = lambda * S_p is complete for exact realisations: a party-local
 * Clifford discards the other qubits of a party).  Only inter-party edges enter the cuts W[X, X^c].
 *
 * Energy E = sum_m |rank_GF2 W[X_m, X_m^c] - t_m|.  Moves: flip one inter-party edge (only the 32 masks separating
 * the two parties change), and with probability p_neu a neutral move: local complementation at a vertex c (toggle
 * the edges between neighbours of c in different parties), which keeps every cut rank.  Restarts from random graphs
 * (or from the start graph given, with a random half of its edges re-drawn after the first restart).
 *
 * Input (stdin), one job per line:  id lambda n0..n6 t1..t63 S [N*N start bits]
 * Arguments: seconds seed [T0 T1 sweeps p_neu dump_best]
 * Output:    id lambda FOUND N hexrow_0 ... hexrow_{N-1}      (row i: bit j = edge ij)
 *            id lambda BEST E restarts steps seconds [N hexrows]   (rows with dump_best = 1)
 * The search proves nothing by itself: the Python wrapper re-verifies every reported graph.
 * Build: gcc -O3 -march=native -fopenmp -o s34_anneal7 s34_anneal7.c -lm
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

#define MAXN 64
#define NS 7
#define NM 64

typedef struct { int id, lam, n[NS], t[NM], has_start; uint64_t start[MAXN]; } Job;

static inline uint64_t rng_next(uint64_t *s) {
    uint64_t x = *s; x ^= x >> 12; x ^= x << 25; x ^= x >> 27; *s = x;
    return x * 0x2545F4914F6CDD1DULL;
}
static inline double rng_unit(uint64_t *s) { return (rng_next(s) >> 11) * (1.0 / 9007199254740992.0); }
static double now_s(void) { struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts); return ts.tv_sec + ts.tv_nsec * 1e-9; }

/* rank over GF(2) of the rows W[i] & cmask, i in rows (zero columns do not change the rank) */
static inline int rank_cut(const uint64_t *W, uint64_t rmask, uint64_t cmask) {
    uint64_t basis[64]; int r = 0;
    uint64_t lead_used = 0;
    while (rmask) {
        int i = __builtin_ctzll(rmask); rmask &= rmask - 1;
        uint64_t v = W[i] & cmask;
        while (v) {
            int lead = 63 - __builtin_clzll(v);
            if (!(lead_used >> lead & 1)) { basis[lead] = v; lead_used |= 1ULL << lead; r++; break; }
            v ^= basis[lead];
        }
    }
    return r;
}

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: %s seconds seed [T0 T1 sweeps p_neu dump_best]\n", argv[0]); return 2; }
    double budget = atof(argv[1]); uint64_t seed0 = strtoull(argv[2], 0, 10);
    double T0 = argc > 3 ? atof(argv[3]) : 1.0, T1 = argc > 4 ? atof(argv[4]) : 0.05;
    double sweeps = argc > 5 ? atof(argv[5]) : 100.0, p_neu = argc > 6 ? atof(argv[6]) : 0.2;
    int dump_best = argc > 7 ? atoi(argv[7]) : 0;
    int cap = 64, nj = 0; Job *jobs = malloc(cap * sizeof(Job));
    while (1) {
        Job J; memset(&J, 0, sizeof J);
        if (scanf("%d %d", &J.id, &J.lam) != 2) break;
        int N = 0;
        for (int p = 0; p < NS; p++) { if (scanf("%d", &J.n[p]) != 1) return 3; N += J.n[p]; }
        for (int m = 1; m < NM; m++) if (scanf("%d", &J.t[m]) != 1) return 3;
        if (scanf("%d", &J.has_start) != 1) return 3;
        if (J.has_start) for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) {
            int b; if (scanf("%d", &b) != 1) return 3;
            if (b & 1) J.start[i] |= 1ULL << j;
        }
        if (nj == cap) { cap *= 2; jobs = realloc(jobs, cap * sizeof(Job)); }
        jobs[nj++] = J;
    }
    #pragma omp parallel for schedule(dynamic, 1)
    for (int jj = 0; jj < nj; jj++) {
        Job *J = &jobs[jj];
        int N = 0, party[MAXN];
        for (int p = 0; p < NS; p++) for (int k = 0; k < J->n[p]; k++) { if (N < MAXN) party[N] = p; N++; }
        char *line = malloc(1 << 15); int len;
        if (N > MAXN || N == 0) {
            snprintf(line, 1 << 15, "%d %d SKIP N=%d\n", J->id, J->lam, N);
            #pragma omp critical
            { fputs(line, stdout); fflush(stdout); }
            free(line); continue;
        }
        uint64_t pmask[NS] = {0}, rm[NM], cm[NM], all = (N == 64) ? ~0ULL : ((1ULL << N) - 1);
        for (int i = 0; i < N; i++) pmask[party[i]] |= 1ULL << i;
        for (int m = 1; m < NM; m++) {
            rm[m] = 0;
            for (int p = 0; p < 6; p++) if (m >> p & 1) rm[m] |= pmask[p];
            cm[m] = all & ~rm[m];
        }
        int sep[NS][NS][NM], nsep[NS][NS];
        for (int p = 0; p < NS; p++) for (int q = 0; q < NS; q++) {
            nsep[p][q] = 0;
            for (int m = 1; m < NM; m++) {
                int ip = p < 6 ? (m >> p & 1) : 0, iq = q < 6 ? (m >> q & 1) : 0;
                if (ip != iq) sep[p][q][nsep[p][q]++] = m;
            }
        }
        int npairs = 0, (*pairs)[2] = malloc(sizeof(int[2]) * N * N);
        for (int i = 0; i < N; i++) for (int j = i + 1; j < N; j++) if (party[i] != party[j]) { pairs[npairs][0] = i; pairs[npairs][1] = j; npairs++; }
        uint64_t rs = seed0 ^ (0x9E3779B97F4A7C15ULL * (uint64_t)(J->id + 1)) ^ ((uint64_t)J->lam << 40);
        for (int w = 0; w < 8; w++) rng_next(&rs);
        uint64_t W[MAXN], Best[MAXN];
        int rk[NM], bestE = 1 << 30, found = 0; long restarts = 0, steps = 0;
        memset(Best, 0, sizeof Best);
        double t_start = now_s();
        long iters = (long)(sweeps * (npairs > 0 ? npairs : 1));
        while (!found && now_s() - t_start < budget) {
            memset(W, 0, sizeof W);
            for (int k = 0; k < npairs; k++) {
                int i = pairs[k][0], j = pairs[k][1], b;
                if (J->has_start) { b = (J->start[i] >> j) & 1; if (restarts > 0 && (rng_next(&rs) & 1)) b = rng_next(&rs) & 1; }
                else b = rng_next(&rs) & 1;
                if (b) { W[i] |= 1ULL << j; W[j] |= 1ULL << i; }
            }
            restarts++;
            int E = 0;
            for (int m = 1; m < NM; m++) { rk[m] = rank_cut(W, rm[m], cm[m]); E += abs(rk[m] - J->t[m]); }
            if (E < bestE) { bestE = E; memcpy(Best, W, sizeof W); }
            double lr = log(T1 / T0);
            for (long it = 0; it < iters && E > 0; it++) {
                if ((it & 0x3FF) == 0 && now_s() - t_start >= budget) break;
                double T = T0 * exp(lr * (double)it / (double)iters);
                if (p_neu > 0 && rng_unit(&rs) < p_neu) {          /* local complementation at c: neutral */
                    int c = (int)(rng_next(&rs) % (uint64_t)N);
                    uint64_t nb = W[c];
                    uint64_t tmp = nb;
                    while (tmp) {
                        int u = __builtin_ctzll(tmp); tmp &= tmp - 1;
                        uint64_t other = nb & ~pmask[party[u]] & ~((2ULL << u) - 1);   /* u' > u, other party */
                        W[u] ^= other;
                        uint64_t o2 = other;
                        while (o2) { int v = __builtin_ctzll(o2); o2 &= o2 - 1; W[v] ^= 1ULL << u; }
                    }
                    steps++;
                    continue;
                }
                int k = (int)(rng_next(&rs) % (uint64_t)npairs), u = pairs[k][0], v = pairs[k][1];
                W[u] ^= 1ULL << v; W[v] ^= 1ULL << u;
                int p = party[u], q = party[v], ns = nsep[p][q], nr[NM], dE = 0;
                for (int s = 0; s < ns; s++) {
                    int m = sep[p][q][s];
                    nr[s] = rank_cut(W, rm[m], cm[m]);
                    dE += abs(nr[s] - J->t[m]) - abs(rk[m] - J->t[m]);
                }
                if (dE <= 0 || rng_unit(&rs) < exp(-dE / T)) {
                    for (int s = 0; s < ns; s++) rk[sep[p][q][s]] = nr[s];
                    E += dE;
                    if (E < bestE) { bestE = E; memcpy(Best, W, sizeof W); }
                } else { W[u] ^= 1ULL << v; W[v] ^= 1ULL << u; }
                steps++;
            }
            if (E == 0) found = 1;
        }
        double secs = now_s() - t_start;
        uint64_t *OUT = found ? W : Best;
        if (found) len = snprintf(line, 1 << 15, "%d %d FOUND %d", J->id, J->lam, N);
        else len = snprintf(line, 1 << 15, "%d %d BEST %d %ld %ld %.3f", J->id, J->lam, bestE, restarts, steps, secs);
        if (found || dump_best) {
            if (!found) len += snprintf(line + len, (1 << 15) - len, " %d", N);
            for (int i = 0; i < N; i++) len += snprintf(line + len, (1 << 15) - len, " %llx", (unsigned long long)OUT[i]);
        }
        snprintf(line + len, (1 << 15) - len, "\n");
        #pragma omp critical
        { fputs(line, stdout); fflush(stdout); }
        free(line); free(pairs);
    }
    free(jobs);
    return 0;
}
