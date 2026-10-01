/* s16: targeted stabilizer (graph-state) search for five-party entropy rays.
 *
 * Target: an integer vector t over the 31 nonempty subsets X of the visible
 * parties A..E (mask order, bit0 = A), scaled by lambda; the purifier F is the
 * sixth party.  We search for a qubit graph state on N = sum_p n_p qubits,
 * party p owning n_p = lambda * S_p consecutive qubits, whose cut-rank entropies
 * rank_GF2 Gamma[X, X^c] equal t_X for all 31 X.
 *
 * Why this search space is complete for a fixed lambda (see README / report):
 *  - every stabilizer state is local-Clifford equivalent to a graph state and
 *    entropies are LU-invariant;
 *  - for a stabilizer state, the marginal of each party p is maximally mixed on
 *    a 2^{S_p}-dimensional stabilizer subspace, so a party-local Clifford moves
 *    |p| - S_p of its qubits into a product |0> state that can be discarded:
 *    WLOG party p has exactly S_p qubits;
 *  - within-party edges never enter any cut Gamma[X, X^c] (X is a union of
 *    whole parties), so only inter-party bits are free.
 *
 * Method: simulated annealing on the inter-party bits, energy
 * E = sum_X |rank_X - t_X|; a flip of (i,j) only changes the 16 masks that
 * separate party(i) from party(j).  Any E = 0 configuration is printed with its
 * adjacency rows; the Python wrapper re-verifies it with an independent
 * implementation before it is accepted as a certificate.
 *
 * Input (stdin), one target per line:
 *   id lambda n0 n1 n2 n3 n4 n5 t1 t2 ... t31
 * Arguments: seconds_per_target seed [T0 T1 sweeps_per_restart [dump_best]]
 * Output (stdout), one line per target:
 *   id lambda FOUND N hexrow0 ... hexrow(N-1)
 *   id lambda BEST E restarts steps            (dump_best = 0, the default)
 *   id lambda BEST E restarts steps N hexrow0 ... hexrow(N-1)   (dump_best = 1:
 *       the best configuration seen, for residual analysis; R29)
 * Build: gcc -O3 -march=native -fopenmp -o s16_anneal s16_anneal.c -lm
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

typedef struct { int id, lam, n[6], t[32]; } Job;

static inline uint64_t rng_next(uint64_t *s) {
    uint64_t x = *s; x ^= x >> 12; x ^= x << 25; x ^= x >> 27; *s = x;
    return x * 0x2545F4914F6CDD1DULL;
}
static inline double rng_unit(uint64_t *s) { return (rng_next(s) >> 11) * (1.0 / 9007199254740992.0); }

static inline int rank_rows(const uint64_t *adj, uint64_t rows, uint64_t cmask) {
    uint64_t basis[64]; uint64_t have = 0; int r = 0;
    while (rows) {
        int i = __builtin_ctzll(rows); rows &= rows - 1;
        uint64_t v = adj[i] & cmask;
        while (v) {
            int p = 63 - __builtin_clzll(v);
            if (have >> p & 1) v ^= basis[p];
            else { basis[p] = v; have |= 1ULL << p; r++; break; }
        }
    }
    return r;
}

static double now_s(void) {
    struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: %s seconds_per_target seed [T0 T1 sweeps]\n", argv[0]); return 2; }
    double budget = atof(argv[1]);
    uint64_t seed0 = strtoull(argv[2], 0, 10);
    double T0 = argc > 3 ? atof(argv[3]) : 1.5, T1 = argc > 4 ? atof(argv[4]) : 0.08;
    double sweeps = argc > 5 ? atof(argv[5]) : 3000.0;
    int dump_best = argc > 6 ? atoi(argv[6]) : 0;

    int cap = 1024, nj = 0; Job *jobs = malloc(cap * sizeof(Job));
    while (1) {
        Job J; if (scanf("%d %d", &J.id, &J.lam) != 2) break;
        for (int p = 0; p < 6; p++) if (scanf("%d", &J.n[p]) != 1) return 3;
        J.t[0] = 0;
        for (int m = 1; m < 32; m++) if (scanf("%d", &J.t[m]) != 1) return 3;
        if (nj == cap) { cap *= 2; jobs = realloc(jobs, cap * sizeof(Job)); }
        jobs[nj++] = J;
    }

    /* masks separating parties p<q (F = party 5 is never inside X) */
    int sep[6][6][32], nsep[6][6];
    for (int p = 0; p < 6; p++) for (int q = 0; q < 6; q++) {
        nsep[p][q] = 0;
        for (int m = 1; m < 32; m++) {
            int ip = (p < 5) ? (m >> p & 1) : 0, iq = (q < 5) ? (m >> q & 1) : 0;
            if (ip != iq) sep[p][q][nsep[p][q]++] = m;
        }
    }

    #pragma omp parallel for schedule(dynamic, 1)
    for (int jj = 0; jj < nj; jj++) {
        Job *J = &jobs[jj];
        int N = 0, off[7], party[64];
        for (int p = 0; p < 6; p++) { off[p] = N; for (int k = 0; k < J->n[p]; k++) party[N + k] = p; N += J->n[p]; }
        off[6] = N;
        char line[8192]; int len = 0;
        if (N > 64 || N == 0) {
            snprintf(line, sizeof line, "%d %d SKIP N=%d\n", J->id, J->lam, N);
            #pragma omp critical
            { fputs(line, stdout); fflush(stdout); }
            continue;
        }
        uint64_t all = (N == 64) ? ~0ULL : ((1ULL << N) - 1), pm[6], qm[32], cm[32];
        for (int p = 0; p < 6; p++) pm[p] = (J->n[p] == 64) ? ~0ULL : (((1ULL << J->n[p]) - 1) << off[p]);
        for (int m = 0; m < 32; m++) { qm[m] = 0; for (int p = 0; p < 5; p++) if (m >> p & 1) qm[m] |= pm[p]; cm[m] = all & ~qm[m]; }
        int npairs = 0; int (*pairs)[2] = malloc(sizeof(int[2]) * N * N);
        for (int i = 0; i < N; i++) for (int j = i + 1; j < N; j++) if (party[i] != party[j]) { pairs[npairs][0] = i; pairs[npairs][1] = j; npairs++; }
        uint64_t rs = seed0 ^ (0x9E3779B97F4A7C15ULL * (uint64_t)(J->id + 1)) ^ ((uint64_t)J->lam << 40);
        for (int w = 0; w < 8; w++) rng_next(&rs);
        uint64_t adj[64], best_adj[64]; int rk[32], bestE = 1 << 30; long restarts = 0, steps = 0;
        double t_start = now_s(); int found = 0;
        long iters = (long)(sweeps * (npairs > 0 ? npairs : 1));
        while (!found && now_s() - t_start < budget) {
            restarts++;
            memset(adj, 0, sizeof adj);
            for (int k = 0; k < npairs; k++) if (rng_next(&rs) & 1) { int i = pairs[k][0], j = pairs[k][1]; adj[i] ^= 1ULL << j; adj[j] ^= 1ULL << i; }
            int E = 0;
            for (int m = 1; m < 32; m++) { rk[m] = rank_rows(adj, qm[m], cm[m]); E += abs(rk[m] - J->t[m]); }
            double lr = log(T1 / T0);
            for (long it = 0; it < iters && E > 0; it++) {
                if ((it & 0xFFFF) == 0 && now_s() - t_start >= budget) break;
                double T = T0 * exp(lr * (double)it / (double)iters);
                int k = (int)(rng_next(&rs) % (uint64_t)npairs), i = pairs[k][0], j = pairs[k][1];
                int p = party[i], q = party[j];
                adj[i] ^= 1ULL << j; adj[j] ^= 1ULL << i;
                int nr[16], dE = 0, ns = nsep[p][q];
                for (int s = 0; s < ns; s++) {
                    int m = sep[p][q][s];
                    nr[s] = rank_rows(adj, qm[m], cm[m]);
                    dE += abs(nr[s] - J->t[m]) - abs(rk[m] - J->t[m]);
                }
                if (dE <= 0 || rng_unit(&rs) < exp(-dE / T)) {
                    for (int s = 0; s < ns; s++) rk[sep[p][q][s]] = nr[s];
                    E += dE;
                    if (E < bestE) { bestE = E; memcpy(best_adj, adj, sizeof adj); }
                } else { adj[i] ^= 1ULL << j; adj[j] ^= 1ULL << i; }
                steps++;
            }
            if (E == 0) { found = 1; memcpy(best_adj, adj, sizeof adj); bestE = 0; }
        }
        free(pairs);
        if (found) {
            len = snprintf(line, sizeof line, "%d %d FOUND %d", J->id, J->lam, N);
            for (int i = 0; i < N; i++) len += snprintf(line + len, sizeof line - len, " %llx", (unsigned long long)best_adj[i]);
            snprintf(line + len, sizeof line - len, "\n");
        } else if (dump_best && bestE < (1 << 30)) {
            len = snprintf(line, sizeof line, "%d %d BEST %d %ld %ld %d", J->id, J->lam, bestE, restarts, steps, N);
            for (int i = 0; i < N; i++) len += snprintf(line + len, sizeof line - len, " %llx", (unsigned long long)best_adj[i]);
            snprintf(line + len, sizeof line - len, "\n");
        } else snprintf(line, sizeof line, "%d %d BEST %d %ld %ld\n", J->id, J->lam, bestE, restarts, steps);
        #pragma omp critical
        { fputs(line, stdout); fflush(stdout); }
    }
    free(jobs);
    return 0;
}
