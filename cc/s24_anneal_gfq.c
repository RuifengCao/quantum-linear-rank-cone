/* s24: weighted graph-state annealing over a finite field GF(q), q in {2, 3, 4, 5, 7, 8, 9} (R31).
 *
 * Same search as cc/s21_anneal_gfp.c (party p owns n_p = lambda * S_p qudits, only inter-party
 * weights enter the 31 cuts, energy = sum |rank W[X, X^c] - target_X|), with the field given by
 * addition / multiplication tables, so that prime powers work too: GF(4) = GF(2)[a]/(a^2 + a + 1),
 * GF(8) = GF(2)[a]/(a^3 + a + 1), GF(9) = GF(3)[a]/(a^2 + 1); element sum_i c_i a^i is encoded as
 * sum_i c_i p^i.  A hit over GF(p^k) is turned into a GF(p) certificate at k * lambda by
 * scripts/s24_galois_search.py (field reduction); every certificate is re-verified there.
 *
 * Input (stdin), one target per line:  id lambda n0..n5 t1..t31
 * Arguments: q seconds_per_target seed [T0 T1 sweeps_per_restart [dump_best]]
 * Output: id lambda FOUND N row0 ... row(N-1)        rows as strings of N symbols 0..q-1
 *         id lambda BEST E restarts steps
 * Build:  gcc -O3 -march=native -fopenmp -o s24_anneal_gfq s24_anneal_gfq.c -lm
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
typedef struct { int id, lam, n[6], t[32]; } Job;

static inline uint64_t rng_next(uint64_t *s) {
    uint64_t x = *s; x ^= x >> 12; x ^= x << 25; x ^= x >> 27; *s = x;
    return x * 0x2545F4914F6CDD1DULL;
}
static inline double rng_unit(uint64_t *s) { return (rng_next(s) >> 11) * (1.0 / 9007199254740992.0); }
static double now_s(void) { struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts); return ts.tv_sec + ts.tv_nsec * 1e-9; }

static int P = 3, INV[16], ADDT[16][16], SUBT[16][16], MULT[16][16];
static void build_field(int q) {
    int p = 0, k = 0, mod = 0;              /* q = p^k, modulus polynomial coefficients (low to high, without leading 1) */
    if (q == 2 || q == 3 || q == 5 || q == 7) { p = q; k = 1; }
    else if (q == 4) { p = 2; k = 2; mod = 1 + 1 * 2; }      /* x^2 = x + 1      -> coeffs [1,1] */
    else if (q == 8) { p = 2; k = 3; mod = 1 + 1 * 2; }      /* x^3 = x + 1      -> coeffs [1,1,0] */
    else if (q == 9) { p = 3; k = 2; mod = 2; }              /* x^2 = -1 = 2     -> coeffs [2,0] */
    else { fprintf(stderr, "q must be 2, 3, 4, 5, 7, 8 or 9\n"); exit(2); }
    int modc[4] = {0, 0, 0, 0};
    if (k > 1) { int m = mod; for (int i = 0; i < k; i++) { modc[i] = m % p; m /= p; } }
    for (int a = 0; a < q; a++) for (int b = 0; b < q; b++) {
        int da[4], db[4], s[4], d[4]; int x = a, y = b;
        for (int i = 0; i < k; i++) { da[i] = x % p; x /= p; db[i] = y % p; y /= p; }
        for (int i = 0; i < k; i++) { s[i] = (da[i] + db[i]) % p; d[i] = (da[i] - db[i] + p) % p; }
        int prod[8] = {0};
        for (int i = 0; i < k; i++) for (int j = 0; j < k; j++) prod[i + j] = (prod[i + j] + da[i] * db[j]) % p;
        for (int deg = 2 * k - 2; deg >= k; deg--) {   /* reduce x^deg = x^(deg-k) * (sum modc[i] x^i) */
            int c = prod[deg]; prod[deg] = 0;
            for (int i = 0; i < k; i++) prod[deg - k + i] = (prod[deg - k + i] + c * modc[i]) % p;
        }
        int sa = 0, sd = 0, sm = 0, pw = 1;
        for (int i = 0; i < k; i++) { sa += s[i] * pw; sd += d[i] * pw; sm += prod[i] * pw; pw *= p; }
        ADDT[a][b] = sa; SUBT[a][b] = sd; MULT[a][b] = sm;
    }
    for (int a = 1; a < q; a++) for (int b = 1; b < q; b++) if (MULT[a][b] == 1) INV[a] = b;
}
/* rank over GF(q) of W[rows, cols] (P holds q) */
static int rankp(uint8_t W[MAXN][MAXN], const int *ri, int nr, const int *ci, int nc) {
    uint8_t M[MAXN][MAXN];
    for (int a = 0; a < nr; a++) for (int b = 0; b < nc; b++) M[a][b] = W[ri[a]][ci[b]];
    int r = 0;
    for (int c = 0; c < nc && r < nr; c++) {
        int piv = -1;
        for (int a = r; a < nr; a++) if (M[a][c]) { piv = a; break; }
        if (piv < 0) continue;
        if (piv != r) for (int b = 0; b < nc; b++) { uint8_t t = M[r][b]; M[r][b] = M[piv][b]; M[piv][b] = t; }
        if (M[r][c] != 1) { int iv = INV[M[r][c]]; for (int b = 0; b < nc; b++) M[r][b] = (uint8_t)MULT[M[r][b]][iv]; }
        for (int a = 0; a < nr; a++) if (a != r && M[a][c]) {
            int f = M[a][c];
            for (int b = 0; b < nc; b++) M[a][b] = (uint8_t)SUBT[M[a][b]][MULT[f][M[r][b]]];
        }
        r++;
    }
    return r;
}

int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage: %s q seconds seed [T0 T1 sweeps [dump_best]]\n", argv[0]); return 2; }
    P = atoi(argv[1]);
    build_field(P);
    argv++; argc--;
    double budget = atof(argv[1]); uint64_t seed0 = strtoull(argv[2], 0, 10);
    double T0 = argc > 3 ? atof(argv[3]) : 1.5, T1 = argc > 4 ? atof(argv[4]) : 0.08;
    double sweeps = argc > 5 ? atof(argv[5]) : 3000.0;
    int dump_best = argc > 6 ? atoi(argv[6]) : 0;
    int cap = 256, nj = 0; Job *jobs = malloc(cap * sizeof(Job));
    while (1) {
        Job J; if (scanf("%d %d", &J.id, &J.lam) != 2) break;
        for (int p = 0; p < 6; p++) if (scanf("%d", &J.n[p]) != 1) return 3;
        J.t[0] = 0;
        for (int m = 1; m < 32; m++) if (scanf("%d", &J.t[m]) != 1) return 3;
        if (nj == cap) { cap *= 2; jobs = realloc(jobs, cap * sizeof(Job)); }
        jobs[nj++] = J;
    }
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
        int N = 0, party[MAXN];
        for (int p = 0; p < 6; p++) for (int k = 0; k < J->n[p] && N < MAXN + 1; k++) { if (N < MAXN) party[N] = p; N++; }
        char line[8192]; int len = 0;
        if (N > MAXN || N == 0) {
            snprintf(line, sizeof line, "%d %d SKIP N=%d\n", J->id, J->lam, N);
            #pragma omp critical
            { fputs(line, stdout); fflush(stdout); }
            continue;
        }
        /* index lists of X and X^c for every mask */
        int xi[32][MAXN], xn[32], yi[32][MAXN], yn[32];
        for (int m = 1; m < 32; m++) {
            xn[m] = yn[m] = 0;
            for (int i = 0; i < N; i++) {
                int in = party[i] < 5 && ((m >> party[i]) & 1);
                if (in) xi[m][xn[m]++] = i; else yi[m][yn[m]++] = i;
            }
        }
        int npairs = 0, (*pairs)[2] = malloc(sizeof(int[2]) * N * N);
        for (int i = 0; i < N; i++) for (int j = i + 1; j < N; j++) if (party[i] != party[j]) { pairs[npairs][0] = i; pairs[npairs][1] = j; npairs++; }
        uint64_t rs = seed0 ^ (0x9E3779B97F4A7C15ULL * (uint64_t)(J->id + 1)) ^ ((uint64_t)J->lam << 40);
        for (int w = 0; w < 8; w++) rng_next(&rs);
        uint8_t W[MAXN][MAXN], B[MAXN][MAXN];
        int rk[32], bestE = 1 << 30, found = 0; long restarts = 0, steps = 0;
        double t_start = now_s();
        long iters = (long)(sweeps * (npairs > 0 ? npairs : 1));
        while (!found && now_s() - t_start < budget) {
            restarts++;
            memset(W, 0, sizeof W);
            for (int k = 0; k < npairs; k++) { int i = pairs[k][0], j = pairs[k][1]; uint8_t v = (uint8_t)(rng_next(&rs) % (uint64_t)P); W[i][j] = W[j][i] = v; }
            int E = 0;
            for (int m = 1; m < 32; m++) { rk[m] = rankp(W, xi[m], xn[m], yi[m], yn[m]); E += abs(rk[m] - J->t[m]); }
            if (E < bestE) { bestE = E; memcpy(B, W, sizeof W); }
            double lr = log(T1 / T0);
            for (long it = 0; it < iters && E > 0; it++) {
                if ((it & 0x3FFF) == 0 && now_s() - t_start >= budget) break;
                double T = T0 * exp(lr * (double)it / (double)iters);
                int k = (int)(rng_next(&rs) % (uint64_t)npairs), i = pairs[k][0], j = pairs[k][1];
                uint8_t old = W[i][j], nw = (uint8_t)((old + 1 + rng_next(&rs) % (uint64_t)(P - 1)) % P);
                W[i][j] = W[j][i] = nw;
                int p = party[i], q = party[j], ns = nsep[p][q], nr[16], dE = 0;
                for (int s = 0; s < ns; s++) {
                    int m = sep[p][q][s];
                    nr[s] = rankp(W, xi[m], xn[m], yi[m], yn[m]);
                    dE += abs(nr[s] - J->t[m]) - abs(rk[m] - J->t[m]);
                }
                if (dE <= 0 || rng_unit(&rs) < exp(-dE / T)) {
                    for (int s = 0; s < ns; s++) rk[sep[p][q][s]] = nr[s];
                    E += dE;
                    if (E < bestE) { bestE = E; memcpy(B, W, sizeof W); }
                } else { W[i][j] = W[j][i] = old; }
                steps++;
            }
            if (E == 0) { found = 1; memcpy(B, W, sizeof W); bestE = 0; }
        }
        free(pairs);
        if (found || dump_best) {
            len = found ? snprintf(line, sizeof line, "%d %d FOUND %d", J->id, J->lam, N)
                        : snprintf(line, sizeof line, "%d %d BEST %d %ld %ld %d", J->id, J->lam, bestE, restarts, steps, N);
            for (int i = 0; i < N; i++) {
                len += snprintf(line + len, sizeof line - len, " ");
                for (int j = 0; j < N; j++) len += snprintf(line + len, sizeof line - len, "%d", B[i][j]);
            }
            snprintf(line + len, sizeof line - len, "\n");
        } else snprintf(line, sizeof line, "%d %d BEST %d %ld %ld\n", J->id, J->lam, bestE, restarts, steps);
        #pragma omp critical
        { fputs(line, stdout); fflush(stdout); }
    }
    free(jobs);
    return 0;
}
