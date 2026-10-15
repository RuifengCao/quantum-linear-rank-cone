/* s37: adversarial search for Theorem 12 (docs/R32_reduction_and_local_dimension.md, section 4.11).  Simulated
 * annealing over qubit graph states of seven systems (parties 0..5 visible, system 6 the purifier), system p owning
 * n_p consecutive qubits, no edges inside a system (party-local controlled-Z gates remove them).  It minimises a
 * linear form F of the entropy vector, given on the 63 masks of visible parties (bit i = party i; the purifier side
 * is folded in by purity).  Entropy of a cut = rank over GF(2) of the adjacency block W[X, X^c].
 *
 * Optional pull towards a target vector t (63 integers): objective = F(S) + mu * sum_m |S_m - t_m|.  The minimum of
 * F itself over all visited graphs is tracked separately and its graph is printed.
 *
 * Input (stdin):  F_1 .. F_63   [mu t_1 .. t_63]
 * Arguments:      n0 .. n6 seconds seed [T0]
 * Output:         MIN <min F seen> OBJ <best objective> RESTARTS <r> STEPS <s>
 *                 GRAPH <N> hexrow_0 ... hexrow_{N-1}      (the graph attaining MIN; row i: bit j = edge ij)
 * The search proves nothing: the Python wrapper (s34_six_party_survey.transfer_anneal) recomputes the entropies of the
 * printed graph exactly and checks every step of the proof on it.
 * Build: gcc -O3 -march=native -o s37_minform7 s37_minform7.c -lm
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define MAXN 64
static int N, owner[MAXN], sizes[7];
static long long Fc[64];
static int Tg[64];
static double mu = 0.0;
static uint64_t adj[MAXN], memb[64];
static uint64_t rng_s;

static inline uint64_t rnd(void) { rng_s ^= rng_s << 13; rng_s ^= rng_s >> 7; rng_s ^= rng_s << 17; return rng_s; }
static inline double urand(void) { return (double)(rnd() >> 11) / 9007199254740992.0; }

static int rank_cut(int m) {
    uint64_t rows[MAXN];
    int nr = 0;
    uint64_t all = (N == 64) ? ~0ULL : ((1ULL << N) - 1);
    uint64_t colmask = ~memb[m] & all;
    for (int v = 0; v < N; v++)
        if (memb[m] >> v & 1) rows[nr++] = adj[v] & colmask;
    int r = 0;
    for (int i = 0; i < nr; i++) {
        uint64_t x = rows[i];
        if (!x) continue;
        uint64_t low = x & (~x + 1);
        r++;
        for (int j = i + 1; j < nr; j++)
            if (rows[j] & low) rows[j] ^= x;
    }
    return r;
}

static int S[64];
static long long fval, objv;

static void full_eval(void) {
    fval = 0; long long dist = 0;
    for (int m = 1; m < 64; m++) { S[m] = rank_cut(m); fval += Fc[m] * S[m]; dist += llabs((long long)(S[m] - Tg[m])); }
    objv = fval + (long long)llround(mu * (double)dist);
}

/* masks whose cut separates systems p and q (q may be the purifier 6) */
static int sep[7][7][64], nsep[7][7];

int main(int argc, char **argv) {
    if (argc < 10) { fprintf(stderr, "usage: s37_minform7 n0..n6 seconds seed [T0] < F (63 ints) [mu t (63 ints)]\n"); return 1; }
    for (int p = 0; p < 7; p++) sizes[p] = atoi(argv[1 + p]);
    double secs = atof(argv[8]);
    rng_s = strtoull(argv[9], 0, 10) * 2654435761ULL + 88172645463325252ULL;
    double T0 = argc > 10 ? atof(argv[10]) : 3.0;
    for (int m = 1; m < 64; m++) if (scanf("%lld", &Fc[m]) != 1) return 2;
    if (scanf("%lf", &mu) == 1) { for (int m = 1; m < 64; m++) if (scanf("%d", &Tg[m]) != 1) return 3; }
    else mu = 0.0;
    N = 0;
    for (int p = 0; p < 7; p++) for (int j = 0; j < sizes[p]; j++) owner[N++] = p;
    if (N > MAXN) return 4;
    for (int m = 0; m < 64; m++) { memb[m] = 0; for (int v = 0; v < N; v++) if (owner[v] < 6 && (m >> owner[v] & 1)) memb[m] |= 1ULL << v; }
    for (int p = 0; p < 7; p++) for (int q = 0; q < 7; q++) {
        nsep[p][q] = 0;
        for (int m = 1; m < 64; m++) {
            int ip = p < 6 ? (m >> p & 1) : 0, iq = q < 6 ? (m >> q & 1) : 0;
            if (ip != iq) sep[p][q][nsep[p][q]++] = m;
        }
    }
    long long minF = 1LL << 60, bestobj = 1LL << 60;
    uint64_t minadj[MAXN];
    memset(minadj, 0, sizeof minadj);
    clock_t t0 = clock();
    long long steps = 0;
    int restarts = 0;
    while ((double)(clock() - t0) / CLOCKS_PER_SEC < secs) {
        double dens = 0.1 + 0.8 * urand();
        memset(adj, 0, sizeof adj);
        for (int u = 0; u < N; u++)
            for (int v = u + 1; v < N; v++)
                if (owner[u] != owner[v] && urand() < dens) { adj[u] |= 1ULL << v; adj[v] |= 1ULL << u; }
        full_eval();
        if (fval < minF) { minF = fval; memcpy(minadj, adj, sizeof adj); }
        if (objv < bestobj) bestobj = objv;
        double T = T0;
        for (long it = 0; it < 200000; it++) {
            int u = (int)(rnd() % (uint64_t)N), v = (int)(rnd() % (uint64_t)N);
            if (owner[u] == owner[v]) continue;
            steps++;
            adj[u] ^= 1ULL << v; adj[v] ^= 1ULL << u;
            int p = owner[u], q = owner[v];
            long long df = 0, dd = 0;
            int newS[64];
            for (int k = 0; k < nsep[p][q]; k++) {
                int m = sep[p][q][k];
                int r = rank_cut(m);
                newS[k] = r;
                df += Fc[m] * (r - S[m]);
                dd += llabs((long long)(r - Tg[m])) - llabs((long long)(S[m] - Tg[m]));
            }
            long long dE = df + (long long)llround(mu * (double)dd);
            if (dE <= 0 || urand() < exp(-(double)dE / T)) {
                for (int k = 0; k < nsep[p][q]; k++) S[sep[p][q][k]] = newS[k];
                fval += df; objv += dE;
                if (fval < minF) { minF = fval; memcpy(minadj, adj, sizeof adj); }
                if (objv < bestobj) bestobj = objv;
            } else { adj[u] ^= 1ULL << v; adj[v] ^= 1ULL << u; }
            T *= 0.99995; if (T < 0.05) T = 0.05;
        }
        restarts++;
    }
    printf("MIN %lld OBJ %lld RESTARTS %d STEPS %lld\n", minF, bestobj, restarts, steps);
    printf("GRAPH %d", N);
    for (int i = 0; i < N; i++) printf(" %llx", (unsigned long long)minadj[i]);
    printf("\n");
    return 0;
}
