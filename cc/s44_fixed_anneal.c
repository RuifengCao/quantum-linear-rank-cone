/* s44_fixed_anneal (R56): GF(2) search for a subspace arrangement with prescribed ranks, some parties FIXED.
 *
 * Driver: scripts/s44_fixed_search.py, which fixes parties whose subspaces the rank data force into a normal form
 * (three parties meeting pairwise in one common subspace whose sum is direct) and verifies every hit independently.
 * Six parties (bits 0..5).  Party i owns c[i] vectors of length N (N <= 64), its subspace is their span.
 * Target: rank of the union over every nonempty X equals h[X].  Fixed parties keep their vectors; the others are
 * annealed.  Energy = sum_X |rank(X) - h(X)|.
 *
 * stdin:  N  seed  seconds
 *         h[1] .. h[63]
 *         for i = 0..5:  c[i]  fixed(0/1)  [c[i] hex vectors if fixed]
 * stdout: "FOUND <seconds>" and then the vectors of all parties (hex, party order), or "BEST <energy>".
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>

typedef uint64_t u64;
static int N, h[64], cnt[6], fixedp[6];
static u64 vec[6][64];
static int rk[64];
static u64 rng_s;
static inline u64 rnd(void) { rng_s ^= rng_s << 13; rng_s ^= rng_s >> 7; rng_s ^= rng_s << 17; return rng_s; }
static double now(void) { struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return t.tv_sec + 1e-9 * t.tv_nsec; }

static int rank_of(int X) {
    u64 basis[64]; memset(basis, 0, sizeof(u64) * N);
    int r = 0;
    for (int i = 0; i < 6; i++) if (X >> i & 1)
        for (int j = 0; j < cnt[i]; j++) {
            u64 v = vec[i][j];
            for (int b = N - 1; b >= 0 && v; b--) if (v >> b & 1) {
                if (!basis[b]) { basis[b] = v; r++; v = 0; break; }
                v ^= basis[b];
            }
        }
    return r;
}

static int energy_all(void) {
    int e = 0;
    for (int X = 1; X < 64; X++) { rk[X] = rank_of(X); e += abs(rk[X] - h[X]); }
    return e;
}

int main(void) {
    u64 seed; double secs;
    if (scanf("%d %llu %lf", &N, (unsigned long long *)&seed, &secs) != 3) return 2;
    for (int X = 1; X < 64; X++) if (scanf("%d", &h[X]) != 1) return 2;
    for (int i = 0; i < 6; i++) {
        if (scanf("%d %d", &cnt[i], &fixedp[i]) != 2) return 2;
        if (fixedp[i]) for (int j = 0; j < cnt[i]; j++) { unsigned long long x; if (scanf("%llx", &x) != 1) return 2; vec[i][j] = x; }
    }
    rng_s = seed * 0x9E3779B97F4A7C15ULL + 1;
    u64 mask = (N == 64) ? ~0ULL : ((1ULL << N) - 1);
    int var[64][2], nv = 0;
    for (int i = 0; i < 6; i++) if (!fixedp[i]) for (int j = 0; j < cnt[i]; j++) { var[nv][0] = i; var[nv][1] = j; nv++; }
    double t0 = now(); int best = 1 << 30; long long moves = 0;
    while (now() - t0 < secs) {
        for (int i = 0; i < 6; i++) if (!fixedp[i]) for (int j = 0; j < cnt[i]; j++) vec[i][j] = rnd() & mask;
        int e = energy_all();
        double T = 2.0;
        for (long long it = 0; it < 4000000 && now() - t0 < secs; it++) {
            T = 2.0 * pow(0.02 / 2.0, (double)(it % 400000) / 400000.0);
            int k = rnd() % nv, p = var[k][0], j = var[k][1];
            u64 old = vec[p][j], nw;
            int mv = rnd() % 3;
            if (mv == 0) nw = old ^ (1ULL << (rnd() % N));
            else if (mv == 1) { int q = rnd() % 6; nw = cnt[q] ? old ^ vec[q][rnd() % cnt[q]] : old; }
            else nw = rnd() & mask;
            if (nw == old || !nw) continue;
            vec[p][j] = nw;
            int de = 0, nr[64];
            for (int X = 1; X < 64; X++) if (X >> p & 1) { nr[X] = rank_of(X); de += abs(nr[X] - h[X]) - abs(rk[X] - h[X]); }
            moves++;
            if (de <= 0 || (double)(rnd() >> 11) * (1.0 / 9007199254740992.0) < exp(-de / T)) {
                for (int X = 1; X < 64; X++) if (X >> p & 1) rk[X] = nr[X];
                e += de;
                if (e < best) best = e;
                if (e == 0) {
                    printf("FOUND %.1f\n", now() - t0);
                    for (int i = 0; i < 6; i++) { for (int jj = 0; jj < cnt[i]; jj++) printf("%llx ", (unsigned long long)vec[i][jj]); printf("\n"); }
                    fflush(stdout); return 0;
                }
            } else vec[p][j] = old;
        }
    }
    printf("BEST %d moves %lld\n", best, moves);
    return 1;
}
