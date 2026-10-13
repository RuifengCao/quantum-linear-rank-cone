/* s32: unimodular realisations of known extreme rays of Stab_5 (round R38, roadmap item S5(a)).
 *
 * A realisation of kappa * r is unimodular if it is an integer symmetric matrix W with zero diagonal on
 * N = kappa * (r(A) + ... + r(F)) vertices, party q owning kappa * r(q) consecutive vertices (purifier F last),
 * such that every visible cut block W[Q_X, Q_{N \ X}] (X a nonempty set of parties among A..E) has Smith invariants
 * in {0, 1} and rational rank kappa * r(X); equivalently, rank kappa * r(X) over GF(p) for EVERY prime p.  Reduced
 * mod d, such a W is a graph state of Z_d qudits with entropy vector kappa * r in units of log d, for every d >= 2
 * (paper, Appendix A 1 and Appendix G 2).
 *
 * This program only searches; scripts/s32_unimodular.py certifies every reported matrix exactly (integer Smith
 * forms), independently of this code.  It anneals over the inter-party entries so that the cut ranks are right over
 * GF(p) for every p in an active set of primes (default 2, 3, 5, 7 and 2^31 - 1, which stands in for Q).  When the
 * active set is satisfied, the ranks are also checked for every prime below 1000 and for 2^61 - 1; a prime that fails
 * joins the active set and the search goes on.  A larger prime can still divide a Smith invariant, which is why the
 * exact certification exists.
 *
 * Two search spaces:
 *   plain  entries of W in {-a, ..., a} (a = amax), one variable per inter-party pair of vertices;
 *   ring   restriction of scalars from a quadratic order R = Z[theta], theta^2 = s theta + t.  The job describes
 *          the ray at multiplicity kappa_R; the search variables are entries v = x + y theta of a symmetric matrix V
 *          over R (x, y in {-a..a}), and W is V with every entry replaced by the symmetric 2 x 2 integer block
 *          G M(v), where M(v) is multiplication by v in the basis {1, theta} (columns: v * 1, v * theta) and
 *          G = [[0, 1], [1, s]] is the Gram matrix of the unimodular trace form (x, y) -> Tr(x y / (2 theta - s)).
 *          Each R-vertex becomes two integer vertices of the same party, so W realises (at best)
 *          kappa = 2 kappa_R; targets are doubled internally.
 *          G is unimodular and the cut blocks of W are (I (x) G) M(V_X), so W is unimodular whenever every V_X has
 *          unit Smith invariants over R (R a principal ideal domain); the certification does not rely on this.
 *
 * Input (stdin), one job per line:  id kappa n0 .. n5 t1 .. t31 S [start]
 *   t_m: target rank of the cut with party mask m (already multiplied by kappa; in ring mode the R-level values).
 *   S = 1: a start matrix follows, N*N integers row-major (plain mode: a matrix W, e.g. a qubit certificate;
 *   ring mode: the x-parts of V, y = 0).  Restarts after the first randomise the signs of its nonzero entries.
 *   S = 0: random starts.
 * Arguments: seconds seed [amax T0 T1 sweeps p_signflip dump_best primes extra ring p_neutral]
 *   primes: initial active set, comma-separated (default 2,3,5,7,2147483647); extra = 0 skips the extra checks
 *   (FOUND then only means: right ranks over the active primes); ring = "s,t" selects ring mode (default: plain);
 *   p_neutral: probability of a neutral move per step (default 0.2).  Neutral moves keep every cut block in the
 *   same Smith class, so they never change the residual; they move the search to another representative:
 *     basis change  V[u'][.] += e V[u][.] outside the party of u, u' (u, u' in one party, e = +-1), i.e. W -> E^T W E
 *                   with E party-local and unimodular;
 *     local complementation at c  V[u][u'] -= a V[u][c] V[u'][c] for neighbours u, u' of c in different parties
 *                   (a = +-1), the symplectic map (x_c, z_c) -> (x_c + a z_c, z_c) on c followed by party-local
 *                   gates that clear the diagonal and the intra-party entries;
 *   both over Z (plain) or over R (ring; there the integer image is a local symplectic map too).  A neutral move is
 *   made only if every new entry stays in the alphabet.
 * Output, one line per job (kappa = the multiplicity of the printed integer matrix):
 *   id kappa FOUND N seconds restarts steps primes row_0 ... row_{N-1}
 *   id kappa BEST E restarts steps seconds primes [N row_0 ... row_{N-1}]   (rows only with dump_best = 1)
 *   (rows: N comma-separated integers; primes: the final active set; E: residual over it)
 * Build: gcc -O3 -march=native -fopenmp -o s32_unimodular s32_unimodular.c -lm
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
#define MAXP 32
#define BIG31 2147483647ULL                    /* 2^31 - 1 */
#define BIG61 2305843009213693951ULL           /* 2^61 - 1 */

typedef struct { int id, lam, n[6], t[32], has_start; int8_t *start; } Job;

static inline uint64_t rng_next(uint64_t *s) {
    uint64_t x = *s; x ^= x >> 12; x ^= x << 25; x ^= x >> 27; *s = x;
    return x * 0x2545F4914F6CDD1DULL;
}
static inline double rng_unit(uint64_t *s) { return (rng_next(s) >> 11) * (1.0 / 9007199254740992.0); }
static double now_s(void) { struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts); return ts.tv_sec + ts.tv_nsec * 1e-9; }

/* ---------------------------------------------------------------- ranks of W[ri, ci] over GF(p) */
static int rank2(int8_t W[MAXN][MAXN], const int *ri, int nr, const int *ci, int nc) {
    uint64_t basis[MAXN]; int r = 0;
    memset(basis, 0, sizeof(uint64_t) * (size_t)nc);
    for (int a = 0; a < nr; a++) {
        uint64_t v = 0;
        const int8_t *row = W[ri[a]];
        for (int b = 0; b < nc; b++) v |= (uint64_t)(row[ci[b]] & 1) << b;
        while (v) {
            int lead = 63 - __builtin_clzll(v);
            if (!basis[lead]) { basis[lead] = v; r++; break; }
            v ^= basis[lead];
        }
    }
    return r;
}

/* small primes: the modulus is a compile-time constant inside each instance, so % is cheap */
#define DEF_RANK_SMALL(P)                                                                                        \
static int rank_p##P(int8_t W[MAXN][MAXN], const int *ri, int nr, const int *ci, int nc) {                      \
    uint8_t M[MAXN][MAXN]; uint8_t inv[P];                                                                       \
    for (int x = 1; x < P; x++) for (int y = 1; y < P; y++) if ((x * y) % P == 1) inv[x] = (uint8_t)y;            \
    for (int a = 0; a < nr; a++) { const int8_t *row = W[ri[a]];                                                 \
        for (int b = 0; b < nc; b++) { int v = row[ci[b]] % P; M[a][b] = (uint8_t)(v < 0 ? v + P : v); } }       \
    int r = 0;                                                                                                   \
    for (int c = 0; c < nc && r < nr; c++) {                                                                     \
        int piv = -1;                                                                                            \
        for (int a = r; a < nr; a++) if (M[a][c]) { piv = a; break; }                                            \
        if (piv < 0) continue;                                                                                   \
        if (piv != r) for (int b = c; b < nc; b++) { uint8_t t = M[r][b]; M[r][b] = M[piv][b]; M[piv][b] = t; }   \
        int iv = inv[M[r][c]];                                                                                   \
        for (int a = r + 1; a < nr; a++) if (M[a][c]) {                                                          \
            int f = (P - (M[a][c] * iv) % P) % P;                                                                \
            for (int b = c; b < nc; b++) M[a][b] = (uint8_t)((M[a][b] + f * M[r][b]) % P);                       \
        }                                                                                                        \
        r++;                                                                                                     \
    }                                                                                                            \
    return r;                                                                                                    \
}
DEF_RANK_SMALL(3) DEF_RANK_SMALL(5) DEF_RANK_SMALL(7) DEF_RANK_SMALL(11) DEF_RANK_SMALL(13)

static uint64_t mulmod_any(uint64_t a, uint64_t b, uint64_t p) { return (uint64_t)((unsigned __int128)a * b % p); }
static uint64_t powmod_any(uint64_t a, uint64_t e, uint64_t p) {
    uint64_t r = 1; a %= p;
    while (e) { if (e & 1) r = mulmod_any(r, a, p); a = mulmod_any(a, a, p); e >>= 1; }
    return r;
}
/* any prime p < 2^62 (slow path: extra checks) */
static int rank_any(int8_t W[MAXN][MAXN], const int *ri, int nr, const int *ci, int nc, uint64_t p) {
    uint64_t M[MAXN][MAXN];
    for (int a = 0; a < nr; a++) { const int8_t *row = W[ri[a]];
        for (int b = 0; b < nc; b++) { int64_t v = row[ci[b]]; M[a][b] = v < 0 ? p - (uint64_t)(-v) % p : (uint64_t)v % p; } }
    int r = 0;
    for (int c = 0; c < nc && r < nr; c++) {
        int piv = -1;
        for (int a = r; a < nr; a++) if (M[a][c]) { piv = a; break; }
        if (piv < 0) continue;
        if (piv != r) for (int b = c; b < nc; b++) { uint64_t t = M[r][b]; M[r][b] = M[piv][b]; M[piv][b] = t; }
        uint64_t iv = powmod_any(M[r][c], p - 2, p);
        for (int a = r + 1; a < nr; a++) if (M[a][c]) {
            uint64_t f = p - mulmod_any(M[a][c], iv, p);
            for (int b = c; b < nc; b++) M[a][b] = (M[a][b] + mulmod_any(f, M[r][b], p)) % p;
        }
        r++;
    }
    return r;
}
/* the large prime 2^31 - 1 (stand-in for the rational rank): fast reduction */
static inline uint64_t red31(uint64_t x) { x = (x & BIG31) + (x >> 31); x = (x & BIG31) + (x >> 31); return x >= BIG31 ? x - BIG31 : x; }
static int rank_big31(int8_t W[MAXN][MAXN], const int *ri, int nr, const int *ci, int nc) {
    uint64_t M[MAXN][MAXN];
    for (int a = 0; a < nr; a++) { const int8_t *row = W[ri[a]];
        for (int b = 0; b < nc; b++) { int v = row[ci[b]]; M[a][b] = v < 0 ? BIG31 - (uint64_t)(-v) : (uint64_t)v; } }
    int r = 0;
    for (int c = 0; c < nc && r < nr; c++) {
        int piv = -1;
        for (int a = r; a < nr; a++) if (M[a][c]) { piv = a; break; }
        if (piv < 0) continue;
        if (piv != r) for (int b = c; b < nc; b++) { uint64_t t = M[r][b]; M[r][b] = M[piv][b]; M[piv][b] = t; }
        uint64_t iv = powmod_any(M[r][c], BIG31 - 2, BIG31);
        for (int a = r + 1; a < nr; a++) if (M[a][c]) {
            uint64_t f = BIG31 - red31(M[a][c] * iv);
            for (int b = c; b < nc; b++) M[a][b] = red31(M[a][b] + f * M[r][b]);
        }
        r++;
    }
    return r;
}

static int rank_dispatch(uint64_t p, int8_t W[MAXN][MAXN], const int *ri, int nr, const int *ci, int nc) {
    if (nr == 0 || nc == 0) return 0;
    switch (p) {
        case 2: return rank2(W, ri, nr, ci, nc);
        case 3: return rank_p3(W, ri, nr, ci, nc);
        case 5: return rank_p5(W, ri, nr, ci, nc);
        case 7: return rank_p7(W, ri, nr, ci, nc);
        case 11: return rank_p11(W, ri, nr, ci, nc);
        case 13: return rank_p13(W, ri, nr, ci, nc);
        case BIG31: return rank_big31(W, ri, nr, ci, nc);
        default: return rank_any(W, ri, nr, ci, nc, p);
    }
}

/* extra checks once the active primes are satisfied: every prime below EXTRA_BOUND, then 2^31 - 1 and 2^61 - 1 */
#define EXTRA_BOUND 1000
static uint64_t EXTRA[200]; static int NEXTRA = 0;
static void build_extra(void) {
    for (int n = 2; n < EXTRA_BOUND; n++) {
        int pr = 1;
        for (int d = 2; d * d <= n; d++) if (n % d == 0) { pr = 0; break; }
        if (pr) EXTRA[NEXTRA++] = (uint64_t)n;
    }
    EXTRA[NEXTRA++] = BIG31; EXTRA[NEXTRA++] = BIG61;
}

/* ring mode: theta^2 = RS * theta + RT; G = [[0, 1], [1, RS]];  G M(x + y theta) =
 *   [[ y,          x + RS y              ],
 *    [ x + RS y,   RT y + RS (x + RS y)  ]]                                                                     */
static int RING = 0, RS = 0, RT = 0;
static void ring_block(int x, int y, int blk[2][2]) {
    int m01 = x + RS * y;
    blk[0][0] = y; blk[0][1] = m01; blk[1][0] = m01; blk[1][1] = RT * y + RS * m01;
}

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: %s seconds seed [amax T0 T1 sweeps p_signflip dump_best primes extra ring]\n", argv[0]); return 2; }
    build_extra();
    double budget = atof(argv[1]); uint64_t seed0 = strtoull(argv[2], 0, 10);
    int amax = argc > 3 ? atoi(argv[3]) : 1;
    double T0 = argc > 4 ? atof(argv[4]) : 1.0, T1 = argc > 5 ? atof(argv[5]) : 0.05;
    double sweeps = argc > 6 ? atof(argv[6]) : 50.0;
    double p_sign = argc > 7 ? atof(argv[7]) : 0.6;
    int dump_best = argc > 8 ? atoi(argv[8]) : 0;
    uint64_t act0[MAXP]; int na0 = 0;
    {
        const char *pl = argc > 9 ? argv[9] : "2,3,5,7,2147483647";
        char buf[512]; strncpy(buf, pl, sizeof buf - 1); buf[sizeof buf - 1] = 0;
        for (char *tok = strtok(buf, ","); tok && na0 < MAXP; tok = strtok(0, ",")) act0[na0++] = strtoull(tok, 0, 10);
    }
    int do_extra = argc > 10 ? atoi(argv[10]) : 1;
    if (argc > 11 && strcmp(argv[11], "0") != 0) {
        if (sscanf(argv[11], "%d,%d", &RS, &RT) != 2) { fprintf(stderr, "ring must be s,t\n"); return 2; }
        RING = 1;
    }
    double p_neu = argc > 12 ? atof(argv[12]) : 0.2;
    if (amax < 1 || amax > 3) { fprintf(stderr, "amax must be 1..3\n"); return 2; }
    int cap = 256, nj = 0; Job *jobs = malloc(cap * sizeof(Job));
    while (1) {
        Job J; memset(&J, 0, sizeof J);
        if (scanf("%d %d", &J.id, &J.lam) != 2) break;
        int N = 0;
        for (int p = 0; p < 6; p++) { if (scanf("%d", &J.n[p]) != 1) return 3; N += J.n[p]; }
        for (int m = 1; m < 32; m++) if (scanf("%d", &J.t[m]) != 1) return 3;
        if (scanf("%d", &J.has_start) != 1) return 3;
        if (J.has_start) {
            J.start = malloc((size_t)N * N);
            for (int k = 0; k < N * N; k++) { int v; if (scanf("%d", &v) != 1) return 3; J.start[k] = (int8_t)v; }
        }
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
    const int mult = RING ? 2 : 1;                     /* integer vertices per search vertex */
    #pragma omp parallel for schedule(dynamic, 1)
    for (int jj = 0; jj < nj; jj++) {
        Job *J = &jobs[jj];
        int NV = 0, vparty[MAXN];                      /* search vertices */
        for (int p = 0; p < 6; p++) for (int k = 0; k < J->n[p]; k++) { if (NV < MAXN) vparty[NV] = p; NV++; }
        int N = mult * NV, party[MAXN], tgt[32];
        char *line = malloc(1 << 16); int len = 0;
        if (N > MAXN || NV == 0) {
            snprintf(line, 1 << 16, "%d %d SKIP N=%d\n", J->id, mult * J->lam, N);
            #pragma omp critical
            { fputs(line, stdout); fflush(stdout); }
            free(line); continue;
        }
        for (int i = 0; i < N; i++) party[i] = vparty[i / mult];
        for (int m = 1; m < 32; m++) tgt[m] = mult * J->t[m];
        int xi[32][MAXN], xn[32], yi[32][MAXN], yn[32];
        for (int m = 1; m < 32; m++) {
            xn[m] = yn[m] = 0;
            for (int i = 0; i < N; i++) {
                int in = party[i] < 5 && ((m >> party[i]) & 1);
                if (in) xi[m][xn[m]++] = i; else yi[m][yn[m]++] = i;
            }
        }
        int npairs = 0, (*pairs)[2] = malloc(sizeof(int[2]) * NV * NV);
        for (int i = 0; i < NV; i++) for (int j = i + 1; j < NV; j++) if (vparty[i] != vparty[j]) { pairs[npairs][0] = i; pairs[npairs][1] = j; npairs++; }
        uint64_t rs = seed0 ^ (0x9E3779B97F4A7C15ULL * (uint64_t)(J->id + 1)) ^ ((uint64_t)J->lam << 40);
        for (int w = 0; w < 8; w++) rng_next(&rs);
        uint64_t act[MAXP]; int na = 0;
        for (int a = 0; a < na0; a++) act[na++] = act0[a];
        int8_t W[MAXN][MAXN], Bst[MAXN][MAXN];
        int8_t VX[MAXN][MAXN], VY[MAXN][MAXN];         /* search variables: plain uses VX only */
        int rk[MAXP][32], bestE = 1 << 30, found = 0; long restarts = 0, steps = 0;
        memset(Bst, 0, sizeof Bst);
        double t_start = now_s();
        long iters = (long)(sweeps * (npairs > 0 ? npairs : 1));
#define APPLY(u, v) do {                                                                                     \
            if (!RING) { W[u][v] = W[v][u] = VX[u][v]; }                                                    \
            else { int blk_[2][2]; ring_block(VX[u][v], VY[u][v], blk_);                                    \
                   for (int e_ = 0; e_ < 2; e_++) for (int f_ = 0; f_ < 2; f_++) {                           \
                       W[2 * (u) + e_][2 * (v) + f_] = (int8_t)blk_[e_][f_];                                \
                       W[2 * (v) + f_][2 * (u) + e_] = (int8_t)blk_[e_][f_]; } }                            \
        } while (0)
        while (!found && now_s() - t_start < budget) {
            /* ---- start state */
            memset(W, 0, sizeof W); memset(VX, 0, sizeof VX); memset(VY, 0, sizeof VY);
            for (int k = 0; k < npairs; k++) {
                int i = pairs[k][0], j = pairs[k][1], x, y = 0;
                if (J->has_start) {
                    x = J->start[i * NV + j];
                    if (restarts > 0 && x) x = (rng_next(&rs) & 1) ? x : -x;
                } else {
                    x = (int)(rng_next(&rs) % (uint64_t)(2 * amax + 1)) - amax;
                    if (RING) y = (int)(rng_next(&rs) % (uint64_t)(2 * amax + 1)) - amax;
                }
                VX[i][j] = VX[j][i] = (int8_t)x; VY[i][j] = VY[j][i] = (int8_t)y;
                APPLY(i, j);
            }
            restarts++;
            int E = 0;
            for (int a = 0; a < na; a++) for (int m = 1; m < 32; m++) {
                rk[a][m] = rank_dispatch(act[a], W, xi[m], xn[m], yi[m], yn[m]);
                E += abs(rk[a][m] - tgt[m]);
            }
            if (E < bestE) { bestE = E; if (dump_best) memcpy(Bst, W, sizeof W); }
            double lr = log(T1 / T0);
            for (long it = 0; it < iters; it++) {
                if (E == 0) {
                    /* ---- extra checks: every prime below 1000 and two large primes */
                    int bad = -1;
                    for (int e = 0; e < NEXTRA && bad < 0 && do_extra; e++) {
                        int inact = 0;
                        for (int a = 0; a < na; a++) if (act[a] == EXTRA[e]) inact = 1;
                        if (inact) continue;
                        for (int m = 1; m < 32; m++)
                            if (rank_dispatch(EXTRA[e], W, xi[m], xn[m], yi[m], yn[m]) != tgt[m]) { bad = e; break; }
                    }
                    if (bad < 0) { found = 1; break; }
                    if (na < MAXP) {
                        act[na] = EXTRA[bad];
                        for (int m = 1; m < 32; m++) {
                            rk[na][m] = rank_dispatch(act[na], W, xi[m], xn[m], yi[m], yn[m]);
                            E += abs(rk[na][m] - tgt[m]);
                        }
                        na++;
                    } else break;
                }
                if ((it & 0x3FF) == 0 && now_s() - t_start >= budget) break;
                double T = T0 * exp(lr * (double)it / (double)iters);
                if (p_neu > 0 && rng_unit(&rs) < p_neu) {
                    /* ---- neutral move: the residual and every rank stay exactly the same */
                    int nchg = 0, cu[MAXN * MAXN / 2], cv[MAXN * MAXN / 2], cx[MAXN * MAXN / 2], cy[MAXN * MAXN / 2], ok = 1;
                    if (rng_unit(&rs) < 0.5) {                       /* basis change inside one party */
                        int a1 = (int)(rng_next(&rs) % (uint64_t)NV), a2 = (int)(rng_next(&rs) % (uint64_t)NV);
                        int e = (rng_next(&rs) & 1) ? 1 : -1;
                        if (a1 == a2 || vparty[a1] != vparty[a2]) ok = 0;
                        for (int w = 0; w < NV && ok; w++) {
                            if (vparty[w] == vparty[a1]) continue;
                            int nx2 = VX[a2][w] + e * VX[a1][w], ny2 = VY[a2][w] + e * VY[a1][w];
                            if (nx2 == VX[a2][w] && ny2 == VY[a2][w]) continue;
                            if (abs(nx2) > amax || abs(ny2) > amax) { ok = 0; break; }
                            cu[nchg] = a2; cv[nchg] = w; cx[nchg] = nx2; cy[nchg] = ny2; nchg++;
                        }
                    } else {                                         /* local complementation at c */
                        int c = (int)(rng_next(&rs) % (uint64_t)NV), a = (rng_next(&rs) & 1) ? 1 : -1;
                        for (int i1 = 0; i1 < NV && ok; i1++) {
                            if (vparty[i1] == vparty[c] || (!VX[i1][c] && !VY[i1][c])) continue;
                            for (int i2 = i1 + 1; i2 < NV; i2++) {
                                if (vparty[i2] == vparty[c] || vparty[i2] == vparty[i1] || (!VX[i2][c] && !VY[i2][c])) continue;
                                /* product (x1 + y1 th)(x2 + y2 th) = x1 x2 + RT y1 y2 + (x1 y2 + x2 y1 + RS y1 y2) th */
                                int x1 = VX[i1][c], y1 = VY[i1][c], x2 = VX[i2][c], y2 = VY[i2][c];
                                int px = x1 * x2 + RT * y1 * y2, py = x1 * y2 + x2 * y1 + RS * y1 * y2;
                                int nx2 = VX[i1][i2] - a * px, ny2 = VY[i1][i2] - a * py;
                                if (abs(nx2) > amax || abs(ny2) > amax) { ok = 0; break; }
                                cu[nchg] = i1; cv[nchg] = i2; cx[nchg] = nx2; cy[nchg] = ny2; nchg++;
                            }
                        }
                    }
                    if (ok && nchg) {
                        for (int z = 0; z < nchg; z++) {
                            int i1 = cu[z], i2 = cv[z];
                            VX[i1][i2] = VX[i2][i1] = (int8_t)cx[z]; VY[i1][i2] = VY[i2][i1] = (int8_t)cy[z];
                            APPLY(i1, i2);
                        }
                    }
                    steps++;
                    continue;
                }
                int k = (int)(rng_next(&rs) % (uint64_t)npairs), u = pairs[k][0], v = pairs[k][1];
                int ox = VX[u][v], oy = VY[u][v], nx, ny;
                int nonzero = ox != 0 || oy != 0;
                if (nonzero && rng_unit(&rs) < p_sign) { nx = -ox; ny = -oy; }
                else {
                    do {
                        nx = (int)(rng_next(&rs) % (uint64_t)(2 * amax + 1)) - amax;
                        ny = RING ? (int)(rng_next(&rs) % (uint64_t)(2 * amax + 1)) - amax : 0;
                    } while (nx == ox && ny == oy);
                }
                VX[u][v] = VX[v][u] = (int8_t)nx; VY[u][v] = VY[v][u] = (int8_t)ny;
                APPLY(u, v);
                int p = vparty[u], q = vparty[v], ns = nsep[p][q], nr[MAXP][16], dE = 0;
                /* the parity of W is unchanged by a sign flip in plain mode */
                int same2 = !RING && (((ox ^ nx) & 1) == 0);
                for (int a = 0; a < na; a++) {
                    for (int s = 0; s < ns; s++) {
                        int m = sep[p][q][s];
                        nr[a][s] = (act[a] == 2 && same2) ? rk[a][m] : rank_dispatch(act[a], W, xi[m], xn[m], yi[m], yn[m]);
                        dE += abs(nr[a][s] - tgt[m]) - abs(rk[a][m] - tgt[m]);
                    }
                }
                if (dE <= 0 || rng_unit(&rs) < exp(-dE / T)) {
                    for (int a = 0; a < na; a++) for (int s = 0; s < ns; s++) rk[a][sep[p][q][s]] = nr[a][s];
                    E += dE;
                    if (E < bestE) { bestE = E; if (dump_best) memcpy(Bst, W, sizeof W); }
                } else {
                    VX[u][v] = VX[v][u] = (int8_t)ox; VY[u][v] = VY[v][u] = (int8_t)oy;
                    APPLY(u, v);
                }
                steps++;
            }
        }
#undef APPLY
        free(pairs);
        double secs = now_s() - t_start;
        char ps[512]; int pl = 0; ps[0] = 0;
        for (int a = 0; a < na; a++) pl += snprintf(ps + pl, sizeof ps - pl, a ? ",%llu" : "%llu", (unsigned long long)act[a]);
        int8_t (*OUT)[MAXN] = found ? W : Bst;
        if (found || dump_best) {
            if (found) len = snprintf(line, 1 << 16, "%d %d FOUND %d %.3f %ld %ld %s", J->id, mult * J->lam, N, secs, restarts, steps, ps);
            else len = snprintf(line, 1 << 16, "%d %d BEST %d %ld %ld %.3f %s %d", J->id, mult * J->lam, bestE, restarts, steps, secs, ps, N);
            for (int i = 0; i < N; i++) {
                len += snprintf(line + len, (1 << 16) - len, " ");
                for (int j = 0; j < N; j++) len += snprintf(line + len, (1 << 16) - len, j ? ",%d" : "%d", OUT[i][j]);
            }
            snprintf(line + len, (1 << 16) - len, "\n");
        } else snprintf(line, 1 << 16, "%d %d BEST %d %ld %ld %.3f %s\n", J->id, mult * J->lam, bestE, restarts, steps, secs, ps);
        #pragma omp critical
        { fputs(line, stdout); fflush(stdout); }
        free(line);
    }
    for (int k = 0; k < nj; k++) free(jobs[k].start);
    free(jobs);
    return 0;
}
