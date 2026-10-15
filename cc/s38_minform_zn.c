/* s38: adversarial search for Theorem 14 (docs/R32_reduction_and_local_dimension.md, section 4.11; R44).  Simulated
 * annealing over CSS states of codes over Z_N, N = 2^k (1 <= k <= 6): an r x m matrix over Z_N whose m columns are the
 * qudits, system p (0..5 visible, 6 the purifier) owning n_p consecutive columns.  By Theorem 1(d), step 3,
 *     S(T) = r(T) + r(E - T) - r(E),   r(T) = log2 |row space of the columns of T|,
 * and by Theorem 1(d) the CSS states generate the cone Stab^(Z_N), so a linear form that is negative somewhere on the
 * cone is negative at some CSS state.  The search minimises a linear form F of the entropy vector, given on the 63
 * masks of visible parties (bit i = party i; the purifier side folded in by purity).
 *
 * Optional pull towards a target vector t (63 integers): objective = F(S) + mu * sum_m |S_m - t_m|, masks with t_m < 0
 * left out (so a target on the seven single systems alone keeps the search away from states with small entropies).
 * The minimum of F itself over all visited codes is tracked separately and its code is printed.
 *
 * Input (stdin):  F_1 .. F_63   [mu t_1 .. t_63]
 * Arguments:      k r n0 .. n6 seconds seed [T0]
 * Output:         MIN <min F seen> OBJ <best objective> RESTARTS <r> STEPS <s>
 *                 CODE <r> <m> e_00 e_01 ... (the code attaining MIN, row-major, entries in 0..N-1)
 * Exhaustive mode: arguments k r n0 .. n6 lifts s nforms, stdin nforms forms (63 integers each) and a base code B
 * (r x m integers): every code B + 2^s P, P a 0/1 matrix (r m <= 30); output LIFTS <count>, then per form
 * FORM <i> MIN <min> ZEROS <codes at 0, first form only> CODE <r> <m> <a minimising code>.
 * The search proves nothing: the Python wrapper (s35_two_power_transfer.anneal) rebuilds the printed code, recomputes
 * its entropies exactly and checks every step of the proof on it.
 * Build: gcc -O3 -march=native -o s38_minform_zn s38_minform_zn.c -lm
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define MAXM 48
#define MAXR 16
static int K, NN, MASKV, R, M, owner[MAXM], sizes[7];
static int code[MAXR][MAXM];
static long long Fc[64];
static int Tg[64];
static double mu = 0.0;
static int cols_of[128][MAXM], ncols_of[128];
static uint64_t rng_s;

static inline uint64_t rnd(void) { rng_s ^= rng_s << 13; rng_s ^= rng_s >> 7; rng_s ^= rng_s << 17; return rng_s; }
static inline double urand(void) { return (double)(rnd() >> 11) / 9007199254740992.0; }

static inline int inv_odd(int u) {                 /* inverse of an odd u modulo 2^K (Newton) */
    unsigned x = (unsigned)u;
    for (int i = 0; i < 5; i++) x *= 2u - (unsigned)u * x;
    return (int)(x & (unsigned)MASKV);
}

/* log2 of the size of the row space of the columns of mask T (pivot of least 2-adic valuation, row elimination) */
static int logsize(int T) {
    int nc = ncols_of[T];
    if (!nc) return 0;
    int A[MAXR][MAXM];
    for (int i = 0; i < R; i++) for (int j = 0; j < nc; j++) A[i][j] = code[i][cols_of[T][j]];
    int total = 0, row = 0;
    while (row < R) {
        int bv = K, bi = -1, bj = -1;
        for (int i = row; i < R && bv > 0; i++)
            for (int j = 0; j < nc; j++) {
                int x = A[i][j];
                if (x) { int v = __builtin_ctz(x); if (v < bv) { bv = v; bi = i; bj = j; if (!v) break; } }
            }
        if (bi < 0) break;
        if (bi != row) for (int j = 0; j < nc; j++) { int t = A[row][j]; A[row][j] = A[bi][j]; A[bi][j] = t; }
        int uinv = inv_odd(A[row][bj] >> bv);
        for (int i = row + 1; i < R; i++) {
            int x = A[i][bj];
            if (!x) continue;
            int f = ((x >> bv) * uinv) & MASKV;
            for (int j = 0; j < nc; j++) A[i][j] = (A[i][j] - f * A[row][j]) & MASKV;
        }
        total += K - bv;
        row++;
    }
    return total;
}

static int S[64];
static long long fval, objv;

static void full_eval(void) {
    int r[128];
    for (int T = 1; T < 128; T++) r[T] = logsize(T);
    fval = 0; long long dist = 0;
    for (int m = 1; m < 64; m++) {
        S[m] = r[m] + r[127 ^ m] - r[127];
        fval += Fc[m] * S[m];
        if (Tg[m] >= 0) dist += llabs((long long)(S[m] - Tg[m]));
    }
    objv = fval + (long long)llround(mu * (double)dist);
}

static int rand_entry(void) {
    double u = urand();
    if (u < 0.25) return 0;
    if (u < 0.45) return 1;
    int x = (int)(rnd() & (uint64_t)MASKV);
    if (u < 0.75) x = (x << (int)(rnd() % (uint64_t)K)) & MASKV;    /* 2-adic scaling */
    return x;
}

/* Exhaustive mode: every code B + 2^s P with P a 0/1 matrix (B read from stdin after the forms); for each of the
 * nforms forms (63 coefficients each) the minimum, the number of codes attaining 0 (first form) and a minimising
 * code.  B = 0, s = 0 gives all 0/1 codes. */
static int lifts(int s, int nforms) {
    static long long Fs[8][64];
    if (nforms < 1 || nforms > 8) return 6;
    for (int f = 0; f < nforms; f++) for (int m = 1; m < 64; m++) if (scanf("%lld", &Fs[f][m]) != 1) return 2;
    int B[MAXR][MAXM];
    for (int i = 0; i < R; i++) for (int c = 0; c < M; c++) if (scanf("%d", &B[i][c]) != 1) return 3;
    int nbits = R * M;
    if (nbits > 30) return 7;
    long long minv[8], zeros = 0;
    static int arg[8][MAXR][MAXM];
    for (int f = 0; f < nforms; f++) minv[f] = 1LL << 60;
    for (long long P = 0; P < (1LL << nbits); P++) {
        for (int i = 0; i < R; i++) for (int c = 0; c < M; c++)
            code[i][c] = (B[i][c] + ((int)((P >> (i * M + c)) & 1) << s)) & MASKV;
        int r[128];
        for (int T = 1; T < 128; T++) r[T] = logsize(T);
        for (int f = 0; f < nforms; f++) {
            long long v = 0;
            for (int m = 1; m < 64; m++) v += Fs[f][m] * (r[m] + r[127 ^ m] - r[127]);
            if (v < minv[f]) { minv[f] = v; memcpy(arg[f], code, sizeof code); }
            if (f == 0 && v == 0) zeros++;
        }
    }
    printf("LIFTS %lld\n", 1LL << nbits);
    for (int f = 0; f < nforms; f++) {
        printf("FORM %d MIN %lld ZEROS %lld CODE %d %d", f, minv[f], f == 0 ? zeros : -1, R, M);
        for (int i = 0; i < R; i++) for (int c = 0; c < M; c++) printf(" %d", arg[f][i][c]);
        printf("\n");
    }
    return 0;
}

int main(int argc, char **argv) {
    if (argc < 12) {
        fprintf(stderr, "usage: s38_minform_zn k r n0..n6 seconds seed [T0] < F (63 ints) [mu t (63 ints)]\n"
                        "       s38_minform_zn k r n0..n6 lifts s nforms < F_1 .. F_nforms (63 ints each) B (r*m ints)\n");
        return 1;
    }
    K = atoi(argv[1]); R = atoi(argv[2]);
    if (K < 1 || K > 6 || R < 1 || R > MAXR) return 5;
    NN = 1 << K; MASKV = NN - 1;
    for (int p = 0; p < 7; p++) sizes[p] = atoi(argv[3 + p]);
    M = 0;
    for (int p = 0; p < 7; p++) for (int j = 0; j < sizes[p]; j++) { if (M >= MAXM) return 4; owner[M++] = p; }
    for (int T = 0; T < 128; T++) {
        ncols_of[T] = 0;
        for (int c = 0; c < M; c++) if (T >> owner[c] & 1) cols_of[T][ncols_of[T]++] = c;
    }
    if (!strcmp(argv[10], "lifts")) return lifts(atoi(argv[11]), atoi(argv[12]));
    double secs = atof(argv[10]);
    rng_s = strtoull(argv[11], 0, 10) * 2654435761ULL + 88172645463325252ULL;
    double T0 = argc > 12 ? atof(argv[12]) : 3.0;
    for (int m = 1; m < 64; m++) if (scanf("%lld", &Fc[m]) != 1) return 2;
    if (scanf("%lf", &mu) == 1) { for (int m = 1; m < 64; m++) if (scanf("%d", &Tg[m]) != 1) return 3; }
    else mu = 0.0;
    long long minF = 1LL << 60, bestobj = 1LL << 60;
    static int mincode[MAXR][MAXM];
    memset(mincode, 0, sizeof mincode);
    clock_t t0 = clock();
    long long steps = 0;
    int restarts = 0;
    while ((double)(clock() - t0) / CLOCKS_PER_SEC < secs) {
        for (int i = 0; i < R; i++) {
            int sh = (int)(rnd() % (uint64_t)K);
            for (int c = 0; c < M; c++) code[i][c] = (rand_entry() << (urand() < 0.3 ? sh : 0)) & MASKV;
        }
        full_eval();
        if (fval < minF) { minF = fval; memcpy(mincode, code, sizeof code); }
        if (objv < bestobj) bestobj = objv;
        double T = T0;
        for (long it = 0; it < 100000; it++) {
            if ((it & 1023) == 0 && (double)(clock() - t0) / CLOCKS_PER_SEC >= secs) break;
            int i = (int)(rnd() % (uint64_t)R), c = (int)(rnd() % (uint64_t)M);
            int oldcol[MAXR];
            for (int ii = 0; ii < R; ii++) oldcol[ii] = code[ii][c];
            double u = urand();
            if (u < 0.2 && M > 2) {                                    /* column c := a col c1 + b col c2 */
                int c1 = (int)(rnd() % (uint64_t)M), c2 = (int)(rnd() % (uint64_t)M);
                if (c1 == c || c2 == c) continue;
                int a = rand_entry() | (urand() < 0.5), b = rand_entry();
                for (int ii = 0; ii < R; ii++) code[ii][c] = (a * code[ii][c1] + b * code[ii][c2]) & MASKV;
            } else if (u < 0.35) {                                     /* add a multiple of another row */
                int i2 = (int)(rnd() % (uint64_t)R);
                if (i2 == i) continue;
                int f = (int)(rnd() & (uint64_t)MASKV);
                code[i][c] = (code[i][c] + f * code[i2][c]) & MASKV;
            } else code[i][c] = rand_entry();
            int same = 1;
            for (int ii = 0; ii < R; ii++) same &= code[ii][c] == oldcol[ii];
            if (same) continue;
            steps++;
            int Sold[64]; long long fold = fval, oold = objv;
            memcpy(Sold, S, sizeof S);
            full_eval();
            long long dE = objv - oold;
            if (dE <= 0 || urand() < exp(-(double)dE / T)) {
                if (fval < minF) { minF = fval; memcpy(mincode, code, sizeof code); }
                if (objv < bestobj) bestobj = objv;
            } else {
                for (int ii = 0; ii < R; ii++) code[ii][c] = oldcol[ii];
                memcpy(S, Sold, sizeof S); fval = fold; objv = oold;
            }
            T *= 0.9999; if (T < 0.05) T = 0.05;
        }
        restarts++;
    }
    printf("MIN %lld OBJ %lld RESTARTS %d STEPS %lld\n", minF, bestobj, restarts, steps);
    printf("CODE %d %d", R, M);
    for (int i = 0; i < R; i++) for (int c = 0; c < M; c++) printf(" %d", mincode[i][c]);
    printf("\n");
    return 0;
}
