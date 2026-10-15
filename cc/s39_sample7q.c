/* s39: one qudit per system over other fields (R45; docs/R32_reduction_and_local_dimension.md, section 4.10, "Other
 * fields").  Random weighted graphs over GF(q) on seven vertices (0..5 visible, 6 the purifier), q = p^k <= 64 with
 * table arithmetic (q in 2 3 4 5 7 8 9 11 13 16 25 27 49); the cut-rank vector (63 cuts, bit i of m = system i) of each
 * graph is the entropy vector, in units of log q, of a GF(q)-linear stabilizer state with one qudit per system.
 * Every distinct vector is written once to <out>.bin (63 bytes) with the 21 weights of one graph realising it to
 * <out>.w (edges (0,1), (0,2), ..., (5,6)).  Weights: uniform nonzero, all 1, all +-1, or drawn from a small random
 * set of nonzero elements; edge density uniform in [0.05, 1].
 * Field elements are integers 0..q-1 whose base-p digits are the coefficients of 1, x, x^2, ... modulo a fixed monic
 * irreducible: GF(4) x^2+x+1, GF(8) x^3+x+1, GF(16) x^4+x+1, GF(9) x^2+1, GF(27) x^3+2x+2, GF(25) x^2-2, GF(49) x^2+1.
 * The sampler proves nothing about completeness; scripts/s36_one_qudit_fields.py re-checks every stored certificate
 * with a separate implementation of the fields.
 * Usage: s39_sample7q q samples seed out_prefix
 * Build: gcc -O3 -march=native -o s39_sample7q s39_sample7q.c
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

static int Q, PCH, DEG;
static uint8_t ADD[64][64], SUB[64][64], MUL[64][64], INV[64];
static uint64_t rs;
static inline uint64_t rnd(void) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; }

/* elements: integers 0..q-1, base-p digits = coefficients of 1, x, x^2, ...; modulus x^k = -(c_{k-1} x^{k-1} + ... + c_0) */
static int modpoly[4];   /* c_0..c_{k-1} of the monic irreducible x^k + c_{k-1} x^{k-1} + ... + c_0 */
static void digits(int a, int *d) { for (int i = 0; i < DEG; i++) { d[i] = a % PCH; a /= PCH; } }
static int undigits(const int *d) { int a = 0; for (int i = DEG - 1; i >= 0; i--) a = a * PCH + d[i]; return a; }
static int fmul(int a, int b) {
    int da[4], db[4], prod[8] = {0};
    digits(a, da); digits(b, db);
    for (int i = 0; i < DEG; i++) for (int j = 0; j < DEG; j++) prod[i + j] = (prod[i + j] + da[i] * db[j]) % PCH;
    for (int t = 2 * DEG - 2; t >= DEG; t--) {          /* x^t = x^{t-k} x^k = -x^{t-k} sum c_i x^i */
        int c = prod[t];
        if (!c) continue;
        prod[t] = 0;
        for (int i = 0; i < DEG; i++) prod[t - DEG + i] = ((prod[t - DEG + i] - c * modpoly[i]) % PCH + PCH) % PCH;
    }
    return undigits(prod);
}
static int setup(int q) {
    Q = q;
    int ps[] = {2, 3, 5, 7, 11, 13};
    PCH = 0;
    for (int i = 0; i < 6; i++) { int p = ps[i], x = q, k = 0; while (x % p == 0) { x /= p; k++; } if (x == 1 && k) { PCH = p; DEG = k; } }
    if (!PCH) return 0;
    memset(modpoly, 0, sizeof modpoly);
    if (DEG == 1) modpoly[0] = 0;                        /* x (unused) */
    else if (q == 4) { modpoly[0] = 1; modpoly[1] = 1; }  /* x^2 + x + 1 */
    else if (q == 8) { modpoly[0] = 1; modpoly[1] = 1; }  /* x^3 + x + 1 */
    else if (q == 16) { modpoly[0] = 1; modpoly[1] = 1; } /* x^4 + x + 1 */
    else if (q == 9) { modpoly[0] = 1; }                  /* x^2 + 1 */
    else if (q == 27) { modpoly[0] = 2; modpoly[1] = 2; } /* x^3 + 2x + 2 */
    else if (q == 25) { modpoly[0] = 3; }                 /* x^2 + 3 = x^2 - 2 */
    else if (q == 49) { modpoly[0] = 1; }                 /* x^2 + 1 */
    else return 0;
    for (int a = 0; a < q; a++) for (int b = 0; b < q; b++) {
        int da[4], db[4], s[4], d[4];
        digits(a, da); digits(b, db);
        for (int i = 0; i < DEG; i++) { s[i] = (da[i] + db[i]) % PCH; d[i] = (da[i] - db[i] + PCH) % PCH; }
        ADD[a][b] = undigits(s); SUB[a][b] = undigits(d);
        MUL[a][b] = DEG == 1 ? a * b % PCH : fmul(a, b);
    }
    for (int a = 1; a < q; a++) for (int b = 1; b < q; b++) if (MUL[a][b] == 1) INV[a] = b;
    for (int a = 1; a < q; a++) if (!INV[a]) return 0;    /* not a field: wrong modulus */
    return 1;
}

static int rank_mat(int R, int C, uint8_t M[4][7]) {
    int r = 0;
    for (int c = 0; c < C && r < R; c++) {
        int piv = -1;
        for (int i = r; i < R; i++) if (M[i][c]) { piv = i; break; }
        if (piv < 0) continue;
        if (piv != r) for (int j = 0; j < C; j++) { uint8_t t = M[r][j]; M[r][j] = M[piv][j]; M[piv][j] = t; }
        int iv = INV[M[r][c]];
        for (int i = 0; i < R; i++) if (i != r && M[i][c]) {
            int f = MUL[M[i][c]][iv];
            for (int j = 0; j < C; j++) M[i][j] = SUB[M[i][j]][MUL[f][M[r][j]]];
        }
        r++;
    }
    return r;
}

#define HBITS 20
#define HSIZE (1u << HBITS)
typedef struct { uint8_t v[63]; uint8_t w[21]; uint8_t used; } entry;
static entry *seen;
static uint64_t hashv(const uint8_t *v) {
    uint64_t h = 1469598103934665603ULL;
    for (int i = 0; i < 63; i++) { h ^= v[i]; h *= 1099511628211ULL; }
    return h;
}

int main(int argc, char **argv) {
    if (argc < 5) { fprintf(stderr, "usage: s39_sample7q q samples seed out_prefix\n"); return 1; }
    if (!setup(atoi(argv[1]))) { fprintf(stderr, "unsupported q\n"); return 2; }
    long long N = atoll(argv[2]);
    rs = strtoull(argv[3], 0, 10) * 2654435761ULL + 88172645463325252ULL;
    seen = calloc(HSIZE, sizeof(entry));
    long long distinct = 0;
    uint8_t A[7][7];
    for (long long s = 0; s < N; s++) {
        double dens = 0.05 + 0.95 * (double)(rnd() % 1000000) / 1e6;
        int mode = (int)(rnd() % 4);
        int nvals = 2 + (int)(rnd() % (uint64_t)(Q - 1));     /* mode 3: weights from a random small set */
        uint8_t vals[64];
        for (int i = 0; i < nvals; i++) vals[i] = 1 + (uint8_t)(rnd() % (uint64_t)(Q - 1));
        memset(A, 0, sizeof A);
        for (int u = 0; u < 7; u++) for (int v = u + 1; v < 7; v++) {
            if ((double)(rnd() % 1000000) / 1e6 >= dens) continue;
            uint8_t w;
            if (mode == 1) w = 1;
            else if (mode == 2) w = (rnd() & 1) ? 1 : SUB[0][1];
            else if (mode == 3) w = vals[rnd() % (uint64_t)nvals];
            else w = 1 + (uint8_t)(rnd() % (uint64_t)(Q - 1));
            A[u][v] = A[v][u] = w;
        }
        uint8_t S[63];
        for (int m = 1; m < 64; m++) {
            int X[7], Y[7], nx = 0, ny = 0;
            for (int v = 0; v < 7; v++) { if (v < 6 && (m >> v & 1)) X[nx++] = v; else Y[ny++] = v; }
            int *rows = nx <= ny ? X : Y, *cols = nx <= ny ? Y : X;
            int R = nx <= ny ? nx : ny, C = nx <= ny ? ny : nx;
            uint8_t M[4][7];
            for (int i = 0; i < R; i++) for (int j = 0; j < C; j++) M[i][j] = A[rows[i]][cols[j]];
            S[m - 1] = (uint8_t)rank_mat(R, C, M);
        }
        uint64_t h = hashv(S) & (HSIZE - 1);
        while (seen[h].used && memcmp(seen[h].v, S, 63)) h = (h + 1) & (HSIZE - 1);
        if (!seen[h].used) {
            seen[h].used = 1; memcpy(seen[h].v, S, 63);
            int k = 0;
            for (int u = 0; u < 7; u++) for (int v = u + 1; v < 7; v++) seen[h].w[k++] = A[u][v];
            distinct++;
        }
    }
    char fn[512];
    snprintf(fn, sizeof fn, "%s.bin", argv[4]); FILE *fv = fopen(fn, "wb");
    snprintf(fn, sizeof fn, "%s.w", argv[4]); FILE *fw = fopen(fn, "wb");
    for (uint64_t h = 0; h < HSIZE; h++) if (seen[h].used) { fwrite(seen[h].v, 1, 63, fv); fwrite(seen[h].w, 1, 21, fw); }
    fclose(fv); fclose(fw);
    printf("q %d SAMPLES %lld DISTINCT %lld\n", Q, N, distinct);
    return 0;
}
