/* s7 tool A: enumerate ALL F_3-weighted graphs on 6 vertices (3^15 = 14,348,907),
 * compute the 31 purified cut-rank entropy coordinates, write raw int8 vectors.
 * Qudit (qutrit) graph states: S(S) = rank_{F_3} A[S, S^c]  (cut-rank, purity built in).
 * Output layer is CONDITIONAL on the ledger's stabilizer definition (see README).
 * Build: gcc -O3 -fopenmp -o fp_cutrank_enum fp_cutrank_enum.c
 * Run:   ./fp_cutrank_enum out_prefix    -> out_prefix_t{tid}.bin  (31 bytes per vector)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif
static const int P = 3;
static int rank_f3(int A[6][6], int rmask, int cmask) {
    int R[6][6], nr = 0, nc = 0, rows[6], cols[6];
    for (int i = 0; i < 6; i++) if (rmask >> i & 1) rows[nr++] = i;
    for (int j = 0; j < 6; j++) if (cmask >> j & 1) cols[nc++] = j;
    for (int i = 0; i < nr; i++) for (int j = 0; j < nc; j++) R[i][j] = A[rows[i]][cols[j]];
    int r = 0;
    for (int c = 0; c < nc && r < nr; c++) {
        int piv = -1;
        for (int i = r; i < nr; i++) if (R[i][c]) { piv = i; break; }
        if (piv < 0) continue;
        for (int j = 0; j < nc; j++) { int t = R[r][j]; R[r][j] = R[piv][j]; R[piv][j] = t; }
        int inv = (R[r][c] == 1) ? 1 : 2;               /* inverses in F_3 */
        for (int j = 0; j < nc; j++) R[r][j] = R[r][j] * inv % P;
        for (int i = 0; i < nr; i++) if (i != r && R[i][c]) {
            int f = R[i][c];
            for (int j = 0; j < nc; j++) R[i][j] = (R[i][j] - f * R[r][j] % P + 3 * P) % P;
        }
        r++;
    }
    return r;
}
int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: %s out_prefix [limit]\n", argv[0]); return 1; }
    long total = 1; for (int i = 0; i < 15; i++) total *= 3;   /* 14,348,907 */
    if (argc > 2) { long lim = atol(argv[2]); if (lim > 0 && lim < total) total = lim; }
    long done = 0;
    int reps[31], nrep = 0;                                    /* masks over 6 without bit5 */
    for (int m = 1; m < 63; m++) if (!(m >> 5 & 1)) reps[nrep++] = m;
#pragma omp parallel
    {
#ifdef _OPENMP
        int tid = omp_get_thread_num();
#else
        int tid = 0;
#endif
        char name[512]; snprintf(name, sizeof name, "%s_t%03d.bin", argv[1], tid);
        FILE *f = fopen(name, "wb");
        unsigned char buf[31];
#pragma omp for schedule(dynamic, 4096)
        for (long code = 0; code < total; code++) {
            int A[6][6]; memset(A, 0, sizeof A);
            long x = code; int e = 0;
            for (int i = 0; i < 6; i++) for (int j = i + 1; j < 6; j++, e++) {
                int w = x % 3; x /= 3;
                A[i][j] = A[j][i] = w;
            }
            for (int k = 0; k < nrep; k++) {
                int m = reps[k];
                buf[m - 1] = (unsigned char) rank_f3(A, m, 63 ^ m);
            }
            fwrite(buf, 1, 31, f);
            long d;
#ifdef _OPENMP
#pragma omp atomic capture
            d = ++done;
#else
            d = ++done;
#endif
            if ((d & ((1L << 20) - 1)) == 0)
                fprintf(stderr, "fp_cutrank: %ld/%ld (%.0f%%)\n", d, total, 100.0 * d / total);
        }
        fclose(f);
    }
    return 0;
}
