/* Exhaustive check for Example 1 of the paper: no weighted graph state of six qutrits (one qutrit per system,
   weights in GF(3)) has the entropy vector q2.  By the normal form for prime local dimension and the local-Clifford
   equivalence of qudit stabilizer states with weighted graph states, this covers every stabilizer state of six
   qutrits with one qutrit per system.

   Input (stdin): the 31 entries of q2 in mask order (mask m over the five parties A..E, bit 0 = A).
   Output: the number of weighted graphs enumerated (3^15) and the number whose cut-rank vector over GF(3)
   equals the input vector.  The entropy of a set X of parties is rank_GF(3) W[X, complement of X]. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int rank3(int rows[6][6], int nr, int nc) {
    int M[6][6];
    for (int i = 0; i < nr; i++)
        for (int j = 0; j < nc; j++) M[i][j] = rows[i][j] % 3;
    int r = 0;
    for (int c = 0; c < nc && r < nr; c++) {
        int piv = -1;
        for (int i = r; i < nr; i++)
            if (M[i][c]) { piv = i; break; }
        if (piv < 0) continue;
        for (int j = 0; j < nc; j++) { int t = M[r][j]; M[r][j] = M[piv][j]; M[piv][j] = t; }
        int inv = (M[r][c] == 1) ? 1 : 2;           /* inverses in GF(3): 1 -> 1, 2 -> 2 */
        for (int j = 0; j < nc; j++) M[r][j] = (M[r][j] * inv) % 3;
        for (int i = 0; i < nr; i++)
            if (i != r && M[i][c]) {
                int f = M[i][c];
                for (int j = 0; j < nc; j++) M[i][j] = ((M[i][j] - f * M[r][j]) % 3 + 3) % 3;
            }
        r++;
    }
    return r;
}

int main(void) {
    int q2[32];
    for (int m = 1; m < 32; m++)
        if (scanf("%d", &q2[m]) != 1) { fprintf(stderr, "need 31 integers\n"); return 2; }
    int pa[15], pb[15], k = 0;
    for (int a = 0; a < 6; a++)
        for (int b = a + 1; b < 6; b++) { pa[k] = a; pb[k] = b; k++; }
    /* order the masks so that cheap tests (singletons, pairs) come first */
    int order[31], no = 0;
    for (int size = 1; size <= 5; size++)
        for (int m = 1; m < 32; m++)
            if (__builtin_popcount(m) == size) order[no++] = m;
    long hits = 0, total = 0;
    int W[6][6];
    long ncode = 1;
    for (int i = 0; i < 15; i++) ncode *= 3;
    for (long code = 0; code < ncode; code++) {
        long c = code;
        memset(W, 0, sizeof W);
        for (int e = 0; e < 15; e++) { int w = (int)(c % 3); c /= 3; W[pa[e]][pb[e]] = w; W[pb[e]][pa[e]] = w; }
        total++;
        int ok = 1;
        for (int t = 0; t < 31 && ok; t++) {
            int m = order[t];
            int X[6], nx = 0, Y[6], ny = 0;
            for (int v = 0; v < 6; v++) {
                if (v < 5 && (m >> v & 1)) X[nx++] = v; else Y[ny++] = v;
            }
            int rows[6][6];
            for (int i = 0; i < nx; i++)
                for (int j = 0; j < ny; j++) rows[i][j] = W[X[i]][Y[j]];
            if (rank3(rows, nx, ny) != q2[m]) ok = 0;
        }
        if (ok) hits++;
    }
    printf("%ld %ld\n", total, hits);
    return 0;
}
