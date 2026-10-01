/* s24: exhaustive check for one qudit per party (R31).
 * Enumerates every symmetric 6 x 6 matrix with zero diagonal over GF(q), q prime (q^15 matrices), i.e.
 * every weighted graph state with one qudit for each of the parties A..E and the purifier F, and counts
 * those whose cut ranks rank W[X, X^c] (X a nonempty subset of A..E) equal the 31 targets.  For a fixed
 * multiplier lambda with lambda * S_p = 1 for every party this covers all stabilizer states (normal form
 * and graph-state lemmas, docs/A2_stabilizer_vs_QLR.md), so "hits 0" proves that lambda * r is not the
 * entropy vector of any stabilizer state of local dimension q with one qudit per party.
 * Input (stdin): t1..t31.  Argument: q (2, 3, 5 or 7; 5 and 7 take very long).
 * Output: "q=Q total T hits H best_nonzero_residual E" (+ up to three hits as 15 upper-triangle entries).
 * Build: gcc -O3 -march=native -o s24_exhaustive6 s24_exhaustive6.c
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
static int q, INV[16];
static int rk(int G[6][6], int X) {
    int ri[6], ci[6], nr = 0, nc = 0;
    for (int v = 0; v < 6; v++) { if (v < 5 && (X >> v & 1)) ri[nr++] = v; else ci[nc++] = v; }
    int M[6][6];
    for (int a = 0; a < nr; a++) for (int b = 0; b < nc; b++) M[a][b] = G[ri[a]][ci[b]];
    int r = 0;
    for (int c = 0; c < nc && r < nr; c++) {
        int p = -1; for (int a = r; a < nr; a++) if (M[a][c]) { p = a; break; }
        if (p < 0) continue;
        for (int b = 0; b < nc; b++) { int t = M[r][b]; M[r][b] = M[p][b]; M[p][b] = t; }
        int iv = INV[M[r][c]]; for (int b = 0; b < nc; b++) M[r][b] = (M[r][b] * iv) % q;
        for (int a = 0; a < nr; a++) if (a != r && M[a][c]) { int f = M[a][c]; for (int b = 0; b < nc; b++) M[a][b] = ((M[a][b] - f * M[r][b]) % q + q) % q; }
        r++;
    }
    return r;
}
int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: %s q < targets\n", argv[0]); return 2; }
    q = atoi(argv[1]); int t[32];
    if (q != 2 && q != 3 && q != 5 && q != 7) { fprintf(stderr, "q must be a prime <= 7\n"); return 2; }
    for (int a = 1; a < q; a++) for (int b = 1; b < q; b++) if (a * b % q == 1) INV[a] = b;
    for (int m = 1; m < 32; m++) if (scanf("%d", &t[m]) != 1) return 3;
    int pairs[15][2], np = 0;
    for (int i = 0; i < 6; i++) for (int j = i + 1; j < 6; j++) { pairs[np][0] = i; pairs[np][1] = j; np++; }
    long long total = 1; for (int k = 0; k < 15; k++) total *= q;
    long long hits = 0, best = 99; int G[6][6];
    int order[31], no = 0;            /* small cuts first: cheap rejection */
    for (int s = 1; s <= 5; s++) for (int m = 1; m < 32; m++) if (__builtin_popcount(m) == s) order[no++] = m;
    for (long long code = 0; code < total; code++) {
        long long c = code; memset(G, 0, sizeof G);
        for (int k = 0; k < 15; k++) { int v = (int)(c % q); c /= q; G[pairs[k][0]][pairs[k][1]] = G[pairs[k][1]][pairs[k][0]] = v; }
        int ok = 1; long long E = 0;
        for (int k = 0; k < 31; k++) { int m = order[k]; int d = rk(G, m) - t[m]; if (d) { ok = 0; E += abs(d); if (E >= best) break; } }
        if (ok) { hits++; if (hits <= 3) { printf("HIT"); for (int k = 0; k < 15; k++) printf(" %d", G[pairs[k][0]][pairs[k][1]]); printf("\n"); } }
        else if (E < best) best = E;
    }
    printf("q=%d total %lld hits %lld best_nonzero_residual %lld\n", q, total, hits, best);
    return 0;
}
