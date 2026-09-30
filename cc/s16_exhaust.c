/* s16_exhaust: exhaustive decision of exact realisability (fixed lambda) for small N.
 *
 * Same normal form as s16_anneal.c: party p owns n_p qubits (n_p = lambda*S_p),
 * only inter-party adjacency bits are free.  All 2^B configurations of the B
 * inter-party bits are enumerated in Gray-code order; each step flips one bit
 * and recomputes the 16 cut ranks that separate the two parties involved.
 * Output: FOUND + adjacency rows (first hit), or NONE after the full sweep --
 * a NONE is a complete (exhaustive) proof that lambda * r is not the entropy
 * vector of any qubit stabilizer state (given the normal-form lemma).
 *
 * Input (stdin): id lambda n0..n5 t1..t31   (one line per target)
 * Output: id lambda FOUND N rows... | id lambda NONE B | id lambda SKIP B
 * Build: gcc -O3 -march=native -fopenmp -o s16_exhaust s16_exhaust.c
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

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

typedef struct { int id, lam, n[6], t[32]; } Job;

int main(int argc, char **argv) {
    int maxB = argc > 1 ? atoi(argv[1]) : 34;
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
        int N = 0, off[6], party[64];
        for (int p = 0; p < 6; p++) { off[p] = N; for (int k = 0; k < J->n[p]; k++) party[N + k] = p; N += J->n[p]; }
        uint64_t all = (N >= 64) ? ~0ULL : ((1ULL << N) - 1), pm[6], qm[32], cm[32];
        for (int p = 0; p < 6; p++) pm[p] = (((1ULL << J->n[p]) - 1) << off[p]);
        for (int m = 0; m < 32; m++) { qm[m] = 0; for (int p = 0; p < 5; p++) if (m >> p & 1) qm[m] |= pm[p]; cm[m] = all & ~qm[m]; }
        int B = 0, pi[2048], pj[2048];
        for (int i = 0; i < N; i++) for (int j = i + 1; j < N; j++) if (party[i] != party[j]) { pi[B] = i; pj[B] = j; B++; }
        char line[8192]; int len;
        if (B > maxB) {
            snprintf(line, sizeof line, "%d %d SKIP %d\n", J->id, J->lam, B);
            #pragma omp critical
            { fputs(line, stdout); fflush(stdout); }
            continue;
        }
        uint64_t adj[64]; memset(adj, 0, sizeof adj);
        int rk[32], bad = 0;
        for (int m = 1; m < 32; m++) { rk[m] = rank_rows(adj, qm[m], cm[m]); bad += (rk[m] != J->t[m]); }
        int found = (bad == 0);
        uint64_t total = (B == 64) ? ~0ULL : (1ULL << B);
        for (uint64_t g = 1; g < total && !found; g++) {
            int k = __builtin_ctzll(g);           /* Gray code: bit k flips at step g */
            int i = pi[k], j = pj[k], p = party[i], q = party[j];
            adj[i] ^= 1ULL << j; adj[j] ^= 1ULL << i;
            for (int s = 0; s < nsep[p][q]; s++) {
                int m = sep[p][q][s];
                int nr = rank_rows(adj, qm[m], cm[m]);
                bad += (nr != J->t[m]) - (rk[m] != J->t[m]);
                rk[m] = nr;
            }
            if (bad == 0) found = 1;
        }
        if (found) {
            len = snprintf(line, sizeof line, "%d %d FOUND %d", J->id, J->lam, N);
            for (int i = 0; i < N; i++) len += snprintf(line + len, sizeof line - len, " %llx", (unsigned long long)adj[i]);
            snprintf(line + len, sizeof line - len, "\n");
        } else snprintf(line, sizeof line, "%d %d NONE %d\n", J->id, J->lam, B);
        #pragma omp critical
        { fputs(line, stdout); fflush(stdout); }
    }
    return 0;
}
