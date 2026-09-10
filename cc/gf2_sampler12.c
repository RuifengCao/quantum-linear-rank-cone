/* s7 tool B: sample GF(2) graph states on 12 vertices, parties = 6 blocks of 2
 * (5 visible + purifier block).  S(mask) = rank_{F_2} of the block cut submatrix.
 * Build: gcc -O3 -fopenmp -o gf2_sampler12 gf2_sampler12.c
 * Run:   ./gf2_sampler12 N seed out_prefix
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
static int rank_rows(uint16_t *rows, int nr) {
    int r = 0;
    for (int c = 0; c < 12 && r < nr; c++) {
        int piv = -1;
        for (int i = r; i < nr; i++) if (rows[i] >> c & 1) { piv = i; break; }
        if (piv < 0) continue;
        uint16_t t = rows[r]; rows[r] = rows[piv]; rows[piv] = t;
        for (int i = 0; i < nr; i++) if (i != r && (rows[i] >> c & 1)) rows[i] ^= rows[r];
        r++;
    }
    return r;
}
int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage: %s N seed out_prefix\n", argv[0]); return 1; }
    long N = atol(argv[1]);
    uint64_t s = strtoull(argv[2], 0, 10);
    char name[512]; snprintf(name, sizeof name, "%s.bin", argv[3]);
    FILE *f = fopen(name, "wb");
    for (long it = 0; it < N; it++) {
        uint16_t adj[12] = {0};
        for (int i = 0; i < 12; i++) for (int j = i + 1; j < 12; j++) {
            s ^= s << 13; s ^= s >> 7; s ^= s << 17;          /* xorshift */
            if (s & 1) { adj[i] |= (uint16_t)1 << j; adj[j] |= (uint16_t)1 << i; }
        }
        unsigned char buf[31];
        for (int m = 1; m < 63; m++) {
            if (m >> 5 & 1) continue;
            uint16_t vmask = 0, sub[12]; int nr = 0;
            for (int b = 0; b < 6; b++) if (m >> b & 1) vmask |= (uint16_t)3 << (2 * b);
            uint16_t cmask = (uint16_t)((~vmask) & 0xFFF);
            for (int v = 0; v < 12; v++) if (vmask >> v & 1) sub[nr++] = adj[v] & cmask;
            buf[m - 1] = (unsigned char) rank_rows(sub, nr);
        }
        fwrite(buf, 1, 31, f);
        if ((it + 1) % 1000000 == 0)
            fprintf(stderr, "gf2_sampler: %ld/%ld\n", it + 1, N);
    }
    fclose(f);
    return 0;
}
