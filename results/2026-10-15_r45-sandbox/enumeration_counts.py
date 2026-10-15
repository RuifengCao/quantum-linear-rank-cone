"""Count the prefixes and seven-vertex graphs that cc/s40_enum7q.c examines, from support patterns only."""
import itertools, sys

def counts(q, reduced=True):
    nf = lambda s: 1 if s == 0 else (q - 1) ** (bin(s).count('1') - 1)   # normal forms with support s
    pre = gra = 0
    for s1, s2, s3, s4, s5 in itertools.product(range(2), range(4), range(8), range(16), range(32)):
        S = [0, s1, s2, s3, s4, s5]                     # support of the weights of vertex k to 0..k-1
        deg = [0] * 6
        for k in range(1, 6):
            for j in range(k):
                if S[k] >> j & 1:
                    deg[k] += 1; deg[j] += 1
        if reduced and any(deg[i] < deg[i + 1] for i in range(5)):
            continue
        w = 1
        for k in range(1, 6):
            w *= nf(S[k])
        pre += w
        dmax = max(deg)
        for s6 in range(64):
            z = bin(s6).count('1')
            if reduced and (z < dmax or any(s6 >> j & 1 and deg[j] >= z for j in range(6))):
                continue
            gra += w * nf(s6)
    return pre, gra

for q in map(int, sys.argv[1:]):
    print(q, 'reduced', counts(q), 'all', counts(q, False))
