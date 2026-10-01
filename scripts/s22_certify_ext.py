#!/usr/bin/env python3
"""s22: exact certificates for CI-extension infeasibilities with extra rows on the extension (R29).

For every ray that s17_ci_test.py (run with --extension ingleton) reported infeasible at a pair
(X, Y), compute a Farkas certificate over Shannon elemental inequalities AND conditional Ingleton
instances on the 7-element extension A..F,Z (s17_certify.certify(..., extension='ingleton')),
translate it into an inequality F'(S) >= 0 on pure-state entropy vectors, and re-verify it with the
independent checker epr1kit.stabcert.verify_exclusion_ext, which rebuilds every Ingleton row from
its descriptor ('I', a, b, c, d, K).  Valid for every stabilizer state of every prime dimension:
the extension by Z = A_X cap A_Y is a subspace arrangement, and subspace arrangements satisfy
Shannon and (conditional) Ingleton inequalities over any field.

  python3 scripts/s22_certify_ext.py --rays rays.npy --ci ci.json --extension ingleton --out certs.json
"""
import argparse, json, os, sys, time
from fractions import Fraction
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..'))
from s17_certify import certify, s_form, extension_rows
from epr1kit import stabcert


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--rays', required=True)
    ap.add_argument('--ci', required=True, help='JSON written by s17_ci_test.py --out')
    ap.add_argument('--extension', default='ingleton', choices=['shannon', 'ingleton'])
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    R = np.load(a.rays).astype(np.int64)
    _, desc = extension_rows(a.extension)
    certs, failed = [], []
    t0 = time.time()
    for rec in json.load(open(a.ci)):
        if not rec['infeasible_pairs']:
            continue
        i = int(rec['index']); X, Y = (int(t) for t in rec['infeasible_pairs'][0])
        out = certify(R[i], X, Y, extension=a.extension)
        if out is None:
            failed.append(i); print(f'ray {i}: ({X},{Y}) -- no exact certificate'); continue
        F_int, y_int, val = out
        F_S = [int(c) for c in s_form(F_int)]
        terms = [(desc[k], Fraction(v)) for k, v in sorted(y_int.items())]
        ok, value = stabcert.verify_exclusion_ext(R[i], X, Y, terms, F_S=F_S)
        if not ok:
            failed.append(i); print(f'ray {i}: ({X},{Y}) -- certificate REJECTED by the independent checker'); continue
        n_ing = sum(1 for d, _ in terms if d[0] == 'I')
        certs.append({'index': i, 'ray': [int(x) for x in R[i]], 'X': X, 'Y': Y, 'extension': a.extension,
                      'F_S': F_S, "F_S_at_ray": int(sum(c * x for c, x in zip(F_S, R[i]))),
                      'terms': [[list(d), f'{v.numerator}/{v.denominator}'] for d, v in terms]})
        print(f'ray {i}: ({X},{Y}) certified; {len(terms)} terms ({n_ing} Ingleton), F\'(r) = {certs[-1]["F_S_at_ray"]}', flush=True)
    json.dump({'extension': a.extension, 'certificates': certs, 'failed': failed}, open(a.out, 'w'), indent=1)
    print(f'{len(certs)} certified, {len(failed)} failed ({time.time() - t0:.0f}s) -> {a.out}')


if __name__ == '__main__':
    main()
