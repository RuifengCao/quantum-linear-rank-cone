import sys, math, time
from fractions import Fraction
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'scripts'))
import s34_six_party_survey as s
t0 = time.time()
ok = True
print('d   2^f12 (exact)   gcd(d,2)/d   f12 (bits)')
for d in range(2, 31):
    q = s.seven_point_f12(d)
    cf = Fraction(math.gcd(d, 2), d)
    ok &= q == cf
    print(f'{d:2d}  {str(q):>12s}  {str(cf):>12s}  {math.log2(q):+.4f}')
print('Smith invariants other than 1:', s.seven_point_smith())
print('all equal to the closed form:', ok, f'({time.time() - t0:.0f} s)')
