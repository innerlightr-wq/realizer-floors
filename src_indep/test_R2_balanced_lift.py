#!/usr/bin/env python3
"""R2 corrected form: the Lemma-H denominator is REPRESENTATIVE-DEPENDENT.
Canonical datum = the balanced lift N* = -C (mod 2^S), |N*| <= 2^{S-1},
giving D = 3^{n - v3(N*)}.  Verified independently (src_indep route).  T1,
with [T2] medians.  Family: u,v,t in 1..12, qa,qb in {2,3,4} -> 15552."""
import math, os, statistics as st, sys
from collections import Counter
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from indep_realizer import word, C_affine, least, QS

def v3(x):
    k = 0
    while x and x % 3 == 0:
        x //= 3; k += 1
    return k

cnt = Counter(); n = 0
one_sided_a = one_sided_b = 0
med0, med1 = [], []
for u in range(1, 13):
    for v in range(1, 13):
        for t in range(1, 13):
            for qa in QS:
                for qb in QS:
                    w = word(u, qa, v, qb, t)
                    S = sum(w); C = C_affine(w); nn = len(w)
                    # balanced lift of -C mod 2^S into (-2^{S-1}, 2^{S-1}]
                    Ns = (-C) % (1 << S)
                    if Ns > (1 << (S - 1)):
                        Ns -= (1 << S)
                    k = v3(abs(Ns))
                    cnt[k] += 1; n += 1
                    one_sided_a += (v3(C) > 0)                 # -C/3^n rep
                    one_sided_b += (v3((1 << S) - C) > 0)      # (2^S-C)/3^n rep
                    share = math.log2(least(w)) / (S + 1)
                    (med0 if k == 0 else med1).append(share)

print(f"family size                              : {n}   (claimed 15552)")
print(f"one-sided rep -C  : v3>0 count           : {one_sided_a}   (claimed 0)")
print(f"one-sided rep 2^S-C : v3>0 count         : {one_sided_b}   (claimed 0)")
print(f"balanced lift N*  : v3>0 count           : {n-cnt[0]}   (claimed 5263)")
print(f"                    fraction             : {(n-cnt[0])/n:.4f}   (claimed ~1/3)")
print()
print("v3(N*) distribution vs geometric(2/3) prediction:")
print(f"  {'k':>3s}{'observed':>10s}{'claimed':>9s}{'geom(2/3)':>11s}")
claimed = {0:10289,1:3478,2:1199,3:417,4:123,5:32,6:9,7:3,8:1,9:1}
for k in range(0, max(cnt)+1):
    pred = n * (2/3) * (1/3)**k
    print(f"  {k:3d}{cnt[k]:10d}{claimed.get(k,0):9d}{pred:11.1f}")
print(f"  max v3 on grid = {max(cnt)}   (claimed 9)")
print()
print(f"dangerous-window probe, median log2 m / (S+1)  [T2]:")
print(f"  v3(N*) = 0 : {st.median(med0):.4f}   (claimed 0.9615)   n={len(med0)}")
print(f"  v3(N*) > 0 : {st.median(med1):.4f}   (claimed 0.9620)   n={len(med1)}")
print(f"  difference : {abs(st.median(med0)-st.median(med1)):.5f}"
      f"  -> denominator reduction is uncorrelated with realizer smallness")
