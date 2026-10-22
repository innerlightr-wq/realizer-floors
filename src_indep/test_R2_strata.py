#!/usr/bin/env python3
"""R2 in the corrected form: stratify gcd(2^S - C_w, 3^n) by which of
u,v,t vanish.  Tests the claim '3-coprime exactly on u,v,t >= 1' and
'structural reduction only at t=0'.  T1 exact."""
import os, sys
from math import gcd
from collections import Counter
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from indep_realizer import word, C_affine, least, QS

def v3(x):
    k = 0
    while x and x % 3 == 0:
        x //= 3; k += 1
    return k

strata = {}
R = 14
for u in range(R):
    for v in range(R):
        for t in range(R):
            for qa in QS:
                for qb in QS:
                    w = word(u, qa, v, qb, t)
                    n, S = len(w), sum(w)
                    N = (1 << S) - C_affine(w)
                    g = v3(N) if N else -1
                    key = (u >= 1, v >= 1, t >= 1)
                    s = strata.setdefault(key, Counter())
                    s[g] += 1

print("stratum (u>=1, v>=1, t>=1) -> distribution of v3(2^S - C_w)")
print(f"{'u>=1':>6s}{'v>=1':>6s}{'t>=1':>6s}{'n':>9s}   v3 histogram")
for key in sorted(strata, reverse=True):
    c = strata[key]
    tot = sum(c.values())
    hist = " ".join(f"{k}:{c[k]}" for k in sorted(c))
    print(f"{str(key[0]):>6s}{str(key[1]):>6s}{str(key[2]):>6s}{tot:>9d}   {hist}")

print()
interior = strata[(True, True, True)]
print(f"INTERIOR u,v,t >= 1 : {sum(interior.values())} configs, "
      f"v3 values present = {sorted(interior)}")
print(f"  -> '3-coprime exactly on u,v,t >= 1' : "
      f"{'CONFIRMED' if set(interior)=={0} else 'REFUTED'}")

# is the reduction confined to t=0?
t0 = sum(sum(v for k, v in strata[key].items() if k > 0)
         for key in strata if not key[2])
tpos = sum(sum(v for k, v in strata[key].items() if k > 0)
           for key in strata if key[2])
print(f"  reductions (v3>0) with t=0  : {t0}")
print(f"  reductions (v3>0) with t>=1 : {tpos}")
print(f"  -> 'reduction only at t=0' : {'CONFIRMED' if tpos==0 else 'REFUTED'}")

# the x=1 cycle word witness
print("\nclaimed witness = the x=1 cycle word.  x=1 under T: 1 -> (3+1)/4 = 1,")
print("so the cycle word is the single letter (2); as a two-defect word with")
print("u=v=t=0 that is qa=qb=2 -> w = (2,2):")
w = (2, 2)
N = (1 << sum(w)) - C_affine(w)
print(f"   w={w}  n={len(w)}  S={sum(w)}  C_w={C_affine(w)}  "
      f"2^S-C_w={N}  v3={v3(N)}  3^n={3**len(w)}")
print(f"   gcd = {gcd(N, 3**len(w))}, so D = 3^n/gcd = {3**len(w)//gcd(N,3**len(w))}")
print(f"   least realizer m(w) = {least(w)}  -> orbit: ", end="")
x = least(w)
for _ in range(2):
    y = 3*x+1; k=0
    while y%2==0: y//=2; k+=1
    print(f"{x} -", end=""); x=y
print(f"> {x}   (x=1 fixed point reached: {x==1})")
