#!/usr/bin/env python3
"""PHASE A1 -- test R4 (beta* = 1/2) and R5 (bridge witnesses)."""
import math, os, sys
from fractions import Fraction as F
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from indep_realizer import least, word, F_c, run_floor, QS

print("R4  uniform reduction: thresholds on the imported one-defect beta")
print("    balanced u=v=t=k, qa=qb=2.  P = 1^k 2 1^k, S_P = 2k+2.")
print("    run-based null (leading run) = k+1 bits.\n")
print(f"    {'k':>5s}{'trunc thresh':>16s}{'adv thresh':>14s}"
      f"{'(k+1)/(2k+2)':>15s}")
for k in (1, 2, 5, 10, 50, 500, 5000):
    # truncation branch: beta*(2k+2) > k+1  <=>  beta > (k+1)/(2k+2)
    tr = F(k + 1, 2 * k + 2)
    # advance branch: (k+2) + beta(2k+2) - alpha(k+1) > k+1
    #              <=> beta > (alpha(k+1) - 1)/(2k+2)
    al = math.log2(3)
    ad = (al * (k + 1) - 1) / (2 * k + 2)
    print(f"    {k:5d}{float(tr):16.6f}{ad:14.6f}{float(tr):15.6f}")
print(f"\n    truncation threshold -> 1/2 EXACTLY (it is (k+1)/(2k+2) = 1/2 for")
print(f"    every k, not just in the limit).  advance threshold -> alpha/2 =")
print(f"    {math.log2(3)/2:.6f}  (the brief's '~0.78').")
print("    => R4's beta* = 1/2 CONFIRMED, and it is exact at every k.")

# empirical beta for one-defect words
print("\n    empirical beta for one-defect words 1^k 2 1^k:")
for k in (5, 20, 60, 150):
    P = (1,) * k + (2,) + (1,) * k
    m = least(P); S = sum(P)
    print(f"      k={k:4d}  S_P={S:5d}  bl(m(P))={m.bit_length():5d}"
          f"  beta_emp={m.bit_length()/S:.4f}")

print("\n" + "=" * 70)
print("R5  bridge witnesses")
print("=" * 70)
print("    claimed: merge sup delta = 0.321 at (u,qa,v,qb,t)=(22,4,2,2,6), |w|=32")
u, qa, v, qb, t = 22, 4, 2, 2, 6
w = word(u, qa, v, qb, t)
n = len(w)
print(f"    |w| = {n}  (claimed 32)  -> {'MATCH' if n == 32 else 'MISMATCH'}")
mw = least(w)
# merge images, as defined in my package
for lbl, img in (("merge keep qa",   (1,)*u + (qa,) + (1,)*(v+t)),
                 ("merge additive",  (1,)*u + (qa+qb-1,) + (1,)*(v+t))):
    mv = least(img)
    d = (mv.bit_length() - mw.bit_length()) / n
    print(f"    {lbl:16s} bl(m(w))={mw.bit_length():3d} "
          f"bl(m(v))={mv.bit_length():3d}  delta={d:+.4f}")

print("\n    sup over the whole grid of (bl(m(v))-bl(m(w)))/|w|, merge keep:")
best = None
for U in range(26):
    for V in range(26):
        for T in range(26):
            for A in QS:
                for B in QS:
                    ww = word(U, A, V, B, T)
                    nn = len(ww)
                    if nn < 6:
                        continue
                    mm = least(ww)
                    im = (1,)*U + (A,) + (1,)*(V+T)
                    d = (least(im).bit_length() - mm.bit_length()) / nn
                    if best is None or d > best[0]:
                        best = (d, U, A, V, B, T, nn)
print(f"      sup delta = {best[0]:+.4f} at (u,qa,v,qb,t)="
      f"({best[1]},{best[2]},{best[3]},{best[4]},{best[5]}), |w|={best[6]}")
print(f"      claimed 0.321 -> {'reproduced' if abs(best[0]-0.321)<0.02 else 'NOT reproduced'}")
