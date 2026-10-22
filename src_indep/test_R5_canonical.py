#!/usr/bin/env python3
"""R5 with the CANONICAL merge definition supplied in the reconciliation
message: image = m2(u, qa+qb, v+t), i.e. the two-run word 1^u (qa+qb) 1^{v+t}.
This differs from BOTH variants I tested on 2026-10-22 (qa, and qa+qb-1).
T1 exact."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from indep_realizer import least, word, QS


def m2(u, q, v):
    return least((1,) * u + (q,) + (1,) * v)


print("1. the claimed R5 witness (u,qa,v,qb,t) = (22,4,2,2,6)")
u, qa, v, qb, t = 22, 4, 2, 2, 6
w = word(u, qa, v, qb, t); n = len(w); mw = least(w)
print(f"   |w| = {n}, bl(m(w)) = {mw.bit_length()}")
for lbl, img in (("canonical  m2(u, qa+qb, v+t)", (u, qa + qb, v + t)),
                 ("keep       m2(u, qa,    v+t)", (u, qa, v + t)),
                 ("additive   m2(u, qa+qb-1,v+t)", (u, qa + qb - 1, v + t))):
    mv = m2(*img)
    d = (mv.bit_length() - mw.bit_length()) / n
    print(f"   {lbl}: bl(m(v)) = {mv.bit_length():3d}  delta = {d:+.4f}"
          f"{'   <-- claimed 0.321' if abs(d-0.321)<0.005 else ''}")

import math
print("\n1b. same witness, TRUE log2 instead of bit length  [T2]")
for lbl, img in (("canonical", (u, qa+qb, v+t)), ("keep", (u, qa, v+t)),
                 ("additive", (u, qa+qb-1, v+t))):
    mv = m2(*img)
    d = (math.log2(mv) - math.log2(mw)) / n
    print(f"   {lbl:10s} delta = {d:+.4f}"
          f"{'   <-- MATCHES claimed 0.321' if abs(d-0.321)<0.001 else ''}")

print("\n2. sup delta, canonical definition, BOTH grids")
for gname, rng, step in (("small full grid (u,v,t <= 25)", range(0, 26), 1),
                         ("sizes >= 32 (sampled, stride 7)", range(0, 120, 7), 7)):
    best = None
    for U in rng:
        for V in rng:
            for T in rng:
                nn = U + V + T + 2
                if gname.startswith("sizes") and nn < 32:
                    continue
                if nn < 6:
                    continue
                for A in QS:
                    for B in QS:
                        ww = word(U, A, V, B, T)
                        mm = least(ww)
                        mv = m2(U, A + B, V + T)
                        d = (mv.bit_length() - mm.bit_length()) / nn
                        if best is None or d > best[0]:
                            best = (d, U, A, V, B, T, nn)
    print(f"   {gname}")
    print(f"      sup delta = {best[0]:+.4f} at (u,qa,v,qb,t)="
          f"({best[1]},{best[2]},{best[3]},{best[4]},{best[5]}), |w|={best[6]}")
