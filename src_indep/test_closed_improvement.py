#!/usr/bin/env python3
"""Does the closed-episode prefix floor strengthen Theorem 5?  T1 exact."""
import math, os, statistics as st, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from indep_realizer import least, klass, word, C_affine, QS
from test_closed_episode import m_closed

def branches(u, qa, v, qb, t):
    P  = (1,)*u + (qa,) + (1,)*v
    Pp = (1,)*u + (qa,)
    Q  = (1,)*v + (qb,) + (1,)*t
    f_open   = least(P)
    f_closed = m_closed(P)
    Sp, Cp = u + qa, C_affine(Pp)
    num = (1 << Sp) * least(Q) - Cp
    f_adv = max(1, -((-num)//3**(u+1))) if num > 0 else 1
    return f_open, f_closed, f_adv

print("1. Fc' = max(m_closed(P), advance)  vs  Fc = max(m(P), advance)")
imp = eq = tot = 0; gains = []
for u in range(12):
    for v in range(12):
        for t in range(12):
            for qa in QS:
                for qb in QS:
                    fo, fc_, fa = branches(u,qa,v,qb,t)
                    old, new = max(fo,fa), max(fc_,fa)
                    m = least(word(u,qa,v,qb,t))
                    assert new <= m, (u,qa,v,qb,t)
                    tot += 1
                    if new > old:
                        imp += 1; gains.append(math.log2(new/old))
                    else: eq += 1
print(f"   configs = {tot}")
print(f"   Fc' > Fc : {imp}  ({100.0*imp/tot:.1f}%)   Fc' = Fc : {eq}"
      f"   Fc' < Fc : 0 (Fc' is never worse)")
if gains:
    print(f"   gain in bits when it improves: mean {st.fmean(gains):.3f}  "
          f"median {st.median(gains):.3f}  max {max(gains):.3f}")

print("\n2. which branch binds, and does the gain survive the max?")
print("   (the closed gain is only visible when truncation binds)")
tb = sum(1 for u in range(12) for v in range(12) for t in range(12)
         for qa in QS for qb in QS
         if branches(u,qa,v,qb,t)[1] >= branches(u,qa,v,qb,t)[2])
print(f"   truncation binds (closed) in {tb}/{tot} = {100.0*tb/tot:.1f}%")

print("\n3. balanced regime medians, |w| >= 62  (bits per |w|)")
rows = []
for T in range(60,181,3):
    b = T//3; u=v=b; t=T-2*b; n=T+2
    if n < 62: continue
    for qa in QS:
        for qb in QS:
            fo, fcl, fa = branches(u,qa,v,qb,t)
            m = least(word(u,qa,v,qb,t))
            rows.append((max(fo,fa).bit_length()/n, max(fcl,fa).bit_length()/n,
                         m.bit_length()/n))
med = lambda i: st.median(r[i] for r in rows)
print(f"   n={len(rows)}   Fc/|w| = {med(0):.4f}   Fc'/|w| = {med(1):.4f}"
      f"   empirical = {med(2):.4f}")
print(f"   improvement in the headline balanced median: "
      f"{med(1)-med(0):+.4f} per unit |w|")

print("\n4. per-regime at |w| ~ 122 (bits)")
T=120
SH={"balanced":(T//3,T//3,T-2*(T//3)),"leading-heavy":(T-2,1,1),
    "middle-heavy":(1,T-2,1),"trailing-heavy":(1,1,T-2)}
print(f"   {'regime':16s}{'Fc':>7s}{'Fc_prime':>10s}{'gain':>7s}{'emp':>7s}")
for lbl,(u,v,t) in SH.items():
    fo,fcl,fa = branches(u,2,v,2,t); m=least(word(u,2,v,2,t))
    print(f"   {lbl:16s}{max(fo,fa).bit_length():7d}"
          f"{max(fcl,fa).bit_length():10d}"
          f"{max(fcl,fa).bit_length()-max(fo,fa).bit_length():7d}"
          f"{m.bit_length():7d}")

print("\n5. sup log2(r/m) on a LARGER grid (does it grow?)")
for lim in (10, 14, 18):
    worst = 0.0
    for U in range(lim+1):
        for V in range(1, lim+1):
            for Q in range(2, 7):
                W = (1,)*U + (Q,) + (1,)*V
                worst = max(worst, math.log2(m_closed(W)/least(W)))
    print(f"   u,v <= {lim:2d}: sup log2(r/m) = {worst:.4f}")
