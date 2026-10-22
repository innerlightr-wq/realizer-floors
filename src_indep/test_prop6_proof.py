#!/usr/bin/env python3
"""Verify each step of the Proposition 6 proof exposition.  T1 exact."""
import math, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from indep_realizer import least, word, C_affine
from test_closed_episode import m_closed

print("S1. -15 = 1 (mod 8):", (-15) % 8 == 1)
print("S2. 3^e = 1 (mod 8) iff e even:",
      all((pow(3,e,8)==1) == (e%2==0) for e in range(1,40)))
print("S3. <9> = {x = 1 mod 8} in (Z/2^k)*, k=4..12:")
ok3 = True
for k in range(4, 13):
    M = 1 << k
    gen9 = {pow(9, j, M) for j in range(M)}
    one8 = {x for x in range(1, M, 2) if x % 8 == 1}
    ok3 &= (gen9 == one8)
print("    ", ok3, " (so every residue =1 mod 8 is an even power of 3)")

print("S4. solutions e of 3^e = -15 (mod 2^{v+6}) form a progression, all even:")
for v in (1,2,3,4,5):
    M = 1 << (v+6)
    sols = [e for e in range(1, 4*M) if pow(3,e,M) == (-15) % M]
    step = sols[1]-sols[0] if len(sols) > 1 else None
    print(f"    v={v}: first={sols[0]} step={step} ord(3)=2^{v+4}={1<<(v+4)}"
          f" all even={all(e%2==0 for e in sols)}")

print("S5. A=1 solves the OPEN congruence mod 2^{v+5} but FAILS the CLOSED one")
print("    mod 2^{v+6}  (q=5, n=q+v):")
for v in (1,2,3,4):
    M6 = 1 << (v+6); M5 = 1 << (v+5)
    e = next(x for x in range(1,4*M6) if pow(3,x,M6)==(-15)%M6)
    open_ok   = pow(3,e,M5) == (1 - (1<<4)) % M5
    closed_ok = pow(3,e,M6) == (1 - (1<<4) + (1<<(5+v))) % M6
    print(f"    v={v} e={e}: open(A=1)={open_ok}  closed(A=1)={closed_ok}"
          f"  -> m_closed is the OTHER lift: {not closed_ok}")

print("S6. advance branch is dominated: ratio advance/truncation -> 0")
print("    advance ~ 2^{e+4} m(Q)/3^e, truncation = 2^e-1, so ratio ~ 16 m(Q)/3^e")
for v in (1,2):
    M = 1 << (v+6)
    e0 = next(x for x in range(1,4*M) if pow(3,x,M)==(-15)%M)
    for e in (e0, e0+(1<<(v+4)), e0+2*(1<<(v+4))):
        u = e-1
        P  = (1,)*u + (5,) + (1,)*v
        Pp = (1,)*u + (5,)
        Q  = (1,)*v + (2,) + (1,)*1
        Sp, Cp = u+5, C_affine(Pp)
        num = (1 << Sp)*least(Q) - Cp
        adv = max(1, -((-num)//3**(u+1))) if num>0 else 1
        tr  = least(P)
        print(f"    v={v} e={e:3d}: trunc={tr.bit_length():4d}b "
              f"advance={adv.bit_length():4d}b  dominated={adv<tr}"
              f"  m(P)=2^e-1: {tr==(1<<e)-1}"
              f"  m_closed=2^e-1+2^{{e+v+5}}: {m_closed(P)==(1<<e)-1+(1<<(e+v+5))}")
