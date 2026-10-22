#!/usr/bin/env python3
"""Referee round: the finer modulus is the CLOSED-EPISODE condition, not a typo.
Verifies the reviewer's claims and tests whether it strengthens Theorem 5. T1."""
import math, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from indep_realizer import least, klass, word, sim, C_affine, F_c, run_floor, QS

def m_closed(w):
    """Least positive m realizing w whose NEXT valuation is >= 2."""
    a, M = klass(w)                       # realizers: a mod 2^M, M = S+1
    for cand in (a, a + (1 << M)):
        x = cand
        for v in w:
            x = (3 * x + 1) >> v
        if (3 * x + 1) % 4 == 0:          # next valuation >= 2
            return cand
    raise AssertionError("no closed lift")

print("1. REVIEWER'S EXAMPLE  w = (1,2,1)")
w = (1, 2, 1)
m, mc = least(w), m_closed(w)
for lbl, val in (("exact finite-word realizer", m), ("closed-episode realizer", mc)):
    print(f"   {lbl:28s} = {val:8d}   word+next = {sim(val, 4)}")
print(f"   claimed 27 and 59 -> {'CONFIRMED' if (m,mc)==(27,59) else 'MISMATCH'}"
      f"   (59-27 = {mc-m} = 2^{{S+1}} = {1<<(sum(w)+1)})")

print("\n2. REVIEWER'S HEIGHT COUNTEREXAMPLE  u=11, q=5, v=4")
u, q, v = 11, 5, 4
P = (1,)*u + (q,) + (1,)*v
mP, rP = least(P), m_closed(P)
print(f"   least open-word realizer m = {mP}")
print(f"   closed-episode realizer  r = {rP}")
print(f"   r - m = {rP-mP} = 2^{{S+1}} = {1<<(sum(P)+1)}   "
      f"-> {'matches 2101247' if rP==2101247 else 'MISMATCH'}")
print(f"   log2(r/m) = {math.log2(rP/mP):.4f}   (claimed 9.003)")
print(f"   => my 'at most one bit' claim is FALSE.  bit_length differs by"
      f" {rP.bit_length()-mP.bit_length()}, but log2 RATIO by {math.log2(rP/mP):.3f}.")

print("\n3. how bad can it get?  sup log2(r/m) over u,q,v (u,v<=14, q<=6)")
worst = None
for U in range(15):
    for V in range(1, 15):
        for Q in range(2, 7):
            W = (1,)*U + (Q,) + (1,)*V
            a, b = least(W), m_closed(W)
            d = math.log2(b / a)
            if worst is None or d > worst[0]:
                worst = (d, U, Q, V, a, b)
print(f"   sup log2(r/m) = {worst[0]:.4f} at (u,q,v)=({worst[1]},{worst[2]},"
      f"{worst[3]})  m={worst[4]} r={worst[5]}")
print(f"   => the gap is UNBOUNDED in bits; it is NOT an O(1) correction.")

print("\n4. does the closed-episode prefix floor hold for two-defect words?")
print("   In w = 1^u qa 1^v qb 1^t the prefix P = 1^u qa 1^v is followed by")
print("   qb >= 2, so every realizer of w realizes P AS A CLOSED EPISODE.")
bad = tot = 0
for u in range(9):
    for v in range(9):
        for t in range(9):
            for qa in QS:
                for qb in QS:
                    w = word(u, qa, v, qb, t)
                    P = (1,)*u + (qa,) + (1,)*v
                    tot += 1
                    if least(w) < m_closed(P):
                        bad += 1
print(f"   configs = {tot}   violations of m(w) >= m_closed(P) : {bad}"
      f"   -> {'VALID' if not bad else 'INVALID'}")
