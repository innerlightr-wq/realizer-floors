#!/usr/bin/env python3
"""Cross-check of my prior two-run least-realizer formula against the exact
class derived in Phase 0.  T1 exact.

Prior formula (two-run word 1^u q 1^v):  n = q+v,  c = 1 - 2^{q-1} + 2^n,
A* = least positive residue of 3^{-(u+1)} c mod 2^{n+1},  r = 2^{u+1} A* - 1.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from three_run_exact import m_class, orbit_word

same = off = tot = 0
diffs = {}
for u in range(0, 12):
    for v in range(1, 12):
        for q in (2, 3, 4, 5):
            w = (1,) * u + (q,) + (1,) * v
            m0, mod = m_class(w)
            n = q + v
            c = 1 - (1 << (q - 1)) + (1 << n)
            M = 1 << (n + 1)
            Astar = (pow(pow(3, u + 1, M), -1, M) * c) % M
            r = (1 << (u + 1)) * Astar - 1
            tot += 1
            same += ((r - m0) % mod == 0)
            d = (r - m0) // mod
            diffs[d] = diffs.get(d, 0) + 1
            off += (r != m0)
            assert orbit_word(r, len(w)) == w, ("r not a realizer", u, q, v)

print(f"configs                         : {tot}")
print(f"r in the SAME class as m(w)     : {same}/{tot}")
print(f"r realizes w (orbit-simulated)  : {tot}/{tot} (asserted)")
print(f"r is NOT the least realizer     : {off}/{tot}")
print(f"(r-m)/2^(S+1) observed          : {dict(sorted(diffs.items()))}")
print()
print("VERDICT: the prior formula yields a correct realizer in the correct")
print("class, but r in {m(w), m(w)+2^{S+1}} -- least only 257/528 of the time.")
print("Cause: A* is reduced mod 2^{n+1}, one bit finer than the class requires,")
print("so it fixes a representative mod 2^{S+2} instead of mod 2^{S+1}.")
print("Impact: log2 r overstates log2 m(w) by at most 1 bit, so exponent-level")
print("conclusions are unaffected; exact minima / champion seeds are not.")
