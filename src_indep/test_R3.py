#!/usr/bin/env python3
"""PHASE A1 -- test R3 (suffix-combination floor F_c) independently."""
import os, statistics, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from indep_realizer import least, word, F_c, run_floor, QS

# ---- 1. is F_c a valid floor? -------------------------------------------
print("1. VALIDITY  F_c(w) <= m(w) ?")
bad = tot = 0
worst = None
for u in range(13):
    for v in range(13):
        for t in range(13):
            for qa in QS:
                for qb in QS:
                    w = word(u, qa, v, qb, t)
                    m = least(w)
                    fc, f1, f2 = F_c(u, qa, v, qb, t)
                    tot += 1
                    if fc > m:
                        bad += 1
                        if worst is None:
                            worst = (u, qa, v, qb, t, m, fc, f1, f2)
print(f"   configs={tot}  violations={bad}  -> "
      f"{'VALID FLOOR' if not bad else 'NOT A FLOOR'}")
if worst:
    print(f"   first violation: {worst[:5]}  m={worst[5]} F_c={worst[6]} "
          f"(trunc={worst[7]} adv={worst[8]})")

# ---- 2. does F_c beat the run-based floor my package used? --------------
print("\n2. F_c vs the RUN-based floor used in my 2026-10-22 package")
beat = tie = lose = 0
for u in range(13):
    for v in range(13):
        for t in range(13):
            for qa in QS:
                for qb in QS:
                    fc, _, _ = F_c(u, qa, v, qb, t)
                    rf = run_floor(u, qa, v, qb, t)
                    if fc > rf: beat += 1
                    elif fc == rf: tie += 1
                    else: lose += 1
print(f"   F_c > run-floor : {beat}")
print(f"   F_c = run-floor : {tie}")
print(f"   F_c < run-floor : {lose}")
print(f"   -> F_c dominates: {'YES' if lose==0 else 'NO'}"
      f"   ({100.0*beat/(beat+tie+lose):.1f}% strictly better)")

# ---- 3. balanced-regime medians, |w| >= 62 ------------------------------
print("\n3. BALANCED regime (u=v=t), |w| >= 62: medians of bits/|w|")
rows = []
for T in range(60, 181, 3):
    b = T // 3
    u = v = b; t = T - 2 * b
    n = T + 2
    if n < 62:
        continue
    for qa in QS:
        for qb in QS:
            m = least(word(u, qa, v, qb, t))
            fc, f1, f2 = F_c(u, qa, v, qb, t)
            rf = run_floor(u, qa, v, qb, t)
            rows.append((fc.bit_length() / n, rf.bit_length() / n,
                         m.bit_length() / n, f1.bit_length() / n,
                         f2.bit_length() / n))
med = lambda i: statistics.median(r[i] for r in rows)
print(f"   n_configs = {len(rows)}")
print(f"   F_c / |w|        = {med(0):.3f}      [brief claims 0.671]")
print(f"   run-floor / |w|  = {med(1):.3f}      [brief claims 0.337]")
print(f"   empirical / |w|  = {med(2):.3f}      [brief claims 1.015]")
print(f"     (truncation branch alone = {med(3):.3f}, "
      f"advance branch alone = {med(4):.3f})")

# ---- 4. which branch binds? --------------------------------------------
print("\n4. which branch of F_c binds, by regime (|w| ~ 122)")
T = 120
SH = {"balanced":(T//3,T//3,T-2*(T//3)), "leading-heavy":(T-2,1,1),
      "middle-heavy":(1,T-2,1), "trailing-heavy":(1,1,T-2)}
print(f"   {'regime':18s}{'trunc':>8s}{'adv':>8s}{'F_c':>8s}"
      f"{'runF':>8s}{'emp':>8s}  binds")
for lbl,(u,v,t) in SH.items():
    fc,f1,f2 = F_c(u,2,v,2,t)
    rf = run_floor(u,2,v,2,t); m = least(word(u,2,v,2,t))
    n = u+v+t+2
    print(f"   {lbl:18s}{f1.bit_length():8d}{f2.bit_length():8d}"
          f"{fc.bit_length():8d}{rf.bit_length():8d}{m.bit_length():8d}"
          f"  {'truncation' if f1>=f2 else 'advance'}")
