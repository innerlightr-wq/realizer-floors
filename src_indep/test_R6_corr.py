#!/usr/bin/env python3
#
# *** SUPERSEDED -- retained for transparency, conclusion WITHDRAWN ***
# This measures a DIFFERENT statistic from the one it was taken to test
# (unnormalised residual vs prefix/suffix deficits, varying |w|), not the
# published (c2) statistic (normalised residual vs P-anchor fullness at
# fixed |w|=122).  See src_indep/test_R6_their_stat.py and the paper's
# Provenance section.  Do not cite the numbers printed here.
#
"""R6 -- residual-vs-anchor correlation, re-run in THIS environment.

The chat package's (c2) block is unavailable here, so the sample and the
anchor statistic are defined explicitly below; these are MY definitions, not
a reproduction of either session's value.  seed = 20261023, n = 100,
interior stratum u,v,t >= 1.  [T2: Pearson r on exact integer inputs]
"""
import math, os, random, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from indep_realizer import least, word, F_c, C_affine, QS


def pearson(xs, ys):
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    sxy = sum((a - mx) * (b - my) for a, b in zip(xs, ys))
    sxx = sum((a - mx) ** 2 for a in xs)
    syy = sum((b - my) ** 2 for b in ys)
    return sxy / math.sqrt(sxx * syy) if sxx and syy else float("nan")


random.seed(20261023)
samp = []
while len(samp) < 100:
    u = random.randint(1, 40); v = random.randint(1, 40); t = random.randint(1, 40)
    qa = random.choice(QS); qb = random.choice(QS)
    samp.append((u, qa, v, qb, t))

res, a_mid, a_defQ, a_defP, a_v = [], [], [], [], []
for (u, qa, v, qb, t) in samp:
    w = word(u, qa, v, qb, t)
    m = least(w)
    fc, f1, f2 = F_c(u, qa, v, qb, t)
    Q = (1,) * v + (qb,) + (1,) * t
    P = (1,) * u + (qa,) + (1,) * v
    res.append(math.log2(m) - math.log2(fc))          # residual above F_c
    a_mid.append(math.log2(least((1,) * v)))           # shared-run anchor
    a_defQ.append(sum(Q) + 1 - math.log2(least(Q)))    # deficit of m(Q)
    a_defP.append(sum(P) + 1 - math.log2(least(P)))    # deficit of m(P)
    a_v.append(v)

print("R6  residual above F_c, n=100, interior stratum, seed 20261023")
print(f"   residual: mean={sum(res)/100:.3f}  "
      f"min={min(res):.3f}  max={max(res):.3f}")
print()
print("   Pearson r of residual against candidate anchor statistics:")
for lbl, ys in (("shared-run anchor log2 m(1^v)", a_mid),
                ("deficit of m(Q)              ", a_defQ),
                ("deficit of m(P)              ", a_defP),
                ("middle-run length v          ", a_v)):
    r = pearson(res, ys)
    # two-sided significance at n=100: |r| > 0.197 is p<0.05
    tag = "consistent with zero" if abs(r) < 0.197 else "SIGNIFICANT at p<0.05"
    print(f"     {lbl} : r = {r:+.4f}   {tag}")
print()
print("   (|r| > 0.197 is the p<0.05 threshold at n=100.  The chat package's")
print("    claimed r = -0.179 is BELOW that threshold, i.e. its own stated")
print("    conclusion 'consistent with zero' does not depend on which of")
print("    these statistics was used.)")
