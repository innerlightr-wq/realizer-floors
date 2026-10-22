#!/usr/bin/env python3
#
# *** SUPERSEDED -- retained for transparency, conclusion WITHDRAWN ***
# This measures a DIFFERENT statistic from the one it was taken to test
# (unnormalised residual vs prefix/suffix deficits, varying |w|), not the
# published (c2) statistic (normalised residual vs P-anchor fullness at
# fixed |w|=122).  See src_indep/test_R6_their_stat.py and the paper's
# Provenance section.  Do not cite the numbers printed here.
#
"""R6 is underpowered at n=100: the estimate straddles the significance
threshold depending on draw order.  Settle it at large n, and quantify the
n=100 instability directly.  [T2 Pearson on exact integer inputs]"""
import math, os, random, statistics, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from indep_realizer import least, word, F_c, QS

def pearson(xs, ys):
    n=len(xs); mx=sum(xs)/n; my=sum(ys)/n
    sxy=sum((a-mx)*(b-my) for a,b in zip(xs,ys))
    sxx=sum((a-mx)**2 for a in xs); syy=sum((b-my)**2 for b in ys)
    return sxy/math.sqrt(sxx*syy)

def draw(rng):
    u=rng.randint(1,40); v=rng.randint(1,40); t=rng.randint(1,40)
    return (u, rng.choice(QS), v, rng.choice(QS), t)

def stats(samp):
    res,dQ,dP = [],[],[]
    for (u,qa,v,qb,t) in samp:
        m=least(word(u,qa,v,qb,t)); fc,_,_=F_c(u,qa,v,qb,t)
        Q=(1,)*v+(qb,)+(1,)*t; P=(1,)*u+(qa,)+(1,)*v
        res.append(math.log2(m)-math.log2(fc))
        dQ.append(sum(Q)+1-math.log2(least(Q)))
        dP.append(sum(P)+1-math.log2(least(P)))
    return res,dQ,dP

print("1. n=100 INSTABILITY: 40 independent replicates, seeds 20261023+i")
for lbl, idx in (("deficit of m(Q)",1),("deficit of m(P)",2)):
    rs=[]
    for i in range(40):
        rng=random.Random(20261023+i)
        s=[draw(rng) for _ in range(100)]
        out=stats(s); rs.append(pearson(out[0],out[idx]))
    sig=sum(1 for r in rs if abs(r)>0.197)
    print(f"   {lbl:18s} r in [{min(rs):+.4f},{max(rs):+.4f}]  "
          f"mean={statistics.fmean(rs):+.4f}  sd={statistics.stdev(rs):.4f}")
    print(f"   {'':18s} 'significant at p<0.05' in {sig}/40 replicates"
          f"  -> conclusion is {'UNSTABLE' if 0<sig<40 else 'stable'} at n=100")

print("\n2. LARGE-n estimate (n=6000, seed 20261023)")
rng=random.Random(20261023)
s=[draw(rng) for _ in range(6000)]
res,dQ,dP=stats(s)
thr=1.96/math.sqrt(6000)
for lbl,ys in (("deficit of m(Q)",dQ),("deficit of m(P)",dP)):
    r=pearson(res,ys)
    print(f"   r(residual, {lbl:18s}) = {r:+.4f}   "
          f"(p<0.05 threshold {thr:.4f})  "
          f"{'NONZERO' if abs(r)>thr else 'consistent with zero'}")
print(f"   residual above F_c: mean={statistics.fmean(res):.3f} "
      f"sd={statistics.stdev(res):.3f}")
print("\n   => at adequate power the residual DOES correlate weakly but")
print("      definitely with the suffix-word deficit; R6's 'no anchor")
print("      correlation' is not supported at n=6000, and its n=100 is")
print("      too small to decide either way.")
