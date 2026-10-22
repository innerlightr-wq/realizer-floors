#!/usr/bin/env python3
"""Single authoritative script regenerating EVERY number in the consolidated
note.  T1 exact integers except where marked [T2].  Standard library only.

Run: python3 src/verify_note_numbers.py
"""
import math, os, random, sys
from math import gcd
from collections import Counter
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))), "src_indep"))
from indep_realizer import (least, klass, word, F_c, run_floor, C_affine,
                            C_R1, sim, QS)

def v3(x):
    k = 0
    while x and x % 3 == 0:
        x //= 3; k += 1
    return k

def m2(u, q, v):
    return least((1,) * u + (q,) + (1,) * v)

OUT = []
def say(s=""):
    print(s); OUT.append(s)

say("=" * 72)
say("N1  modulus convention: realizer class is one class mod 2^{S+1}")
say("=" * 72)
bad = scans = 0
for u in range(3):
    for v in range(3):
        for t in range(3):
            for qa in QS:
                for qb in QS:
                    w = word(u, qa, v, qb, t); S = sum(w)
                    if S > 15: continue
                    a, M = klass(w)
                    found = [m for m in range(1, (1 << M) + 1, 2)
                             if sim(m, len(w)) == w]
                    scans += 1
                    bad += (found != [a] or M != S + 1)
say(f"    full odd scans = {scans}   failures = {bad}")

say()
say("=" * 72)
say("N2  R1 three-run closed form")
say("=" * 72)
n1 = b1 = 0
for u in range(13):
    for v in range(13):
        for t in range(13):
            for qa in QS:
                for qb in QS:
                    n1 += 1
                    b1 += (C_R1(u,qa,v,qb,t) != C_affine(word(u,qa,v,qb,t)))
say(f"    configs = {n1}   mismatches = {b1}")

say()
say("=" * 72)
say("N3  R2 corrected: denominator D = 3^{n - v3(2^S - C_w)}, by stratum")
say("=" * 72)
strata = {}
for u in range(14):
    for v in range(14):
        for t in range(14):
            for qa in QS:
                for qb in QS:
                    w = word(u,qa,v,qb,t); S = sum(w)
                    N = (1 << S) - C_affine(w)
                    strata.setdefault((u>=1,v>=1,t>=1), Counter())[v3(N)] += 1
say(f"    {'u>=1':>6s}{'v>=1':>6s}{'t>=1':>6s}{'n':>8s}   v3 histogram")
for key in sorted(strata, reverse=True):
    c = strata[key]
    say(f"    {str(key[0]):>6s}{str(key[1]):>6s}{str(key[2]):>6s}"
        f"{sum(c.values()):>8d}   " + " ".join(f"{k}:{c[k]}" for k in sorted(c)))
interior = strata[(True,True,True)]
tpos = sum(sum(v for k,v in strata[k2].items() if k>0)
           for k2 in strata if k2[2])
t0 = sum(sum(v for k,v in strata[k2].items() if k>0)
         for k2 in strata if not k2[2])
say(f"    interior u,v,t>=1 : v3 values = {sorted(interior)}  "
    f"(3-coprime: {set(interior)=={0}})")
say(f"    reductions v3>0 : t=0 -> {t0}   t>=1 -> {tpos}")
w0 = (2,2)
say(f"    boundary witness w=(2,2) [the x=1 cycle word]: m={least(w0)}, "
    f"2^S-C={(1<<sum(w0))-C_affine(w0)}, v3={v3((1<<sum(w0))-C_affine(w0))}, "
    f"D={3**len(w0)//gcd((1<<sum(w0))-C_affine(w0),3**len(w0))}")

say()
say("=" * 72)
say("N4  R3 suffix-combination floor F_c")
say("=" * 72)
viol = beat = tie = lose = tot = 0
for u in range(13):
    for v in range(13):
        for t in range(13):
            for qa in QS:
                for qb in QS:
                    m = least(word(u,qa,v,qb,t))
                    fc,_,_ = F_c(u,qa,v,qb,t); rf = run_floor(u,qa,v,qb,t)
                    tot += 1; viol += (fc > m)
                    if fc > rf: beat += 1
                    elif fc == rf: tie += 1
                    else: lose += 1
say(f"    configs = {tot}   F_c > m(w) violations = {viol}")
say(f"    vs all-ones-run floor:  better {beat}   tie {tie}   worse {lose}")
say(f"    strictly better on {100.0*beat/tot:.1f}%")
import statistics
rows = []
for T in range(60,181,3):
    b = T//3; u=v=b; t=T-2*b; n=T+2
    if n < 62: continue
    for qa in QS:
        for qb in QS:
            m = least(word(u,qa,v,qb,t)); fc,f1,f2 = F_c(u,qa,v,qb,t)
            rf = run_floor(u,qa,v,qb,t)
            rows.append((fc.bit_length()/n, rf.bit_length()/n,
                         m.bit_length()/n))
med = lambda i: statistics.median(r[i] for r in rows)
say(f"    balanced |w|>=62, n={len(rows)}: F_c/|w|={med(0):.3f}  "
    f"run-floor/|w|={med(1):.3f}  empirical/|w|={med(2):.3f}")

say()
say("=" * 72)
say("N5  R4 uniform reduction thresholds")
say("=" * 72)
al = math.log2(3)
say(f"    {'k':>6s}{'truncation':>14s}{'advance':>12s}")
for k in (1,2,5,10,50,500,5000):
    say(f"    {k:6d}{(k+1)/(2*k+2):14.6f}{(al*(k+1)-1)/(2*k+2):12.6f}")
say(f"    truncation threshold = 1/2 exactly for every k")
say(f"    advance threshold -> alpha/2 = {al/2:.6f}")
say("    empirical beta for 1^k 2 1^k:")
for k in (5,20,60,150):
    P = (1,)*k+(2,)+(1,)*k
    say(f"      k={k:4d}  S_P={sum(P):4d}  beta_emp={least(P).bit_length()/sum(P):.4f}")

say()
say("=" * 72)
say("N6  R5 merge bridge, canonical image m2(u, qa+qb, v+t)")
say("=" * 72)
u,qa,v,qb,t = 22,4,2,2,6
w = word(u,qa,v,qb,t); n=len(w); mw=least(w)
say(f"    witness (22,4,2,2,6): |w|={n}")
for lbl,img in (("canonical m2(u,qa+qb,v+t)",(u,qa+qb,v+t)),
                ("keep      m2(u,qa,v+t)   ",(u,qa,v+t)),
                ("additive  m2(u,qa+qb-1,..)",(u,qa+qb-1,v+t))):
    mv = m2(*img)
    say(f"      {lbl}: bits {mv.bit_length():3d}  "
        f"delta_bits={(mv.bit_length()-mw.bit_length())/n:+.4f}  "
        f"delta_log2={(math.log2(mv)-math.log2(mw))/n:+.4f}   [T2 log2]")
for gl, rng, interior_only in (("small full grid u,v,t<=25", range(26), False),
                               ("small grid, INTERIOR u,v,t>=1", range(1,26), True),
                               ("|w|>=32 sampled stride 7", range(0,120,7), False)):
    best = None
    for U in rng:
        for V in rng:
            for T2_ in rng:
                nn = U+V+T2_+2
                if gl.startswith("|w|>=32") and nn < 32: continue
                if nn < 6: continue
                for A in QS:
                    for B in QS:
                        mm = least(word(U,A,V,B,T2_))
                        mv = m2(U,A+B,V+T2_)
                        d = (math.log2(mv)-math.log2(mm))/nn
                        if best is None or d > best[0]:
                            best = (d,U,A,V,B,T2_,nn)
    say(f"    sup delta_log2, {gl}: {best[0]:+.4f} at "
        f"({best[1]},{best[2]},{best[3]},{best[4]},{best[5]}), |w|={best[6]}")

say()
say("=" * 72)
say("N7  R6 residual-vs-anchor correlation  [T2]  seed 20261023, n=100")
say("=" * 72)
def pearson(xs,ys):
    n=len(xs); mx=sum(xs)/n; my=sum(ys)/n
    sxy=sum((a-mx)*(b-my) for a,b in zip(xs,ys))
    sxx=sum((a-mx)**2 for a in xs); syy=sum((b-my)**2 for b in ys)
    return sxy/math.sqrt(sxx*syy)
random.seed(20261023)
samp=[(random.randint(1,40),random.choice(QS),random.randint(1,40),
       random.choice(QS),random.randint(1,40)) for _ in range(100)]
res,dQ,dP,vv=[],[],[],[]
for (u,qa,v,qb,t) in samp:
    m=least(word(u,qa,v,qb,t)); fc,_,_=F_c(u,qa,v,qb,t)
    Q=(1,)*v+(qb,)+(1,)*t; P=(1,)*u+(qa,)+(1,)*v
    res.append(math.log2(m)-math.log2(fc))
    dQ.append(sum(Q)+1-math.log2(least(Q)))
    dP.append(sum(P)+1-math.log2(least(P))); vv.append(v)
say(f"    residual above F_c: mean={sum(res)/100:.3f} "
    f"min={min(res):.3f} max={max(res):.3f}")
for lbl,ys in (("deficit of m(Q)",dQ),("deficit of m(P)",dP),
               ("middle-run length v",vv)):
    r=pearson(res,ys)
    say(f"      r(residual, {lbl:22s}) = {r:+.4f}  "
        f"{'consistent with zero' if abs(r)<0.197 else 'significant p<0.05'}")
say("    (|r|>0.197 is the p<0.05 two-sided threshold at n=100)")

say()
say("=" * 72)
say("N8  two-run import: notation reconciliation check")
say("=" * 72)
say("    TWO_RUN_FLOOR.md states: m realizes 1^u q 1^v IFF m+1=2^{u+1}A,")
say("    A odd, 3^{u+1}A = c (mod 2^{n+1}), c = 1-2^{q-1}+2^n, n=q+v.")
eq = lt = tot8 = 0
for u in range(12):
    for v in range(1,12):
        for q in (2,3,4,5):
            w=(1,)*u+(q,)+(1,)*v
            m0,_=klass(w); n_=q+v; c=1-(1<<(q-1))+(1<<n_); M=1<<(n_+1)
            A=(pow(pow(3,u+1,M),-1,M)*c)%M
            r=(1<<(u+1))*A-1
            tot8+=1; eq+=(r==m0); lt+=(r!=m0)
say(f"    configs = {tot8}   r == least realizer: {eq}   r != least: {lt}")
say(f"    => the mod-2^{{n+1}} criterion is ONE BIT FINER than the realizer")
say(f"       class (which is mod 2^{{S+1}}, i.e. A mod 2^n).  It selects one")
say(f"       of the two realizers in each class mod 2^{{S+2}}, so 'least' holds")
say(f"       in {eq}/{tot8} cases only.  NOTATION MUST BE RECONCILED before §2.")

print()
print(f"[verify_note_numbers] {len(OUT)} lines emitted")
