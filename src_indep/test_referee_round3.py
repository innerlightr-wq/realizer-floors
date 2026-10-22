#!/usr/bin/env python3
"""Referee round 3: verify (A) unbounded gain, (B) denominator truncation,
(C) Hypothesis beta domain.  T1 exact."""
import math, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from indep_realizer import least, klass, word, sim, C_affine
from test_closed_episode import m_closed

def Fc_pair(u, qa, v, qb, t):
    P  = (1,)*u + (qa,) + (1,)*v
    Pp = (1,)*u + (qa,)
    Q  = (1,)*v + (qb,) + (1,)*t
    fo, fcl = least(P), m_closed(P)
    Sp, Cp = u+qa, C_affine(Pp)
    num = (1 << Sp)*least(Q) - Cp
    fa = max(1, -((-num)//3**(u+1))) if num > 0 else 1
    return max(fo,fa), max(fcl,fa), fo, fcl, fa

print("="*70); print("A. the reviewer's two-defect example  w = 1^11 5 1^4 2 1")
u,qa,v,qb,t = 11,5,4,2,1
Fc, Fcp, fo, fcl, fa = Fc_pair(u,qa,v,qb,t)
m = least(word(u,qa,v,qb,t))
print(f"   F_c(w)   = {Fc:>12d}   (claimed 4095)")
print(f"   F'_c(w)  = {Fcp:>12d}   (claimed 2101247)")
print(f"   m(w)     = {m:>12d}   (claimed 14684159)")
print(f"   gain     = {math.log2(Fcp/Fc):.4f} bits   (claimed 9.003)")
print(f"   branches: trunc_open={fo} trunc_closed={fcl} advance={fa}")
print(f"   both floors valid: {Fc<=m and Fcp<=m}")

print("\n"+"="*70); print("A2. the infinite family  w = 1^{e-1} 5 1^v 2 1,  3^e = -15 mod 2^{v+6}")
print(f"   {'v':>3s}{'e':>6s}{'F_c':>14s}{'F_c_prime':>18s}{'gain bits':>11s}{'> v+5?':>8s}")
ok = True
for v in range(1, 8):
    M = 1 << (v+6)
    tgt = (-15) % M
    e = None
    for cand in range(2, 4*M, 2):               # even exponents
        if pow(3, cand, M) == tgt:
            e = cand; break
    if e is None:
        print(f"   {v:3d}     no even e found mod 2^{v+6}"); ok=False; continue
    # need e large enough that advance < truncation; step up by ord(3) if needed
    for bump in range(0, 40):
        ee = e + bump*(M//2 if M//2 else 1)
        if ee < 3: continue
        U = ee-1
        Fc, Fcp, fo, fcl, fa = Fc_pair(U,5,v,2,1)
        if fo > fa:                              # truncation binds, as required
            break
    gain = math.log2(Fcp/Fc)
    trunc_pred = (1 << ee) - 1
    closed_pred = (1 << ee) - 1 + (1 << (ee+v+5))
    match = (fo == trunc_pred) and (fcl == closed_pred)
    print(f"   {v:3d}{ee:6d}{Fc:>14d}{Fcp:>18d}{gain:>11.4f}{str(gain> v+5):>8s}"
          f"   formulas {'match' if match else 'DIFFER'}")
    ok = ok and (gain > v+5)
print(f"   => gain > v+5 for every v tested: {ok}")
print(f"   => the improvement is NOT uniformly O(1) bits.")

print("\n"+"="*70); print("B. denominator at the boundary:  w = (6,2)")
w = (6,2); n=len(w); S=sum(w); C=C_affine(w); N=(1<<S)-C
def v3(x):
    k=0
    while x and x%3==0: x//=3; k+=1
    return k
from math import gcd
print(f"   n={n} S={S} C_w={C}  2^S-C_w={N} = 3^{v3(N)} * {N//3**v3(N)}")
print(f"   v3 = {v3(N)} > n = {n}")
print(f"   naive 3^(n-v3) = 3^{n-v3(N)}  -> NOT an integer denominator")
print(f"   true reduced denominator = {3**n//gcd(N,3**n)}   (anchor {N}/{3**n} = {N//gcd(N,3**n)}/{3**n//gcd(N,3**n)})")
print(f"   => D = 3^max(0, n-v3) = 3^{max(0,n-v3(N))} = {3**max(0,n-v3(N))}  CORRECT")
# how often does v3 exceed n?
bad=tot=0
for U in range(7):
    for V in range(7):
        for T in range(7):
            for A in range(2,8):
                for B in range(2,8):
                    ww=word(U,A,V,B,T); NN=(1<<sum(ww))-C_affine(ww)
                    tot+=1
                    if v3(NN) > len(ww): bad+=1
print(f"   v3 > n on {bad}/{tot} configs (qa,qb up to 7) -> truncation is needed")

print("\n"+"="*70); print("C. Hypothesis beta: admissible domain")
w14=(1,4); m14=least(w14)
print(f"   w=(1,4) [U=1,q=4,V=0]: m={m14}, S={sum(w14)}, "
      f"log2(m)/S = {math.log2(m14)/sum(w14):.4f}  < 1/2 -> beta>1/2 FALSE if V=0 allowed")
print(f"   orbit check: {sim(m14,2)}")
worst=None
for U in range(0,26):
    for V in range(0,26):
        for Q in range(2,8):
            P=(1,)*U+(Q,)+(1,)*V
            r=math.log2(least(P))/sum(P)
            if worst is None or r<worst[0]: worst=(r,U,Q,V)
print(f"   min log2 m / S over U,V in 0..25, q in 2..7 : {worst[0]:.4f} "
      f"at (U,q,V)=({worst[1]},{worst[2]},{worst[3]})")
wi=None
for U in range(1,26):
    for V in range(1,26):
        for Q in range(2,8):
            P=(1,)*U+(Q,)+(1,)*V
            r=math.log2(least(P))/sum(P)
            if wi is None or r<wi[0]: wi=(r,U,Q,V)
print(f"   INTERIOR U,V>=1 : min = {wi[0]:.4f} at (U,q,V)=({wi[1]},{wi[2]},{wi[3]})"
      f"   {'> 1/2 OK' if wi[0]>0.5 else '<= 1/2 PROBLEM'}")
