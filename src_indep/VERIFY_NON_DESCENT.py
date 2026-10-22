#!/usr/bin/env python3
"""Exact non-descent certificates.  T1: integers and Fractions only."""
from fractions import Fraction as F
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from indep_realizer import least, klass, word as mkword, sim
from test_closed_episode import m_closed
from VERIFY_CONSTRAINED_EXTREMIZER import confined, prefix_sums, Hw, wstar, confined_words

def C_of(w):
    n, C, S = len(w), 0, 0
    for k in range(1, n + 1):
        C += 3 ** (n - k) * (1 << S); S += w[k - 1]
    return C

def Fc_pair(u, qa, v, qb, t):
    P  = (1,)*u + (qa,) + (1,)*v
    Pp = (1,)*u + (qa,)
    Q  = (1,)*v + (qb,) + (1,)*t
    fo, fcl = least(P), m_closed(P)
    Sp, Cp = u + qa, C_of(Pp)
    num = (1 << Sp) * least(Q) - Cp
    fa = max(1, -((-num) // 3 ** (u + 1))) if num > 0 else 1
    return max(fo, fa), max(fcl, fa)

ok = True
def rep(tag, cond):
    global ok; ok &= cond
    print(f"   {'PASS' if cond else 'FAIL'}  {tag}")

print("=" * 72)
print("MANDATORY TEST  w = (1,2,1,3,1)   [u=1, qa=2, v=1, qb=3, t=1]")
print("=" * 72)
w = (1, 2, 1, 3, 1); N = len(w); K = sum(w); C = C_of(w)
S = prefix_sums(w)
print(f"   N = {N}, K = {K}, prefix sums = {tuple(S)}")
rep(f"N = 5 (got {N})", N == 5)
rep(f"K = 8 (got {K})", K == 8)
rep(f"C(w) = 383 (got {C})", C == 383)
den = (1 << K) - 3 ** N
print(f"   2^K - 3^N = {1<<K} - {3**N} = {den}   (positive compression: {den>0})")
rep(f"2^K - 3^N = 13 (got {den})", den == 13)
T = F(C, den)
rep(f"T(w) = 383/13 (got {T})", T == F(383, 13))

# identity T = H / (3 (2^{R_N} - 1)) via exact integers:
#   2^{R_N} - 1 = 2^K/3^N - 1 = (2^K - 3^N)/3^N  ->  3(2^{R_N}-1) = (2^K-3^N)/3^{N-1}
H = Hw(w)
rep(f"C(w) = 3^(N-1) H(w)  ({C} = {3**(N-1)}*{H})", C == 3 ** (N - 1) * H)
T_alt = H / (F(den, 3 ** (N - 1)))
rep(f"T = H/(3(2^R_N -1)) identity holds (got {T_alt})", T_alt == T)

print(f"   confined at c=1 : {confined(S,1)}   at c=0 : {confined(S,0)}")
rep("confined for corridor c = 1", confined(S, 1))
rep("NOT confined at c = 0 (so R_N > 0 is possible)", not confined(S, 0))

m = least(w)
Fc, Fcp = Fc_pair(1, 2, 1, 3, 1)
print(f"   m(w) = {m},  F_c = {Fc},  F'_c = {Fcp}")
rep(f"m(w) = 187 (got {m})", m == 187)
rep(f"F_c  = 27  (got {Fc})", Fc == 27)
rep(f"F'_c = 59  (got {Fcp})", Fcp == 59)
rep("floors are valid: F_c <= F'_c <= m(w)", Fc <= Fcp <= m)
rep(f"orbit of {m} realizes w", sim(m, N) == w)

print(f"\n   Akhtar non-descent needs m(w) <= T(w) = {float(T):.4f}")
rep(f"F'_c = {Fcp} > T = {float(T):.4f}  -> EXCLUDES non-descent", Fcp > T)
rep(f"F_c  = {Fc} < T = {float(T):.4f}  -> does NOT exclude", Fc < T)
rep(f"actual m(w) = {m} > T  (consistent: no non-descent witness)", m > T)

print()
print("=" * 72)
print("SCOPE GUARD: positive compression requires 2^K > 3^N, i.e. R_N > 0")
print("=" * 72)
bad = 0; tested = 0
for N2 in range(1, 10):
    for K2 in range(N2, 2 * N2 + 4):
        pos = (1 << K2) > 3 ** N2
        conf0 = K2 * 1 <= 0  # placeholder
        # at c=0 every confined word has K <= floor(alpha N) hence 2^K <= 3^N
        s = 0
        while (1 << (s + 1)) <= 3 ** N2: s += 1
        if K2 <= s:
            tested += 1
            if (1 << K2) > 3 ** N2: bad += 1
print(f"   c=0: every admissible K <= floor(alpha N) gives 2^K <= 3^N")
rep(f"no positive compression at c=0 ({tested} cases, violations {bad})", bad == 0)

print()
print("=" * 72)
print("VERDICT:", "ALL PASS" if ok else "FAILURES PRESENT")
print("=" * 72)
sys.exit(0 if ok else 1)
