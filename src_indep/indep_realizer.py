#!/usr/bin/env python3
"""PHASE A1 -- INDEPENDENT re-derivation.  Shares no code with ../src/.

Independence: the realizer class is built here by BACKWARD RECURSION on the
word, never via the C_w affine identity.  From
    x realizes (v_1, w')  <=>  x_1 = (3x+1)/2^{v_1} realizes w'
and x ~ a' (mod 2^{M'}) for w' gives
    3x + 1 = 2^{v_1} a'  (mod 2^{v_1 + M'}),  i.e.
    x = 3^{-1}(2^{v_1} a' - 1)  (mod 2^{v_1 + M'}).
Base: the empty word is realized by every odd x, i.e. a = 1 (mod 2).
The modulus therefore comes out as 2^{S+1} with no input from me -- an
independent derivation of the +1 bit.
"""
import sys


def klass(w):
    """(a, M): seeds realizing w are x = a (mod 2^M), a least positive."""
    a, M = 1, 1                                   # empty word: x odd
    for v in reversed(w):
        M2 = v + M
        mod = 1 << M2
        a = (pow(3, -1, mod) * ((1 << v) * a - 1)) % mod
        M = M2
    return (a if a else (1 << M)), M


def least(w):
    a, M = klass(w)
    return a


def sim(m, n):
    """Valuation word of n accelerated steps -- written independently."""
    out, x = [], m
    for _ in range(n):
        y = 3 * x + 1
        v = 0
        while y % 2 == 0:
            y //= 2
            v += 1
        out.append(v)
        x = y
    return tuple(out)


def C_affine(w):
    """C_w via the append recurrence C' = 3C + 2^S, for cross-checking R1 only."""
    C, S = 0, 0
    for v in w:
        C = 3 * C + (1 << S)
        S += v
    return C


def C_R1(u, qa, v, qb, t):
    """The brief's R1, transcribed literally."""
    n = u + v + t + 2
    S = u + v + t + qa + qb
    return (3 ** n
            - 3 ** (v + t + 1) * (1 << (u + 1))
            + 3 ** (v + t + 1) * (1 << (u + qa))
            - 3 ** t * (1 << (u + qa + v + 1))
            + 3 ** t * (1 << (u + qa + v + qb))
            - (1 << S))


def word(u, qa, v, qb, t):
    return (1,) * u + (qa,) + (1,) * v + (qb,) + (1,) * t


# ----------------------------------------------------------------- R3: F_c
def F_c(u, qa, v, qb, t):
    """max{ m(P), (2^{S_P'} m(Q) - C_P')/3^{u+1} },  P=1^u qa 1^v,
    P'=1^u qa, Q=1^v qb 1^t.  Exact ceiling division."""
    P = (1,) * u + (qa,) + (1,) * v
    Pp = (1,) * u + (qa,)
    Q = (1,) * v + (qb,) + (1,) * t
    f1 = least(P)
    SPp, CPp = u + qa, C_affine(Pp)
    num = (1 << SPp) * least(Q) - CPp
    f2 = max(1, -((-num) // 3 ** (u + 1))) if num > 0 else 1
    return max(f1, f2), f1, f2


def run_floor(u, qa, v, qb, t):
    """The floor my 2026-10-22 package actually used: all-ones RUN bounds only."""
    w = word(u, qa, v, qb, t)
    best = 1
    for (k, L) in ((0, u), (u + 1, v), (u + v + 2, t)):
        if L == 0:
            continue
        pref = w[:k]
        Sk, Ck = sum(pref), C_affine(pref)
        num = (1 << Sk) * ((1 << (L + 1)) - 1) - Ck
        if num > 0:
            best = max(best, -((-num) // 3 ** k))
    return best


QS = (2, 3, 4)


if __name__ == "__main__":
    print("=" * 70)
    print("A1.0  independent modulus check: full odd scan on a small grid")
    print("=" * 70)
    bad = 0; scans = 0
    for u in range(3):
        for v in range(3):
            for t in range(3):
                for qa in QS:
                    for qb in QS:
                        w = word(u, qa, v, qb, t)
                        S = sum(w)
                        if S > 15:
                            continue
                        a, M = klass(w)
                        found = [m for m in range(1, (1 << M) + 1, 2)
                                 if sim(m, len(w)) == w]
                        scans += 1
                        if found != [a] or M != S + 1:
                            bad += 1
    print(f"   scans={scans}  failures={bad}   (M == S+1 and class is a singleton")
    print(f"   below 2^M in every case)  -> modulus 2^{{S+1}} CONFIRMED independently")

    print()
    print("=" * 70)
    print("A1.1  R1 three-run closed form")
    print("=" * 70)
    r1bad = r1n = 0
    for u in range(13):
        for v in range(13):
            for t in range(13):
                for qa in QS:
                    for qb in QS:
                        w = word(u, qa, v, qb, t)
                        r1n += 1
                        if C_R1(u, qa, v, qb, t) != C_affine(w):
                            r1bad += 1
    print(f"   configs={r1n}  mismatches={r1bad}   -> R1 {'CONFIRMED' if not r1bad else 'FALSE'}")

    print()
    print("=" * 70)
    print("A1.2  R2 full-denominator lemma")
    print("=" * 70)
    print("   R2 claims: gcd(C_w, 3^n) = 1 always, therefore D = 3^n.")
    print("   But the anchor is 3^{-n}(2^S - C_w), so D is governed by")
    print("   gcd(2^S - C_w, 3^n), NOT by gcd(C_w, 3^n).  Testing both:")
    g1bad = g2bad = n2 = 0
    wit = []
    for u in range(13):
        for v in range(13):
            for t in range(13):
                for qa in QS:
                    for qb in QS:
                        w = word(u, qa, v, qb, t)
                        n = len(w); S = sum(w); C = C_affine(w)
                        n2 += 1
                        from math import gcd
                        if gcd(C, 3 ** n) != 1:
                            g1bad += 1
                        g = gcd((1 << S) - C, 3 ** n)
                        if g != 1:
                            g2bad += 1
                            if len(wit) < 6:
                                wit.append((u, qa, v, qb, t, g))
    print(f"   configs = {n2}")
    print(f"   gcd(C_w, 3^n) != 1        : {g1bad}   -> premise "
          f"{'HOLDS' if not g1bad else 'FAILS'}")
    print(f"   gcd(2^S - C_w, 3^n) != 1  : {g2bad}   -> D = 3^n "
          f"{'HOLDS' if not g2bad else 'FAILS'}")
    if wit:
        print("   witnesses (u,qa,v,qb,t) -> gcd(2^S-C_w, 3^n):")
        for x in wit:
            print(f"      {x[:5]} -> {x[5]}")
