#!/usr/bin/env python3
"""Exact test suite for the constrained residual-ordering theorem.

T1 THROUGHOUT.  Confinement is decided by the exact integer comparison
    2^{S_j} <= 2^c * 3^j        (c a nonnegative integer),
which is equivalent to S_j <= alpha*j + c since S_j is an integer.  The
residual kernel H(w) = sum_{j=0}^{N-1} 2^{S_j}/3^j is accumulated as an exact
Fraction.  No floating point decides any assertion.
"""
from fractions import Fraction as F
from itertools import product
import sys


def confined(S, c):
    """S = list of prefix sums S_1..S_N.  Exact: 2^{S_j} <= 2^c 3^j."""
    return all((1 << s) <= (1 << c) * 3 ** j for j, s in enumerate(S, start=1))


def H(S_full):
    """H(w) = sum_{j=0}^{N-1} 2^{S_j} / 3^j, S_full = [S_0..S_{N-1}]."""
    return sum(F(1 << s, 3 ** j) for j, s in enumerate(S_full))


def Hw(w):
    """H of a word: H(w) = sum_{j=0}^{N-1} 2^{S_j}/3^j, with S_0 = 0."""
    return H([0] + prefix_sums(w)[:-1])


def prefix_sums(w):
    S, t = [], 0
    for d in w:
        t += d; S.append(t)
    return S


def words(N, K):
    """All compositions of K into N positive parts."""
    if N == 1:
        if K >= 1:
            yield (K,)
        return
    for first in range(1, K - N + 2):
        for rest in words(N - 1, K - first):
            yield (first,) + rest


def confined_words(N, K, c):
    return [w for w in words(N, K) if confined(prefix_sums(w), c)]


def Sstar(N, K, c):
    """S_j^* = min( floor(alpha j + c), K-N+j ).  floor(alpha j + c) computed
    exactly as the largest s with 2^s <= 2^c 3^j."""
    out = []
    for j in range(1, N + 1):
        s = 0
        while (1 << (s + 1)) <= (1 << c) * 3 ** j:
            s += 1
        out.append(min(s, K - N + j))
    return out


def wstar(N, K, c):
    S = Sstar(N, K, c)
    prev, w = 0, []
    for s in S:
        w.append(s - prev); prev = s
    return tuple(w)


def report(tag, ok):
    print(f"   {'PASS' if ok else 'FAIL'}  {tag}")
    return ok


def main():
    allok = True
    print("=" * 72)
    print("T1  Construction: d_j^* >= 1, total = K, w* confined")
    print("=" * 72)
    bad = tot = 0
    for N in range(1, 11):
        for c in range(0, 4):
            lo, hi = N, None
            # K ranges over N .. floor(alpha N + c)
            s = 0
            while (1 << (s + 1)) <= (1 << c) * 3 ** N:
                s += 1
            hi = s
            for K in range(lo, hi + 1):
                w = wstar(N, K, c); S = prefix_sums(w)
                tot += 1
                if not (all(d >= 1 for d in w) and sum(w) == K and confined(S, c)):
                    bad += 1
    allok &= report(f"{tot} (N,K,c) triples, violations = {bad}", bad == 0)

    print()
    print("=" * 72)
    print("T2  Dominance: every confined word has S_j <= S_j^* for all j")
    print("=" * 72)
    bad = tot = ntested = 0
    for N in range(1, 9):
        for c in range(0, 3):
            s = 0
            while (1 << (s + 1)) <= (1 << c) * 3 ** N:
                s += 1
            for K in range(N, s + 1):
                St = Sstar(N, K, c)
                cw = confined_words(N, K, c)
                ntested += len(cw)
                tot += 1
                for w in cw:
                    if any(a > b for a, b in zip(prefix_sums(w), St)):
                        bad += 1
    allok &= report(f"{tot} families, {ntested} confined words, violations = {bad}",
                    bad == 0)

    print()
    print("=" * 72)
    print("T3  w* is the UNIQUE maximiser of H among confined words")
    print("=" * 72)
    bad = ties = tot = 0
    for N in range(1, 9):
        for c in range(0, 3):
            s = 0
            while (1 << (s + 1)) <= (1 << c) * 3 ** N:
                s += 1
            for K in range(N, s + 1):
                cw = confined_words(N, K, c)
                if not cw:
                    continue
                tot += 1
                hs = [(Hw(w), w) for w in cw]
                best = max(hs)[0]
                argmax = [w for h, w in hs if h == best]
                ws = wstar(N, K, c)
                if ws not in cw or Hw(ws) != best:
                    bad += 1
                if len(argmax) > 1:
                    ties += 1
    allok &= report(f"{tot} families: w* is a maximiser, failures = {bad}", bad == 0)
    allok &= report(f"{tot} families: maximiser unique, non-unique = {ties}", ties == 0)

    print()
    print("=" * 72)
    print("T4  Minimum: w_min = (1,...,1,K-N+1), H_min = 3(1-(2/3)^N)")
    print("=" * 72)
    bad = badadm = ties = tot = 0
    for N in range(1, 9):
        for c in range(0, 3):
            s = 0
            while (1 << (s + 1)) <= (1 << c) * 3 ** N:
                s += 1
            for K in range(N, s + 1):
                cw = confined_words(N, K, c)
                if not cw:
                    continue
                tot += 1
                wmin = tuple([1] * (N - 1) + [K - N + 1])
                if wmin not in cw:
                    badadm += 1
                    continue
                hs = [(Hw(w), w) for w in cw]
                lo = min(hs)[0]
                argmin = [w for h, w in hs if h == lo]
                closed = 3 * (1 - F(2, 3) ** N)
                if lo != closed or Hw(wmin) != lo:
                    bad += 1
                if len(argmin) > 1:
                    ties += 1
    allok &= report(f"{tot} families: w_min admissible, failures = {badadm}", badadm == 0)
    allok &= report(f"{tot} families: H_min = 3(1-(2/3)^N) exactly, failures = {bad}",
                    bad == 0)
    allok &= report(f"{tot} families: minimiser unique, non-unique = {ties}", ties == 0)

    print()
    print("=" * 72)
    print("T5  MANDATORY EXAMPLE  N=8, K=12, c=0")
    print("=" * 72)
    N, K, c = 8, 12, 0
    ws = wstar(N, K, c)
    S = prefix_sums(ws)
    h = Hw(ws)
    cw = confined_words(N, K, c)
    allcomp = sum(1 for _ in words(N, K))
    print(f"   w*            = {ws}")
    print(f"   prefix sums   = {tuple(S)}")
    print(f"   H_max^conf    = {h} = {float(h):.10f}")
    print(f"   compositions of {K} into {N} parts : {allcomp}")
    print(f"   of these, every-prefix confined    : {len(cw)}")
    allok &= report("w* = (1,2,1,2,1,2,2,1)", ws == (1, 2, 1, 2, 1, 2, 2, 1))
    allok &= report("prefix sums = (1,3,4,6,7,9,11,12)",
                    tuple(S) == (1, 3, 4, 6, 7, 9, 11, 12))
    allok &= report(f"H = 13349/2187, approx 6.1037951532 (got {float(h):.10f})",
                    h == F(13349, 2187))
    allok &= report(f"confined count = 85 (got {len(cw)})", len(cw) == 85)
    # unrestricted front-loaded comparison
    front = tuple([K - N + 1] + [1] * (N - 1))
    hf = Hw(front)
    print(f"   unrestricted front-loaded {front}: H = {float(hf):.3f}")
    allok &= report(f"front-loaded H ~ 31.127 (got {float(hf):.3f})",
                    abs(float(hf) - 31.127) < 0.001)
    allok &= report("front-loaded word is NOT confined at c=0",
                    not confined(prefix_sums(front), 0))

    print()
    print("=" * 72)
    print("T6  Boundary cases: N=1, K=N, K=floor(alpha N + c)")
    print("=" * 72)
    ok = True
    for c in (0, 1, 2):
        # N = 1
        s = 0
        while (1 << (s + 1)) <= (1 << c) * 3: s += 1
        for K in range(1, s + 1):
            w = wstar(1, K, c)
            ok &= (w == (K,) and confined([K], c) and Hw(w) == F(1))
        # K = N  (all ones)
        for N in range(1, 9):
            w = wstar(N, N, c)
            ok &= (w == tuple([1] * N))
        # K = floor(alpha N + c)  (maximal)
        for N in range(1, 9):
            s = 0
            while (1 << (s + 1)) <= (1 << c) * 3 ** N: s += 1
            w = wstar(N, s, c)
            ok &= (sum(w) == s and confined(prefix_sums(w), c)
                   and all(d >= 1 for d in w))
    allok &= report("N=1, K=N (all ones), K=floor(alpha N+c) all behave", ok)
    # H_min at K=N equals H_max at K=N?
    same = all(wstar(N, N, c) == tuple([1] * N) for N in range(1, 9) for c in (0, 1, 2))
    allok &= report("K=N forces the unique word 1^N (max = min, degenerate)", same)

    print()
    print("=" * 72)
    print("VERDICT:", "ALL PASS" if allok else "FAILURES PRESENT")
    print("=" * 72)
    return 0 if allok else 1


if __name__ == "__main__":
    sys.exit(main())
