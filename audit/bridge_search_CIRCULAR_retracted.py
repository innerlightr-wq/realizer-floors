#!/usr/bin/env python3
#
# *** RETRACTED 2026-10-22 -- KEPT ONLY SO THE RETRACTION IS REPRODUCIBLE ***
# This script computes delta by substituting the MEASURED m(w).  But (ii)
# rearranges to  m(w) >= (m(v)/(A 2^{delta|w|}))^{1/q}, which IS the bridge
# conclusion.  So its "margins" assume a lower bound on m(w) of the same
# strength as the thing being proved.  Superseded by bridge_search.py.
#
"""Phase 2 -- bridge search in anchor coordinates.  T1 exact throughout.

A bridge w -> v needs
   (i)  |v| >= a|w| - b
   (ii) m(v) <= A m(w)^q 2^{delta|w|}
plus a PROVED floor m(v) >= K 2^{c|v|} on the image family, yielding
   m(w) >= (K 2^{-cb}/A)^{1/q} 2^{(ca-delta)|w|/q}.

Measured with q=1, A=1 (the only form any candidate here supports), so in bits
   yield(w) = [proved floor bits on v]  -  max(0, bl(m(v)) - bl(m(w)))
and the bridge is INTERESTING iff yield > [proved floor bits on w], uniformly
over a nonempty regime.  The proved floor on BOTH sides is the mechanism floor
F (Phase 1 certified it valid and certified the brief's alpha-floor invalid, so
F is the only admissible comparator).

delta must be UNIFORM over the family, so the uniform column takes
max(bl(m_v)-bl(m_w)) and min(floor bits) over the regime -- the worst case,
which is what a theorem would need.

PREREGISTERED KILL (declared before execution): a candidate is FAILED in a
regime if  uniform_yield <= uniform_comparator  there.  BRIDGE FOUND requires
uniform_yield > uniform_comparator on >= 100 configurations of some regime.
"""
import csv, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from three_run_exact import (three_run_word, m_three_run_closed, m_class,
                             C_of, S_of)
from null_model import mech_floor

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(HERE, "data")


def mech_floor_word(w):
    """Mechanism floor for an arbitrary word: max over maximal all-ones runs."""
    n = len(w)
    best, k = 1, 0
    while k < n:
        if w[k] == 1:
            L = 0
            while k + L < n and w[k + L] == 1:
                L += 1
            pref = w[:k]
            Sk, Ck = S_of(pref), C_of(pref)
            num = (1 << Sk) * ((1 << (L + 1)) - 1) - Ck
            if num > 0:
                best = max(best, -((-num) // 3 ** k))
            k += L
        else:
            k += 1
    return best


def mw(w):
    m, _ = m_class(w)
    return m


# ------------------------------------------------------------- the candidates
def cand_a1(u, q, v, qp, t):
    "(a) defect merge, keep first height: 1^u q 1^v q' 1^t -> 1^u q 1^{v+t}"
    return three_run_word(u, q, v + t, 1, 0)[:u + 1 + v + t]


def cand_a2(u, q, v, qp, t):
    "(a) defect merge, additive height: -> 1^u (q+q'-1) 1^{v+t}"
    return (1,) * u + (q + qp - 1,) + (1,) * (v + t)


def cand_b1(u, q, v, qp, t):
    "(b) delete SECOND defect, splice: -> 1^u q 1^{v+t}"
    return (1,) * u + (q,) + (1,) * (v + t)


def cand_b2(u, q, v, qp, t):
    "(b) delete FIRST defect, splice: -> 1^{u+v} q' 1^t"
    return (1,) * (u + v) + (qp,) + (1,) * t


def cand_c(u, q, v, qp, t):
    "(c) suffix bridge: the shifted suffix word 1^v q' 1^t"
    return (1,) * v + (qp,) + (1,) * t


def cand_d1(u, q, v, qp, t):
    "(d) reversal: 1^t q' 1^v q 1^u"
    return (1,) * t + (qp,) + (1,) * v + (q,) + (1,) * u


def cand_d2(u, q, v, qp, t):
    "(d) control: prefix truncation at the end of the middle run (the NULL model)"
    return (1,) * u + (q,) + (1,) * v


CANDS = [("a1_merge_keep", cand_a1), ("a2_merge_additive", cand_a2),
         ("b1_del_second", cand_b1), ("b2_del_first", cand_b2),
         ("c_suffix", cand_c), ("d1_reversal", cand_d1),
         ("d2_prefix_control", cand_d2)]

REGIMES = {
    "balanced":      lambda T: [(T // 3, T // 3, T - 2 * (T // 3))],
    "leading_heavy": lambda T: [(T - 2, 1, 1)] if T >= 4 else [],
    "middle_heavy":  lambda T: [(1, T - 2, 1)] if T >= 4 else [],
    "trailing_heavy": lambda T: [(1, 1, T - 2)] if T >= 4 else [],
}


def main():
    os.makedirs(DATA, exist_ok=True)
    QS = (2, 3, 4)
    rows = []
    for name, fn in CANDS:
        for rname, shapes in REGIMES.items():
            agg = {}
            for T in range(6, 91, 2):
                for (u, v, t) in shapes(T):
                    for q in QS:
                        for qp in QS:
                            w = three_run_word(u, q, v, qp, t)
                            vv = fn(u, q, v, qp, t)
                            if len(vv) == 0:
                                continue
                            m_w, m_v = mw(w), mw(vv)
                            Fw = mech_floor_word(w).bit_length()
                            Fv = mech_floor_word(vv).bit_length()
                            dbits = max(0, m_v.bit_length() - m_w.bit_length())
                            # drift correction for the shifted-suffix candidate:
                            # m >= (2^{S_k} m(v) - C_k)/3^k, k=u+1
                            if name == "c_suffix":
                                k = u + 1
                                Sk, Ck = S_of(w[:k]), C_of(w[:k])
                                num = (1 << Sk) * mech_floor_word(vv) - Ck
                                Fv = (max(1, -((-num) // 3 ** k))).bit_length()
                                dbits = 0
                            yld = Fv - dbits
                            a = agg.setdefault(T, dict(
                                n=u + v + t + 2, yld=[], cmp=[], d=[], Fv=[]))
                            a["yld"].append(yld); a["cmp"].append(Fw)
                            a["d"].append(dbits); a["Fv"].append(Fv)
            for T, a in sorted(agg.items()):
                uy = min(a["yld"])            # uniform (worst-case) yield
                uc = min(a["cmp"])            # what the null model already gives
                rows.append(dict(candidate=name, regime=rname, n=a["n"],
                                 uniform_yield_bits=uy,
                                 null_floor_bits=uc,
                                 margin=uy - uc,
                                 max_delta_bits=max(a["d"]),
                                 min_image_floor_bits=min(a["Fv"]),
                                 interesting=int(uy > uc)))
    with open(os.path.join(DATA, "phase2_candidates.csv"), "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=list(rows[0]))
        wr.writeheader(); wr.writerows(rows)

    print("=== Phase 2: uniform yield vs the null model it must beat ===")
    print("candidate            regime          max_margin  n_interesting/n  verdict")
    for name, _ in CANDS:
        for rname in REGIMES:
            sub = [r for r in rows if r["candidate"] == name
                   and r["regime"] == rname]
            if not sub:
                continue
            best = max(r["margin"] for r in sub)
            ni = sum(r["interesting"] for r in sub)
            verdict = "FAILED" if ni == 0 else f"CHECK({ni})"
            print(f"{name:20s} {rname:15s} {best:+10d}  {ni:5d}/{len(sub):<5d}"
                  f"      {verdict}")
    tot_int = sum(r["interesting"] for r in rows)
    print(f"\nconfigurations with uniform_yield > null floor: {tot_int}"
          f" / {len(rows)}   (BRIDGE FOUND needs >= 100)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
