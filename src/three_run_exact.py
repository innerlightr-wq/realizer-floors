#!/usr/bin/env python3
"""Phase 0 -- exact least-realizer formula for accelerated 3x+1 valuation words.

T1: exact integer arithmetic only (Python int).  No floats anywhere.

Conventions
-----------
Accelerated Syracuse map on odd x:  T(x) = (3x+1)/2^{v(x)},  v(x) = v_2(3x+1).
Word w = (v_1,...,v_n), S_k = v_1+...+v_k, S = S_n, |w| = n (letter count).
Affine identity:  2^{S} x_n = 3^n m + C_w,   C_w = sum_{k=1}^n 3^{n-k} 2^{S_{k-1}}.
Composition:      C_{uv} = 3^{|v|} C_u + 2^{S_u} C_v.
"""
import csv, os, sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(HERE, "data")


# ---------------------------------------------------------------- basic objects
def orbit_word(m, n):
    """Valuation word of the first n accelerated steps from odd seed m."""
    w, x = [], m
    for _ in range(n):
        y = 3 * x + 1
        v = (y & -y).bit_length() - 1
        w.append(v)
        x = y >> v
    return tuple(w)


def S_of(w):
    return sum(w)


def C_of(w):
    """C_w = sum_{k=1}^n 3^{n-k} 2^{S_{k-1}}, exact."""
    n, C, S = len(w), 0, 0
    for k in range(1, n + 1):
        C += 3 ** (n - k) * (1 << S)
        S += w[k - 1]
    return C


def C_compose(wu, wv):
    """Composition law check: C_{uv} = 3^{|v|} C_u + 2^{S_u} C_v."""
    return 3 ** len(wv) * C_of(wu) + (1 << S_of(wu)) * C_of(wv)


# ------------------------------------------------- the realizer class, derived
# m realizes w exactly  <=>  3^n m + C_w == 2^S  (mod 2^{S+1}).
#
# Proof of the two ingredients (both re-verified numerically below):
#  (1) 2^{S_n} | 3^n m + C_n  implies  2^{S_k} | 3^k m + C_k for every k<n,
#      because C_n = 3 C_{n-1} + 2^{S_{n-1}} gives
#      3^n m + C_n = 3(3^{n-1} m + C_{n-1}) + 2^{S_{n-1}}, and S_n >= S_{n-1};
#      hence the top divisibility alone encodes every "v_2 >= v_k" condition,
#      and x_k is automatically odd for k < n.
#  (2) the single remaining condition is x_n odd, i.e. one extra bit; so the
#      modulus is 2^{S+1}, NOT 2^S.
def m_class(w):
    """(m0, mod): realizer set = {m = m0 mod `mod`}; m0 is the least positive."""
    n, S, C = len(w), S_of(w), C_of(w)
    mod = 1 << (S + 1)
    inv = pow(pow(3, n, mod), -1, mod)
    m0 = (inv * ((1 << S) - C)) % mod
    return (mod if m0 == 0 else m0), mod


def anchor(w):
    """Lemma-H anchor N/D in lowest terms, with D a power of 3, plus S."""
    n, S, C = len(w), S_of(w), C_of(w)
    N, D = (1 << S) - C, 3 ** n
    while N % 3 == 0 and D % 3 == 0:
        N //= 3
        D //= 3
    return N, D, S


# ------------------------------------------------------- three-run closed forms
def three_run_word(u, q, v, qp, t):
    return (1,) * u + (q,) + (1,) * v + (qp,) + (1,) * t


def C_three_run_closed(u, q, v, qp, t):
    """C_w = 3^n + 3^{v+t+1}(2^{u+q}-2^{u+1}) + 3^t(2^{u+q+v+qp}-2^{u+q+v+1}) - 2^S."""
    n, S = u + v + t + 2, u + q + v + qp + t
    return (3 ** n
            + 3 ** (v + t + 1) * ((1 << (u + q)) - (1 << (u + 1)))
            + 3 ** t * ((1 << (u + q + v + qp)) - (1 << (u + q + v + 1)))
            - (1 << S))


def m_three_run_closed(u, q, v, qp, t):
    """m(w)+1 = -3^{-(u+1)}2^{u+1}(2^{q-1}-1) - 3^{-(u+v+2)}2^{u+q+v+1}(2^{qp-1}-1).

    Both defect terms carry 2^{u+1} with an odd cofactor, so 2^{u+1} || m+1
    exactly -- the leading run of u ones is what forces m+1 to be 2-adically deep.
    """
    S = u + q + v + qp + t
    mod = 1 << (S + 1)
    i1 = pow(pow(3, u + 1, mod), -1, mod)
    i2 = pow(pow(3, u + v + 2, mod), -1, mod)
    t1 = i1 * ((1 << (u + 1)) * ((1 << (q - 1)) - 1)) % mod
    t2 = i2 * ((1 << (u + q + v + 1)) * ((1 << (qp - 1)) - 1)) % mod
    m0 = (-1 - t1 - t2) % mod
    return (mod if m0 == 0 else m0), mod


# ----------------------------------------------------------------- verification
def exhaustive_class(w):
    """Every odd m in [1, 2^{S+1}] that realizes w, by direct orbit simulation."""
    n, S = len(w), S_of(w)
    return [m for m in range(1, (1 << (S + 1)) + 1, 2) if orbit_word(m, n) == w]


def main():
    os.makedirs(DATA, exist_ok=True)
    QS = (2, 3, 4)
    RMAX = 12

    # ---- A. closed form vs general formula vs direct simulation -------------
    rows, bad = [], 0
    for u in range(RMAX + 1):
        for v in range(RMAX + 1):
            for t in range(RMAX + 1):
                for q in QS:
                    for qp in QS:
                        w = three_run_word(u, q, v, qp, t)
                        n, S = len(w), S_of(w)
                        C_gen = C_of(w)
                        C_cl = C_three_run_closed(u, q, v, qp, t)
                        m_gen, mod = m_class(w)
                        m_cl, mod_cl = m_three_run_closed(u, q, v, qp, t)
                        realizes = orbit_word(m_cl, n) == w
                        nextone = orbit_word(m_cl + mod, n) == w
                        ok = (C_gen == C_cl and m_gen == m_cl
                              and mod == mod_cl and realizes and nextone)
                        bad += (not ok)
                        rows.append(dict(u=u, q=q, v=v, qp=qp, t=t, n=n, S=S,
                                         C=C_gen, m=m_cl, mod=mod,
                                         C_closed_agrees=int(C_gen == C_cl),
                                         m_closed_agrees=int(m_gen == m_cl),
                                         realizes=int(realizes),
                                         class_stable=int(nextone),
                                         all_ok=int(ok)))
    with open(os.path.join(DATA, "three_run_verify.csv"), "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=list(rows[0]))
        wr.writeheader()
        wr.writerows(rows)
    print(f"[A] three-run configs checked : {len(rows)}")
    print(f"[A] failures                  : {bad}")

    # ---- B. exhaustive modulus / minimality certification ------------------
    # Scans every odd m below 2^{S+1}; certifies BOTH that the realizer set is
    # exactly one class and that the modulus is 2^{S+1} (not 2^S, not 2^{S+2}).
    exh, ebad = [], 0
    for u in range(0, 5):
        for v in range(0, 5):
            for t in range(0, 5):
                for q in QS:
                    for qp in QS:
                        S = u + q + v + qp + t
                        if S > 17:
                            continue
                        w = three_run_word(u, q, v, qp, t)
                        found = exhaustive_class(w)
                        m0, mod = m_class(w)
                        pred = [m for m in range(m0, (1 << (S + 1)) + 1, mod)]
                        ok = (found == pred) and len(found) >= 1
                        ebad += (not ok)
                        exh.append(dict(u=u, q=q, v=v, qp=qp, t=t, S=S,
                                        n_realizers_below_2Sp1=len(found),
                                        least_found=found[0] if found else -1,
                                        least_pred=m0, modulus=mod,
                                        set_matches=int(found == pred)))
    with open(os.path.join(DATA, "three_run_exhaustive.csv"), "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=list(exh[0]))
        wr.writeheader()
        wr.writerows(exh)
    print(f"[B] exhaustive scans          : {len(exh)}")
    print(f"[B] failures                  : {ebad}")

    # ---- C. composition law + prefix nesting -------------------------------
    cbad = pbad = 0
    for u in range(0, 4):
        for v in range(0, 4):
            for t in range(0, 4):
                for q in QS:
                    for qp in QS:
                        w = three_run_word(u, q, v, qp, t)
                        cut = u + 1
                        cbad += (C_of(w) != C_compose(w[:cut], w[cut:]))
                        # prefix nesting: the full class refines every prefix class
                        m0, mod = m_class(w)
                        for k in range(1, len(w) + 1):
                            pm, pmod = m_class(w[:k])
                            if (m0 - pm) % pmod != 0:
                                pbad += 1
    print(f"[C] composition-law failures  : {cbad}")
    print(f"[C] prefix-nesting failures   : {pbad}")

    # ---- D. the two reference identities -----------------------------------
    d1 = all(m_class((1,) * k)[0] == (1 << (k + 1)) - 1 for k in range(1, 40))
    print(f"[D] m(1^k) == 2^{{k+1}}-1 for k<40 : {d1}")

    ok = (bad == 0 and ebad == 0 and cbad == 0 and pbad == 0 and d1)
    print("\nPHASE 0:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
