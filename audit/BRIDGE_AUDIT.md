# Two-Defect Bridge Audit — accelerated 3x+1 realizer floors

**Round:** falsification-first, 2026-10-22. **Working dir:** `~/scratch/two-defect-bridge/`
(the brief said `./two-defect-bridge/`; the session cwd was the `number-theory`
repo, so the round was run outside it to avoid polluting an unrelated project).

## VERDICT: **STOP**

No bridge exists in the class specified, and the reason is structural rather
than quantitative. The round also **refutes the null model as stated in the
brief** and **corrects the realizer modulus**, both by exact computation.

Everything below is exact integer arithmetic (Python `int`). `mpmath` is not
installed on this machine, so the brief's T2 float stage was replaced by exact
integers — strictly stronger. **No floating-point value decides any assertion
in this report.** The one estimate that is not exact (Phase 3 constant sizing)
is labelled T2/T3 and is parametric, because its constant was not verified at
source this round.

---

## 0. Corrections to the fixed context (both forced by computation)

### 0.1 The realizer class modulus is `2^{S+1}`, not `2^S`

The brief states the realizing seeds form one class mod `2^{S_w}`. It is
`2^{S_w+1}`. Exactly:

> `m` realizes `w` ⟺ `3^n m + C_w ≡ 2^S (mod 2^{S+1})`.

Two ingredients, both verified: (a) the single top divisibility
`2^{S_n} | 3^n m + C_n` already implies `2^{S_k} | 3^k m + C_k` for every
`k < n` — from `C_n = 3C_{n-1} + 2^{S_{n-1}}` — so it encodes every
"`v_2 ≥ v_k`" condition and makes each `x_k` (`k<n`) automatically odd;
(b) the one remaining condition is `x_n` odd, worth exactly one more bit.

Certified by exhaustive scan of **every** odd `m < 2^{S+1}` on 1104
configurations (`data/three_run_exhaustive.csv`): the realizer set is always
one class, and the modulus is `2^{S+1}` — not `2^S`, not `2^{S+2}`. **0
failures.** Reference identity `m(1^k) = 2^{k+1} − 1` holds for all `k < 40`.

**Cross-check against my prior two-run formula** (`src/two_run_crosscheck.py`,
528 configs). That formula — `A* = 3^{-(u+1)}c mod 2^{n+1}`, `r = 2^{u+1}A*-1`,
`n=q+v`, `c=1-2^{q-1}+2^n` — yields a **correct realizer in the correct class**
in 528/528 cases (orbit-simulated). But `A*` is reduced mod `2^{n+1}`, one bit
finer than the class requires, so it fixes a representative mod `2^{S+2}`:
`r ∈ {m(w), m(w)+2^{S+1}}`, and `r` is the **least** realizer only **257/528**
of the time. `log₂ r` therefore overstates `log₂ m(w)` by at most 1 bit —
exponent-level conclusions are unaffected, exact minima and champion seeds are
not. This is a representative-selection issue, not a wrong modulus.

### 0.2 The stated gap floor is **false**, and α belongs elsewhere

> `log₂ m(w) ≥ α·(longest defect-free run) − O(h)` — **refuted.**

It fails on **83125 / 357399 = 23.3%** of two-defect words with `|w| ≤ 62`
(`src/check_floors.py`). Worst case `w = 1^60 2 2`: `bl(m) = 61` but
`⌊α·60⌋ = 95` — **short by 34 bits**, unreachable by any `O(h)` with `h = 2`.
The failure is structural, not marginal: `m(w) < 2^{S+1}` always, while
`α·maxrun` can exceed `S+1` outright (it does so on 6055/34676 scatter rows;
smallest witness `u=0,q=2,v=0,q'=2,t=14`, where `α·14 = 22 > S+1 = 19`).

Where α actually belongs (`src/why_alpha_fails.py`). Under the conjugation
`y = x+1` an all-ones run is `y ↦ 3y/2`, so a run of `L` ones at `x_k` forces
`2^{L+1} | x_k+1`, i.e. **`x_k ≥ 2^{L+1} − 1` — coefficient 1**; and
`x_{k+L} = 3^L(x_k+1)/2^L − 1 ≥ 2·3^L − 1` — **coefficient α, on the orbit
maximum**. Transporting back to the seed via `m = (2^{S_k}x_k − C_k)/3^k`
cancels the `αL` against `−αL` exactly, leaving `L + 1 + (S_k − αk)`.
Measured: at `u=20`, `bl(seed) = 21 = u+1` while `bl(orbit max) = 33` against
`⌊α·20⌋+1 = 32`. **α is an orbit-maximum coefficient; on the seed it is 1.**

### 0.3 The corrected null model (used as the comparator throughout)

> `F(w) = max over all-ones runs of ⌈(2^{S_k}(2^{L+1} − 1) − C_k)/3^k⌉`

for a run of length `L` beginning after `k` letters. **0 violations in 357399
configurations.** This is the only admissible comparator: scoring a bridge
against the refuted α-floor would let it "win" against a bound that is not a
bound, which is exactly how the first Phase 2 run produced spurious margins.

### 0.4 Lemma H is vacuous on this family

`m ≥ 2^{S−1}/D` is non-trivial on only **220 / 111069** configs with `|w| ≤ 42`,
all at tiny `n`. Reason: `D = 3^{n−v₃(2^S−C_w)}` and two-defect words have
`S = n + h` with `h ≤ 6`, so `log₂(2^{S−1}/D) ≈ n(1−α) + h − 1 < 0` for
`n ≳ 10`. Denominator control needs `S ≫ αn` (high mean valuation); this family
has mean valuation `1 + h/n → 1`. Lemma H is a tool for a different family.

---

## 1. Phase 0 — exact three-run formula  [T1]

With `w = 1^u q 1^v q' 1^t`, `n = u+v+t+2`, `S = u+q+v+q'+t`:

```
C_w = 3^n + 3^{v+t+1}(2^{u+q} − 2^{u+1}) + 3^t(2^{u+q+v+q'} − 2^{u+q+v+1}) − 2^S

m(w) ≡ −1 − 3^{−(u+1)}·2^{u+1}(2^{q−1} − 1)
          − 3^{−(u+v+2)}·2^{u+q+v+1}(2^{q'−1} − 1)        (mod 2^{S+1})
```

`m(w)` is the least positive residue. Both defect terms carry `2^{u+1}` with an
odd cofactor, so **`2^{u+1} ‖ m+1` exactly** — the leading run is what forces
2-adic depth. Degenerate check: at `q=q'=1` the formula collapses to
`C = 3^n − 2^n`, the all-ones value.

**Verified** on all `u,v,t ≤ 12`, `q,q' ∈ {2,3,4}` (**19773** configs): closed
form vs general `C_w`, vs the class computation, vs direct orbit simulation of
the predicted seed, plus class stability at `m+2^{S+1}`. **0 failures.**
Composition law and prefix-nesting: 0 failures. (`data/three_run_verify.csv`.)

---

## 2. Phase 1 — null-model scatter  [T1, kill condition preregistered]

Kill condition was written into `src/null_model.py` before execution.

| regime | \|w\| | S+1 | bl(m) | null F | brief F | null/\|w\| | true/\|w\| | brief floor valid? |
|---|---|---|---|---|---|---|---|---|
| balanced (u=v=t) | 30 | 33 | 29 | 10 | 15 | 0.333 | 0.967 | yes |
| balanced | 60 | 63 | 63 | 20 | 31 | 0.333 | 1.050 | yes |
| balanced | 120 | 123 | 123 | 40 | 63 | 0.333 | 1.025 | yes |
| leading-heavy (v=t=1) | 30 | 33 | 33 | 27 | 41 | 0.900 | 1.100 | **no, short 8** |
| leading-heavy | 60 | 63 | 62 | 57 | 88 | 0.950 | 1.033 | **no, short 26** |
| leading-heavy | 120 | 123 | 123 | 117 | 183 | 0.975 | 1.025 | **no, short 60** |
| middle-heavy (u=t=1) | 120 | 123 | 121 | 117 | 183 | 0.975 | 1.008 | **no, short 62** |
| trailing-heavy (u=v=1) | 120 | 123 | 122 | 117 | 183 | 0.975 | 1.017 | **no, short 61** |

**Result: daylight exists** against the valid floor, so the preregistered rule
says *proceed* (`KILL(mech) = False`). But the daylight is sharply localised:

- **Single-run-heavy regimes: the null model is nearly tight** — rate 0.90–0.975
  against a truth of ~1.0, i.e. within 1–6 bits. Orbit-advance + truncation is
  essentially the whole story there.
- **Balanced regime: the entire gap.** 0.333 vs ~1.0, because the longest run is
  only `n/3` while the realizer still occupies `S+1` bits.

**What the daylight is.** The deficit `S+1 − bl(m)` has median **0–1** and
maximum ≈ `log₂(family size)` (14 at `|w| = 60` over 15930 configs; adversarial
minimum = maximal deficit). That is the pseudo-random/equidistribution profile,
and it reproduces the realizer-deficit law already established in my earlier
work (deficit ≈ geometric(1/2), mean ≈ 1). **The Phase 1 daylight is not new
information**, and the adversarial floor sits at `S+1 − log₂(#family)`, which is
an equidistribution statement, not a provable floor.

---

## 3. Phase 2 — bridge search  [T1]

### 3.1 Retraction (logged as it happened)

The **first** Phase 2 run reported margins up to **+30 bits** and 434
"interesting" configurations — over the BRIDGE-FOUND threshold. **Retracted.**
It computed δ by substituting the *measured* `m(w)`. But (ii) rearranges to

```
m(v) ≤ A·m(w)^q·2^{δ|w|}   ⟺   m(w) ≥ ( m(v) / (A·2^{δ|w|}) )^{1/q}
```

which is *literally the bridge's conclusion* with `K·2^{c|v|}` replaced by
`m(v)`. The run used the quantity being bounded to bound it. Script preserved
as `src/bridge_search_CIRCULAR_retracted.py` and output as
`data/phase2_candidates_CIRCULAR.csv` so the retraction is reproducible.

### 3.2 The Tautology Lemma (what replaced it)

> **Lemma.** Condition (ii) is never weaker than the floor it is meant to
> produce. Hence a bridge is non-circular only if (ii) is established by a
> mechanism that does not reference the *size* of `m(w)`.

In this framework exactly one such mechanism exists: a realizer of `w` also
realizes every contiguous **factor** of `w` — directly (prefix: `m(v) ≤ m(w)`)
or after orbit advance (factor at offset `k`:
`m(v) ≤ x_k = (3^k m(w) + C_k)/2^{S_k}`). So the admissibility test is
combinatorial, not statistical: **is the image word a contiguous factor of `w`?**

### 3.3 Candidate table

`a=1, b=1` for every length-preserving candidate. "margin" = provable floor
bits on `w` via the candidate, minus the null-model floor bits.

| candidate | image `v` | factor of `w`? | A, q, δ | max margin | where it died |
|---|---|---|---|---|---|
| (a1) merge, keep height | `1^u q 1^{v+t}` | **0 / 1548** | none exist | — | no non-circular (ii) |
| (a2) merge, additive height | `1^u (q+q'−1) 1^{v+t}` | **0 / 1548** | none exist | — | no non-circular (ii) |
| (b1) delete 2nd defect | `1^u q 1^{v+t}` | **0 / 1548** | none exist | — | no non-circular (ii) |
| (b2) delete 1st defect | `1^{u+v} q' 1^t` | **0 / 1548** | none exist | — | no non-circular (ii) |
| (c) suffix bridge *(main event)* | `1^v q' 1^t` | 1548 / 1548 | `A=1, q=1, δ\|w\|=α(u+1)−(u+q)` | **+0** | = null model exactly |
| (d1) reversal | `1^t q' 1^v q 1^u` | 174 / 1548 (palindromic only) | `A=1,q=1,δ=0` | **+0** | = null model |
| (d2) prefix control | `1^u q 1^v` | 1548 / 1548 | `A=1,q=1,δ=0` | **+0** | is the null model |

**Configurations with a non-circular (ii) and positive margin: 0.**
(BRIDGE FOUND needed ≥ 100.) Margin ranges where a factor relation exists:
`[−88, +0]` for (c), `[−89, +0]` for (d2), `[+0,+0]` for (d1) — **never positive.**

The merge and deletion images carry a genuine *potential* gain (≈ `0.415·L` bits
in the balanced regime, since `1^u q 1^{2L}` has a longer run than `w` does),
which is what the circular run was detecting. It is unobtainable: for these, (ii)
with `q=1` is **equivalent to an upper bound on the realizer deficit
`S+1 − log₂ m(w)`** — i.e. to the floor itself.

### 3.4 Candidate (c) resolved by an identity, not a measurement

The brief's main event asks whether the two one-defect subword anchors constrain
each other mod a common power of 2. They do not, and the reason is exact
(`src/anchor_checks.py`, 6561 configs, **0 deviations in every column**):

- the suffix congruence pulled back through `3^{u+1}m = 2^{u+q}x_{u+1} − C_{u+1}`
  has modulus `2^{S_k}·2^{S_s+1} = 2^{S+1}` — **it already *is* the full
  realizer congruence**;
- the prefix congruence is a **coarsening** of it (`2^{S_p+1}` divides `2^{S+1}`).

So the conjunction of the two subword conditions equals the suffix condition
alone. There is no independent information to combine through the shared run
`1^v`, and the best available bound is `max{prefix floor, drift-discounted
suffix floor}` — the null model, verbatim.

---

## 4. Phase 3 — effective calibration  [T1 structure, T2/T3 constants]

Triggered: Phase 2 found nothing, Phase 1 showed daylight.

**The reduction (exact, 0 violations in 19773 configs).** With
`A = (m+1)/2^{u+1}` odd,

```
Λ = 3^{u+v+2}·A + Z,   Z = 3^{v+1}(2^{q−1}−1) + 2^{q+v}(2^{q'−1}−1)
ord_2(Λ) ≥ M := q+v+q'+t = S−u
```

i.e. `ord_2( 3^{u+v+2}·A·(−Z)^{−1} − 1 ) ≥ M`: a linear form in **three**
logarithms, heights `(3, A, Z)`.

**Why two defects are categorically worse than one.** In the single-defect case
the third height is `2^{q−1}−1 ≤ 15`, so `log A₃ = O(1)`. Here
**`log₂ Z = Θ(v)`** — measured exactly equal to `⌊α(v+1)⌋+1` at `v = 5,10,20,40,80`
— and **`Z` is never {2,3}-smooth (19773/19773)**, so it cannot be absorbed into
the 2- and 3-adic data. A Matveev/Yu-class bound
`M ≤ C₃·log A₁·log A₂·log A₃·log b'` then gives
`log A ≥ M/(C₃·log3·log Z·log b')`, so the recoverable exponent fraction carries
an extra `1/v`.

Recoverable exponent fraction (parametric in `C₃`; 1.0 = full exponent, and the
null model already gives ≈ 1/3 in the balanced regime for free):

| regime | \|w\| | C₃=1e6 | C₃=1e9 | C₃=1e12 | C₃=1e18 |
|---|---|---|---|---|---|
| middle run bounded (v=2) | 10⁸ | 8.9e−09 | 8.9e−12 | 8.9e−15 | 8.9e−21 |
| balanced (u=v=t) | 100 | 5.8e−09 | 5.8e−12 | 5.8e−15 | 5.8e−21 |
| balanced | 10⁴ | 2.7e−11 | 2.7e−14 | 2.7e−17 | 2.7e−23 |
| balanced | 10⁸ | 1.3e−15 | 1.3e−18 | 1.3e−21 | 1.3e−27 |
| middle-run heavy (v=n−4) | 10⁸ | 6.5e−16 | 6.5e−19 | 6.5e−22 | 6.5e−28 |

In the balanced and middle-heavy regimes the fraction **tends to 0 as `|w|`
grows, for every fixed `C₃`** — the method loses the linear-in-`S` exponent
altogether, not merely constant factors. Only `v = O(1)` retains single-defect
quality, and there the word is a single-defect word plus a bounded perturbation.

**Two-log baseline, for scale.** My prior audited two-run work (Bugeaud–Laurent
1996 Théorème 4 / Corollaire 2 at `p=2`) gave
`n ≤ 267.4724·max{(q−1)ln2, ln A*, ln2}·M(u)²`,
`M(u) = max{ln(u+1)+0.889179, 5.545178}`, crossover `S = 32464`. **The
three-log constant `C₃` was NOT verified at source in this round**, which is
why the table above is parametric rather than asserting a crossover.

**The equivalent incomplete-window statement.** `m(w)` small ⟺ the residue
`A ≡ −3^{−(u+1)}(2^{q−1}−1) − 3^{−(u+v+2)}2^{q+v}(2^{q'−1}−1) (mod 2^{S−u})`
is small. Writing `ξ = 3^{−(u+1)}` in `(Z/2^{S−u})*`, the two coefficients are
`ξ` and `ξ·3^{−(v+1)}`. So the remaining gap is equivalent to a **pointwise**
bound on

> `#{(a,b) in a box : 3^{−a} and 3^{−a−b} lie in prescribed dyadic windows
> simultaneously}`

— a **two-dimensional** incomplete-window problem for the cyclic group generated
by 3 in `(Z/2^M)*`. Its one-dimensional case is exactly the single-defect
obstruction (ternary digits of `2^{−n}` mod `3^{u+1}`), for which
character/complete-sum methods were already ruled out (NO TRANSFER, 2026-09-21).
The two-dimensional case is strictly stronger, so that closure carries over
a fortiori — consistent with the brief's instruction not to re-propose them.

---

## 5. Retraction log

1. **Phase 2 first run — RETRACTED.** Margins to +30 bits, 434 "interesting"
   configs. Cause: δ computed from the measured `m(w)`, making (ii) and the
   conclusion the same inequality. Replaced by the factor test. Artefacts
   preserved for reproducibility.
2. **Phase 1 `rand_median_deficit` column — WITHDRAWN.** Computed `S+1` from the
   balanced config while the median came from a sample with different `q,q'`,
   hence the impossible value `−1`. Superseded by `src/check_floors.py`, which
   takes the deficit per record (median 0–1, max 14).
3. **Brief's gap floor — REFUTED** (§0.2), not retracted by me but recorded here
   because the round's kill condition was phrased against it.
4. **My own earlier two-run modulus — CORRECTED** (§0.1): it carried one
   redundant bit (`2^{n+1}` where `2^{S+1}` in `m` is the true class).

---

## 6. What this round established (written to be citable)

1. **The least-realizer class of a valuation word has modulus `2^{S+1}`**, via
   `3^n m + C_w ≡ 2^S (mod 2^{S+1})`; the top divisibility implies all prefix
   divisibilities, and the extra bit is `x_n` odd. Exhaustively certified.
2. **Exact closed form for three-run (two-defect) words** — `C_w` and `m(w)`
   above — with `2^{u+1} ‖ m+1`. Verified 19773 configs against direct simulation.
3. **The floor `log₂ m(w) ≥ α·(longest run) − O(h)` is false**, failing on 23.3%
   of two-defect words with shortfalls to 34 bits. **α is a coefficient on the
   orbit maximum; on the seed the coefficient is exactly 1**, because the drift
   factor `3^k/2^{S_k}` cancels it identically. The corrected floor `F(w)` has
   0 violations in 357399 configurations.
4. **Lemma H (denominator control) is vacuous on the two-defect family**, since
   its mean valuation tends to 1 while the lemma needs `S ≫ αn`.
5. **Tautology Lemma:** in the bridge formalism (i)+(ii), condition (ii) is
   logically a lower bound on `m(w)` of the same strength as the conclusion.
   A bridge is therefore non-circular only when (ii) comes from a factor
   relation, and the factor relations are exactly truncation and orbit-advance.
   **The bridge class collapses onto its own null model.** This extends, to the
   bridge formalism, the prefix-locality obstruction I had already recorded
   (a least realizer always realizes every prefix); it is not a new mechanism,
   but it is the first time the collapse is proved for the *transformation*
   class rather than for invariants.
6. **The two one-defect subword anchors of a two-defect word do not constrain
   each other.** The pulled-back suffix congruence already has the full modulus
   `2^{S+1}` and the prefix congruence is a coarsening (exact, 0/6561). The
   "combine two one-defect floors through the shared middle run" route is empty.
7. **Two defects are categorically harder than one for effective methods**, not
   just quantitatively: the third logarithm's height satisfies `log₂ Z = Θ(v)`
   (`= ⌊α(v+1)⌋+1`) and `Z` is never {2,3}-smooth, so the recoverable exponent
   fraction acquires a factor `1/v` and tends to 0 for every fixed constant,
   except when the middle run is bounded.
8. **The residual gap is a two-dimensional pointwise incomplete-window problem**
   for `⟨3⟩ ≤ (Z/2^M)*`, whose one-dimensional case is the already-closed
   single-defect obstruction.
9. **The Phase 1 daylight is the known realizer deficit**, not new information:
   median deficit 0–1, adversarial deficit `≈ log₂(family size)`. The null model
   is within 1–6 bits of the truth whenever one run dominates; the whole gap is
   in the balanced regime.

## 7. Scope and non-claims

- **Nothing here bears on the Collatz conjecture**, and nothing excludes
  divergent orbits. Every statement is about least realizers of finite
  valuation words.
- **No formal verification.** No Lean, no proof assistant. "Certified" above
  means exhaustive or exact finite computation, nothing more.
- **Validity / source verification / finite computation / effectiveness /
  novelty are kept separate.** Items 1–8 of §6 are proofs plus exact finite
  checks. Item 7's *constant* is an estimate (T2/T3) with the three-log constant
  unverified at source. Items 5–6 restate or extend obstructions I had already
  recorded; no priority is claimed for them, and absence from a literature
  search would not be evidence of novelty — none was performed this round.
- **No stronger conjecture is proposed**, and no claim is promoted beyond what
  the computation supports.

## 8. Reproducing

```
python3 src/three_run_exact.py        # Phase 0  -> data/three_run_{verify,exhaustive}.csv
python3 src/check_floors.py           # validates both candidate null models
python3 src/why_alpha_fails.py        # locates alpha on the orbit max
python3 src/null_model.py             # Phase 1  -> data/phase1_{summary,scatter}.csv
python3 src/bridge_search.py          # Phase 2  -> data/phase2_candidates.csv
python3 src/anchor_checks.py          # candidate (c) identity + Lemma H
python3 src/effective_three_logs.py   # Phase 3
python3 src/regime_table.py           # the table in section 2
python3 src/two_run_crosscheck.py     # prior two-run formula cross-check
```

Every number printed in this report is regenerated by one of these scripts.


---

## Appendix A — Phase A re-derivation (2026-10-22, later session)

An independent re-implementation (`src_indep/`, backward-recursion route,
no shared code) **confirms** the three-run closed form (0/19773) and
re-derives the `2^{S+1}` modulus from scratch. It also **corrects this
report on one substantive point**:

> **The null model used above is needlessly weak.** Sections 0.3, 2 and 3
> truncate only to *all-ones runs* (`m(1^L) = 2^{L+1}-1`). Truncating instead
> to the full **one-defect prefix** `P = 1^u q 1^v` and using its exact least
> realizer `m(P)` — still pure truncation, the same mechanism — gives a floor
> `F_c` that is **never worse and strictly better on 97.4%** of 19773
> configurations. In the balanced regime the median rate is **0.683, not the
> 0.333 reported in the section 2 table**. The Phase 1 "daylight" is therefore
> **0.683 vs 1.034**, roughly half the gap claimed above.

This does **not** overturn the STOP verdict: `F_c` is still a factor-relation
floor, so the Tautology Lemma (section 3.2) applies to it unchanged and no
bridge beats it. What changes is the *quantitative* characterisation of the
null model, and with it the "margin +0" readings in section 3.3 — those
compared a run-based floor against a run-based floor on both sides.

Full detail, plus the status of four further claimed results (R2 false as
stated; R5 not reproduced; R6 unavailable), in `notes/PHASE_A_REDERIVATION.md`.

### Appendix A.1 — reconciliation outcome (later same day)

Provenance was supplied for R3–R6 (a second, chat-session implementation of the
same prompt). Its package is **absent from this machine**, so its report and
`RECONCILIATION.md` are cited by title only, marked not verified. Outcome:

| item | outcome |
|---|---|
| R1 | confirmed, 0/19773 |
| R2 | **my Phase A objection withdrawn**: premise and conclusion both hold on the interior `u,v,t ≥ 1` (v₃ = 0 on all 19773); reductions confined to `t=0` (1176 vs 0); witness = the `x=1` cycle word `w=(2,2)`, `D=1` |
| R3 | confirmed; supersedes this report's null model (see Appendix A) |
| R4 | confirmed and sharpened: truncation threshold exactly 1/2 at every `k`; advance threshold exactly `α/2 = 0.792481` |
| R5 | **reconciled**: canonical image `m₂(u,q_a+q_b,v+t)` with δ by true `log₂` gives +0.3209 vs claimed 0.321 |
| R6 | **not carried**: underpowered at `n=100` (r ∈ [−0.117,+0.350] over 40 replicates); at `n=6000` r = +0.063/+0.092, nonzero |

Consolidated note: `notes/consolidated.md`, `papers/realizer-floors/main.pdf`
(5pp). Authoritative numbers: `src/verify_note_numbers.py` → `data/note_numbers.txt`.
