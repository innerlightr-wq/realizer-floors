# Exact Realizer Floors for One- and Two-Defect Valuation Words

Code and data for the technical note *Exact Realizer Floors for One- and
Two-Defect Valuation Words in the Accelerated 3x+1 Map* (E. De Jesús, 2026).

**The paper proves nothing about the Collatz conjecture** and excludes no
divergent orbits. Every statement concerns least realizers of finite valuation
words. There is no formal verification: "verified" means exhaustive or exact
finite computation on a stated grid.

The paper is in [`paper/`](paper/). Prose companion: [`consolidated.md`](consolidated.md).

## Setup

Python **3.14** (3.10+ should work), **standard library only** — no `numpy`,
`sympy` or `mpmath`. Nothing to install.

```sh
git clone https://github.com/innerlightr-wq/realizer-floors
cd realizer-floors
python3 src/verify_note_numbers.py       # regenerates every number in the paper
```

All arithmetic is exact `int`. The only floating point is `math.log2` applied to
exact integers, used for reporting ratios (marked **T2** in the paper); no
inequality is decided by a float.

## Conventions

Accelerated Syracuse map on odd `x`: `T(x) = (3x+1)/2^{v}` with `v = v₂(3x+1)`.
A valuation word is `w = (v₁,…,v_n)`, `S = Σvₖ`, `|w| = n`, `α = log₂3`.
`m(w)` is the least positive realizer. A two-defect word is
`w = 1^u q_a 1^v q_b 1^t` with `q_a, q_b ≥ 2`.

## Where each result in the paper comes from

| Paper | Claim | Script |
|---|---|---|
| Prop. 1 | realizer class is one class mod `2^{S+1}` | `src_indep/indep_realizer.py` (full odd scan, 243 words) |
| Prop. 2 | open vs closed-episode congruences | `src_indep/test_closed_episode.py` |
| Thm. 3 | two-defect closed form for `C_w` | `src_indep/indep_realizer.py` (0/19773 mismatches) |
| Lem. 4 | denominator `D = 3^{max(0, n−v₃)}` | `src_indep/test_R2_strata.py` |
| Lem. 4 | balanced lift, geometric(2/3) | `src_indep/test_R2_balanced_lift.py` |
| Thm. 5 (a),(b) | `F_c` validity and dominance | `src_indep/test_R3.py` |
| Thm. 5 (c) | closed-episode branch, 35.0% strict gain | `src_indep/test_closed_improvement.py` |
| **Prop. 6** | **unbounded gain, infinite family** | `src_indep/test_prop6_proof.py` |
| §6 | `β* = 1/2` exact; advance threshold `α/2` | `src_indep/test_R4_R5.py` |
| §7 | merge/deletion witnesses | `src_indep/test_R5_canonical.py` |
| §3 | two-run import cross-check (257/528) | `src/two_run_crosscheck.py` |
| all | single authoritative regeneration | `src/verify_note_numbers.py` → `data/note_numbers.txt` |

`src_indep/` is an **independent second implementation**: it builds the realizer
class by backward recursion on the word rather than through the affine identity,
so the modulus `2^{S+1}` falls out with no input from the author. The two
implementations agree on every shared quantity.

## Retractions and superseded files

Kept deliberately, so the record is checkable:

- `audit/bridge_search_CIRCULAR_retracted.py` — estimated its transfer constant
  from the measured `m(w)`, which is circular. See the paper's Provenance.
- `src_indep/test_R6_power.py`, `src_indep/test_R6_corr.py` — **conclusions
  withdrawn**; they measure a different statistic from the one they were taken
  to test. Do not cite their output.
- `audit/BRIDGE_AUDIT.md` — the separate two-defect *bridge* audit (verdict
  STOP) that preceded this note.

## Third-party inputs

`src_indep/test_R6_their_stat.py` re-runs a correlation defined in a separate
package that is **not redistributed here**. It exits with instructions if that
package is absent; point `CHAT_PKG` at it to run. Reported result: `r = −0.179`
(n=100), reproduced in this environment, with `r = −0.1673` at n=20000.

## Prior work

Proposition 1 is **not new here**: it was published by the same author as an
erratum to the Terras–Everett residue class
([10.5281/zenodo.21889704](https://doi.org/10.5281/zenodo.21889704), 11 Aug
2026), and appears independently in Akhtar
([10.5281/zenodo.22015248](https://doi.org/10.5281/zenodo.22015248), 19 Aug
2026). The affine identity is Böhm–Sontacchi (1978), Prop. 4. §9 of the paper
gives the full audit with verification status on each citation. The prior-art
check is **not exhaustive**, and absence from it is not evidence of novelty.

## Licence

Code (`src/`, `src_indep/`, `audit/*.py`): **MIT** — see `LICENSE`.
Paper and prose (`paper/`, `consolidated.md`, `audit/BRIDGE_AUDIT.md`):
**CC BY 4.0** — see `LICENSE-PAPER`.

## Citing

Cite the Zenodo deposit of the note. A DOI will be added here once deposited.
