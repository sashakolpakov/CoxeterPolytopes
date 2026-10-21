# Independent source review: finite-volume tail scalar kernel

Date: 2026-10-20. Source-only mathematical review; no Lean, scientific,
TeX, test, or certificate job was executed by this reviewer. No canonical
source was edited.

Reviewed bindings:

* `formal/CoxeterProofs/FiniteVolumeTailKernel.lean`:
  `b8f03fd42fd60248379280d07949f596c21dfd6f16adc44c99060574b4c39daa`.
* `docs/finite_volume_tail_kernel_assumptions_20261020.md`:
  `fbc3bedadb9a7af5ec0998ed06311c2eaa041908e8e778f1245f0a16c02c02a7`.
* `manuscripts/finite_volume57.tex`:
  `c300db69c58211ffdf69c22970b289325224ad284d9e63a3c28d88821fca82a0`.

**Source-level mathematical composition: PASS.** Actual compilation and
axiom-print status are separate. This review identifies no reason to change
the theorem statements or arithmetic before the bounded compile/audit.

## Literal finite caps and indexing

The `u8` definition at Lean lines 19–30 exactly equals eight times every
finite cap printed in the paper's tables at 1963–1971, 2231–2240 and
2499–2507:

|Dimension|Paper cap U|Lean u8|
|---|---:|---:|
|59|783/2|3132|
|60|404|3232|
|61|420|3360|
|62|434|3472|
|63|449|3592|
|64|464|3712|
|65|479|3832|
|66|495|3960|
|67|504|4032|
|68|512|4096|
|69|544|4352|
|70–77|8n|64n|
|78–116|9n|72n|
|117–126|10n|80n|
|127–136|11n|88n|
|137–146|12n|96n|
|147–156|13n|104n|
|157–163|14n|112n|

Combining 78–116 is valid: the manuscript uses 9n on 78–80, 81–86,
87–91, 92–97, 98–108, and 109–116. The branch lengths are
`11+8+39+10+10+10+10+7=105`. There is no uncovered or duplicated dimension.

`beta8` is exactly the paper's parity formula multiplied by eight
(paper 1846–1849; Lean 33–34). Natural subtraction is ordinary subtraction
on the admitted range. `Fin 105` runs from i=0 through i=104, hence from
n=59 through n=163. The transfer theorem's `i=n-59` uses both range
hypotheses and reconstructs n exactly (Lean 38–47). It does not assume
the desired cap comparison. Ordinary kernel `decide` is the proposed
proof of the complete finite proposition; actual acceptance still
requires its successful compilation.

## Separate costs, denominator clearing, and payment subtraction

`sum_le_cap` derives a finite-list sum bound from each individual cap.
`cusp_sum_comparison` derives the comparison of the separate cusp-cost
and curvature sums from the comparison of each actual pair. Thus
`payment_contradiction` does not consume a supplied aggregate maximum.
Combining its hypotheses yields
`beta*finite.length + cuspCurvatures + correction
 <= finite.sum + cuspCosts
 <= beta*finite.length + cuspCurvatures`, contradicting positive correction.
It works with an empty cusp list and only weak cusp comparison, as the
paper's strict global correction already suffices (Lean 49–100).

The extra natural scale is necessary and correctly placed. Integrality
of 8U does not imply integrality of eight times each charge. After choosing
a positive integer s clearing all remaining denominators, the actual
common scale is 8s: finite cap `s*u8 n`, shape coefficient `s*beta8 n`, and
all list entries, deficit and correction in that same scale. Multiplication
by s preserves the cap comparison. The theorem permits s=0 as a formally
stronger statement; it does not assert that zero is a valid geometric
encoding. The assumptions document explicitly requires a positive
application scale (lines 20–28).

The Nat encoding also requires nonnegative encoded curvature and deficit.
Those are application obligations: cusp domination with nonnegative
charge gives nonnegative curvature, and the positive shape lower bound
gives nonnegative deficit in the hypothetical configuration. The formal
module makes no claim to encode arbitrary negative values of D5.

## Lorentz arithmetic and the boundary

The scalar meanings agree with paper 2597–2644:
`2p=MW-2H`, `H<=2M`, `9T=2N-2`, with N=n-1.

For the near range, `T-W<=4` implies `2N<=9W+38`. Multiplying by
nonnegative M and using `18H<=36M` gives
`18p=9MW-18H >= 2M(N-37)`, hence `9p>=M(N-37)`.
For the far range, the separate budgets give C0+C1<=4. Therefore
`3(T-W)<=4(C0+C1)<=16`, which implies `2N<=9W+50`, and the same Lorentz
substitution gives `9p>=M(N-43)` (Lean 113–137). Neither desired pairing
bound is assumed.

At N>=163 and M>=1, the near result gives `9p>=126`, hence p>=14.
In the far result, `21<=20M` implies `126<=120M`; nonnegativity of
`M(N-163)` gives `M(N-43)>=120M>=126`. Again p>=14, contradicting the
strict input p<14. In particular, **20M=21 is sufficient**, including
N=163, where the far lower bound reaches exactly 14 (Lean 140–165).

The final theorem makes an exhaustive order split on `20M<=21`. Its
complement implies `21<=20M`; equality goes through the near branch and
is also independently covered by the far theorem. The near upper ratio
enters through the conditional near-loss input, not through an unused
extra premise in the standalone scalar inequality (Lean 169–181).

## Scope and remaining application inputs

No conclusion `n<=58`, narrowness contradiction, or cap comparison is
smuggled in as a premise of the corresponding derived theorem. The generic
payment theorem intentionally takes `U<=beta`; its bounded successor
derives the literal comparison. The near/far theorems derive the pairing
bound and contradiction from the primitive scalar hypotheses.

The assumptions document correctly leaves actual simultaneous payment,
identification of the literal finite caps, cusp envelopes/concentration,
positive shape correction, framing-wall existence, Lorentz/common-height
interpretation, near loss<=4, and the four-thirds estimate as application
inputs. It does not claim to prove those difficult geometric/table facts.
The compact-facet alternative is absent. Consequently this is a derived
arithmetic kernel for the bounded tail and the unbounded **cusp branch**,
not a complete proof of the tail n<=58 or the final n<=57 theorem.

There are nine literal axiom-print commands, matching the document's
description. This source review does not report their outputs or upgrade
the draft's execution status.
