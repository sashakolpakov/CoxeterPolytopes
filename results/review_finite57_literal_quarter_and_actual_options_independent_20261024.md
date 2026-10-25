# Independent source review: quarter and actual-option bridges

24 October 2026. Reviewer `/root/audit_finite57_algorithms`; implementation
author `/root/audit_sharp_ra`. I did not author these two supplemental
drivers. **No blocking defect found at their stated upper-bound contracts.**
Their actual execution remains pending root's resource admission. No import,
syntax check, test, worker, optimizing recurrence, Lean or TeX ran in this
review. The existing arithmetic draft and source ledger remained frozen.

## Exact source boundary

- `scripts/check_finite57_literal_quarter_binding_20261024.py`, 213 lines,
  SHA256 `3c4cbd7bc029cd189092b0989853ecf1a8de71b6cf8af2382348c130f0d4b594`.
- `scripts/check_finite57_literal_actual_options_20261024.py`, 119 lines,
  SHA256 `f76b84cdd1852c8ea544ab1becd95dd8c2050a088b4776aace03f5e09d15e494`.
- Common driver SHA256
  `b3780a3978fd3ac9bf5517b831437f2b304640209641a3c00b83872d79ce4d0f`;
  literal library SHA256
  `f6eba7af77af50b821e9e99c589388c1352f9ba5a99052aaeee64393f5ada1e0`.
- Contract:
  `results/finite57_literal_quarter_and_actual_option_contract_20261024.md`,
  SHA256 `73200e1ca26cd396cf2eec3eebdc054b70eb14dfab04a766df0238e4762662ed`.

I read every line of both drivers and the contract. The adapter, coefficient,
schema and frozen-input checks are recorded in the separate common-driver
review `review_finite57_literal_binding_driver_independent_20261024.md`.
I also read the relevant original producer definitions for `increasing_family`,
the conditional refinement tables and the complete actual `complement` loop.
Those producers were read as sources and were never imported.

## Quarter domain and the exact count intervals

The domain requires a fixed rank10 representative and residual rank at
least seven. This is necessary and sufficient to fit at least one A1,
A2 and A4 after the fixed components. The common domain already enumerates
every rank9--17 count, every A12 split and every rank5 count. It follows
that `expected_domain` covers every potentially triggering representative
class rather than choosing them from a saved maximizer list.

For residual rank R, the explicit count loops are exactly

    1 <= a1 <= R-6,
    1 <= a4 <= floor((R-a1-2)/4),
    1 <= a2 <= floor((R-a1-4*a4)/2).

Thus `count_options` sums the correct number of possible positive a2 values
for every a1/a4 pair. There is no missing `+1`: the a2 interval begins at
one. These are complete **actual** A1/A2/A4 counts, not merely reserved
lower counts. The index requires all 539 class keys and all 37985 raw
coarse options; equal set sizes prevent duplicate rows.

The coarse threshold is applied to every saved coarse row, giving exactly
14 selected classes. Actual fixed assignments are rebuilt by taking every
combinations-with-replacement product in each representative fiber, then
filtering for an actual A10. This retains mixed repeated rank10 types and
all other actual fixed types. Rank10 cannot be hidden in the residual
alphabet. The 77 selected multisets are therefore the complete actual
quarter fibers at the stated common-domain premises.

The index explicitly distinguishes raw fine options from reachable fine
options. It does **not** derive 23376 reachable options merely from a saved
scalar: that figure must equal the sum of the four executed fine-shard
counters. This is correctly stated in the contract.

## Coarse quarter inequality

Each coarse option uses the literal barred fixed/item coefficients and
the saved no-A2 residual upper. The explicit contextual terms are correct:

- `36*a2*(a2-1)` when at least two A1 are present covers the two-A2 guard
  with the actual A10 supplied by the quarter family.
- Under the rank9 representative guard, `64*a2*k5` and `40*a1*k5`
  cover the A2/D5 and second A1/D5 increments. Actual A9 and D5 form
  subsets of these representative classes. The first A1/D5 3/8 is
  already in the barred cross coefficient and is not added twice.
- `64*a2*a4` is the complete new A2/A4 increment: two A2 marks times
  two A4 endpoints times 1/4, for each component pair.
- `22*count(A17)` corrects the special barred base in its proper currency.

Every positive a1/a2/a4 choice is visited and compared with the frozen row
bound. The old no-A2 table is an upper relaxation: its formal residual may
contain extra A1 or A4. The proof uses the decomposition obtained by removing
**all actual counts** of those types first. It does not assert that an
undercharged formal reconstruction has the full literal context cost.
The common driver supplies current barred transition domination; inherited
saved-table upper validity remains essential. Coarse equality or actual
attainment is neither required nor asserted.

## Fine quarter inequality and the residual-context issue

Each selected actual fixed multiset is joined to its exact ancestor key
and row. The code independently reconstructs every ancestor fixed fiber,
not only the A10 subset. Its recorded residual alphabet must be exactly
the common full alphabet, with no duplicate type or weight omissions.
Every current literal ordinary item is checked **at most** its old recorded
weight, and residual ordinary interactions with A1, A2 and A4 are separately
zero. This is a domination comparison; it does not confuse the two epochs'
actual coefficients with equalities.

The remaining contextual proof is valid for the actual decomposition.
All actual A4 have been removed into their explicit count, so the actual
remainder has no A1, A2 or A4. The syntactic support check excludes every
other residual type capable of changing the context polynomial. A3 may
remain, but its only relevant D14 presence guard is already activated by
the mandatory A1. It follows that literal Psi is constant on every actual
remainder. Ordinary residual interactions vanish, and the coefficient
domination gives

    literal actual cost <= literal explicit-prefix cost + saved_upper[left].

The code checks that expression for **every reachable positive count triple**
and bounds it by its frozen refined row bound, itself at most 2861/8<377.
A saved table is allowed to contain extra A4 only as an upper-domain
relaxation. No literal-context identity is asserted for that formal extra
A4. The contract and return scope preserve this important distinction.

The four fine shards must total 77 fixed rows and 23376 reachable options.
A `None` residual value is a valid reason to skip only at the separately
assumed complete saved reachability/upper-table boundary.

## Every explicit option in the actual complement

The four `options` shards address the gap identified in the common driver.
They bind every table's exact allowed type set and literal item coefficient,
then enumerate all actual A1/A2 counts. The D14 guard is exhausted by
exactly these modes:

- A1 or A2 present, or no fixed D14: no A3 reservation is needed.
- Otherwise exclude A3 entirely, or reserve one actual A3 and leave the
  remaining A3 copies in the table alphabet.

The quarter exclusion is also exact: when A1, A2 and actual A10 are all
present, the complement contains no A4 and its table must exclude A4.
The shared complete fixed-fiber domain ensures actual A10 is accounted
for among the fixed components. No arbitrary representative is silently
substituted into the literal guard.

For each option the entire prefix is evaluated directly using the literal
function. It includes fixed bases, fixed-fixed cross costs, all A1/A2
ordinary terms, the reserved A3, and all Psi terms. Consequently the proof
does not need equality of independently exported native w1/w2 or fixed-cost
fields. Any such historical undercharge would be caught if it made a new
literal option exceed the recorded maximum. Residual item equality, the
current residual/residual zero premise and the flag support argument give
the required upper bound after this exact prefix.

Every required table key must exist. Duplicate table keys are rejected,
and `used==set(tabs)` rejects uncovered or unused stored flags. A key is
counted as used before a `None` reachability skip, consistently with the
producer's table-construction loop. The enumeration matches the complete
native loop's mathematical domain without importing its implementation.

All options are checked as inequalities against each frozen row maximum;
the older full literal attainer equality remains a separate lower-bound
check in the common driver. Combined with complete domains and valid
residual uppers, these are the two different directions needed for exact
maximal values. The supplements alone do not establish table validity.

Root must compare the four actual counters with 39 classes, 748 fibers and
96528 reachable options. The printed `complete_four_shard_totals` dictionary
is an expected contract, not itself an assertion about completed executions.

## Parsing, termination, resources and acceptance limits

Only the pinned common driver is imported, followed by its pinned literal
library and bounded adapters. All pins are checked before and after.
The exact JSON schemas match the frozen producer definitions: quarter rows
use `k5` and `key[9]`; selected refinements use `k5,a12`; ancestor classes
use `fixed_A5_count,fixed_A12_count`. Rational-string row bounds are
converted through Fraction; scaled saved table values remain integers.

The source loops are finite count/fiber/flag/rank loops, with no optimizing
recurrence. The 24 MiB selected-value and 40 MiB file boundaries and the
frozen-hash parsing limitations are unchanged. No general arbitrary-JSON
or heap-cap claim follows. Actual runtime and memory fit remain unmeasured.

The supported public entry rejects optimized Python and preserves the
bounded runner, parent/shared-lock checks, one thread, inner alarm, outer
59-second/1 GiB monitor, resource admission and abort floors. The new drivers'
own `assert __debug__` does not authorize an optimized direct route. Root
alone owns any actual launch after the unchanged 6 GiB memory admission.

These sources are suitable for resource-admitted **scoped bridge** batches.
Neither this review nor the existence of the source files supplies any
539/77/23376 or full-option runtime PASS. Local terminal-receipt bindings,
the validity of inherited computed tables, ideal geometric coverage,
layer59--65/higher-tail provenance and whole-manuscript verification remain
separate. No source, certificate, canonical manuscript or bound was changed.
