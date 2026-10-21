# The missing B1 face payments in the global bound-20 proposal

The [exact screen](../scripts/fresh_sqrt5_global_b1_first_extensions_20261016.py)
and [saved data](fresh_sqrt5_global_b1_first_extensions_20261016.json) evaluate
the frozen [global weighting proposal](fresh_sqrt5_global_discount_20261016.md)
on all 745 inherited baseline-one minimal cores and their 3,820 retained first
white extensions. Generation and a fresh saved replay passed.

Every common spectator context at total corner rank 21 is included: 3,067,805
cases, yielding 97,711 distinct exact payment profiles. Of the 4,565 rows,
918 receive less than one in at least one context. These comprise 716 finite
first extensions, 122 symbolic dotted first extensions, 61 extensible minimal
cores and 19 already rejected sealed nonsingular minimal cores. No sealed
singular minimal core fails. The minimum payment is exactly one minus the
adjacent-edge discount, approximately 0.770363341946.

The screen deliberately includes partial or rejected cores. Its 918 failures
are therefore not 918 actual compact faces. A nonsingular core needs further
touching roots before it can be a complete quadrilateral star. Appending
orthogonal roots alone does not fix that rank obstruction.

The separate [finite first-extension census](fresh_sqrt5_global_b1_finite_first_20261016.md)
retains all singular stopping choices. The [B6 witness](fresh_sqrt5_global_b6_discount_witness_20261016.md)
then verifies an actual compact local quadrilateral, with admissible arithmetic
rank-21 spectator data, whose payment is below one. This establishes a real
local coverage failure of the restricted global proposal, while making no
claim of a compact ambient realization.

The [dotted completion screen](fresh_sqrt5_global_b1_dotted_first_completions_20261016.md)
addresses the 122 underpaid symbolic first extensions separately. Further
touching extensions and distance domains must also be covered. No dimension
bound follows from these screens. The older compact arithmetic bound-21
candidate uses a different rule and remains unchanged.

```sh
python3 scripts/fresh_sqrt5_global_b1_first_extensions_20261016.py --verify
```
