# Complete Q(√5) coverage with capped hub counts

The new internally reviewed computer-assisted proof candidate excludes
dimension 23 for compact arithmetic Coxeter polytopes with ground field
Q(√5). Together with the frozen exclusions of dimensions 24–29 and the
published compact ceiling 29, this gives **dimension at most 22** for that
field. Generation and saved replay of the
[binder](../scripts/fresh_sqrt5_hub_count_coverage_20261016.py) and its
[bound-22 certificate](fresh_sqrt5_hub_count_coverage_20261016_bound22.json)
pass. External mathematical validation and priority remain open.

The [new rule and separate review](fresh_sqrt5_hub_count_weights_20261016_review.md)
preserve all fixed short-path and H4 payments. Its additional bonuses use
the entire A1 count and the I2(5) count capped at two. Its exact vertex
maximum is 485211/10000 = 48.5211. The denominator in a cross angle remains
the actual source count, irrespective of the capped coefficient context.

## Direct coverage under the changed rule

The binder imports the frozen complete arithmetic classifications and
literal reduction maps from the
[previous coverage](fresh_sqrt5_dimension23_coverage_20261015.md). It
reevaluates every older sealed branch, singular first extension, singular
continuation, dotted stopping state, and complete finite-queue singular
state with the new rule. No previous numerical payment is assumed to
transfer. All 2,661 zero-touch rows remain in the output; 1,764 fit dimension
23 and receive 52,029 common spectator checks. The other 897 are explicitly
excluded at this target dimension by corner rank. The binder verifies
471 literal maps and 76 recursive manifests with 872 hash references.

A k-node core has external rank 25−k. Every common external pair `(e,f)`
of A1 and I2(5) counts is tested when `e+2f≤25−k` and the leftover rank is
zero or at least two. A2 and A3 components suffice to realize every leftover
rank at least two without changing either count. Each corner adds its own
core counts before applying the cap or either denominator. This is a
conservative enumeration; no arithmetic realization of each forest is
assumed.

Structural completeness is inherited from the independently reviewed
finite and dotted queues, their complete one-white domains and the earlier
conditional-sealing proofs. A first touching white root either creates a
persistent baseline payment or follows a retained branch. Orthogonal white
components are precisely the spectators above. All singular stopping
alternatives are checked, including expanded states. Nonsingular no-touch
states are rank-excluded. The justified 31-node limit covers every target
dimension at most 29, including 23.

## The strict contradiction

The 76 triangle markings still receive at least two. Every quadrilateral
receives at least one, and all larger-face weights are nonnegative. Thus
the average vertex weight is at least `5 f₂/f₀ − binom(23,2)`. The symmetric
nonnegative h-vector of a simple compact polytope gives `f₂/f₀>121/2`;
strictness follows from `h₀=h₂₃=1` and the larger extreme paired coefficient.
The average vertex weight is consequently strictly greater than
`5·121/2−253=99/2`. This contradicts the exact upper bound 48.5211.

The certificate explicitly combines this new exclusion with the complete
previous exclusions of dimensions 24–29. Other arithmetic fields and the
noncompact case are handled in the
[all-field composition](fresh_arithmetic_hub_dimension_composition_20261016.md).

```sh
python3 scripts/fresh_sqrt5_hub_count_coverage_20261016.py --dimension 23 --hub-cap 2
python3 scripts/fresh_sqrt5_hub_count_coverage_20261016.py --dimension 23 --hub-cap 2 --verify
```
