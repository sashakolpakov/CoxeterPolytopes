# The integer-cap inner gap survives endpoint replacement

Changing the old endpoint radii does not rescue the radius-two integer-gap
cap family. Every cap in that family which is compatible with the central
radius-two cap misses `x=3z_000/20` by ball power at least `1/100`.
This result no longer assumes any particular endpoint layer.

The [exact checker](../scripts/fresh_h20_integer_inner_gap_20261014.py) and
[certificate](fresh_h20_integer_inner_gap_20261014.json) pass generation and
saved replay. The [independent review](fresh_h20_integer_inner_gap_20261014_review.md)
reconstructs all archived recipes with separate rational matrix arithmetic.

## Complete exact bound

The family consists of radius-two caps with nonnegative integer seed gaps
whose sums in the three affine factors are `7,7,8`. The factor composition
counts are `1716,1716,6435`. Their centers are uniquely determined by those
gaps. Central compatibility is the exact condition

```text
center norm squared = 8 or 12, or center norm squared >= 16.
```

Indeed the normalized central pairing is rational, so a finite angle can
only be right or `pi/3`; the other permitted relation is parallel or dotted.

At `x=3z_000/20`, every factor distance is nonnegative and has denominator
dividing 11,200. The checker enumerates every cap whose total power is at
most `1/100`, including equality. There are 651 such recipes. Exactly two
are central-compatible, both with power exactly `1/100`. Every omitted
recipe has larger power by the safe factor-distance filter.

Consequently the point

```text
E = (3z_000/20, 881/400, 1),       E² = 1/200
```

is strict on the central cap, all 22 seed walls and every compatible cap in
the entire integer-gap family. Its central power is `41/100` and its minimum
seed gap is `17/20`. Even allowing a cap of any radius at most 11 at **every**
old vertex leaves power at least `2061/100` there. The common timelike point
`(0,200,1)` is strict on these families. Every finite selection therefore
has infinite volume, even when the mutual cap-compatibility requirements
are relaxed. Cyclic symmetry transports the same statement to all 392
vertex directions.

This does not exclude other cap radii, noninteger gaps or nonvertex endpoint
caps. In particular, it is not a general obstruction to repairing the
central neighborhood.

## Scouts that led to the bound

The [segment scout](../scripts/fresh_h20_integer_segment_20261014.py)
enumerates all 8,066 integer-gap caps whose closed radius-two balls meet
`[0,z_000]`. It checks endpoint radius squares 1 through 21 with exact
integer contact gates. No tested bank starts a covering path beyond the
central ball on that segment. This is recorded as a numerical interval/DFS
result in its [archive](fresh_h20_integer_segment_20261014.json); the exact
inner-gap theorem above supplies the stronger, endpoint-independent reason.

The subsequent [noninteger inner-cap scout](../scripts/fresh_h20_inner_contact_20261014.py)
uses 19 prescribed seed contacts, mostly `pi/3`, solves the central-contact
quadratic, and transports every accepted profile to all 392 vertex
directions. It aims to cover the 392 probes `z/5` jointly.

| bounded profile bank | base recipes | transported occurrences | best greedy probe coverage |
|---|---:|---:|---:|
| central right; at most two seed labels changed to 4/5/6 | 1,255 | 491,960 | 4 of 392 |
| central labels 2/3/4/5/6 or parallel; at most one changed seed label | 344 | 134,848 | 3 of 392 |

Each bank used eight deterministic-seeded restarts. In these runs, the first
selected cap left an uncovered probe without any compatible remaining
candidate. This is **greedy failure, not infeasibility of either bank**.
The [first archive](fresh_h20_inner_contact_20261014_free.json) and
[second archive](fresh_h20_inner_contact_20261014_free_m2_3_4_5_6_0_d1.json)
retain the complete base recipes, selections, bounds and source hashes.
Neither scout constitutes a new exact chamber or a finite-volume certificate.

```sh
python3 scripts/fresh_h20_integer_inner_gap_20261014.py --verify
python3 scripts/fresh_h20_integer_inner_gap_20261014_review.py --verify
```
