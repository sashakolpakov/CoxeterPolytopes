# Integer-gap cap obstruction: independent review

**PASS.** The [primary checker](../scripts/fresh_h20_balanced_integer_obstruction_20261014.py) and [certificate](fresh_h20_balanced_integer_obstruction_20261014.json) passed a fresh saved replay. A separate [Fraction matrix/witness checker](../scripts/fresh_h20_balanced_integer_obstruction_20261014_review.py) and [record](fresh_h20_balanced_integer_obstruction_20261014_review.json) passed generation and replay. No mathematical issue was found.

The independent checker inverts the literal unit A6/A7 Grams, reconstructs all 1,110 archived cutter recipes, verifies their norms and probe powers, and checks 218,115 endpoint relations, including every first incompatible endpoint witness. It imports no primary arithmetic. Completeness of the composition search was reviewed at source level and exercised by the fresh primary replay; the independent script is not a second exhaustive enumeration.

The all-angle gate is valid. For a radius-two cap against an endpoint sphere of radius `sqrt(28)`, the normalized negative pairing is

```text
c = (distance²−32)/(8sqrt(7)).
```

Its square is rational. If `c=cos(pi/m)`, then `2cos(2pi/m)=4c²−2` is a rational algebraic integer in `[-2,2]`. For finite `m≥2` this permits only `c²=0,1/4,1/2,3/4`. The three nonzero possibilities require a rational square equal to `7/4`, `7/2`, or `21/4`, impossible. Thus the only finite angle is right; the remaining compatible cases are parallel or dotted. The producer correctly checks the sign before squaring the parallel/dotted threshold. Against the radius-two central cap, the cosine itself is rational, so finite contacts are right or `pi/3`, giving center square 8 or 12; parallel/dotted starts at 16.

Every nonnegative integer seed-gap vector with factor sums `7,7,8` occurs in the weak-composition enumeration. Each factor contributes a nonnegative squared distance to the target, so discarding a factor whose contribution already exceeds radius squared four cannot lose a cutter or tangent. The remaining Cartesian product is checked exactly, including equality. There are 555 strict cutters at each probe; all 555 are baseline-compatible at `2z/5`, while none is baseline-compatible at `3z/5`. There are no tangent recipes at either probe.

The rational Gram inverse gives denominator dividing 28 for every center norm and center–vertex product. Consequently all powers at `3z/5` lie in `(1/700)Z`. Excluding nonpositive power therefore gives a uniform lower bound `1/700` for the **entire** baseline-compatible family. The lift

```text
E = ((3/5)z_000, 98783/2800, 1),       E²=1/1400
```

is strict on all these caps and on the retained 415-wall baseline. A common timelike strict point exists at `(0,200,1)` for every such radius-two cap, since its center norm is nonnegative. Thus every finite selection has infinite volume, even without checking mutual cap compatibility.

The theorem concerns fixed radius two, nonnegative integer seed gaps, and the retained complete radius-`sqrt(28)` endpoint layer. It does not exclude noninteger gaps, other radii, or replacement of endpoint facets.

```sh
python scripts/fresh_h20_balanced_integer_obstruction_20261014.py --verify
python scripts/fresh_h20_balanced_integer_obstruction_20261014_review.py --verify
```
