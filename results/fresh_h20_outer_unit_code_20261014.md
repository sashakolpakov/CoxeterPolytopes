# Outer code layers repair one escape but leave another

The exact construction retains the weighted 39-wall core and the entire binary radius-one orbit. Three new outer layers give **7,545,383 actual facets**. They cut the preceding escape, but the new chamber still has an explicit strict spacelike escape. Even with the heavy and golden coordinates fixed at one outer corner, their balls do not cover the unit ten-cube.

The [producer](../scripts/fresh_h20_outer_unit_code_20261014.py) and [certificate](fresh_h20_outer_unit_code_20261014.json) use exact rational arithmetic in Q(φ), with compressed sign-displacement counts rather than a matrix of millions of normals. The retained binary orbit and weighted core are bound to their [frozen certificate](fresh_h20_weighted_binary_20261014.json).

The [independent review](fresh_h20_outer_unit_code_20261014_review.md) passed generation and saved replay using separate Fraction arithmetic in the basis `1,√5`. It reconstructs all six pair-count tables by generating polynomials and checks the seed, central and binary bounds, preserved endpoints, old escape cut and new escape. The primary saved replay also passed.

## Three complete layers

Coordinates have spatial weights `1` ten times, `2` once, and `φ²=1+φ` eight times. The box is `[-1,1]^19`. All outer-layer heavy and golden coordinates are fully signed, independently. A sphere with center `c` and radius `R` has normal `(c,(c²−R²)/2,1)` of norm `R²` in the form `x²−2DU`.

Choose the ten distinct parity-check columns

`(0,1,2,3,4,5,8,10,12,15)` in F₂⁴.

Their xor is zero. The even zero-patterns with zero syndrome form a binary `[10,5,4]` code C of size 32. Its weight enumerator is `1+15z⁴+15z⁶+z¹⁰`. The 192 odd zero-patterns whose syndromes are not selected columns are exactly the odd patterns at distance three from C. Their enumerator is `60z³+72z⁵+60z⁷`.

| Layer | Unit coordinates | Radius squared | Number |
|---|---|---:|---:|
| Code | zero on a pattern of C; independently ±1 elsewhere | 2 | 1,139,200 |
| Holes | zero on one of the 192 odd patterns; ±1 elsewhere | 1 | 5,357,568 |
| Small | independently ±1/2 | 1/2 | 524,288 |
| Retained binary | the frozen parity-dependent triples | 1 | 524,288 |
| Retained core | 38 vertical sides and central radius-two sphere | — | 39 |

The zero codeword includes all original vertex-centered radius-√2 caps. Thus the new code layer contains that preceding layer; it does not discard any of its facets.

For two ternary unit centers with zero-patterns `x,y`, their unit squared distance is `|x△y|+4r`, where `r` counts opposite signs in their common nonzero positions. Equal-outer-sign pair relations are exhaustive:

| Pair | Distance values and relations |
|---|---|
| Code/code | 4 right, 6 π/3, 8 parallel, larger dotted |
| Hole/hole | 2 right, 4 parallel, larger dotted |
| Code/hole | 3 right, 5 π/4, larger dotted |
| Small/code | `5/2+2r`: right at r=0, parallel at r=1, then dotted |
| Small/hole | `5/2+2r`: π/4 at r=0, then dotted |
| Small/small | number of changed unit signs: 1 right, 2 parallel, then dotted |

Changing an outer sign adds 8 or `4φ²` to the squared distance. The sole additional parallel case is code/code with identical unit coordinates and only the heavy sign changed. All other such pairs are dotted. The certificate counts every unordered same-layer pair and every cross-layer pair with exact multiplicities.

Code zeros meet both unit sides at π/4; hole zeros meet both in parallel. Small centers meet the nearer unit side at π/4. Nonzero full coordinates give a right contact on one side and a dotted relation on the other. All three layers are dotted to the central cap.

They are also **strictly dotted to every retained binary cap**. The squared-distance lower bounds are `17/4+2φ` for code/hole centers and `13/4+2φ` for small centers, strictly larger than the respective external-tangency thresholds. A zero unit coordinate only increases the first bound at a binary cap's full coordinates. This argument permits both binary triples and all signs, so it does not require any unproved transitivity assumption.

The timelike point `(0,100,1)` is strict for every wall. Projection to any one wall, using all nonpositive pairings, proves that every listed wall is an actual facet. The original `2^19` finite endpoints `(z,(z²+2)/2,1)` remain actual finite vertices: their nineteen orthogonal vertical walls and original vertex cap are incident, and every other new wall is strict. This is a preservation statement, not a complete new vertex or cusp census.

## Repaired and remaining gaps

At the old escape's spatial point, with all unit coordinates zero and all outer coordinates `99/100`, the code cap with zero-pattern `1023` has power

`−1999/1000+φ/1250 < 0`.

Thus it cuts the old strict spacelike escape as well as its spatial ideal point.

For a new uncovered point, use the dual-code support `{1,3,5,9}`. Set unit coordinate 1 to zero, coordinates 3,5,9 to one, and the other six to one half; set every outer coordinate to one. The specified zero-pattern on those four coordinates has odd parity, while every codeword has even parity there. Consequently the minimum code-ball power is 1/2. Hole-ball and small-ball powers also have minimum 1/2. This is an exact failure of coverage of the full unit cube at fixed outer signs, not just a sampled gap.

Now move the three unit ones and all nine outer coordinates to `99/100`, leaving the zero and halves unchanged. Call the resulting spatial point x. The exact minimum powers are:

| Family | Power at x |
|---|---|
| Code | `4813/10000+φ/1250` |
| Holes | `5013/10000+φ/1250` |
| Small | `4713/10000+φ/1250` |
| Binary, uniform lower bound | `26313/10000+(2401/1250)φ` |
| Central | `102413/10000+(9801/1250)φ` |

All exceed 1/4, and every vertical side is strict. Hence

`X=(x,(x²−1/4)/2,1)`

is a future strict spacelike point of norm 1/4 for the entire chamber. Its exact coordinates and the minimum-attaining code/hole patterns are saved. This proves infinite volume.

There is a stronger reason to stop adding these particular outer layers: they are all dotted to the binary orbit. The separate [dotted-only completion theorem](fresh_h20_weighted_binary_20261014_dotted_completion_review.md) rules out a finite-sided finite-volume completion obtained solely by additions dotted to every binary facet. The present package supplies a concrete instance; it makes no obstruction claim about new walls that meet binary facets at finite angles or about coherent replacements.

Replay:

```sh
python scripts/fresh_h20_outer_unit_code_20261014.py --verify
```
