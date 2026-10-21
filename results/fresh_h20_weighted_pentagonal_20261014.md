# Pentagonal-radius caps and a four-family covering obstruction

The radius `R=4/phi` produces a new, exactly compatible mixed-radius motif: two such caps and two radius-one half-vertex caps cut all four probes of a golden square, with right and `pi/5` mutual contacts. The resulting 43-wall chamber has an exact strict escape.

The global compatibility test is stronger. Retaining the weighted-box 39-wall core, **no finite coherent union of the following four families has finite volume**:

- arbitrary seed-compatible radius-two caps;
- arbitrary seed-compatible radius-`4/phi` caps;
- radius-one caps centered at literal half-vertices;
- arbitrary-radius caps centered at literal box vertices.

Every finite seed label, parallel contact and dotted contact is allowed. At least **35,328 half-vertex probes are not strictly cut** by any such selection. The only possible boundary contacts at those probes admit a common outward perturbation, so at least one strict escape always remains. Arbitrary radius-one centers, other nonvertex radii, and changes to the retained core are outside this theorem.

The [exact checker](../scripts/fresh_h20_weighted_pentagonal_20261014.py) and [certificate](fresh_h20_weighted_pentagonal_20261014.json) pass generation and saved replay. A separate [Q(sqrt(5)) arithmetic and proof review](fresh_h20_weighted_pentagonal_20261014_review.md) also passes, including the new pair classes, all distance-two pairs, the motif and escape, and the covering count. The certificate depends on the frozen [weighted-box construction](fresh_h20_seed_deformation_20261014.md) and [all-label mixed covering gate](fresh_h20_weighted_mixed_probe_20261014.md).

## Complete center and cutter formulas

Write `phi=(1+sqrt(5))/2`. The box has ten unit inradii, one `sqrt(2)` inradius, and eight golden inradii `phi`. Its central sphere has radius two. A half-vertex probe has half of each signed box coordinate.

For `R=4/phi=4phi-4`, the smaller opposite-wall gap is at most the coordinate inradius. In a unit coordinate it must be zero, since `R/2>1`. In a heavy or golden coordinate it can only be zero or `R/2`; the next finite value, `R/sqrt(2)`, is already too large. The heavy `R/2` possibility fails because its complementary normalized gap lies strictly between `cos(pi/3)` and `cos(pi/4)`. Exact rational radical enclosures establish the strict interval.

Consequently every possible center has full signed unit and heavy coordinates. Each golden coordinate is either full, `±phi`, or small, `±(2-phi)`. This is an all-label classification, not a label-bounded scout.

If `k` golden coordinates are full, then

`|c|²=52-24phi+k*(4phi-4)`,

`central cosine=(phi-1)/2+k/4`.

The central gate rejects `k=0,1`; `k=2` gives `pi/5`, and every `k>=3` is dotted. The ball power at the matching half-probe is

`21-14phi+k*(4phi-6)`.

It is negative exactly for `k=2,3`. Any opposite coordinate sign makes the power positive, so a cutter has the probe's full sign pattern. There are exactly `C(8,2)+C(8,3)=84` cutters per probe. The two powers are `9-6phi` and `3-2phi`, respectively.

Their seed contacts are also explicit. A full unit coordinate gives right and `pi/5` contacts; a full heavy or golden coordinate gives right and dotted contacts; a small golden coordinate gives `pi/3` and `pi/5` contacts. All cap-cap and central normalized cosines lie in `Q(phi)`, so their only finite values are `0,1/2,phi/2`. The same field-degree argument used in the precursor makes the ensuing pair tests complete for all labels.

## A genuine mixed-radius motif

Fix the eleven unit/heavy signs and all but two golden signs. On golden patterns `00,10,01,11`, place:

| Probe pattern | Cap |
|---|---|
| `00` | radius one at its half-vertex |
| `10` | radius `4/phi`, full golden coordinates `{0,2,3}` |
| `01` | radius `4/phi`, full golden coordinates `{4,5,6}` |
| `11` | radius one at its half-vertex |

All four cut their own probes strictly. Their six mutual relations are four right, one `pi/5`, and one dotted. The checker verifies all seed and central relations and a common strict timelike point `(0,100,1)`. Projection onto each wall establishes that the core plus these caps has 43 actual facets.

This proves a useful finite-contact connection between the two radii. It is not a claim of finite-vertex or global covering closure. A literal escape has normalized spatial coordinates zero in the eleven unit/heavy positions and `99/100` in all eight golden positions. With `D=(|x|²-1)/2` and `U=1`, its Lorentz square is one and every wall inequality is strict. The certificate saves every new-cap power at this point.

## Why the checkerboard cannot extend globally

Consider two probes differing in one golden bit. The exact pair census is:

| Cutter families | Complete pair count | Compatible cases |
|---|---:|---|
| `R5/R5` | 7,056 | none |
| `R5/R2` | 7,392 | none |
| `R2/R2` | 7,744 | 1,260 right |
| `R1/R5` | 84 | 35 right, 21 at `pi/5` |
| `R1/R2` | 88 | 11 dotted |
| `R1/R1` | 1 | none |

Here `R5` means `4/phi`, and `R1` means the cap centered at the probe itself. The `R1/R5` cases require `k=3`; their relation is `pi/5` when the changed coordinate is full and right otherwise. The `R1/R2` cases require the radius-two cutter's unique full golden coordinate to be the changed bit. The inherited `R2/R2` rule likewise designates exactly one endpoint by its full golden coordinate.

The next-distance test supplies the key global restriction. For two `k=3` radius-`4/phi` cutters at probes differing in two golden bits `S`, compatibility is possible only at a right angle. Their full-coordinate triples `A,B` must be disjoint, and exactly one bit of `S` must lie in `A union B`. The checker verifies all 3,136 profile pairs: 240 right, 2,896 forbidden. An exact coordinate-distance identity transports this statement to every choice of `S`.

Thus three pairwise-distance-two probes cannot all receive compatible `k=3` cutters: their full-coordinate triples would be three disjoint three-element subsets of an eight-element set. This already rejects a complete `R1/R5` parity checkerboard on a three-cube. It also gives a uniform degree bound without assuming a checkerboard arrangement.

## The mixed covering bound and boundary cases

In each golden `Q8` fiber, let `A` be the probes strictly cut by any selected `R1`, `R2` or `R5` cap; choose one cutter for each probe. Every induced edge is of one of two kinds:

- An edge incident to a radius-two cutter is assigned to the unique endpoint designating that golden bit. Each radius-two cutter accepts at most one edge.
- An `R1/R5` edge is assigned to its radius-one endpoint. Its radius-`4/phi` neighbors are pairwise at golden distance two, so there are at most two of them.

Hence `e(A)<=2*n_R1+n_R2<=2|A|`. Cube degree counting gives `e(A)>=8|A|-1024`, so `|A|<=170`. Across 2,048 fibers these families strictly cut at most 348,160 probes.

The preceding exact vertex-radius gate shows that every seed-compatible vertex-centered radius is either four or has square at most eight. The latter caps miss all half-probes strictly. At most 512 radius-four caps are mutually compatible, and each covers 275 half-probes. Therefore at most another 140,800 probes are covered, leaving at least

`524288-348160-140800=35328`

probes not strictly cut.

There are no radius-two, radius-`4/phi`, or vertex-cap boundary exceptions. For radius two, setting the golden coefficient of a hypothetical zero power to zero leaves a strictly positive rational coefficient. For radius `4/phi`, the rational coefficient of every half-probe power is odd, so the power cannot vanish. The vertex cases were excluded exactly in the precursor.

A radius-one half-vertex cap can meet such a probe on its boundary only when its center differs in exactly one unit sign. At a missed probe, displace its ten unit coordinates by `epsilon` times their current signed half-values. The power of every tangent radius-one sphere increases by

`epsilon+(5/2)*epsilon² > 0`.

For sufficiently small positive `epsilon`, all other strict inequalities persist and the point stays inside the vertical box. A further small spacelike Lorentz perturbation then gives a strict escape. This handles boundary-only coverage and proves infinite volume for every finite coherent selection from the stated four-family union.

```sh
python scripts/fresh_h20_weighted_pentagonal_20261014.py --verify
```

The next construction must leave at least one of these cap templates or change the core. The new finite-contact motif is valid, but its local success does not overcome the global compatibility bound.
