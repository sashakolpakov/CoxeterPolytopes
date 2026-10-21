# A parity-controlled H20 cap orbit and its exact obstructions

The construction passes the exact Coxeter pair and actual-facet checks, but
**has infinite volume**. It strictly cuts every one of the 524,288 old
half-vertex probes; an explicit different escape survives. A separate local
certificate shows that none of the new binary caps contains a finite vertex.
Thus this is a tested construction mechanism, not a new finite-volume example.

The [producer](../scripts/fresh_h20_weighted_binary_20261014.py) and
[certificate](fresh_h20_weighted_binary_20261014.json) have a separate
[arithmetic review](fresh_h20_weighted_binary_20261014_review.md).
The [local producer](../scripts/fresh_h20_weighted_binary_20261014_local.py)
and [local certificate](fresh_h20_weighted_binary_20261014_local.json) have a
separate [graph and arithmetic review](fresh_h20_weighted_binary_20261014_local_review.md).
All primary generation, saved replay, and independent replay checks pass.

## Exact construction

Put `phi=(1+sqrt(5))/2`. Use normalized spatial coordinates with squared
weights `1` ten times, `2` once, and `phi^2=1+phi` eight times. The Lorentz
form is `|x|^2-2DU`, in signature `(20,1)`. Retain the 38 walls of the
box `|x_i|<=1` and the central sphere of radius two. A sphere with center
`c` and radius `R` has outward normal `(c,(|c|^2-R^2)/2,1)`; the retained
side is the exterior of its ball on the ideal boundary.

Let `A={0,1,2}` and `B={0,3,4}` among the ten unit coordinates. For each
of the `2^19` sign patterns, use the triple A if the sign parity among the
unit and golden coordinates is even, and B if it is odd. The heavy sign
does not enter parity. Put the three selected unit coordinates at signed
one and every other coordinate at signed one-half; use radius one.

For equal types the squared center distance is

```text
4 h_full + h_other + 2 h_heavy + (1+phi) h_golden,
```

with even unit-plus-golden difference parity. For different types it is

```text
1 + 4 h_intersection + 2 h_symmetric_difference
  + h_neither + 2 h_heavy + (1+phi) h_golden,
```

with odd parity. The 828 exact compressed rows, with binomial
multiplicities, cover all pairs by symmetry. Distinct pairs have squared
distance two, three, `2+phi`, or at least four: the relations are right,
`pi/3`, `pi/5`, parallel, or dotted. Each cap meets the seed in three right,
seven `pi/3`, one `pi/4`, and eight `pi/5` relations; its other seed
relations and central relation are dotted.

The matching half-vertex probe has ball power `-1/4`. This explains the
strict probe coverage and distinguishes these centers from the previously
excluded literal-half-vertex radius-one family.

Also include all `2^19` box-vertex-centered spheres of radius `sqrt(2)`.
Their mutual and seed relations are Coxeter, and every cross-pair with a
binary cap is dotted, since the squared distance is at least
`17/4+2phi > (1+sqrt(2))^2`.

The total is 1,048,615 distinct actual facets. The common timelike point
`(0,100,1)` is strict on all walls. Projecting it orthogonally onto any
one wall stays timelike and strict on every other wall, because all
distinct wall pairings are nonpositive. This proves the actual-facet
claim without enumerating a million-wall matrix. The 524,288 box-vertex
endpoints at height `sqrt(2)` survive, with positive-definite incident
diagram `20 A1`. This is not a complete finite-vertex census.

## Why the chamber is not finite-volume

Take all ten unit coordinates zero and all nine heavy/golden coordinates
equal to `99/100`. Its squared norm is
`(99/100)^2 (10+8phi)`. At this point the minimum binary, endpoint, and
central ball powers respectively are

```text
19/4 + (49/100)^2 (10+8phi) - 1,
10   + ( 1/100)^2 (10+8phi) - 2,
(99/100)^2 (10+8phi) - 4.
```

Each exceeds one. The point is strictly inside every coordinate bound.
Its norm-one Lorentz lift, with `U=1` and `D=(|x|^2-1)/2`, is consequently
strict for every wall and spacelike. This is an exact global escape and
proves infinite volume; passing the half-vertex probes was insufficient.

There is a second obstruction to improving this construction locally.
One reference binary facet has exactly 95 finite binary neighbors and
19 finite seed neighbors. The complete finite-pair graph on these 114
neighbors has 2,644 edges. A positive-definite incident diagram can
contain at most one `pi/5` neighbor of the reference root: two give a
negative three-by-three principal determinant. The 73 saved proper
colorings cover the zero-or-one possibilities and bound every such
clique, including the reference root, by 19. A finite H20 vertex requires
20 independent incident normals, so none lies on this facet. Explicit
sign and coordinate symmetries extend the conclusion to every binary cap.

This local certificate concerns the frozen chamber. Additional walls
that meet its binary facets require a new analysis. The global escape
above independently settles the current chamber's volume.

A subsequent [dotted-only completion theorem](fresh_h20_weighted_binary_20261014_dotted_completion_review.md)
also excludes every finite acute-compatible completion whose added walls
are all dotted to even one retained binary facet. Such walls leave that
facet unchanged. Its actual `pi/5` ridges would require finite vertices
in a finite-volume Coxeter polytope, contradicting the local certificate.
New walls must meet every retained binary facet, or the construction must
delete or replace the obstructed facets.

```sh
python3 scripts/fresh_h20_weighted_binary_20261014.py --verify
python3 scripts/fresh_h20_weighted_binary_20261014_review.py --verify
python3 scripts/fresh_h20_weighted_binary_20261014_local.py --verify
python3 scripts/fresh_h20_weighted_binary_20261014_local_review.py --verify
```
