# Independent review of the second completed H20 cusp

**PASS.** The [review source](../scripts/fresh_h20_second_cusp_completion_20261016_review.py)
and [certificate](fresh_h20_second_cusp_completion_20261016_review.json)
independently reconstruct the [151-wall chamber](fresh_h20_second_cusp_completion_20261016.md).
Generation and a fresh saved replay both pass. The primary report's
pending-review statement is historical; this separate report supplies the
review without changing that frozen source.

Starting from the independently reconstructed 146-wall geometry, the five
new normals are negative root sums plus positive multiples of the same
null vector nu2. Their heights are 3,3,3,5,6. The last two roots retain the
native squared lengths of walls 142 and 143. All 740 new pairs are checked:
102 right, six pi/3, two parallel and 630 dotted. A strict timelike center
and its facet projections certify all 151 facets are actual.

The complete 25-node incident diagram is

    affine A10 + affine A4 + affine A2 + 3 affine A1.

Exact rational matrices give six positive null kernels and total Euclidean
rank 19. Each kernel maps, by the physical root recipes, to a positive
multiple of the same future null ray. The incident physical normals span
nu2-perpendicular with rank 20. Therefore the horospherical link is a
bounded product of affine simplices. Together with the four preserved
old cusps, this certifies at least five proper cusps.

Every new wall is dotted to each of the 126 old walls outside this cusp
star. Any finite vertex using it would therefore require a positive
definite rank-20 block inside a positive semidefinite Gram matrix of rank
19, which is impossible. At least 171 old finite vertices remain, with no
new ones. The exact strict escape E5 has square 31500061/5000000000 and
negative pairing with every wall, so the chamber still has infinite volume.

There is also an unrestricted local obstruction. Any additional acute
spacelike normal r satisfies r.nu2<0. Positive affine relations first imply
r.nu2<=0. Equality would force r perpendicular to all 25 incident normals,
and hence proportional to the null vector nu2, contradicting spacelikeness.
This excludes further walls through the completed cusp, while allowing
walls away from it.

The review inherits the old isolated parameter and old chamber acuteness,
and does not claim a complete ambient cone, a complete cusp census, a
finite-volume example, or a new dimension. Equalities are justified by
exact recipes; interval overlaps serve only as consistency checks.

```sh
python3 scripts/fresh_h20_second_cusp_completion_20261016_review.py --verify
```
