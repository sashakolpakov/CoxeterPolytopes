# Third-vertex costs in the joint rank13 ADE triangle census

The complete [two-selected-vertex triangle census](fresh_sqrt5_rank13_high_ade_triangles_20261017.md)
has24,523 local occurrences. Reconstructing all three literal vertex forests
gives360 distinct third-vertex forests and3,558 distinct marked triangle
constraints, after identifying the two selected endpoints. These are complete
necessary local lists within the inherited215-type ADE selection, not assertions
of simultaneous realization.

The marked pair at the third vertex is the framing bond. There are3,411
occurrences with third component I2(5),10,278 with H3,8,131 with H4, and2,703
with I2(10). The established baseline charge is one at the first, second and
fourth kinds, but **two at the label-five bond of H4**. A triangle of the last
kind therefore needs no additional charge at its two selected ADE vertices.
For every other kind, one half at each selected vertex supplies the remaining
unit.

Counting a root pair only once gives these sharper capacities:

| Selected vertex forest | Previous possible triangle pairs | Pairs that still need one half |
|---|---:|---:|
| D13 |9|5|
| A13 |10|6|
| D5+2D4 |19|11|
| A1+4A3 |72|24|

The capacity decreases for214 of the215 selected forests. A pair admitting
both an H4 and a non-H4 third corner is retained in the half-charge set;
only H4-only pairs are removed. Thus this reduction does not assume which
one of several possible local triangles is actually present.

The [source](../scripts/fresh_sqrt5_rank13_triangle_third_costs_20261017.py)
reconstructs both selected endpoints and the full third forest from the
literal old graph, saved attachment and deleted root. It verifies their types,
marks, baseline charges and all occurrence counts. Its [saved certificate](fresh_sqrt5_rank13_triangle_third_costs_20261017.json)
includes all marked constraints, actual third-forest baseline totals, and every
selected vertex's refined corner set. Generation and a fresh saved replay pass.
This composes inherited classification engines; it is not an independently
implemented reproduction of the original triangle enumeration.

The accounting change matters: one cannot keep paying only one at every
nonspecial corner while using the reduced selected charges. A simple valid
upper rule pays one at every nonspecial root pair and an extra one at each
H4 label-five bond. A rank13 vertex then costs at most
`binom(13,2)+floor(13/4)=81`. The saved third-forest information permits
finer allocation than this coarse bound.

If a separate complete quadrilateral analysis gives the set J of possible
all-selected quadrilateral pairs, and I is the refined half-charge set here,
the selected-vertex cost is at most
`|I|/2+|J minus I|/4`. A root pair determines one two-face, so the triangle
and quadrilateral charges are not added at the same corner. This package does
not supply J or a joint dimension bound.

```sh
python3 scripts/fresh_sqrt5_rank13_triangle_third_costs_20261017.py
python3 scripts/fresh_sqrt5_rank13_triangle_third_costs_20261017.py --verify
```
