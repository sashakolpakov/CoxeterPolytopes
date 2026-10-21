# Independent review of the marked short-face incidence bound

**PASS:** on the selected narrow facet, let `H` count the two marked
baseline-38 vertex types and `L` all other vertices. Then
`34H<155L`, so their fraction is less than `155/189`. The inherited baseline
average is strictly below `7148/189`, approximately 37.8201. This remains
above the required 35 and gives no new dimension ceiling.

The [review source](../scripts/fresh_sqrt5_marked_short_faces_20261016_review.py)
and [certificate](fresh_sqrt5_marked_short_faces_20261016_review.json) import
neither the [primary producer](../scripts/fresh_sqrt5_marked_short_faces_20261016.py)
nor the transition producer. They use the independent Golden arithmetic,
matrix operations and distinct-color tree codes from the earlier marked
flow review. The finite connected-diagram catalogue is an explicitly
inherited, independently reviewed dependency.

## Complete transitions are reconstructed

The review re-enumerates all 54,040 subsets of up to three attachment roots
over all 40 literal ports of the two maximal marked types. Tree codes retain
the facet root and exactly distinguish the two targets. The full 66 accepted
attachment list agrees with the primary. The common positive Gram and its
Schur data are independently reconstructed for every port. A unit-bond path
from the old framing root to a new attachment makes the cyclic product
equal to plus or minus the framing product, establishing its integrality.

The scalar check uses an enlarged box `-5<=a<=1`, `-7<=b<=-1` for
`g=a+b*phi`, justified by minimality 14 and conjugate definiteness. Every
admissible framing scalar and attachment agrees with the complete transition
certificate. Thus the short-face review does not silently rely on the
earlier review of only 20 selected edge witnesses.

## The fourth pairing is free

Both maximal marked types and all cross transitions are included. The
review reconstructs all 3,028 pairs of edge choices. For 3,010 pairs the
Schur determinant is strictly positive at both ends of `[-2,0]`; its
concavity makes it positive throughout. Both residual diagonals are
negative, so every such putative finite closure has two negative directions
and cannot occur in the ambient Lorentz space.

For the remaining 18 pairs, a finite Coxeter contact has normalized cycle
in `Q(sqrt5)`. The real cyclotomic degree condition gives `phi(m)<=4`.
The possible labels are initially `2,3,4,5,6,8,10,12`: the prime factors
of `m` can only be 2, 3 and 5, with the displayed bounded exponents.
Labels 8 and 12 would require square roots of 2 and 3 in `Q(sqrt5)` and
are excluded. For every remaining pair the cross projection `c` is nonzero.
The corank equation

`ad-(x-c)^2=0`

and `x^2 in Q(sqrt5)` imply `x=(x^2+c^2-ad)/(2c)` is in that field.
Products for labels 4 and 6 are excluded by rational square tests; label 10
has squared product `2+phi` of norm 5, which cannot be the square of a field
element. Only products zero, minus one and minus phi remain. The resulting
54 exact tests have just one rank root, giving the reported quadrilateral
at the end corner of `A1+D20` with all four vertices maximal.

Direct Gram multiplication verifies its root dependence, four future finite
rays and every wall inequality. The fourth vertex is checked with a marked
tree code. The conclusion allows an arbitrary finite fourth vertex, so a
quadrilateral with three maximal vertices must have four; there is no
unexamined nonmaximal fourth-corner option.

## Triangles and the incidence argument

An all-maximal triangle would use the same new wall for two continuations
from one maximal vertex. At a released root its product would have to be
both a simply laced endpoint contact (zero or minus one) and a maximal-edge
framing product (minus phi or dotted). This is impossible.

For a triangle with exactly two maximal vertices, the maximal edge must
have framing minus phi; dotted framing is incompatible with the third
finite vertex. The reviewer exhausts all such edges and all other released
roots, using its own finite graph codes. Exactly three directed cases
remain, matching the primary. It directly verifies all nine finite future
rays and recomputes the third vertices' baseline 35. At each maximal vertex,
these triangles and any all-maximal quadrilateral use the same single end
corner.

Writing `T2` for triangles with two maximal vertices and `Q4` for
all-maximal quadrilaterals, simplicity therefore gives

`2T2+4Q4<=H`.

Every other triangle and quadrilateral has at least two nonmaximal vertices.
Its positive face excess (two or one, respectively) can be assigned with at
most one unit to each participating nonmaximal corner. There are 190
two-face corners per vertex of the 20-dimensional facet. Consequently

`2T+Q<=H+190L`.

For the compact simple facet, symmetry and nonnegativity of the h-vector
give `f2>45f0`. Since the total number of polygon corners is `190f0`,

`35(H+L)<5f2-190(H+L)<=2T+Q<=H+190L`.

This is precisely `34H<155L`. The independently checked baseline has cost
38 only at the two specified marked types and integer cost at most 37
elsewhere, so its average is less than `37+155/189=7148/189`.
No assumption that the average number of polygon sides is below four is
used. The existence and minimality of the selected facet remain the
published narrow-facet dependency already recorded in the
[central-E8 review](fresh_sqrt5_narrow_e8_ports_20261016_review.md).

Generation and fresh saved replay passed. This is a conditional restriction
on the selected facet of a hypothetical compact arithmetic 21-polytope;
it is neither a whole-polytope frequency claim nor a finite-volume
realization of the saved local faces.

```sh
python3 scripts/fresh_sqrt5_marked_short_faces_20261016_review.py
python3 scripts/fresh_sqrt5_marked_short_faces_20261016_review.py --verify
```
