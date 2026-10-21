# A conditional reduction of the missing face coverage

This reduces the coverage problem for a proposed compact arithmetic21 weight
rule. It does **not** certify the proposed rule or exclude dimension21.
Assume the distance-two discount is exactly zero, the adjacent discount is
0 <= delta <= 2/3, and every added compensation or cross payment is
nonnegative. The adjacent discount applies only to a3-edge in an entire
spherical component of rank at least6. These hypotheses must be checked
on the eventual exact coefficients; a numerical vector is insufficient.

The baseline is the graph-distance1/2 indicator with the two entire-H4
corrections. Its complete triangle theorem and extension monotonicity are
inherited from the frozen [triangle package](fresh_sqrt5_triangle_weights_20261015.md).
The [source](../scripts/fresh_sqrt5_discounted_face_domain_reduction_20261016.py)
and [certificate](fresh_sqrt5_discounted_face_domain_reduction_20261016.json)
check all68 spherical component types of rank at most21. Every discounted
pair has baseline exactly1, including eligible B,D,E components. All
branching vertices of these finite trees have degree3 and only3-labelled
incident bonds.

## The only payments that can need new checks

Write B for the integral baseline total and d for the number of discounted
corners of a face. Its payment is at least B-d*delta >= B-2d/3. All corner
weights remain nonnegative.

| Face | Remaining potentially insufficient cases (B,d) |
|---|---|
| Triangle, required payment2 | (2,1), (2,2), (3,2) |
| Quadrilateral, required payment1 | (0,0), (1,1), (2,2) |

The baseline-zero quadrilateral case already belongs to the inherited
coverage domain, but its complete coverage still must be rebound when
the free payments change. The genuinely new paid-quadrilateral cases are
B=1 with its paid corner discounted, and B=2 with both paid corners
discounted. Triangles with B>=4 and quadrilaterals with B>=3 are safe.
Faces with at least five sides need only nonnegative weights in this
particular combinatorial argument.

There cannot be three discounted corners of a compact triangle. Normalize
the three black normals to square2. Three3-bonds give the affine-A2 Gram,
and z=n0+n1+n2 is null. For a strict future interior H, the product H.z
is negative, so z is nonzero and future null. Acuteness gives z.n<=0
for every other wall. Its ray is therefore an ideal point of the chamber,
contradicting compactness. This argument does not need arithmeticity.

## Complete local contact choices for triangle white extensions

The same positive-null argument excludes any affine-D4 principal diagram:
its null combination has coefficient2 at the center and1 at each of the
four leaves. This gives the following degree bound for a white node in a
compact triangle star.

Every diagram obtained by removing one black node is spherical. A white
node of degree at least5 would still have degree at least4 in a corner,
which is impossible. A degree4 white node must be joined to all three
black nodes and one white node; otherwise some corner retains degree4.
Each corner then forces the three retained bonds to have label3. Any edge
between two of the four neighbors would survive in a corner and create
a cycle, contradicting spherical finiteness. The five nodes consequently
form an induced affine-D4 diagram, again impossible. Thus every white
node has degree at most3 in the full local star.

For connected growth, a new white node therefore has these contact choices:

* Degree1: label3,4,5,6 or10.
* Degree2: ordered labels (3,3), (3,4), (4,3), (3,5) or (5,3).
* Degree3 with at least one white neighbor: all three labels are3.
* Degree3 with all three black neighbors: all labels3, or exactly one
  label4 or5 and the other two labels3.

For degree2, a corner contains both neighbors. The three-node principal
Gram is positive definite only when the two squared bond products sum
to less than4; the source checks exactly the five listed ordered choices.
For degree3 with a white neighbor, a corner retains the whole star, so
the spherical branch classification forces all3-bonds. With three black
neighbors, every pair of bonds must pass the same three-node test; the
source checks exactly the seven ordered triples above.

These are exhaustive local branching choices, not sufficient acceptance
conditions. Every retained graph still needs all three spherical corners,
the physical and conjugate Gram tests, exact isomorphism control and the
rank bound. Orthogonal components and the marked facet root must be
accounted for when evaluating the payments.

The baseline is nondecreasing under allowed white extensions, and the
number of3-edges between the original black nodes never changes. Hence
one may prune B>=4, B=3 with fewer than two black3-edges, and B=2 with
no black3-edge. One may **not** prune merely because a current corner
component has rank below6: further white nodes can make it eligible.

Generation and a fresh saved replay of the finite algebra and catalogue
checks pass. The domain arguments above are mathematical proof premises
in addition to those checks. No complete graph-growth run or face-payment
enumeration is asserted by this package.

```sh
python3 scripts/fresh_sqrt5_discounted_face_domain_reduction_20261016.py --verify
```
