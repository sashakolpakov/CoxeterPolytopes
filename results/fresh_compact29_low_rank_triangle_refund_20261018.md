# A valid triangle refund at the compact rank29 target

18 October 2026. **New conditional face-charge rule; no dimension improvement.**
For the unconditional Vinberg weighting, a marked pair in a finite component
of rank2,3 or4 has weight2, while every E7 pair has weight1. These are the
weights in the [frozen compact29 reduction](fresh_compact_dimension_bound_20261015.md),
not the distance-at-most-two weights used in the arithmetic lane.

At every actual vertex of full type T=E7+E7+E7+E8, delete the unit on an E7
angle if its two-face is a triangle having a corner of original weight2.
Make no other change. This is a valid global face rule. The weight2 corner
cannot be a modified T angle, since T has no component of rank at most4.
It remains and by itself pays the triangle's required two. This also handles
several eligible deletions on one triangle. No quadrilateral changes, and
all corner weights remain nonnegative. Thus all old face inequalities survive.

Let R(v) count these refunded E7 triangle angles at a T vertex. Its new
total is exactly87-R(v). The rule gives the desired upper84 at a particular
target vertex precisely when R(v)>=3. A proof of this inequality at every
T vertex is not provided below. The other2,125 extremal forest types also
remain separate obligations for a universal compact29 exclusion.

## A designated triangle forces six refunds

Use the E7 numbering0-1-2-3-4-5,2-6 from the
[reviewed triangle-star theorem](fresh_compact29_three_triangle_star_constraint_20261018.md).
Its complete actual-face classification says that any designated triangle
at pairs(2,4),(2,5),(3,5) has the same third wall f: its only old product
is f.r5=-phi, where phi=(1+sqrt5)/2; every root has square2. It is orthogonal
to the other two E7 factors and the E8.

In fact replacing any old node j in{2,3,4,5} by f gives a finite full corner.
The changing rank7 part of those four diagrams is respectively

|j|Changing corner forest|Component containing f|
|---:|---|---|
|2|A2+H4+A1|H4 on{3,4,5,f}|
|3|A4+H3|H3 on{4,5,f}|
|4|D5+I2(5)|I2(5) on{5,f}|
|5|E6+A1|isolated f|

The inverse-Z-Gram projection lemma gives actual vertices. For each pair
j<k in that four-set, f and the two old side walls form the actual triangle
through the old vertex. At its corner replacing j, both k and f lie in the
displayed rank-at-most4 component, and that angle has weight2. At the corner
replacing k, j and f are disconnected and its weight is zero. The old E7
angle has weight1. Every one of these six triangles therefore has original
payment1+2+0=3 and qualifies for the refund.

For completeness, the five finite vertices are the entire actual simplex
four-face supported by old roots0,1,6 and the untouched components. The
positive-definite replacement sets realize its simplex edges, which exhaust
degree four at each vertex; convexity/connectedness excludes extra vertices.
This also proves persistence of the six triangles against every other wall.

Consequently R(v)>=6 whenever any designated triangle occurs. The reviewed
cross-factor incompatibility allows at most one such E7 factor. No quad
consumer needs classification for these six deletions: the changed angles
are on actual triangles, and their retained weight2 corners pay them.

## An inactive one-paid-corner quadrilateral, and its refund

The absence of designated triangles does not justify deleting their quad
units. Here is an explicit local counterconfiguration to that inference.
In one E7 choose the designated pair(3,5), with white support
W={0,1,2,4,6}=A4+A1. Put

    Y=4cos^2(pi/7),       X=(3Y-4)/(5Y-15).

Add norm-two walls c,d whose only old contacts are c.r3=-sqrt(X) and
d.r5=-sqrt(Y), and set c.d=0. The latter product is forced by the common
finite-vertex realization: the inverse E7 entries are

    J33=15/2,   J55=3/2,   J35=5/2,

and the zero-product condition is exactly

    (5X-3)(Y-3)=5.

The number Y is the root in(3,4) of Y^3-5Y^2+6Y-1=0. Its enclosure
81/25<Y<13/4 follows by substitution and the positive derivative on(3,4).
Hence X>23/5>4: c.r3 is strictly dotted and d.r5 is label7.

All four adjacent white-plus-side corners are positive definite. The three
new marked pairs are disconnected, in A4+A2+A1, E6+A1, and A4+3A1,
respectively. Their original payments are zero, while the old E7 corner
pays one. The white Schur block on the old sides has diagonals3/10,3/2
and product-1/2. Its full four-side extension has Lorentz signature(2,1,0)
by the displayed corank equation, with both opposite pairs indefinite.
It defines a compact Lambert quadrilateral. The acute projection lemma
makes its four finite corners actual in the local chamber.

No principal affine diagram occurs. A subset containing c and root3 has
a dotted negative minor; otherwise c is isolated. A connected subset
containing d and root5 contains a label7, which no affine Coxeter diagram
has; otherwise d is isolated. The remaining old E7 subdiagrams are finite.

This quad forces a different actual triangle at the old pair(4,5), with
third wall d. Its other corner diagrams are E6+A1 and D5+I2(7). Their
marked payments are0 and2. Thus its three payments are1,0,2, and refunding
the old E7 unit is safe. In particular E7(4,5) itself has weight1, not2.
The finite corners realize a closed three-edge cycle, proving the triangle
is actual; the same projection argument prevents further walls cutting it.

Three such local quads, one per E7 factor of T, coexist in one acute local
chamber. The new height squares are(15/2)X-2 and(3/2)Y-2, both greater than2.
Across distinct old factors the projections are orthogonal, so every mutual
new-wall product is negative and strictly dotted. A small negative old-dual
perturbation of the common vertex is a strict interior point. All finite
corners persist by projection, and cross-factor dotted pairs cannot belong
to an affine principal diagram. No designated triangle wall is among these
walls. The three forced(4,5) triangles give R(v)>=3, so the refund rule
saves the required three units at this local configuration.

This is not a compact or finite-volume chamber construction. It proves
that an inactive quad can consume its entire old unit, and that its forced
triangle can instead fund a legitimate refund. It does not prove that an
arbitrary inactive T vertex has this Lambert carrier in each factor, or
that every alternative quad consumer forces enough distinct refunds.

The remaining mathematical target is R(v)>=3, or another compatible saving,
at every T vertex with no designated triangle. That requires exhaustive
actual-quad or incidence/multiplicity control; the at-most-one active-factor
theorem alone is insufficient. No new catalogue computation, optimizer,
Lean process, or verification run was performed for this note. All frozen
inputs are unchanged. The new charge and local algebra have not yet had an
independent review.
