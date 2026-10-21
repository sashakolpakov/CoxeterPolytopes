# Triangle payments with finite distance9 and ideal distance6

This local theorem gives every actual triangular two-face total corner payment
at least2 under the following rule. A finite corner pays1 if the marked sides
have graph distance at most9 in its full vertex diagram. An ideal corner pays1
if its varying affine factor has at most9 nodes, or the marks have distance at
most6, or their distance is at most9 and an endpoint is within6 of a branch
vertex or an endpoint of a nonsimple bond. A corner pays at most1 under this
rule. There are no cross-component coupons in this triangle theorem.

The [exact checker](../scripts/fresh_ideal6_triangle_payments_20261017.py)
and [data](fresh_ideal6_triangle_payments_20261017.json) enumerate all920
marked triples in the inherited minimal indefinite diagrams of orders8,9,10.
The complete classification is an explicit input from the
[high-rank census](fresh_high_rank_quasilanner_census_20261017.md).
The actual supporting-white and cusp-factor geometry is inherited from the
[factor framework](fresh_cusp_factor_payment_framework_20261017.md).
This package alone gives no new ambient dimension bound.

## Reduction to the finite census

Let W be the elliptic supporting diagram of an actual triangular two-face,
and let its sides be a,b,c. Its full star is indefinite. Every W together with
two sides is a finite or ideal corner and hence has positive-semidefinite
Gram matrix. Consequently any minimal indefinite subdiagram contains all
three sides. Such a diagram has order at most10 by the stated classification.
For order at most7, a Steiner tree on the three marked nodes has at most6
edges. At least two of its marked-to-marked paths avoid the third mark.
Each gives a paid corner, including at an ideal vertex.

For orders8,9,10, remove the third mark for each corner and find its marked
path and branch/bond features. There are224,336,360 triples respectively.
Only nine triples fail the pure distance6 test. The feature rule pays at
least two corners in all920 triples except diagram25 with marks(1,2,9).
Every saved feature witness persists under additional supporting whites:
an existing branch or nonsimple bond cannot disappear, and an existing
path cannot become longer. At an ideal corner the connected marked core
lies in its varying affine factor. At a finite corner every saved path has
length at most9. Thus these are witnesses in the full corner, not assumptions
that a small finite subdiagram remains a whole finite component.

## The one exceptional triple

Diagram25 is simply laced, with edges

    0-1, 0-2, 0-4, 2-3, 4-5, 5-6, 6-7, 7-8, 8-9.

Deleting9 leaves affine E8 on nodes0,...,8. This is present at corner(1,2).
The exact checker verifies its positive kernel, rank8 and a positive-definite
eight-node principal matrix. Every extra supporting white is orthogonal to
this entire affine diagram: positive semidefiniteness at corner(1,2) forces
its pairing with the positive kernel to vanish, and all individual pairings
are nonpositive. Corner(1,2) already pays through its path1-0-2.

At corner(1,9), deleting2 leaves A8+A1, with marked path

    1-0-4-5-6-7-8-9

and isolated node3. The checker verifies this labelled path and positive
definiteness exactly. The preceding orthogonality forces node1 to remain
a simple-bond leaf and its neighbor0 to retain degree2 in the full corner.
No classical irreducible affine diagram has a simple-bond leaf whose
neighbor has degree2: affine A has no such leaf, affine C has nonsimple
end bonds, and the simple leaves of affine B and D meet branch vertices.
The remaining irreducible affine diagrams have at most9 nodes. This uses
the standard finite/affine classification already assumed by the census.

Therefore, if corner(1,9) is finite, its length7 path pays1. If it is ideal,
its varying affine factor must be exceptional and has at most9 nodes,
so the whole-factor rule pays1. Either way the triangle receives its second
payment, with arbitrary extra supporting whites included.

## Verification scope

The checker uses standard-library rational arithmetic and no project arithmetic
imports. It verifies all920 triples, every path/feature witness, and the literal
exceptional Gram calculations. It does not redo the inherited classification,
formalize the geometric argument, or prove a quadrilateral payment rule.
Generation and a separate saved replay both terminated exit0/PASS. A separate
agent read the complete source and written argument and found no identified
gap in the kernel, fixed-leaf, feature-persistence or lower-order steps. That
review did not rerun the920 arithmetic cases or the inherited classification.
No external refereeing or proof-assistant validation is claimed. This source,
data and report are final/frozen.

    python3 scripts/fresh_ideal6_triangle_payments_20261017.py
    python3 scripts/fresh_ideal6_triangle_payments_20261017.py --verify
