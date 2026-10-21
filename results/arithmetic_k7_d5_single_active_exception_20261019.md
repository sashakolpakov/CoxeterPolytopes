# A single active D5 exception over the cubic-seven field

19 October 2026. **New written charge theorem candidate, awaiting
independent review.** No enumeration, scientific job, Lean, TeX or build
was run. All older source-bound evidence is unchanged.

## 1. The inherited complete quadrilateral theorem

Let P be a compact arithmetic hyperbolic Coxeter polytope with actual
ground field K7=Q(cos(2pi/7)). We use the already independently reviewed
complete quadrilateral theorem in

* `fresh_compact_global_bound_20261015.md`, SHA256
  `0f7d1a4516a011c52ac6106721f820802d47b2cbd7a36840b20506be961dc59a`;
* its separate full arithmetic/coverage review
  `fresh_compact_global_bound_20261015_review.md`.

Its native finite labels are2,3,4,6,7,14. Every Lanner triangle
has a7 or14 bond; that bond is an entire rank-two component at
its finite corner. A triangle therefore keeps weight2 under the
distance-at-most-two baseline plus the entire-I2(7/14) extra unit.

The complete all-white quadrilateral theorem says that a quad
with no marked distance-at-most-two corner is the following sealed
exception. Put k=cos(pi/7). There are four side roots a,c,b,d,
with opposite products a.c=b.d=-k, and three white roots u,v,w.
All other nonzero core products are simple:

    c.u=b.v=d.w=u.v=u.w=-1/2.                     (1)

There are no other core contacts. The two exceptional finite
corners (c,b) and (c,d) are D5; their marked roots are the long
leaf and a short leaf, at graph distance3. The other two corners
are A4+A1 and have baseline zero.

Two sealing statements from the full theorem are essential here.
All black-to-white neighbors have been exhausted in constructing
this core. Eliminating arbitrary remaining attached white branches
gives diagonal parameters p,q,r<=1; the exact singularity argument
forces p=q=r=1, so there is NO nonorthogonal remaining branch.
Consequently every remaining white component is orthogonal to the
ENTIRE seven-root core, not merely to its five-root D5 corner.
This supplies orthogonality of the new walls to the whole vertex
spectator. We do not infer it just from D5 being a whole component.

## 2. New activation rule

Start from the K7 charge:

* one unit at marked distance at most2;
* an extra unit for an entire marked I2(7) or I2(14) component;
* one half at each of the TWO distance-three long-leaf/short-leaf
  pairs of an entire D5 component.

Change only the last line. A D5 half-angle is ACTIVE exactly when
its actual two-face is a quadrilateral with no marked short corner.
Retain its1/2 if active and give it zero otherwise. This is a
physical face selector, not a replacement of the proved exception
classification by a numerical constraint.

All polygon payments are preserved. Triangles keep their unchanged
entire high-I2 weight2. A quad with a short corner keeps its baseline
unit, independently of the D5 additions. A quad without one is
exactly(1), so both its D5 half-angles are active and pay1 together.
For k-gons with k>=5, nonnegativity suffices. Thus activation is
simultaneously safe on EVERY actual face.

## 3. Capacity: at most one active half-angle per D5 component

Fix an actual finite vertex p, normalized by p^2=-1, and one of
its entire D5 components. Name its fork u, its short leaves b,e,
its long near node v and long leaf l. Its two bonus pairs are
(b,l) and(e,l). Suppose both are active.

At the exception for(b,l), the core description(1) supplies an
additional ACTUAL outward unit wall r whose product is-k with b
and zero with EVERY other old vertex root. At the exception for
(e,l), it supplies an actual wall r' whose product is-k with e
and zero with every other old root. The respective complete
spectator orthogonalities follow from the sealing argument in
Section1. In particular these are contact rows for the entire
rank-n finite vertex, not merely a selected D5 subdiagram.

For the unit D5 Gram inverse the two short-leaf entries are

    (G^-1)_bb=(G^-1)_ee=5/2,
    (G^-1)_be=3/2.                               (2)

These follow directly from the standard D5 root coordinates, or
from solving the simple star with arms1,1,2. Spectators are
orthogonal, so they do not change(2). Let R,R' be the projections
of r,r' to the whole old positive vertex span p-perp. Then

    R^2=R'^2=q=(5/2)k^2>1,
    R.R'=(3/2)k^2.                               (3)

Both extra facets are strict at p. Hence r.p<0 and r'.p<0.
Writing

    r=R+lambda*p, r'=R'+lambda'*p,

fixes lambda,lambda'>0. Unit normalization and(3) give
lambda=lambda'=sqrt(q-1), with NO sign ambiguity. Therefore

    r.r'=(3/2)k^2-(q-1)=1-k^2>0.                (4)

The contact rows are different, so r,r' are distinct walls.
Equation(4) contradicts the nonpositive products of actual outward
walls of an acute Coxeter chamber. Thus at most ONE of the two
half-angles is active at this component.

The proof uses the full actual finite vertex to fix the timelike
complement and the signs. A formal D5 Gram pattern without whole
spectator sealing and actual outward orientations would not prove
this capacity statement.

## 4. Improved component budget and complete rank20 frontier

The distance-two count of D5 is8. Its old two half-bonuses raised
its cost to9; activation and capacity give the new upper

    D5 cost <=17/2.                               (5)

This theorem applies separately to every whole D5 component at a
vertex, in every dimension where it occurs. It does not require
choosing a favorable active state or deleting a half-angle on an
exceptional face. The face selector is determined first, and(4)
then bounds its actual multiplicity.

For clarity, the COMPLETE rank20 forest frontier for THIS K7
charge, before importing ANY later eight-field refunds, is:

|Entire forest|Vertex upper|
|---|---:|
|A20,B20|37|
|D20|38|
|D5+D15|73/2|
|D5+A15,D5+B15|71/2|
|D_r+D_(20-r), r=4,6,7,8,9,10|36|
|D12+E8,D13+E7,D14+E6|36|
|D19+A1|36|
|D18+I2(7),D18+I2(14)|36|

Every other rank20 finite forest costs at most35. Here is the
complete reduction, requiring no forest enumeration. A component's
deficit from twice its rank is at least2, except that D5 has
deficit at least3/2. Paths have deficit3, A1 deficit2, and the
entire high-I2 components have deficit2. There are no H components
because label5 is not native; F4 is a path with deficit3.

With at least four components the deficit is at least6, so the
cost is at most34. With three components, zero or one D5 gives
deficit at least11/2. Two D5 factors leave rank10, whose smallest
possible deficit is2 (at D10), giving cost at most35. Three D5
factors alone have rank15 rather than20. With two components,
a D5 leaves rank15, allowing exactly A15,B15,D15. Without D5,
cost above35 requires two deficit2 components, and the table
lists every rank20 possibility. With one component the only
finite types are A20,B20,D20. This proves completeness.

In particular the old extra heavy forests4D5 and2D5+D10 are now
controlled, with upper34 and35 respectively. D5+A15/B15 and
D5+D15 remain at35.5 and36.5. These are explicit next obligations,
not omitted cases or excluded actual vertices.

## 5. Frontier and limits

K7 is a concrete next field because the ONLY exceptional component
addition in its complete quad rule was the D5 half-pair correction,
and the new activation removes its larger multiplicities. In
contrast the existing Q(sqrt2) rule also has F4,E6,E7,E8 additions,
and its separately source-bound silver opposed-fork star survives
the present local Lorentz and arithmetic tests. That star prevents
an unbalanced three-unit D20 deletion from being justified by those
tests alone; it is not a completed compact20-polytope.

The new theorem here is the ACTIVE D5 rule and capacity(4), with
the precise component/forest consequences(5). We have NOT extended
any later eight-field A/B/D refund to K7. The two spectator families
in the table still need real consumer proofs, and label7/14 contacts
must be included whenever such an extension is attempted. No new
dimension ceiling, example, global decomposition or novelty claim
is made. The existing K7 ceiling remains20 until the remaining
rank20 vertex costs are all brought to at most35 with valid
triangle and quadrilateral payments.
