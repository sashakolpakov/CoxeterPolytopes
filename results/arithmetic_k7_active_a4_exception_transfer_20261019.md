# Moving the cubic-seven exception charge from D5 to A4+A1

19 October 2026. **New written charge theorem candidate, awaiting full
independent review.** Derived and checked by batch_guard while comparing
the remaining arithmetic fields. No scientific job, enumeration, Lean,
TeX or build ran. Earlier source-bound reports, including the D5
one-active theorem, are unchanged.

## 1. Inherited geometric input and exact core

Let P be a compact arithmetic Coxeter polytope over the actual ground
field K7=Q(cos(2pi/7)). Use the independently reviewed complete
quadrilateral classification in

    fresh_compact_global_bound_20261015.md,
    SHA256 0f7d1a4516a011c52ac6106721f820802d47b2cbd7a36840b20506be961dc59a,

and its separate full coverage/arithmetic review. Put k=cos(pi/7).
The old proof gives k>9/10, hence

    4/5<k^2<1.                                    (1)

Every triangular two-face keeps weight2 at an entire I2(7/14)
component. Every quadrilateral either has a marked pair at graph
distance at most2, or has precisely this sealed core, with four
cyclic sides a,b,c,d and three white roots u,v,w:

    a.c=b.d=-k,
    c.u=b.v=d.w=u.v=u.w=-1/2.                    (2)

All other core products are zero. Its corners(a,b) and(d,a)
are A4+A1; its corners(b,c) and(c,d) are D5. At(a,b), the
entire A4 is b-v-u-w and the entire A1 is a. At(d,a), the
entire A4 is d-w-u-v and the entire A1 is again a.

ENTIRE is an essential proved qualification. The inherited proof
exhausts all black-to-white neighbors, eliminates arbitrary other
attached white branches and uses singularity to force all three
reduced white diagonals p=q=r=1. Thus there is no attached branch.
Every remaining actual white component is orthogonal to ALL seven
core roots. At either A4 corner this makes both the A4 path and
the isolated A1 root whole connected components of the actual
finite vertex. It also makes each new core wall orthogonal to
the complete old spectator. We do not infer that stronger statement
just from the finite corner diagram alone.

## 2. Charge rule and simultaneous face payment

Retain the two original nonnegative terms:

* one unit for marked graph distance at most2;
* one extra unit for an entire marked I2(7) or I2(14) component.

Remove all D5 half-bonuses. Instead add1/2 at a marked pair
consisting of an endpoint of an entire A4 and the root of an
entire A1, IF its actual two-face is a quadrilateral with no
marked distance-at-most-two corner. Give no new bonus otherwise.
Call such a coupon ACTIVE. Its selector uses the actual two-face,
whose existence and incident marked pair are fixed; it does not
choose a favorable state to minimize a vertex budget.

Every face is paid simultaneously. A triangle keeps its whole
high-I2 weight2. A quad with a short marked corner keeps its
original baseline unit. A quad without a short corner is exactly
(2), so its TWO A4+A1 coupons are active and pay1 together.
These are the two corners incident to side a, not the D5 corners.
All weights remain nonnegative, so every k-gon with k>=5 has
at least its required payment5-k.

The remaining issue is the total number of active coupons at one
vertex. Counting every abstract endpoint/isolate pair would be
wrong; the next section uses actual extra wall geometry.

## 3. At most one active coupon for each entire A4 component

Fix an actual finite vertex p with p^2=-1. Fix one of its entire
A4 components, write its path as

    b-v-u-w,

and number these positions1,2,3,4. An active coupon at endpoint b
and an old entire isolate h supplies, by(2), an extra ACTUAL unit
wall r with contacts

    r.h=-k, r.u=-1/2,
    r perpendicular to every other old vertex root.             (3)

The corresponding coupon at endpoint w supplies a wall touching
its isolate by-k and v by-1/2. Complete spectator orthogonality
in(3) follows from the sealed-core theorem in Section1, including
all other isolates and all other connected vertex components.

For the unit A4 Gram G,

    (G^-1)_22=(G^-1)_33=12/5,
    (G^-1)_23=8/5.                                (4)

These follow from (G^-1)_ij=2*min(i,j)*(5-max(i,j))/5. Since
the isolate is orthogonal to the whole A4, the projection R of
any r in(3) against the FULL old finite vertex has square

    q=R^2=k^2+3/5>1.                              (5)

At the actual old vertex r is a distinct facet, so r.p<0.
Its decomposition is r=R+lambda*p, with

    lambda=sqrt(q-1)=sqrt(k^2-2/5)>0.

Thus two such actual walls have the SAME positive timelike
coefficient. This fixes the sign of their timelike cross product;
one cannot choose it independently to satisfy acuteness.

Suppose two distinct coupons for this fixed A4 are active, and
take their actual walls r,r'. There are exactly three cases.

* SAME endpoint, DIFFERENT isolates: the positive projection
  cross product is3/5. Therefore

      r.r'=3/5-(k^2-2/5)=1-k^2>0.

* OPPOSITE endpoints, SAME isolate: the positive cross product
  is k^2+2/5. Therefore

      r.r'=k^2+2/5-(k^2-2/5)=4/5>0.

* OPPOSITE endpoints, DIFFERENT isolates: the positive cross
  product is2/5. Therefore

      r.r'=2/5-(k^2-2/5)=4/5-k^2 in(-1/5,0).      (6)

The first two contradict acuteness of distinct outward walls.
The last contradicts the actual Coxeter gap: a negative product
strictly between-1 and0 gives an actual ridge by the two-normal
positive-definite projection argument. The nonright Coxeter angle
then has product-cos(pi/m)<=-1/2,m>=3. No actual product can
lie in(-1/2,0). Equation(1) puts(6) strictly in that forbidden
interval.

For completeness, the projection argument starts from a strict
interior timelike point and projects orthogonally to the two wall
normals. Their positive-definite Z-Gram has nonnegative inverse,
so the required coefficients add nonnegative multiples of the two
outward normals. Every other wall product remains strictly
negative by acuteness, and the projected vector is timelike.
The intersection is consequently an actual ridge. This justifies
using the gap for the additional walls in(6), rather than treating
arbitrary formal acute Gram entries as Coxeter labels.

The contact rows of two distinct coupons differ, so the two walls
cannot coincide. For one endpoint and one isolate there is only
one marked angle, hence only one actual two-face, so there is no
fourth case involving two different coupons at the same pair.

We have proved the capacity statement

    active coupons at each entire A4 component <=1.             (7)

This includes arbitrary numbers of entire A1 spectators and both
A4 endpoints. It is not a per-isolate bound or a selected-support
inequality.

## 4. Vertex budgets in arbitrary rank

An A4 component has original distance-two cost5. By(7), all new
coupons assigned to it contribute at most1/2, even with many
isolates. Thus its complete charge is at most11/2. D5 has no
new bonus and has cost8. Every other component retains its
original distance-two/high-I2 cost.

Consequently EVERY connected component of rank r has charge at
most2r-2. In particular the A4 upper11/2 is below6, while A1
contributes0. The coupon normalization is accounted for once per
A4 component, not once per abstract endpoint/isolate pair.

For a rank20 vertex with at least three components this gives
total charge at most40-2*3=34. If it has one or two components,
it has NO active coupon: a coupon requires BOTH a whole A4 and
a whole A1, whose ranks sum to5, so a rank20 vertex containing
them must have at least one further component. Thus every one-
or two-component rank20 cost is exactly its uncorrected baseline.

The complete forests with rank20 charge possibly above35 are
therefore precisely:

|Entire forest|Upper|
|---|---:|
|A20,B20|37|
|D20|38|
|D_r+D_(20-r),r=4,...,10|36|
|D12+E8,D13+E7,D14+E6|36|
|D19+A1|36|
|D18+I2(7),D18+I2(14)|36|

To see completeness, a path of rank at least2 has deficit3
from twice its rank; A1,D/E and the entire high-I2 components
have deficit2. In a two-component forest, both deficits must
be2 to exceed35. The displayed D/E pairs, D19+A1 and the
two high-I2 spectators exhaust that rank equation. The only
connected rank20 finite types are A20,B20,D20. No forest
enumeration or field-specific height cutoff is used here.

In particular D5+A15/B15 now have charge at most35, D5+D15
has charge36, 2D5+D10 has charge34 and4D5 has charge32.
The 35.5/36.5 obstructions from retaining an active D5 half
have disappeared. This stronger transfer supersedes that budget
while preserving the earlier one-active-D5 theorem as evidence.

## 5. What is still required for a dimension19 result over K7

The new theorem is the valid ACTIVATED A4+A1 transfer and its
complete per-component capacity, not a new dimension ceiling.
The corrected charge still has A20,B20,D20 above35 and the
displayed two-component cost36 families. Later eight-field
refunds have not been silently extended to K7.

Their single-half-source local branches are promising: F/T
native cos(pi/7) already fail the tiny full-star conjugate window,
the B finite-square table includes labels7 and14, and the
universal dotted-long unit proofs have no field-list restriction.
But the complete old-source audit, added native e-label7/14
coverage, all simultaneous consumers and common-white interfaces
remain obligations before importing any such rule. The active
coupon's complete white graph is A3 plus its orthogonal spectator;
this must be tracked when new deletions are composed.

The inherited complete K7 core theorem, actual finite-vertex
simplicity, arithmetic finite/Lanner restrictions, acute-wall face
correspondence and the actual Coxeter gap are the explicit inputs.
The new capacity uses only(1)--(5) and full outward orientation.
It does not assume a universal facet decomposition, a completed
rank20 example, a new enumeration fact, or published priority.
