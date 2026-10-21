# Deep finite distance-eight payments can be deleted at feature depth8

Start with the complete simultaneous payment rule of the
[dimension61 proof](fresh_rank15_rank16_dimension61_20261018.md).
Set to zero every remaining payment at a distance8 marked pair in a
whole finite A/B/D component when both marks have distance at least8
from every endpoint, branch and endpoint of a nonsimple bond. Retain
every other finite weight, every ideal weight and every component-pair
coupon. Then **every triangle receives at least2 and every quadrilateral
at least1**.

The change lowers the feature threshold of the completed
[depth11 zero theorem](fresh_depth11_finite8_zero_transfer_20261018.md).
Its main new step is an estimate for the Schur diagonal after subtracting
the new root's squared norm, rather than separate bounds on the projection
and correlation. A second argument excludes the minimal cores that could
escape the shorter retained path. The two small cycle boundary cases are
handled separately; their profile domains differ from the long tail.

The [source](../scripts/fresh_depth8_finite8_zero_transfer_20261018.py)
and [certificate](fresh_depth8_finite8_zero_transfer_20261018.json) check
the exact matrix and scalar inequalities. This is a local payment theorem.
It does not itself give a complete finite-forest maximum, improve the
ideal envelopes or prove a dimension bound below61.

## The actual-face inputs and the shorter path

Use an actual quadrilateral star W+a+b+c+d with opposite black pairs
{a,c} and {b,d}, normalized to squared root norm2. W and W plus any one
black are positive definite. W plus either opposite pair is indefinite.
These are the established actual-face inputs of the
[universal corridor theorem](fresh_universal_deep_finite8_corridor_20261018.md).
Every distinct-root product is nonpositive, and every nonzero native
contact magnitude is at least1. Arbitrarily large finite labels, parallel
products and dotted products are included.

Suppose W+a+b has the stated whole finite component and marked pair.
Feature depth8 supplies two simple outward arms of eight whites, in
addition to the eight-edge a-to-b path. A feature exactly at distance8
is allowed: all bonds before it are still simple. Keep six outward
whites on either side. The resulting21 roots form an induced A21 with

    a at6, b at14, indexed from0.

Its other19 roots are supporting whites. The two additional whites on
each outward arm remain available for actual-edge estimates.

Assume there is no path of length at most7 between the opposite black
pairs whose internal roots are white. Then c.d=c.b=d.a=0. The own-opposite
products c.a and d.b remain unrestricted. For a norm-two A_n path write

    J_n[i,j]=(min(i,j)+1)*(n-max(i,j))/(n+1).

In either full eight-white outward arm, the inverse diagonals at distances
j=3,4,5,6 from the old mark are respectively

    2, 20/9, 20/9, 2.

A contact of magnitude at least1 already has projection square at least2.
The actual edge containing that white arm and c is strictly positive
definite, so equality is also excluded. Positivity of the path inverse
means that extra contacts cannot reduce this obstruction. This excludes
those four outer positions on either arm, and likewise for d.

A c-contact at an inner position at least two steps from a supplies a
path to b of length at most7. A c-contact within two steps beyond b does
the same. Hence the only possible retained contact positions are

    c: 4,5,6,7;       d: 13,14,15,16.

No bound on their contact magnitudes is imposed. Contacts outside A21
remain allowed throughout.

## Nonzero retained contact, including a disconnected deletion

Choose a minimal indefinite principal K inside W+a+c. It contains both
a and c, since deleting either leaves a principal part of a strictly
positive-definite actual edge. K is connected and has order at most10.
The published Chein/Vinberg order bound, including the bridge for dotted
order2 diagrams and disconnected positive-semidefinite deletions, is
explicitly recorded and reviewed in the bound
[universal corridor report](fresh_universal_deep_finite8_corridor_20261018.md).
This proof uses that all-label order ceiling, not a bounded search.

Let B be the connected component through a in K minus c. It has at most
nine roots and lies in the same whole finite tree as a. Feature depth8
forces it to be a simple path: reaching a possible feature at distance8
already uses all nine available nodes, leaving no node for branching or
a bond beyond the feature. It cannot pass through b, which is absent
from K. Connectedness of K forces c to contact B.

Suppose c has no retained A21 contact. All of B toward b lies inside the
retained path. On the outward side, the only possible contacts missed by
that path are distances7 and8. The full eight-white edge arm gives

    J8[6,6]=14/9, J8[7,7]=8/9, J8[6,7]=7/9.

At distance7 a nonsimple contact has projection at least28/9>2, so only
a simple bond is possible. At distance8 a label at least5 has squared
native magnitude at least(3+sqrt5)/2>5/2, giving projection greater than
20/9>2. Thus only labels3 and4 are possible there. Contacts at both
positions have projection at least(14+8+14)/9=4 and are impossible.

There are now five exhaustive cases.

* If the contact is only at distance8, B already has nine roots. There
  is no further node in K beyond c. Appending c with bond3 or4 gives
  finite A10 or B10, contradicting indefiniteness.
* If the contact is only at distance7 and B has nine roots, the ninth
  B-node is either one inward step from a or outward distance8. The
  simple c-contact then gives finite A10 or D10, respectively.
* If the contact is only at distance7 and B has eight roots, it is
  exactly a and the seven outward whites. With no further white in K,
  appending c gives finite A9.
* In the same eight-root case, K can instead contain one further white
  z outside B. The definition of B makes z orthogonal to B, and
  connectedness forces z to contact c. The strictly positive-definite
  diagram K minus a is an A8 path (the seven outward whites and c)
  with z attached to its end. Its inverse endpoint8/9 excludes a
  c-z label at least5. A bond3 or4 makes the full K finite A10 or B10.

The second bullet has two alternatives, giving five cases in the
certificate. In particular **K minus c need not be connected**; the
extra z case is essential and is not silently discarded. The source
checks positive pivots for every finite outcome. For B10 it uses the
explicit rational congruence that scales its terminal root by sqrt2.

Every case contradicts minimal indefiniteness. Therefore c has a nonzero
retained contact. Reflection proves the same for d. The argument permits
arbitrary supporting spectators and does not assume that the two minimal
indefinite cores intersect.

## A sharper uniform Schur inequality

Let x,y be the nonnegative contact-magnitude vectors against A21 and put
J=J21, q=x^T J x, r=y^T J y and s=x^T J y. Each vector has a coordinate
at least1. On the left set u_i=i+1, on the right set u_i=21-i, and put

    alpha=17/352, U=sum_left x_i*u_i, V=sum_right y_i*u_i.

On either four-position block the matrix

    C[i,j]=J[i,j]-alpha*u_i*u_j

is entrywise nonnegative and every diagonal is at least2. Thus one
nonzero native coordinate gives x^T C x>=2, and similarly for y. The
cross block of the path inverse factors exactly. Consequently

    q-2 >= alpha*U^2 >0,
    r-2 >= alpha*V^2 >0,
    s = U*V/22.

Since c.d=0, eliminating A21 leaves the Schur matrix

    [2-q, -s; -s, 2-r].

Its diagonal entries are negative and its determinant satisfies

    (q-2)*(r-2)-s^2
      >= (alpha^2-1/22^2)*U^2*V^2
       = (3/11264)*U^2*V^2
      >= 1875/11264 >0,

where U,V>=5. It is negative definite, giving two negative directions
in an actual Lorentz Gram matrix. This is impossible. The full-star
short7 alternative therefore exists at feature depth8.

The source checks the A8/A21 inverses, all positive leading determinants,
32 residual block entries, eight diagonal lower bounds and16 cross
factorizations. The bound applies to all nonnegative native magnitudes;
it is not inferred from a finite contact sample.

## The cycle interface has two distinct boundary domains

The inherited [cycle interface](fresh_deep_finite8_cycle_interface_20261017.md)
and [finite8 replacement](fresh_finite8_cycle_replacement_20261017.md)
give the complete single-root profile classification. Endpoint4 has a
nonsimple bond at its new mark. Second-node3 has a fixed branch next to
its new mark. Neither can give a changed deep finite8 pair. Endpoint3
on an arc of length8 is the only affected single-contact profile.

Its opposite old endpoint stays a whole-component leaf, at distance
complementary_arc minus1 from the old mark. Thus the outputs for arcs
(6,8),(7,8),(8,8) remain nondeep and unchanged. Only P=8,Q>=9 remains.
If the endpoint3 root on the8-arc has no spectator neighbor, that new
root is itself a leaf and its pair is also unchanged. An unpaid changed
pair therefore forces an actual spectator neighbor.

At Q=9 a second neighbor requirement is inherited: endpoint3 on the
9-arc otherwise gives a whole-leaf9 payment. The complete residual lists
are as follows. Subscript0 denotes the8-arc, and subscript1 the Q-arc.

|Q|Residual profiles|Profiles requiring a spectator neighbor|
|---:|---|---|
|9|isolated, end3_0, end3_1, both_ends3|end3_0, end3_1|
|10|isolated, end3_0, end3_1, end4_1, both_ends3|end3_0|
|>=11|isolated, end3_0, end3_1, end4_1, next3_1, both_ends3|end3_0|

At Q=9, end4_1 is a sealed whole-leaf9 unit and next3_1 a nondeep
distance8 unit. At Q=10, next3_1 is a sealed whole-leaf9 unit. Their
whole finite components have ranks17 or18, so the separate large-D
leaf9 reductions starting at rank20 do not affect them. These profiles
must be discharged before taking the residual domain; retaining all six
tail profiles at Q=9 or10 would give false nonpositive lower bounds.
For Q>=11 all three long-arc distances Q,Q,Q-1 are at least10, so the
six-profile list is constant.

The spectator argument is the inherited unequal-height one. Adjacent
black contacts or a common spectator contact already pay through paths
of length1 or2. Otherwise the spectator supports x,y are disjoint. For
one required neighbor i of c, the norm-two spectator Gram G gives

    (x/h-y/k)^T G^(-1)(x/h-y/k) >= x_i^2/(2h^2) >=1/(2h^2).

Use the c-side if it requires a neighbor, otherwise the d-side. The side
is fixed for the entire reciprocal-height rectangle, even if both require
neighbors. Two coordinate bounds are never added.

With A=1/8, B=1/Q, u=1/h, v=1/k, L=u*l1+v*m1 and M=u*l2+v*m2, the
inherited normalized determinant is

    1-2*(A*L+B*M)-(A+B)*(c0*u^2+d0*v^2)-A*B*(L-M)^2,

plus the nonnegative spectator term. The chosen Cauchy bound reduces
that endpoint3 norm constant from1 to1/2. Corrected constants and moments
are nonnegative. The lower bound is concave on the reciprocal-height
rectangle and nonincreasing in B. Four rectangle vertices suffice for
each ordered profile pair; Q=11 bounds the entire constant-list tail.

|Region|Ordered pairs|Exact rectangle values|Strict minimum|
|---|---:|---:|---:|
|Q=9|16|64|1/18|
|Q=10|25|100|1/10|
|Q>=11|36|144|1/11|

All308 values are positive, contradicting the exact zero determinant of
an unpaid actual star. No equality case remains. The source independently
normalizes the eight inherited profile records in Q(sqrt2); the endpoint4
square roots cancel against their reciprocal heights. This scalar check
does not purport to reprove the inherited geometric classification.

The other cycle regions are unchanged. The classical-tree short7 route
uses finite paths at most7, persistent ideal features or sealed whole A7;
the small-cycle routes have fixed nearby leaves or sealed small components.
They never consume a newly deleted deep distance8 payment.

## Composition with the current dimension61 rule

The precise baseline is the completed dimension61 rule, including all
rank15/rank16 changes, the A17 restoration, D17 cost120, D13-zero coupon
increments and actual E6/D11 coupon cost4. No historical cheaper ideal
rule is silently substituted for it.

First handle the baseline's independent actual-corner priorities: C11
endpoint zero, B12 center/terminal zero, D13 center zero, cycle11 distance5
half, and the coupled order11 B/C/D halves. Their proofs are bound in the
[sparse composition](fresh_sparse_coupon_order11_13_dimension62_20261018.md),
[order11 composition](fresh_order11_order17_dimension63_20261018.md) and
[current ideal closure](fresh_order14_ambient62_ideal_margin_closure_20261018.md).
Their outputs are paths of length at most7, retained whole affine center
units, typed small whole-component coupons, whole finite components of
rank at most14, whole A9/A18 endpoint8 pairs, or whole D18 leaf9 pairs.
Endpoint8 pairs are nondeep, and the other outputs miss the new deletion.
Distinct retained halves still add to1. These priorities remain independent
of the generic short8 route.

In the remaining branch, consider the actual corners of the face. If
none has a newly deleted finite pair, all its weights agree pointwise
with the proved dimension61 rule and its coverage is unchanged. If a
newly deleted corner occurs, the new corridor gives a full-star path of
length at most7. At a finite corner it is an unchanged unit. At an ideal
corner the established elementary, small-tree and protected-feature units
or the complete cycle/tree replacements apply. The priority types above
are absent on the whole face, so no alternative is asked to spend one
of their removed units. The cycle and tree interfaces just proved supply
the same retained payments under the depth8 rule.

The later order17 D zeros and order14 B/C or order16 B deletions have
distance10 or11. They cannot remove a short7 output. The predecessor's
18 small classical tree targets have their short8 roles explicitly
checked, with the later C11/B12/D13 exceptions already discharged by
priority. For uses of the generic short8 route inside those later proofs,
the new corridor replaces a deep finite8 output by short7, and the new
Q=9/Q=10/tail cycle analysis replaces the old depth11 cycle interface.
No invocation asks the deleted distance8 pair itself to pay.

Every recent finite rank15/rank16/D17 decrease has distance at least10;
these and all existing coupons keep their exact current values. The new
deep pairs require components larger than those small ranks in any event.
The existing finite7/ideal5 triangle theorem, including protected ideal
6/7 and small whole-factor interfaces and the cycle11 triangle exclusion,
uses no finite distance8 payment. Hence triangles still receive2 and
quadrilaterals1 under the single simultaneous rule stated above.

## Verification scope

Final generation a91fd9 and separate fresh replay9f838f both terminated
exit0/PASS, obtaining corridor margin1875/11264 and the three cycle
minima1/18,1/10,1/11. Two peer agents completed full source/report
mathematical reads. Both checked the five missed-retained-contact cases,
the disconnected K-minus-c extra white, the strict A8 equality exclusion,
the shifted Schur bound, both small cycle domains, the fixed-side Cauchy
correction and the constant-profile tail. They also reread the current
ideal priorities and protected finite outputs. Both passed the simultaneous
current61 composition and triangle interface without a mathematical
correction. A prose typo was corrected. Neither reviewer reran the path
matrices,308 scalar values or inherited classification/census.
The source imports no project arithmetic module. Its exact checks certify
the finite matrix/profile inequalities and bind the completed inputs;
the all-label minimal-core order bound, profile classification, factor
fibers and actual-corner geometry remain explicit written theorem inputs.
No historical geometric census, forest optimization or ideal convolution
is repeated. No new dimension theorem, finite-volume example, external
refereeing or proof-assistant formalization is asserted. Direct source
bindings, Python syntax and local report links pass. Source, certificate
and report are final/frozen after these completed checks and reviews.

    python3 scripts/fresh_depth8_finite8_zero_transfer_20261018.py
    python3 scripts/fresh_depth8_finite8_zero_transfer_20261018.py --verify
