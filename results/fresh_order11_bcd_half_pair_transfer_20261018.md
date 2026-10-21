# Coupled half payments for the order11 affine B, C and D factors

For whole varying affine factors of order11, halve the two terminal/fork
pairs in B10, the one terminal pair in C10, and the four opposite-fork
leaf pairs in D10. Their total pair costs become54,109/2,53 respectively.
Every affected quadrilateral still receives at least1, provided the finite
rule includes a half between a whole A1 and the doubled terminal of a
separate whole B10. This adds exactly1/2 per such component pair. No other
finite coupon is added. The proof allows arbitrary supporting spectators
and all native labels.

Together with the frozen [cycle11 theorem](fresh_cycle11_half_pair_transfer_20261018.md),
the four order11 costs are

|Family|Total pair cost|
|---|---:|
|A cycle|99/2|
|B|54|
|C|109/2|
|D|53|

The all-family normalized cost at affine rank10 is therefore109/22.
The [source](../scripts/fresh_order11_bcd_half_pair_transfer_20261018.py)
and [certificate](fresh_order11_bcd_half_pair_transfer_20261018.json)
check the exact identities and the complete finite boundary cases below.
The finite coupon budget and global dimension comparison are separate
obligations. No new dimension bound is asserted by this packet.

## Intrinsic marked pairs and actual-star setup

All roots have physical squared norm2. A simple bond has contact magnitude1
and a doubled bond magnitude sqrt2. The following canonical descriptions
specify automorphism-invariant pair sets; numbering only records the checks.

* C has path0,...,10, with doubled bonds0--1 and9--10. Halve(0,10),
  the two doubled terminal leaves, at distance10.
* B has spine0,...,8, doubled bond7--8, and fork leaves9,10 at0.
  Halve(8,9),(8,10), each joining the doubled terminal to a fork leaf,
  at distance9.
* D has spine0,...,6, fork leaves7,8 at0 and9,10 at6. Halve all four
  pairs between opposite forks, at distance8.

Write an actual quadrilateral star as W+a+b+c+d, with opposite pairs
{a,c},{b,d}. At the chosen affected ideal corner F is the whole varying
affine factor containing a,b. Its remaining supporting diagram Z is
positive definite and orthogonal to F. These are the inherited
[actual factor-fiber facts](fresh_cusp_factor_payment_framework_20261017.md).
Every actual edge W+one black is strictly positive definite; every actual
corner is positive semidefinite. The positive kernel nu of F gives strictly
positive heights h=−c·nu and k=−d·nu.

Any adjacent black contact, or a shared nonzero spectator coordinate of
c and d, gives a marked path of length1 or2 at another actual corner.
Those pairs retain units in the simultaneous rule, so they already finish
the payment. Henceforth adjacent black products vanish, and the native
nonnegative contact vectors x,y to Z have disjoint supports. This is a
branch of the proof, not an assumption about all quadrilaterals.

## Complete all-label end profiles

Deleting a marked leaf from any of these three factors leaves a connected
finite B10 or D10. Its simple long endpoint lies next to the omitted leaf
(or is the old mate at a fork). The only possible positive-semidefinite
one-root extensions, when restricted to this old ten-root component, are

1. no contact;
2. a simple contact at that long endpoint;
3. a doubled contact at that endpoint;
4. a simple contact at the next node.

Here is the complete classification argument. An extension with two old
contacts would contain a cycle. A connected finite Coxeter diagram is a
tree; the only connected affine diagrams with a cycle are affine A, which
cannot contain the old doubled bond or fork. Exceptional finite or affine
diagrams are too small to contain B10 or D10. For one contact, the native
inverse diagonals of either old component are1 at the long endpoint,2 at
the next node, and greater than2 everywhere else. The Schur condition
q² Jii<=2 therefore permits only q=1 or sqrt2 at the endpoint, and only
q=1 at the next node. This includes equality and imposes no label cutoff:
the next native magnitude after sqrt2 is phi, whose square exceeds2.
Adding spectator roots cannot repair a non-PSD principal extension.

The opposite product with the omitted old mark has magnitude t>=0 and
is not restricted by that corner extension. It can be zero, an ordinary
native magnitude, parallel or dotted. A positive nonzero native magnitude
is at least1; if it exceeds1 it is at least sqrt2. The isolated profile
requires t>0 because its height is positive. Call the four profiles
isolated, end3, end4, next3. At a fork end, end3/end4 attach to the mate
and next3 attaches to the fork center.

Choose kernel coefficients1 on the spine,1/sqrt2 at doubled terminal
leaves, and1/2 at short fork leaves. Eliminating the nine supporting roots
in F leaves the scaled old-mark matrix

    (1/10) [[1,−1],[−1,1]].

For each new root write its height h, moment l, and residual constant c0
as follows; its end cost is g=2l/h+c0/h².

|End|Profile|h|l|c0|
|---|---|---|---|---:|
|Heavy|isolated|t/sqrt2|0|2|
|Heavy|end3|t/sqrt2+1|1|1|
|Heavy|end4|t/sqrt2+sqrt2|sqrt2|0|
|Heavy|next3|t/sqrt2+1|2|0|
|Fork|isolated|t/2|0|2|
|Fork|end3|(t+1)/2|1|1|
|Fork|end4|(t+sqrt2)/2|sqrt2|0|
|Fork|next3|1+t/2|2|0|

The exact face Schur determinant, which vanishes because the actual
quadrilateral quotient has dimension3, gives

    0 = 1−(g_c+g_d)/10
          +(x/h−y/k)^T G_Z^(-1)(x/h−y/k)/10.             (1)

The source checks every physical Gram, the24 end rows,48 cross rows and
the symbolic determinant identity. The spectator quadratic form is
nonnegative; no limit on its rank or number of components was used.

## Affine alternatives are distinct second halves

Every end4 or next3 profile gives a whole affine order11 factor at the
adjacent corner replacing the omitted old mark by the new root. Its
positive kernel, together with corner PSD and nonpositive products,
forces every additional supporting root orthogonal to that factor.
Thus its type and marked roles remain valid in the whole actual corner.

|Original end|end4 target|next3 target|
|---|---|---|
|C heavy|C terminal pair|B terminal/fork pair|
|B heavy|B terminal/fork pair|D opposite-fork pair|
|B fork|C terminal pair|B terminal/fork pair|
|D fork|B terminal/fork pair|D opposite-fork pair|

All targets are exactly pairs retained at1/2 by this theorem. The original
corner retains1/2 and this different actual corner contributes1/2.
This is a simultaneous two-half payment, with no recursive appeal to
another transfer theorem. The source checks all12 oriented instances,
positive kernels, whole factor types, marked leaf roles and distances.

Only isolated and end3 profiles remain. A heavy isolated cost is at most4;
a heavy end3 cost is at most3. At a fork the maximum8 occurs exactly for
isolated t=1 or end3 t=0. In either case the height is1/2, and the root has
exactly one unit contact to the selected old leaf or its mate. Call these
two profiles extreme. Every other remaining fork profile has cost at most4:
an isolated nonextreme t is at least sqrt2, while an end3 nonextreme t is
at least1 and has cost at most3. These bounds follow directly from the
positive-denominator formulas, which decrease with t.

For C the total cost is at most8, contradicting(1). For B an extreme fork
is necessary. For D at least one end must be extreme. This leaves only
one extreme plus a nonextreme, or D with two extremes.

## One extreme end: coordinate Cauchy and whole finite roles

Let the extreme end have height1/2 and contact vector x, and let the other
end be nonextreme. Equation(1) requires

    E=(2x−y/k)^T G_Z^(-1)(2x−y/k)=g_other−2,  0<=E<=2.

For a nonzero x-coordinate i, disjoint supports and Gram-metric Cauchy give
E>=4x_i²/2>=2. If E<2, x is empty. At equality E=2, necessarily x_i=1
and equality in Cauchy forces

    2x−y/k = G_Z e_i.                                  (2)

The sign pattern in(2) makes i the only x-coordinate. If i has a neighbor
j in Z, then y_j>0, giving the c–i–j–d marked path of length3. That
unchanged full unit finishes the proof. Otherwise i is a whole isolated
A1 component, and(2) forces all remaining contacts of both roots to vanish.
Consequently the selected extreme root has either no spectator contact,
or exactly a unit contact to a whole A1, with every other supporting
component orthogonal to the selected corner components.

Choose the adjacent corner containing that extreme root and the far old
mark. If the extreme root attached to the omitted selected old leaf, it
lies in a separate whole A1 or A2, while the far old mark lies in the
whole B10 doubled terminal (B case), or a whole D10 short leaf (D case).
If it attached to the old mate, the marks are in one whole finite B11/B12
or D11/D12 component. Thus the selected corner supplies at least1/2:

* A1/B10 doubled terminal uses the new half in this theorem.
* A2/B10 doubled terminal uses the existing
  [typed coupon](fresh_finite8_d12_b10_coupon_budget_20261017.md).
* A1/A2 versus D10 short leaf uses the existing
  [small-component half](fresh_large_finite_d_leaf_half_transfer_20261017.md).
* The joined finite component has a retained full unit.

The original ideal half plus this distinct corner's payment totals at
least1. This argument explicitly determines whole components, so it does
not mistake an induced A1/A2 or B10/D10 subdiagram for a coupon role.

## Two D extremes: the complete spectator energy3/2 split

Both heights now equal1/2 and both end costs equal8. Equation(1) becomes

    (x−y)^T G_Z^(-1)(x−y)=3/2.                         (3)

If either support is empty, select that extreme root and use the
spectator-free corner just proved. Otherwise both supports are nonempty.
Both separate augmented diagrams Z+c and Z+d are strict PD principal
subdiagrams of actual edges. Therefore each root contacts at most one
node in any connected spectator component, since a second contact would
create a cycle in a finite diagram.

If a component has contacts from both sides at distance<=3, it gives a
c–d path of length<=5. Such an ideal or finite corner contributes at least
1/2 in the simultaneous rule: only the cycle11 distance5 orbit can have
half rather than a full unit. The original half then finishes the payment.

Suppose a mixed component has contact nodes at distance>=4. Apply the
all-rank inverse calculation and positive quadratic remainder in the
[five-quarters spectator proof](fresh_spectator_five_quarters_transfer_20261017.md).
We use its intermediate inequalities, not its conclusion for energy5/4.
Writing the inverse block as [[a,c],[c,b]], the contribution is
a p²+b q²−2cpq, with p,q>=1 and a,b<2 by strict separate extension.
The complete possibilities with energy<=3/2 are as follows.

|Whole component and roles|a|b|c|Unit energy|
|---|---:|---:|---:|---:|
|A5 endpoints|5/6|5/6|1/6|4/3|
|A6 endpoints|6/7|6/7|1/7|10/7|
|A7 endpoints|7/8|7/8|1/8|3/2|
|D6 long/short leaves|1|3/2|1/2|3/2|
|E6 eligible endpoints|4/3|4/3|2/3|4/3|

For completeness, in A_r put N=r+1 and let d be the contact distance.
The unit energy is d(N−d)/N>=2d/(d+2). Thus d>=7 is impossible. For
d=4,5,6, the condition energy<=3/2 respectively gives N<=32/5,50/7,8;
together with N>=d+2 these force N=6,7,8 and the endpoint pairs above.
B has no far pair of allowed contact positions. In D only a long/short
pair in D6 or D7 is possible, with energy r/4, so D7 is too large. E6
gives the displayed block; E7 has at most one allowable node and E8 none.
F,H and I2 have diameter<=3. These are all ranks and native labels.

In each displayed block, either nonunit contact is at least sqrt2.
The source checks that the gradients are nonnegative at(sqrt2,1) and
(1,sqrt2), and the quadratic remainder is positive semidefinite. The
minimum on either corresponding quadrant is therefore attained there,
and both values exceed3/2. All surviving contacts are unit.

Every active component contributes at least1/2. For a one-sided component
this follows from Jii>=1/2 and a native magnitude>=1, or from the complete
[spectrum below one](fresh_spectator_energy_below_one_20261017.md).
For a mixed component with disjoint contact nodes i,j, Gram-metric Cauchy
applied to e_i−e_j gives

    E >= (p+q)²/(4−2Gij) > 1/2,

because the positive-definite two-node block has Gij>−2. This lower bound
is used to count additional active components, not added coordinatewise
inside one component.

The unit A5,A6,E6 energies leave less than1/2 before the total3/2, so no
other component can supplement them; alone they do not satisfy(3).
Only a whole A7 endpoint pair or whole D6 long/short pair remains, with
all other spectator contacts zero. These are whole components of Z by
construction, not smaller subdiagrams chosen inside a larger spectator.

For A7 there are four choices of the two extreme attachments (selected
old leaf or mate). If at least one attaches to its selected leaf, the
c–d corner is finite with marked distance8 in whole A9 or A18 and has
an endpoint mark. If both attach to mates, an adjacent finite whole D18
corner has a marked leaf at distance9. Each is a retained full unit.
For D6 the two orientations and four old attachment choices give eight
literal cases. The strict edge/corner tests exclude six; each remaining
case has a finite whole E8 marked path of length6. The source checks all12
configurations with exact rational inertia, every actual edge, all four
corners, component roles and finite marked distances. Other untouched
spectator components are orthogonal direct summands and change none of
these checks. The D18 leaf9 unit is retained: the large-D leaf halves begin
at order20. The endpoint distance8 pairs cannot be depth11 pairs.

It remains to treat(3) with no mixed component. Every active component is
one-sided with energy>=1/2, so at most three are active. With three, all
have energy1/2; the below-one classification makes them whole A1 with unit
contacts. Both roots have support, so one of them contacts exactly one A1.
With two active components, each side has exactly one; the smaller energy
is at most3/4. The same complete classification makes that side a whole
A_r with a unit endpoint contact, r<=3. Select this extreme root.

If it attaches to its omitted old leaf, the adjacent corner has whole
A_(r+1), of rank<=4, versus a separate whole D10 short leaf, earning the
existing half. If it attaches to the mate, the adjacent marked component
is whole finite D_(11+r), of rank<=14, earning a unit. The source checks
all these cases for r=0,1,2,3, as well as the four B cases above:12 exact
finite corners with their positive pivots and whole-component roles.
This completes every support pattern in(3).

## Triangles and simultaneous later interfaces

All changed B/C/D pairs have distance at least8. The inherited
[finite7/ideal5 triangle theorem](fresh_finite7_ideal5_triangle_payments_20261017.md)
uses a finite marked path at most7, an ideal path at most5, a protected
ideal path6/7 near a feature, or a whole affine factor of order at most9.
None of those units is changed here. The distinct cycle11 distance5
halves do not affect triangles, by the frozen cycle11 triangle exclusion.
The new finite coupon only adds weight. Thus triangle payments remain
valid without a new triangle census or coupon.

For quadrilaterals, if no actual corner has a changed B/C/D pair, its old
coverage is unchanged (with any separate cycle11 change handled by that
theorem). If a changed B/C/D corner exists, apply the complete proof above
at that actual corner. Its immediate affine alternatives contribute a
distinct half directly. No step asks a second changed B/C/D pair to supply
a full unit, and no branch uses a refined order17 ideal orbit.

The finite outputs are whole components of rank<=14, whole A9/A18 endpoint
distance8, whole D18 leaf9, whole E8 path6, or the precisely typed small
A/D10 and A1/A2/B10 halves. They survive the depth11 zero rule on deep
finite distance8 pairs and the D16/D17 exchanges. In particular the proof
does not consume a generic whole-D16 or whole-D17 long-pair unit.
The cycle11 proof uses finite paths<=6 or ideal paths<=5; the latter
still contribute at least half after every change here. Conversely, the
short<=5 branch above also requires only half from a cycle11 corner.

The separate order17 half-coupon successor can safely use a short8
alternative: even if it is one of the newly changed D10 pairs, the
alternative corner retains1/2 and the original changed order17 corner
retains1/2. This observation is an interface statement, not a new proof
of that successor's other geometry. All new weighted sets are invariant
under the automorphism groups of their whole native factors.

## Verification scope

Generation35489d and separate fresh saved replay5544c5 both terminated
exit0/PASS. The exact source
checks24 end identities,48 cross identities,12 whole affine alternatives,
12 small finite corners, the complete five far inverse blocks and12 sealed
far spectator configurations. The all-rank arguments are the classification,
inverse formulas, positivity and support reasoning written above; bounded
matrix checks are not being used as a rank or label cutoff. No previous
selected-map census is needed for this universal local theorem.

The parent completed a full source/report mathematical read, passing the
complete B10/D10 one-root classification, physical R10 identities,
nonrecursive distinct affine halves, equality-Cauchy whole A1 sealing,
complete energy3/2 far blocks and nonunit exclusions, literal A7/D6
boundary cases, unmixed whole A_r selection, and triangle/depth11/finite
exchange/cycle11 interfaces. A second agent independently completed a
full source/report mathematical and interface read of the same arguments,
including the positive-quadratic remainder on each nonunit quadrant and
the order17 half interface. Neither review reran the matrix, profile,
scalar or inherited census arithmetic. Both found no concrete gap in
the stated universal local theorem.

Direct audit60a6fd passed all25 source bindings, all9 report links and
source syntax. Source and data were unchanged between generation and
fresh replay; only this verification account was finalized afterward.
Source, certificate and report are final/frozen. External refereeing,
formal verification, a finite budget or a new ambient bound are not
claimed by this local packet.

    python3 scripts/fresh_order11_bcd_half_pair_transfer_20261018.py
    python3 scripts/fresh_order11_bcd_half_pair_transfer_20261018.py --verify
