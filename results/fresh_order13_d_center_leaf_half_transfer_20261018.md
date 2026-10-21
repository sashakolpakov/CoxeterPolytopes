# Order13 D71 using existing whole-component halves

In a whole affine D12 factor of order13, lower the four center-to-opposite-
short-leaf pairs from1 to1/2. Number the spine0 through8, with leaves9,10
at0 and11,12 at8. The changed orbit is

    (0,11), (0,12), (8,9), (8,10).

Every pair has distance9. The factor total drops73 to **71**. No finite
weight or coupon changes. The earlier center pair(0,8) remains zero.
Every triangle still receives at least2 and every quadrilateral at
least1 under the simultaneous rule described below. The B and C families
at this order still cost73, so this theorem alone does not improve the
all-family envelope or the completed dimension bound60.

The [producer](../scripts/fresh_order13_d_center_leaf_half_transfer_20261018.py)
and [certificate](fresh_order13_d_center_leaf_half_transfer_20261018.json)
bind the [current fractional rule](fresh_d17_fractional_ambient60_frontier_20261018.md),
the [earlier center-zero theorem](fresh_order13_d_center_pair_zero_transfer_20261018.md),
the [complete center/leaf profile argument](fresh_affine_d14_long_pair_transfer_20261017.md),
and the [arbitrary-rank energy-one theorem](fresh_spectator_energy_one_boundary_20261017.md).
The argument is universal in the supporting rank and native labels; it
does not infer coverage from the selected overlap maps.

## Actual face and complete profiles

Write an actual quadrilateral star as W+a+b+c+d with opposite pairs
a,c and b,d. Normalize physical roots to square2. The supporting-white
diagram and each white-plus-one-black edge diagram are strictly positive
definite; adjacent corners are positive semidefinite. At the selected
ideal corner the varying whole affine factor F contains a,b, and its
orthogonal complement in that corner is a finite spectator diagram Z.
These are the inherited factor-fiber hypotheses, including nonsimple
ideal vertices.

By diagram automorphisms it suffices to take a=0 and b=11. Adjacent black
contacts give a retained distance1 unit; a common spectator contact of
c,d gives a distance2 unit. Otherwise their adjacent products vanish,
their nonnegative Z-contact vectors x,y have disjoint supports, and the
opposite products remain arbitrary native products, including parallel
and dotted products. Each nonzero magnitude is at least1; above1 the
next possibility is sqrt2. Both heights relative to the selected cusp
are strictly positive.

Give the spine cusp coefficient1 and each short leaf coefficient1/2.
The resistance between a,b is R=8+2=10. The finite principal diagram
F minus the center splits as two A1 roots and D10. The finite diagram
F minus b is D12. For a further root, either large D component permits
only isolation, a simple or doubled bond to its long endpoint, or a
simple bond to its second node. The inverse diagonals at those nodes
are1 and2; all other diagonals exceed2. A second contact creates a cycle
in a component with a fork. The finite/affine classification excludes
such a cycle and excludes exceptional types containing D10 or D12.
The strict white-plus-one-black edge condition is used where the corner
alone could have equality. This is an all-label classification, not a
maximum-label test.

For the center end, after allowing its two old isolated leaves, the
complete symbolic profiles are:

|Profile|Contacts|c0|Moment l|Height h|
|---|---|---|---|---|
|outer|a:t,9:u,10:v|2−(u²+v²)/2|0|t+(u+v)/2|
|inner3|a:t,1:1,9:u,10:v|1−(u²+v²)/2|1|t+1+(u+v)/2|
|inner4|a:t,1:sqrt2|0|sqrt2|t+sqrt2|
|next3|a:t,2:1|0|2|t+1|

The inner3 outer contacts are exactly none, one by3, one by4, or both
by3. These follow from the same finite/affine classification after the
large block has extended to D11. The outer profile admits every native
u,v allowed by the strict edge inequality u²+v²<4. The inner4 and next3
profiles are handled by the immediate actual-corner priority below.

For the leaf end, let t be its opposite product and12 its mate:

|Profile|Contacts|d0|Moment m|Height k|
|---|---|---|---|---|
|isolated|b:t|2|0|t/2|
|mate3|b:t,12:1|1|1|(t+1)/2|
|mate4|b:t,12:sqrt2|0|sqrt2|(t+sqrt2)/2|
|center3|b:t,8:1|0|2|1+t/2|

The producer checks the two principal inverse lists, all eight end
identities and all16 separated cross identities directly at R=10.
For costs g=2l/h+c0/h² and f=2m/k+d0/k², eliminating F and Z gives

    det(S)/(h²k²) = 1−(g+f)/10 + E/10,
    E = (x/h−y/k)^T G_Z^(-1)(x/h−y/k) ≥0.

The actual four-black Schur determinant is zero. The formal identity is
checked symbolically, including the spectator terms and unequal heights.

## The two immediate affine cases use current weights

At center inner4, the adjacent corner contains a whole affine B10 factor
of order11 with marks at its doubled terminal and opposite short leaf,
distance9. Its current literal weight is zero. We do **not** use an old
whole-order cutoff here. Apply the completed independent
[BC11 zero theorem](fresh_order11_bc_zero_transfer_20261018.md) first to
any actual changed B/C11 corner; that theorem supplies a full payment
using its own retained short paths and distinct finite/ideal alternatives.
Its outputs do not involve an order13 center/leaf pair.

At center next3, the adjacent whole affine D10 factor has its two marks
at opposite short leaves, distance8. Its literal weight is still1/2 in
the [coupled order11 rule](fresh_order11_bcd_half_pair_transfer_20261018.md).
It is a distinct corner from the selected D13 half, so the two halves
give1 directly. A positive kernel seals each immediate affine component
against every other acute root in that corner. The checker reconstructs
both whole components, their kernels, all canonical marked images and
their exact current weights.

## Cost bounds and the unequal-height branch

After the immediate cases, the center cost is at most6. Its unique
extreme is a simple new leaf c on either9 or10, with t=0 and h=1/2.
Every other center profile has g≤3. More precisely, a nonextreme center
with g≥2 has only these possibilities, up to exchanging9 and10:

|Profile|Height|Cost|
|---|---|---|
|outer isolated, t=1|1|2|
|outer one doubled contact, t=0|1/sqrt2|2|
|inner3 with no outer contact, t=0|1|3|

To justify this discrete list, the outer isolated cost is2/t². A single
outer contact with t=0 has cost8/u²−2; the native gap gives the extreme6
at u=1, cost2 at sqrt2, and smaller values for all larger labels. With
t≥1 and a nonzero outer contact its cost is at most2/3. Two outer
contacts give at most1. For inner3, any nonzero t gives at most5/4;
one outer simple contact gives at most14/9, a doubled contact less,
and two simple contacts at most1. Thus the displayed list covers all
nonextreme profiles with g≥2, and in particular **h≤1** there.

The leaf cost is at most8. The two extremes are a simple d attached to
the selected old leaf11 or to its mate12; both have k=1/2. Every other
leaf profile has f≤4, by the monotone table expressions and native gap.

If neither side is extreme, g+f≤7 and the normalized determinant is at
least3/10, impossible. If the center is extreme but the leaf is not,
g+f≤10. Equality is necessary and E=0. Positive definiteness and
disjoint supports then force x=y=0. Choose the adjacent corner that
contains c,b and omits a,d. It has whole A2 marked at c and whole D10
marked at the short leaf b. Their existing half plus the selected half
pays1. Orthogonal Z cannot enlarge or join these components.

Suppose instead that the leaf is extreme and the center is not. The
determinant equation gives E=g+8−10=g−2, so g≥2, h≤1 and E≤1.
The entire leaf-side vector y must vanish. If y is nonzero and x=0,
coordinate Cauchy gives E≥(1/k)²/2=2. If both supports are nonempty,
choose i in the support of y and j in the disjoint support of x.
Cauchy with e_i−e_j gives

    E ≥ (y_i/k+x_j/h)²/(4−2(G_Z)_ij) > 9/8 > 1.

Here the numerator is at least(2+1)², and strict positive definiteness
of the two-node spectator principal gives (G_Z)_ij>−2 and hence the
strict denominator bound<8. Both alternatives contradict E≤1.

With y=0, choose the adjacent corner containing a,d and omitting b,c.
It is independent of all contacts of c. If d attaches to b, it is a
whole A1 against the fork center of whole D12, with the existing half.
If d attaches to the mate, both marks lie in whole finite D13 and have
its unchanged unit. In either case the selected half gives total≥1.
Every remaining supporting white is in Z and orthogonal to F and d,
so the stated components are whole without a spectator-rank bound.

## Both extreme and the complete energy-one boundary

Now h=k=1/2 and g+f=14. The determinant identity yields

    D=(x−y)^T G_Z^(-1)(x−y)=1.

If one vector vanishes, choose the corresponding extreme root's corner
as above: center side gives whole A2/D10-shortleaf half; leaf side gives
whole A1/D12-center half or finite D13 unit. The other vector can have
arbitrary contacts and rank because its root is absent from this corner.

If both vectors are nonzero, the inherited complete energy-one theorem
applies: Z+c and Z+d are strict positive-definite principal subdiagrams
of the actual edges. With disjoint supports, either the contacts are
unit contacts on two distinct whole A1 components, or they are at the
two endpoints of A3 or the two short leaves of D4 through D7. The latter
cases give an actual c-to-d path of length4 and an unchanged unit.
All other spectator components have zero contacts.

For the two separate A1 case, choose the center-root corner c,b. Its
whole small component is A3 with c its middle, and the large component
is whole D10 with b a short leaf. The pre-existing small-component
short-leaf coupon pays1/2 for every mark of A3, so it pays this middle
mark. Together with the original half the total is1. No new A3-middle
coupon is introduced: that new role in the earlier center-zero theorem
concerned the **D10 fork center**, whereas the present large mark is a
**short leaf** already covered by the old whole-small rule.

The checker constructs the four center-side finite representatives
(two old leaves, with or without the isolated spectator), the two
leaf-side representatives, and the five mixed energy-one diagrams.
These corroborate the exact component roles. Arbitrary-rank completeness
comes from the inherited energy-one theorem and the analytic one-sided
argument, not from enumerating these representatives.

## Simultaneous payment interfaces

The only finite payments used are the existing whole A2/A3 against
whole D10 short-leaf halves, whole A1 against whole D12 fork-center half,
and the whole finite D13 unit. No newly reduced finite A15/B15/D15,
A16/D16, A17 or D17 pair is used. The separate new D16/A16 local work
therefore does not alter these outputs. No finite cost needs a new
optimization for this theorem alone.

Apply BC11's independent actual-corner priority first. If such a corner
does not occur, the inner4 branch is absent. For any changed actual D13
center/leaf corner the complete argument above then obtains1 from its
retained half and another existing half/unit, or directly from a short
distance1,2 or4 unit. No result invokes the proposed D13 theorem to pay
its own output. On faces with no changed actual corner the prior rule
agrees pointwise. This also covers protected-overlap and immediate uses
of the four pairs without a new finite selected-map census.

The old D13 center-zero theorem uses only order11 center/terminal units,
short4 or whole D10-center coupons. The current B/C completion proofs,
order15/17 transfers and finite zero rules retain their independent
short8→7→6 routes. Their universal tree/cycle replacements use smaller
ideal paths or sealed finite roles, while this theorem handles any
newly encountered D13 center/leaf corner directly. All its short outputs
have length at most4, so no deep finite7/8 zero is consumed.

Triangle payments are unchanged: all altered pairs have distance9,
beyond the elementary finite7/ideal5 triangle interface. The earlier
protected small-factor and fixed-leaf triangle exceptions do not use
these pairs. All short8 literal D13 weights are unchanged, including the
already-zero center pair at distance8.

## Verification and scope

Final generationdb9c64 and separate fresh saved replay83e0ba both
terminated exit0/PASS. The producer checks
exact inverse/profile/cross identities, the spectator determinant,
native boundary values, current immediate roles and finite components.
It binds the complete inherited profile and energy-one theorems.

Two independent full source/report mathematical reviews passed. The root
review included the inverse-diagonal guards and inherited energy-one and
below-one arguments, the complete g≥2 profile list, unequal-height
Cauchy, current immediate roles and six whole finite representatives.
The construction review additionally reread the large-D one-root profile
proof and energy-one checker; it checked that vanishing of the entire
extreme-side vector suffices despite arbitrary contacts of the absent
opposite root, and distinguished the existing short-leaf coupon from
the separate fork-center coupon. Both checked current priorities and
triangle protection. Neither repeated the formula, profile, inertia,
classification or historical census arithmetic. Direct immutable
bindings, report links and Python syntax passed after review. The
source, certificate and report are final/frozen.

No historical all-label census, finite forest optimization or global
dimension proof is rerun. No external or formal verification is claimed.

    python3 scripts/fresh_order13_d_center_leaf_half_transfer_20261018.py --verify
