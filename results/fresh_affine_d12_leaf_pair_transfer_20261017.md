# Four opposite-fork-leaf transfers in affine D12

Let an actual quadrilateral of a finite-sided finite-volume hyperbolic Coxeter
polytope have a varying affine D12 factor F of order13 at an ideal corner.
If its two marks are short leaves at opposite forks, its ideal unit may be
removed after adding the full finite whole-A1/whole-D12-short-leaf coupon.
The quadrilateral is then paid by a retained short path, a whole finite
component of order at most17, this new full coupon, or two inherited halves.
All four cross-fork-leaf pairs are covered. There is no arithmetic, simple
cusp, facet-count or spectator-rank restriction.

The [source](../scripts/fresh_affine_d12_leaf_pair_transfer_20261017.py) and
[certificate](fresh_affine_d12_leaf_pair_transfer_20261017.json) verify the
exact Schur identities and boundary configurations. The additional coupon
fits the [joint finite8 budget](fresh_finite8_d12_b10_coupon_budget_20261017.md)
through ambient dimension75. This is a universal face-local transfer, so it
also covers the immediate affine-D12 alternatives in earlier proofs, not
only the saved preferred-core contexts. No new ambient bound is asserted here.

## Actual face and complete profiles

Use the inherited [factor-fiber setup](fresh_cusp_factor_payment_framework_20261017.md)
for a face star W+a+b+c+d, with opposite pairs a,c and b,d. The white diagram
and each white-plus-one-black edge diagram are strictly positive definite.
Every adjacent corner is positive semidefinite. At the chosen ideal corner,
W+a+b splits as F plus a positive-definite spectator diagram Z orthogonal
to F. The actual four-black Schur complement has determinant zero. These
properties hold also at a nonsimple ideal vertex.

Adjacent black contact already gives a path of length1 and a retained unit.
Likewise a shared spectator contact of c and d gives a path of length2.
Assume neither occurs. Thus b.c=a.d=c.d=0, and the nonnegative native
contact vectors x,y of c,d with Z have disjoint supports. The opposite
products a.c and b.d remain unrestricted native products, including
parallel and dotted cases. Every nonzero native magnitude is at least1;
if it exceeds1, it is at least sqrt2.

Number the spine of F as0--1--...--8, with short leaves9,10 at0 and11,12
at8. The marks are a in{9,10}, b in{11,12}. Give spine roots cusp
coefficient1 and short leaves coefficient1/2. The scaled affine Gram is a
Laplacian, and its resistance between a and b is8+2+2=12.

Deleting a short leaf leaves finite D12, with its long endpoint the mate
of that leaf and its second node the old fork center. At this rank its
complete one-root positive-semidefinite extensions are: an isolated root,
a simple bond to the long endpoint, a doubled bond to that endpoint, or a
simple bond to its second node. This is the same finite/affine
classification used in the [larger-D transfer](fresh_affine_d14_long_pair_transfer_20261017.md).
Exceptional finite or affine types cannot contain D12. Two contacts would
create a cycle together with an existing fork and cannot give a finite or
affine component. Thus this is an all-label classification, not a search
with a maximum label or supporting rank.

Put t=-a.c, and call the mate L and the old center0. Eliminating F minus
a,b gives the following complete left-end data; the right end is reflected.

|Contacts in F minus a|Residual constant c0|Moment l|Height h|Cost g=2l/h+c0/h²|
|---|---:|---:|---:|---|
|none|2|0|t/2|8/t²|
|L by3|1|1|(t+1)/2|4/(t+1)+4/(t+1)²|
|L by4|0|sqrt2|(t+sqrt2)/2|4sqrt2/(t+sqrt2)|
|0 by3|0|2|1+t/2|4/(1+t/2)|

The residual square is c0+l²/12, the contact with the scaled far mark is
-l/12, and the cross projection of the separated end profiles is lm/12.
The checker verifies these statements for all32 end profiles and64
ordered cross-profile pairs and verifies the formal determinant identity
including spectators:

    det(S)/(h²k²)=1-(g_c+g_d)/12
                 +(x/h-y/k)^T G_Z^(-1)(x/h-y/k)/12.

Each end has g<=8. Equality occurs exactly when c is simply attached to
one of the two old short leaves: either t=1 in the isolated profile or
t=0 in mate3. These two extreme profiles both have h=1/2. All other
profiles have g<=4. This follows from the native gap and the strictly
decreasing expressions in the table, with isolated t>=sqrt2,
nonextreme mate3 t>=1, and the last two profiles t>=0.
In the isolated profile h>0 forces t>0: an opposite wall does not contain
the chosen cusp. Thus the native lower bound t>=1 applies there.

## A nonextreme end: the zero boundary seals the spectators

If at least one end is nonextreme, g_c+g_d<=12. If it is smaller, the
normalized determinant is strictly positive even before the nonnegative
spectator correction, which is impossible. Equality requires one extreme
end of cost8 and one nonextreme end of cost4. The latter is precisely:
isolated with t=sqrt2; mate4 with t=0; or center3 with t=0.

The spectator correction must vanish. Positive definiteness gives x/h=y/k;
their disjoint supports then force x=y=0. Thus every spectator is orthogonal
to F and both new roots, and cannot enlarge a displayed corner component.

Choose the adjacent corner containing the extreme new root and the old
mark at the opposite end, excluding the nonextreme new root. If the extreme
root was attached to the selected old leaf, that selected leaf is absent:
the new root is a whole A1 and the other mark is a short leaf of whole D12.
This receives the new full coupon. If it was attached to the mate instead,
the long arm extends to whole finite D13, and both marks receive the old
whole-component unit. Adding the orthogonal positive-definite Z leaves
this actual corner finite.

There are48 literal cases: four old mark pairs, two choices of the extreme
profile, three nonextreme equality profiles, and two orientations. The
checker reconstructs every Gram, its signature(13,1,1), white and edge
positivity, all corner semidefiniteness, and the chosen whole finite
components. Doubled bonds are made rational by scaling only the relevant
new root by sqrt2; its diagonal is4 and its single bond entry is-2.
The chosen corner excludes that root, so its type and mark checks use
the unscaled simple-bond Gram directly.

## Both ends extreme: an arbitrary one-sided spectator needs no classification

Both new roots are simple leaves attached to old fork leaves and h=k=1/2.
The determinant identity becomes

    0=-1/3+D/3,  D=(x-y)^T G_Z^(-1)(x-y),

so D=1. This unscaled equality is used only after the heights have been
proved equal to1/2.

If y=0, choose the corner containing d and a, which excludes c. All of Z
is orthogonal to F and d, regardless of x, its support, or the ranks and
types of its components. Exactly as above, this corner has whole A1 and
D12 with the required short-leaf mark if d attaches to selected b; if it
attaches to the mate, it has whole D13 containing both marks. It is finite
and receives a unit. The case x=0 is reflected. No classification of the
nonzero one-sided contact vector is needed. The32 literal choices in the
checker verify this whole-component statement; their arbitrary spectator
block is handled by the stated orthogonality argument.

## Two nonempty supports: a short path or two sealed halves

The separate augmented diagrams Z+c and Z+d are strictly positive definite,
as principal subdiagrams of the actual edge diagrams. Therefore the
[complete energy-one boundary theorem](fresh_spectator_energy_one_boundary_20261017.md)
applies. With disjoint nonempty supports and D=1, the only possibilities
are unit contacts on two distinct whole A1 components, the two endpoints
of a whole A3, or the short fork leaves of a whole D4 through D7. Every
unmentioned spectator component has zero contacts.

For A3 or D4 through D7, the two contacted nodes have distance2. Adding c
and d gives a marked path of length4 at their adjacent actual corner.
This is retained by both the finite8 and ideal5 rules, independently of
any other supporting roots.

For two distinct A1 components, the literal core consists of F, c,d and
one separate simple spectator leaf on each of c,d. Every other spectator
is orthogonal to this17-root core. All16 choices of old marks and attachment
leaves have signature(15,1,1); their three new corners are finite. If the
new root attaches to the selected old leaf, the corresponding adjacent
corner has a whole A2 against a short leaf of whole D12 and receives the
inherited half. An attachment to the mate instead extends the long arm
to D14 and receives the old whole-component unit. The checker gives four
cases with two halves, eight with a unit and one half, and four with three
units. Each receives at least one, with no use of the new full A1 coupon.

This exhausts all possible spectator behavior. Every use of a finite coupon
has established whole actual components, rather than only a finite
subdiagram inside an ideal corner.

## Payment interface and verification

The four removed canonical pairs are(9,11),(9,12),(10,11),(10,12), each
of graph distance10. The frozen order13 diagnostic identifies exactly
these four as the difference between the78-pair D12 rule and its74 basic
pairs. All basic and elementary-triangle payments are therefore retained.
Since the transfer is universal, it covers every induced selected-core use
and the `inner_next3` immediate alternative of the older affine-D14 proof.
It does not need to enumerate those provenance cases again.

Generation and a fresh separate-process saved replay both terminated
exit0/PASS on all32 end identities,64 cross identities,48 nonextreme
configurations,32 arbitrary one-sided corner checks and16 two-A1
configurations with48 finite corners. A separate full internal mathematical
and source review passed the nonextreme sealing, arbitrary one-sided
shortcut, height reduction before the energy-one theorem, strict edge
hypotheses and whole finite component arguments. It did not repeat the96
literal configurations/corner checks. Direct source bindings, syntax and
report links pass. Source, data and report are final/frozen. No external
refereeing, proof-assistant formalization or new finite-volume construction
is claimed.

    python3 scripts/fresh_affine_d12_leaf_pair_transfer_20261017.py
    python3 scripts/fresh_affine_d12_leaf_pair_transfer_20261017.py --verify
