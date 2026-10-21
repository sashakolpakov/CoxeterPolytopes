# Four opposite-leaf transfers in affine D11

For an actual quadrilateral with a varying affine D11 factor of order12
at an ideal corner, each of the four opposite short-leaf pairs can be
paid elsewhere after upgrading the whole-A1/whole-D11-short-leaf finite
coupon from one half to one. The affine cost becomes62 instead of66.
The other payments are retained paths of length at most5, whole finite
components of order at most17, and inherited finite half-coupons.
There is no arithmeticity, simple-cusp or spectator-rank hypothesis.

The [source](../scripts/fresh_affine_d11_leaf_pair_transfer_20261017.py)
and [certificate](fresh_affine_d11_leaf_pair_transfer_20261017.json)
check the Schur identities and literal finite-corner configurations.
This is a universal local transfer. Its new coupon must be included in
any global finite budget; the A/B/C factors of order12 and any ambient
dimension consequence remain separate obligations.

## Actual star and end profiles

Use the inherited [factor-fiber setup](fresh_cusp_factor_payment_framework_20261017.md)
for W+a+b+c+d, with opposite pairs a,c and b,d. All roots have square2;
the white and actual-edge diagrams are positive definite, corner diagrams
are positive semidefinite, and the actual four-black Schur complement
has determinant zero. At the chosen cusp, W+a+b splits as the varying
affine factor F and a positive-definite spectator diagram Z orthogonal
to F. Opposite walls have positive cusp height.

An adjacent black contact already pays at distance1. A shared spectator
neighbor of c,d pays at distance2. Otherwise the native contact vectors
x,y with Z have disjoint supports, and b.c=a.d=c.d=0. The opposite
products a.c and b.d remain unrestricted native products.

Number the spine0--...--7, with leaves8,9 at0 and10,11 at7. Give the
spine roots cusp coefficient1 and the leaves1/2. The effective resistance
between opposite short leaves is7+2+2=11. Deleting a leaf leaves finite
D11. The complete positive-semidefinite one-root extension list is
isolated, simple attachment to its long endpoint (the old mate), doubled
attachment there, or simple attachment to the next node (the old center).
This follows from the finite/affine classification as in the
[preceding D12 transfer](fresh_affine_d12_leaf_pair_transfer_20261017.md):
exceptional types cannot contain D11, and multiple contacts would make
a cycle together with a surviving fork.

For native opposite magnitude t, the end data are:

|Profile|c0|Moment l|Height h|Cost g=2l/h+c0/h²|
|---|---:|---:|---|---|
|isolated|2|0|t/2|8/t²|
|mate3|1|1|(t+1)/2|4/(t+1)+4/(t+1)²|
|mate4|0|sqrt2|(t+sqrt2)/2|4sqrt2/(t+sqrt2)|
|center3|0|2|1+t/2|4/(1+t/2)|

The checker directly eliminates the white part for all32 end profiles
and64 ordered cross profiles. It verifies the resistance block, residual
square c0+l²/11, far contact l/11, cross term lm/11 and determinant identity

    det(S)/(h²k²)=1-(g_c+g_d)/11
                 +(x/h-y/k)^T G_Z^(-1)(x/h-y/k)/11.

Every end has cost at most8. Equality occurs exactly when the new root
is a simple leaf attached to either old short leaf, with height1/2:
isolated t=1 or mate3 t=0. All other profiles have cost at most4.
The separate diagrams Z+c and Z+d are positive definite as principal
subdiagrams of actual edges.

## A nonextreme end forces the extreme end's spectator contacts to vanish

Two nonextreme ends have total cost at most8<11 and cannot have zero
determinant. If one end is extreme, the other needs cost at least3.
Strict decrease of the displayed expressions and the native angle gaps
leave precisely these five nonextreme choices:

|Profile and opposite magnitude|Height|Cost|Required weighted spectator energy E|
|---|---|---|---|
|isolated, sqrt2|sqrt2/2|4|1|
|isolated, (1+sqrt5)/2|(1+sqrt5)/4|12-4sqrt5|9-4sqrt5|
|mate3, 1|1|3|0|
|mate4, 0|sqrt2/2|4|1|
|center3, 0|1|4|1|

For isolated roots, positive height excludes t=0; labels at least6
give cost at most8/3. For nonextreme mate3, t>=1 and equality3 requires
t=1. For mate4 or center3, any nonzero t>=1 gives cost below3. Thus
this list covers arbitrary native finite labels, parallel and dotted
opposite products.

Let x be the extreme end's contacts and y the nonextreme end's contacts.
The scaled extreme contacts are2x. The nonextreme height is at most1,
and E lies between0 and1. If x is nonempty and y is empty, inverse
nonnegativity and diagonal lower bound1/2 give E>=4/2=2, impossible.
If both supports are nonempty, take one node from each. Gram-metric
Cauchy and disjoint supports give

    E >= (2*x_i+y_j/h)^2/(4-2*Gij) > 9/8 > 1.

Here both native magnitudes are at least1, h<=1, and Gij>-2 in a finite
diagram. This is again impossible. Therefore x=0. When E=0, positive
definiteness also forces y=0.

Choose the adjacent corner containing the extreme new root and the far
old mark, excluding the nonextreme new root. Every spectator is orthogonal
to this chosen corner's displayed F-based components. If the extreme
root attaches to the selected old leaf, that leaf is absent: the new
root is a whole A1 and the far mark is a short leaf of whole D11. This
gets the upgraded unit coupon. If it attaches to the old mate, the
component is whole finite D12 and receives the inherited unit.
In particular the center3 case never requires repayment at a second
deleted affine pair.

## Both ends extreme

Both heights are1/2 and both costs8. The determinant identity now gives

    (x-y)^T G_Z^(-1)(x-y)=5/4.

If one support is empty, use the adjacent corner containing the
spectator-free new root and the far old mark. The same whole A1/D11
unit or whole D12 unit pays, regardless of the other contact vector.
The source checks all32 literal choices of old marks, leaf attachments
and excluded black. These checks concern the selected finite corner;
they do not claim that the displayed core without spectators is itself
a realizable full face star.

If both supports are nonempty, the
[five-quarters spectator theorem](fresh_spectator_five_quarters_transfer_20261017.md)
gives either a path of length at most5 between c,d or exactly one whole
A1 and one whole A3, contacted simply at leaves from opposite sides.
The short path retains its finite8 or ideal5 unit.

For A1+A3, the checker constructs all32 literal configurations, including
which new black contacts A1. Their18-root Grams have signature(16,1,1),
their white diagrams and all four edge diagrams are positive definite,
the old corner is affine, and all three new corners are finite. Every
other spectator is orthogonal to this core. If a new root attaches to
the selected old leaf, its adjacent corner has a whole A2 or A4 against
a short leaf of whole D11, receiving an inherited half. An attachment
to the mate instead joins the marks in a whole finite component of
order at most17 and receives a unit. The saved corner totals are at
least one in every case; the new full A1 coupon is not used here.

## Payment interface and validation

The removed pairs are(8,10),(8,11),(9,10),(9,11), all at distance9.
All other62 pairs remain at unit weight. The transfer applies to every
actual quadrilateral using such a pair, including low-core, selected-core,
protected and immediate-factor uses. The short paths and finite payments
used in its conclusion do not call back to a removed ideal pair.
Elementary triangle payments are not altered by this local theorem;
their retained inclusion is a separate global-composition check.

Generation and fresh separate-process saved replay both terminated
exit0/PASS on all32 end identities,64 cross identities, five
nonextreme boundary profiles,32 arbitrary one-sided corner checks and
32 A1+A3 configurations with96 finite corners. A separate full internal
source and mathematical review passed the complete native profile list,
scaled Cauchy argument, whole finite-component choice without recursion,
height-before-energy reduction and the sealed A1+A3 corner payments.
The reviewer also checked compatibility with the finite A17 endpoint
deletion: these new finite corners have at most16 roots. That review
did not duplicate the32/64 arithmetic calculations. Source, certificate
and report are final/frozen.
No ambient dimension improvement or new finite-volume example is
asserted by this package alone.
