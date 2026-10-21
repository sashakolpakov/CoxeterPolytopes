# Halving the affine D12 fork-center pair

An actual ideal varying factor of type affine D12 has order13. Its two fork
centers are distance8 apart. Their pair may be reduced from1 to1/2 after
adding the following finite payments between separate whole components:

* either endpoint of a whole A2 against the fork center of a whole D10:
  weight1/2, total increment1 per component pair;
* the unique middle of a whole A3 against the fork center of a whole D10:
  weight1/2, total increment1/2 per component pair.

All other weights are retained. The affine D12 total falls from74 to147/2,
with normalized cost147/26. Every affected quadrilateral has a distinct
unit-paid corner, or a distinct new half-coupon corner which combines with
the retained original half. This is a universal local theorem, allowing
arbitrary native labels and arbitrary supporting spectators. Its finite
budget and any ambient dimension conclusion are separate obligations.

The [source](../scripts/fresh_order13_d_center_pair_half_transfer_20261018.py)
and [certificate](fresh_order13_d_center_pair_half_transfer_20261018.json)
verify the exact formulas, affine alternatives and finite component roles.

## Actual star and complete center profiles

Write the actual quadrilateral star as W+a+b+c+d, with opposite pairs
{a,c},{b,d}. All roots have square2; W and each actual edge W plus one
black root are positive definite, adjacent corners are positive
semidefinite, and the four-root residual Gram after eliminating W has
determinant zero. At the selected ideal corner W+a+b, the
[factor-fiber theorem](fresh_cusp_factor_payment_framework_20261017.md)
gives the whole varying affine factor F containing a,b and a positive
definite spectator diagram Z orthogonal to F. This formulation applies
also to nonsimple ambient cusps.

Use spine0--1--...--8, leaves9,10 at0, and leaves11,12 at8; a=0,b=8.
The positive kernel has coefficient1 on the spine and1/2 on the four
leaves. Its products with the opposite roots give positive heights h,k.
If b.c, a.d or c.d is negative, a distinct actual corner has a retained
distance1 unit. A common nonzero spectator contact of c,d gives a retained
distance2 unit. Hence assume these three black products vanish and the
nonnegative spectator vectors x,y have disjoint supports. Opposite products
a.c and b.d remain unrestricted: their magnitudes are zero or at least1;
if larger than1 they are at least sqrt2. This includes parallel and dotted
opposite sides.

Deleting a leaves a finite D10 block B on1,...,8,11,12, and the two
isolated old leaves9,10. In B inverse diagonals at nodes1,2 are1,2; every
other diagonal exceeds2 (the fork leaves have5/2). A single nonzero
contact therefore can only be simple or doubled at1, or simple at2.
There cannot be two contacts to B: they create a cycle in the augmented
connected diagram already containing the remote fork. Finite and affine
Coxeter classification permits no such diagram. Exceptional types are
too small to contain D10. Thus the complete positive-semidefinite
principal extensions B+c have precisely these four profiles, including
the isolated one. This is a classification argument, not a bounded label
search. In the doubled-at1 and simple-at2 cases, a positive affine kernel
forces every other root of the actual corner to be orthogonal to that
component. Consequently c has no contact with the old isolated leaves or
with Z in those cases.

Put t=-a.c and let u,v be the contact magnitudes to9,10. Eliminating
F minus a,b gives the following exact data, with R=8:

|Profile|c0|l|h|
|---|---:|---:|---|
|outer only|2-(u²+v²)/2|0|t+(u+v)/2|
|simple contact at1|1-(u²+v²)/2|1|t+1+(u+v)/2|
|doubled contact at1|0|sqrt2|t+sqrt2|
|simple contact at2|0|2|t+1|

Here the residual c square before Z elimination is c0+l²/8, and its
cross term with the far scaled mark is -l/8. Reflection gives the right
end data. The supports are separated, and every left/right projection
cross product is lm/8. The source checks all four profiles and all16
ordered cross identities directly in the physical norm-two Gram.

The last two profiles already give distinct retained **unit** corners.
They are whole affine B10 and D10 factors of order11, respectively, with
the old far mark a fork center and c a terminal. Their marked distances
are8 and7. The current
[B/C/D11 rule](fresh_order11_bcd_half_pair_transfer_20261018.md)
only halves B11 terminal/fork-leaf distance9 and D11 opposite fork-leaf
distance8 pairs; these center pairs remain units. The source checks the
four reflected literal affine components, positive kernels, marked roles
and the exact retained weights. The C11 endpoint-zero change does not
affect either output.

We may therefore discard the last two profiles after paying a distinct
unit. For the simple-at1 profile, the principal corner B+c plus its two
old leaf contacts is positive semidefinite, so0<=c0<=1. Since h>=1,
its cost g=2l/h+c0/h² is at most3. No discrete outer-label list is needed.
For the outer-only profile the strict actual edge gives u²+v²<4.
The following cases exhaust every native contact and opposite product:

* u=v=0: positive height forces t>=1, so g<=2;
* exactly one outer contact, t>=1: c0<=3/2 and h>=3/2, so g<=2/3;
* exactly one outer contact, t=0: g=8/u²-2, equal to6 only for a simple
  contact u=1, and at most2 otherwise because u>=sqrt2;
* two nonzero outer contacts: u,v>=1, hence c0<=1,h>=1 and g<=1.

Thus every unpaid end has g<=6; the unique extreme is a simple new leaf
on one old outer fork leaf, with t=0 and height1/2. Every nonextreme has
g<=3. These statements hold with all additional spectator contacts:
those only strengthen the positive-semidefinite restrictions already used.

## Determinant and arbitrary spectators

Write J=G_Z^(-1). Eliminating F minus a,b and then Z gives the exact
four-root determinant identity

    0 = 1-(g_c+g_d)/8 + E/8,
    E=(x/h-y/k)^T J (x/h-y/k) >= 0.

The source verifies this identity symbolically, including the signs of
the off-diagonal spectator contribution. Its use precedes any unscaled
native-contact estimate. Two nonextreme ends would give g_c+g_d<=6 and
a positive gap at least1/4, which is impossible.

Suppose exactly one end, say c, is extreme. Then h=1/2 and
E=g_d-2<=1. If some x_i is nonzero, disjoint supports give y_i=0.
Gram-metric coordinate Cauchy yields

    E >= (x_i/h-y_i/k)²/G_ii = 2*x_i² >= 2,

a contradiction. Hence x vanishes on the entire Z, including any
connected tails. In the actual adjacent corner W+b+c, deleting a
separates the old fork leaf and c as a **whole A2**, while the old far
component is a **whole D10** marked at its fork center b. The old other
fork leaf is an isolated A1. Every remaining spectator is orthogonal to
these components, since it was orthogonal to F and x=0. The new A2/D10
half plus the original ideal half gives1. No restriction on y or on
components incident only to the opposite root d is necessary: d is
absent from this chosen corner.

If both ends are extreme, h=k=1/2 and g_c+g_d=12, so

    (x-y)^T J (x-y)=1.

Both separate spectator extensions Z+c and Z+d are strictly positive
definite, being principal subdiagrams of actual edges. Thus the complete
[energy-one spectator theorem](fresh_spectator_energy_one_boundary_20261017.md)
applies with its native magnitudes and disjoint support hypotheses.

If either support is empty, choose that extreme root and obtain the same
whole A2/D10 half as above. This does not require classification of the
other support. If both are nonempty, the theorem gives exactly:

* one simple contact on each of two distinct whole A1 components;
* simple contacts at the two endpoints of one whole A3;
* simple contacts at the twin short leaves of one whole D_r,4<=r<=7.

In the last two cases the c-to-d path has length4, giving a distinct
retained unit. The theorem's whole-component and no-other-contact
conclusions are essential: they follow from equality in Gram-metric
Cauchy, the native contact gap, and the strict separate edge inequalities,
not from enumerating spectator ranks. Its written arbitrary-rank proof
is inherited here; the source corroborates the five mixed finite types.

In the two-whole-A1 case, c is simply attached to one old fork leaf and
to its isolated spectator. In W+b+c these three roots form a **whole A3**,
with c its unique middle. The far old block remains whole D10 with mark b
at the fork center. Every other spectator has zero c contact and is
orthogonal to this component; the other contacted A1 is isolated in this
corner because d is absent. Thus the narrower A3-middle/D10-center half
is sufficient. Combined with the retained ideal half it pays1. We do not
charge the two A3 endpoints and do not require both symmetric coupon
corners. Eight reflected literal corners in the certificate verify the
whole A2/A3 and D10 roles and finite signatures.

## Simultaneous rule and interfaces

The pair of fork centers is an intrinsic singleton orbit. The completed
[dimension63 rule](fresh_order11_order17_dimension63_20261018.md)
inherits affine D12 cost74 from the
[whole-leaf transfer](fresh_affine_d12_leaf_pair_transfer_20261017.md).
The [dimension64 composition](fresh_deep_finite8_dimension64_20261018.md)
records the exact paid set; the source binds these completed theorems,
checks the center unit, and changes only that weight to1/2. All ideal
distance at most7 pairs remain literal units.

For every actual quadrilateral, first ask whether an actual corner has
the changed D13 center pair. If so, apply this theorem to that corner.
It uses its retained half and a distinct coupon half, or a distinct
retained unit. If no such actual corner exists, every actual weight is
the inherited weight, apart from finite coupon increases, and the prior
whole-face proof applies unchanged. This priority is independent of which
corner an earlier proof had selected; no provenance of a selected witness
is substituted for an actual-corner statement.

The theorem is compatible with the independent
[C11 endpoint-zero rule](fresh_order11_c_endpoint_zero_transfer_20261018.md):
its only ideal outputs are the retained B/D11 center units or a path of
length at most4. The B12 center/terminal deletion likewise cannot change
these outputs. Short finite paths at most4 are unaffected by deep finite8
zeros, the D16/D17 exchanges, or the proposed sparse finite refinements.
The coupons concern whole D10, not whole D16 or D17. No existing finite
same-component distance8 unit is consumed. These observations concern
local compatibility, not an assertion that future finite refinements have
already been globally budgeted.

The [finite7/ideal5 triangle theorem](fresh_finite7_ideal5_triangle_payments_20261017.md)
uses finite distances at most7, ideal distances at most5, protected ideal
leaf distances6/7, or whole ideal factors of order at most9. It therefore
does not use this order13 distance8 center pair. Triangle minimum2 is
unchanged; quadrilateral minimum1 is proved above.

The only additional finite costs are exactly1 per whole A2/D10 component
pair and1/2 per whole A3/D10 component pair, with the stated middle and
fork-center roles. A separate complete finite-forest bound must include
them. No new global dimension bound follows from this packet alone.

## Verification scope

Generation04dc72 and a separate fresh saved replay6a0f6d terminated
exit0/PASS, checking two finite inverse
diagonal lists, four exact end profiles,16 cross identities, the formal
determinant identity and scalar endpoint values, four retained whole
affine unit corners, eight whole finite coupon corners, and five mixed
energy-one finite types. The profile classification and arbitrary-rank
spectator theorem are written mathematical inputs, not consequences of
the finite checks. Direct audit ed80c9 passes26 source bindings,10 local
links and Python syntax.

A separate full producer/report mathematical read passes the complete
D10 attachment classification, positive-kernel sealing of the affine
B10/D10 unit alternatives, all-native extreme/nonextreme bounds,
coordinate-Cauchy whole A2 conclusion without restricting the opposite
spectator vector, energy-one whole A3-middle alternative, determinant
signs, finite roles and actual-corner/triangle/current-zero interfaces.
That reviewer did not rerun formulas, inertia, profiles or spectator
classification; the arbitrary-rank energy-one theorem remains an explicit
inherited input. Parent separately completed a full source/report
mathematical read of the complete D10 attachment cases, positive-kernel
sealing, cost bounds, determinant sign, coordinate Cauchy, whole A2 and
A3-middle outputs and changed-corner compatibility; no formula or inertia
replay was duplicated. Source, data and report are final/frozen. No external
refereeing or proof-assistant formalization is claimed.

    python3 scripts/fresh_order13_d_center_pair_half_transfer_20261018.py
    python3 scripts/fresh_order13_d_center_pair_half_transfer_20261018.py --verify
