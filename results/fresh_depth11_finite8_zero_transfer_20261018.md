# Deep finite distance-eight payments can be deleted at feature depth11

Retain the current complete finite and ideal payment rule, except assign
weight zero to a marked distance8 pair in a whole finite A/B/D component
when both marks are at distance at least11 from every endpoint, branch
and endpoint of a nonsimple bond. Every other weight and coupon remains.
Then **every triangle receives at least2 and every quadrilateral at least1**.

This improves the completed depth17 half rule in two ways: the required
feature distance falls to11, and the affected pair receives zero rather
than one half. The [source](../scripts/fresh_depth11_finite8_zero_transfer_20261018.py)
and [certificate](fresh_depth11_finite8_zero_transfer_20261018.json) check the
new exact matrix and cycle inequalities. A separate finite budget is needed
for subsequent counting improvements. No bound below64 is claimed here.
Generation, fresh replay and a separate full source/report mathematical
review have passed. The verification scope is recorded below.

## A shorter forced corridor

Use the actual quadrilateral star W+a+b+c+d, with opposite pairs {a,c}
and {b,d}, and roots of squared norm2. The white diagram and W plus any
one black are positive definite, while W+a+c and W+b+d are indefinite.
These are the same actual-face inputs as in the completed
[universal depth17 theorem](fresh_universal_deep_finite8_corridor_20261018.md).

Suppose W+a+b is a finite corner with the stated distance8/depth11 pair.
Its whole finite tree contains two simple outward arms of11 whites. Keep
nine whites on either side of the eight-edge a-to-b path. These27 roots
form an induced A27, with a at position9 and b at17, indexed from zero.
The other25 roots are supporting whites. The two extra whites on each
outer arm remain available for the actual-edge estimate.

Assume that the full star has no path of length at most7 joining the two
opposite black pairs with only whites inside. In particular c.d, c.b and
d.a vanish. The own opposite contacts c.a and d.b are unrestricted.
For an A_n path its norm-two inverse is

    J_n[i,j]=(min(i,j)+1)*(n-max(i,j))/(n+1).

Any retained outer contact at distance j=3,...,9 from its old mark gives,
in an actual edge, projection square at least

    J_11[j-1,j-1]=j*(12-j)/12 >=9/4>2.

Every nonzero native contact has magnitude at least1, and all inverse
entries are positive, so further contacts cannot lower this estimate.
It contradicts strict edge positivity. Contacts farther along the middle
path, or within two steps beyond the opposite old mark, supply a path of
length at most7 and are excluded by the current assumption. Thus c and d
can contact the retained A27 only at

    c: 7,8,9,10;       d: 16,17,18,19.

Both contact vectors are nonzero. Indeed choose a minimal indefinite
principal K inside W+a+c. It contains both a,c. The published all-label
order bound is10, with the exact Chein/Vinberg classification bridge recorded
in the depth17 theorem. The connected component B through a in K minus c
has at most9 roots and radius at most8 from a. Since B lies in the actual
finite tree, feature depth11 makes it a simple subpath. It cannot pass
through b, which is absent from K. Nine outward whites and the seven inner
whites therefore contain B. Connectedness of K forces c to contact B.
The reflected argument applies to d. No core-overlap or finite-label
cutoff has entered.

## Exact uniform Schur inequality

Let x,y be the two nonnegative contact-magnitude vectors, J=J_27, and
q=x^T J x, r=y^T J y, s=x^T J y. Both vectors have a coordinate at least1,
so positive inverse entries give

    q,r >=40/7.

Put U=sum x_i*(i+1) on the left block and V=sum y_j*(27-j) on the right.
The cross block has s=U*V/28. On either same-side block,

    J[i,k] >=(17/308)*u_i*u_k,

where u_i=i+1 on the left and u_i=27-i on the right. Consequently

    s^2 <=(121/289)*q*r.

Eliminating A27 from its Gram with c,d leaves

    [2-q, -s; -s, 2-r].

Both diagonals are negative. Writing q=40/7+u and r=40/7+v gives the
following lower bound for its determinant:

    36/289 +(382/289)*(u+v) +(168/289)*u*v >0.

The Schur matrix is negative definite, contradicting the single negative
direction available in the actual Lorentz space. Hence the full-star
short7 alternative exists. The source checks the A11/A27 inverse identities,
all positive leading determinants, seven edge gates,32 same-block bounds,
16 cross factorizations and the exact positive determinant coefficients.
The estimate applies to arbitrary nonnegative native magnitudes and arbitrary
supporting spectators. It is not inferred from a sampled contact census.

## The cycle interface at the shorter depth

The completed [cycle interface](fresh_deep_finite8_cycle_interface_20261017.md)
identifies every finite distance8 output that could change. Endpoint4 has
a nonsimple bond at its new mark, and second-node3 has a nearby branch;
neither is deep. An endpoint3 output has a fixed whole-component leaf at
distance complementary_arc minus1 from the old mark. Thus cycles with arcs
(6,8), (7,8), or (8,Q) for8<=Q<=11 have unchanged outputs. There are six
such endpoint cases.

Only the universal family with arcs P=8,Q>=12 remains. A changed endpoint3
output forces an actual spectator neighbor of its new root. Otherwise that
root is itself a whole-component leaf. The predecessor's fixed-side
unequal-height Cauchy inequality then subtracts1/2 from that profile's norm
constant in the determinant lower bound. If both profiles require neighbors,
use one fixed side throughout; no two coordinate bounds are added.

For every Q>=12 the residual list consists of the same six profiles:
isolated, endpoint3 on the8-arc, endpoint3 or endpoint4 or second-node3 on
the Q-arc, and both-end3. The long-arc distances are Q,Q,Q-1, all at least11,
so no change of classification occurs in this tail. The corrected norm
constants and moments are nonnegative. The determinant lower bound is
concave in reciprocal heights and nonincreasing in1/Q. It suffices to check
Q=12 and the four height-rectangle vertices for every ordered profile pair.

All36 pairs and144 exact rational endpoint values are strictly positive.
The minimum is

    5/16-2/12=7/48.

It occurs only for the two endpoint3 profiles on the8-arc. This contradicts
the zero determinant of an actual quadrilateral star in the unpaid branch.
There is no newly surviving equality case. The source uses normalized
rational moments; for the endpoint4 profile the square root cancels its
reciprocal height exactly, as in the frozen independent scalar checker.

The argument never uses the weight of the changed distance8 pair. It
supplies retained alternatives or excludes the remaining configuration.
It therefore permits weight zero. The earlier half was not an indispensable
term in a two-corner payment sum. All other cycle ranges retain their proofs,
including whole-leaf9 and separate-component coupons. Small-cycle outputs
have a nearby fixed leaf or sealed small components. The classical-tree
transfer outputs finite paths of length at most7, persistent ideal features,
or the sealed whole A7 equality. These are all unchanged.

## Complete simultaneous local coverage

Use the actual-corner case split of the completed
[dimension64 composition](fresh_deep_finite8_dimension64_20261018.md).
If a quadrilateral has no newly changed finite corner, all its actual
weights equal the proved current values and its total is at least1.
Otherwise the new corridor supplies an actual adjacent-corner path of
length at most7. Finite short7 pairs remain units. In ideal corners the
current small-tree short7 units, elementary ideal units and complete tree
and cycle transfers remain valid by the preceding interface argument.
This supplies total at least1 without adding any part of the deleted
corner's old weight. No cyclic reuse of a changed finite payment occurs.

The existing24-factor inclusion certificate is unchanged because all ideal
weights remain unchanged. The completed finite7 triangle theorem uses no
finite distance8 unit, so every triangle still receives at least2. All
current cross-component coupons and all remaining same-component weights
keep their values. This proves the stated complete local rule.

## Verification status

Primary generation35f9a0 and separate fresh replay29c166 terminated
exit0/PASS. A separate full source/report mathematical review passed the
forced corridor, complete contact supports, nonzero-vector argument,
uniform Schur obstruction, constant-profile cycle tail and fixed-side
Cauchy correction. It checked that the actual-corner split never consumes
the deleted half and that all triangle, ideal and coupon interfaces remain.
No unresolved step was identified. The reviewer did not repeat the144
endpoint values or any inherited geometric census. The source imports no project
arithmetic module; its finite identities do not replace the stated inherited
classification and actual-face proofs. The short7 small-factor inclusion
check is inherited from the frozen dimension64 binder. No new finite-forest
optimization, ideal-envelope improvement, dimension theorem, polytope,
external refereeing or proof-assistant formalization is asserted. All direct
source bindings, syntax and local report links pass. Source, certificate
and report are final/frozen.

    python3 scripts/fresh_depth11_finite8_zero_transfer_20261018.py
    python3 scripts/fresh_depth11_finite8_zero_transfer_20261018.py --verify
