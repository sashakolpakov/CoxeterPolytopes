# The cycle interface survives deep finite-distance-eight halves

The inherited ideal-five cycle and classical-tree replacements remain valid
if unit payments at graph distance8 in whole finite A/B/D components are
reduced to halves when **both marks have distance at least t from every
endpoint, branch and endpoint of a nonsimple bond**, for either t=17 or
t=20. All other current payments, including whole-leaf9 payments and separate-component coupons,
retain their existing values.

This is a local interface theorem. It does not prove that the proposed
finite rule pays every quadrilateral, give its finite-forest maximum, or
improve the unrestricted dimension ceiling65. The
[source](../scripts/fresh_deep_finite8_cycle_interface_20261017.py) and
[certificate](fresh_deep_finite8_cycle_interface_20261017.json) extend the
[finite8 cycle replacement](fresh_finite8_cycle_replacement_20261017.md).
For each depth the additional calculation has36 ordered profile pairs
and144 exact rectangle vertices. The strict minima are29/144 at depth17
and73/336 at depth20. This local statement is independent of whether a
separate corridor theorem has been proved at either depth; it does not
itself assert a depth17 high-overlap classification.

## Exactly which old finite units can change

At an actual ideal cycle corner, write its varying affine factor as F and
its orthogonal positive-definite spectator diagram as Z. Let the old marked
arcs have lengths P,Q, with P<=Q. An adjacent black root c has one of the
complete inherited profiles: isolated, endpoint3, endpoint4, second-node3,
or both-end3. The classification and actual factor-fiber argument allow
arbitrarily many spectators and arbitrary native labels; they are not a
bounded list of full vertex diagrams.

For a single attachment at position j=1 or2 on an arc of length L, the
adjacent finite marked distance is L-j+1. The endpoint of the old path
opposite the attachment, and its degree-two neighbor, stay fixed because
Z is orthogonal to F. This whole-component simple leaf excludes an affine
completion in the relevant ranks. The adjacent corner is finite even
after restoring all spectators.

An endpoint4 profile has a doubled bond at c. A second-node3 profile has
an existing fork one step from c. Neither gives a changed deep pair,
regardless of any additional contacts. Their stronger inverse-diagonal
sealing, when used by the predecessor to pay a whole-leaf9 pair, is also
unchanged. Thus only endpoint3 with L=8 can lose a unit.

In this endpoint3 case the fixed opposite old endpoint is a leaf of the
whole actual component. Its distance from the other marked root is
exactly **complementary_arc-1**. Spectators cannot remove the leaf or
shorten this distance in a finite Coxeter tree. Therefore:

|Cycle arcs|Endpoint3 attachment|Why the distance8 pair is unchanged|
|---|---|---|
|(6,8)|8-arc|Old mark is5 from a whole-component leaf|
|(7,8)|8-arc|Old mark is6 from a whole-component leaf|
|(8,Q), 8<=Q<=t|8-arc|Old mark is Q-1<t from that leaf|
|(8,Q), Q>=t+1|8-arc|The sole newly residual family|

These are all possibilities in the current cycle replacements, whose
marked distances are6 through9. The source records the15 finite endpoint
cases at depth20 and12 at depth17 in the first three rows. At(8,8),
reflection treats either8-arc.
Cycle12 and cycle13 single-profile payments have distances at most6 and7;
the cycle14 distance8 profile is the first row above.

Consequently a full-star path of length6 or7 at an ideal cycle corner
still has its old replacement even when that ideal pair itself has
weight zero. The small-cycle equality outputs are sealed small whole
components or unchanged separate coupons. In the classical-tree branch,
the immediate finite outputs have distances at most7; the ideal outputs
have a persistent feature one step from a mark. Its sole additional
P=6 equality is the sealed whole A7 pair. These outputs cannot be changed
deep finite-distance8 pairs. This uses the complete
[ideal-five local transfer](fresh_finite8_ideal5_local_replacement_20261017.md)
and its cycle12/13/14 successors, rather than treating every short ideal
pair as a literal unit.

## The new endpoint requires an actual spectator neighbor

It remains to treat P=8,Q>=t+1, where t is17 or20. The corner produced by
endpoint3 on the8-arc starts with a finite A_(P+Q) path whose new endpoint is c. If c has
no spectator neighbor, it is a leaf of the whole component and is not
deep. A changed pair therefore forces at least one actual spectator
neighbor of c. This is an assertion about the whole actual component,
not its displayed A path. Indirect extension cannot begin without such
a direct neighbor.

Normalize every root to squared norm2. Write h,k>0 for the native heights
of the two opposite black roots c,d against the positive kernel of F,
and x,y>=0 for their contact-magnitude vectors with Z. Every nonzero
native magnitude is at least1. A shared nonzero coordinate gives an
actual c-white-d path of length2, which retains a full payment whether
that corner is finite or ideal. We can therefore suppose the supports
are disjoint.

Any direct adjacent-black contact already supplies an unchanged
distance1 payment. In the unpaid branch, the adjacent-black products
that vanish in the inherited determinant formula therefore still vanish.

If c is the newly residual endpoint3 profile, choose one spectator
neighbor i. Then x_i>=1,y_i=0. With G the positive-definite norm-two Gram
matrix of Z and z=x/h-y/k, the Gram-metric Cauchy inequality gives

    z^T G^(-1) z >= (z^T e_i)^2/(e_i^T G e_i)
                     = x_i^2/(2h^2) >= 1/(2h^2).

The corresponding bound for d is1/(2k^2). This does not require h=k,
does not restrict other spectator contacts, and does not use a rank or
label cutoff. If both profiles require neighbors, use the c-side bound
throughout the entire rectangle. The two coordinate bounds are not added,
and the chosen side does not change between rectangle vertices.

## Complete constant-profile tail and exact strict gap

Put A=1/8, B=1/Q, u=1/h, v=1/k, and use the inherited nonnegative
profile moments l1,l2,m1,m2 and norm constants c0,d0. Let
L=u*l1+v*m1 and M=u*l2+v*m2. The inherited normalized determinant before
spectator correction is

    F=1-2(A*L+B*M)-(A+B)(c0*u^2+d0*v^2)-A*B*(L-M)^2.

Eliminating Z adds (A+B) z^T G^(-1) z. The required-neighbor bound
therefore replaces c0=1 by1/2, or d0=1 by1/2 for the fixed d-side choice.
All modified norm constants remain nonnegative. The source checks the
corrected formula against the raw formula plus this explicit correction
at every saved vertex.

For every Q>=t+1 the complete residual profile list is the same:

    isolated, end3_0, end3_1, end4_1, next3_1, both_ends3.

The index0 denotes the8-arc. The old end4_0 and next3_0 profiles still
pay through their nondeep distance8 and distance7 corners. The new
end3_0 profile requires the spectator neighbor. On the Q-arc, the
three single profiles have distances Q,Q,Q-1, so they remain in the
old unrestricted list for all Q>=18. Isolated and both-end3 remain
unrestricted. Their inherited minimum heights are respectively
1,1,1,sqrt(2),1,2, independent of Q. Thus this is a genuine constant-list
tail, not an extrapolation from a changing finite sample.

With corrected constants nonnegative, F is concave in(u,v) on the
closed reciprocal-height rectangle. Its minimum occurs at a vertex.
It is also nonincreasing in B because the derivative is

    -2*M-(c0*u^2+d0*v^2)-A*(L-M)^2 <= 0.

Since B<=1/(t+1), it suffices to take Q=t+1 and the four rectangle vertices
for each of the36 ordered profile pairs. The source checks both tails:

|Feature depth t|Tail Q>=|Ordered pairs|Exact vertices|Strict minimum|
|---:|---:|---:|---:|---:|
|17|18|36|144|29/144|
|20|21|36|144|73/336|

All288 values are strictly positive. Both minima are attained for
end3_0/end3_0 with u=v=1.
A strictly positive lower bound contradicts the exact zero determinant
required by the actual quadrilateral star. No equality or additional
spectator case remains in this tail.

All other cycle ranges and the classical-tree branch retain their old
payments and proofs. Existing finite rank16/17 refinements preserve
the short8 and whole-leaf9 interfaces used here. The large finite-D
leaf9 half transfer also creates no new cycle dependency: a second-node3
distance9 output requires an arc of length10 and P<=9, hence a sealed
whole D component of order at most19. Its two changed large-D leaf
orbits start at order20. Existing separate coupons are unchanged.

## Verification scope

Final two-depth generation672f44 and separate fresh replay0277d7 both
terminated exit0/PASS on all27 protected endpoint cases,72 ordered
pairs and288 exact rectangle vertices. The earlier depth20-only run is
superseded. The producer deliberately shares the frozen exact arithmetic,
profile and determinant implementations.

Two other agents completed full source/report mathematical reads,
including the inherited cycle profile/classification and fixed-side
Cauchy interfaces. Their final parameter-delta reads passed both depths,
whole-component leaf/neighbor reasoning, complete constant-tail domains,
concavity, reciprocal-Q monotonicity and unchanged tree/small-cycle/leaf9
outputs. They did not re-enumerate the inherited all-label geometry.
One also checked the two minimum values algebraically at u=v=1,
without repeating the288-value producer computation.

A separate [saved scalar checker](../scripts/check_deep_finite8_cycle_saved_20261017.py)
terminated exit0/PASS, chunk576e4b. It imports no project module and uses
explicit normalized endpoint moments with rational Fraction arithmetic;
the sqrt(2) endpoint4 contact cancels its reciprocal height. It reconstructs
all288 raw values, corrections and lower bounds, all72 pair minima and
corrected constants, the27 endpoint cases and both complete scalar
domains. It independently obtains the minimum5/16-2/Q at Q=18 and21.
This is an independent scalar reconstruction, not an independent proof
of the geometric profile classification or the full payment system.
The complete helper has also received a separate source/mathematical read.

All direct input bindings and report links, and both Python syntax checks,
pass. No complete changed payment system, finite budget, external
refereeing or formal verification is claimed. Source, data and report
are final/frozen after the recorded reviews.

    python3 scripts/fresh_deep_finite8_cycle_interface_20261017.py
    python3 scripts/fresh_deep_finite8_cycle_interface_20261017.py --verify
