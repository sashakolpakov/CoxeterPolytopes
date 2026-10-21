# Removing the order12 B fork-center / doubled-terminal unit

For a whole varying affine B11 factor of order12, remove the pair joining
its fork center to its doubled terminal. Every affected quadrilateral is
paid elsewhere after adding a half-payment between either marked endpoint
of a whole finite A2 and the doubled terminal of a separate whole finite
B9. The order12 B cost falls from64 to63. The new finite component-pair
cost is1. A joint finite-forest bound remains a separate obligation.

This is a universal actual-face argument, with arbitrary supporting roots
and native labels. It does not restrict to the eight preferred-core
contexts in which the singleton occurs. The [source](../scripts/fresh_order12_b_center_terminal_transfer_20261018.py)
and [certificate](fresh_order12_b_center_terminal_transfer_20261018.json)
check its exact Schur identities, finite inverse diagonals and boundary
corners. No new ambient bound or finite-volume example is asserted.

## Complete endpoint profiles

Write an actual quadrilateral star as W+a+b+c+d, with opposite pairs
{a,c},{b,d}. The original ideal corner has whole varying factor F and
positive-definite spectator diagram Z, orthogonal to F. Actual edges
are strictly positive definite and actual corners positive semidefinite,
as in the inherited factor-fiber framework.

Number F by the spine0--...--9, with doubled bond8--9, and short leaves
10,11 at0. The removed marks are a=0,b=9. Normalize physical root squares
to2 and cusp coefficients to1 on0,...,8, to1/sqrt2 at9, and to1/2 at10,11.
The marked resistance is9. The singleton is invariant under the fork-leaf
swap and has graph distance9.

An adjacent black contact already supplies a retained unit at a different
corner. A common spectator contact of c and d similarly supplies a
distance2 unit. Otherwise their native nonnegative spectator contact
vectors x,y have disjoint supports, and all adjacent black contacts vanish.
Opposite products remain unrestricted: each magnitude is0 or at least1,
and a value greater than1 is at least sqrt2, including parallel/dotted cases.

Deleting a leaves finite B9 on1,...,9 and two isolated A1 leaves. The
inverse diagonals of B9 are1,2,3,...,8,9/2. Its complete one-root PSD
extension therefore has no contact, a simple or doubled contact at1,
or a simple contact at2. Two contacts would create a cycle containing
the old doubled bond, which no finite or affine Coxeter diagram permits.
Exceptional types are too small. The inverse diagonal test covers every
native label, including equality; it is not a finite label search.

If there is a doubled contact at1 or a simple contact at2, the resulting
affine component has a positive kernel. All additional contacts to the two
old A1 leaves and to spectators must vanish. If the contact at1 is simple,
Schur elimination of the whole B9 block requires u²+v²<=2 for the two
outer contact magnitudes u,v. Thus the possible outer contacts are none,
one simple, one doubled, or both simple. If there is no inner contact,
the actual edge requires u²+v²<4; no label enumeration is needed below.

Deleting b leaves finite D11. Its long endpoint8 has inverse diagonal1,
its next node7 has diagonal2, and every other inverse diagonal exceeds2.
Its four complete profiles are no contact, simple or doubled contact at8,
and simple contact at7. The same tree/classification and Schur arguments
exclude every other profile.

Put t=-a.c on the center side, or t=-b.d on the terminal side. The exact
endpoint data are:

|End and profile|c0|Moment l|Height h|
|---|---:|---:|---|
|Center, outer contacts only|2−(u²+v²)/2|0|t+(u+v)/2|
|Center, simple contact at1|1−(u²+v²)/2|1|t+1+(u+v)/2|
|Center, doubled contact at1|0|sqrt2|t+sqrt2|
|Center, simple contact at2|0|2|t+1|
|Terminal, isolated|2|0|t/sqrt2|
|Terminal, simple contact at8|1|1|t/sqrt2+1|
|Terminal, doubled contact at8|0|sqrt2|t/sqrt2+sqrt2|
|Terminal, simple contact at7|0|2|t/sqrt2+1|

Each endpoint cost is g=2l/h+c0/h². Eliminating the ten old supporting
roots and then Z gives the exact identity

    0 = 1−(g_c+g_d)/9 + E/9,
    E = (x/h−y/k)^T G_Z^(-1)(x/h−y/k) >=0.

The source checks all eight endpoint rows, all sixteen separated cross
products, the marked Schur matrix and the formal determinant identity.

## A unique center extreme and three terminal possibilities

The center cost is6 only for t=0 and a single unit outer contact. This
means c is simply attached to one of the old fork leaves10,11 and has
height1/2. Every other center profile costs at most4:

* Without inner or outer contacts, positive height gives t>=1 and g<=2.
* For one outer contact, t=0 gives g=8/u²−2, at most2 unless u=1.
  For t>=1, c0<=3/2 and h>=3/2 give g<=2/3.
* For two nonzero outer contacts, c0<=1 and h>=1 give g<=1.
* Inner simple profiles cost at most3 because 0<=c0<=1 and h>=1.
  Inner doubled and next-node profiles cost at most2 and4 respectively.

Every terminal profile costs at most4. A nonextreme center would therefore
give determinant at least1/9, even before E, which is impossible. The
center must be the stated extreme. The determinant then requires g_d>=3.
The native gap and strictly decreasing terminal formulas leave only:

|Terminal profile|Opposite magnitude|Cost|Height|E|
|---|---:|---:|---:|---:|
|isolated|1|4|1/sqrt2|1|
|simple at8|0|3|1|0|
|simple at7|0|4|1|1|

For the next-node profile, the actual corner W+a+d contains the whole
affine D11 factor of order12. Its positive kernel seals the factor. The
marks are its fork center and an opposite short leaf at distance8, an
unchanged unit under the [completed payment rule](fresh_order11_order17_dimension63_20261018.md).
This directly pays the face. The literal next-node template checks only
that affine corner; it is not asserted to realize an entire face star
without further spectators.

## The remaining spectator block is forced exactly

In the other two cases h=1/2 and E<=1. If x has a nonzero coordinate i,
disjoint support gives (x/h−y/k)_i=2x_i. Gram-metric Cauchy gives

    E >= (2x_i)²/G_Z(ii) >=2,

a contradiction. Thus x=0. For the simple-at8 case E=0 forces y=0 too.
All spectators are orthogonal to the entire displayed core. The actual
corner W+a+d has a whole finite D12 with marks at the fork center and
long endpoint, distance9. Its retained unit pays the face.

For the isolated case k=1/sqrt2 and E=1 give y^T G_Z^(-1)y=1/2. A nonzero
coordinate has native magnitude at least1. Equality in Cauchy forces
y_i=1 and y=(1/2)G_Z e_i. Nonnegative contacts and nonpositive off-diagonal
Gram entries force i to be an isolated whole A1 in Z, with every other
contact zero. Thus exactly one spectator e attaches simply to d, while
every remaining supporting root is orthogonal to F,c,d,e. Indirect tails
and bridges are excluded by this equality, not merely omitted from a search.

The two relevant finite corners are consequently whole components:

|Corner marks|Marked components and roles|Payment|
|---|---|---:|
|b,c|B9 doubled terminal; A2 endpoint|new1/2|
|a,d|D11 fork center; A2 endpoint|existing1/2|

The existing half is part of the [D11/D13 coupon theorem](fresh_affine_d13_long_pair_transfer_20261017.md).
These two distinct corners total1 without the removed ideal unit. The
third finite corner has separate A2 marked components and is not needed.
Both choices of the old fork leaf give the same roles.

The source reconstructs both isolated-spectator cases and both simple-at8
cases, their full signatures, strict white and edge diagrams, all corners
and actual finite component roles. Physical doubled bonds are handled by
a recorded rational congruence; scaled diagonal4 roots are not treated as
physical norm4 roots. The two next-node affine corners are checked separately.

## Simultaneous rule and finite cost

The [order12 B/C theorem](fresh_order12_bc_unit_coupon_transfer_20261017.md)
already leaves exactly64 B pairs. Deleting only(0,9) leaves63. All short8
and elementary triangle pairs remain. No triangle payment changes.

For quadrilaterals, give priority to any actual corner with the newly
removed B singleton and apply this universal theorem. Its outputs are
distance1/2 units, the whole affine D11 center/leaf distance8 unit, the
whole finite D12 center/endpoint distance9 unit, or the two stated whole
finite halves. None is changed by the order11 C endpoint-zero proposal,
the existing deep finite zeros or the finite D16/D17 exchange. If there
is no affected actual B corner, all its old payments agree with the
completed rule. This avoids an unproved provenance assumption for earlier
short-path or immediate alternatives that happened to use the singleton.

The new half has two eligible A2 marks and one B9 doubled terminal, adding
exactly1 per whole A2/B9 component pair. B9 can no longer be silently
replaced by the old rank9 representative in a finite budget; its new
interaction must be counted explicitly. The local theorem by itself
supplies neither that maximum nor a new dimension bound.

## Verification

Generation5661ec and a separate fresh replay494d67 both terminated exit0/PASS.
A separate full mathematical/source review passed the complete principal
extension domains, all-native cost bounds, R9 identity, Cauchy and isolated-A1
sealing, physical gauge interpretation, all whole-component payment roles
and simultaneous actual-corner/triangle interfaces. The reviewer did not
repeat the matrix calculations or perform a finite-forest optimization.
Finite/affine classification, the actual
factor-fiber setup, and the previously retained weights remain inherited
mathematical inputs. No geometric core census is required for this
universal proof, and no external review or proof-assistant formalization
is asserted. Source, certificate and report are final/frozen after these checks.
