# Finite distance7 at feature depth8: complete profiles and a sealed equality

The local rule deletes finite distance7 payments in whole A/B/D
components when both marks are at least eight steps from every endpoint,
branch and endpoint of a nonsimple bond. All other weights and coupons
are unchanged. The argument below strengthens the completed
[depth10 rule](fresh_depth10_finite7_zero_transfer_20261018.md).
Exact generation, a separate fresh replay and two independent internal
mathematical composition reviews pass; no new finite forest maximum
or dimension bound is asserted by this local packet.

The [producer](../scripts/fresh_depth8_finite7_profile_transfer_20261018.py)
binds the current
[order13 rule](fresh_order13_sparse_ambient60_frontier_20261018.md).
The proof uses a complete finite contact classification with arbitrary
opposite-wall products; it does not extrapolate sampled labels.

## Actual star and retained path

Use the inherited actual quadrilateral star W+a+b+c+d, with opposite
pairs a,c and b,d. Roots have square2. W and each W plus one black are
strictly positive definite; adjacent corners are positive semidefinite,
and opposite pairs with W are indefinite. The full Gram has at most
one negative direction. Assume no adjacent black pair has a white-
interior path of length at most6. Then c.b=d.a=c.d=0.

The newly changed finite corner supplies the retained path A24, numbered
0 through23, with a=8 and b=15. There are eight outward supporting whites
on either side. A feature exactly eight steps away may be at an old
endpoint; all bonds inside this path are simple. Extra supporting whites
and their contacts outside the retained path remain arbitrary.

On an eight-white arm the inverse diagonals at outward distances3 through6
are2,20/9,20/9,2. A nonzero contact at any of these sites violates strict
actual-edge positivity. Together with the no-short6 assumption this leaves
c's retained white contacts only at

    0,1,6,7,9,22,23.

Its contact at the opposite black a is independent, of magnitude t.
The d profile is the reflection i→23−i, with opposite product u at b.
Each of t,u is0 or at least1, including arbitrary parallel/dotted products.

## Complete seventeen white-contact profiles

The retained W principal splits as A8+A6+A8. Its smallest inverse
diagonal is6/7, and its inverse entries are nonnegative. Thus three
nonzero white contacts would already have projection at least18/7>2.
There are at most two. A label at least5 has contact square at least
(3+sqrt5)/2>5/2, giving projection greater than15/7>2. Each nonzero
white contact is consequently1 or sqrt2.

Exhaust all99 such vectors on the seven positions. Impose exact strict
positivity on W+c and semidefiniteness on the adjacent retained corner
(A24 minus a)+c. Both are single-root Schur tests against positive
path blocks. Exactly these seventeen profiles remain (unwritten
coefficients are1):

    empty;
    {23}, {23:sqrt2}, {22};
    {9}, {9,23}, {9:sqrt2};
    {7}, {7,23}, {7,9}, {7:sqrt2};
    {6}, {1};
    {0}, {0,23}, {0,9}, {0:sqrt2}.

These are necessary profiles, not claimed realizable chambers. The strict
and weak inequalities are kept distinct, and no unknown Gram entry is
set to zero. Outside contacts can only impose further conditions.

## A central contact is forced, including disconnected deletion

At least one c contact at6,7,a=8,9 is nonzero. To prove this, choose a
minimal indefinite principal K in W+a+c. The inherited complete theorem
bounds its order by10, and it contains a,c and is connected. The component
B through a in K minus c has at most9 roots. Feature depth8 makes this
component a simple retained path: reaching an old endpoint requires all
nine roots, leaving no room for a branch beyond that endpoint. B cannot
pass through the absent b. Connectedness forces c to contact B.

Suppose t=0 and no contact6,7,9 occurs. Rightmost contacts22/23 cannot
reach B. If c contacts0, the whole nine-root path0..8 is required, so K
has exactly those roots and c. The seventeen-profile list makes its
only contact in K a simple or doubled bond at0: K is finite A10 or B10.

The only other possibility is c's single unit contact at1. The path1..8
and c already use nine roots. With no extra root K is A9. At most one
extra root is possible. If it joins B in the selected finite component,
feature depth8 permits only extension to0 or9; K is D10 or A10. If it
is in a separate component of K minus c, connectedness requires a
contact to c. The strict edge contains the seven whites1..7 and c,
an A8 with c its endpoint. Its inverse diagonal8/9 forces the new
bond to be3 or4, so K is again A10 or B10. No further extra component
fits the order bound. Every alternative contradicts indefiniteness.
This argument does not assume K minus c connected. Reflection proves
the d central-contact statement at14,b=15,16,17.

## Complete Schur certificate and its unique equality

Let J be the exact A24 inverse. Extend each white profile by x_a=t and
y_b=u. Put q=x^T J x, r=y^T J y and s=x^T J y. The central-contact
statement and positive inverse entries imply q>2 and r>2. The full
principal Schur block is

    [[2−q, −s], [−s, 2−r]].

For every pair of the seventeen profiles, discard pairs already giving
a white-interior c-to-d path of length at most6. Split each opposite
product into t=0 or t=1+T, and u=0 or u=1+U, with T,U≥0. In the zero
branches require the central white contact just proved. Exactly581
branches remain. The determinant polynomial

    (q−2)(r−2)−s²

has nonnegative exact Q(sqrt2) coefficients in every branch. In580 its
constant coefficient is strictly positive, with exact minimum47/25.
Thus it gives two negative
directions, impossible in Lorentz space. The only remaining branch has
t=u=0 and x=e9,y=e14. Its determinant is zero, with q=r=6 and s=4.
This is an exhaustive coefficient certificate for unbounded T,U, not
a numerical parameter grid. The candidate finite profiles were themselves
exhaustively obtained from the native gap and strict edge bound.

## Equality seals the finite path and supplies a retained unit

The full A24+c+d Gram has inertia(24,1,1) and radical

    k=(J(e9−e14),1,−1),   k_0=1/5, k_23=−1/5.

Its represented vector is zero in actual Lorentz space: the span contains
a timelike vector, whose orthogonal complement is positive definite.
Every extra white therefore annihilates this radical.

If the actual selected finite component extends the old path, take its
first extra root z. It can attach only at endpoint0 or23. An interior
attachment creates a branch within seven steps of a or b; a second
attachment creates a cycle, both forbidden by the finite component and
the feature-depth assumption. By reflection take endpoint0, with bond
magnitude p. Projection onto A24 gives p²·24/25<2, hence p=1 or sqrt2.

For p=1, the eight outward whites and z form A9, with inverse diagonal
9/10 at z. Strict actual-edge positivity restricts x_z=−c.z and
y_z=−d.z to0,1,sqrt2. Radical annihilation requires

    y_z−x_z=1/5,

which none of the nine pairs satisfies. For p=sqrt2 the same white
principal is B9 with inverse diagonal9/2 at z; both contacts must vanish,
contradicting radical annihilation with nonzero endpoint contribution.
Thus there is no first attachment: the selected component is whole A24.

Every other supporting white is orthogonal to A24. Its radical relation
now says x_z=y_z; a nonzero value gives the adjacent c-z-d path of
length2. In the no-short6 branch both vanish. All extra whites are
therefore orthogonal to the entire core. In the corner with marks c,d,
the middle path c−9−...−14−d is a whole A8, with both marks leaves.
Its distance7 unit is retained: its feature depth is zero. It is distinct
from the deleted selected corner and pays the face without any new coupon.

## Current payments, triangles and counting scope

Outside equality, the proof forces a full-star short6 corner. The
completed depth11/depth10 arguments show that every current ideal
short6 alternative uses retained short paths, finite7/8 paths protected
by a feature within five steps, a fixed marked leaf, or separate coupons.
Reducing the threshold to8 preserves these outputs. The later order13/14
changes have ideal distances at least9; rank15–17 finite changes have
distances at least10. None changes these short6 replacements.

Every one of the14 elementary triangle distance7 payments has a feature
within five steps; the exceptional distance7 alternate has a fixed
whole-component marked leaf. These retain their units. Thus the new
actual-corner priority preserves triangle payment2 and quadrilateral
payment1 under the current simultaneous rule.

The new deleted pair counts are max(0,n−23), max(0,n−24) and
max(0,n−25) for A_n,B_n,D_n. All finite compensation is zero. A complete
new forest count, including every existing coupon and restoration, is a
separate obligation. No dimension59 conclusion is asserted here.

## Verification status

Generation session6477, terminal chunk93f980, and separate fresh replay
session88473, terminal chunk51cfa0, pass exit0. Both reconstruct all99
candidate white vectors,17 admitted profiles and581 coefficient branches.
All coefficient signs are determined in Q(sqrt2) by rational square
comparisons, without floating-point positivity tests. The producer also
checks central-contact finite boundaries, the full radical, every
endpoint-extension contact and the whole A8 alternative. Direct graph
counts corroborate the all-rank spine formulas in236 A/B/D graphs through
rank80. No finite-forest optimization or ambient comparison is performed.

The construction peer completed a full producer/report mathematical read,
including the disconnected minimal-core cases, exact polynomial domain,
Lorentz radical relation, equality sealing, current short6 and triangle
interfaces, and all-rank counting formulas. This was a source and proof
review, not an execution or arithmetic replay. The ideal-lane peer also
completed a full producer/report mathematical composition read, checking
the complete native reduction, disconnected-core argument, exact sign
and monotonicity logic, radical sealing and the distinct whole A8 unit.
That review additionally checked compatibility with the new order12 and
order16 transfers. It did not execute or repeat profile, polynomial,
matrix, graph or historical-census computations. The minimal-core classification
and current short6/triangle geometry remain explicit inherited inputs.
No external referee or proof-assistant validation is claimed.

    python3 scripts/fresh_depth8_finite7_profile_transfer_20261018.py --verify
