# Half payments on the affine11-cycle distance-five pairs

In the current payment rule, replace each of the11 distance5 unit pairs
in a whole varying affine11-cycle by a half. Its total pair cost decreases
from55 to99/2. Every affected quadrilateral still receives at least1, and
no triangle can contain such a marked cycle pair. No finite payment or
coupon is added or changed. Arbitrary supporting spectator diagrams and
all native labels are allowed.

The [source](../scripts/fresh_cycle11_half_pair_transfer_20261018.py) and
[certificate](fresh_cycle11_half_pair_transfer_20261018.json) corroborate
the proof below with exact matrices and scalar identities. This is a
local theorem. Other order11 affine families still have cost55, so the
all-family envelope does not yet improve. No new dimension bound or
finite-volume example is asserted.

## Complete actual-corner profiles

Write an actual quadrilateral star as W+a+b+c+d, with opposite pairs
{a,c},{b,d}. At the chosen ideal corner its actual varying affine factor
F is the11-cycle, marked at distance5; the remaining spectator diagram Z
is positive definite and orthogonal to F. These are the inherited
[actual factor-fiber facts](fresh_cusp_factor_payment_framework_20261017.md).
All roots have squared norm2. Both actual edge diagrams W+c and W+d
are strictly positive definite. Heights h=−c·nu and k=−d·nu at the
positive kernel nu of F are strictly positive.

An adjacent black contact or a common nonzero spectator coordinate of c,d
gives a marked path of length1 or2 in another actual corner. Such pairs
retain their full units. Hence suppose these alternatives are absent,
so the adjacent-black products vanish and the native nonnegative contact
vectors x,y to Z have disjoint supports.

Deleting either mark from F leaves A10. The complete positive-semidefinite
one-root extensions are isolated; one3 or4 bond at either endpoint; one3
bond at either second node; or3 bonds at both endpoints. This is the same
all-label path-extension classification used in the
[cycle12 transfer](fresh_cycle12_ideal5_transfer_20261017.md) and its
[weak-profile predecessor](fresh_distance8_ideal_corner_replacement_20261017.md).
No finite or affine exceptional diagram contains A10. The classical
classification gives precisely the listed possibilities, with no imposed
native-label cutoff. Adding spectators cannot repair an excluded principal
extension.

Each of the six single-contact profiles gives an adjacent marked path
of length at most6. The old endpoint opposite the attachment is a simple
leaf with a degree-two neighbor. Both degrees are fixed after all
spectators are restored because Z is orthogonal to F. No classical affine
component has that configuration; exceptional affine components are too
small. Thus the whole connected component containing the marks is finite,
and all other components are positive-definite principal parts of Z.
The adjacent corner is finite and its distance<=6 unit is unchanged.
The source checks all six physical-norm2 core matrices, strict positive
definiteness, marked distances and fixed leaf/neighbor roles.

If either profile has both endpoint3 contacts, its adjacent corner contains
a full affine11-cycle, obtained by replacing the deleted old mark by the
new black root. The positive kernel and actual corner positive
semidefiniteness force every further root to be orthogonal to this cycle:
its nonpositive products have zero positive-kernel sum. Hence it is the
whole varying affine factor, with the two marks again at distance5. This
different actual corner receives a half, as does the original corner.
Their two halves already total1. This argument does not require the other
opposite root's profile or height to be fixed.

Only two isolated profiles remain. They contact F only at their respective
opposite old marks. Their positive heights are therefore native contact
magnitudes:1, sqrt2, phi=(1+sqrt5)/2, or at least sqrt3. Parallel and dotted
opposite products are included in the last range.

## The three isolated height cases

The old arcs have lengths5 and6, so alpha=1/5+1/6=11/30. The inherited
Schur identity, reverified symbolically in the source, is

    0 = 1−2alpha/h²−2alpha/k²
          +alpha*(x/h−y/k)^T G_Z^(-1)(x/h−y/k).

The final term is nonnegative. If both heights are at least sqrt2, the
first three terms are at least4/15. If one height is1 and the other at
least sqrt3, they are at least1/45. Thus the only unordered height pairs
that remain are(1,1),(1,sqrt2),(1,phi).

For any nonzero x-coordinate i, disjoint supports and Gram-metric Cauchy
give

    (x/h−y/k)^T G_Z^(-1)(x/h−y/k) >= x_i²/(2h²) >=1/(2h²).

There is the corresponding bound for a y-coordinate. These bounds are
used one coordinate at a time; they are not added.

At(1,sqrt2), the required normalized spectator energy is3/11<1/2.
Consequently the height1-side support is empty. The other side would have
unscaled one-sided energy6/11. The complete
[strict spectrum below one](fresh_spectator_energy_below_one_20261017.md)
then forces a unit contact at a leaf of whole A_r with r/(r+1)=6/11,
so r=6/5, impossible.

At(1,phi), the required energy is25/11−sqrt5. It is positive but strictly
less than1/(2phi²), the smaller of the two possible nonempty-coordinate
bounds. Both supports must therefore be empty, which would give energy0.
This is again impossible. Reflection treats the reversed ordered heights.

At(1,1), the required unscaled energy is14/11. The next argument shows
that this necessarily supplies a marked c–d path of length at most5.

## Spectator energy14/11 forces a short marked path

Both separate augmented diagrams Z+c and Z+d are strictly positive
definite. Every finite Coxeter diagram is a forest, so each new root can
contact at most one node in any connected component of Z. If both roots
contact one component, the all-rank inverse argument in the
[five-quarters spectator theorem](fresh_spectator_five_quarters_transfer_20261017.md)
applies to that component. We use its independently stated intermediate
bound, not its conclusion specialized to energy5/4: contact nodes at
distance>=4 contribute at least4/3 to the difference energy. The proof
allows arbitrary native magnitudes, uses the strict inverse-diagonal
restrictions, and covers all classical ranks and exceptional types.
Every other component contributes nonnegative energy. Since14/11<4/3,
a mixed component must instead give contact distance<=3, hence a c–d
path of length<=5 with only spectator whites inside.

Suppose there is no mixed component. Every contacted component then
contributes a one-sided positive energy. Each is at least1/2, so there
are at most two active components. We need the following complete
one-component spectrum up to14/11:

    r/(r+1), r>=1 integral; 1; 6/5; 5/4.

Here is its all-rank derivation. There is only one contact within a component
by the forest argument. Its magnitude is1 or sqrt2: a native magnitude
at least phi already gives energy>=phi²/2>14/11. For magnitude sqrt2,
the inverse diagonal must be<=7/11<1. The inherited below-one subtree
argument forces a whole A_r leaf, and r/(r+1)<=7/11 forces r<=7/4.
Only A1 remains, with energy1.

For a unit contact on A_r, write N=r+1 and the inverse diagonal as
i(N−i)/N. A leaf gives r/(r+1). An interior position with both i,N−i>=3
has diagonal>=3/2, too large. The remaining positions have one of these
integers equal2, giving2−4/N; the energy bound forces N<=11/2. Thus
N=4 or5, giving1 or6/5.

For B_r, the all-rank inverse diagonals are the integers1,...,r−1 and
r/2, so only1 fits. For D_r they are the spine integers and r/4 at
the short leaves, leaving1 and5/4, the latter at D5. Exact exceptional
Grams show every diagonal of E6,E7,E8,F4,H3,H4 exceeds14/11. Finally,
for I2(m), m>=5, the diagonal2/(4−theta²) increases with the native bond
theta and already exceeds14/11 at theta=phi. The cases m=3,4 are A2,B2.
This proves the spectrum without a spectator-rank or label cutoff.

One active component cannot have energy14/11, since it is none of those
values. With two active components, each is at most14/11−1/2=17/22<1.
Each is therefore r/(r+1) with r<=3. The six unordered sums from
{1/2,2/3,3/4} are checked exactly; none equals14/11. No unmixed-support
case remains. The required short c–d path follows.

At its actual corner, a finite marked path<=5 retains a full unit. An
ideal marked path<=5 also retains a full unit, except possibly another
distance5 pair in an affine11-cycle; that one retains a half. In every
case this distinct corner contributes at least1/2. Together with the
original chosen half it gives the required total1. This completes the
quadrilateral proof. It never assumes that the original changed pair
still has a full unit.

## No triangle has an affected corner

For an actual triangle, write its third black root as c and put the
old affine11-cycle marks at a=0,b=5. Both other corner diagrams are
positive semidefinite. The A10 extension classification applied to F\a
allows c to contact only cycle nodes{1,2,9,10}; applied to F\b it allows
only{3,4,6,7}. These sets are disjoint. The old marks themselves are
interior positions in the opposite deleted paths, so contacts to a or b
are also forbidden. Hence c is orthogonal to all of F.

Its height at the positive kernel of F would then be zero, whereas the
third side of the actual triangle does not contain this ideal vertex
and has strictly positive height. This contradiction excludes every
triangle with the affected corner. Spectators cannot repair an excluded
principal PSD extension, so the conclusion has arbitrary-support scope.
Rotation/reflection treats all11 distance5 pairs. Triangle coverage is
therefore unchanged without importing an extra triangle coupon.

## Simultaneous-rule compatibility and verification

All alternative finite paths have length<=6 and all other terminal marked
paths have length<=5. Therefore this proof is compatible both with the
completed depth17 half rule and with a prospective depth11 zero rule on
finite distance8 pairs. No such distance8 unit or half is consumed here.
The other retained short5 ideal units and all coupons are unchanged.
The new cycle weights are invariant under its full dihedral automorphism
group. The11 half pairs and44 unit pairs have total99/2; the other
order11 families and their current envelope55 are not altered.

Generation95bef3 and separate fresh saved replay97502e both terminated
exit0/PASS. Two other agents completed full source/report mathematical
reads, including the inherited strict below-one spectrum and the
five-quarters proof's intermediate4/3 bound. They passed the complete
A10 extension profiles, arbitrary-rank/magnitude spectator classification,
native-height gates, whole-cycle positive-kernel sealing, distinct-half
payments, fixed-leaf finite outcomes and triangle contact-set exclusion.
Neither review duplicated the matrix/scalar arithmetic or any historical
census. No gap was identified in the stated local scope.
The source imports no project arithmetic helper. It
checks six physical single-profile matrices, a positive-kernel whole
cycle, the formal isolated determinant and all native-height gates,
the exact finite exceptional inverse diagonals, all six two-component
energy sums, the triangle contact-set exclusion and the complete weighted
pair partition. Classical all-rank inverse formulas, the finite/affine
classification, actual factor-fiber theorem and mixed-component4/3 bound
remain explicit inherited proof inputs. The finite exceptional checks
are complete cases; no finite numerical rank sample is used as an
arbitrary-rank proof. All15 direct bindings, seven local report links and
Python syntax pass. Source, data and report are final/frozen after the
recorded reviews. All reviews are internal; external refereeing and
proof-assistant formalization are absent.

    python3 scripts/fresh_cycle11_half_pair_transfer_20261018.py
    python3 scripts/fresh_cycle11_half_pair_transfer_20261018.py --verify
