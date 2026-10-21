# Feature depth10 suffices for the finite7 zero transfer

The [depth11 theorem](fresh_depth11_finite7_zero_transfer_20261018.md)
can be strengthened to feature depth10. In the completed dimension60
rule, delete finite distance7 payments in whole A/B/D components when
both marks are at least ten steps from every endpoint, branch and endpoint
of a nonsimple bond. Every triangle still receives at least2 and every
quadrilateral at least1. All other payments remain unchanged.

The [source](../scripts/fresh_depth10_finite7_zero_transfer_20261018.py)
and [certificate](fresh_depth10_finite7_zero_transfer_20261018.json) check
the new equality argument below. Final generation, fresh replay and two
full internal mathematical source/report reviews pass.

## The shorter retained path

Use the predecessor's actual quadrilateral star W+a+b+c+d, with opposite
pairs a,c and b,d. Retain eight outward whites on each side of the
distance7 pair a,b. They form A24, with a at8 and b at15. The full
outward arms have ten whites. Their A10 inverse diagonals at distances
3 through8 are

    24/11, 28/11, 30/11, 30/11, 28/11, 24/11.

All exceed2. Exactly the predecessor's strict actual-edge and
no-short6 arguments restrict c's retained contacts to6,7,8,9 and d's
to14,15,16,17. Both vectors are nonzero: a minimal indefinite core
has order at most10, so its component through a after deleting c
has radius at most8 and fits in the retained path. Feature depth10
keeps that component simple. Disconnected deletion is still allowed.

Put J=J24 and alpha=1/25, with factors u_i=i+1 on the left and
u_i=24-i on the right. Every residual entry

    C[i,j]=J[i,j]-u_i*u_j/25

is strictly positive. Its left diagonal is

    77/25, 72/25, 63/25, 2,

and its right diagonal is the reflection. Thus, for nonzero native
contact vectors x,y, with U=sum x_i*u_i and V=sum y_i*u_i,

    q-2 >= U^2/25 >0,
    r-2 >= V^2/25 >0,
    s=U*V/25.

The Schur determinant is at least zero. A strictly positive determinant
again gives two negative directions and is impossible. Unlike the
depth11 argument, equality must now be resolved explicitly.

## The unique equality contradicts a genuine supporting white

Equality requires x^T C x=y^T C y=2. Every nonzero native coordinate
is at least1; all residual entries are positive, and only one diagonal
on each side equals2. Hence

    x=e9, y=e14,
    q=r=6, s=4.

The Schur block is [[-4,-4],[-4,-4]]. The full principal Gram on A24+c+d
has inertia(24,1,1), and its radical is

    k=(J(e9-e14), 1, -1).

In particular k's old endpoint0 coefficient is1/5. The source verifies
every entry of the full Gram-kernel identity, not just its two-dimensional
Schur counterpart.

This radical is an actual linear relation among the roots. Indeed the
principal span contains a timelike vector. Any vector in its radical
is orthogonal to that timelike vector and has norm zero; the timelike
orthogonal complement in the ambient Lorentz space is positive definite.
The radical vector therefore maps to zero. It pairs to zero with every
other actual root, including roots outside this principal diagram.

Let z be the next outward supporting white on the left, at distance9
from a. It exists because the feature depth is at least10. Inside the
retained old path it has precisely the simple contact with endpoint0.
Let x_z=-c.z and y_z=-d.z. Pairing z with k forces

    -1/5 - x_z + y_z = 0.

Both x_z and y_z must be0 or1. To see this, use the full ten-white arm
in either actual edge. Its inverse diagonal at z is18/11. Any nonzero
nonsimple native contact has squared magnitude at least2, giving
projection norm at least36/11>2, contrary to strict edge positivity.
Additional arm contacts only increase this lower bound. The only
remaining finite native labels are2 and3, with magnitudes0 and1.

The required difference y_z-x_z=1/5 is therefore impossible. All four
allowed contact pairs are checked exactly. This excludes the unique
zero-determinant case without limiting spectator rank, other contacts
or the ambient Coxeter labels.

The full-star short6 alternative consequently holds at feature depth10.

## Retained payments and finite counting consequence

The predecessor's actual-corner case split applies unchanged: a newly
zero finite7 corner forces a short6 actual corner; its current ideal
cycle/tree alternatives supply a full payment. If no newly zero corner
occurs, the completed dimension60 rule applies pointwise.

All finite7 outputs used in that short6 argument have a fixed leaf five
steps from a mark or a branch one step from it. The existing finite8
outputs have their separately established depth8 protection. Small
sealed components and separate-component coupons are unchanged. Thus
reducing the new feature threshold from11 to10 loses no such alternative.

For triangles, all14 elementary distance7 payments have a persistent
feature within five steps of a mark. The sole exceptional distance7
alternate has a marked whole-component leaf. Their protections remain
strictly below10. The source checks these interfaces in the completed
predecessor certificate; it does not repeat the1,025-triple enumeration.

The deleted pair counts in A_n,B_n,D_n are respectively
max(0,n-27), max(0,n-28), max(0,n-29). This follows by placing the first
mark at least ten steps from the left endpoint and the second at least
ten from the right endpoint, doubled-bond endpoint, or fork center.
Both marks lie on the undecorated spine, seven steps apart.

This still leaves the complete ambient60 finite maximum at413:
the [current frontier](fresh_ambient60_payment_frontier_20261018.md)
gives upper413, deleting nonnegative payments cannot increase it,
and the actual formal forest2D17+D26 loses no pair and still attains413.
This conclusion needs no new DP. It is an obstruction of the payment
rule, not a claimed geometric realization. Curvature remains405; a
dimension59 theorem still needs further finite and ideal improvements.

## Verification scope

Generation `2c69dd` and separate fresh replay `039503` both terminated
exit0/PASS. Two peer agents completed full source/report mathematical
delta reads after their earlier conceptual reviews. Both checked the
unique equality, full radical and its actual Lorentz interpretation,
the genuine ninth white, complete zero/simple contact alternatives,
disconnected minimal-core interface, retained short6/triangle outputs,
all-rank pair counts and unchanged413 counting conclusion. Neither
repeated matrix, scalar, triangle-census, DP or classification arithmetic.
These are internal agent reviews, not external mathematician refereeing.
The producer checks A10/A24 inverses, strict arm bounds, every residual
block entry, unique equality positions, exact Schur data, the full
principal radical and all allowed contacts with the ninth white.
It shares the predecessor's simple rational path-inverse routine.
Minimal-core, actual-face and current short6/triangle arguments remain
explicit inherited theorem inputs. There is no new geometric census,
finite optimization, ideal convolution, external refereeing or formal
proof-assistant validation. No new finite-volume example is established.

    python3 scripts/fresh_depth10_finite7_zero_transfer_20261018.py
    python3 scripts/fresh_depth10_finite7_zero_transfer_20261018.py --verify
