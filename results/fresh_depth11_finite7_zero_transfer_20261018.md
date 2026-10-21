# Deep finite distance-seven payments can be deleted

Start with the complete simultaneous rule of the
[dimension60 proof](fresh_order11_order16_dimension60_20261018.md).
Set to zero precisely those finite distance7 pairs in whole A/B/D
components for which both marks have graph distance at least11 from
every endpoint, branch and endpoint of a nonsimple bond. Every triangle
still receives at least2 and every quadrilateral at least1. All other
finite weights, all ideal weights and all separate-component coupons
remain as in the completed rule, including the existing deep finite8 zeros.

The [source](../scripts/fresh_depth11_finite7_zero_transfer_20261018.py)
and [certificate](fresh_depth11_finite7_zero_transfer_20261018.json) check
the matrix and triangle interfaces below. Generation, fresh replay and
two full internal mathematical reviews pass. This is a local theorem;
it supplies no new dimension bound by itself.

## A retained A26 path forces a shorter actual corner

Use the inherited actual quadrilateral star W+a+b+c+d, with opposite
pairs a,c and b,d. All roots have squared norm2. W and W plus any one
black root are strictly positive definite; W+a+c and W+b+d are indefinite.
Every nonzero native contact magnitude is at least1. The full actual
Gram has at most one negative direction. These are the arbitrary-support
inputs of the [universal corridor theorem](fresh_universal_deep_finite8_corridor_20261018.md).

Suppose a,b form a newly zero finite pair. Their distance is7, and
feature depth11 supplies two outward simple arms of eleven whites.
Keep nine on each side. The retained principal diagram is A26, with
a at index9 and b at index16. A feature exactly eleven steps from a
mark is allowed: the eleven-white path up to it still has simple bonds.
The entire retained path lies inside the actual finite component.

Assume there is no white-interior path of length at most6 between an
adjacent black pair. Then c.d=c.b=d.a=0. The opposite products c.a
and d.b remain unrestricted. Write

    J_n[i,j]=(min(i,j)+1)*(n-max(i,j))/(n+1)

for the inverse norm-two A_n Gram. In an eleven-white outward arm,
the diagonal entries at outward distances3 through9 are

    9/4, 8/3, 35/12, 3, 35/12, 8/3, 9/4.

Each exceeds2. A nonzero c- or d-contact at one of these positions
would already project with squared norm greater than2 onto that white
arm, contradicting the positive-definite actual edge. Additional
contacts cannot reduce the bound because the path inverse is positive.

The no-short-path assumption excludes every remaining retained
c-contact except positions7,8,9,10: an inward contact at least two
steps from a gives a path to b of length at most6; contacts within
two outward steps of b do so as well. The reflected d-support is
15,16,17,18. Contacts outside the retained path remain unrestricted.

Both retained contact vectors are nonzero. Indeed choose a minimal
indefinite principal K inside W+a+c. It contains a,c, is connected,
and has at most10 roots by the inherited all-label minimal-core theorem,
including its disconnected-semidefinite-deletion bridge. The component
B through a in K minus c has at most9 roots and lies in the original
whole finite tree. Feature depth11 makes B a simple path, entirely
within eight steps of a; it cannot pass through the absent b. Hence
B lies in the retained A26. Connectedness forces c to contact B.
This argument does not assume K minus c is connected. Reflection gives
the same conclusion for d.

## A uniform negative-definite Schur block

Let x,y be the nonnegative retained contact vectors, J=J26, and put
q=x^T J x, r=y^T J y, s=x^T J y. On the left support set u_i=i+1;
on the right set u_i=26-i. Define

    alpha=122/3267,
    U=sum_left x_i*u_i, V=sum_right y_i*u_i.

On either four-position block the matrix

    C[i,j]=J[i,j]-alpha*u_i*u_j

is entrywise nonnegative and every diagonal is at least2. Each contact
vector has a native coordinate at least1. Thus x^T C x>=2 and y^T C y>=2.
The cross block factors exactly, giving

    q-2 >= alpha*U^2 >0,
    r-2 >= alpha*V^2 >0,
    s = U*V/27.

Since c.d=0, elimination of A26 leaves the Schur block

    [2-q, -s; -s, 2-r].

Its diagonals are negative and its determinant satisfies

    (q-2)*(r-2)-s^2
      >= (alpha^2-1/27^2)*U^2*V^2
       = U^2*V^2/43923
      >= 4096/43923 >0,

because U,V>=8. It is negative definite, producing two negative
directions in an actual Lorentz Gram matrix. This contradiction proves
the full-star short6 alternative. There is no spectator-rank, contact
label, contact-magnitude or supporting-extension cutoff.

The source checks the A11/A26 inverses and positive leading determinants,
all seven arm bounds,32 residual entries, eight diagonal bounds and16
cross factorizations. These finite checks verify formulas used for every
allowed nonnegative contact vector; they do not enumerate that infinite
domain.

## Short6 ideal corners still supply a complete payment

At a finite corner, the resulting path of length at most6 retains a
unit. At an ideal corner, factor-fiber geometry places its marks in one
varying affine factor. The current rule retains every distance-at-most5
unit except the distance5 halves on the11-cycle. Handle that exception
first using its [independent theorem](fresh_cycle11_half_pair_transfer_20261018.md):
its alternative finite paths have length at most6, or a distinct retained
half combines with the original half. No distance7 finite payment is used.

Small exceptional affine factors retain their whole-factor units.
Every current classical-tree change has distance at least8, so the
small-tree short6 units and radius-five protected-feature units are
unchanged. For a remaining unpaid distance6 corner, the only universal
cases are a cycle of order at least12 or a feature-free classical tree.

For cycles, the inherited profile classification gives adjacent finite
single-contact corners. A single endpoint3 attachment on an arc of
length L has marked distance L, and the opposite old endpoint remains
a whole-component leaf at distance other_arc minus1 from the other mark.
An endpoint4 attachment has a nonsimple bond at its new mark; a
second-node3 attachment has a branch one step from it. The latter two
can never produce a newly deep pair.

The [cycle12 proof](fresh_cycle12_ideal5_transfer_20261017.md) has only
finite single-contact paths of length at most6. The
[cycle13 proof](fresh_cycle13_ideal5_transfer_20261017.md) can give
distance7 on the longer arc, but its old mark is five steps from the
fixed leaf. The [cycle14 proof](fresh_cycle14_ideal5_transfer_20261017.md)
at arcs6/8 can give distance7 only through a second-node attachment,
which has the nearby branch. Its distance8 endpoint output is protected
by a fixed leaf five steps from the old mark, as already required by the
current deep8 rule. Thus every one of these payments is retained.

For cycles of order at least15 with shorter arc6, the other arc is at
least9. A single-contact distance7 output cannot occur: the short-arc
outputs have distances5 or6, while the long-arc outputs have distance
at least8. The possible distance8 output is a second-node profile on
the9-arc and has a fixed branch. The inherited
[P=6 local proof](fresh_finite8_ideal5_local_replacement_20261017.md)
also uses retained leaf9 outputs, small sealed components, and whole
A1/A17 coupons. These all keep their current weights. Its height and
spectator arguments are therefore unchanged; no new cycle tail is
assumed from sampled ranks.

In the same report's radius-five classical-tree proof at P=6, an inner
second-node3 profile supplies a path of length5. Every other non-equality
case has a strict determinant gap. The unique equality seals a whole
A7 component with marked distance6. These outputs are also unchanged.

Now consider any actual quadrilateral. If a newly zero finite7 corner
occurs, the corridor and the short6 alternatives just established supply
a complete payment under the new rule. If no such corner occurs, every
actual corner weight agrees pointwise with the completed dimension60
rule. Its existing coverage applies. This actual-corner case split
avoids recursively asking a deleted finite7 pair to pay.

## Triangle coverage

The source freshly reconstructs all1,025 high-core marked triples of
the [finite7 triangle theorem](fresh_finite7_ideal5_triangle_payments_20261017.md)
and checks exact agreement with its saved record. Every one of its14
elementary distance7 payments has a persistent branch or nonsimple-bond
feature within five steps of a marked endpoint. Restoring supporting
roots cannot erase that feature. At a finite corner such a pair cannot
be newly deep; if the actual distance shortens, it is also retained.

The four exceptional alternate lengths are6,6,6,7. In the last case,
marked triple(1,2,9), the alternate pair(1,9) has mark1 protected as a
leaf of the whole actual component by the inherited affine E8 kernel
argument. Its depth is zero, so its finite unit remains. The ideal
exceptional-factor payment is unchanged. Low cores already give two
paths of length at most5 by the unchanged all-label Steiner argument.
Thus every triangle retains two distinct units.

## Verification and remaining scope

Generation `2d1eb8` and separate fresh saved replay `a2aecd` both
terminated exit0/PASS. Two peer agents completed full mathematical
producer/report reads. They also reread the arbitrary-support minimal-core
bridge, including disconnected deletion, the current rule and relevant
short6 cycle/tree proofs, and the inherited triangle/E8 interfaces. Both
passed the support restrictions, nonzero-contact argument, shifted Schur
bound, noncircular actual-corner case split and retained triangle payments.
Neither repeated matrix, marked-triple, scalar, DP or classification
arithmetic. These are internal agent reviews. The producer uses its own rational path formulas
and reuses the inherited triangle graph/kernel routines. The minimal-core
order bound, actual-corner/factor-fiber geometry, complete affine extension
classification, and old short6 cycle/tree arguments remain explicit
written theorem inputs. No historical core classification, finite DP,
ideal convolution, external refereeing or proof-assistant check is repeated.

This changes the generic finite rule in sufficiently large components.
It does not remove the ambient60 obstruction2D17+D26, whose ranks are
too small for a depth11 distance7 pair. The completed dimension bound
therefore remains60 pending further local improvements and counting.
No new finite-volume polytope is constructed.

    python3 scripts/fresh_depth11_finite7_zero_transfer_20261018.py
    python3 scripts/fresh_depth11_finite7_zero_transfer_20261018.py --verify
