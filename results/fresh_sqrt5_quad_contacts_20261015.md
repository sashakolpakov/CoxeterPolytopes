# H4 payments and complete one-white contact domains over Q(√5)

The [exact checker](../scripts/fresh_sqrt5_quad_contacts_20261015.py) and
[certificate](fresh_sqrt5_quad_contacts_20261015.json) refine the 396 classes
left by the [principal-core filters](fresh_sqrt5_quad_filters_20261015.md).
They produce a finite reduction and two all-size consequences described
below. The remaining extensions are unresolved; **no new dimension bound
is claimed** by this package.

Use the weight rule from the [triangle package](fresh_sqrt5_triangle_weights_20261015.md):
base weight 1 for marked distance at most two, and an extra unit on the
endpoint pair of an entire H4 component. The initial union census already
assumes no marked distance-one or distance-two corner. All uses of “paid”
below mean that a quadrilateral already receives its required total of at
least one from this H4 rule.

## Reductions before adding a white node

Exactly 84 surviving core classes have an H4 component whose marked nodes
are its endpoints: 76 finite-label cores and eight dotted families. Such an
H4 cannot extend to any larger spherical connected diagram, so it remains
an entire component of the full finite corner. Its unit payment is valid
independently of all external white nodes and of the dotted distance.

Two further dotted classes, IDs 32 and 33, contain explicit affine
five-node principal diagrams. Each witness has exact inertia (4,0,1),
verified after the positive radical rescaling. Compactness excludes them.
These witnesses avoid the unresolved dotted edge.

Thus 286 finite-label and 24 dotted classes require further attention.
The finite-label classes consist of 112 singular and 174 nonsingular cores.

## Complete single-node domains

For each remaining core, append one white vertex with every possible
contact in the complete alphabet 2,3,4,5,6,10. Every pair containing a white
vertex belongs to a finite corner, so this alphabet covers all possibilities.
All four corner matrices must be positive definite, and no new marked path
of length at most two may appear. Adding a nonright bond cannot repair a
failure of these conditions, which justifies the monotone pruning.

There are 4,811 nonempty contact profiles after these structural tests:

| Outcome | Profiles |
|---|---:|
| Creates a sufficient entire-H4 payment | 62 |
| Contains an explicitly certified affine principal subdiagram | 508 |
| Fails the conjugate-core kernel condition | 817 |
| Cyclic product outside the ground field | 798 |
| Conjugate positive semidefiniteness fails | 98 |
| Forces two physical negative directions | 1,701 |
| Finite-label profile passing the necessary tests | 458 |
| Dotted-distance profile still unresolved | 369 |

Every saved affine witness is checked as positive semidefinite of corank
one. The search for these witnesses is deliberately used only in the
positive direction: absence of a found witness is not asserted to prove
absence of every affine subdiagram. For a singular conjugate core, any
new row of a positive-semidefinite extension must annihilate its nullspace.
The checker saves the exact kernel and every nonzero obstructing product.
Remaining finite-label profiles receive the full physical/conjugate Schur
tests from the preceding package.

Profiles already paid by H4 need no further arithmetic feasibility test
for this weight argument. Similarly, a row left unresolved in a dotted
family is not asserted to be arithmetically realizable.

## Two conclusions allowing arbitrarily many white nodes

All 112 singular finite-label cores have **no accepted unpaid nonempty
contact profile**. Consequently, in an actual quadrilateral that is still
unpaid, every external white vertex is orthogonal to every node of such a
core. Its affected finite corner components are entire. Arbitrarily large
orthogonal white components may remain. This assertion is conditional on
remaining unpaid: 19 of these cores have contact profiles that would instead
produce a sufficient H4 payment, so unconditional “no attachment” would
be too strong.

The same statement holds for 44 nonsingular finite-label cores. They cannot
be the principal core of an unpaid complete quadrilateral star. If all
external white roots are orthogonal to that nonsingular core, their Gram
block is positive definite because they lie in the spherical white support.
The entire star is then nonsingular. But a quadrilateral star with |W|+4
roots lies in the |W|+3 dimensional normal space of its two-face, and must
be singular. This contradiction excludes those 44 classes **from the
unpaid case**, without a cutoff on the external white support.

The remaining 130 nonsingular finite-label cores have 458 accepted profiles
between them. The 24 dotted families have 369 unresolved profiles. Multiple
new white vertices may interact, and their individual tests do not settle
the full star rank, conjugate positivity or dotted integrality. Those are
the next obligations. The 112 singular cores still require a suitable
weight assignment, or stronger global exclusions, to obtain a dimension
bound.

## Verification and scope

The exact one-node enumeration is a necessary-condition calculation, not
an enumeration of ambient polytopes. Its all-size consequences use the
first-contact argument and the full-star rank condition explicitly.
The source manifest binds the complete union census and the independent
core arithmetic package. Independent contact-domain review is being
prepared separately. All files from the preceding Q(√2) checkpoint remain
unchanged.

Replay: `python3 scripts/fresh_sqrt5_quad_contacts_20261015.py --verify`.
