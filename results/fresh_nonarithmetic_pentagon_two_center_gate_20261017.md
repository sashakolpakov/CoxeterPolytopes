# Exact two-center obstruction for the changed pentagonal cavity

After deleting12,47,52 from the [exact63-wall control](fresh_nonarithmetic_release49_62_replacement_20261017.md),
the following specified compact disk closure is impossible: retain the
five boundary sides [7,8,65,64,30] and all other old walls; insert at most
two new actual facets of S6; keep the five boundary vertices and edges
distinct, with no old chord or old-only triangle inside the disk; and
make every endpoint of this patch finite and simple. The result is
independent of the new finite Coxeter labels and their number field.

It does not exclude three or more new facets, additional wall deletion
or deformation, a patch meeting additional old sides and changing this
boundary disk, or arbitrary finite-volume completion. In particular,
there is no new finite-volume example here.

## Complete reduction to mandatory caps

S6 has type H4+A2. Every additional wall that defines an actual finite
facet of its H3 face is right to H4 and has one of seven A2 profiles:
zero, or a single bond3,4,5 to either A2 root. At every required endpoint,
the five-node diagram A2 plus two old sides and the new cap must be
spherical. Complete finite-diagram classification supplies the domain.
A finite bond m>=6 can occur only in an isolated I2(m) component, so
symbolic label6 covers every such m, without a numerical cutoff.

For two interior vertices of a triangulated disk, the new vertices must
be joined: otherwise an old chord separates them, contrary to the
stated topology. Each receives one consecutive arc of the old boundary.
The arc sizes are1+4 or2+3. For2+3 there are five placements. For1+4,
the owner of four old edges meets all five old sides. A single interior
vertex also meets all five. These latter cases contain the mandatory
four-contact frame [64,30,7,8] excluded below. This reduction invokes
only the stipulated disk topology, not a theorem that every possible
completion must have that topology.

Four of the five2+3 placements require the respective four-contact
frames

    [65,64,30,7], [64,30,7,8], [30,7,8,65], [7,8,65,64].

Their finite-corner domains contain12,4,4,12 low recipes and2,1,1,2
high-label tails. Every one has zero A2 profile. Let J be the exact
inverse of its projected four-frame Gram and c the contact covector.
The frame spans S6-perpendicular, so the unique proposed normal has
square c^t J c. All32 low recipes differ strictly from2.

Each of the six tails has one nonzero contact −t, t=2cos(pi/m), m>=6.
Its square is a t^2 with a>0. Five are uniformly below2 or above2 on
the enlarged closed real interval [sqrt3,2]. The remaining tail in
[65,64,30,7] is nonright only
to64. Its product with retained wall8 is a strictly positive constant
times t, so every positive t is obtuse. No finite-angle recognition or
field extension is used for this rejection.

The remaining placement requires a cap meeting both edges30–7 and7–8.
Both endpoints contain an H4 component, forcing the sole recipe: zero
A2 profile and rightness to30,7,8. Their orthogonal complement inside
S6-perpendicular is one spacelike line. Projecting retained wall53
onto it gives a nonzero vector k with k^2>0. The two norm2 possibilities
are ±sqrt(2/k^2)k. One has a strictly positive product with retained
wall4; the other with retained wall53 (also with64,65). Both are excluded
by acuteness, without choosing a common-center orientation.

Thus every one-/two-center compact disk in the stated class has an
impossible mandatory wall. The next closure of this particular cavity
must use at least three new actual facets or change the retained
adjacency/topology further. This is a necessary structural step, not a
construction of such a larger patch.

## Exact computation and scope of verification

The [checker](../scripts/fresh_nonarithmetic_pentagon_two_center_gate_20261017.py)
and [certificate](fresh_nonarithmetic_pentagon_two_center_gate_20261017.json)
use the frozen outward rational-interval physical engine, exact inherited
literal contacts, and the frozen finite-diagram classification. The
four-frame inverses,32 strict norm inequalities, six complete real tail
checks and both all-right branch blockers are newly computed. Each
decision has separated rational endpoints. The [numerical scout](fresh_nonarithmetic_release12_47_52_pentagon_20261017.md)
is bound for provenance but supplies no acceptance sign.

All walls used by a blocker are still retained after12,47,52 and the
earlier49,62 deletions. In particular, deleting the old rigid12–52 bond
does not evade this new gate. Exact deletion cost remains seven old
finite vertices and four compact two-faces; no new chamber is promoted.

Generation64042 and fresh saved replay94989 both exited0 with exactly
the saved certificate. Source syntax, all direct source manifests and
local report links pass. Source, data and report are final and frozen.
This is one primary exact implementation with shared frozen arithmetic,
not an independently reimplemented interval audit.
