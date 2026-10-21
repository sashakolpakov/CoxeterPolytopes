# Affine-carrier exclusion with the order-nine prohibition

**Final internally reviewed proof.** Let an actual quadrilateral
two-face have four finite vertices and zero payment at every corner under
Vinberg's finite-vertex rule. Suppose the ambient Coxeter diagram has no
connected affine subdiagram of order at most9. Then both indefinite
complements H−U and H−V in its private/shared decomposition are
affine-free. Simplicity of other ideal vertices is not assumed.

This strengthens the cutoff of the frozen
[private-tree proof](fresh_private_affine_tree_exclusion_20261017.md),
[cycle-carrier proof](fresh_affine_cycle_carrier_exclusion_20261017.md), and
[tree-carrier proof](fresh_affine_tree_carrier_exclusion_20261017.md).
The new ingredient is an affine-E7 leaf witness at shorter classical
orders. No predecessor is changed. This note alone does not prove
payment at an ideal corner or an ambient dimension bound.

## Actual roles and the only changes to the cutoff

Use opposite black pairs {u1,u2}, {v1,v2}; all other roots are white.
Every W+ui+vj is positive definite. The actual star has one indefinite
component H, negative index one and corank one. The frozen
[finite-quad audit](fresh_simple_finite_quad_vinberg_transfer_20261017.md)
proves these facts from compactness of this particular face, without
assuming the ambient polytope compact. Its private/shared decomposition
and published finite-weight lemmas are retained.

A connected affine subdiagram of H−U must contain both v-roots and
otherwise only white roots. By the order-nine prohibition it has order
m>=10, and is classical affine A/B/C/D. The exceptional affine diagrams
all have order at most9. The symmetric assertion holds for H−V.

The former lower bound m>=18 was used in the private/carrier arguments
only for the following purposes. Each has a replacement valid at m>=10:

| Use | Sufficient replacement |
|---|---|
| Exclude exceptional affine types | m>=10 |
| Find a small affine diagram after a single leaf attachment | The E7/E8 dichotomy below |
| Separate opposite end decorations in the resistance classification | R=m−1>=9, whereas every two-contact path must have R<8 |
| Exclude exceptional finite one-node extensions of a grounded B/D tree | Its order m−1>=9 already exceeds the finite exceptional maximum8 |
| Identify a continued two-step terminal profile as affine B/D | An exact graph identity at every m, not a size exclusion |
| Exclude a finite extension at the doubled end of B_m | m>=10 suffices; F4 has only four nodes |

In particular the old resistance-role table, including its two
white-contact exceptions, does not gain a new row at orders10–17.

## A small affine witness after every simple or doubled leaf

Normalize classical affine tree roots to square2. At a main node choose
the simple path between the two decorations, stopping before a doubled
terminal bond, or ending at one selected fork leaf. Its total number of
simple edges is m−3 in each of B/C/D. Let L,R be the available simple
lengths on the two sides of an attached simple leaf.

If L,R>=2, then L+R>=7. If both are at least3, the new leaf and two
length-three arms give affine E7, with arms (1,3,3). Otherwise the shorter
arm has length2 and the longer has length at least5, giving affine E8
with arms (1,2,5). These are induced subdiagrams: the old diagram is a
tree, and unused fork leaves and longer continuations are omitted.
Their orders are8 and9 respectively.

If L or R is below2, the same explicit near-decoration witnesses as in
the private-tree proof apply: the doubled three-arm star, its one-edge
extension, the F4 affine path at a doubled terminal, and the D4
four-leaf star at a fork center. At a fork leaf, the original fork gives
arms1 and2 after the new leaf is added; the remaining simple direction
has at least m−4>=6 edges, so a length-five arm gives affine E8.

For a doubled new leaf, a main node with two simple neighbors gives the
affine doubled three-arm star. A node incident to a doubled edge gives
the 4–4 affine C2 path. At a fork leaf, retain its center, its sibling,
the first main-path node, the contacted leaf and the new doubled leaf.
The five-node affine B4 witness has positive kernel
(sqrt2,sqrt2,sqrt2/2,sqrt2/2,1), in the order
(center,contacted leaf,sibling,first main,new leaf).

Thus every simple or doubled single leaf on a classical affine tree of
order at least10 contains an affine subdiagram of order at most9.
The same assertion for an affine A cycle uses arms (1,2,5); since m>=10,
the unused cycle arc separates the two selected endpoints, so the
selected tree is induced. A doubled leaf uses the doubled star.

These are all-rank constructions. The
[checker](../scripts/fresh_affine_carriers_order10_20261017.py) records
648 literal tree-leaf witnesses at orders10–17, sixteen cycle witnesses,
and their exact positive kernels as corroboration. It does not infer
the all-rank statement from those probes.

## Private diagrams at the new cutoff

For a private affine U, a shared neighbor z can contact only its two
u-roots. The actual shared path and the zero-payment hypothesis reduce
nonzero labels to3 or4, exactly as in the frozen proof: a5 bond puts the
retained marked pair in a paid H3/H4 component, and a higher finite bond
cannot admit the required connected spherical continuation.

One contact is excluded by the preceding leaf lemma. For two contacts,
the same grounded Green diagonal is nu_i² R_ij, with nu equal to1 on
main nodes, 1/sqrt2 at doubled terminals and1/2 at fork leaves. The
strict spherical inequalities are unchanged. Opposite ends now have
R>=9; even their smallest possible coefficient gives R/4>=9/4>2.
All other roles give the old adjacent-main affine A2, sibling-leaf
affine A3, or the unique remaining terminal triangle: a doubled terminal
a and its adjacent main b, with z.b=−1 and z.a equal to−1 or−sqrt2.

The finite diagram U−a has order at least9 and is B or D. Its only
connected finite one-node extension prolongs the long endpoint b by3.
Two shared prolongations give affine B/D (or a negative perturbation
of its positive kernel), contradicting the actual corner. Thus the
connector z is unique. The block U+z has inverse zz-entry zero, and its
Schur complement on the remaining roots X is exactly G_X. Actual
negative index one and corank one force a unique affine component of X;
both v-roots belong to it by the two finite-corner deletions.

That component also has order at least10. The same finite-extension and
single-leaf arguments force it to be precisely the private V, attached
through the two v-roots. Deleting the marked doubled terminals in U,V
joins two long finite B/D arms through z, giving affine C/B/D inside a
finite corner. This is the same all-rank graph contradiction as before.
Private affine cycles are excluded by the frozen cycle argument, whose
long-A deletion rank and induced E8 requirements now hold at m>=10.
Hence neither private diagram is affine under the order-nine prohibition.

## Carriers inside the larger complements

For an affine tree A in H−U, an external root with both marked-deletion
corners spherical has at most two A-contacts: any contact-pair path
must contain both v-marks, and three tree paths cannot all contain the
same two distinct vertices. For an opposite u-root, or a shared root
reached from a u-root, the same paid H3/H4 argument reduces labels to3/4.
Single contacts are excluded by the new leaf lemma.

For two contact nodes x,y at resistance distance L, grounding at a mark
at position s gives exactly

` q(s)=t_x² nu_x² s+t_y² nu_y²(L−s)<2. `

Both coefficients are at least1/4, so L<8. Since opposite decorations
have distance at least9, the complete frozen908-case role reduction
applies unchanged. Every non-small-affine survivor contacts a marked
v-root. Therefore the opposite u-roots, which are right to both v-roots,
are right to A.

Take the first connector z on a shortest u-to-A path. Its predecessor w
is right to A. The fork-switch exception is excluded by the four-node
affine doubled star in one finite corner. The other white-contact
exception starts x−4−p−3−y, with marks p,y and z simple to x,y. Deleting
p and adjoining w−3−z produces affine B/D with the unchanged remote
decoration. This is a principal diagram of an actual finite corner, so
its size need not be small. The sixteen new literal boundary-order
checks confirm this graph identity at m=10,...,17.

The surviving first connector is the terminal triangle. A−a has rank
at least9, so every further outside neighbor touching a surviving root
would be a second simple prolongation of b, again impossible in a finite
corner. Remaining sole contacts to a are excluded by the small leaf
witnesses or the paid high-bond component along the actual shared path.
The connector is unique, all A-whites are private, and A is exactly V.
The preceding private-diagram argument excludes it.

For an affine A-cycle carrier, the two old arcs between any two external
contacts force those contacts to be precisely the two v-marks. A single
contact is excluded by the new leaf lemma. The long A_(m−1) deletions
allow only endpoint3/4 and next-to-end3 attachments; small affine cycles
leave only a doubled-terminal B_m extension. Its shared-path predecessor
cannot attach there in a finite diagram. This repeats the frozen
cycle-carrier proof with every rank requirement explicitly lowered.

Thus H−U and H−V contain no affine diagrams. This result removes their
large-carrier obstruction using only the order-nine prohibition. The
later published finite-quadrilateral case audit and the small-order
ideal-corner payments are separate obligations before any new dimension
theorem can be composed.

## Verification scope

Generation passes the648 new tree witnesses,16 cycle witnesses,16
boundary continuation identities and the new9/4 resistance margin.
The old908-role table is bound as frozen data and is not recomputed.
The new proof uses the published finite/affine classification, inherited
actual-face and Schur arguments, and explicit all-rank constructions.
It is not an independent replay of its predecessor proofs. Generation
and a fresh saved replay pass. A separate root mathematical/source read
passes the simple-path E7/E8 dichotomy, grounded rank and resistance
thresholds, actual-path continuation and private closure. That read does
not duplicate the old proof/table arithmetic or the new literal probes.
Syntax, direct source hashes and local links pass. Source, data and report
are final and frozen; this is internal review, not external refereeing or
proof-assistant formalization.
