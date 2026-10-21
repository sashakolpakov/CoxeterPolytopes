# A full corridor settles the deep-pair question throughout the high domain

For every corner of the complete frozen 6,703-context high-overlap domain,
the following statement holds. If the actual whole finite A/B/D component
has its two marked roots at distance 8, and both marks have distance at least
20 from every endpoint, branch and endpoint of a nonsimple bond, then the
full actual supporting star has an opposite-black path of length at most 7.

This is a **short-path alternative**. Only in the branch with no such path
anywhere in the full supporting star does it exclude the deep corridor.
The absence of a path in the saved old core alone is insufficient. There
is no new pair weight, complete payment system or dimension bound here.

The [source](../scripts/fresh_full_a43_corridor_high_domain_20261017.py) and
[exact certificate](fresh_full_a43_corridor_high_domain_20261017.json)
cover all 26,812 corners of the
[inherited high domain](fresh_distance7_high_supporting_frontier_20261017.md),
including all four corners of the 251 formerly unresolved contexts and all
earlier structural outputs. No old principal matrix is required to be
nonsingular or Lorentzian. The proof instead uses a forced positive definite
path in the proposed actual completion.

## The forced path and every possible old embedding

Normalize all actual wall roots to squared norm 2. A feature of a finite
A/B/D diagram is a vertex of degree other than 2, or either endpoint of a
bond whose label exceeds 3. Depth at least 20 at each mark, together with
marked distance 8, gives two outward simple arms of length 20 and the simple
path joining the marks. All vertices other than the two marks are supporting
whites. The two omitted black roots do not belong to this actual corner.

Keep the first 17 outward whites on each arm. Together with the marks and
the seven whites between them, these form an induced A43 diagram. Number
its vertices 0 through 42; the left and right marks are at 17 and 25.
The unused three outward whites on each side are retained as geometric
inputs for the A20 contact estimates below, but are not needed in the final
matrix. A feature at outward distance exactly 20 is harmless: all bonds
inside the displayed A20 arm are simple.

The 342 old native cores are reconstructed from their frozen diagrams and
positive diagonal gauges. They have order at most 19, so any displayed
corner has at most 17 roots. Their gauges are 1 or 2; dividing by the square
roots of the gauges puts every relevant product in Q(sqrt(2)). The source
checks all four left/right choices in every old context and independently
reconstructs the displayed component partitions.

A displayed component containing a mark must be a simple path. An existing
branch or nonsimple bond remains a feature in the actual component and is
at distance at most 16 from the mark. An existing cycle cannot occur inside
a finite A/B/D tree. The checker saves one such witness for each rejected
nonsimple marked block.

If the two marks are in the same displayed path, that path must already
give distance 8: a tree cannot acquire a second, shorter route. Its embedding
in A43 is then fixed once the two marks are placed at 17 and 25. If the
marks lie in different displayed paths, orient each path from its outward
endpoint toward the eventual joining segment. Let t_L,t_R be the inward
distances from the marks to the chosen endpoints. The connecting segment
must contain exactly

    k = 7 - t_L - t_R

additional whites, with k at least 1. The two displayed components are
orthogonal, so no zero-white connection is possible. Every orientation
of both paths is enumerated, with singleton orientations counted once.
All retained embeddings fit inside A43; all old pair products are checked
against the corresponding A43 entries. This also covers paths joined through
old spectator roots: those roots may occupy positions not selected in the
reduced old embedding. Their identities and extra contacts are not assumed
away. The Gram tests use a necessary subset of actual roots only.

## All-label restrictions leave three white slots per omitted black

Write b for the omitted left black root and d for the omitted right black
root. The full-star no-short7 hypothesis implies b·d=0 and forbids any
contact from b to the right mark, or from d to the left mark. Their products
with their own marks are fixed native products and remain in the calculation.

For a norm-2 A20 arm, the inverse is positive and its diagonal at outward
distance j is j(21-j)/21. The omitted black and all these whites belong to
an actual edge diagram, which is strictly positive definite. Hence:

* At outward distance 1, only nonzero contacts of labels 3 or 4 are possible.
  A label at least 5 has squared magnitude greater than 5/2, and
  (20/21)(5/2)>2.
* At distance 2, only label 3 is possible, since (38/21)·2>2.
* At distances 3 through 18, even a unit contact gives projection greater
  than 2, so no contact is possible.

Additional nonpositive products cannot lower these projection bounds,
because the inverse is positive. Strictness excludes equality as well.
These statements have no finite label cutoff.

In the inner A7 segment, b contacting a white at distance j at least 2 from
the left mark gives a route to the right mark of length 9-j, at most 7.
Only the first inner white can therefore meet b in the no-short7 branch.
Its A7 inverse diagonal is 7/8, which excludes label 5 and every larger
label by the same squared-magnitude comparison. Contacts to the first two
outward whites beyond the opposite mark also give short paths; more distant
contacts are already excluded by the other A20 arm. The argument for d is
the reflection of this argument.

Thus b can contact, apart from its own mark at 17, only positions 16,15,18,
with allowed nonzero labels respectively (3 or 4), 3, (3 or 4). Root d has
the reflected slots 26,27,24. The source first checks every fixed old contact.
A forbidden one supplies either a recorded short path, an A20 projection
obstruction, or a strict edge label obstruction. Four orientations are
rejected by the last reason; the remaining fixed-contact exclusions use
the former two reasons.

One omitted black cannot contact both consecutive outer whites. The positive
A2 inverse is (1/3)[[2,1],[1,2]], and two unit contacts already have projection
square 2, contrary to strict edge positivity. The remaining slot products
are enumerated exactly as 0,1,sqrt(2), respecting every old fixed contact.
They are an overestimate of realizable completions, which suffices for an
obstruction.

## Two necessary Gram tests

For a norm-2 simple path of order n, the inverse entries, using indices from
0, are

    J_n[i,j] = (min(i,j)+1)(n-max(i,j))/(n+1).

The source verifies this inverse formula by multiplication for n=2,7,17,20,43.
Deleting the two marked roots from A43 gives three independent white blocks
A17+A7+A17. Let J_W be their inverse. For the positive contact vectors x,y
of b,d with the A43 path, strict positive definiteness of each actual edge
requires

    2 - x_W^T J_W x_W > 0,     2 - y_W^T J_W y_W > 0.

Every profile violating either strict inequality is impossible. Otherwise
eliminate the full positive definite A43 block. Because b·d=0, the remaining
two-by-two Schur matrix is

    [ 2 - x^T J_43 x       -x^T J_43 y     ]
    [   -x^T J_43 y      2 - y^T J_43 y    ].

In every surviving enumerated profile its first diagonal is strictly negative
and its determinant strictly positive. It therefore has two negative
directions. By congruence the full Gram matrix would have two negative
directions, impossible for actual vectors in a real Lorentz space. This
argument requires no assumption on the old core's signature. Adding further
actual roots or spectator contacts cannot repair a principal Gram matrix
with two negative directions.

All signs and products use the
[exact Q(sqrt(2)) helper](../scripts/fresh_quadratic_sqrt2_exact_20261017.py).
No floating-point tolerance or optimizer is used.

## Complete counts and scope

|Corner outcome|Corners|
|---|---:|
|Nonsimple marked block incompatible with depth 20|16,627|
|Connected marks at a fixed distance other than 8|2,187|
|No admissible distance-8 joining orientation|14|
|Complete corridor/contact calculation|7,984|
|Total|26,812|

The last row of cases has 13,644 complete path orientations. Fixed contacts
exclude 8,337 orientations. The remaining 5,307 use 379 distinct contact
templates, with every template's multiplicity saved. There are 4,983 distinct
joint contact profiles. With orientation multiplicities restored, 57,582
profiles fail strict actual-edge positivity and 50,137 give two negative
directions. None remains unresolved within this deep-corridor hypothesis.

This settles the proposed deep-pair alternative throughout the **high-overlap
domain**, including the formerly unresolved 251 contexts. It does not settle
those contexts for arbitrary pair shapes, and does not complete the finite-7
payment theorem. Low/disjoint cases, the actual retained pair maps and all
affected universal payment interfaces must still be checked before changing
a weight or importing a finite budget. The unrestricted dimension ceiling
remains 65. No new finite-volume polytope is constructed.

## Verification

The first full-domain pilot used direct, uncompressed contact enumeration
and found the same orientation-weighted outcomes. The authoritative producer
compresses identical fixed-slot domains and saves every one of the 4,983
distinct exact profile checks. Generation 81249 terminated exit0/PASS
(terminal chunk c22b97). Separate fresh saved replay 96652 terminated
exit0/PASS (terminal chunk 8076b4). An earlier producer run exposed a missing explicit
fixed-label rejection for four orientations; it was corrected before this
successful generation. The final calculation records these rejections.

Two other agents completed full producer/report mathematical reads. They
checked complete connected/disconnected embeddings, the old size bound and
persistent-feature argument, physical gauge normalization, all-label slot
completeness, arbitrary spectator contacts, strict white-edge positivity
and the Lorentz index-two contradiction. Neither repeated the 4,983-profile
arithmetic or the inherited overlap census.

A separate [saved-template checker](../scripts/check_full_a43_corridor_saved_20261017.py)
completed exit0/PASS (terminal chunk 8cba59). It imports no project module,
uses generic rational Gauss-Jordan inversion instead of the path-inverse
formula, and decides Q(sqrt(2)) signs with rational isolating intervals
instead of the producer's squared-norm sign comparison. It reconstructs
all 379 saved slot domains, checks all 4,983 edge/Schur profiles exactly,
and reconciles the 13,644 orientation multiplicities. A second agent read
the complete helper and confirmed this scope. The helper does not independently
reconstruct the old native cores, fixed-contact exclusion witnesses, path
embeddings, all-label geometry or original high-overlap census.

The inherited all-label high-overlap census is a frozen proof input; this
packet reconstructs its native matrices, corner partitions and complete new
corridor/contact domain, but does not repeat the original overlap census.

    python3 scripts/fresh_full_a43_corridor_high_domain_20261017.py
    python3 scripts/fresh_full_a43_corridor_high_domain_20261017.py --verify
