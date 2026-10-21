# A complete E7 face and obstruction-guided compact search

**2026-10-14. No new finite-volume hyperbolic Coxeter polytope or dimension is established.**

The new H20 branch completes a fifteen-facet H13 face and repairs both of
its ambient cusp links. Its final chamber has 87 actual facets and 44
verified finite vertices on that face. An exact global escape still proves
infinite volume. The compact H9 search checks 74 further realizable
fifteen-facet incidence patterns; necessary Coxeter conditions exclude all
of them. The sources and previous archives remain preserved.

## Literature restrictions remain in force

The user's dimension concern is correct for noncompact arithmetic groups:
they occur only in dimensions at most 19 or exactly 21. In particular,
there is no noncompact arithmetic H20 or H22+ search to pursue. Borcherds'
[Section 8, Example 5 and the 1998 remark](https://math.berkeley.edu/~reb/papers/al/al.pdf)
describe the rank-22 example and Esselmann's restriction in rank at least
21. Lattice rank is one greater than hyperbolic dimension. The precise
dimension range and the distinct compact/arithmetic upper bound of 29 are
also stated in [Belolipetsky, Section 1](https://arxiv.org/pdf/1506.03111).
The [existing literature note](arithmetic_dimension_bounds_and_effective_high_dimension_search.md)
records the original references and recent frontier checks.

[Allcock's published edge-walking algorithm](https://web.ma.utexas.edu/users/allcock/research/edgewalking3.pdf)
is explicitly a lattice algorithm over the integers. Its edge traversal
and binary quadratic reduction can be more efficient than Vinberg's
distance-ordered wall enumeration, but it is not an unrestricted
nonarithmetic combinatorial generator. The present geometric construction
and Gale search therefore use separate angle, realization, and global
volume checks. The dimension targets remain nonarithmetic finite-volume
H20 and compact H9.

## A matching complete H13 source

The previous E7 face exposed an `A11+A2` corner but no suitable source was
then known in the search. Matching the rational quadratic form to a
reflective lattice did not supply the required corner. This turn finds it
as an E8 face of the known Borcherds H21 chamber, at archived corner orbit
243. Its original spatial corner is `A11(norm2)+A2(norm4)`; normalizing the
A2 walls and the timelike direction adjoins sqrt2. A geometric isometry
need not preserve the chosen rational lattice.

The [complete source](fresh_h13_a11a2_completion_20261014.md) has fifteen
actual facets, 44 finite vertices, and two bounded intrinsic cusp links.
Its finite volume is independently established in the extracted frame:
two positive inequalities bound all thirteen corner-gap coordinates, and
all 105 thirteen-wall active subsets are checked. All vertices are causal
and both cusp quotients have rank 12. The primary extraction also classifies
all 210 old H21 wall projections.

The [independent review](fresh_h13_a11a2_completion_20261014_review.md)
uses another complete enumeration: a vertex of the positive orthant cut
by two inequalities has at most two positive coordinates. It checks the
origin, thirteen axes, and 78 coordinate pairs, reconstructing the same
46 vertices. It also checks all wall angles and the normalized transfer
Gram. Intrinsic completeness does not depend on assuming completeness of
the original H21 projection census.

A separate [integral-overlattice gate](fresh_h13_a11a2_lattice_gate_20261014.md)
explains why the earlier fixed-lattice route was restrictive. The form
`A11+A2+<-2>` has only one proper integral overlattice, an odd index-two
extension; both contain the norm-six mirror bisecting the A2 chamber.
There is no proper even or index-three integral extension. This excludes
neither restricted reflection subgroups nor the successful geometric
transfer over a larger field. Bounded exploratory lattice and diagonal
prefixes are saved without reflectivity or completeness claims.

## Backtracking and closing the two ambient cusp links

All four corner orientations of the source and both necessary E7 lifts
fail against the fixed previous 87-wall chamber. The right source caps
have an explicit acuteness or forbidden-angle obstruction at old wall 26;
the nonright lift has a forbidden normalized cosine of `1/4` with an A11
wall. These sixteen transfer tests are exact.

The [replacement](fresh_h20_e7_intrinsic_caps_20261014.md) removes old walls
`17,26,75,76` and adds the two right source caps. It has 85 actual facets,
a strict interior point, and the complete H13 face. The 44 finite source
vertices are actual ambient H20 vertices. Its two ideal rays initially
have a remaining spherical E7 factor in the ambient link, with an
explicit spacelike recession. Intrinsic cusp boundedness alone is not
sufficient for ambient finite volume.

The [independent 85-wall review](fresh_h20_e7_intrinsic_caps_20261014_review.md)
reconstructs both caps, all 3,570 wall pairs, all facet witnesses, and an
explicit isometry of every source vertex. The 63 other nonincident wall
projections are future timelike, so they impose no extra boundaries on
the complete hyperbolic face.

The [two-cusp extension](fresh_h20_e7_intrinsic_caps_20261014_cusps.md)
then adds `r_i=-omega+2nu_i`, where `nu_i` are the two intrinsic ideal rays
and `omega` is the E7 fundamental weight of norm 2 at old wall 15.
Both added normals have norm 2. Their affine E7 factors complete the
ambient cusp quotient to rank 19. All 171 new pairs pass: 39 right,
two at pi/3, and 130 divergent. The final chamber has 87 facets and
preserves all 44 verified finite vertices on this face. It is a different
87-wall chamber from the predecessor: four old walls were replaced by
two source caps and two cusp fillers.

The [independent cusp review](fresh_h20_e7_intrinsic_caps_20261014_cusps_review.md)
reconstructs the two caps and all 171 new pairs, checks all 87 facet
witnesses and all 44 preserved vertices, and derives the affine blocks
directly from the full physical incident Gram matrices. Their sizes are
`9+5+8` and `13+8`; positive kernel relations and quotient rank 19 certify
both full cusp links. Primary and independent saved replays pass.

The original exact global escape survives. Its squared norm remains

`1336028658157091827 / 250000000000000000 > 0`.

It is strict against every final wall. Thus the new chamber still has
infinite volume. Neither a full census of ambient vertices nor global
closure is claimed. The old completed H12 face and its cusp links are
not assumed to survive the four deletions.

A [fresh four-wall witness](fresh_h20_e7_intrinsic_caps_20261014_cusps_nonquasi.md)
re-establishes the conditional non-quasi-arithmetic obstruction after the
deletions. Current walls 19 and 20 have Cartan cycle `12+8sqrt2`, while
walls 22 and 23 give the rational unit-normal Gram block
`[[1,-3],[-3,1]]`, of determinant `-8`. Any finite-volume completion
retaining these four facets therefore fails the necessary conjugate
positive-semidefiniteness condition. The current infinite-volume chamber
is not claimed to be a lattice.

## An all-parameter obstruction for retaining the old walls

The [new exact gate](fresh_h20_e7_global_escape_gate_20261014.md) proves
that every compatible additional wall which supplies an actual facet on
the selected E7 face retains the predecessor's global escape. It covers
both necessary lift types and every Coxeter angle label and distance.
It is not a finite-height or finite-label search.

For right lifts, three old-wall inequalities give a positive rational
dual certificate. For nonright lifts, two old-wall inequalities and the
Coxeter coefficient gaps force `d_12=0` and `d_0>=sqrt2`. A strictly
convex quadratic then bounds the projected norm below by a value greater
than `1/2`, contradicting the required norm of an escape-cutting lift.
The full finite-corner inverse proves the needed time orientation for
any valid added wall; no separate apex-preservation assumption is made.

The [independent arithmetic review](fresh_h20_e7_global_escape_gate_20261014_review.md)
reconstructs every coefficient and sign using sparse multiplication and
outward rational radical intervals, separately from the producer's
recursive quadratic arithmetic. The gate excludes only the stated face
facet class while retaining all old walls. Other walls, including null
projections used for cusp repair, and replacements remain outside it.
The [separate proof review](fresh_h20_e7_global_escape_gate_20261014_proof_review.md)
checks completeness of the two lift types, the automatic positive time
coefficient, both Coxeter gaps, and the constrained quadratic minimum.

## Compact H9: two different guided proposal mechanisms

The [guided search report](fresh_h9_obstruction_guided_gale_20261014.md)
records 12,000 obstruction-targeted minor-boundary proposals, followed by
2,000 balanced positive-circuit carrier attempts. Together they produce
74 new realizable nerves outside the old 1,063 keys and the standard
single-face-truncation family. Seventy fail literal local obstructions
and four fail the two-Lanner-block gate. No full label assignment reaches
metric solving, and no saved case remains unresolved.

Separate integer incidence replay checks 222,222 five-minors, 370,370
six-circuits, 1,110 links, and 1,226 explicit guidance witnesses. It
verifies applications of inherited exact obstruction proofs; it does
not claim a fresh re-enumeration of those label libraries. A separate
625-recipe 5-to-5 flip diagnostic explains the search bottleneck but
makes no polytopality claim for the unrealized resulting complexes.

The next compact search should change the incidence model rather than
extend the same walk. Nonuniform Gale arrangements or controlled paths
through larger missing faces require fresh exact realization checks.
There is no general fifteen-facet impossibility conclusion here.

## Verification and remaining work

All primary mathematical packages have saved replay commands in their
individual reports. Independent checks cover the intrinsic source,
ambient replacement, both full cusp links, and universal E7 face-facet
gate. The guided Gale
certificate binds every saved case and terminal search state.

The [checkpoint auditor](../scripts/audit_e7_face_and_guided_search_checkpoint_20261014.py)
checks recursive source hashes, nested case records, and local report
links while preserving the preceding checkpoint's frozen mathematical
sources and reports. Its [saved result](e7_face_transfer_and_obstruction_guided_search_2026_10_14_audit.json)
is an integrity certificate, not an additional mathematical proof.

The target remains a new finite-volume example. Further H20 construction
must remove the remaining global escape; completing another face without
a global closure check is insufficient. The current explicit source,
repaired cusp links, and obstruction certificates narrow that next search.
