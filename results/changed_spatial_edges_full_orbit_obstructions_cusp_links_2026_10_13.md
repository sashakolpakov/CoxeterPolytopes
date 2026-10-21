# D9 nonreflectivity, complete symmetry tests, and new cusp repairs

**No new compact or finite-volume Coxeter example has been established.**
Complete edge paths now prove the quartic D9 lattice globally nonreflective.
Geometric constructions reach 102 facets in H9 and complete a full orbit of
241,920 new cusp links in H22. Exact escapes still prove infinite volume in
every chamber constructed in this round. Two whole cap banks and one
continuous perpendicular-wall family are also excluded with exact proofs.

| Search | Verified advance | Remaining limitation |
|---|---|---|
| Quartic arithmetic H9, D9 lattice | Eleven complete edges and an exact infinite-order symmetry of the full chamber | The entire lattice is nonreflective; further mirror enumeration here is unnecessary |
| Geometric H9 | Six transported caps cut the old escape and create three finite vertices | One escape survives the entire six-cap symmetry bank, including every compatible subset or replacement |
| Geometric H20 | Complete all-label obstruction on a prescribed positive plane; whole polar-reflected bank excluded | The 176-facet chamber remains infinite-volume; unrestricted completion is not excluded |
| Geometric H22 | Three compatible whole-orbit placements, each with 5,397,548 actual facets and 241,920 repaired links | Even the placement that cuts the prior remote escape has a new exact escape |

This supersedes the geometric stages in the
[preceding checkpoint](a9_nonreflectivity_compact_prism_face_orbits_2026_10_13.md).
Its full-host exclusions for the earlier quartic F4 and A9 lattices remain in
force. A different integral lattice does not by itself prove a different
reflection group or commensurability class.

## D9: complete local edge walking gives a global obstruction

The new host is

```text
K=Q(T),  T=2cos(pi/8),  T^4-4T^2+2=0,
L= <-(T-1)> + D9,
D9={y in O_K^9 : sum(y_i) in 2O_K}.
```

Its full zero-height reflection chamber is B9, because individual coordinate
sign changes preserve D9. Every possible primitive reflective norm is included
in a thirteen-class pool dividing eight. This differs from the smaller norm
pool for the earlier unimodular lattices.

Nine complete steps reach a second B9 corner through corners of type
`Bk+B(9-k)`. A separate initial A8 diagonal edge reaches another B9 corner.
An eleventh step on the initial A2+B6 edge reaches B3+B6. The complete censuses
examine 212 finite profiles and 1,364 relative-norm orbit points. They certify
the first mirror along each selected edge, including simultaneous contacts
and the full endpoint stabilizer. The resulting twenty walls pass all 190
pair tests and reach twelve complete finite corners.

The search reduces each selected edge to relative norm equations and enumerates
the full relevant unit orbits over its finite interval. Certified PARI ideal
class and fundamental-unit computations are an explicit dependency. A separate
Fraction implementation repeats the profiles, orbit normalization and interval
censuses; it does not independently implement PARI's number-field algorithms.

The new edge uses a different relative extension, `K(sqrt(3(T-1)))`, with its
own certified full unit and absolute-norm data. Its complete census has 28
profiles and twenty interval points, including two simultaneous roots. Exact
cone relations select the norm-four facet defining the full B3+B6 corner.
It cuts the preceding finite chamber's escape, although another strict escape
of square `2537/1250-T>0` survives the twenty walls.

The stronger result compares this corner with stage three of the chain. An
exact form isometry maps all nine oriented primitive incident normals and
the future corner line from one to the other. Both it and its inverse preserve
the actual D9 module. Its trace is

```text
-7378736 - 3993344*T + 12596300*T^2 + 6817068*T^3 > 10.
```

Four of the complete edges, covering 83 profiles and 74 interval points,
connect these corners to the same initial chamber.
The transport therefore preserves the full reflection chamber, not just a
finite wall set. The trace proves infinite order; discreteness then forces
infinite chamber volume. This closes the D9 host, including every finite-volume
reflection subgroup using its mirrors. There is no need to exhaust its
remaining edge queue.

The [new edge and global conclusion](fresh_h9_changed_spatial_next_20261013.md),
[global theorem](fresh_h9_d9_symmetry_transfer_20261013.md),
[exact transport certificate](fresh_h9_d9_symmetry_transfer_20261013.json), and
[root source review and fresh replay](fresh_h9_d9_symmetry_transfer_20261013_root_review.md)
record the premises. The earlier [ten-edge report](fresh_h9_changed_spatial_20261013.md)
and [independent review](fresh_h9_changed_spatial_20261013_review.md) remain
supporting records. Coordinate boosts that failed D9 integrality were not used
as chamber symmetries. No claim is made about every lattice in this rational
quadratic space or all compact H9 polytopes.

## H9: a useful construction and a full finite-bank exclusion

Transporting six old caps by an exact seed isometry extends the 96-facet chamber
to 102 facets. The transported walls cut the preceding escape and create a
separate component with three finite vertices, two finite edges and 23 open
rays. Every new pair passes, and every possible new positive-semidefinite
subdiagram is positive definite. A new strict spacelike escape remains.

The complete weighted sixteen-seed symmetry group has order 2,304. Its six
source caps give 13,824 literal images, before compatibility screening or
deduplication. One exact strict spacelike point survives all these images and
all 102 old walls. A common strict timelike center is also verified. Hence any
selection from this entire dictionary still has infinite volume, even after
arbitrary old-wall deletions or replacements within the dictionary.

This closes a specific symmetry construction, including the replacement
candidates found by the numerical scout. It does not exclude different source
normals or arbitrary H9 completions. The
[construction and full-bank proof](fresh_h9_prism_completion_20261013.md),
[independent source review](fresh_h9_prism_completion_20261013_source_review.md),
and saved-result replays pass.

## H20: two different construction routes excluded

The last facet of the existing 176-facet chamber has eighteen finite neighbors.
All proposed normals perpendicular both to it and to those neighbors lie in a
positive two-plane. Projection of an old wall bounds its allowed angle labels
by seven possibilities, `2,...,8`. Both directions for each label are tested
rigorously; all fourteen branches fail a retained-wall condition. This is an
all-real-direction, all-integer-label result on that prescribed plane.

A second proposal reflects all 175 predecessor normals across the old gap's
polar hyperplane. Exact pairing signs show that the current escape survives
every reflected normal. Therefore every selection from the entire reflected
bank and existing walls remains infinite-volume, including deletions of the
last facet. The argument covers the whole bank without relying on the numerical
compatibility screen.

Both obstructions are scoped. They leave nonright neighbor relations, different
polar hyperplanes and more general replacements open. The
[proof and search scopes](fresh_h20_facet_bridge_20261013.md) and
[independent source review and replay](fresh_h20_facet_bridge_20261013_review.md)
record the exact bounds and dependencies.

## H22: all specified new cusp links repaired, but volume still infinite

At angle label six, the preceding cap orbit has 241,920 signed deficient cusps.
Each link contains seven spherical components and an affine A1, leaving two
missing tangent directions. Seven negative highest-root caps turn the
spherical components affine; four more caps bound the two remaining directions.
Their orthogonal components have type

```text
affine A5 + affine A4 + 2 affine A3 + 6 affine A1.
```

Positive affine relations and total quotient rank 21 prove these links bounded.
Exact equivariant tangent frames transport the repair around the complete
120,960-element old-coordinate cusp orbit and both height signs. All pair
relations pass, and acuteness with a strict center proves that the resulting
5,397,548 inequalities are distinct actual facets.

Three placements are separately certified. A closer placement creates new
permitted finite angles with original walls and cuts the remote escape from
the first two placements. Another exact rational strict escape then remains,
with square `3511216341/10^10`. Thus every constructed chamber still has infinite
volume, uniformly in the old center label `m>=3`. These repairs fix the new
label at `k=6`; they do not extend the previous all-`k>=6` theorem to a completed
finite-volume family.

The closer singleton caps also introduce additional parallel pairs. Complete
incidence checks at two resulting null rays find only eighteen incident walls,
quotient rank seventeen and four independent unbounded directions in each
21-dimensional horosphere. The
[exact cusp diagnosis and root review](fresh_h22_lower_new_cusp_20261013_review.md)
therefore identify a new local defect as well as the global escape. Completing
the named cusp orbit has not completed every new end. The
[construction report](fresh_h22_new_cusp_completion_20261013.md) gives every
placement and exact escape. The
[independent review](fresh_h22_new_cusp_completion_20261013_review.md) includes
fresh large-orbit replays and separate Fraction checks of all local link data,
old-wall relations, and the target cut.

A further [exact all-real-shift census](fresh_h22_lower_incident_family_20261013.md)
excludes finite repairs from the same eleven transported affine-ray families
while the lowered cap `c_4` and all original 210 facets are retained. Any new
wall incident at its exposed cusp must satisfy two linear equations. The full
signed orbit gives 1,126,349 literal rational solutions, with all degeneracies and
orientations handled. After separating redundant old-wall copies, five literal
solutions representing three new normals pass the incident-wall angle gate.
Every one has an exact positive pairing with another retained original facet.
No additional link constraint is available from these families at any real
shift. A repair therefore needs different normal families or changed retained
facets; this is not unrestricted H22 nonexistence.
The [independent family review](fresh_h22_lower_incident_family_20261013_review.md)
and both the reviewer's and root task's fresh exact replays pass.

## Consequences for the algorithmic search

The [Vinberg–Allcock audit](vinberg_allcock_search_strategy_20261013.md) remains
the guide: Vinberg enumerates walls by distance; Allcock traverses edges and
uses lattice arithmetic to certify the next corner. A combinatorial adaptation
must retain metric variables and prove compatibility and global closure.
The arithmetic work here has complete selected edge steps over a number field;
it does not claim that Allcock's published integer-lattice theorem already
supplies that extension.

The arithmetic search must now change the closed D9 host. The geometric search
needs new H9 normals outside the excluded six-source bank, H20
relations outside the prescribed perpendicular plane, and the remaining ends
of the repaired H22 chambers. Larger counts within the excluded dictionaries
cannot change their exact infinite-volume conclusions.

A concrete arithmetic candidate for the next continuation is
`< -3(T-1) > + D9`. It remains admissible at the four real embeddings. Its
determinant differs from the just-closed host by the nonsquare 3: the existing
[field certificate](fresh_h9_glued_host_20261013_structure.json) proves that 3
is inert, so its valuation is odd. Thus this changes the determinant
squareclass, rather than merely choosing another integral lattice in the same
diagonal form. Its coordinate-axis norm field is `K(sqrt(3(T-1)))`, just
certified above, while its initial A2+B6 edge uses
`K(sqrt(9(T-1)))=K(sqrt(T-1))`, also already certified. The new lattice's full
primitive norm pool, including the prime above 3, must be established before
any first-edge claim. This is an untested next host, not an additional result
of this checkpoint.
