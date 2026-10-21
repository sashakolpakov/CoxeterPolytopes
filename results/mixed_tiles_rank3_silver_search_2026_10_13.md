# Mixed tiles, targeted cuts, and silver lattices — 13 October 2026

**No new compact or finite-volume Coxeter dimension has been established.**
This continues the [polar-cut and gluing search](polar_blocks_h20_search_2026_10_13.md),
targeting compact H9 and finite-volume H20. The strongest positive
geometric result is an exactly certified 767-facet finite-volume H20
polytope; an explicit forbidden ridge proves that it is not Coxeter.
New H9 cones have up to 934 walls, but all saved candidates have exact
spacelike escapes and therefore infinite volume.

Each exclusion below has an explicit scope. None excludes a dimension.

## Golden H9: select directions, then search reflection roots

Freeing the scalar/E8 ratio of the earlier 372-wall cone's escape gives
a rank-three slice with form `<-2phi>+<14>+<1>`. A complete stated box
tests 3,847,150 parameter pairs and finds seven reflective cuts. Their
block-swap images give seven two-wall orbits, with four maximal
compatible selections. All four resulting 382-wall cones have exact
spacelike escapes; every finite angle passes the complete Coxeter gate.

For one surviving escape, two nonnegative projection identities show
that **no compatible real wall in either tested rank-three slice can
cut it**, at any bound. This is a slice obstruction, not full-host
nonreflectivity. See the [rank-three report](fresh_rank3_h9_report_20261013.md)
for the independent Fraction replay and the literal identities.

The next screen tests all 28 rank-four spaces spanned by time, the
scalar coordinate, and two standard E8 fundamental weights. Exact
Farkas certificates exclude 27 of them for the chosen escape. The
remaining pair is `omega4,omega7`, using zero-based indices, with Gram
`[[12,9],[9,8]]`. It contains a geometric cutter of norm `16-9phi`
pairing positively with the escape as `9-11phi/2`. That cutter is not
an integral reflection root. The
[28-slice certificate](fresh_rank4_h9_slice_screen_20261013.json) is
replayed by [this script](../scripts/fresh_rank4_h9_slice_screen_20261013.py).

The binary line through that cutter gives
`(82+107phi)t^2-2phi*x^2=q`, with `q=1,2`.
No integral solution occurs in the complete stated range `0<t<=10^6`.
The [independent bounded replay](fresh_rank3_h9_newplane_audit_20261013.json)
also verifies all 28 slice certificates. This failure concerns that
line, not the full rank-four space.

Freeing all three spatial coefficients succeeds:

`r=(x,t*omega4+u*omega7,z)`.

The box `0<t<=300`, `0<u<=600`, `0<=z<=600`, with exact conjugate bounds,
tests 48,383,940 triples. It finds **335 reflective cuts**, 86 of norm
one and 249 of norm two. Each is compatible with all 382 existing walls
and cuts the previous escape. An independent arithmetic audit replays
all 127,970 cross-pair angle checks, reflection coefficients, strict
interior, bounds, and cut inequalities.

- [Free-ratio search](../scripts/fresh_rank4_h9_free_search_20261013.py)
- [335-wall archive](fresh_rank4_h9_free_search_20261013.json)
- [Independent cut audit](fresh_rank3_h9_free_audit_20261013.json)

Block-swap images give a **670-wall dictionary**, with 340 incompatible
pairs. Twelve compatible selections add 472–552 walls to the original
382, giving cones with **854–934 walls**. Every one has a saved exact
spacelike extreme ray. The
[completion archive](fresh_rank4_h9_completion_20261013.json) records
each selection; the [dictionary](fresh_rank4_h9_dictionary_20261013.json)
records every normal and conflict. These counts are not finite-volume
polytope discoveries.

There is a stronger conclusion than twelve failed selections. One
exact spacelike vector satisfies **all 1,052 old and new inequalities**,
so no subset of this dictionary closes volume. The
[independent audit](fresh_rank4_h9_selection_audit_20261013.md) also
checks all 224,115 possible new-root pairs, all 340 conflicts, and
every saved selection and escape.

More strongly, a surviving escape of norm
`901480/358801-2767458phi/1794005>0` has exact positive-cone projection
identities onto both `S=span(time,omega4,omega7,scalar)` and its
block-swap image `U(S)`. The identities use the original 382 walls,
with eight strictly positive coefficients. Thus **no finite collection
of additional original-compatible real walls in S or U(S), at any
bound, can complete this cone to finite volume**. The
[literal certificate](fresh_rank3_h9_rank4_obstruction_20261013.json)
and [independent replay](../scripts/fresh_rank3_h9_rank4_obstruction.py)
close this pair of slices without further enumeration. They do not
exclude directions outside those slices or a construction that changes
the original walls. The [combined slice report](fresh_rank3_h9_newplane_report_20261013.md)
distinguishes this obstruction from the earlier 27-slice screen, which
used a different escape.

## Silver H9: complete norm alphabets and duplicate groups

Over `Z[sqrt2]`, two equal-norm F4 blocks give hosts
`<-alpha>+F4+F4+<1>`. For `alpha=1+sqrt2`, `sqrt2`, and `2+2sqrt2`,
complete candidate norm alphabets and an exact height-30 search give
92, 40, and 92 walls. All have exact spacelike escapes. The one-block
calibration closes as a known compact H4 simplex.

Exact transports and a discriminant argument identify the first and
third hosts' **full integral reflection groups** with the older silver
K3 carrier's full group. This is not an integral isometry of the
three lattices. The new search includes additional norm classes omitted
from the older norm-two search; it does not establish a new group.
The middle host also lies in an earlier diagonal quadratic space,
though equality of its integral group is not proved.

The old norm-two stabilizer argument cannot simply transfer: its
involution is reflection in a norm-`2+sqrt2` mirror already among the
92 walls. That mirror misses their escape. Full-group reflectivity
remains open. See the [construction report](fresh_silver_f4_h9_report_20261013.md)
and [independent audit](fresh_silver_f4_independent_report_20261013.md).

## H20: mixed tiles and cuts that preserve finite volume

The [mixed construction report](fresh_h20_mixed_20261013.md) certifies
a convex finite-volume H20 polytope assembled from four short and four
long facets of the known Borcherds H21 chamber. Eight full seams,
6,136 nonnegative support certificates, and 28 separating planes
prove the geometry. It has **767 true facets**, including 415
projectively nonrational planes over `Q(sqrt2)`.

An exact strict timelike ridge point proves that two actual facets
meet with squared cosine `1/9`. This rejects the Coxeter condition
without confusing arbitrary Gram pairs with actual ridges.

Four rational planes repair selected cusp factors using affine E7
and D4 alcoves but leave earlier global escapes. Outward angle
rotations produce 24 compatible selections, all with exact escapes.
See the [cusp and rotation report](fresh_h20_cusp_repair_20261013.md).

Combining old and mixed planes with a new rotation gives 1,142 distinct
planes. Exact escape-hitting clauses and pair conflicts exclude every
finite-volume Coxeter subset of this dictionary. Adaptive rotations
add 120 further planes, including inward cuts; an exact SAT obstruction
also excludes the literal 1,262-plane dictionary. These conclusions
do not exclude other mirrors or assemblies.

Replacing `r` by `v=r+t*s`, with `t<0`, while retaining `s`, makes the
old inequality a positive combination of new inequalities. Containment
preserves finite volume, while a common strict timelike point ensures
full dimension. The first cut, `v=3r0-sqrt2*r74`, repairs one angle to
`pi/8`. Twelve successive cuts preserve finite volume but retain an
exact forbidden ridge. A simultaneous block projection changes an
alternating-projection loop to a new quartic-field plane with angle
`pi/5`; it too retains a certified forbidden ridge. The
[independent inward/block audit](fresh_h20_cusp_repair_inward_20261013.md)
checks both containment and actual ridge failure.

## Evidence and limits

Every linked lane includes exact data and replay commands. Numerical
optimization proposes supports or extreme rays; exact arithmetic
checks the certificates used for conclusions. Independent arithmetic
audits and focused regression suites pass. The H20 support, dictionary,
containment, and actual-ridge replays pass.

The free-ratio producer explicitly restricts its accepted command-line
bounds to the audited `300,600,600` box or smaller, preserving its
integer-array overflow guarantees. The six square, sign, polynomial,
and forged-certificate regression tests pass, including rejection of
a corrupted projection identity.

Finite wall lists with escapes do not prove nonreflectivity of the
full arithmetic groups. Finite dictionaries do not exclude additional
planes. The finite-volume H20 constructions fail the angle gate.
No world-first polytope or new-dimensional example is claimed.
