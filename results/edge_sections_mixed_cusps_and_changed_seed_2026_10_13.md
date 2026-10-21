# Edge sections, mixed cusps, and a changed compact seed

**No new compact or finite-volume Coxeter example, and no new dimension,
is established.** This round adds continuous geometric exclusions and a
different exact H9 starting chamber. The previously certified 29-wall H9
chamber remains infinite-volume; its saved 66-vertex/455-edge census is
unchanged.

## Dimension and algorithm restrictions remain in force

Noncompact arithmetic examples have dimension at most 19 or exactly 21.
The 21 bound does not apply to compact arithmetic examples; the general
arithmetic and compact upper bound is 29. Allcock's published algorithm
works with integral lattices. The geometric edge organization used here
does not make an unrestricted nonarithmetic search an implementation of
that lattice algorithm. See the [literature audit](literature_frontier_recheck_20261013.md)
for the primary references and exact scopes.

## All sections crossing the all-right edges

The [new exact edge theorem](fresh_h20_allright_edge_gate_20261013.md)
excludes every real H20 section crossing any of the 42 open edges incident
to the source's two all-right vertices. It covers generic transverse
sections, not only hyperplanes through old vertices or cusps. An actual
second vertex has a projected D4 block forcing a forbidden angle.

The standard-library checker verifies 42 first edges, 420 auxiliary
cusp-to-cusp edges, and 8,862 incidence-avoidance cases. A separate
Fraction reconstruction imports no producer or repository geometry helper
and passes. The proof's facet-survival and transversality steps were also
reviewed separately.

The argument supplies a useful broader reduction: at any new finite
section vertex lying in an old edge, a Coxeter section normal can pair
nontrivially with at most three of the edge's 20 normals. A three-entry
support must form an old two-edge star, with a cancellation. Thus each
possible support places the section normal in at most five linear
dimensions. This leaves a concrete construction problem at the other
source-edge orbits; it does not exclude them.

## Mixed cusp factors are excluded for all real coefficients

The [mixed-factor theorem](fresh_h20_retained_link_20261013.md) excludes
every section through any of the 18 cusp representatives whose tangent
normal is active in at least two irreducible affine factors. Arbitrary
real coefficients and arbitrary displacement along the null ray are
included. The proof first handles actual facets and redundancy, then
uses Euclidean simplex geometry and exact inverse-metric bounds.

The finite data consist of 53 affine factors, 411 positive-definite
one-wall deletions, 247 connected deletions, and 2,067 attachment profiles.
The producer's saved replay and a
[separate exact inverse reconstruction](fresh_h20_retained_link_20261013_review.json)
pass. Single-factor cusp directions remain
open, apart from the particular planes excluded in preceding work.

## All-real subdivision obstruction in the 59 twin-leaf planes

The [new subdivision theorem](fresh_h20_irrational_subdivision_20261013.md)
excludes every finite face-to-face subdivision of every full-dimensional
section in the 59 previously tested planes into Coxeter polytopes.
Irrational parameters and arbitrary angle labels are included. This is
stronger than the preceding forbidden-angle test for the section itself.

In 57 planes an actual spherical link triangle cannot be tiled by
spherical Coxeter triangles: side-length incompatibility and a noninteger
area ratio rule out the only possible corner angles. The other two planes
have two actual ridges whose putative rational angles force a nonzero
algebraic integer to have all conjugates in [0,1), which is impossible.
The exact certificate checks 118 orientations, including 114 links and
eight exceptional ridge certificates, using 122 face interpolations.
The saved producer replay and a
[separate Fraction geometry reconstruction](fresh_h20_irrational_subdivision_20261013_review.json)
pass, and both mathematical arguments were reviewed separately. Other
planes and constructions that are not finite
face-to-face subdivisions remain open.

## Constructive generic-edge sections with H4 vertex links

A [different bounded construction](fresh_h20_h4_edge_sections_20261013.md)
uses the sparse-pairing reduction at edges with an A4 component. Changing
one endpoint pairing replaces this local A4 diagram by H4. The section
crosses an ordinary open edge, so it has an actual finite Coxeter vertex
and inherits finite volume from the ambient chamber.

The pilot samples nine distinct endpoint-orbit pairs, both A4 endpoints,
three interior positions, and both radical signs: 108 numerical sections.
Each has a numerical defect at an actual finite vertex reached in one
section-edge step. Nine midpoint sections were then reconstructed exactly
over Q(sqrt(5)) with one positive radical. All nine have a separately
certified actual obtuse ridge, checked against all 210 inequalities. Fresh
saved replay and [source/mathematical review](fresh_h20_h4_edge_sections_20261013_review.md)
pass. This rejects those nine exact cases, not all positions or all such
edge sections. No global Coxeter candidate survived this bounded pilot.

## A different compact H5 face for the H9 construction

The [changed compact-face seed](fresh_h9_compactification_20261013.md)
glues the two arithmetic Makarov prisms P0 and P2 along their common
orthogonal simplex facet. The resulting seven-facet compact H5 face is
non-quasi-arithmetic, verified directly by its field and Gram data.
Adding an orthogonal A4 chamber gives an exact eleven-facet H9 Coxeter
chamber with ten finite vertices. This uses different face arithmetic
and transverse geometry from the preceding P0-plus-H4 construction.
The mixed H5 face is not claimed as a new published example.

The first exact cap family has 810 profiles: 60 have nonpositive old
scale, 438 have a forbidden old-wall angle, and 312 are individually
compatible with the seed. All 312 share a strict exact spacelike escape,
so no subset completes finite volume. Their mutual compatibility is not
assumed. A wider 24,010-proposal numerical scout also retains a gap;
the exact exclusion does not extend to that wider bank. Saved exact
replay and a separate source/proof review pass.

## Continuing from this checkpoint

For H20, use other source-edge orbits with the sparse-pairing reduction,
single-factor cusp directions outside excluded planes, or different
ambient assemblies. For compact H9, use different centers or contacts
on the new mixed face, or split incidences in the remaining 15/16-facet
noncompact sources before solving a compact Gram realization. Repeating
the excluded finite cap families cannot succeed.

The [preceding checkpoint](complete_frontier_polar_cover_and_cusp_sections_2026_10_13.md)
retains the exact 29-wall construction, complete frontier, polar-cover
exclusion, and earlier section proofs. Local Coxeter geometry and exact
facet counts remain intermediate results until global finite volume is
proved.
