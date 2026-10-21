# Compact dodecahedral/double spectrum: an intrinsic cutting theorem

20 October 2026. New written research proof, awaiting independent review.
No calculation, census, optimization, Lean run or typesetting run was used.

## Scope and result

Let D be the regular compact right-angled hyperbolic dodecahedron, let L
be its edge length, and let D2 be the reflection double of D across one
pentagonal face. Let T be the regular right-angled hyperbolic 120-cell.

**Theorem proposed for review.** Suppose P is a compact right-angled
hyperbolic four-polytope and every intrinsic facet is isometric to D or
D2. Then P admits cuts along pairwise disjoint actual hyperbolic
dodecahedra into finitely many copies of T. If N is its number of D2
facets and m is the number of pieces, then

    N = 12(m-1),       Vol(P) = m Vol(T).

In particular, within this explicitly specified intrinsic facet class,
the minimum is T and the next distinct volume is 2 Vol(T), attained only
by the reflection double of T across a facet. The result does not say
that an arbitrary compact right-angled four-polytope has this facet
spectrum, or that its decomposition graph must be a chain.

This is a structural reduction, not the missing universal minimum or
second-volume theorem. To apply it universally at low volume one must
still exclude every other facet type. The existing cutoffs of 46
vertices at or below Vol(T), and 96 strictly below 2 Vol(T), do not do
that.

## 1. Foundations and precisely imported facts

Use the Lorentz model of H4 with signature (-,+,+,+,+). Actual wall
normals are outward unit spacelike vectors. Distinct actual wall
products are nonpositive. Compact right-angled polytopes are simple:
an edge has exactly three incident facets, a ridge exactly two, and a
vertex exactly four. The orthogonal-projection argument in the current
compact paper proves that mutually orthogonal actual walls have their
expected actual intersection.

The only regular metric objects used below are the given D and T.
Their existence, and the regularity/symmetries inherent in those
objects, are inputs. No four-dimensional classification is used.
All-D recognition is proved by continuation, as in the preserved
zero-excess proof. Standard convex facts used are connectedness of a
compact polytope's boundary and facet adjacency graph, and uniqueness
of a face as the intersection of two given facets.

## 2. The canonical pentagon in D2

Construct D2 as D together with its reflected copy across the gluing
pentagon S. A regular pentagon not meeting S remains a regular
pentagonal face of D2. The five pentagons adjacent to S pair with their
reflections and merge into five right-angled hexagons. There are
twelve remaining pentagons. Thus D2 has 17 faces and 30 vertices.

Each of those five hexagons has cyclic side lengths

    (2L, L, L, 2L, L, L).

Indeed the two sides meeting either endpoint of the removed gluing
edge become a single straight side of length 2L; all other edges
retain length L. The two long sides are opposite. Consequently the
segment joining their midpoints is intrinsically defined. It is exactly
the former gluing edge, of length L. The five such segments form S.

This shows that S does not depend on a chosen presentation of D2 as a
double: the long sides and their midpoints recover every edge of S.
The five vertices of S are the midpoints of the five long edges of
D2. Each long edge is shared by two hexagonal faces. Every other
edge has length L. Reflection in the plane of S interchanges the two
halves and fixes each hexagonal supporting plane; that plane meets S
orthogonally. The pentagon S is regular and right-angled.

## 3. Canonical seams glue without a choice

Join D2 facets of P when they share a hexagonal ridge. A hexagon cannot
border a D facet. In either adjacent D2 facet it has the same two long
sides, so its two descriptions give exactly the same midpoint segment.
The canonical pentagons therefore glue along their actual geometric
edges, not merely abstract identifications.

Consider a vertex of a canonical pentagon in a facet F. It is the
midpoint of an actual edge e of F, hence of P. Let the other two parent
facets containing e be G and H. The faces F intersect G and F intersect H
are the two hexagons at the long edge e. Therefore G and H are D2.
Viewed in G, the same edge e is a long edge of its hexagonal face
G intersect F. A long edge of D2 is incident to two hexagons, so its
other face, G intersect H, is hexagonal too. The same holds in H.
All three canonical pentagon vertices are the midpoint of e.

Thus every seam edge has two incident pentagons and every seam vertex
has three. Locally the seams form a closed two-manifold, and every
hexagonal-adjacency component supplies one connected such surface.
This local assertion uses actual simplicity, not a guessed face graph.

## 4. All seams in a component lie in one hyperplane

For a D2 facet F write n_F for its actual unit normal and let V_F be the
linear Lorentz three-space of its canonical pentagon plane. Its five
hexagonal neighbors have actual normals equal to the intrinsic unit
normals of those hexagonal faces inside the wall of F: right angles
give orthogonality to n_F. Since the hexagonal planes are orthogonal
to the pentagon plane, all five normals lie in V_F. They span V_F,
being the side normals of a compact right-angled pentagon.

Set U_F = span(n_F, V_F). This is a nondegenerate Lorentz four-space.

Suppose F and G share a hexagon. Its canonical seam has two endpoints.
Let H and I be the third facets at the two corresponding long edges.
They are distinct D2 facets by Section 3. In the pentagon of F, their
side normals are the two sides adjacent to the side indexed by G.
They are nonadjacent sides of that pentagon, so their supporting lines
are ultraparallel. Therefore n_H and n_I span a nondegenerate Lorentz
two-plane. Both are orthogonal to n_F and n_G, while n_F and n_G are
orthogonal unit vectors. Hence

    span(n_F,n_G,n_H,n_I)

is a nondegenerate Lorentz four-space. It equals U_F because n_G,n_H,n_I
already span V_F. Interchanging F and G proves that it also equals U_G.
Connectedness makes U_F the same space U throughout the component.

Let H0 = H4 intersect U. It is a genuine hyperbolic three-hyperplane;
its normal line U-perpendicular is spacelike. Every component wall has
normal in U, so H0 is orthogonal to every such wall. Furthermore

    U intersect n_F-perpendicular = V_F,

so H0 intersects F exactly in its canonical pentagon, including the
whole polygon. This is the geometric recognition step that a merely
closed abstract surface would not supply.

## 5. The actual section is a whole dodecahedron

Let Q = P intersect H0. It has nonempty three-dimensional interior.
For example, a relative interior point of one canonical pentagon is
strict on every parent inequality except its own facet inequality;
moving slightly inward within H0 makes that last inequality strict too.
Compactness of P makes Q compact.

Every component pentagon is an actual whole facet of Q, since
P intersect H0 intersect F equals that pentagon. Their union is closed
in the boundary of Q. It is also open there. At a pentagon interior
there is one parent wall; at an edge interior there are exactly its
two parent walls; at a pentagon vertex there are exactly the three
parent walls of the actual edge of P whose midpoint it is. Those
walls all belong to the same component, by Section 3. No outside wall
passes through any of these points. These three cases exhaust the
pentagon surface, and show that it provides a full boundary
neighborhood at each point. Since the boundary of Q is connected,
these pentagons are its entire boundary.

The section Q is right-angled because its facet normals are the actual
orthogonal parent normals. It is a compact simple three-polyhedron
all of whose faces are pentagons. If q is its number of faces, incidence
and Euler give 5q=2e, 3v=2e and v-e+q=2, whence q=12.

It is isometric to D without a further classification input. Match one
regular pentagonal face to a face of D. The adjacent supporting plane
is determined by the common edge and its right angle to the matched
plane. A regular right-angled pentagon with that full edge on the
specified side is unique. Thus the adjacent whole pentagon agrees;
continue across the connected face graph. The finite matched set is
closed under neighbors in both polyhedra, so all faces and oriented
halfspaces agree. In particular the chosen component contains exactly
twelve D2 facets.

## 6. The cut introduces no geometric defects

An outside facet does not meet Q at any point of its boundary: its
relative face, edge and vertex strata have just been identified with
actual parent strata having respectively one, two and three component
facets. It cannot meet the interior of Q either. Its supporting
functional is nonpositive on Q; vanishing at an interior point would
force it to vanish on H0 and hence also on every canonical pentagon,
contradicting the face-interior assertion.

Thus no outside facet meets H0 inside P, and its connected facet lies
strictly on one side. An interior point of Q is strict on every parent
inequality, so it is an interior point of P and both halves have
nonempty interior. H0 cuts each of the twelve D2 facets into its two
original D copies and cuts no other facet. Its new facet Q is D, and
all new dihedral angles are right, because H0 is orthogonal to all
twelve cut wall planes. The two pieces are compact and convex.

The D/D2 spectrum is therefore preserved, and the sum of their D2
facet counts is N-12. Repetition terminates. An outside component stays
on one side of this section: hexagon-adjacent outside facets cannot
lie on opposite sides, since their actual shared ridge would have to
meet the section. Its own canonical section, the convex hull of its
pentagon surface, stays on that side as well. The new and old sections
are consequently disjoint; all earlier cut facets remain uncut later.

## 7. Final recognition, tree, and numerical consequences

A final piece has only D facets. Place one on a facet of T. Across any
whole pentagon the adjacent wall is uniquely determined by
orthogonality. Within that wall, the regular dodecahedron containing
the whole pentagon on the specified side is unique: any symmetry of
the regular pentagon extends to a symmetry of D. Continue across
connected facet graphs in both the piece and T. Their whole facets
and oriented halfspaces agree, so the piece is T. This is the same
geometric continuation as the independently reviewed W=0 result;
there is no inference of whole-polytope equality just from counts.

Each cut removes twelve D2 facets from the total and adds one piece.
Hence N=12(m-1). The geometric sections are disjoint, so the adjacency
graph of the T pieces is a tree. No assertion that it is a path is
made. Volume additivity gives Vol(P)=m Vol(T). When m=2 there is one
common whole dodecahedral facet; placing the two copies on its opposite
sides gives exactly the reflection double. The symmetry group of T is
transitive on its facets, so its isometry type is unique.

For completeness the resulting face data are

    number of D2 facets = 12(m-1),
    number of D facets  = 94m+26,
    total facets       = 106m+14,
    vertices           = 560m+40,
    pentagon excess W  = 30(m-1).

One may obtain the facet counts by reversing the cuts: each gluing
removes two D facets and combines twelve D pairs into D2 facets.
Alternatively sum the 20 or 30 intrinsic vertices over facets and
use four incidences at each parent vertex. The volume formula
Vol(P)=pi^2(f+W+16)/12 then agrees with m*34pi^2/3.

## 8. Why this is not the universal answer

The hypothesis specifies intrinsic metric facet types, not just that
the ridges are pentagons and hexagons. For instance the Loebell L6 has
14 faces and defect two; D2 has 17 faces and defect five. They are
different polyhedra. An L6 facet is not admitted here. Nor does an
arbitrary facet with twelve pentagons and five hexagons automatically
have the stated D2 metric without a separate recognition proof.

The all-wall sign, actual parent strata, Lorentz span propagation and
whole-section argument explicitly prevent the usual recognition gaps.
We do not assume all compact polytopes arise by gluing T copies. The
universal task remains to force the required spectrum from a low
volume assumption, or to find a different argument proving f+W>=120.

## 9. Relation to old and current material

The construction of polytopes by assembling dodecahedra and 120-cells
is published. A primary example is Kolpakov--Martelli--Tschantz,
*Some hyperbolic three-manifolds that bound geometrically*, Section 3:
https://arpi.unipi.it/bitstream/11568/765892.8/5/1311.2993.pdf .
That is a forward construction, not used here as a converse theorem.
The converse above is derived from the intrinsic D/D2 spectrum.
No worldwide novelty claim is made; a full literature comparison is
still pending independent mathematical review of this note.

Kolpakov's thesis Section 4.4 explicitly poses the universal compact
minimum-volume and minimum-facet questions:
https://homeweb.unifr.ch/kellerha/pub/KolpakovThesis_Final.pdf .
The present note does not close those unrestricted questions.

Source interfaces read in this pass (not modified):

* manuscripts/compact_ra4_volume_constraints.tex,
  SHA256 afe807672f3a6bbf4bb2dd95d2861b60391d7c2d7c879c6d6382e749e7ad976f.
* results/compact_ra4_neighbor_excess_and_zero_excess_rigidity_20261019.md,
  SHA256 9c670a7b3c25dda48997c3f70a6d1be9be88bcdefb475b7710f5011068aebcbe.
* manuscripts/ideal_ra4_volume_and_structure.tex,
  SHA256 e503ad2cb8edafd827f735a62568a3372407e9dfdaf3f50d09b05fb3bc75f025.

The ideal paper's opposite-cusp distinctness, canonical section,
all-wall width equations and legitimate-cut proof were reread as an
analogy. No geometric recognition flaw was found in those steps; that
is an independent written audit scope, not a substitute for their
proofs and not a claim of new machine verification. The present compact
argument avoids transferring the ideal cusp-width equations: it uses
the different four-normal span identity proved in Section 4.
