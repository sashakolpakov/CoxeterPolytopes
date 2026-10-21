# Root written review of the positive ideal-four cusp defect

19 October2026. **PASS for the written theorem V(Q)≥3V(R)+3**, relative
to the explicitly inherited acute-normal, actual-face projection and
entirely ideal box-link facts. No computation was used in this review.
This does not prove an ideal dimension exclusion.

Reviewed complete source:
[positive cusp-defect proof](right_angled_ideal4_cusp_defect_research_20261019.md),
SHA256 `a50bb5ecc0947870c969098445618de8504c9ade119bd0e4706ccd010769f336`.
The discarded stronger+6 attempt plays no part in the reviewed proof.

## Actual facets and the exact defect

With unit normals, n+h_u is a future null vector at the actual cusp u,
orthogonal to n. At an edge uv of R the other two facet normals a,b
are orthogonal to h_u and h_v. If the latter were orthogonal too,
four actual mutually orthogonal walls would force a finite vertex.
Thus h_u.h_v≤−1 and −t_u.t_v≥2.

At equality, the three sides n,h_u,h_v bound an ideal triangle in the
actual a,b section. Its third null vertex h_u+h_v is feasible against
every other wall; positive combinations and acuteness give that directly.
Convexity of the three feasible ideal vertices proves that no additional
wall truncates the triangle. Conversely an outside vertex with two
R-neighbors gives the two endpoints of one cusp-box edge in the
capacity-two proof. Their common vertical walls and n give that same
triangle. The R-edge determines its two-face and third vertex uniquely.
Hence q=N1/2 is an equality, not an assumed maximum multiplicity.

## Spinor normalization and the rank argument

The formula T(a,b)=(a²+b²,a²−b²,2ab) gives the exact Lorentz product
−2det(v,w)². Consecutive determinant magnitudes are therefore at least1.
For ordered projective lifts the wrap signs v0=−vk and v(k+1)=−v1
are consistent. Sending the middle lift to(1,0) gives the displayed
two unit side normals and the coefficient gamma with no missing factor2.
The same opposite pair contributes the same coefficient on its two
polygonal faces, so the four corner coefficients at a cusp are paired.

The recurrence matrix has two independent coordinate kernel vectors.
Deleting any exceptional diagonal leaves a path; signs of its edge
entries do not affect the absolute-value lower bound. Diagonal≥2 and
edge magnitude≤1 give a strictly positive quadratic form, including
the two endpoint terms. This contradicts nullity≥2 if at most one
coefficient is<2. The argument applies to triangles as well: the wrap
edge remains distinct and deletion leaves a two-vertex path.

For a triangle with all three outgoing edge-faces triangular, the three
consecutive determinants equal1 and the three coefficients equal1.
R has at least8 triangles by its degree4/Euler polygon identity. At
most2q triangles meet exceptional edges, so the stated lower bound on
small corners follows. Pairing at every vertex gives
G≥2+max(0,4−q).

## Good opposite facets and empty-neighbor vertices

If h_u were orthogonal to another h_v then t_u.t_v=−1. The two vertices
cannot be adjacent in R. Nonadjacent vertices share at most one face,
so an opposite pair through u misses v. Its two products with t_v
are strictly negative and hence at most−1 by the actual-wall product
gap. Its coefficient at u must then be≥2, contradicting goodness.
Thus a good opposite facet is orthogonal to no other opposite facet.

The ideal three-facet h_u has at least6 vertices; u has its4 outgoing
R-edges inside it. An additional vertex w is not adjacent to u, even
in the ambient polytope: an edge whose endpoints lie on a face lies
in that face. If w had a neighbor v on R, the edge wv would lie on
h_v. At w the two opposite facets h_u,h_v must be parallel faces of
the rectangular box, because they cannot be orthogonal. In upper
half-space R's sphere is tangent to both parallel vertical planes.
The known box corner v then forces the other tangency point u to be
the opposite box corner. The geodesic u-w is feasible by convexity;
the common vertical walls make it an actual edge. This contradicts
the choice of w. Consequently w has no R-neighbor.

Any three distinct facets through such a vertex include an orthogonal
pair in its three-box link. Therefore at most2 of the good h_u can
contain it, proving 2N0≥G. Combining this with the previous bound
gives 2N0+q≥6 and N0+q≥3. In particular N0≥1 is a separate structural
consequence. No octahedral rigidity or invalid sharing bound is used.

## Numerical consequences and limits

For the existing positive-kappa carrier with m≥11, the new inequality
gives V(Q)≥36. The separate injection still gives F(Q)≥25. Ambient
injection yields F5≥62, and strict Euler incidence gives V5≥189.
Another injection gives F6≥252; the unchanged transport endpoint
99f4≤99F+136V and2f4≥62F give68V≥1485F. At F252 the integer
ceiling is5504. The updated ideal-face fractions are<2/23 and<1/15.
These steps use the previously reviewed small-facet bridge and the
published minimum24 and NK/Dufour premises with their stated scope.
They neither construct an ideal five-/six-polytope nor exclude one.

The new conditional Lean successor source formalizes the integer
counting from these individual interfaces. Its compilation, if executed,
must be recorded separately; this written review is not a kernel audit.
