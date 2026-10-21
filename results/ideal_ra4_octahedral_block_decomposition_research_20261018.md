# Decomposition with octahedral and doubled-octahedral facets

18 October2026. New direct written proof, with separate independent
mathematical reviews of the decomposition and its path-rigidity
strengthening. Manuscript compilation is pending. This is a RESTRICTED theorem. It does
not assert that every entirely ideal right-angled H4 polytope has the
facet spectrum assumed below. No four-dimensional census is used.

**Internally checked theorem.** Suppose every facet of a finite-volume entirely
ideal right-angled hyperbolic4-polytope P is a regular ideal octahedron
O or the double D of O along a triangular face. Then P can be cut into
regular ideal24-cells by disjoint regular ideal octahedral sections.
If P has N doubled facets, N=8k for an integer k>=0, the decomposition
has k+1 cells, and Vol(P)=(k+1)*4pi^2/3. The adjacency graph of the cells
across the cuts is a path, with successive cuts opposite in each middle
cell. An existing polytope in this class is determined up to isometry
by its number of cells. No every-length existence assertion is needed.

The proof uses standard Lorentz wall geometry and the explicit O and D
models. It does not use a classification of arbitrary4-polytopes or a
preexisting universal24-cell decomposition assertion. The existing
second-volume argument supplies useful already-written local calculations,
which are repeated as needed below independently of the defect24 budget.

## Cusp boxes and the internal triangle

D has three nonsquare cusp rectangles, each a2-by-1 rectangle formed
by joining two unit squares. At these cusps its two quadrilateral sides
are the two long sides. Its six other cusp rectangles are square and
have exactly one incident quadrilateral side. O has only square cusp
rectangles and no quadrilateral faces. D has three quadrilateral faces;
each pair meets at one of its nonsquare cusps. The convex hull of its
three nonsquare cusps is its internal gluing triangle T. Each side of
T is the interior diagonal of one quadrilateral face.

All face rectangles of an ambient cusp box have aspect1 or2. Therefore
a noncubic box has lengths2:1:1 or2:2:1. The latter is impossible:
the four nonsquare D faces would require both long sides quadrilateral,
so each square end face would have four quadrilateral sides at its
cusp. O has none and D at a square cusp has exactly one. Thus every
noncubic box has lengths2:1:1. Its four nonsquare D faces form the belt,
every belt adjacency is a quadrilateral ridge, and the two ends are O:
all four end incidences are triangular, excluding a square cusp of D.

The triangle T meets a2-by-1 cusp rectangle along its transverse middle
segment. The four such segments in an ambient belt meet at the centers
of their shared long edges and form one simple closed square on the
box surface. This explicitly describes the local link at the ideal
point after compactification. There is no dangling end or branching.

## A finite closed surface for each component

Join two D facets when they share a quadrilateral. Glue their internal
triangles along the common internal diagonal. The facet interiors are
disjoint; a quadrilateral lies in exactly two facets; and its diagonal
is in its relative interior. Thus triangles intersect along precisely
these edges or their actual ideal endpoints. At finite points the result
is a2-manifold. At each ideal point, the four-segment belt circle just
described is its link, so adjoining that ideal point makes it a disk.
We obtain a finite closed2-manifold in the Klein compactification of
the boundary of P. Its connected components correspond exactly to the
components of the quadrilateral adjacency graph: all four triangles at
an ideal point already belong to the same belt cycle in that graph.

In one component let n,e,v count its triangular faces, edges and vertices.
Each triangle has three sides, every edge lies in two triangles, and
every vertex lies in four triangles. Hence

    e=3n/2, v=3n/4, chi=v-e+n=n/4.                 (1)

For completeness, chi<=2 can be proved without assuming the component
is a sphere or importing a classification of surfaces. Choose a
spanning tree A in its one-skeleton. The dual graph remains connected
after deleting the edges dual to A. Otherwise the boundary between a
nonempty proper set of triangles and its complement is a nonempty
Eulerian subgraph contained in A: the circular link of every vertex
makes the number of transitions even. A forest has no nonempty
Eulerian subgraph, a contradiction. Choose a dual spanning tree using
the remaining edges. The two disjoint trees use v-1 and n-1 edges,
respectively, so e>=v+n-2 and chi<=2. Equation(1) gives n<=8 and n
divisible by4.

There are at least seven triangles in the component. For any D facet F,
its three quadrilateral neighbors and the three facets opposite F at
its nonsquare cusps are six distinct D facets. Ridge neighbors differ
because intersections of convex facets are faces. The cusp opposites
differ from every ridge neighbor and from one another because their
walls are asymptotic to F's wall at exactly the specified cusp. All
six belong to the same component: cusp opposites lie across the same
four-facet quadrilateral belt. Consequently7<=n<=8 and n is divisible
by4, so n=8. Each component has six ideal vertices and twelve edges.

## Reconstructing its hyperplane

For this component make a graph J with the six ideal vertices as nodes
and its twelve edges as edges. It is4-regular. It is simple: parallel
edges would come from two distinct quadrilateral ridges with the same
two nonsquare cusps. Their incident D facets include at least three
distinct belt facets. Two of these are opposite at one cusp, so their
walls are asymptotic at that unique ideal point and cannot contain the
second cusp. Therefore J=K6 minus a perfect matching, the octahedral
graph. The three vertices of a D triangle give a triangle in J.
Distinct D facets give distinct triangles because the triangle interior
is in the relative interior of that facet. J has eight triangles, so
all occur. Label the three opposite pairs of J by1,2,3 and signs.
The eight D facets are indexed by s in{+1,-1}^3; quadrilateral
adjacency changes exactly one sign.

Let n_s be their outward unit normals in R^{4,1}. At a belt cusp the
two opposite unit-normal sums are equal because the two short widths
of its2:1:1 box are equal. In upper-half-space, a vertical wall x_i=a
has unit normal(a,e_i,a), and its opposite x_i=b has unit normal
(-b,-e_i,-b); the sum is(a-b)xi for the common cusp null vector xi.
This verifies the width scaling and sign directly. The six belt
relations are precisely the cube-face parallelogram equations. Thus

    n_s=a+s1*b1+s2*b2+s3*b3.

Expanding all eight unit equations kills the nonconstant sign characters:
a.bi=0,bi.bj=0(i!=j),a^2+sum bi^2=1. Right angles along cube edges give
1-2bi^2=0. Hence a^2=-1/2,bi^2=1/2. These independent vectors span a
Lorentz4-space U, whose spacelike orthogonal line defines a hyperplane
H in H4. All eight walls are orthogonal to H.

Their halfspaces in H cut out the regular ideal octahedron
S={|x1|+|x2|+|x3|<=t}. Its six null rays are a+/-bi, exactly the six
actual cusps of this component by the positive opposite-pair sums.
Convexity gives S subset P, while the eight inequalities give
P intersect H subset S. Therefore P intersect H=S. Each face of S is
the internal triangle of the corresponding D facet; each edge interior
is in the interior of a quadrilateral ridge of P.

## The cut and induction

Any facet outside these eight cannot meet S at a finite point. A
supporting functional vanishing in the relative interior of a face of
S vanishes on that face. A triangular face interior belongs to just its
own D facet, and an edge interior belongs to just its two incident D
facets; neither can belong to another facet. The same supporting-face
argument excludes the interior of S. Its vertices are ideal. Thus any
outside facet meeting both sides of H would give a finite intersection
with S, a contradiction. This applies equally to outside O and outside D.

The center of S is strict on every original supporting inequality and
is an interior point of P. Both halves therefore have interior. The
eight cut D facets split into regular octahedra, and the new section
facet is S. Every other facet is uncut. In Klein coordinates a vertex
of a half on H must be a vertex of S, by extremality; a vertex away
from H is an old vertex. Thus no finite vertex is introduced. New
dihedral angles involve S and the eight D walls and are right angles;
all other angles are inherited. Both halves again have only O or D
facets and together contain exactly eight fewer D facets.

Repeat. After k=N/8 cuts there are k+1 pieces, each with only O facets.
Here is also an explicit disjointness check. No outside D facet crosses
the first H. Two quadrilateral-adjacent outside D facets cannot lie on
opposite sides: their common quadrilateral would then lie in H and give
a forbidden finite intersection with S. Consequently every other D
component lies on one side, and its octahedral section, the convex hull
of its six cusps, lies on that side as well. It has no ideal vertex on
the first section: a noncubic cusp belongs to one belt component and
its two ends are O. A finite boundary contact is likewise excluded by
the facet-interior/ridge-interior argument. Thus the sections are
pairwise disjoint, including their ideal endpoints. Each new cut splits
one piece into two; the adjacency graph of terminal pieces is a tree.

## Direct recognition of the all-octahedral pieces

One may use Kolpakov's unique-minimum theorem here, but there is also
a direct propagation argument. Place one regular octahedral facet F
of a piece on a facet of a model regular24-cell C. A regular ideal
octahedron on a given side of one of its triangular face planes is
determined by that face. Indeed an isometry between two such regular
octahedra can have its restriction to the face canceled by a facet
stabilizer symmetry; an ambient H3 isometry fixing the three face
vertices is the identity on that side (the other possibility reflects
the face plane). Therefore every facet of the piece adjacent to F
coincides as a whole polyhedron with the corresponding facet of C:
right angles first determine its supporting wall from the intrinsic
triangle of F, then the shared triangle determines the octahedron.

Apply the same argument to each newly matched facet. The facet adjacency
graph of C is connected, so all24 model facets, with their full ideal
vertices, are actual facets of the piece. Their halfspaces imply the
piece is contained in C. Conversely its closed convex hull contains
all vertices of C and hence C. Thus the piece equals C. This uses
regular-octahedron face symmetry and existence of the regular24-cell,
not a census or classification of4-polytopes.

Volume additivity now gives(k+1)*4pi^2/3. Equivalently the defect identity
gives W=3N/2=12k and Vol=4pi^2/3+W*pi^2/9. These agree.

## What is not proved

The Oct/D facet assumption is essential. The argument does not exclude
other facets at larger volume, prove that every ideal RA4 polytope is
arithmetic, or establish a universal volume divisibility theorem. In
the stated class the proved argument above gives a finite tree of
24-cells. The normal and convexity arguments remain
written mathematics until their geometric interfaces are formalized.

## The path-rigidity strengthening

In one24-cell block use outward unit normals(1,r), with r a D4 root:
exactly two spatial entries are+/-1 and the others zero. Distinct
cut facets have disjoint closures by the section theorem, so their
root dot product is negative: either-1 or-2. If r1.r2=-1, then r1+r2
is again a D4 root, and its facet C is orthogonal to both cut walls.
C is not itself a cut facet because adjacent cuts are impossible.
The facet of the original P containing C would therefore contain
three regular ideal octahedra: C and its reflected neighbors across
the two cuts. Their relative interiors are disjoint. This exceeds
the volume of either O or D, contradicting the facet hypothesis.
Thus any two cut facets in a block are opposite, with r2=-r1.
There cannot be three pairwise opposite facets. The dual tree has
maximum degree two and is therefore a path. Moreover the next cut
in each internal cell is uniquely opposite the previous one. Full
facet symmetry extends to the24-cell, so gluing twists do not change
the isometry type. This determines the polytope
uniquely by its number of cells within this restricted class. This
last paragraph does not assert existence for every path length.

Separate mathematical reviews:
`review_ideal_ra4_octahedral_block_decomposition_mathematical_20261018.md`
and `review_ideal_ra4_octahedral_block_path_rigidity_20261018.md`.
The complete illustrated manuscript is
`manuscripts/ideal_ra4_octahedral_decomposition.tex`. Neither written
review is an external human review, certificate execution, or full
Lean verification of the geometry.
