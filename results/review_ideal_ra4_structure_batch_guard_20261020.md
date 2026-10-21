# Adversarial structural audit: ideal RA4 octahedral chains

20 October 2026. **No mathematical gap found in the reviewed structural
chain or its application at the second volume.** This is an internal
source audit, not external refereeing, an execution result, or a formal
verification. Earlier project PASS records were not used as proof. I
read the current argument and reconstructed the geometric implications;
I have previously reviewed related internal material. No primary source
was edited and no checker, solver, TeX or Lean job was run.

Reviewed source: `manuscripts/ideal_ra4_volume_and_structure_publication.tex`,
SHA256 `e0c2ca24cec24ffcd736e201d3c7e06d75b28ef7be8093decd321ae3ac54bda1`
(1102 lines). The focus is the equality-spectrum reduction and Sections
8--16, rather than merely the discrete volume calculation.

## Equality really supplies the metric facet hypothesis

At total facet defect 24, the positive-facet argument does not assume
the desired decomposition. The defect-two exclusion counts positive
triangular neighbors in addition to quadrilateral neighbors and cusp
opposites; these are distinct actual facets. The minimum-four branch
uses the matching of the four quadrilaterals to bound each neighbor's
available positive neighbors. The minimum-at-least-five branch produces
four mutually adjacent facets, and the timelike orthogonal projection
makes their intersection an actual finite vertex. Thus those are
geometric contradictions, not abstract graph prohibitions.

For seven positive facets, the absence of a positive triangular neighbor
of any D9 excludes defect four. The remaining defect-six facet has two
alternating quadrilaterals at each vertex. Its embedded medial graph,
not just its abstract degree sequence, recovers the cuboctahedral face
lattice. The explicit model has nonsquare links at every cusp, giving
the injection contradiction. Consequently exactly eight positive facets
are genuine triangular-face doubles of the regular ideal octahedron;
the other facets are regular ideal octahedra. This precedes use of the
chain theorem and supplies its intrinsic metric hypothesis.

## Closed surfaces and eight-triangle components

I checked the possible cusp identifications in lines 432--637. A
noncubic box can only be `2:1:1`: the proposed `2:2:1` case would require
four quadrilateral sides on a square end, whereas an allowed end has
zero or one. The four nonsquare faces of a `2:1:1` box are exactly the
doubled belt facets, and the ends are ordinary octahedra.

The canonical diagonal in a shared quadrilateral is the same from
both sides: its two endpoints are the nonsquare cusps, which the common
ambient box identifies intrinsically. At such a cusp the four middle
segments make one simple square on the box. Thus the compactified link
is a circle, not two sheets or a branch. At finite diagonal points there
are exactly two facet triangles, since the point is inside an actual
ridge. Triangle interiors lie in distinct facet interiors. These facts
establish the claimed embedded closed two-manifold, including the cusp
points. A cubic cusp contributes no triangle vertex.

The Euler bound does not presume orientability or spherical topology.
The primal-tree/dual-tree proof gives `chi <= 2`; the incidence counts
give `chi = n/4`. The seven distinct facets consisting of a doubled
facet, its three ridge neighbors and its three opposite-cusp facets
belong to this same component. Hence `7 <= n <= 8` and `4 | n`, forcing
eight triangles and six vertices. Parallel edges in the six-cusp graph
would put two different cusps on an opposite asymptotic pair of belt
walls. This is impossible. The octahedral graph and all eight distinct
triangle cycles therefore follow without an embedding-classification
assumption.

## One actual hyperplane, and legitimate disjoint cuts

The equal-width equation in lines 638--733 uses **unit** outward wall
normals on a single common horosphere. Opposite pairs bound the two
equal short directions, so the sums agree with their orientations.
The six cube equations eliminate every mixed sign character. Unit
length and adjacent orthogonality then give the stated orthogonal
Lorentz basis with squares `-1/2, 1/2, 1/2, 1/2`. The future orientation
is fixed by the positive null sums. Its span consequently defines a
genuine hyperbolic hyperplane, not just coplanar combinatorial data.

The six null rays are the six actual parent cusps. Their convex hull
lies in P, while the eight parent inequalities give the reverse
inclusion in H. Hence the entire section is the regular ideal
octahedron S. Its edge interiors are inside parent quadrilateral
ridges, not parent edges. This distinction is essential to the
no-outside-facet proof and is correctly used.

Every finite point of S is covered by the interior/face-interior/
edge-interior cases in lines 734--788. An outside wall cannot vanish
in any case. Cutting therefore leaves every outside facet whole,
splits each chosen double along its genuine gluing triangle, and
introduces no finite vertex: a new vertex on H must be extreme in S.
The section center is strict in every parent inequality, so both
pieces are four-dimensional; their finite volume follows by containment.
New ridge adjacencies are precisely the eight orthogonal cut-wall
adjacencies. Outside facets may still meet S at old ideal cusps;
the phrase that S "meets only" the cut walls is read in this
ridge-adjacency context, consistently with the explicit finite-point
statement.

The disjointness argument survives ideal boundary scrutiny. Another
connected belt component is on one side of H, as are the convex hulls
of its cusps. There is no extra null point in `H intersect closure(P)`
outside the six vertices of S: joining such a point to a finite interior
point of S would give a segment in `P intersect H`, forcing the null
point into `closure(S)`. The latter has only its six listed ideal
points. At each of those cusps there is only one belt component.
Different sections thus have disjoint closures, not just disjoint
finite interiors. Each cut removes exactly eight doubled facets in
the total count over its two pieces, so iteration terminates.

## Recognition and the exclusion of branching

The all-octahedral recognition in lines 789--832 propagates **whole**
facets. A shared triangular face and orthogonality fix the neighboring
wall. A regular ideal octahedron on a specified side of that triangle
is unique: its face symmetries extend, and fixing the three ideal
vertices pointwise leaves only identity or reflection in their plane.
Propagation over the model's connected facet graph gives both the
supporting-halfspace containment and the reverse convex-hull
containment. No arbitrary four-polytope classification is needed.

Disjoint separating cuts produce a tree; the tree alone would not
imply a path. Lines 833--912 supply the additional argument. Two
nonopposite disjoint cut facets in a cell have spatial root product
`-1`. Their sum is a D4 root giving a cell facet orthogonal to both.
This facet and its two reflected neighbors are three full octahedra
in one original supporting facet. None can be an internal cut, since
it meets an existing cut along a triangle and all cut closures are
disjoint. Their relative interiors cannot overlap: the common inward
normal collar would make the four-dimensional cell interiors overlap.
Their total intrinsic volume exceeds both permitted original facet
volumes. Thus incident cuts must be opposite, so degree is at most
two and the tree is a path. The cell's flag transitivity extends all
facet symmetries, justifying the reflection placement and removing a
gluing-twist ambiguity.

At second volume the independently established eight doubled facets
give exactly one cut and two cells. The chain theorem therefore
identifies the equality model without assuming universal volume
quantization or universal decomposition. The attaining double's
finite exact incidence outputs remain recorded model checks; I did
not rerun them. Their role does not replace any of the preceding
universal geometric steps.

## External boundary and limits

The imported published facts are Kolpakov's volume/incidence formulas,
the antiprism bound with equality, the admissible nine-vertex entry
and the ten-vertex triangles/quadrilaterals entry with its matching,
and three-dimensional ideal right-angled rigidity. I checked that
the cited small-facet scope is the one stated in Kolpakov, Section
4.3, Case D, Figures 6--7 (printed pp. 1951--1952), rather than an
unrestricted four-dimensional classification:
[published primary paper](https://msp.org/agt/2012/12-4/agt-v12-n4-p03-s.pdf).
Their historical classification proofs are not rederived here.
The acute-wall sign property and elementary convex/topological facts
are also used as stated. The fixed model counts and coordinate
checks are not newly executed results of this audit.

This finding is limited to the manuscript's stated theorems. In
particular it gives no decomposition or divisibility statement for
arbitrary ideal right-angled four-polytopes without the Oct/D facet
hypothesis, and does not assert existence at every chain length.
