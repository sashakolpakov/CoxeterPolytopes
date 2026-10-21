# Adversarial read of the second-volume equality argument

18 October 2026. Scope: a fresh mathematical read of the complete equality
report and the inherited gap argument, directed specifically at possible
failure of opposite-cusp injection, N=7 exclusion, unit-normal width
identities, and the actual entirely ideal right-angled cut. No proof or
geometry script was executed, no Lean job was run, and no frozen source
was changed. Review counts are not a premise of any argument below.

**Outcome: no actual flaw found in these steps.** This is a statement of
what this adversarial read found, not a certificate that no flaw exists.
The external small-facet classification and geometric rigidity remain
explicit dependencies; they are not internally reclassified here.

## Opposite-cusp injection is genuinely injective

For a cusp x, take x to infinity in upper-half-space. Its incident facet
walls are vertical Euclidean hyperplanes and its link is a rectangular
box. The facet F and the facet F_x opposite its link face are parallel
vertical walls x_1=b and x_1=a, with a>b. Their closures have exactly one
common ideal point, infinity. Thus the same F_x cannot be opposite F at
another cusp, and cannot share a ridge with F. Different nontriangular
ridge neighbors are also distinct because the intersection of two faces
of a convex polytope is a face. Hence F, its q nontriangular ridge
neighbors, and its s nonsquare-cusp opposites are distinct positive facets.
The non-minimal inequality 24>=d(F)+m(q+s) follows. No uniqueness of
opposites across different choices of F is needed or claimed.

## N=7 really forces the excluded cuboctahedral case

For each D9 facet under N=7, the one facet, its three quadrilateral
neighbors, and its three nonsquare-cusp opposites exhaust the positive
facets. Therefore no triangular ridge of a D9 can have a positive
neighbor. The excess defect above seven copies of three is exactly three:
the only partitions are6;5+4;4+4+4.

A defect-four facet has four quadrilaterals paired by edge adjacency.
If one of its quadrilateral neighbors were D9, the neighbor across the
paired quadrilateral would cut an edge-adjacent face in that D9. That face
cannot be a quadrilateral (D9 quadrilaterals do not share edges), and
cannot be a triangle (the neighbor is positive). All four neighbors must
therefore have defect greater than three, but there are at most two other
such facets. This rules out every partition involving defect four.

The remaining defect-six facet has six quadrilaterals, no two sharing an
edge. Every one of its twelve four-valent vertices has exactly two
alternating quadrilaterals: there are24 incidences and at most two per
vertex. Joining quadrilateral centers through these vertices gives a
simple four-regular six-vertex embedded graph, hence the octahedral graph.
Its eight triangular complementary regions are distinct three-cycles;
these exhaust the graph's eight three-cycles. Thus its medial recovery
has the full cuboctahedral face lattice, rather than merely a matching
abstract graph. The displayed cuboctahedral normal model gives twelve
nonsquare cusp rectangles. Rigidity then transfers that metric conclusion
to the facet, and injection requires19 positive facets, contradicting7.
I found no missing N=7 partition or tacit four-dimensional census here.

## The unit-normal width equalities have the correct scaling

There is a direct upper-half-space normalization check. With cusp null
vector xi=(1,0,...,0,1), the vertical wall x_i=a has unit normal
(a,e_i,a). The opposite wall x_i=b has outward unit normal
(-b,-e_i,-b). Their sum is (a-b)xi. Thus the coefficient of the null
vector is proportional to the width, not its inverse, with positive
sign for the outward pair.

At a2:1:1 cusp the four D9 belt facets bound the two equal short
directions. Equal widths therefore give equal opposite-pair sums with
no free rescaling of individual normals. These are all six cube-face
parallelogram identities. Their affine form n_s=a+sum(s_i*b_i), followed
by all unit equations and cube-edge orthogonalities, gives

a^2=-1/2, b_i^2=1/2, a.b_i=b_i.b_j=0 (i!=j).

The positive opposite-pair sums fix the future orientation. The resulting
Lorentz four-space U defines a genuine hyperbolic hyperplane H, and all
D9 supporting walls are orthogonal to it. The six octahedral ideal
vertices in H are exactly the actual six nonsquare cusps. Convexity
therefore supplies the inclusion S subset P; the eight D9 inequalities
supply P intersect H subset S. No formal arrangement is being substituted
for the actual section P intersect H=S.

## H produces actual right-angled entirely ideal halves

The face interiors of S lie in the relative interiors of D9 facets, and
its edge interiors are internal diagonals of actual quadrilateral ridges.
Consequently no ordinary octahedral facet can contain a finite point of S:
a supporting functional vanishing at a relative-interior point vanishes
on that minimal face, contradicting the listed facet/ridge interiors.
If an ordinary octahedral facet crossed H, a finite segment in that facet
would cross S, contradicting this fact. Thus only D9 facets are cut, each
along its genuine internal gluing triangle into two regular octahedra.

One useful explicit interior check is this: the center of S satisfies all
eight D9 inequalities strictly. No other facet meets S at any finite
point, so it satisfies all other supporting inequalities strictly too.
It is an interior point of P, and both halves have nonempty interior.

In the Klein model, a vertex of P intersect H^+ lying on H must be a
vertex of the section S; otherwise a nontrivial segment in S would
contradict extremality. A vertex away from H is an old vertex of P.
The section has only the six old ideal vertices, so no finite vertex
can be introduced. The same applies to the other side. New dihedral
angles involving S occur only with the cut D9 walls, which are orthogonal
to H; all other angles are inherited. This proves the needed entirely
ideal and right-angled assertions for the actual halves, not just their
facet incidence model. Their facets are all octahedra, hence all ridges
are triangular. The volume identity and published minimum uniqueness
then identify them as24-cells. Extending octahedral facet symmetries by
24-cell flag transitivity supplies the final untwisted reflection gluing.

## Precise remaining dependencies and verification limits

The read uses the standard rectangular ideal link and acute-chamber
intersection facts; the established volume and defect identities;
the unique minimum24-cell theorem; the antiprism bound with equality;
Andreev-type rigidity for the admitted ideal right-angled3-facets; the
complete admissible nine- and ten-vertex facet list; and regular24-cell
flag transitivity. The local D9 geometry is supplied by its explicit
double-octahedron model, not by the four-dimensional equality conclusion.

I independently read the primary classification completeness text in
Kolpakov, pages1951–1952, Section4.3 CaseD, and the opposite-facet counting
paragraph on page1949. The text explicitly says Figures6–7 give the
realizable cases through ten vertices. The four-quadrilateral matching
read from Figure7 remains an external incidence transcription in this
review: I did not obtain a fresh usable figure image and do not claim to
have independently rechecked its picture or the underlying historical
catalogue. This dependence is stated in the equality report itself.
Primary source: [Kolpakov2012](https://msp.org/agt/2012/12-4/agt-v12-n4-p03-s.pdf).

No unproved extra geometric implication was identified within the reviewed
steps once those dependencies are assumed. This does not supply end-to-end
Lean validation, external human review, or any classification of higher
volume garlands/doublings. The supporting local scripts alone would not
prove the universal equality theorem; its counting and geometric reduction
are the written arguments assessed above.
