# Independent written review: octahedral-block decomposition

18 October2026. **PASS**, for the restricted theorem in
[the author's draft](ideal_ra4_octahedral_block_decomposition_research_20261018.md),
including the explicit finite-disjointness argument described below.
This review read the full draft, the full earlier
[second-volume equality argument](ideal_ra4_second_volume_equality_20261018.md),
and the opposite-cusp injection in the frozen
[gap argument](ideal_ra4_second_volume_20261018.md). It did not run a
program, inspect a four-polytope census, compile a manuscript, or execute
Lean. Frozen originals were not changed.

The hypothesis is essential: every facet is a regular ideal octahedron O
or its triangular-face double D9. The result gives a tree of24-cells in
that class. It does not classify all entirely ideal right-angled H4
polytopes, establish universal volume integrality, or prove that all
decompositions are chains.

The following points were checked independently.

* The local rectangle data force every noncubic ambient cusp box to have
  ratios2:1:1. The forbidden2:2:1 case would give four quadrilateral
  incidences at a square end face, allowed by neither O nor D9. In the
  surviving box, exactly four nonsquare D9 faces form one quadridge belt;
  the two ends are O. Each internal triangle cuts its rectangle along
  the segment transverse to the long direction at its midpoint. The four
  segments form one embedded square circle, so the compactified triangle
  complex has a disc neighborhood at the cusp, with no crossing or branch.
* At finite points triangle interiors lie in distinct facet interiors,
  and their edge interiors lie in distinct quadrilateral ridge interiors.
  Thus no additional finite identifications are being assumed away. The
  complex is an actual closed surface. A cusp cannot join separate graph
  components, since its four triangles already lie on one belt cycle.
* For a component with n triangles, each edge has two incident triangles
  and each vertex four, giving chi=n/4. The primal-tree/dual-tree proof of
  chi<=2 is valid without orientability or a surface-classification input.
  A dual cut would give a nonempty even-degree subgraph of a forest.
* The seven-facet injection remains inside the chosen component. The three
  box-opposite facets at a D9's nonsquare cusps are the opposite members
  of those same four-facet belt cycles, hence quadridge-connected to it.
  Their walls are asymptotic to the chosen facet at distinct cusps, making
  them mutually distinct and different from the three ridge neighbors.
  Consequently7<=n<=8 and4 divides n, so n=8.
* The component's six-cusp graph is simple: two distinct ridges with the
  same two cusp endpoints would put a pair of opposite belt facets
  through both endpoints, contradicting their unique asymptotic point.
  Its simple4-regular graph is therefore the octahedral graph. Distinct
  internal triangles cannot repeat a cusp triple, since their finite
  interiors lie in distinct supporting facets. This gives the cube
  labeling and all six parallelogram equations for outward unit normals.
* Equal short box widths imply equal opposite-normal sums. The cube
  equations and the eight norm equations give the Lorentz4-space with
  a^2=-1/2 and b_i^2=1/2. The resulting regular ideal octahedron S is
  contained in P by its six actual cusp vertices and convexity; the eight
  wall inequalities give P cap H=S in the opposite direction. This is an
  actual global section, not merely a formal local wall configuration.
* A facet outside the eight cannot touch S finitely: a supporting
  functional's zero set on S is a face, while all triangular face
  interiors and edge interiors already belong respectively to their
  designated facet and their two designated ridge facets. Hence outside
  facets cannot cross H. The new pieces have no finite vertices, and
  their new dihedral angles are right angles.
* For clarity the finite-disjointness implication must remain explicit.
  Every other D9 component lies on one side of H: its facets are uncut,
  and two quadridge neighbors on opposite sides would put a finite point
  of their shared ridge in S, violating the preceding no-contact fact.
  Its cusps cannot be vertices of S because each noncubic cusp has one
  belt component. Their convex hull, its later section, is therefore
  disjoint from S, even at ideal points. Equivalently in recursive pieces
  an earlier S is an outside O facet, to which the same no-contact lemma
  applies. The tree conclusion follows by successive separating cuts.
* The direct all-O recognition is sound. A shared ideal triangle and the
  right-angle condition determine the neighboring supporting wall, and a
  regular octahedron on that wall is determined by its triangle and side.
  Octahedral facet symmetry cancels any triangle automorphism. Matching
  propagates through the connected facet graph of the model24-cell. Its
 24 actual walls give containment in the model, and its actual cusp
  vertices give the reverse containment. No universal decomposition or
  published four-dimensional classification is hidden in that step.

Thus N=8k, k+1 cells, and Vol=(k+1)4pi^2/3 follow within the stated facet
class. The defect identity independently agrees: sum d=3N, W=3N/2=12k.
This is a written mathematical review, not end-to-end kernel verification
or an assessment of worldwide priority.
