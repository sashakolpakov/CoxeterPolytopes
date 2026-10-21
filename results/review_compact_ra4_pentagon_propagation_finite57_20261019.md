# Independent review of compact RA4 pentagon propagation

19 October2026. **Written mathematical review PASS**, conditional on the
explicit inherited acute/capacity lemmas and the named compact3D metric
uniqueness/classification inputs. No job, checker body, numerical program,
census, Lean or build ran. No reviewed source or pinned evidence was edited.

Read the complete regular-pentagon propagation source, the earlier W=5
note, the neighbor/zero-excess note and root review, and the full original
perpendicular-capacity proof and its root review/scope audit. This review
does not infer a stronger result from an agent message.

The principal reviewed source is
`compact_ra4_regular_pentagon_propagation_20261019.md`, SHA256
`0536b40623cc942c529865a3e56ac067f5a89ba46bed2ad341fb2dc88aec69b7`.
The earlier W=5 source is pinned at
`255a16112211991134cc9d9e97010d5dc4b1f995c285f7ea6a02976c0927f5cd`;
the neighbor/zero-excess source at
`9c670a7b3c25dda48997c3f70a6d1be9be88bcdefb475b7710f5011068aebcbe`;
the capacity source at
`8e6db87b001b59a16c1b9b55f03fa9373630bcfb20c1697796f4fd88fbaff8ff`.
Its conclusions
are W=0 or W>=9, and the stronger local implication that any degree-two
positive-excess facet forces W>=15. The possible forthcoming full-adjacency
minimum-five/small-graph extension is NOT reviewed by this note.

## The degree-two implication checks

The weighted degree at F is exactly d(F). A k-gonal ridge forces
d(F)>=2(k-5) by the inherited3D capacity theorem. With precisely two
nonpentagonal ridges, their weights are therefore equal and equality
in the capacity theorem applies. The full closing proof supplies the
Loebell type, including TWO disjoint bases and exactly2k lateral
pentagons, each incident to an entire base edge. No unproved3D catalogue
is being substituted for this equality argument.

The compact marked Andreev uniqueness clause realizes the one-step
combinatorial rotation as an isometry. Its restriction to a convex base
is the one-step orientation-preserving rotation, so the base is a regular
right-angled k-gon. The center/vertex/edge-midpoint right triangle gives
cosh(L_k/2)=sqrt(2)cos(pi/k) and cosh(L_k)=1+2cos(2pi/k). Hence for k>=6
its edge length exceeds L_5, whose cosh is phi.

Each lateral pentagon consequently differs METRICALLY from the regular
right-angled pentagon in a zero-excess dodecahedral facet. Its neighbor
must be positive, even though the shared ridge is itself pentagonal.
All2k lateral neighbors and both base neighbors are distinct: distinct
facets of a convex polytope cannot share two separate ridges. Summing
their excesses and F's against the exact identity sum d=2W gives

    2W>=2(k-5)+4(k-5)+4k=10k-30.

There is no disjoint-charge assumption: facets, rather than individual
ridge charges, are summed. Thus W>=5k-15>=15 is valid.

## Four equal edges really force a pentagon

A directed four-edge chain with all lengths L_5 and coherent right-angle
inward turns has a unique realization up to isometry. Matching it to
the regular pentagon fixes both endpoint perpendicular supporting lines;
they are the same fifth line. A full-dimensional convex polygon has a
single connected intersection face with that supporting line. The segment
between the two endpoint vertices must therefore be the fifth edge.
It cannot be split by genuine additional right-angle vertices or replaced
by a detour. The inward-turn hypothesis excludes an oppositely oriented
chain. Thus the lemma applies as stated, without a hexagon formula.

For a nonpentagonal actual ridge, a zero-excess third facet supplies an
entire regular-pentagon edge, hence length L_5. If exactly two positive
carriers are adjacent, all remaining k-2>=4 edges form one consecutive
zero-carrier block, contradicting the lemma. A clique of four ACTUAL
facet normals gives their actual common finite vertex by the reviewed
orthogonal-projection lemma, so the required adjacency is geometric.

## The small full-graph argument checks

For 0<W<=8 the degree-two implication and the previous degree-one
exclusion give minimum degree three in G. Since its unweighted edge count
is at most W, there are at most five positive facets. The proof correctly
passes to FULL adjacency H, which includes pentagonal ridges; an induced
square of G alone would not be a contradiction.

With four vertices H=K4. The two remaining positive facets give the only
two positive edge-carriers of any selected nonpentagonal ridge, and their
four-normal clique makes them adjacent. With five vertices, the complement
of H is a matching. Zero missing edges violates positive index four;
two missing edges give an induced square, whose compact normal Gram
has two negative directions. The sole remaining H=K5 minus uv forces
all three edges incident to u into G. Choosing ua gives exactly b,c as
the positive ridge carriers; v is disjoint from u. Their clique makes
the two carriers adjacent, so the polygon lemma again applies. All cases
are covered and W=0 or W>=9 follows.

The strict nonedge normal gap is used with COMPACTNESS. A merely acute
right-angled noncompact graph would need a separate treatment of null
pairs. The inherited capacity packet proves that strict gap and the
actual-face projection rather than assuming abstract flag incidence.

## Earlier W=5 metric note and result boundary

The earlier C5 argument also checks on its own. With W=5 the graph is a
unit-weight five-cycle, whose facets are L6. The one-step symmetry from
marked metric uniqueness puts the two belt radii at sqrt(2), time
coordinates cosh(h) and axial coordinates +/-sinh(h). The two cross-belt
orthogonality equations give their angle offset pi/6 and base-normal
product -sqrt(3). Thus the five ambient normals have Gram
I-sqrt(3)A(C5), with three negative eigenvalues because sqrt(3)>phi.
This contradicts Lorentz index one. Its stated limit to the specific
unit C5 component is correct; it does not exclude every abstract cycle.

The W=0 regular120-cell recognition separately relies on the stated
metric dodecahedral classification and propagates ALL dodecahedral
facets. The present proof does not propagate an arbitrary lone facet.
Neither W>=9 nor the degree-two W>=15 implication proves f+W>=120,
universal minimum volume, a second-volume gap, a global tessellation,
or a new dimension ceiling. Those boundaries remain essential.
