# Raw357 forces three extra finite vertices

An actual raw357 quadrilateral in the compact arithmetic H13/Q(sqrt5)
model forces three further finite vertices. Their complete native frame
types are2A2+H4, A1+2A2+H3 and3A2+I2(5), together with the original
rank5 spectator forest. All three lie on an actual compact4-simplex.
This supplies a new incidence mechanism; it supplies no capacity,
exported inequality, fixed-mass improvement or H13 exclusion yet.

The [producer](../scripts/fresh_raw357_extra_vertex_incidence_pilot_20261018.py)
and [saved certificate](fresh_raw357_extra_vertex_incidence_pilot_20261018.json)
bind the [new varying-mass model](fresh_all_marker_families_full_model_20261018.md)
and its [support inventory](fresh_all_marker_families_model_support_20261018.json).
That model has certified lower bound12.339709954130539... . Its nine
unambiguous positive raw357 classes have total multiplier mass
83284683/10000000000. This mass selects the pilot; none of the geometric
statements depends on it.

## Complete core and source inventory

Raw357 is canonical core318 in the frozen catalogue. It has ten walls,
physical inertia(8,1,1), conjugate inertia(9,0,1), and spectator rank5.
Its non-right native edges are:

* dotted0–2, with squared norm2 contact14+6sqrt5;
* label5 on1–3;
* label3 on2–3,4–5,3–6,5–6,2–7,7–8,0–9.

Every other pair is right. The exact native automorphism group is trivial.
The catalogue contains one quadrilateral partition, with opposite pairs
{1,3} and{0,2}; the aggregate demand is1.

All45 deletions of two core roots are checked exactly. Seven leave a
positive-definite8-root support:

    {0,1}, {0,3}, {1,2}, {2,3}, {2,4}, {2,5}, {2,6}.

The first four are precisely the original quad vertices. The last three
are extra vertices. This is subtraction from the complete actual source
quad's vertex inventory, rather than from one arbitrarily pointed view.
Deleting a spectator instead cannot create another finite vertex within
these15 walls: at least nine core roots remain, whereas the complete core
has positive index8. Interlacing therefore excludes a positive-definite
full13-root support. The45 core tests exhaust this within-star domain.

The complete inherited spectator catalogue gives47 unambiguous source
classes, all with exactly one source context and aggregate demand1.
All47 are retained, including the38 with zero mass in the current dual.
Each has three target occurrences, giving141 pointed target domains and
138 distinct target forest types. The nine positive classes are70400,
70404,70406,70408,70414,70420,70422,70426 and70445. The last has explicit
provenance in the model; its unique catalogue context is checked directly.

## Why the vertices and simplex are actual

Assume the full star occurs in a realizing compact Coxeter polytope.
Append all five spectator roots to any of the seven positive-definite
core supports. This gives13 linearly independent spacelike facet normals
with acute Gram G. Let h be a strict future timelike interior vector,
let b_i=h·n_i<0, and put c=G^{-1}b. A positive-definite acute Gram has
entrywise nonnegative inverse, hence c_i≤0. Then

    x=h−sum_i c_i n_i

is future timelike and orthogonal to the13 selected walls. For every
other actual facet normal m, acuteness gives

    x·m=h·m−sum_i c_i(n_i·m)≤h·m<0.

Consequently the selected walls define an actual finite vertex, and no
other facet is active there. Arbitrary additional facets cannot remove
any of these seven vertices. In particular, every such vertex is simple.

Fix the spectators and core walls0,7,8,9. The remaining side walls
1,3,4,5,6 have a residual Gram of inertia(4,1,0), with every entry of
its exact inverse strictly negative. They bound a compact4-simplex.
Its five vertices are the supports omitting{2,j}, for
j in{1,3,4,5,6}; all five have just been proved actual. The actual face
lies inside the simplex and contains all five vertices. Convexity gives
equality, so no unlisted facet truncates it.

At any of the three extra vertices, four simplex-side walls are incident.
Each pair of these four ports lies on an actual triangular2-face.
An ordinary quadrilateral cannot occupy that same port pair at the
simple target vertex. This gives six necessary port exclusions per
marker. It is a compatibility statement and does not add another
triangle payment or increase the original quad demand.

## Full normal recovery and multiplicity

|Omitted core roots|Active finite frame|Formal orientations at a target|
|---|---|---:|
|2,4|2A2+H4|4b(b−1)h4|
|2,5|A1+2A2+H3|4a b(b−1)h3|
|2,6|3A2+I2(5)|16b(b−1)(b−2)p|

Here a,b,h3,h4,p count whole target components A1,A2,H3,H4,I2(5).
The A2 roles are ordered. Each A2 has two native reflections, H3 and
H4 have none, and I2(5) has two; these give the displayed factors4
and16. The producer independently compares the formulas with complete
component embedding counts. Summing across the47 classes gives552,
1144 and8928 orientations, respectively:10624 in total.

These orientation counts are local domain sizes, not occurrence charges.
Each source quad contributes exactly one occurrence at each specified
extra vertex, three in total. The three vertices are distinct because a
wall omitted at one is incident at the other. The original source demand
remains1; it does not become3 or10624.

At each extra target, record both missing walls' covectors against the
entire actual finite frame, with zero contacts on all spectator coordinates.
The exact projection has two strictly positive timelike height squares.
Subtracting the positive product-height root recovers the prescribed
mutual native product. A complete covector fixes the spacelike projection;
norm2 and the common future-side sign fix the remaining coefficient.
Thus the marker uniquely recovers both missing actual wall normals.

The recovered full frame and normals determine the connected nonfinite
ten-wall core and its orthogonal finite spectators. Its trivial native
automorphism group leaves no orientation stabilizer to divide or multiply.
Its unique quadrilateral partition identifies the same source quad.
The finite source vertices and their convex hull identify its actual
face and full star. Passive component permutations do not change it.
This is the reverse recovery needed for once-per-quad incidence; a future
capacity proof must still count simultaneous distinct markers correctly.

The bounded native pilot uses highest-mass source70408, whose spectator
forest is2A1+H3. Its three complete target domains have8,48 and96 maps.
All152 maps check the full native covectors, positive heights, missing-wall
product and duplicate-free unordered recovered-wall keys. The map domain
is checked against complete whole-component embeddings. Other source
classes are included in the47-class/141-domain inventory and the general
recovery argument; this packet does not claim to have numerically replayed
all10624 native maps.

## Workload and next structural question

The current ordinary327/600 domains are reconstructed by chronological
whole-orbit references. Their literal-star counts are kept distinct from
count-vector representatives. The existing raw56/57/76/108/112 marker
inventories are added once each. The unreduced future pair positions are:

|Pair domain|Positions|
|---|---:|
|New marker versus current ordinary star|1743264|
|New marker versus existing marker|41984|
|New marker versus new marker|1610304|

These counts are a workload estimate, not compatible-pair counts. No pair
has been retested and all historical owners and unknowns remain unchanged.
The next useful step is a complete relative-role reduction using ordered
A2 components and their endpoint reflections, followed by exact native
tests and a safe mixed-marker capacity with one vertex allowance.
In particular, the large3A2/I2(5) orientation domain should not trigger
an unreduced million-pair scan.

The six simplex port exclusions and complete two-wall normal recovery
are available for that next step. They do not by themselves give a
nontrivial marker-clique capacity or a valid shared ordinary/marker charge.
No fixed-mass gain is inferred from the positive raw357 mass, and no new
optimization is performed here.

## Verification and limits

Final producer generation33721/8145e3, separate fresh replay41121/8ee999,
and [standalone saved checker](../scripts/check_raw357_incidence_saved_20261018.py)
92616/06aa36 all terminated exit0/PASS. Two earlier producer attempts
stopped while adapting the old-domain count schema, before any certificate
was written: literal count-vector representatives are not native stars.
The final code uses complete corner-orbit star sets and checks their
disjoint contiguous index domain.

The checker imports no producer. Its stdlib Q(sqrt5), ordered
multiquadratic, positive-definiteness and Gauss–Jordan routines are reused
from the independent four-core checker. It independently checks all45
support decisions, a nonsingular9-root minor plus full null vector,
native scaling and automorphisms, the simplex inverse, all three height
identities, complete47/141 inventory and orientation formulas, all152
saved native maps, and the pair-workload arithmetic. Contacted finite
blocks are inverted exactly; untouched passive blocks may be replaced
by diagonal placeholders only after checking zero contacts and complete
component separation. It does not claim to invert irrelevant actual
passive blocks.

Historical core/spectator classification, actual-face existence arguments
and old ordinary ownership completeness remain explicit proof inputs.
The helper does not replay a physical pair census, compatibility graph,
local cover or global matrix. No new capacity or column is exported.
The root completed a full producer/report/helper mathematical read and
reread the inherited scaling, active-graph, forest-type and embedding-count
interfaces. That review checked the actual-vertex projection, complete
within-star support domain, simplex equality, two-wall future-sign
recovery, trivial-automorphism/unique-quad multiplicity, orientation versus
occurrence counts, complete47/141 inventory versus152-map pilot, and the
workload-only scope. It repeated no execution, exact support/native
arithmetic, historical classification or owner census. Final integrity
checks verify direct source hashes, both Python ASTs and local report links.

    python3 scripts/fresh_raw357_extra_vertex_incidence_pilot_20261018.py --verify
    python3 scripts/check_raw357_incidence_saved_20261018.py
