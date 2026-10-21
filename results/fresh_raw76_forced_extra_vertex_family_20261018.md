# Raw76 quadrilaterals force an extra vertex: a bounded coupled-charge pilot

This packet concerns the compact arithmetic H13/Q(sqrt5) model. It adds a
new necessary incidence mechanism: every represented raw76 quadrilateral
has a unique fifth finite vertex outside the quadrilateral. A marker for
that occurrence can be charged together with ordinary corners at the
extra vertex, using a single vertex allowance.

The complete inventory contains **235 unique-context face classes** and
3,204 marker orientations. Physical compatibility is tested for only the
three currently highest-mass classes. Two resulting inequalities strictly
reject the frozen current327 dual. No global optimization, new fixed-mass
repair, new dimension bound, or physical realization is asserted.

|Face|Extra-vertex type|Markers|New pairs|Excluded|Retained unknown|Marker capacity|Best charge minus current allowance|
|---|---:|---:|---:|---:|---:|---:|---:|
|43092|1781|4|326|326|0|1|−4039181/10000000000|
|43012|944|6|1215|1215|0|1|241899/80000000|
|43010|904|12|2274|2268|6|2|23181393/20000000000|

All marker–ordinary-star pairs are excluded in these three domains. At
904, six marker–marker pairs remain allowed. They are unknown outcomes,
not asserted simultaneous physical realizations. Their retention gives
capacity2 rather than1.

The [producer](../scripts/fresh_raw76_forced_extra_vertex_family_20261018.py),
[exact data](fresh_raw76_forced_extra_vertex_family_20261018.json), and
[standalone saved checker](../scripts/check_raw76_forced_extra_vertex_saved_20261018.py)
bind the [current327 global model](fresh_current327_full_model_20261018.md).
The model's exact optimum interval is unchanged. Its4,624 face masses are
used only to select this bounded pilot and evaluate universally valid
charges; the inequalities themselves do not depend on those masses.

## An actual fifth vertex, not merely a positive principal matrix

Canonical core65, raw core76, has seven native normals labelled0 through6.
Its unique quadrilateral partition is `{1,3}|{0,2}`. The quadrilateral fixes
walls `{4,5,6}` together with its rank8 spectator forest. Its four corners
omit `{0,1}`, `{0,3}`, `{1,2}`, and `{2,3}` from the seven core walls.

Scale normal0 by the positive square root of `(7+sqrt5)/2`, leaving the
other six scales equal to1. The resulting Gram is over Q(sqrt5), with
inertia `(5,1,1)`. The source checks all21 five-root principal supports.
Exactly five are positive definite: the four quadrilateral corners and
the additional support

    {0,1,3,4,6}, omitting {2,5}.

Its finite type is `2A1+H3`. With the rank8 spectators it is a positive
definite13-wall support in the actual H13 chamber. It determines an
actual finite vertex even in the presence of arbitrary additional facets.
Here is the argument, also used in the
[frozen1325 incidence proof](fresh_forest1325_cross_vertex_incidence_20261017.md).
Choose a strict future timelike interior point h, write
`b_i=h·n_i<0`, and set `c=G^{-1}b` for the selected positive definite Gram.
Since G is acute, `G^{-1}` is entrywise nonnegative and `c<=0`. The vector

    x = h − sum_i c_i n_i

is perpendicular to the13 selected normals. For every other actual wall
m, acuteness gives

    x·m = h·m − sum_i c_i(n_i·m) <= h·m < 0.

Moreover `x²=h²−bᵀG^{-1}b<0` and `x·h=x²<0`. Thus x is future timelike,
and the13 independent walls leave exactly its projective line. The other
walls are strictly inactive at this vertex. In particular, this vertex
is outside the original quadrilateral because its support omits the
quadrilateral's fixed wall5.

## Recovery and exactly one marker per quadrilateral

At the extra vertex, the active incident walls are `{0,1,3,4,6}`. The
producer enumerates whole-component embeddings of these `2A1+H3` roles
into the complete finite13-wall frame. The two missing normals2 and5
have their complete native contact covectors saved, including zero
contacts with every unused spectator coordinate. Both normals have
square2. Their strictly positive squared timelike heights are checked
from the full inverse finite Gram; the common future-side choice also
reproduces their prescribed mutual product0.

A complete contact covector fixes the spacelike projection of a normal.
Norm2 and the common positive timelike height then fix its remaining
coefficient uniquely. Thus a marker at an actual finite vertex recovers
the two actual missing outward walls, not merely congruent diagram labels.
Together with the13 actual frame walls, it recovers the original15-wall
configuration. The fixed walls `{4,5,6}` and the eight spectators recover
the quadrilateral as their intersection with the actual polytope. Two
different quadrilaterals therefore cannot produce the same marker.

Conversely, each actual source quadrilateral has the one specified fifth
support and hence produces one marker. This counts actual occurrences,
not labelled embeddings. The seven-root native graph has trivial
automorphism, checked explicitly. The dotted versus native5 contacts
distinguish the two A1 roles, and the H3 native5/native3 chain fixes its
roles. The saved unordered two-covector marker keys are all distinct.
Permuting unused equal spectator components changes neither these
covectors nor the actual fixed-wall set and supplies no extra multiplicity.
With a A1 components and h H3 components in the target forest, the complete
orientation domain has `a(a−1)h` markers. These are alternative occurrences
at a vertex; each source quadrilateral still contributes coefficient1.

This is an incidence statement about markers that arise from actual
source quadrilaterals. An abstract covector pair that passes the necessary
tests is not asserted to realize any quadrilateral.

## The triangle-port exclusion and existing base accounting

Fix core walls `{0,4,6}` and the eight spectators. The three remaining
side walls `{1,3,5}` have three pairwise completions among the five saved
positive definite supports. The Schur residual has inertia `(2,1,0)` and
strictly negative inverse. Its three dual rays form a compact simplex.
All three vertices are actual by the strict-center argument above.

The actual fixed-wall face lies inside that three-side simplex. Conversely
all three simplex vertices lie in the actual polytope, so their convex
hull does too. The two sets are equal. Additional facets cannot truncate
this compact triangle. At the fifth vertex its two ports are the images
of core walls1 and3. The unique two-face determined by those ports is this
triangle. An ordinary quadrilateral using that same occupied port pair
therefore cannot coexist with the marker. This agrees with the inherited
[simplex-port theorem](fresh_simplex_disjoint_ports_20261017.md).

The triangle is already present in the
[older prism/base-simplex incidence mechanism](fresh_sqrt5_prism_simplex_base_incidence_20261017.md).
It is not a new independently charged triangle. The new information is
that the forced marker and ordinary Q-stars compete for the **same**
extra-vertex budget. The marker adds a source-quadrilateral coefficient,
with no added triangle demand or separate triangle allowance.

## Complete inventory and three conservative physical domains

The235 classes are all existing canonical core65 face rows with one
unambiguous source context. The source checks their rank8 spectator
forests and any explicit provenance annotation. Merged face rows are not
silently assigned to one provenance. The235 extra-vertex types are
distinct. They have3,204 marker orientations and39,966 within-target
marker–marker pair positions in total. Only16 of the235 classes currently
have positive mass, with total `57192559/5000000000`. These are catalogue
and workload counts, not235 physical compatibility certificates.

The three tested classes have spectators `A3+H3+I2(5)`, `A1+A3+A4`, and
`A1+A2+H3+I2(5)`, respectively. Their extra-vertex types are1781,944,904.
The complete old native graphs are reconstructed chronologically from
the290-class domain, historical32 additions, raw4 additions, and the two
later raw76 patches. All old nodes, ordinary occurrence coefficients and
Q–Q edges are preserved literally, including every historically owned
unknown and the two unknown pairs at forest1763 outside this pilot.

Each newly tested pair contains at least one marker. There are3,815 such
positions, and no old Q–Q outcome is retested or replaced. Every physical
comparison selects a fresh immutable squareclass field from full native
covectors. The strict self heights and common-side cross products are
reconstructed exactly. Positive products, the forbidden interval `(−1,0)`,
or failure of the arithmetic two-cycle field/integrality tests give
necessary exclusions. Identical individual wall covectors denote the same
wall and are not tested against themselves. All six unresolved pairs stay
as compatibility edges. No higher-order physical realizability is claimed.

## A single coupled vertex charge

Let w be a previously proved nonnegative ordinary-corner charge, normalized
to cost1. For each maximal clique of the enlarged conservative graph,
record its old charge x and marker count k. The source completely enumerates
12,804 maximal cliques over the three domains; a guard rejects a truncated
census. It then computes the exact feasible polygon

    lambda >= 0, t >= 0,  lambda*x + t*k <= 1.

Its rational vertices yield mass-independent inequalities. Nonnegativity
of both weights extends every maximal-clique check to all subcliques.
Summing over actual finite vertices gives the old ordinary face coefficients
multiplied by lambda, plus `t` times the source raw76 face count, bounded
by one target-vertex count. Old occurrence coefficients already count
actual corners with their multiplicities. The extra marker has coefficient1
by the recovery argument above. There is no second vertex allowance.

For this pilot every marker–ordinary pair is excluded, so an old charge of
cost1 combines with a marker weight1,1,1/2, respectively, without scaling
the old charge. The best violated inequalities are

    2 a26328 + a26394 + a26410 + (1/2) a26864
      + a28678 + a28760 + a43012 <= v944,

    a26304 + a26362 + a26392 + a26627
      + (1/2) a26862 + (1/2) a43010 + a45821 <= v904.

They exceed the corresponding current exact dual allowances by
`241899/80000000` and `23181393/20000000000`. The first target1781 retains
enough allowance slack that its analogous inequality is not violated.
These are new valid candidate columns, not a reoptimization of the frozen
model. All14 existing multi-owner columns remain in that model unchanged.
No claim is made that this restricted two-parameter charge family contains
every facet of the enlarged graph cone or optimizes a new local demand.

## Verification and limits

Initial bounded generation session31929, terminal chunk9e85d9, passed.
Final source-bound generation session31368, terminal chunke4f19e, and
separate fresh replay session53243, terminal chunk9b3f1b, both passed exit0.
The replay reproduced the complete geometry, pair decisions, graphs and
charges. No local LP or global optimizer is used.

The standalone checker uses only the Python standard library. It implements
its own Q(sqrt5) arithmetic and ordered multiquadratic fields; it verifies
the saved21-support geometry, compact triangle inverse, native symmetry,
235-class catalogue accounting, marker maps and full-frame zero blocks,
strict heights and physical product exclusions. It independently enumerates
the saved graphs' maximal cliques and checks every charge and current-mass
violation. Historical owner chronology and the original catalogue/native
embedding completeness remain producer/dependency obligations. This is
not an independent reclassification of Coxeter diagrams or an independent
proof of the geometric counting argument.

The final hardened checker session90736, terminal chunkd8745d, passed
exit0:21 supports, five positive supports,235 family rows,22 marker maps,
1,317 independent wall expressions,3,809 pair exclusions,12,804 maximal
cliques,20 charge-envelope inequalities and ten stored charge candidates.
All six unknowns are retained. The earlier checker run exposed a saved
frame schema mismatch, which was corrected before its successful replay;
the final version also explicitly checks both saved Gram symmetries.

The root completed a full producer/report/helper mathematical read,
including the strict-center actual-vertex argument, compact triangle,
complete-frame recovery and once-per-quad count, active-role multiplicity,
triangle-port exclusion without another payment, preserved Q–Q outcomes,
nonnegative charge polygon, and the helper's independent field, graph and
charge arithmetic. The construction peer independently completed the same
full source/report/helper mathematical-scope read, additionally rereading
the inherited whole-component embedding/projection and strict-center/
triangle-port interfaces. The peer checked the final finite-frame schema
and Gram-symmetry guards. Neither reviewer executed the package or repeated
physical projections, pair tests, clique enumeration or arithmetic replay.
No external referee or proof-assistant validation is claimed.

The complete235-class inventory is distinguished throughout from the three
tested domains. No untested marker–ordinary compatibility is discarded.
The frozen current327 model and all preceding packets remain unchanged.

    python3 scripts/fresh_raw76_forced_extra_vertex_family_20261018.py --verify
    python3 scripts/check_raw76_forced_extra_vertex_saved_20261018.py
