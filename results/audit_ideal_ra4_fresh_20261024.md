# Fresh audit of ideal right-angled four-volume and structure

24 October 2026. Reviewer: `/root/audit_sharp_ra`.

**The complete written second-volume, unique-double and restricted chain
arguments survive this fresh mathematical audit at their declared
external-input boundary.** I found no omitted defect case, circular
equality argument, unidentified wall in a section, invalid cut, or
unjustified gluing parameter. I also read every selected finite-model
producer/checker and found its claimed mathematical scope consistent
with the written proof.

One external incidence check remains limited: the original publication's
complete small-facet list was retrieved as text, but a fresh visual
transcription of its Figure 7 was unavailable through the PDF image
service. The four-quadrilateral matching therefore remains the explicit
external incidence premise already stated in the manuscript. I checked
the stored cycles and their matching directly; I do not relabel that as
a new visual check of the original figure or a reproof of the historical
classification.

This agent supplied no derivation or edit to the audited frozen sources
before this audit. Previous PASS labels were not used to infer the
deductions below. No manuscript, PDF, worker, certificate or Lean source
was changed. No worker import, mathematical enumeration, numerical or
symbolic program, test, Lean or TeX process ran. All mathematical checks
below are manual source reconstruction. Earlier execution records retain
their original bindings and narrower scopes.

## 1. Sources read and authenticated

I read the ENTIRE 1,102-line publication TeX before consulting its older
audit verdict. I also read the complete gap, equality and restricted
decomposition sources, the path-rigidity source, both model producers,
both independent saved-data checkers, the D4 checker, and both selected
conditional Lean modules. The selected sum/partition lemmas in
`CertificateBounds.lean` were inspected at their definitions.

|Source|Lines|SHA256|
|---|---:|---|
|`manuscripts/ideal_ra4_volume_and_structure_publication.tex`|1102|`e0c2ca24cec24ffcd736e201d3c7e06d75b28ef7be8093decd321ae3ac54bda1`|
|`manuscripts/build/ideal_ra4_volume_and_structure_publication.pdf`|—|`d95d1e29abc01bee7b4cf3d66f2c760871e1f6ac9f725fb4a017312c0ab69c47`|
|`results/ideal_ra4_second_volume_20261018.md`|178|`9dff11c375bddcf62f16f3f01ddabd7f5b5c137756cf02ebff275076e6869c25`|
|`results/ideal_ra4_second_volume_equality_20261018.md`|166|`4d227e9f7f4577f0a61b30a3ebf2790c763feeb4f1acddf7799a83f96fc1d7e2`|
|`results/ideal_ra4_octahedral_block_decomposition_research_20261018.md`|229|`664aa3925d61bef86ac257db946c6280528539156036bdc2dc923053ceda7e81`|
|`results/review_ideal_ra4_octahedral_block_path_rigidity_20261018.md`|85|`f4f57cae61ffe6581a79a0fcc4572e32a7786bf56e4ea00d1811585f3329bf85`|
|`scripts/fresh_ideal_ra4_double24_volume_20261018.py`|180|`16731c71d949e287d4cdce85e0ef98f8e5551beaaef7a6acfdd196c1fb0f5323`|
|`scripts/check_ideal_ra4_second_volume_saved_20261018.py`|196|`05ebb4f12e7413a13c82da8be366fe874933c93c3f6e57e8bc9bb4e932f64c09`|
|`scripts/fresh_ideal_ra4_equality_structure_20261018.py`|225|`298b3722b6aeceecab918f3c137882bf4f672cf94fb7fe646a6ac724d353ef5f`|
|`scripts/check_ideal_ra4_equality_structure_saved_20261018.py`|148|`464fa7261976233865b170bb747d9fa499e498860d1b337d2fc8ae3bfd096b00`|
|`scripts/check_ideal_ra4_octahedral_path_roots_20261018.py`|93|`87dca49df0d666db7f65bb541d2997bb1c625ef9e11f47f2e66ec479dfa1548a`|
|`formal/CoxeterProofs/IdealRA4Gap.lean`|76|`a39656988e358c9d53323f0927e216080cdb89e7251aa9f83f187c4ae7257a24`|
|`formal/CoxeterProofs/ManuscriptIdeal4.lean`|314|`c7726547a68d8617806436138c3e10d2f5dac338229b5b537c91685b163b295f`|

The PDF is hash-bound, not newly visually reviewed. The current proof map
is `docs/ideal_ra4_volume_and_structure_proof_map.md`, SHA256
`2a907af42e075e909ce5c259978fd5c6339e29d1d861d69ac07bf7dfca6d990f`.

## 2. Statements and exact external boundary

The universal gap and equality theorem applies to convex, finite-sided,
finite-volume hyperbolic four-polytopes with EVERY vertex ideal and
EVERY actual dihedral angle right. No arithmeticity is assumed.
The chain theorem separately requires every intrinsic facet to be a
regular ideal octahedron O or its triangular-face reflection double D.
This additional hypothesis is not assumed in the universal equality
proof. It is established there before the chain theorem is invoked.

The original [Kolpakov journal paper](https://msp.org/agt/2012/12-4/agt-v12-n4-p03-p.pdf)
was retrieved directly. The following are external inputs, with these
precise scopes:

|Input|Original location|
|---|---|
|Ideal RA4 incidence identities and volume formula|Lemmas 1–2, printed pp.1943–1944|
|Unique minimum attained by the regular ideal 24-cell|Theorem 3, printed p.1944|
|A k-gonal face in ideal RA3 forces at least 2k vertices, with antiprism equality|Proposition 5, printed pp.1945–1946|
|Complete admissible small octahedrite list, including the unique nine- and ten-vertex entries|Section 4.3, Case D, printed pp.1951–1953, Figures 6–7|
|Three-dimensional metric rigidity|The stated Andreev rigidity input, also invoked there|

The text directly confirms completeness and admissibility of the small
list. Its historical census and the published proofs remain external.
For the ten-vertex entry the EXTRA incidence premise is that its four
quadrilaterals have edge-adjacency graph two disjoint edges. That
transcription is addressed in Section 10 below. No classification of
four-polytopes or universal 24-cell decomposition is imported.

All-wall acuteness is a separate standard geometric input. It was
checked directly in [Vinberg, printed p.42, Section 2.1](https://www.mathnet.ru/php/getFT.phtml?jrnid=rm&option_lang=eng&paperid=2140&what=fullteng)
in the immediately preceding sharp-dimension audit. It is not being
deduced from the finite model checks here. The same finite-point and
Klein-face foundations justify the ideal right-angled three-dimensional
facets and their rectangular cusp links.

## 3. Defect identities and the universal volume gap

The rectangular ambient cusp has eight vertices and twelve edges. Thus
the ambient incidence counts are

\[
e=4v,\quad r=3v+f,\quad\sum_R k(R)=12v,
\quad W=3(v-f).
\]

For a four-valent ideal facet, (e_F=2v_F), (f_F=v_F+2), and the
polygon-edge sum is (4v_F). Consequently

\[
d(F)=\sum_{R\subset F}(k(R)-3)=v_F-6,
\qquad\sum_Fd(F)=2W.
\]

The published formula gives

\[
\operatorname{Vol}(P)=4\pi^2/3+W\pi^2/9,
\qquad W\in3\mathbb Z_{\ge0}.
\]

Only the smaller unit (\pi^2/3) is available universally. No
divisibility by the 24-cell volume has been assumed.

A positive-defect facet contains a k-gon with (k\ge4), so the
antiprism bound gives (d\ge2k-6\ge2). If there is just one
nontriangle, its defect (k-3) contradicts this bound. If there are
exactly two of sizes (k\ge l\ge4), then

\[
(k-3)+(l-3)\ge2(k-3)
\]

forces (l=k) and equality in the antiprism bound. This genuinely
identifies an antiprism; it does not only identify its vertex count.

For the injection, every nontriangle has its distinct positive-defect
ridge neighbor. At a nonsquare cusp the ambient box's opposite face
has the same nonsquare intrinsic rectangle, so belongs to another
positive-defect facet. Opposite wall planes are asymptotic at that
ONE ideal point. Hence an opposite facet cannot repeat at a second
cusp or be a ridge neighbor of the original facet. These are distinct
ACTUAL facets, not repeated cusp incidences. If m is minimum positive
defect, the valid simultaneous payment is

\[
2W\ge d(F)+m(q(F)+s(F)).
\]

The antiprism cusp calculation uses unit normalization consistently.
I expanded the four displayed normals: adjacent products are zero;
each opposite pair has unit product (-1). Their two unit-normal sums
are positive multiples (2u/\mu) and ((h/r+1/\mu)u). With

\(h^2+r^2=1\), (bc=r\), (\mu r=2h\), the ratio is

\[
{2/\mu\over h/r+1/\mu}={2r^2\over1+h^2}
 =2\cos(\pi/k).
\]

It exceeds one for every k at least four, with no truncation on k.
The fixed geometric antiprism and its symmetries are justified by the
published metric equality/rigidity input. D has exactly three nonsquare
cusps and three quadrilaterals; its other six cusps are square. The
nine-vertex complete list identifies it uniquely.

The complete minimum-defect split is:

|Selected facet|Lower bound for (2W)|
|---|---:|
|Antiprism, (m\ge2)|(m(m+9)\ge22)|
|Non-antiprism, (m=2)|Impossible: exactly two quadrilaterals|
|(m=3), unique nine-vertex facet|(3(1+3+3)=21)|
|(m=4), not the already treated pentagonal antiprism|(4(1+4)=20)|
|Non-antiprism, (m\ge5)|(4m\ge20)|

Thus (W\ge10), and divisibility by three forces (W\ge12).
This excludes precisely the prospective intermediate defects 3, 6, 9
without bounding the total number of facets or assuming a four-polytope
census. The gap proof is independent of the later O/D spectrum.

## 4. Equality spectrum: no positive-defect type is omitted

At the second volume (\sum d=24). The acute-intersection calculation
uses (u=p-\sum\langle p,n_i\rangle n_i) for orthogonal unit normals.
Its norm is negative, it is future, and every unselected wall product
is strict. It supplies actual common faces, not only formal pairwise
adjacencies. In dimension four, a four-clique would give a forbidden
finite vertex. In particular adjacency between two neighbors of F
anywhere in P forces their two faces in F to share an actual edge.

The minimum-defect reductions all close:

* An antiprism with (d\ge4) gives (24\ge3d+16\ge28), impossible.
* An A4 facet has eight nonsquare cusps. At each, the box needs another
  nonsquare face pair, forcing a positive triangular neighbor beyond
  its unique quadrilateral neighbor at that cusp. At least three
  distinct triangular neighbors cover the eight cusps. Together with
  F, two quadrilateral neighbors and eight distinct cusp opposites,
  this gives 14 positive facets and defect at least 28.
* For (m\ge5), at most four positive facets exist. Each needs three
  nontriangular neighbors, so exactly four form a clique. The actual
  intersection calculation excludes it.
* For (m=4), at most six positive facets exist. F has four quadrilateral
  neighbors, whose possible mutual adjacencies lie in the stated
  matching. Each needs at least three nontriangular neighbors, forcing
  the sixth facet and exactly three such neighbors. Hence none of the
  four neighbors has defect four; their defects are at least five.
  The full budget is (4+4\cdot5+4=28), impossible.

Therefore m is three and a D facet gives (7\le N\le8). At (N=8),
all defects are three. At (N=7), every D's three quadrilateral
neighbors and three cusp opposites exhaust all other positive facets,
so NO triangular face of D can have a positive neighbor.

The excess above seven defects of three is only three. Its partitions
are one six, one five plus one four, or three fours. A defect-four
facet's paired quadrilaterals cannot have a D neighbor: the paired
positive facet would create either edge-adjacent quadrilaterals or a
positive triangular neighbor in D. Thus all four quadrilateral
neighbors would have defect above three, but there are at most two
other such facets. Every partition involving a four is excluded.

The remaining facet has defect six and twelve vertices. All six of
its nontriangles are quadrilaterals because their neighbors are D.
No two share an edge by the same argument. There are 24
quad–vertex incidences and at most two per four-valent vertex, so
each vertex alternates two quadrilaterals and two triangles.

The graph with six quadrilaterals as vertices and original vertices
as edges is simple and 4-regular. Two facets sharing two vertices
would share an edge, so parallel edges are excluded. The complement
is a perfect matching. The medial embedding is supplied by the
actual alternating cellulation: its eight regions are the original
triangles. They give the eight distinct 3-cycles of the octahedral
graph. Consequently the whole face lattice is cuboctahedral; the
argument does not infer a geometric polyhedron merely from an
abstract graph.

The explicit cuboctahedral model has square normals of square one
and triangle normals of square (1/2). At each cusp their opposite
RAW sums agree; unit normalization multiplies the triangle sum by
\(\sqrt2\). All twelve rectangles are therefore nonsquare. The
injection requires (1+6+12=19) positive facets, contradicting seven.
This completes the independent conclusion: exactly eight D facets,
all other facets O. There is no use here of a future decomposition.

## 5. The restricted surface argument is geometric and closed

Now assume the O/D spectrum without any volume assumption. All ambient
box face aspects are one or two. Three distinct side lengths would
force a ratio four, so a noncube is (2:1:1) or (2:2:1).
The latter would force each square end to have four quadrilateral
sides, whereas O has none and a square cusp of D has exactly one.
Thus only (2:1:1) occurs. Its four nonsquare D faces form one belt,
all belt adjacencies are quadrilaterals, and the ends are O.

This also ensures that two D facets sharing a quadrilateral use the
SAME diagonal: its nonsquare cusp vertices are exactly the ends of
the long direction. Different intrinsic triangle interiors lie in
different actual facet interiors. An internal diagonal lies in the
relative interior of its actual ridge, with exactly two incident
triangles. At a cusp the four transverse middle segments close into
one simple square link. Compactification gives a disk, not branching
or an open end. Components cannot be joined accidentally by an ideal
point, because all four belt facets are already graph-connected.

For one closed component, (e=3n/2), (v=3n/4), and (\chi=n/4).
The primal/dual tree proof of (\chi\le2) is valid also for an
unoriented surface: the boundary of any chosen face set has even
degree at every circular vertex link, so it cannot be a nonempty
subgraph of a tree. A disjoint dual spanning tree then gives the
inequality. It follows that (n\le8) and (4\mid n).
The original D and its three quadrilateral neighbors and three
distinct cusp opposites are all in this SAME belt component. Hence
\(n\ge7\), forcing (n=8,e=12,v=6). No assertion that the surface
was already a sphere was needed.

## 6. Actual hyperplane and complete cut

The six-cusp graph has no parallel edges: two quadrilateral ridges
with the same two cusps would involve at least three belt facets at
one cusp, hence an opposite asymptotic pair containing both cusps.
That is impossible. Simplicity and degree four identify the
octahedral graph. Two actual facets cannot give the same triangle,
whose finite interior would lie in both facet interiors. Thus its
eight 3-cycles are exactly the eight distinguished triangles, with
cube dual.

Equal short widths give equal OPPOSITE UNIT-normal sums on each belt.
The upper-half-space formula makes both scale and sign explicit:
the sum is the positive separation times the fixed cusp null vector.
The six cube equations kill every mixed sign character, so

\[
n_s=a+s_1b_1+s_2b_2+s_3b_3.
\]

Expansion of the eight unit equations gives (a\perp b_i),
\(b_i\perp b_j\), and (a^2+\sum b_i^2=1). Each cube-edge
orthogonality gives (1-2b_i^2=0). Thus (a^2=-1/2) and
\(b_i^2=1/2). These are EXISTING wall vectors, not a claimed
realization from an arbitrary Gram matrix. Their span U is a genuine
Lorentz four-space, with spacelike orthogonal complement defining H.

The six rays (a\pm b_i) are exactly the actual cusp rays, by the
positive opposite sums. Their convex hull is the regular ideal
octahedron S in H. Convexity gives (S\subset P); the eight actual
wall inequalities give (P\cap H\subset S). Therefore (P\cap H=S).
Its triangles are genuine internal D triangles and its edge interiors
are ridge DIAGONALS, not edges of P.

An outside supporting functional is nonpositive on S. Vanishing at a
finite relative-interior point of any face forces vanishing on that
face. This is excluded for triangular interiors by facet interiors,
for edge interiors by the ridge's exactly two facets, and for S's
interior by the same face-interior argument. Hence every outside
facet misses S at finite points and cannot cross H. This checks
EVERY original wall, including outside doubled facets.

The center of S is strict against all original walls, so both halves
have interior. A new vertex on H must be a vertex of S; otherwise a
nontrivial section segment disproves extremality. Vertices away from
H are old vertices. All vertices are therefore ideal. The new facet
meets precisely the eight cut walls, orthogonal to H. The halves are
right-angled and finite-volume by containment. Exactly eight D facets
become O; no other facet is cut and the new facet is O.

The number of D facets is a strictly decreasing nonnegative integer,
which proves termination. An outside belt component lies on one
side of an earlier H. Its section lies there by convexity and misses
the first section's boundary at finite points. Two sections cannot
share a cusp because a noncubic cusp has only one D belt component.
These facts prove disjointness INCLUDING ideal endpoints.

## 7. Recognition, opposition and absence of twists

For an all-octahedral piece, place one facet on a model 24-cell facet.
A shared ideal triangle and orthogonality uniquely determine its
neighboring supporting wall. In that H3 wall the regular octahedron
with that face on a specified side is unique. Every ideal triangle
symmetry extends to an octahedral symmetry; after canceling it an H3
isometry fixes the three ideal points, so is the identity or reflection
in their plane. The specified side selects the identity.

Propagating through the connected model facet graph matches all 24
whole facets, not just wall planes. Their halfspaces put the piece
inside the model; their ideal vertices and convexity give the reverse
containment. This proves recognition without a four-dimensional
classification or an assumed universal minimum realization.

The disjoint cuts give a tree of pieces. In the model 24-cell, two
cut facets with disjoint closures have spatial root product (-1)
or (-2): product one gives a triangular ridge, zero a common cusp.
If it is (-1), the root sum is a D4 root orthogonal in Lorentz
space to both cut walls. Its facet and its two reflected neighbors
lie in one original supporting wall. None can itself be cut because
it meets an existing cut. Their interiors are disjoint: a shared
open three-dimensional patch would give overlapping inward collars
and hence overlapping four-dimensional piece interiors.

That original facet contains at least three octahedral volumes,
contradicting its permitted O/D intrinsic volume. Thus two cuts in
one cell are opposite. There cannot be three pairwise opposite
facets, so the tree has degree at most two and is a path.

Full flag transitivity supplies the final metric point: every
octahedral facet isometry extends to a 24-cell symmetry. Once that
action is canceled, an ambient isometry fixes the full gluing
hyperplane pointwise and is identity or reflection. The side picks
one. The neighbor is exactly the reflection, and an intermediate
next cut is the unique opposite facet. There is no independent
twist parameter.

At the second volume N is eight by Section 4, so there is exactly one
cut and two cells. The argument therefore proves the UNIQUE reflection
double, not merely two cells or their abstract adjacency. The general
chain theorem remains restricted to O/D, without an existence assertion
for each prescribed length.

## 8. Attainment and finite model counts

The reflection-union calculation has the correct sign:

\[
\langle x,s_n(m)\rangle
 =\langle x,m\rangle-2\langle m,n\rangle\langle x,n\rangle\le0
\]

for (x\in C), because both factors in the subtracted product are
nonpositive. Conversely each side of the gluing hyperplane supplies
all inequalities of the corresponding original cell. Thus the new
intersection is exactly the convex union, not an incompletely tested
Gram chamber. Its two interiors are disjoint and volume is additive.

The two finite ideal vertex sets have six common vertices and 42
distinct vertices in total. Their convex union is their convex hull;
strict convexity of the ideal sphere makes each distinct ideal point
extreme. Hence the new vertex set is complete. The facet-intersection
algorithm described next supplies the recorded face lattice. Independently,
the proven gluing counts give f-vector ((42,168,164,38)).

For m cells, no extra nonadjacent cusp identification can occur:
such a common ideal point would lie on each separating cut along
the path, contradicting closure-disjoint sections. The counts

\[
(v,e,r,f)=(18m+6,72m+24,68m+28,14m+10)
\]

and (8(m-1)) doubled facets follow. The formula and volume additivity
agree with the independent defect identity. None extends the
24-cell-volume divisibility to arbitrary ideal RA4 polytopes.

## 9. Algorithm correctness, halting and replay scope

The double producer and independent checker use exact `Fraction`
linear algebra. The producer uses forward pivots, while the checker
uses reverse columns. Each solves all 10,626 four-wall systems of
the original bounded 24-wall intersection. Its origin is strict,
and the inequalities bound every coordinate. At an extreme point,
active wall coefficients must span R4; otherwise a sufficiently small
perturbation in both directions preserves every inactive strict
inequality. Thus some enumerated independent four-wall subsystem
finds every vertex. Conversely a feasible full-rank solution is extreme.
Singular systems are correctly discarded.

Face reconstruction closes the full vertex set under intersection with
facet vertex sets. Each generated nonempty set is a face, and every
face is an intersection of containing facets. With at most (2^{42})
subsets and monotone addition of new sets, the loop halts and is
complete for this MODEL. Affine rank over Q gives its dimension.
The saved checker repeats the construction independently, checks all
ridge/edge/cusp incidences, and pins the complete literal D4 wall set
and gluing wall. The seven antiprism rational-function identities have
nonzero denominators for (c=\cos(\pi/k)), (k\ge3).

The equality producer covers these specifically bounded domains:

* The fixed ten-vertex planar rotation system and its twelve face
  cycles, producing the four-quad matching. The saved checker checks
  its edges, two-sided incidence, embedding consistency and matching.
* All 364 active triples in the fixed rationally rescaled cuboctahedral
  model. The active-wall completeness argument applies in dimension
  three. Its boundedness follows from the coordinate square walls.
  The independent checker guards the entire wall set, vertices and
  every cusp's squared width ratio two.
* All triples of the fifteen possible missing edges on six labels,
  retaining exactly perfect matchings. This is exactly the class of
  simple 4-regular six-vertex graphs. Every complement's eight triangles
  and dual cube labels are checked.
* The six rational cube-face equations, their rank four, the four
  independent constant/linear modes and their annihilation. This
  proves completeness of THAT linear solution space.

The bounded defect-partition recursion generates the four partitions
used in the prose; its length increases at every recursive call.
All loops and rank computations terminate over finite inputs. The
standalone helper imports no producer. It performs its mathematical
checks at module level, so it must NOT be casually imported during
source inspection or planning; I did not import it.

The D4 checker covers all 276 unordered root pairs and their explicit
facet vertex intersections. Its 96 root-sum cases give 1,152 reflected
octahedral vertex incidences. This supports the finite root identities,
not the universal disjoint-cut or inward-collar proof. Every producer
and checker explicitly keeps historical classification replay and
universal geometric proof outside its claimed scope. No finite census
of four-polytopes is hidden in the algorithms.

## 10. Exact external-incidence and execution limitations

The ten-vertex transcription has quadrilateral cycles

    ABCD, ABFE, EHIJ, CGHI.

Their shared EDGES are exactly AB and HI. The other shared vertices
C and E do not create edge adjacency. I checked this implication
directly in the producer, saved checker and stored cycle definitions.
The original text explicitly calls Figures 6–7 the complete admissible
list. A fresh image fetch of Figure 7 failed; identifying this precise
stored drawing with that external figure remains an EXTERNAL
TRANSCRIPTION premise in this audit. Earlier recorded visual reviews
retain their own dates and authorship. I do not claim a new primary
image check or historical census replay.

The exact saved model files remain:

|File|SHA256|
|---|---|
|`results/fresh_ideal_ra4_double24_volume_20261018.json`|`386d9650b8f669e3edcb3cdee56ed53d1c9ea4a90af2b4ba5c8c1b83f9d7bfa3`|
|`results/fresh_ideal_ra4_equality_structure_20261018.json`|`589a7901494481ddbefb8b3cf2a51854c70bd804bd78128386bd650dcf4e19fc`|

Their arithmetic contents were not regenerated. The inherited core
receipt is `results/ideal_ra4_unified_core_verification_final_20261018.json`,
SHA256 `5763d85328eb435254e24f286411ca18dd692f4baf139ca47ce7603ab511243b`.
It binds the 16-page predecessor, not the current publication TeX.
The 17-page publication has the separate successful build record
`results/ideal_ra4_volume_and_structure_publication_build_20261018.json`,
SHA256 `6787fbb05a4c593b76da05ed55b10c711675bccf17d584b6b53183d8b6ba5fd8`.
That PDF build does not create a fresh combined mathematical replay.

## 11. What the conditional Lean modules verify

`IdealRA4Gap` checks the exact integer case partition and its budget
consequence (12\le W), then the volume numerator implication. Its
`FacetCase` is the EXTERNAL geometric/small-facet classification
interface, rather than a formal classification of hyperbolic facets.
`CertificateBounds` proves the selected-list lower sum by induction;
the disjoint actual selected/rest partition must be supplied.

`ManuscriptIdeal4` verifies the eight-facet arithmetic from named
equality stages, the strong recursion on the D-facet count, the number
of cells and volume numerators, and the final conditional unique-double
implication. The geometric source supplies minimum-three and the N7
reduction before its granular equality interface is constructed.

`OctDGeometry.positiveCut` contains the actual cut and decrease by eight;
`zeroRecognition`, `pathRecognition` and `canonicalTwoCellGluing` remain
WRITTEN geometric premises. Their internal proofs were audited above
but are not formalized by accepting those records. `EqualitySpectrum`
can be constructed from the more granular written-stage and model
interfaces; it is not being passed off as a published theorem.
The volume numerator is in units (\pi^2/3). The larger-unit
conversion is supplied explicitly and is not assumed universally.

No concrete real polytope instantiates these records here. Printed
axiom lists do not list ordinary parameters, and successful kernel
audits do not verify their geometric truth. The actual earlier Lean
execution boundary is recorded in the inherited core receipt; this
audit adds no compilation, axiom-audit or end-to-end geometric credit.

## 12. Consolidation recommendation and verdict

Keep this as ONE manuscript. Its current order correctly separates
the universal gap, independent equality-spectrum proof, general
restricted chain theorem, and final equality/attainment application.
The difficult structural steps are written out rather than relegated
to computational certificates. Five figures depict the exact cusp,
surface, graph and cut structures needed by the proof.

No mathematical manuscript repair is required by this audit at the
stated boundaries. For the final reproduction package, keep the
published small-facet/classification input and the Figure 7 incidence
transcription explicit, preserve the source-bound predecessor records,
and distinguish the publication build from a combined runtime pass.
The accepted supporting algorithms already have finite halting and
completeness arguments in the publication appendix. Exploratory
antiprism-spectrum and larger-facet extensions should stay outside
this manuscript unless they supply a separately chosen new theorem.

The conclusions are the internally written universal second-volume
and unique-double theorem, and the SEPARATE restricted O/D chain
theorem. This verdict does not claim universal volume divisibility,
arbitrary ideal-polytope chain decomposition, arithmeticity, external
refereeing, complete geometric Lean verification or a new execution.
Root acceptance is a separate coordinating step.
