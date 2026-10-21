# Independent full-source review: compact RA4 excess at least seventeen

19 October 2026. **PASS, relative to the external geometric inputs stated
in the manuscript.** This was a complete written-source and mathematical
review, with no scientific calculation, graph census, checker, Lean job,
TeX build or PDF inspection. No reviewed or frozen source was changed.

Reviewed all 1188 lines of
`manuscripts/compact_ra4_volume_constraints.tex`, SHA256
`afe807672f3a6bbf4bb2dd95d2861b60391d7c2d7c879c6d6382e749e7ad976f`,
and its current preparation/proof map
`docs/compact_ra4_volume_constraints_preparation_20261019.md`.
The previous full review of source20e and the independent W=16
endpoint-exclusion review remain separate historical records. This
review binds the present W>=17 source, not a successor PDF.

## Inherited geometry, capacity and volume sections

Sections1--6 retain the complete arguments needed later. The two-normal
positive-definite projection adds nonnegative outward-normal multiples
to a strict interior point; every other wall remains strict. The
null-sum case excludes product -1 by compactness, so the compact gap is
0 or strictly below -1. Orthogonal-wall projection proves actual
intersections, rather than using abstract adjacency alone. The
perpendicular-edge opposite has the correct orientation and formula;
the product of two proposed walls on the same perpendicular is positive,
which proves its uniqueness as an actual facet.

The polygon closing lemma includes both the no-extra and one-extra
cases. Its complementary discs and ear count use the exact two outside
neighbors of each B_i; a single extra meeting one A_i instead has only
one outside neighbor. The following no-chord argument identifies the
Loebell equality lattice. The four-dimensional capacity proof correctly
uses singleton footprints when e_i=1 and at most two incidences for
every remaining extra facet. It does not treat the lower bounds e_i>=1
as disjoint contributions without that additional argument.

I rechecked W=v-5f, d(F)=(V-20)/2=M-12, sum d=2W, and the orbifold
Euler/volume identity. The neighbor excess proof sums distinct facets;
shared ridge contributions are legitimate in that exact sum. The
cutoffs46 and96 use the correct thresholds120 and256, evenness of V,
and the explicit endpoint exclusions48 and98. The double's counts
1160 vertices and226 facets give W=30 and volume68*pi^2/3.

Zero-excess recognition assumes every facet is the metric regular
dodecahedron, using the expressly named three-dimensional classification.
Whole-facet propagation fixes both supporting hyperplanes and their
actual regular dodecahedra; connectedness supplies the two inclusions
of facet sets. It does not propagate from one arbitrary regular facet.

## L6 metric and the new degree-two bound

Section7's Loebell symmetry has the correct p_i/q_i indices. Compact
three-dimensional metric uniqueness realizes its reflection, producing
the cyclic lateral-pentagon pattern (b,a,c,c,a). The regular-base
formula and the stated Lorentz proof of the pentagon identity give
cosh(b)=2, cosh(a)=sqrt(2), cosh(c)=sqrt(3). Thus b, the base edge,
is uniquely longest in the actual shared pentagon.

The proposition `prop:degree-two-twenty` is universal, not restricted
to W=16. For k>=7 the earlier estimate already gives W>=20. For
k=6, its fifteen selected facets are distinct and contribute the
baseline30 to sum d. A lateral neighbor of excess2 has two unit
G-edges and is another L6. The unique-longest edge in the shared
pentagon forces its hexagonal-base parent to be the original B_j:
the actual edge has precisely three parent facets by compact
simplicity, and its parent F intersects the lateral neighbor
pentagonally. Hence there are 6-m_j distinct additional unit edges
at B_j, besides F B_j. Its additional excess is at least
max(0,5-m_j); together with the m_j high lateral facets this is
at least five for each row. The two rows/base sets are disjoint,
so sum d>=40 and W>=20. No ridge-disjointness assumption is used.

The three-L-edge polygon lemma uses actual closing walls and their
orthogonal projection; an uncounted wall cannot truncate their
intersection. Consequently eligibility is proved for actual third
facet carriers. The G/H distinction is preserved throughout:
pentagonal positive adjacencies remain in H when excluding induced
squares. The four-node link enumeration and minimum H-degree5
are correct whenever minimum G-degree is at least3.

## Complete small-graph and marked-polygon arguments

Section8's subcubic nine-node lemma includes triangle, bipartite and
shortest-odd-cycle cases. In the tight triangle case all three
outside slots are used as a matching; its three zero-contact vertices
are indeed isolated. In the C5 case every outside vertex has zero
or two cycle contacts, and zero-contact vertices are isolated by
the induced2K2 prohibition. The constructed independent fives are
valid. The reductions for N<=8 cover each possible nontrivial
component. Thus `prop:ten-positive` proves N>=10 as a reusable
statement, without borrowing a previous global W>=15 result.

Section9's pentagon and hexagon lemmas have the correct metric scope.
The two-L pentagon uses a unique common perpendicular. For the
four-L hexagon, the exceptional positions must be opposite, and
the two explicit Lorentz equations give cosh(a)=cosh(d)=phi^3,
so both lengths are 2L.

The five-neighbor exclusion is local to an excess3 facet, under
minimum G-degree3. Its universal-vertex, triangle with/without
the outside edge, and trianglefree forest/C5 cases exhaust the
five-node link. The exclusions of marked degree-two vertices with
adjacent neighbors and degree-three vertices with a neighbor-path
come from actual polygon closure. The final P5/C5 cases use a
genuine pentagonal consumer forced regular by its L-edges, whose
shared 2L edge is impossible. The proof does not silently assume
the other facets have excess3. The example C5[K2] appropriately
warns that the weaker abstract ten-node hypotheses alone do not
suffice.

I checked all cases of `lem:ten-complement`. With two nonadjacent
degree-four vertices u,v, the partition X,A,B,S has sizes
t,4-t,4-t,t, and all other vertices are subcubic. Every A--B
contact is forced by the induced2K2 restriction. This leaves
t=2,3,4. At t=2 or3, an S-edge would have only four spare endpoint
degrees to cover respectively six or five vertices in N(u) union
N(v), so S is independent. The displayed independent fives then
follow. At t=4 an S-edge spends all spare degrees on X and cannot
have another S-edge adjacent to it; a disjoint second edge would
itself induce2K2. Therefore S has an independent triple, to which
u,v can be added. The deletion reductions for zero/one high-degree
vertex or two adjacent ones are also valid.

## Final composition, figures and limits

The final proof is a single nonstaged argument. For 0<W<=16 the
degree-two>=20 proposition gives minimum G-degree3, so 3N<=2E<=32.
The reusable small-graph theorem gives N>=10, hence N=10. Since
sum d<=32 and each d>=3, at most two facets have excess greater
than3. Every other facet has H-degree at least6 by the local
five-neighbor exclusion; the exceptional facets still have
H-degree at least5. The complement therefore satisfies exactly
the ten-node lemma, producing the forbidden K5 in H. This proves
W=0 or W>=17 without an omitted W15/W16 case or circular dependency.

All four vector figures were read against their arguments. Figure4
retains the U--V edge and shows the correct independent five; its
lower explanatory nodes are moved to -2.95/-3.1 while the s_i
labels remain at -1.65/-1.8. Thus the requested downward spacing
is present in source. No claim about rendered spacing or a new
PDF is made. The chain figure correctly depicts three fixed
L-edges and actual closing walls, and captions distinguish
combinatorial drawings from hyperbolic metric models.

The source explicitly retains the acute-wall background, compact
RA3 twelve-face metric classification, compact Andreev uniqueness,
standard convex-polytope facts and orbifold Gauss--Bonnet as
external inputs. Kellerhals--Perren's opposite-facet mechanism
and Inoue's excess identity receive antecedent credit, with the
correct Proposition3.2 reference. No correction is requested.
There is no end-to-end geometric Lean verification, universal
minimum-volume or second-volume theorem, arbitrary120-cell
decomposition, new example or dimension-bound claim in this PASS.
