# Thirty vertices: reductions for an entire regular pentagon

20 October 2026. Paper reasoning only. These are supporting reductions,
not a proof that one regular pentagon forces the dodecahedral double.
No search, worker, build, or accepted-source change was made.

Let Q be a compact right-angled H3 polyhedron with V=30, hence F=17.
Put phi=(1+sqrt(5))/2 and L=arcosh(phi). A regular pentagon below is
an ENTIRE actual regular right-angled pentagonal face of side L.

The inputs are the actual-neighborhood/capacity, three-L closure, and
Andreev metric-recognition conventions in
`manuscripts/compact_ra4_volume_constraints.tex`, SHA256
`afe807672f3a6bbf4bb2dd95d2861b60391d7c2d7c879c6d6382e749e7ad976f`,
and the full-flower theorem
`results/compact_ra4_uniform_dodecahedron_neighbor_block_20261020.md`,
SHA256 `69d900288cd778e6cda38cb6eedb932ef0583eb363649f60727be41dd98b9da3`.
The dual sphere is flag with no induced four-cycle. A K4 is impossible.
Every face has at least five sides and capacity gives max face size7.
If p,h,s count pentagons, hexagons, heptagons, then

    p=12+s,       h+2s=5,       s=0,1,2.                 (1)

## 1. A regular pentagon cannot have exactly one nonpentagonal neighbor

Let v_i be its cyclic vertices, n its outward normal, and s_i>0 its
outgoing lengths. The five actual opposite normals are distinct and

    m_i=sinh(s_i)v_i-cosh(s_i)n,    u_i=coth(s_i)>1.

Adjacent v-products are -phi and nonadjacent ones -phi^2. Thus

    m_i.m_j=(u_i*u_j-cosh dist(v_i,v_j))
                       /sqrt((u_i^2-1)(u_j^2-1)).        (2)

If the face between m_i,m_(i+1) is pentagonal, those two sides meet and
their product is0, so u_i*u_(i+1)=phi. If that face is nonpentagonal,
the two actual sides are nonadjacent and their product is <=-1. An
orthogonal intersection would be the actual vertex closing the five
already fixed sides of that face, making it pentagonal instead.

Suppose just the neighbor between m_4,m_0 is nonpentagonal. The four
zero products force u_0=u_2=u_4=t and u_1=u_3=phi/t. Equation(2) for
m_4,m_0 gives

    (t^2-phi)/(t^2-1)<=-1,       t^2<=phi^2/2.

For the distinct nonconsecutive m_1,m_3, equation(2) is strictly
negative and equals phi^2*(1-t^2)/(phi^2-t^2). The right-angled native
gap therefore forces

    t^2>=2phi^2/(phi^2+1)>phi^2/2,

the contradiction, since phi^2<3.

If there are ZERO nonpentagonal neighbors, the same five zero-product
equations form an odd cycle and force every u_i=sqrt(phi), hence every
s_i=L. All five neighboring pentagons are regular: each has its base
L-edge and two adjacent L-edges. The complete flower theorem then
identifies Q with the reflected double of D at V30.

Consequently a non-double candidate containing a regular pentagon must
give that pentagon at least TWO nonpentagonal neighbors. This is not
an exclusion: the genuine double also has regular pentagons incident
to two hexagons, so one may not assume the originally selected face
has a complete flower.

## 2. The entire heptagon branch has one combinatorial form

Take a heptagonal face A. The capacity neighborhood consists of A,
its seven cyclic neighbors A_i, and the seven distinct outgoing
opposites B_i. Exactly two extra faces C,D remain. Each extra meets
at most one A_i.

In the complementary B-disc, every extra meeting an A_i is inserted
on the boundary between B_(i-1),B_i. If t extras do so, its boundary
has 7+t vertices and its interior has 2-t vertices. Boundary B's have
only two outside A-neighbors and thus disc degree at least3; an inserted
extra has only one outside A-neighbor and thus disc degree at least4.
For t=2 the ear contradiction excludes a triangulation with no interior
vertices. For t=1, a boundary chord cuts off a side without the sole
interior vertex; the ear away from its endpoints is again impossible.
There are therefore no chords, so the interior vertex joins every
one of the eight boundary vertices. Its degree8 contradicts maxface7.
Thus t=0, and all seven A_i are pentagons.

Both C,D now lie inside the B7-disc. If they were nonadjacent, each
would have at least five B-neighbors. Their common B-neighbors must
form a clique: a nonadjacent pair would give an induced four-cycle
with C,D. That clique has size at most2, or it makes a K4 with C.
Their union would consequently contain at least8 of the seven B's.
Hence C,D are adjacent.

Every B-chord has C,D on the same side, since their edge cannot cross
it; an innermost empty side gives the forbidden ear. Thus there are
no B-chords. A triangulated disc with seven boundary and two interior
vertices has17 edges. Removing its seven boundary edges and CD leaves
nine interior-to-boundary incidences. The two common C,D-neighbors
are exactly the third vertices of the two triangles on CD. Each
interior vertex needs at least four B-neighbors. Their consecutive
blocks therefore have sizes4 and5, overlapping at those two endpoints.

This uniquely describes the graph up to exchanging C,D and reversing
the B-cycle. The two junction B's and the degree6 cap are hexagons;
the degree5 cap and other B's are pentagons. There is exactly ONE
heptagon and three hexagons. In particular the two-heptagon case in(1)
cannot occur. No metric assertion about regular pentagons in this
split-cap graph has been established here.

## 3. All pentagons regular is already sufficient

This additional hypothesis forces the double without a fullerene
enumeration. In the induced nonpentagonal-face graph H, a hexagon
needs at least two neighbors: otherwise it has three consecutive
L-edges and closes as a pentagon. When it has exactly two, they must
be opposite in its six-cycle. In particular a degree2 hexagon cannot
belong to a triangle of H. A heptagon needs at least three neighbors.

For s=2, H has three vertices and cannot supply degree3 to a heptagon.
For s=1, H has four vertices. The heptagon joins all three hexagons;
the latter need at least one edge each among themselves. A triangle
would make K4; otherwise they form a path, whose end hexagon has
degree2 in a triangle, impossible.

For s=0, H has five vertices of minimum degree2. It has no triangle:
given a triangle, call the other vertices x,y. If xy is absent, each
has either degree2 inside a triangle or all three triangle neighbors
(a K4). If xy is present, each has one or two triangle neighbors.
An attachment of size1 shared with the other's attachment makes that
vertex degree2 in a triangle; disjoint attachments produce an induced
four-cycle; equal size2 attachments make K4; different size2 attachments
again make an induced four-cycle. These exhaust the possibilities.
The triangle-free graph must therefore be C5: any shorter cycle would
be an induced four-cycle, and a chord of the five-cycle makes a triangle.

The H-cycle separates the dual sphere into two pentagonal caps. The two
H-neighbors of each hexagon are opposite, so each side has two pentagon
neighbors at that hexagon. The five third vertices U_i of successive
H-edges are distinct on each side: a repeated U_i would give a common
pentagon adjacent to both H-neighbors of some hexagon, contradicting
their opposite positions. They form a five-cycle. Inside it every U_i
needs degree at least3, having only two outside H-neighbors. The ear
argument requires at least one interior pentagon per cap. Both caps
therefore use at least six of the twelve pentagons; each uses exactly
six. A chord would again cut off an empty ear-side, so the sole interior
pentagon joins all five U_i. It and its five regular neighbors form
the complete flower. The cited theorem identifies Q with doubleD.

## 4. Remaining obstruction and global scope

What remains unproved is propagation from ONE regular pentagon to all
pentagons or to SOME complete flower. The heptagon reduction does not
exclude accidental regularity in its one explicit graph. In the
pentagon/hexagon branch, the three C30 fullerenes are the external count
in Brinkmann--Goedgebeur--McKay, Table2, printed p15, row(30,17):
https://arxiv.org/pdf/1207.7010v2 . That count alone proves no metric
rigidity, and this note does not classify their regular faces.

The already accepted d>=5 boundary charge gives, with the notation of
the frozen global-charge note,

    Nu>=f-(26/15)W+(13/15)TI-m+(2/3)WIB+E5/5.

Since TI>=2m and f>=65, W<=37 forces Nu>=13/15>0. At or below the
120-cell volume, a connected D-subcomplex in a non120-cell therefore
requires W>=38 and f+W>=103. These statements do not give f+W>=120,
nor do the reductions above identify every boundary facet of excess5.
