# Compact RA4: excluding excess sixteen

19 October 2026. **New written candidate; awaiting independent review.**
No calculation, finite graph census, optimizer, Lean or TeX build ran.
The primary manuscript at SHA256
`20e5b4a3aff508520b7db3c4ae1a8904c0aa2d3299b4164a657bb6432433c49e`
is unchanged. The proposed successor is W=0 or W>=17; this is not a
compact four-dimensional minimum-volume or second-volume theorem.

## 1. Exact inherited interfaces

Use d(F), W, G and H from the reviewed compact RA4 proofs. G has an
edge of weight k-5 for each nonpentagonal ridge, while H is FULL
adjacency on the positive-excess facets. Thus weighted G-degree is
d(F), and sum d(F)=2W. The source packets are:

* `compact_ra4_positive_facet_link_gap_fifteen_20261019.md`, SHA256
  `8201b2f5cc7da5988b6d164b258da8e14981cf94c6b5e980a5a87e9deb70dcd7`;
* `compact_ra4_fifteen_excess_equality_exclusion_research_20261019.md`,
  SHA256
  `237d9ccc7c91062b268124fc78c4fa82dbcf4bbf6b0b6b32c8eae013c713ff87`.

Their geometric inputs include compact 3D metric uniqueness and the
12-face regular-dodecahedron classification. The following proved
interfaces, rather than any four-dimensional classification, are used:

1. Each positive facet has d>=2; a G-edge of weight w forces d>=2w
   at both ends. Thus G has no degree-one vertex.
2. A degree-two G vertex has equal weights w, k=w+5, and its facet
   is a Loebell L_k: two k-gon bases and 2k lateral pentagons. Its
   distinct base neighbors have d>=2(k-5), and its 2k distinct
   lateral neighbors have d>=2. Including itself gives W>=5k-15.
3. If G has minimum degree three, H has minimum degree at least five.
   H has no K5 and no induced C4. The resulting graph conditions
   exclude nine or fewer positive facets.
4. Every induced-2K2-free subcubic graph on nine vertices has an
   independent five. This is the explicit triangle/bipartite/C5 proof,
   not a graph census.
5. W=0 or W>=16; zero identifies the regular right-angled 120-cell.

We also use the following local conclusion extracted from the full
five-neighbor proof in the second packet:

**Local link exclusion.** A positive facet of excess three cannot have
exactly five positive neighbors, whenever G has minimum degree three.

To check its scope, such a facet has three hexagons and twelve
pentagons. Its five positive-neighbor faces give a graph J with three
marked hexagons; all remaining faces border zero-excess dodecahedra.
The proof uses only that local data, the no-K4/no-induced-C4 condition
on J, and the shared L- and 2L-edge lengths. Universal-vertex, triangle,
forest and C5 cases exhaust J. No part assumes that its neighbors or
all other positive facets have excess three. In particular this local
conclusion is not restricted to the global equality W=15.

Here L is the regular right-angled pentagon edge length. The relevant
polygon facts, proved in that packet, are: three consecutive L-edges
force the regular pentagon by projection onto its two actual closing
walls; a right-angled pentagon with two adjacent L-edges is regular;
and a right-angled hexagon with four L-edges has its other two sides
opposite and both of length 2L.

## 2. The lateral pentagon of L6

**Metric lemma.** Each lateral pentagon of the compact right-angled
Loebell L6 has cyclic lengths

    (b,a,c,c,a),   cosh b=2,   cosh a=sqrt(2),   cosh c=sqrt(3).

In particular its edge on its hexagonal base is its UNIQUE longest
edge. This is a marked metric assertion, not just the face lattice.

Write the upper/lower bases U,D and their respective lateral rows as
p_i,q_i, indices modulo six. The incidences are p_i adjacent to
U,p_(i-1),p_(i+1),q_(i-1),q_i, and q_i adjacent to
D,q_(i-1),q_(i+1),p_i,p_(i+1). The combinatorial reflection

    p_i -> p_(-i),     q_i -> q_(-i-1),     U -> U, D -> D

is realized isometrically by compact 3D uniqueness. It fixes p_0,
exchanges its two neighboring p faces and exchanges its two q faces.
Consequently p_0 has lengths (b,a,c,c,a), with b on U. The same holds
for every lateral face by rotation and exchange of the two rows.

The regular-base formula already proved in the inherited packet gives
cosh b=1+2cos(2pi/6)=2. For completeness the right-angled pentagon
identity for cyclic lengths s_0,...,s_4 is

    cosh(s_i)=sinh(s_(i+2))*sinh(s_(i+3)).

One direct Lorentz derivation chooses the common endpoint of sides
s_(i+2),s_(i+3) as (1,0,0), with their outgoing directions the two
coordinate axes. The outward normals of the sides beyond their far
ends are (sinh u,cosh u,0) and (sinh v,0,cosh v). Their product is
-sinh u sinh v. The connecting perpendicular is side s_i, so the
same product is -cosh(s_i), proving the identity.

Applied to the displayed pentagon it gives

    2=sinh(c)^2,        cosh(a)=sinh(c)*sinh(a).

Therefore cosh c=sqrt(3) and coth a=sqrt(2), hence cosh a=sqrt(2).
All lengths are positive, and 2>sqrt(3)>sqrt(2). The longest edge
is therefore exactly b, without a second equal-length alternative.

## 3. All degree-two cases at W=16 are impossible

Assume W=16 and let F have G-degree two. The inherited estimate
5k-15<=16 forces k=6. The original F has d=2. Its two base
neighbors B_1,B_2 and twelve lateral neighbors A_i are distinct and
all have excess at least two. The sum over these fifteen facets is
therefore at least thirty, while the complete sum is thirty-two.
Thus every selected facet has d<=4, and at most two of the twelve
lateral neighbors can have d>2. This remains true if there is an
additional positive facet elsewhere: it consumes the same two units
of remaining excess, rather than invalidating the estimate.

Each base B_j meets six of these lateral neighbors in F. At least
four of the six have d(A_i)=2. Each such A_i is itself L6 by the
degree-two equality interface: d=2 and d>=2w imply precisely two
unit-weight G-edges. Its shared pentagon R=F intersect A_i is a
lateral pentagon in BOTH Loebell facets. Seen from F, the unique
longest edge of R lies in B_j. Seen from A_i, the same unique
longest edge lies in one of A_i's hexagonal bases. These are the
same geometric edge. Compact simplicity gives exactly three parent
facets at its interior: F,A_i,B_j. The parent of A_i's hexagonal
base is consequently B_j, since F intersects A_i pentagonally.
Thus A_i intersect B_j is a hexagon, giving a unit G-edge A_i--B_j.

Together with the original F--B_j hexagon, the four such distinct
A_i force d(B_j)>=5. But d(B_j)<=4 by the exact thirty-two-unit
sum above. This contradiction excludes EVERY degree-two case at
W=16. It does not assume that all lateral neighbors are L6, or that
their shared pentagons are regular.

## 4. A ten-vertex complement lemma

**Graph lemma.** If C has ten vertices, has no induced 2K2, has maximum
degree at most four, and at most two vertices have degree four, then
C has an independent set of size five.

If at most one vertex has degree four, delete it (or any vertex when
there is none). The remaining nine-vertex graph is subcubic and has
an independent five by the inherited lemma. If exactly two degree-four
vertices u,v are adjacent, deleting u again leaves a subcubic graph.

It remains to treat two nonadjacent degree-four vertices. Put

    X=N(u) intersect N(v),          t=|X|,
    A=N(u) minus N(v),              |A|=4-t,
    B=N(v) minus N(u),              |B|=4-t,
    S=the remaining vertices,       |S|=t.

Sets X,A,B,S partition the eight vertices other than u,v. Every
other vertex has degree at most three. Every a in A must be adjacent
to every b in B; otherwise u-a and v-b induce 2K2. If A is nonempty,
this gives degree(a)>=1+(4-t), forcing t>=2. Therefore t=2,3,4.

At t=2, A and B each have two vertices, and all are saturated by
their complete A--B contacts and their u- or v-contact. Thus A is
independent and has no S contacts. Also S is independent: an edge
yz in S would have to contact collectively every vertex of
N(u) union N(v), of size six, to avoid 2K2 with one of the u- or
v-edges. But y,z have at most four spare degrees after yz. Now
{v} union A union S is an independent five.

At t=3, an S-edge would similarly have to contact the union
N(u) union N(v), now of size five, with only four spare degrees.
Thus S is independent and {u,v} union S is an independent five.

At t=4, A=B=empty. Any S-edge yz must collectively contact every
one of the four vertices of X. The endpoints' four spare degrees
are therefore all used on X. In particular neither endpoint has
another S-neighbor. A second S-edge would then induce 2K2 with yz.
Hence S has at most one edge and has an independent triple. Add
u,v to obtain an independent five. This exhausts every case.

The covering condition for an S-edge is literal: S has no u- or
v-contacts, so if neither endpoint contacts x in N(u), the four
vertices y,z,u,x induce 2K2; the same applies to N(v). It is not an
assumption about common-neighbor sizes beyond this induced subgraph.

## 5. The remaining W=16 case

By Section 3, G has no degree-two vertex. It therefore has minimum
degree at least three. If N is its number of positive facets and E
its unweighted edge count, then

    3N<=2E<=2W=32,   so N<=10.

Minimum H degree is at least five. The inherited small-graph argument
excludes N<=9, hence N=10. Each facet has d>=3, and the exact sum32
means at most two facets have d>=4. Every other facet has d=3 and
cannot have H-degree five by the local link exclusion in Section 1.
Consequently H has minimum degree at least five, with at most two
vertices of degree five; all its other vertices have degree at least
six. Its complement C therefore has maximum degree four and at most
two vertices of degree four. It has no induced 2K2 because H has no
induced C4. Section 4 gives an independent five in C, equivalently a
forbidden K5 in H. This is the final contradiction.

Thus W=16 is impossible in both cases. Together with the reviewed
W=0 or W>=16 theorem this yields the proposed successor

    W=0 or W>=17.

The zero case is the regular 120-cell. No new numerical cutoff,
minimum-volume assertion, second-volume classification, garland
decomposition, higher-dimensional nonexistence or novelty claim
is inferred. The current primary source and earlier packets were
not changed; review is required before promotion or integration.
