# Compact RA4: a positive-excess gap of fifteen

19 October2026. **Written theorem candidate, awaiting root and independent
review.** No scientific job, finite census, optimizer, Lean or typesetting
ran. All prior source-bound notes and the compact-capacity manuscript are
unchanged. This is not a compact minimum-volume or second-volume theorem.

## Statement and explicit inputs

For a compact right-angled hyperbolic four-polytope P, put

    W=sum_(ridges R)(number_of_sides(R)-5).

**Claim.** W=0 or W>=15.

The reviewed zero-excess recognition makes W=0 exactly the regular
right-angled120-cell. The claim does not establish f+W>=120.

For a facet F, write

    d(F)=sum_(R contained in F)(number_of_sides(R)-5)
        =number_of_faces(F)-12,
    sum_F d(F)=2W.                                  (1)

Use the reviewed own3D capacity theorem: a compact right-angled
three-polyhedron containing a k-gon has at least2k+2 faces. At equality
it has Loebell combinatorics, two disjoint k-gons and2k pentagons.
Thus a positive facet has d>=2. We also use actual-face orthogonal
projection, the strict compact normal gap, and its consequences that
the full facet graph has no induced square and no clique of size five.
The source is
`compact_ra4_neighbor_excess_and_zero_excess_rigidity_20261019.md`,
SHA256 `9c670a7b3c25dda48997c3f70a6d1be9be88bcdefb475b7710f5011068aebcbe`,
with root review SHA256
`292187dbf4fb5ea327baed6086177ec1b0fac174116ca3c579a6d75bc43f0baa`.

Two external3D metric inputs are explicit: zero-excess facets are regular
right-angled dodecahedra; and the compact nonobtuse Andreev uniqueness
clause realizes every marked combinatorial symmetry of a right-angled
Loebell polyhedron by an isometry. For the latter clause see
[Roeder--Hubbard--Dunbar, Theorem1.4, p832](https://aif.centre-mersenne.org/article/AIF_2007__57_3_825_0.pdf#page=9).
No four-dimensional classification or volume-minimum input is used.

Let G have one vertex per positive-excess facet and one edge of weight
k-5 for each shared nonpentagonal ridge. Its weighted degree is d(F)
and its total edge weight is W. Let H be FULL adjacency among the same
positive facets, including pentagonal shared ridges. Thus G is a spanning
subgraph of H. These two graphs must not be conflated.

## 1. Degree two in G already forces W>=15

An edge of weight w forces d(F)>=2w at either endpoint, by3D capacity.
Consequently G has no degree-one vertex. At a degree-two vertex with
weights w_1,w_2 this gives w_1+w_2>=2max(w_1,w_2), so both weights
are w>=1. Put k=w+5. Equality in3D capacity identifies its facet F as
the Loebell type with two k-gonal bases and2k lateral pentagons.

The combinatorial k-fold rotation of this type preserves each base and
rotates its faces in cyclic order. By the stated metric uniqueness it
is a genuine k-fold isometry of F. Its convex base is therefore a
regular right-angled k-gon. The triangle from its center to a vertex
to the adjacent edge midpoint gives its edge length L_k:

    cosh(L_k/2)=sqrt(2)cos(pi/k),
    cosh(L_k)=1+2cos(2pi/k).                         (2)

If L=L_5, then cosh(L)=phi=(1+sqrt(5))/2 and L_k>L for k>=6.
Every lateral pentagon of F contains a whole edge of its base of length
L_k, and hence is not a regular right-angled pentagon. Its neighbor in P
cannot be a zero-excess regular dodecahedron. All2k lateral neighbors
are therefore positive, of excess at least2. They and the two base
neighbors are distinct: two facets of a convex polytope cannot share
two different ridges. Each base neighbor has excess at least2(k-5).
Summing these distinct facets with F in(1) gives

    2W>=2(k-5)+2*2(k-5)+(2k)*2=10k-30,
    W>=5k-15>=15.                                   (3)

Thus if0<W<15, the ordinary degree of every vertex of G is at least
three. This conclusion does not assume any outside facet is sealed
without proving it: nonregular lateral pentagons force all the counted
outside neighbors to be positive.

## 2. Four-edge closure and positive-neighbor eligibility

A compact right-angled polygon with four consecutive edges of length L
is the regular right-angled pentagon. Match their inward-oriented chain
to four edges of that pentagon; fixed lengths and90-degree turns
determine the chain uniquely. At either end the next supporting line
is the perpendicular through the end vertex. Both are the same fifth
pentagon line. Convexity allows only one edge on that supporting line,
which closes the boundary. An additional intervening edge is impossible.

For any nonpentagonal ridge R=F intersect A, an edge carried by a
zero-excess third facet is an entire edge of a regular pentagon, so has
length L. Hence R cannot have four consecutive zero-excess carriers.
If R has at most one positive-excess carrier, at least k-1>=5 edges
are a single zero-carrier block, impossible. If R has exactly two
positive carriers which are adjacent edges, its remaining k-2>=4 edges
are also such a block, impossible.

Fix a positive facet F. In the graph induced on its H-neighbors, call
a vertex A eligible if its degree is at least three, or its degree
is two and its two neighbors are not adjacent. Every G-neighbor A is
eligible. Indeed its degree in this graph is exactly the number of
positive third-facet carriers on R=F intersect A. A common H-neighbor
gives an actual edge of R by orthogonal projection. Two such carriers
which are adjacent in the graph give an actual vertex of R by the
same projection, so they are adjacent polygon edges. The preceding
zero-block argument proves the eligibility requirement.

**Four-vertex graph fact.** A graph on at most four vertices, with no
K4 and no induced C4, has at most two eligible vertices.

Here is a complete proof. If a vertex has degree three, its other three
vertices induce at most two edges, since three give K4. With zero or
one of those edges only the degree-three vertex can be eligible. With
two edges they form a path: its middle and the original vertex have
degree three and are eligible, while each endpoint has two adjacent
neighbors and is not. If no vertex has degree three, components are
paths or cycles. A triangle has no eligible vertex, C4 is excluded,
and paths on at most four vertices have at most two internal vertices.
This also covers fewer than four vertices.

Assume now0<W<15. Every positive F has at least three G-neighbors,
all eligible. If its H-degree were at most four, its neighbor graph
would contradict this fact: it has no K4 (otherwise with F there is
an ambient K5) and no induced C4. Consequently

    minimum_degree(H)>=5.                           (4)

## 3. The remaining small graph obstruction, proved without a census

Let N be the number of positive facets and E the unweighted number of
edges of G. Under0<W<15, its minimum degree three gives

    3N<=2E<=2W<=28,       N<=9.                     (5)

Let C be the complement of H. The absence of K5 in H means that C
has no independent set of size five. The absence of an induced C4 in
H means that C has no induced pair of disjoint edges (an induced2K2).
By(4), maximum_degree(C)<=N-6.

N<=5 is incompatible with(4), and N=6 would make H=K6. If N=7,
C has maximum degree one. At most one edge is possible, since two
would be an induced2K2. Therefore C has an independent set of size
at least six, a contradiction.

For N=8, C has maximum degree two. At most one component can contain
an edge, since edges in different components induce2K2. Its nontrivial
component is a path on at most four vertices or a cycle on at most
five: a five-vertex path or a cycle of length at least six already has
an induced2K2. In the path cases, adding all isolated vertices gives
an independent set of size at least six. Cycles of lengths3,4,5
give independent sets of sizes6,6,5, respectively. Again impossible.

It remains to prove the following exact nine-vertex lemma.

**Lemma.** Every induced2K2-free graph C on nine vertices with maximum
degree at most three has an independent set of size at least five.

If C contains a triangle T, its vertices have only three degree-slots
outside T in total. Partition the other six vertices into D, those
having a neighbor in T, and S, those having none. Thus |D|<=3 and
|S|>=3. S is independent: any edge in S together with a triangle edge
would induce2K2. If |S|>=4, add any triangle vertex to get five.
Otherwise |S|=|D|=3, and the three T-to-D contacts form a perfect
matching: each triangle vertex has one slot and each D vertex needs
one contact. There are no D-to-S edges, for an edge from d to s
together with the triangle edge avoiding d's unique T-neighbor would
induce2K2. Hence the three S vertices are isolated. A triangle vertex
and a D vertex matched to a DIFFERENT triangle vertex form an independent
pair; adding S gives the required independent five.

Suppose C is trianglefree. If it is bipartite, an independent part has
at least five vertices. Otherwise a shortest odd cycle is induced.
It must have length five, because every induced cycle of length at
least six contains an induced2K2. Denote this cycle by T. Each of its
five vertices has at most one neighbor outside T, so the total number
of contacts from the other four vertices to T is at most five.

Trianglefreeness allows an outside vertex at most two T-neighbors,
which must be nonadjacent. Exactly one T-neighbor is impossible: if
that neighbor is t_i, the edges d-t_i and t_(i+2)-t_(i+3) induce2K2.
Thus every outside vertex has zero or two T-neighbors, and at least
two of the four outside vertices have zero. Each zero-contact vertex
is isolated. For if d-e were an outside edge with d having no T
contact, e's at most two T-neighbors do not cover all cycle edges;
an uncovered cycle edge together with d-e would induce2K2.

If there are at least three such isolated vertices, add two independent
cycle vertices. If there are exactly two, choose another outside
vertex e. Its two T-neighbors are nonadjacent, and deleting them from
C5 leaves three vertices containing an independent pair. That pair,
e, and the two isolated vertices are independent. This proves the
lemma in all cases.

Apply the lemma at N=9, where(4) gives maximum_degree(C)<=3. It
contradicts the absence of an independent five. This exhausts every
N allowed by(5), proving W=0 or W>=15.

## Scope and remaining obstruction

The W=15 endpoint is NOT excluded here. In particular the proof does
not classify graphs on ten or more positive facets when degree-two
vertices are absent. A degree-two facet at W=15 would force equality
in(3), but that separate endpoint analysis is not needed for the theorem.

At or below the120-cell volume, a distinct polytope must therefore have
W>=15 and all facets with at most46 vertices. Strictly below the double,
the known local cutoff remains96. The proof does not establish a global
facet minimum, f+W>=120, a double-volume gap, a garland theorem, or a
new dimension ceiling. It is a paper proof relative to the explicit
inherited lemmas and3D metric inputs above, with no execution or
end-to-end formal verification claim.
