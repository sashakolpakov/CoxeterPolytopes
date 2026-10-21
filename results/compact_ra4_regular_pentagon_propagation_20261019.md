# Compact RA4: regular-pentagon propagation forces excess at least nine

19 October2026. New written proof, awaiting independent review. No
scientific calculation, enumeration, optimizer, build or Lean ran.
Earlier pinned sources and the compact-capacity manuscript are unchanged.
This note strengthens the positive-excess obstruction; it does not settle
the compact minimum-volume or double-volume questions.

## Inputs and statement

Let P be a compact right-angled hyperbolic four-polytope. For each facet F
and ridge R put

    d(F)=sum_(R contained in F)(number_of_sides(R)-5),
    W=sum_R(number_of_sides(R)-5).

Then d(F)=number_of_faces(F)-12 and sum_F d(F)=2W. A zero-excess
facet is the regular right-angled dodecahedron, using the explicit3D
metric classification already recorded in the reviewed compact packet.
Every positive-excess facet has d>=2: a k-gonal ridge with k>=6 forces
at least2k+2 faces and hence d>=2(k-5)>=2.

Use these reviewed own geometric/combinatorial results:

* a compact right-angled3-polyhedron with a k-gon has at least2k+2
  faces; at equality it has Loebell combinatorics, with two disjoint
  k-gons and2k pentagons;
* the positive-excess graph G, with an edge of weight k-5 for each
  nonpentagonal shared ridge, has no degree-one vertex; at a degree-two
  vertex the two weights are equal and the corresponding facet is that
  Loebell polyhedron;
* orthogonal actual outward normals determine an actual face, and the
  full facet graph has no induced square and no clique of size five.

These are in `compact_ra4_neighbor_excess_and_zero_excess_rigidity_20261019.md`,
SHA256 `9c670a7b3c25dda48997c3f70a6d1be9be88bcdefb475b7710f5011068aebcbe`,
and its root review, SHA256
`292187dbf4fb5ea327baed6086177ec1b0fac174116ca3c579a6d75bc43f0baa`.

One further input is the compact3D Andreev metric-uniqueness clause:
the marked combinatorial type and its dihedral angles determine its
compact nonobtuse realization up to isometry. The precise clause is
in Roeder--Hubbard--Dunbar,
[Theorem1.4, p832](https://aif.centre-mersenne.org/article/AIF_2007__57_3_825_0.pdf#page=9).
We use it only to realize the combinatorial rotational symmetry of a
Loebell facet, not to classify four-polytopes.

**Theorem.** Every compact right-angled four-polytope satisfies

    W=0 or W>=9.                                     (1)

The reviewed zero-excess recognition identifies W=0 with the regular
right-angled120-cell. Every distinct polytope therefore has W>=9.

## 1. A degree-two positive facet costs at least fifteen globally

Suppose F has exactly two nonpentagonal ridges. If their edge weights
are w_1,w_2, the3D capacity inequalities give

    d(F)=w_1+w_2>=2 max(w_1,w_2).

Thus w_1=w_2=w>=1. Write k=w+5. Equality in the closing proof gives
the Loebell type with two k-gonal bases and2k lateral pentagons.

The bases are regular right-angled k-gons. Indeed the combinatorial
order-k rotation preserving both bases and rotating the two belts
extends, by the stated metric uniqueness, to an order-k isometry.
Its restriction to either convex base is the one-step rotational
symmetry. Thus all its sides and angles agree. If L_k is its side
length, the right triangle from center to vertex to edge midpoint has
angles pi/k,pi/4,pi/2, giving

    cosh(L_k/2)=sqrt(2)*cos(pi/k),
    cosh(L_k)=1+2*cos(2pi/k).                         (2)

For the regular right-angled pentagon let L=L_5, so cosh(L)=phi,
phi=(1+sqrt(5))/2. Since k>=6, equation(2) gives L_k>L.

Every lateral pentagon of F contains a whole base edge of length L_k.
It is therefore NOT a regular right-angled pentagon. The other facet
of P across this pentagon cannot be a zero-excess facet: all faces of
the latter's regular dodecahedron are regular pentagons of edge length L.
Each of the2k lateral pentagons consequently has a positive-excess
neighbor, of excess at least2.

These2k neighbors are distinct from one another and from the two base
neighbors: two facets of a convex polytope have only one intersection
face, and cannot share two different ridges. The two base neighbors
each have excess at least2(k-5) by3D capacity. Summing these DISTINCT
facets, along with F, against sum d=2W gives

    2W >= 2(k-5)+2*2(k-5)+(2k)*2 = 10k-30,
    W >= 5k-15 >= 15.                               (3)

The argument does not presume that F's pentagonal neighbors have no
further nonpentagonal ridges, or that the selected excesses occur on
disjoint ridges. Only distinct facets are summed in the exact identity.

In particular, if 0<W<15, every vertex of G has ordinary degree at
least three. This is stronger than an abstract degree-two cycle test.

## 2. Four consecutive regular-pentagon edges force closure

**Polygon lemma.** Let R be a compact right-angled hyperbolic polygon.
If it has four consecutive edges of length L, then R is the regular
right-angled pentagon.

Match their ordered, inward-oriented four-edge chain to four consecutive
edges of a regular right-angled pentagon. Such a match is unique up to
isometry: the first directed edge and the successive fixed lengths and
right-angle inward turns determine the entire chain. At either end of
the matched chain, the next supporting line is the unique perpendicular
to its end edge through the corresponding endpoint. In the regular
pentagon those two lines are the SAME fifth supporting line. They are
therefore also the same line for R. A compact convex polygon has only
one edge in one supporting line; the connecting segment is the fifth
edge and closes the boundary. There is no room for an additional edge
between them. This proves the lemma, including its side-orientation
condition. It uses no hexagon-length formula or polygon classification.

**Application.** Take a nonpentagonal ridge R=F intersect F', with
k>=6 sides. Each of its edges is carried by a distinct third facet H.
If H has zero excess, the face F intersect H in H is a regular
pentagon. The shared edge is an entire edge of this pentagon, so has
length L. Consequently R cannot have four consecutive edges carried
by zero-excess facets.

If R has exactly two positive-excess edge carriers H_1,H_2 which are
adjacent in R, its other k-2>=4 edges form a consecutive zero-excess
block, contradicting the lemma. Orthogonality of the four actual normals
of F,F',H_1,H_2 is sufficient to give this adjacency: their actual
common vertex lies in R, by the reviewed projection lemma.

## 3. The small full positive-facet graph cannot close

Suppose toward a contradiction that 0<W<=8. By(3), G has minimum
degree at least three. If E is its unweighted edge count, then E<=W<=8,
so the number N of positive-excess facets is at most floor(2E/3)<=5.
There are at least four such facets.

Let H denote the FULL adjacency graph on these N positive-excess
facets, including pentagonal shared ridges. Thus G is a spanning
subgraph of H, and H also has minimum degree at least three.

If N=4 then H=K4. Choose any edge of G and let R be its nonpentagonal
ridge. The two remaining positive facets both meet R, by orthogonal
projection. They give two adjacent edge carriers of R because all
four normals are orthogonal. There are no other positive facets.
The application in Section2 gives a contradiction.

If N=5, the complement of H has maximum degree at most one. Hence it
is a matching of size zero, one or two. Size zero would make five
mutually orthogonal positive unit normals, impossible in signature(4,1).
Size two would leave an induced square on its four endpoints, also
impossible by the acute-normal index argument. Therefore

    H=K5 minus one edge uv.

Write a,b,c for the other three vertices. Vertex u has exactly three
neighbors in H; its degree in G is at least three, so ua,ub,uc all
belong to G. In particular R=F_u intersect F_a is nonpentagonal.

The only other positive facets meeting R are F_b,F_c: they meet it
by orthogonal projection, and F_v is disjoint from F_u. The two edges
of R carried by F_b,F_c are adjacent since the four normals for
u,a,b,c are orthogonal. All other edges of R have zero-excess
carriers. Section2 again gives a contradiction.

This exhausts N<=5 and proves(1). Importantly, it uses FULL adjacency
H to deal with a possible pentagonal diagonal of G; an induced square
in G alone would not suffice.

## Consequences and limits

Every distinct polytope at or below the120-cell volume must have W>=9
and all facets with at most46 vertices. A counterexample below the
double still has all facets at most96 vertices. These conclusions do
not give f+W>=120 or a universal double-volume gap. There is no
classification of the possible positive graphs with six or more
vertices, no proof of a global120-cell tessellation, and no dimension
improvement. This successor theorem currently has written proof only,
with the two explicit3D metric inputs and inherited acute/capacity
lemmas above; independent review and formalization remain separate.
